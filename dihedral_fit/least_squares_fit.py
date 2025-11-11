#!/usr/bin/env python3
"""
fit_dihedral_gromacs.py

Fit Fourier (cos/sin) series to a target dihedral energy profile (QM - MM_zeroed)
and produce dihedral terms suitable for use in GROMACS (multiplicity, Vn, phase).

Usage:
    python fit_dihedral_gromacs.py --qm qm.txt --mm mm_zeroed.txt --out dihedral.itp --nmax 6

Input file format (text, whitespace or comma separated):
    angle(deg)   energy(kJ/mol)

What the script does:
  - Reads QM and MM_zeroed profiles (angles must match or will be interpolated)
  - Builds target energy: E_target = E_QM - E_MM_zeroed
  - Fits E_target(φ) = C0 + sum_{n=1..nmax} [A_n*cos(nφ) + B_n*sin(nφ)] by linear least-squares
  - Converts each (A_n, B_n) to GROMACS form: term = K_n*(1 + cos(nφ - phase)) where
      K_n = sqrt(A_n^2 + B_n^2)
      phase = atan2(B_n, A_n) (radians, converted to degrees)
    Note: Because GROMACS uses K_n*(1+cos(...)) which adds a constant K_n, we absorb that
    into the fitted C0 (constant offset). The produced K_n and phase reproduce the oscillatory part.
  - Writes a simple .itp-like file listing every fitted term and a short summary.

Dependencies: numpy, scipy (scipy optional; if available will perform a small nonlinear refine)

"""

import argparse
import sys
from math import pi, atan2, degrees
import numpy as np

try:
    from scipy.optimize import least_squares
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


def read_profile(path):
    data = np.loadtxt(path, comments="#", delimiter=None)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"File {path} must have at least two columns: angle energy")
    angles = data[:,0]
    energies = data[:,1]
    return angles, energies


def interp_profile(angles_from, energies_from, angles_to):
    # angles in degrees; we'll normalize to [0,360) then interpolate using wrap-around
    a_from = np.mod(angles_from, 360.0)
    idx = np.argsort(a_from)
    a_sorted = a_from[idx]
    e_sorted = energies_from[idx]
    # append point at +360 to allow interpolation across 360->0
    a_ext = np.concatenate([a_sorted, a_sorted + 360.0])
    e_ext = np.concatenate([e_sorted, e_sorted])
    a_query = np.mod(angles_to, 360.0)
    # linear interpolation
    e_query = np.interp(a_query, a_ext, e_ext)
    return e_query


def build_design_matrix(phis_rad, nmax):
    # columns: constant, cos(n phi) and sin(n phi) for n=1..nmax
    N = len(phis_rad)
    cols = [np.ones(N)]
    for n in range(1, nmax+1):
        cols.append(np.cos(n * phis_rad))
        cols.append(np.sin(n * phis_rad))
    return np.column_stack(cols)


def fit_linear(phis_rad, Etarg, nmax):
    X = build_design_matrix(phis_rad, nmax)
    coeffs, *_ = np.linalg.lstsq(X, Etarg, rcond=None)
    C0 = coeffs[0]
    An = coeffs[1::2]
    Bn = coeffs[2::2]
    return C0, An, Bn


def convert_to_gromacs_terms(An, Bn):
    terms = []
    for i, (a, b) in enumerate(zip(An, Bn), start=1):
        R = np.hypot(a, b)  # amplitude for cos(n phi - phase)
        if R < 1e-12:
            phase_deg = 0.0
        else:
            phase = atan2(b, a)  # radians
            phase_deg = degrees(phase)
        terms.append((i, R, phase_deg))
    return terms


def refine_nonlinear(phis_rad, Etarg, terms_initial):
    # refine amplitudes and phases by fitting sum_k Kk*(1+cos(n phi - phase_k)) + C0
    # Parameter vector: [C0, K1, p1, K2, p2, ...]
    def residuals(p):
        C0 = p[0]
        val = np.full_like(phis_rad, C0)
        for idx, (K, ph_deg) in enumerate(zip(p[1::2], p[2::2])):
            n = idx + 1
            phase = ph_deg
            val += K * (1.0 + np.cos(n * phis_rad - phase))
        return val - Etarg

    # prepare initial parameter vector: convert phase to radians
    p0 = [0.0]
    for n, K, ph_deg in terms_initial:
        p0.append(K)
        p0.append(np.deg2rad(ph_deg))
    p0 = np.array(p0)

    # optionally fix bounds: allow K to be any real (but non-negative better), phase any
    lb = [-np.inf] + [0.0 if i%2==1 else -np.inf for i in range(1, len(p0))]  # attempt K>=0
    ub = [np.inf] + [np.inf]* (len(p0)-1)

    try:
        res = least_squares(residuals, p0, bounds=(lb, ub), verbose=0)
    except Exception as e:
        print("scipy least_squares failed or not available. Skipping nonlinear refine:", e, file=sys.stderr)
        return terms_initial, None

    p_opt = res.x
    C0_opt = p_opt[0]
    terms_opt = []
    for i in range(1, len(p_opt), 2):
        K = p_opt[i]
        ph = p_opt[i+1]
        n = (i+1)//2
        terms_opt.append((n, K, degrees(ph)))
    return terms_opt, C0_opt


def write_itp(outfile, terms, C0, header_comment=None):
    with open(outfile, 'w') as f:
        if header_comment:
            f.write(f"; {header_comment}\n")
        f.write("; Dihedral terms produced by fit_dihedral_gromacs.py\n")
        f.write("; Format (human): multiplicity   K_n[kJ/mol]   phase[deg]   (GROMACS form: K_n*(1+cos(n*phi - phase)))\n")
        f.write("[ dihedraltypes ]\n")
        f.write(";  multiplicity   K_n   phase_deg\n")
        for n, K, phase_deg in terms:
            f.write(f"{n:3d}    {K:12.6f}    {phase_deg:10.4f}\n")
        f.write("\n; Constant offset (C0) fitted:\n")        
        f.write(f"; C0 = {C0:.6f} kJ/mol (this is an overall energy shift, GROMACS dihedral functional form adds +K_n constants for each term)\n")


def main():
    p = argparse.ArgumentParser(description="Fit Fourier series to dihedral energy and produce GROMACS-style terms.")
    p.add_argument('--qm', required=True, help='QM profile file: angle(deg) energy(kJ/mol)')
    p.add_argument('--mm', required=True, help='MM profile with dihedral zeroed: angle energy (same units)')
    p.add_argument('--out', default='dihedral.itp', help='Output itp-like filename')
    p.add_argument('--nmax', type=int, default=6, help='Maximum multiplicity n to fit')
    p.add_argument('--refine', action='store_true', help='Do nonlinear refinement (requires scipy)')
    p.add_argument('--no-offset', action='store_true', help='Force zero constant offset (C0=0)')
    args = p.parse_args()

    # Read profiles (energies in kcal/mol), then convert to kJ/mol
    aq, Eq = read_profile(args.qm)
    am, Em = read_profile(args.mm)
    
    Eq *= 4.184
    Em *= 4.184

    # define common angle grid: use QM angles
    angles = np.mod(aq, 360.0)
    Em_interp = interp_profile(am, Em, angles)

    Etarg = Eq - Em_interp
    phis_rad = np.deg2rad(angles)

    C0, An, Bn = fit_linear(phis_rad, Etarg, args.nmax)

    if args.no_offset:
        C0 = 0.0

    terms = convert_to_gromacs_terms(An, Bn)

    # optionally refine
    if args.refine:
        if not SCIPY_AVAILABLE:
            print("scipy not available; cannot refine. Install scipy to use --refine.", file=sys.stderr)
        else:
            # prepare initial tuple list as (n, K, phase_deg)
            terms_initial = terms
            terms_opt, C0_opt = refine_nonlinear(phis_rad, Etarg, terms_initial)
            if terms_opt is not None:
                terms = terms_opt
                C0 = C0_opt if C0_opt is not None else C0

    # compute fitted energy and RMSD
    fitted = np.zeros_like(Etarg) + C0
    for n, K, ph_deg in terms:
        fitted += K * (1.0 + np.cos(n * phis_rad - np.deg2rad(ph_deg)))
    rmsd = np.sqrt(np.mean((Etarg - fitted)**2))

    header = f"Fitted up to n={args.nmax}, RMSD={rmsd:.6f} kJ/mol"
    write_itp(args.out, terms, C0, header_comment=header)

    print("Wrote:", args.out)
    print(header)
    print("Terms (multiplicity, K[kJ/mol], phase_deg):")
    for n, K, ph in terms:
        print(f"  n={n:2d}    K={K:10.6f} kJ/mol    phase={ph:9.4f} deg")
    print(f"Constant offset C0 = {C0:.6f} kJ/mol")
    print("\nNotes:")
    print(" - GROMACS dihedral functional often used: K_n*(1 + cos(n*phi - phase))")
    print(" - The script reports K_n and phase (degrees). Each term also contributes +K_n constant to the total energy; C0 is the residual offset.")
    print(" - You can paste the produced terms into an .itp and then reference them in your topology. Verify units (this script assumes kJ/mol for energies).")

if __name__ == '__main__':
    main()

