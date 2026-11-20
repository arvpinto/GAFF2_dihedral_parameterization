#!/usr/bin/env python3
"""Fit a cosine series to the difference between a QM and an MM torsional profile.

The target of the fit is

    E_target(phi) = E_QM(phi) - E_MM0(phi)

where E_MM0 is the molecular mechanics profile computed with the torsion being
parameterised removed from the topology. The series has the functional form used
by both GROMACS (proper dihedral, function 1 or 9) and AMBER,

    E_fit(phi) = C0 + sum_n K_n * [1 + cos(n*phi - delta_n)]

so the fitted K_n and delta_n can be written straight into a topology. C0 is an
overall shift that no force field term carries and that cancels in any relative
energy.

The two profiles do not need to share a grid: the MM profile is interpolated
onto the QM angles, periodically.

Examples
--------
  # three terms, phases free, GROMACS and parmed lines for atoms 10 1 4 7
  least_squares_fit.py --qm qm_scan.dat --mm mm_scan.dat --nmax 3 \\
      --indices 10 1 4 7 --out dihedral.itp --plot fit.png

  # only the multiplicities 1 and 3, phases restricted to 0 or 180 degrees
  least_squares_fit.py --qm qm_scan.dat --mm mm_scan.dat \\
      --terms 1 3 --phases symmetric

  # how much each extra multiplicity buys, before choosing nmax
  least_squares_fit.py --qm qm_scan.dat --mm mm_scan.dat --scan 6
"""

import argparse
import sys

import numpy as np

try:
    from scipy.optimize import least_squares as _least_squares
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

KJ_PER_KCAL = 4.184


# --------------------------------------------------------------------------- io

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--qm", required=True, metavar="FILE",
                   help="QM profile: two columns, angle(deg) and energy")
    p.add_argument("--mm", required=True, metavar="FILE",
                   help="MM profile with the torsion removed, same two columns")
    p.add_argument("--units", choices=["kcal", "kJ"], default="kcal",
                   help="units of the energy column of BOTH input files "
                        "(default: %(default)s)")

    g = p.add_mutually_exclusive_group()
    g.add_argument("--nmax", type=int, default=3,
                   help="fit the multiplicities 1..NMAX (default: %(default)s)")
    g.add_argument("--terms", nargs="+", type=int, metavar="N",
                   help="fit only these multiplicities, e.g. --terms 1 2 3 6")

    p.add_argument("--phases", choices=["free", "symmetric"], default="free",
                   help="'free' fits any phase, 'symmetric' restricts it to 0 or 180 "
                        "degrees, which keeps the term symmetric about phi=0 "
                        "(default: %(default)s)")
    p.add_argument("--refine", action="store_true",
                   help="refine the linear solution with a non-linear least squares "
                        "step (needs scipy; changes little when the grid is regular)")
    p.add_argument("--no-offset", action="store_true",
                   help="force the constant term C0 to zero")
    p.add_argument("--drop-duplicates", action="store_true", default=True,
                   help="keep only the first point of each repeated angle (default)")
    p.add_argument("--keep-duplicates", dest="drop_duplicates", action="store_false",
                   help="keep repeated angles, for instance the two ends of a full turn")

    p.add_argument("--scan", type=int, metavar="NMAX",
                   help="only report the RMSD obtained with 1..NMAX terms and stop")
    p.add_argument("--indices", nargs=4, type=int, metavar="I",
                   help="atom indices written in the GROMACS [ dihedrals ] lines")
    p.add_argument("--types", nargs=4, metavar="T",
                   help="atom types written in the GROMACS [ dihedraltypes ] lines, "
                        "e.g. --types oh c3 c3 c3")
    p.add_argument("--funct", type=int, choices=[1, 9], default=9,
                   help="GROMACS dihedral function type; 9 is the multi-term proper "
                        "dihedral (default: %(default)s)")
    p.add_argument("--out", default=None, metavar="FILE",
                   help="write the topology lines to this file")
    p.add_argument("--plot", default=None, metavar="FILE",
                   help="save a figure with the two profiles, the target and the fit")
    p.add_argument("--dpi", type=int, default=600,
                   help="resolution of the figure (default: %(default)s)")
    p.add_argument("--show", action="store_true", help="also display the figure")
    return p.parse_args(argv)


def read_profile(path, to_kj, drop_duplicates):
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise SystemExit(f"{path}: expected two columns, angle(deg) and energy")
    angles, energies = data[:, 0].astype(float), data[:, 1].astype(float) * to_kj
    if drop_duplicates:
        _, keep = np.unique(np.round(np.mod(angles, 360.0), 3), return_index=True)
        if len(keep) != len(angles):
            print(f"note: {len(angles) - len(keep)} repeated angle(s) dropped from "
                  f"{path}", file=sys.stderr)
        keep = np.sort(keep)
        angles, energies = angles[keep], energies[keep]
    order = np.argsort(angles)
    return angles[order], energies[order]


def interpolate(angles_from, energies_from, angles_to):
    """Periodic linear interpolation of a torsional profile onto new angles."""
    a = np.mod(angles_from, 360.0)
    order = np.argsort(a)
    return np.interp(np.mod(angles_to, 360.0), a[order], energies_from[order], period=360.0)


# ------------------------------------------------------------------------- fit

def linear_fit(phi, target, orders, phases, offset=True):
    """Least squares in the cos/sin basis. Returns C0 and [(n, K, delta_deg), ...]."""
    columns = [np.ones_like(phi)] if offset else []
    for n in orders:
        columns.append(np.cos(n * phi))
        if phases == "free":
            columns.append(np.sin(n * phi))
    design = np.column_stack(columns)
    coeff, *_ = np.linalg.lstsq(design, target, rcond=None)

    k = 0
    c0 = 0.0
    if offset:
        c0 = coeff[0]
        k = 1

    terms = []
    for n in orders:
        a = coeff[k]
        k += 1
        if phases == "free":
            b = coeff[k]
            k += 1
            amplitude = float(np.hypot(a, b))
            delta = float(np.degrees(np.arctan2(b, a)))
        else:
            amplitude = float(abs(a))
            delta = 0.0 if a >= 0 else 180.0
        terms.append((n, amplitude, delta))

    # Each term is written as K_n*[1 + cos(n*phi - delta_n)] and so carries a
    # constant K_n of its own; the fitted C0 has to give that back for the two
    # forms to describe the same curve.
    c0 -= sum(amplitude for _n, amplitude, _delta in terms)
    return float(c0), terms


def evaluate(phi, c0, terms):
    value = np.full_like(phi, c0)
    for n, amplitude, delta in terms:
        value += amplitude * (1.0 + np.cos(n * phi - np.deg2rad(delta)))
    return value


def refine(phi, target, c0, terms, phases):
    """Non-linear refinement of the amplitudes and phases."""
    if not HAVE_SCIPY:
        print("note: scipy is not installed, --refine ignored", file=sys.stderr)
        return c0, terms

    orders = [n for n, _, _ in terms]
    free = phases == "free"
    p0 = [c0]
    for _n, amplitude, delta in terms:
        p0.append(amplitude)
        if free:
            p0.append(np.deg2rad(delta))
        else:
            # a signed amplitude covers the phases 0 and 180 degrees
            p0[-1] = amplitude if delta < 90.0 else -amplitude

    def unpack(p):
        out = []
        k = 1
        for n in orders:
            if free:
                out.append((n, p[k], float(np.degrees(p[k + 1]))))
                k += 2
            else:
                out.append((n, abs(p[k]), 0.0 if p[k] >= 0 else 180.0))
                k += 1
        return float(p[0]), out

    def residual(p):
        c, t = unpack(p)
        return evaluate(phi, c, t) - target

    result = _least_squares(residual, np.array(p0, dtype=float))
    return unpack(result.x)


def quality(phi, target, c0, terms):
    residual = evaluate(phi, c0, terms) - target
    return float(np.sqrt(np.mean(residual ** 2))), float(np.max(np.abs(residual)))


# ---------------------------------------------------------------------- output

def topology_lines(terms, args):
    lines = []
    lines.append("; proper dihedral fitted with least_squares_fit.py")
    lines.append("; E = sum_n K_n * [1 + cos(n*phi - delta_n)]")
    if args.types:
        lines.append("[ dihedraltypes ]")
        lines.append(";   i    j    k    l  func      delta(deg)       K(kJ/mol)   n")
        for n, amplitude, delta in terms:
            lines.append(f"  {args.types[0]:>3s}  {args.types[1]:>3s}  {args.types[2]:>3s}"
                         f"  {args.types[3]:>3s}   {args.funct:3d}  {delta:14.5f}"
                         f"  {amplitude:14.6f}  {n:3d}")
        lines.append("")
    if args.indices:
        i, j, k, l = args.indices
        lines.append("[ dihedrals ]")
        lines.append(";   ai     aj     ak     al  func      delta(deg)       K(kJ/mol)   n")
        for n, amplitude, delta in terms:
            lines.append(f"{i:6d} {j:6d} {k:6d} {l:6d}   {args.funct:3d}"
                         f"  {delta:14.5f}  {amplitude:14.6f}  {n:3d}")
        lines.append("")
    if not args.types and not args.indices:
        lines.append("; give --indices and/or --types to get lines ready to paste")
        lines.append(";    n      K(kJ/mol)       delta(deg)")
        for n, amplitude, delta in terms:
            lines.append(f"; {n:4d}  {amplitude:14.6f}  {delta:14.5f}")
    return lines


def parmed_lines(terms, args):
    mask = " ".join(f"@{i}" for i in args.indices) if args.indices else "@i @j @k @l"
    lines = ["# parmed, AMBER convention: the height is in kcal/mol",
             "#   addDihedral <mask1> <mask2> <mask3> <mask4> <K> <n> <delta> <scee> <scnb>"]
    for n, amplitude, delta in terms:
        lines.append(f"addDihedral {mask} {amplitude / KJ_PER_KCAL:.6f} {n} "
                     f"{delta:.4f} 1.2 2.0")
    return lines


def make_plot(path, angles, qm, mm, target, fitted, rmsd, dpi, show):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = np.argsort(angles)
    a = angles[order]
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(7.0, 7.4), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 0.8], "hspace": 0.12},
        facecolor="white",
    )
    for axis in (top, bottom):
        axis.set_facecolor("white")
        axis.grid(True, color="0.88", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_xlim(-180, 180)
        axis.set_xticks(range(-180, 181, 60))

    top.plot(a, (qm - qm.min())[order], "o-", color="black", markersize=3.5,
             linewidth=1.4, label="QM")
    top.plot(a, (mm - mm.min())[order], "s--", color="#d62728", markersize=3.5,
             linewidth=1.4, label="MM, torsion removed")
    top.set_ylabel("relative energy (kJ/mol)")
    top.legend(frameon=False, fontsize=9, loc="lower left",
               bbox_to_anchor=(0.0, 1.01), borderaxespad=0.0, ncol=2)

    bottom.plot(a, (target - target.mean())[order], "o", color="#1f77b4",
                markersize=3.5, label="target, QM - MM")
    bottom.plot(a, (fitted - target.mean())[order], "-", color="#1f77b4",
                linewidth=1.8, label=f"fit, RMSD = {rmsd:.2f} kJ/mol")
    bottom.set_xlabel("dihedral angle (degrees)")
    bottom.set_ylabel("energy (kJ/mol)")
    bottom.legend(frameon=False, fontsize=9, loc="lower left",
                  bbox_to_anchor=(0.0, 1.01), borderaxespad=0.0, ncol=2)

    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    if show:
        plt.show()
    plt.close(fig)


# ------------------------------------------------------------------------ main

def main(argv=None):
    args = parse_args(argv)
    to_kj = KJ_PER_KCAL if args.units == "kcal" else 1.0

    angles, qm = read_profile(args.qm, to_kj, args.drop_duplicates)
    mm_angles, mm_raw = read_profile(args.mm, to_kj, args.drop_duplicates)
    mm = interpolate(mm_angles, mm_raw, angles)

    target = qm - mm
    phi = np.deg2rad(angles)
    offset = not args.no_offset

    if args.scan:
        print(f"{'nmax':>5} {'RMSD(kJ/mol)':>14} {'max dev':>10}")
        for nmax in range(1, args.scan + 1):
            c0, terms = linear_fit(phi, target, range(1, nmax + 1), args.phases, offset)
            if args.refine:
                c0, terms = refine(phi, target, c0, terms, args.phases)
            rmsd, worst = quality(phi, target, c0, terms)
            print(f"{nmax:5d} {rmsd:14.4f} {worst:10.4f}")
        return

    orders = sorted(set(args.terms)) if args.terms else list(range(1, args.nmax + 1))
    if min(orders) < 1:
        raise SystemExit("multiplicities must be 1 or larger")

    c0, terms = linear_fit(phi, target, orders, args.phases, offset)
    if args.refine:
        c0, terms = refine(phi, target, c0, terms, args.phases)
    if args.no_offset:
        c0 = 0.0
    rmsd, worst = quality(phi, target, c0, terms)

    print(f"{len(angles)} points, {len(orders)} term(s), phases {args.phases}")
    print(f"RMSD {rmsd:.4f} kJ/mol ({rmsd / KJ_PER_KCAL:.4f} kcal/mol), "
          f"largest deviation {worst:.4f} kJ/mol")
    print(f"constant offset C0 = {c0:.4f} kJ/mol (not a force field term)")
    print()
    print(f"{'n':>4} {'K(kJ/mol)':>14} {'K(kcal/mol)':>14} {'delta(deg)':>14}")
    for n, amplitude, delta in terms:
        print(f"{n:4d} {amplitude:14.6f} {amplitude / KJ_PER_KCAL:14.6f} {delta:14.4f}")
    print()

    lines = topology_lines(terms, args) + [""] + parmed_lines(terms, args)
    text = "\n".join(lines) + "\n"
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(f"; RMSD {rmsd:.4f} kJ/mol over {len(angles)} points\n")
            fh.write(text)
        print(f"topology lines -> {args.out}")
    else:
        print(text)

    if args.plot:
        make_plot(args.plot, angles, qm, mm, target, evaluate(phi, c0, terms),
                  rmsd, args.dpi, args.show)
        print(f"figure -> {args.plot} ({args.dpi} dpi)")


if __name__ == "__main__":
    main()
