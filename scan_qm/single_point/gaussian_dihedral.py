#!/usr/bin/env python3
"""Build a torsional energy profile from a set of Gaussian output files.

For every output file the script measures a dihedral angle on the last geometry
and reads the last SCF energy, then prints the two columns sorted by angle, with
the energies relative to the lowest point of the set.

Examples
--------
  # profile of the 10-1-4-7 dihedral, in kcal/mol
  gaussian_dihedral.py 10 1 4 7 sp_*.log -o qm_scan.dat

  # the same in kJ/mol, with the energies relative to the first file
  gaussian_dihedral.py 10 1 4 7 sp_*.log --units kJ --reference first

  # read the MP2 energy instead of the SCF one
  gaussian_dihedral.py 10 1 4 7 sp_*.log --energy mp2
"""

import argparse
import math
import re
import sys

HARTREE = {"kcal": 627.509474, "kJ": 2625.499639, "hartree": 1.0}

# Energy lines Gaussian writes, and how to pull the number out of them.
ENERGY_PATTERNS = {
    "scf": (re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)"), None),
    "mp2": (re.compile(r"EUMP2\s*=\s*(-?[\d.]+D[+-]?\d+)"), "fortran"),
    "oniom": (re.compile(r"ONIOM:\s+extrapolated energy\s*=\s*(-?\d+\.\d+)"), None),
    "ccsd(t)": (re.compile(r"CCSD\(T\)=\s*(-?[\d.]+D[+-]?\d+)"), "fortran"),
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("atoms", nargs=4, type=int,
                   help="the four 1-based atom indices that define the dihedral")
    p.add_argument("files", nargs="+", help="Gaussian output files (.log or .out)")
    p.add_argument("-o", "--output", default=None,
                   help="write the profile to this file instead of the screen")
    p.add_argument("--units", choices=["kcal", "kJ", "hartree"], default="kcal",
                   help="energy units of the second column (default: %(default)s)")
    p.add_argument("--energy", choices=sorted(ENERGY_PATTERNS), default="scf",
                   help="energy to read from each file (default: %(default)s)")
    p.add_argument("--orientation", choices=["input", "standard"], default="standard",
                   help="orientation block the geometry is read from (default: %(default)s)")
    p.add_argument("--reference", choices=["min", "first", "none"], default="min",
                   help="energy the profile is measured from (default: %(default)s)")
    p.add_argument("--range", choices=["180", "360"], default="180",
                   help="report angles in (-180,180] or [0,360) (default: %(default)s)")
    p.add_argument("--precision", type=int, default=6,
                   help="decimals of the energy column (default: %(default)s)")
    return p.parse_args(argv)


def to_float(text, style):
    """Gaussian writes some energies in Fortran D notation."""
    return float(text.replace("D", "E")) if style == "fortran" else float(text)


def read_energy(lines, kind):
    pattern, style = ENERGY_PATTERNS[kind]
    value = None
    for line in lines:
        m = pattern.search(line)
        if m:
            value = to_float(m.group(1), style)
    return value


def read_geometry(lines, orientation):
    tag = "Input orientation:" if orientation == "input" else "Standard orientation:"
    start = None
    for i, line in enumerate(lines):
        if tag in line:
            start = i
    if start is None:
        return []
    geom = []
    i = start + 5
    while i < len(lines) and "----" not in lines[i]:
        f = lines[i].split()
        if len(f) >= 6:
            geom.append(tuple(float(v) for v in f[-3:]))
        i += 1
    return geom


def dihedral(p1, p2, p3, p4):
    """Dihedral angle in degrees, in the convention used by Gaussian and GROMACS.

    This is the standard formulation: the two end bonds are projected onto the
    plane normal to the central bond and the signed angle between them is taken.
    A formula that looks similar but returns the opposite sign is easy to write
    by accident, so the result is worth checking once against the value the QM
    program prints for the same four atoms.
    """
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def scale(a, f):
        return (a[0] * f, a[1] * f, a[2] * f)

    b0 = sub(p1, p2)
    b1 = sub(p3, p2)
    b2 = sub(p4, p3)
    b1 = scale(b1, 1.0 / math.sqrt(dot(b1, b1)))
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def wrap(angle, rng):
    if rng == "360":
        return angle % 360.0
    a = (angle + 180.0) % 360.0 - 180.0
    return 180.0 if a == -180.0 else a


def main(argv=None):
    args = parse_args(argv)
    idx = [a - 1 for a in args.atoms]
    if min(idx) < 0:
        raise SystemExit("atom indices are 1-based and must be positive")

    data = []
    for name in args.files:
        with open(name) as fh:
            lines = fh.readlines()
        energy = read_energy(lines, args.energy)
        geom = read_geometry(lines, args.orientation)
        if energy is None:
            print(f"warning: no {args.energy} energy in {name}, skipped", file=sys.stderr)
            continue
        if len(geom) <= max(idx):
            print(f"warning: no geometry (or too few atoms) in {name}, skipped",
                  file=sys.stderr)
            continue
        angle = wrap(dihedral(*(geom[i] for i in idx)), args.range)
        data.append((angle, energy, name))

    if not data:
        raise SystemExit("no usable output file")

    if args.reference == "min":
        ref = min(e for _, e, _ in data)
    elif args.reference == "first":
        ref = data[0][1]
    else:
        ref = 0.0

    factor = HARTREE[args.units]
    data.sort(key=lambda row: row[0])

    unit = {"kcal": "kcal/mol", "kJ": "kJ/mol", "hartree": "hartree"}[args.units]
    out = open(args.output, "w") if args.output else sys.stdout
    try:
        print(f"# dihedral {'-'.join(str(a) for a in args.atoms)}"
              f"   {args.energy.upper()} energy   {len(data)} points", file=out)
        print(f"#{'angle(deg)':>14} {'energy(' + unit + ')':>22}", file=out)
        for angle, energy, _ in data:
            print(f"{angle:15.3f} {(energy - ref) * factor:22.{args.precision}f}", file=out)
    finally:
        if args.output:
            out.close()
            print(f"{len(data)} points -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
