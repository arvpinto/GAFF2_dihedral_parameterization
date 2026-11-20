#!/usr/bin/env python3
"""Build a torsional energy profile from a set of GROMACS single-point runs.

The potential energy is read from the mdrun log of each point. The dihedral
angle is taken from the file name by default (dihe_-63.517.log -> -63.517), or
measured on the corresponding structure file when four atom indices are given,
which is the safer option whenever the name may carry other numbers.

Examples
--------
  # angle from the file name, energies in kcal/mol
  gromacs_dihedral.py dihe_*.log -o mm_scan.dat

  # angle measured on the .gro files that sit next to the logs
  gromacs_dihedral.py dihe_*.log --atoms 10 1 4 7 --coords . -o mm_scan.dat

  # keep the GROMACS units and do not shift the profile
  gromacs_dihedral.py dihe_*.log --units kJ --reference none
"""

import argparse
import math
import os
import re
import sys

KJ_PER_KCAL = 4.184
ENERGY_RE = re.compile(r"Potential Energy\s*=\s*(-?[\d.]+[eEdD]?[+-]?\d*)")
NUMBER_RE = re.compile(r"[-+]?\d+\.\d+|[-+]?\d+")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("files", nargs="+", help="GROMACS log files of the single-point runs")
    p.add_argument("-o", "--output", default=None,
                   help="write the profile to this file instead of the screen")
    p.add_argument("--units", choices=["kcal", "kJ"], default="kcal",
                   help="energy units of the second column (default: %(default)s)")
    p.add_argument("--atoms", nargs=4, type=int, metavar="I",
                   help="measure the dihedral on the structure files instead of "
                        "reading it from the file name (1-based atom indices)")
    p.add_argument("--coords", default=None,
                   help="directory holding the .gro files (default: next to each log)")
    p.add_argument("--coord-ext", default=".gro",
                   help="extension of the structure files (default: %(default)s)")
    p.add_argument("--name-regex", default=r"([-+]?\d+\.?\d*)",
                   help="regular expression whose first group holds the angle in the "
                        "file name (default: %(default)s)")
    p.add_argument("--reference", choices=["min", "first", "none"], default="min",
                   help="energy the profile is measured from (default: %(default)s)")
    p.add_argument("--range", choices=["180", "360"], default="180",
                   help="report angles in (-180,180] or [0,360) (default: %(default)s)")
    p.add_argument("--precision", type=int, default=6,
                   help="decimals of the energy column (default: %(default)s)")
    return p.parse_args(argv)


def read_energy(path):
    """Potential energy in kJ/mol from a GROMACS log, or None."""
    value = None
    with open(path) as fh:
        for line in fh:
            m = ENERGY_RE.search(line)
            if m:
                value = float(m.group(1).replace("d", "e").replace("D", "E"))
    return value


def angle_from_name(path, pattern):
    stem = os.path.basename(path)
    stem = stem[: stem.rfind(".")] if "." in stem else stem
    m = re.search(pattern, stem)
    if not m:
        raise ValueError(f"no angle in the name {stem!r}")
    return float(m.group(1))


def read_gro(path):
    """Coordinates of a .gro file, in nanometres."""
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[1].split()[0])
    coords = []
    for line in lines[2:2 + n]:
        # fixed-width gro record: 20 characters of labels, then 3 x 8 characters
        coords.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    return coords


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


def coord_path(log, args):
    stem = os.path.basename(log)
    stem = stem[: stem.rfind(".")] if "." in stem else stem
    folder = args.coords if args.coords else os.path.dirname(log) or "."
    return os.path.join(folder, stem + args.coord_ext)


def main(argv=None):
    args = parse_args(argv)
    idx = [a - 1 for a in args.atoms] if args.atoms else None

    data = []
    for log in args.files:
        energy = read_energy(log)
        if energy is None:
            print(f"warning: no potential energy in {log}, skipped", file=sys.stderr)
            continue
        try:
            if idx is None:
                angle = angle_from_name(log, args.name_regex)
            else:
                coords = read_gro(coord_path(log, args))
                if len(coords) <= max(idx):
                    raise ValueError(f"only {len(coords)} atoms in {coord_path(log, args)}")
                angle = dihedral(*(coords[i] for i in idx))
        except (OSError, ValueError) as err:
            print(f"warning: {log}: {err}, skipped", file=sys.stderr)
            continue
        data.append((wrap(angle, args.range), energy, log))

    if not data:
        raise SystemExit("no usable log file")

    seen = {}
    for angle, _, log in data:
        seen.setdefault(round(angle, 3), []).append(log)
    for angle, logs in seen.items():
        if len(logs) > 1:
            print(f"warning: {len(logs)} files at {angle} deg ({', '.join(logs)})",
                  file=sys.stderr)

    if args.reference == "min":
        ref = min(e for _, e, _ in data)
    elif args.reference == "first":
        ref = data[0][1]
    else:
        ref = 0.0

    factor = 1.0 if args.units == "kJ" else 1.0 / KJ_PER_KCAL
    unit = "kJ/mol" if args.units == "kJ" else "kcal/mol"
    data.sort(key=lambda row: row[0])

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        source = "measured on the structures" if idx else "read from the file names"
        print(f"# GROMACS potential energy, angles {source}, {len(data)} points", file=out)
        print(f"#{'angle(deg)':>14} {'energy(' + unit + ')':>22}", file=out)
        for angle, energy, _ in data:
            print(f"{angle:15.3f} {(energy - ref) * factor:22.{args.precision}f}", file=out)
    finally:
        if args.output:
            out.close()
            print(f"{len(data)} points -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
