#!/usr/bin/env python3
"""Extract the optimised geometries of a Gaussian relaxed (ModRedundant) scan.

Each converged step of the scan is written to an XYZ trajectory and, optionally,
to one Gaussian input file per step so that single-point energies can be
recomputed at a higher level of theory, or to one GROMACS .gro file per step so
that the same geometries can start the MM side of the scan.

The script is self-contained: it only needs the Gaussian output file.

Examples
--------
  # XYZ trajectory with every converged scan point
  gaussian_scan_geoms.py propanediol_scan.log -o scan_geoms.xyz

  # the same, plus one single-point input per step
  gaussian_scan_geoms.py propanediol_scan.log -o scan_geoms.xyz \\
      --write-inputs --prefix sp --route "#p M062X/6-311++G(d,p)" \\
      --charge 0 --mult 1 --mem 4GB --nproc 8

  # one .gro per step, named after the dihedral 10-1-4-7 of that geometry
  gaussian_scan_geoms.py propanediol_scan.log --drop-last \\
      --write-gro --gro-template propanediol_converted.gro --dihedral 10 1 4 7

  # drop the last point when the scan closes a full turn on the first geometry
  gaussian_scan_geoms.py propanediol_scan.log -o scan_geoms.xyz --drop-last
"""

import argparse
import math
import os
import sys

# Gaussian prints atomic numbers; only the elements likely to appear in an
# organic fragment are needed, with a numeric fallback for anything else.
ELEMENTS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 26: "Fe", 29: "Cu",
    30: "Zn", 35: "Br", 53: "I",
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("logfile", help="Gaussian output file of the relaxed scan")
    p.add_argument("-o", "--output", default="scan_geoms.xyz",
                   help="XYZ trajectory to write (default: %(default)s)")
    p.add_argument("--orientation", choices=["input", "standard"], default="input",
                   help="orientation block to read; 'input' keeps the frame and the "
                        "atom order of the input file (default: %(default)s)")
    p.add_argument("--drop-last", action="store_true",
                   help="discard the last scan point (a 360 degree scan ends on its "
                        "starting geometry, which would be counted twice)")
    p.add_argument("--write-inputs", action="store_true",
                   help="also write one Gaussian input file per scan point")
    p.add_argument("--prefix", default="sp",
                   help="base name of the input files: PREFIX_0001.com, ... "
                        "(default: %(default)s)")
    p.add_argument("--outdir", default=".",
                   help="directory for the input files (default: %(default)s)")
    p.add_argument("--route", default="#p M062X/6-311++G(d,p)",
                   help="route section of the single-point inputs (default: %(default)s)")
    p.add_argument("--charge", type=int, default=0, help="total charge (default: %(default)s)")
    p.add_argument("--mult", type=int, default=1, help="spin multiplicity (default: %(default)s)")
    p.add_argument("--mem", default=None, help="%%mem line of the inputs, e.g. 4GB")
    p.add_argument("--nproc", type=int, default=None, help="%%nprocshared line of the inputs")
    p.add_argument("--write-gro", action="store_true",
                   help="also write one GROMACS .gro per scan point, which needs "
                        "--gro-template and --dihedral")
    p.add_argument("--gro-template", default=None, metavar="FILE",
                   help=".gro file of the same molecule, in the same atom order, "
                        "used for the atom names and the box")
    p.add_argument("--dihedral", nargs=4, type=int, metavar="I",
                   help="four 1-based atom indices; each .gro is named after the "
                        "value this dihedral has in that geometry")
    p.add_argument("--gro-prefix", default="dihe",
                   help="base name of the .gro files: dihe_63.517.gro "
                        "(default: %(default)s)")
    p.add_argument("-q", "--quiet", action="store_true", help="do not print a summary")
    args = p.parse_args(argv)
    if args.write_gro and not (args.gro_template and args.dihedral):
        p.error("--write-gro needs --gro-template and --dihedral")
    return args


def dihedral(p1, p2, p3, p4):
    """Dihedral angle in degrees, in the convention used by Gaussian and GROMACS."""
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

    b0, b1, b2 = sub(p1, p2), sub(p3, p2), sub(p4, p3)
    b1 = scale(b1, 1.0 / math.sqrt(dot(b1, b1)))
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


def write_gro_files(points, args):
    """One .gro per scan point, named after the dihedral of that geometry.

    The template supplies the title, the atom and residue names and the box;
    only the coordinates change, so the atom order of the template has to match
    the order of the Gaussian input.
    """
    with open(args.gro_template) as fh:
        template = fh.readlines()
    natoms = int(template[1].split()[0])
    if natoms != len(points[0][1]):
        raise SystemExit(f"{args.gro_template} has {natoms} atoms, the scan has "
                         f"{len(points[0][1])}")

    box = [float(v) for v in template[2 + natoms].split()[:3]]

    os.makedirs(args.outdir, exist_ok=True)
    idx = [i - 1 for i in args.dihedral]
    written = []
    for _energy, geom in points:
        angle = dihedral(*(geom[i][1:] for i in idx))
        # Gaussian centres its orientations on the origin; put the molecule in
        # the middle of the box instead, so it does not straddle the boundary
        shift = [box[k] / 2.0 - sum(a[k + 1] for a in geom) / (10.0 * len(geom))
                 for k in range(3)]
        name = os.path.join(args.outdir, f"{args.gro_prefix}_{angle:.3f}.gro")
        with open(name, "w") as fh:
            fh.write(template[0])
            fh.write(template[1])
            for line, (_el, x, y, z) in zip(template[2:2 + natoms], geom):
                # keep the first 20 characters (residue and atom names, index)
                fh.write(f"{line[:20]}{x / 10.0 + shift[0]:8.3f}"
                         f"{y / 10.0 + shift[1]:8.3f}{z / 10.0 + shift[2]:8.3f}\n")
            fh.write(template[2 + natoms])
        written.append(name)
    return written


def read_orientation(lines, start, want_input):
    """Read the coordinate block that starts at `start` ("... orientation:")."""
    geom = []
    # header: "Center  Atomic  Atomic  Coordinates", a units line and two rules
    i = start + 5
    while i < len(lines) and "----" not in lines[i]:
        f = lines[i].split()
        if len(f) >= 6:
            z = int(f[1])
            x, y, zc = (float(v) for v in f[-3:])
            geom.append((ELEMENTS.get(z, str(z)), x, y, zc))
        i += 1
    return geom


def parse_scan(path, orientation="input"):
    """Return [(energy_hartree, geometry), ...] for every converged scan step."""
    with open(path) as fh:
        lines = fh.readlines()

    tag = "Input orientation:" if orientation == "input" else "Standard orientation:"
    points = []
    energy = None
    last_geom = None

    for i, line in enumerate(lines):
        if tag in line:
            last_geom = read_orientation(lines, i, orientation == "input")
        elif "SCF Done:" in line:
            f = line.split()
            try:
                energy = float(f[4])
            except (IndexError, ValueError):
                pass
        elif "Optimization completed" in line:
            if energy is None or not last_geom:
                print(f"warning: incomplete scan step near line {i + 1}", file=sys.stderr)
                continue
            points.append((energy, last_geom))

    if not points:
        raise SystemExit(
            f"{path}: no converged scan step found. Is this the output of an "
            "opt=modredundant scan?"
        )
    return points


def write_xyz(path, points):
    with open(path, "w") as fh:
        for n, (energy, geom) in enumerate(points, start=1):
            fh.write(f"{len(geom)}\n")
            fh.write(f" Scan point {n} scf energy: {energy:.9f}\n")
            for el, x, y, z in geom:
                fh.write(f"{el}\t{x:.6f}     {y:.6f}     {z:.6f}\n")


def write_inputs(points, args):
    os.makedirs(args.outdir, exist_ok=True)
    names = []
    for n, (_energy, geom) in enumerate(points, start=1):
        name = os.path.join(args.outdir, f"{args.prefix}_{n:04d}.com")
        with open(name, "w") as fh:
            if args.mem:
                fh.write(f"%mem={args.mem}\n")
            if args.nproc:
                fh.write(f"%nprocshared={args.nproc}\n")
            fh.write(f"{args.route}\n\n")
            fh.write(f"scan point {n}\n\n")
            fh.write(f"{args.charge} {args.mult}\n")
            for el, x, y, z in geom:
                fh.write(f" {el:<2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
            fh.write("\n")
        names.append(name)
    return names


def main(argv=None):
    args = parse_args(argv)
    points = parse_scan(args.logfile, args.orientation)
    if args.drop_last and len(points) > 1:
        points = points[:-1]

    write_xyz(args.output, points)
    written = write_inputs(points, args) if args.write_inputs else []
    gro = write_gro_files(points, args) if args.write_gro else []

    if not args.quiet:
        lo = min(e for e, _ in points)
        hi = max(e for e, _ in points)
        print(f"{len(points)} converged scan points -> {args.output}")
        print(f"SCF energy range: {(hi - lo) * 627.509474:.3f} kcal/mol")
        if written:
            print(f"{len(written)} Gaussian inputs: {written[0]} ... {written[-1]}")
        if gro:
            print(f"{len(gro)} GROMACS structures: {gro[0]} ... {gro[-1]}")


if __name__ == "__main__":
    main()
