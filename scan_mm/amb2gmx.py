#!/usr/bin/env python3
"""Convert an AMBER topology and its coordinates to the GROMACS format.

The conversion is done by parmed, which keeps the parameters as they are and
only rewrites them: AMBER and GROMACS use the same functional form for proper
dihedrals, so the heights become kJ/mol and the multiplicities and phases are
carried over unchanged.

Examples
--------
  amb2gmx.py propanediol.prmtop propanediol.rst7

  amb2gmx.py propanediol.prmtop propanediol.rst7 \\
      --top propanediol.top --gro propanediol.gro --box 10 10 10
"""

import argparse
import os
import sys

try:
    import parmed as pmd
except ImportError:
    raise SystemExit("this script needs parmed: pip install parmed")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("prmtop", help="AMBER topology (.prmtop or .parm7)")
    p.add_argument("coordinates", help="AMBER coordinates (.rst7 or .inpcrd)")
    p.add_argument("--top", default=None,
                   help="GROMACS topology to write (default: <prmtop>_converted.top)")
    p.add_argument("--gro", default=None,
                   help="GROMACS structure to write (default: <coordinates>_converted.gro)")
    p.add_argument("--box", nargs=3, type=float, default=None, metavar=("X", "Y", "Z"),
                   help="rectangular box in nanometres, for a topology with no box")
    p.add_argument("--combine", default="all", choices=["all", "none"],
                   help="write one [ moleculetype ] for the whole system ('all') or "
                        "one per molecule ('none') (default: %(default)s)")
    p.add_argument("--overwrite", action="store_true",
                   help="replace the output files if they already exist")
    return p.parse_args(argv)


def default_name(path, suffix):
    stem = path.rsplit(".", 1)[0]
    return f"{stem}_converted{suffix}"


def main(argv=None):
    args = parse_args(argv)
    top = args.top or default_name(args.prmtop, ".top")
    gro = args.gro or default_name(args.coordinates, ".gro")

    for path in (top, gro):
        if os.path.exists(path) and not args.overwrite:
            raise SystemExit(f"{path} exists; use --overwrite to replace it")

    try:
        structure = pmd.load_file(args.prmtop, xyz=args.coordinates)
    except Exception as err:
        raise SystemExit(f"could not read {args.prmtop} / {args.coordinates}: {err}")

    if args.box:
        structure.box = [args.box[0] * 10.0, args.box[1] * 10.0, args.box[2] * 10.0,
                         90.0, 90.0, 90.0]
    elif structure.box is None:
        print("note: the topology carries no box; gmx editconf or --box can add one",
              file=sys.stderr)

    structure.save(top, format="gromacs", combine=args.combine, overwrite=True)
    structure.save(gro, overwrite=True)
    print(f"{len(structure.atoms)} atoms -> {top}, {gro}")


if __name__ == "__main__":
    main()
