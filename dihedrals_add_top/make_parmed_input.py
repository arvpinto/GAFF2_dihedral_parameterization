#!/usr/bin/env python3
"""Write a parmed input that swaps a torsion for a fitted set of terms.

The script finds every dihedral of a prmtop whose four atom types match the
fragment that was parameterised, and writes the deleteDihedral and addDihedral
commands that replace them. Heights are converted from kJ/mol to the kcal/mol
that AMBER topologies use, which is the step most easily forgotten.

Both files use the same functional form, E = K * [1 + cos(n*phi - delta)], so
the multiplicities and the phases carry over unchanged.

Examples
--------
  # terms taken from the .itp written by least_squares_fit.py
  make_parmed_input.py propanediol_system.prmtop --types oh c3 c3 c3 \\
      --itp ../dihedral_fit/dihedral.itp -o dihedrals_parmed.in \\
      --parmout propanediol_system_newdihe.prmtop

  # terms given on the command line, as K(kJ/mol):n:delta(deg)
  make_parmed_input.py propanediol_system.prmtop --types oh c3 c3 c3 \\
      --term 4.591222:1:70.6111 --term 1.345460:2:-12.9995 \\
      --term 0.295361:3:100.3075 -o dihedrals_parmed.in

  # list the matching dihedrals and stop
  make_parmed_input.py propanediol_system.prmtop --types oh c3 c3 c3 --list
"""

import argparse
import re
import sys

try:
    import parmed as pmd
except ImportError:
    raise SystemExit("this script needs parmed: pip install parmed")

KJ_PER_KCAL = 4.184


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("prmtop", help="AMBER topology of the full system")
    p.add_argument("--types", nargs=4, required=True, metavar="T",
                   help="the four atom types of the torsion, e.g. oh c3 c3 c3")
    p.add_argument("--itp", default=None, metavar="FILE",
                   help="file written by least_squares_fit.py, read for the terms")
    p.add_argument("--term", action="append", default=[], metavar="K:N:DELTA",
                   help="one term as height(kJ/mol):multiplicity:phase(deg); "
                        "repeat the option for every term")
    p.add_argument("--units", choices=["kJ", "kcal"], default="kJ",
                   help="units of the heights that are read in (default: %(default)s)")
    p.add_argument("--scee", type=float, default=1.2,
                   help="1-4 electrostatic scaling of the new terms (default: %(default)s)")
    p.add_argument("--scnb", type=float, default=2.0,
                   help="1-4 van der Waals scaling of the new terms (default: %(default)s)")
    p.add_argument("-o", "--output", default="dihedrals_parmed.in",
                   help="parmed input to write (default: %(default)s)")
    p.add_argument("--parmout", default=None, metavar="FILE",
                   help="topology parmed should write at the end")
    p.add_argument("--list", action="store_true",
                   help="only print the matching dihedrals and stop")
    return p.parse_args(argv)


def read_itp(path):
    """Pick (K, n, delta) out of the GROMACS lines written by the fitting script."""
    terms = []
    pattern = re.compile(
        r"^\s*(?:\S+\s+){4}\s*(?:1|9)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+)\s*$"
    )
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith((";", "#", "[")):
                continue
            m = pattern.match(line)
            if m:
                delta, height, n = float(m.group(1)), float(m.group(2)), int(m.group(3))
                terms.append((height, n, delta))
    if not terms:
        raise SystemExit(f"{path}: no dihedral line found")
    # [ dihedraltypes ] and [ dihedrals ] hold the same terms twice
    unique = []
    for t in terms:
        if t not in unique:
            unique.append(t)
    return unique


def parse_terms(args):
    terms = read_itp(args.itp) if args.itp else []
    for text in args.term:
        parts = text.split(":")
        if len(parts) != 3:
            raise SystemExit(f"--term {text!r}: expected height:multiplicity:phase")
        terms.append((float(parts[0]), int(parts[1]), float(parts[2])))
    if not terms:
        raise SystemExit("give the terms with --itp or with one or more --term")
    factor = 1.0 / KJ_PER_KCAL if args.units == "kJ" else 1.0
    return [(height * factor, n, delta) for height, n, delta in sorted(terms, key=lambda t: t[1])]


def matching_dihedrals(structure, types):
    """Every proper dihedral with these four atom types.

    A dihedral and its reverse describe the same torsion, so both orderings
    count as a match and only one of the two is written out. In a symmetric
    fragment such as 1,3-propanediol this finds two paths per molecule.
    """
    wanted = (tuple(types), tuple(reversed(types)))
    found = []
    for d in structure.dihedrals:
        if d.improper:
            continue
        atoms = (d.atom1, d.atom2, d.atom3, d.atom4)
        seq = tuple(a.type for a in atoms)
        if seq in wanted or tuple(reversed(seq)) in wanted:
            mask = tuple(a.idx + 1 for a in atoms)
            if mask not in found and tuple(reversed(mask)) not in found:
                found.append(mask)
    return found


def main(argv=None):
    args = parse_args(argv)
    terms = parse_terms(args)
    structure = pmd.load_file(args.prmtop)
    paths = matching_dihedrals(structure, args.types)

    if not paths:
        raise SystemExit(f"no dihedral of type {'-'.join(args.types)} in {args.prmtop}")

    if args.list:
        print(f"{len(paths)} dihedral(s) of type {'-'.join(args.types)}")
        for mask in paths:
            print("  " + " ".join(f"@{i}" for i in mask))
        return

    with open(args.output, "w") as fh:
        fh.write(f"# {len(paths)} dihedral(s) of type {'-'.join(args.types)} in "
                 f"{args.prmtop}\n")
        fh.write("# heights in kcal/mol, phases in degrees\n")
        for mask in paths:
            fh.write("deleteDihedral " + " ".join(f"@{i}" for i in mask) + "\n")
        for mask in paths:
            text = " ".join(f"@{i}" for i in mask)
            for height, n, delta in terms:
                fh.write(f"addDihedral {text} {height:.6f} {n} {delta:.4f} "
                         f"{args.scee} {args.scnb}\n")
        if args.parmout:
            fh.write(f"parmout {args.parmout}\n")

    print(f"{len(paths)} dihedral(s) x {len(terms)} term(s) -> {args.output}",
          file=sys.stderr)
    for height, n, delta in terms:
        print(f"  n={n}  K={height:.6f} kcal/mol  delta={delta:.4f} deg", file=sys.stderr)
    print(f"run it with:  parmed {args.prmtop} -i {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
