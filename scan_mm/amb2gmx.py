#!/usr/bin/env python3

import sys
import parmed as pmd

def main():
    if len(sys.argv) != 3:
        print("Usage: python amber2gmx.py <topology.prmtop> <coordinates.inpcrd>")
        sys.exit(1)

    prmtop_file = sys.argv[1]
    inpcrd_file = sys.argv[2]

    try:
        print(f"Loading AMBER files: {prmtop_file}, {inpcrd_file}")
        amber_structure = pmd.load_file(prmtop_file, xyz=inpcrd_file)

        top_file = prmtop_file.rsplit('.', 1)[0] + '_converted.top'
        gro_file = inpcrd_file.rsplit('.', 1)[0] + '_converted.gro'

        print(f"Saving GROMACS files: {top_file}, {gro_file}")
        amber_structure.save(top_file, format='gromacs', combine='all')
        amber_structure.save(gro_file, format='gro')

        print("Conversion complete.")

    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

