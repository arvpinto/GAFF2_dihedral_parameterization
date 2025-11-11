#!/usr/bin/env python3
import sys
import math
import argparse
import numpy as np

HARTREE_TO_KCAL = 627.509474

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Read Gaussian .out files, extract a dihedral and final relative energy (kcal/mol)."
    )
    parser.add_argument(
        "atoms", nargs=4, type=int,
        help="Four atom indices (1-based) defining the dihedral."
    )
    parser.add_argument(
        "files", nargs="+",
        help="Gaussian output files (.out or .log)."
    )
    return parser.parse_args()

def get_final_energy(lines):
    """Return the final SCF energy in Hartree."""
    energy = None
    for line in lines:
        if "SCF Done:" in line:
            parts = line.split()
            try:
                energy = float(parts[4])
            except (IndexError, ValueError):
                continue
    return energy

def get_final_geometry(lines):
    """Extract the final Cartesian coordinates (in Angstroms)."""
    geom = []
    start = None
    for i, line in enumerate(lines):
        if "Standard orientation:" in line:
            start = i
    if start is None:
        return []
    # geometry starts 5 lines after "Standard orientation:"
    start += 5
    for line in lines[start:]:
        if "----" in line:
            break
        parts = line.split()
        if len(parts) >= 6:
            x, y, z = map(float, parts[3:6])
            geom.append((x, y, z))
    return geom

def calc_dihedral(p1, p2, p3, p4):
    """Return dihedral angle in degrees for four 3D points."""
    b1 = np.array(p2) - np.array(p1)
    b2 = np.array(p3) - np.array(p2)
    b3 = np.array(p4) - np.array(p3)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return math.degrees(math.atan2(y, x))

def main():
    args = parse_arguments()
    atom_indices = [a - 1 for a in args.atoms]  # convert to 0-based

    data = []  # list of (dihedral_deg, energy_kcal)

    # --- Read data from all files in given order ---
    for filename in args.files:
        with open(filename, 'r') as f:
            lines = f.readlines()
        energy = get_final_energy(lines)
        geom = get_final_geometry(lines)
        if energy is None or not geom:
            print(f"Warning: {filename} missing data", file=sys.stderr)
            continue
        dihedral = calc_dihedral(*(geom[i] for i in atom_indices))
        # Normalize dihedral angle to range [-180, 180]
        dihedral = ((dihedral + 180) % 360) - 180
        energy_kcal = energy * HARTREE_TO_KCAL
        data.append((dihedral, energy_kcal))

    if not data:
        print("No valid data found.")
        sys.exit(1)

    # --- Compute relative energies ---
    min_energy = min(e for _, e in data)
    rel_data = [(d, e - min_energy) for d, e in data]

    # --- Sort and print results by dihedral angle ---
    rel_data.sort(key=lambda x: x[0])  # sort by dihedral angle

    print(f"{'#Dihedral(deg)':>15} {'RelEnergy(kcal/mol)':>20}")
    for d, e in rel_data:
        print(f"{d:15.3f} {e:20.6f}")

if __name__ == "__main__":
    main()

