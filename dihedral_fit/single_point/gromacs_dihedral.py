#!/usr/bin/env python3
import sys
import re

# GROMACS energies are typically in kJ/mol — conversion factor to kcal/mol
KJ_TO_KCAL = 1 / 4.184

def extract_energy(filename):
    """Extract potential energy (in kcal/mol) from a GROMACS log file."""
    with open(filename, 'r') as f:
        for line in f:
            if "Potential Energy" in line:
                parts = line.split()
                for val in parts:
                    try:
                        energy = float(val)
                        return energy * KJ_TO_KCAL
                    except ValueError:
                        continue
    raise ValueError(f"Could not find potential energy in {filename}")

def extract_angle(filename):
    """Extract floating-point dihedral angle from filename (e.g. dihe45.0.log -> 45.0)."""
    match = re.search(r"([-+]?\d*\.?\d+)", filename)
    if not match:
        raise ValueError(f"No floating-point angle found in filename: {filename}")
    return float(match.group(1))

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py dihe*log")
        sys.exit(1)

    results = []
    for fname in sys.argv[1:]:
        angle = extract_angle(fname)
        energy = extract_energy(fname)
        results.append((angle, energy))

    # Find minimum energy and compute relative energies
    min_energy = min(e for _, e in results)
    results = [(ang, e - min_energy) for ang, e in results]

    # Sort by dihedral angle before printing
    results.sort(key=lambda x: x[0])

    print(f"{'#Dihedral(deg)':>15} {'RelEnergy(kcal/mol)':>20}")
    for ang, e in results:
        print(f"{ang:.2f}\t{e:.6f}")

if __name__ == "__main__":
    main()

