#!/usr/bin/env python3
"""Convert the Dissolve Cu-BTC (HKUST-1) example into a PDB file for EPSR's readpdb.

usage: make_pdb.py DISSOLVE_INPUT OUTPUT_PDB

DISSOLVE_INPUT is experimental_data/cu-btc/input.txt. The species there is one
cubic unit cell of the framework (624 atoms, a = 26.3336 Å). EPSR cannot treat a
periodic network as one molecule, so the framework is split at its Cu-O bonds:

* each benzene-1,3,5-tricarboxylate linker (C9H3O6, 18 atoms, held together by the
  Dissolve C1-O, C1-C2, C2-C3 and C3-H bonds) becomes one residue -> one molecule;
* each Cu2 paddlewheel (two Cu 2.629 Å apart; the next Cu is 8 Å away) becomes one
  residue -> a two-atom molecule. Single Cu residues cannot be used: readpdb (in
  the Fortran and, faithfully, in the Rust port) writes single-atom residues
  without their atom records (see test_problems/.../EPSR-test-ex49).

The Cu-O coordination bonds are then intermolecular and are described by the
Lennard-Jones and Coulomb terms of the reference potential (and the empirical
potential). The Cu pairs are written first, then the linkers, each with its atoms
in the order O, C1, C2, C3, H, so that readpdb numbers the atom types
CU=1, O=2, C1=3, C2=4, C3=5, H=6. The PDB atom name is the Dissolve atom type,
which readpdb uses as the EPSR atom label (its first two columns give the element).
No CONECT records are needed: readpdb finds the bonds from van der Waals radii.
"""

import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

TYPE_ORDER = ["O", "C1", "C2", "C3", "H"]  # order of the atoms within each linker


class Atom(NamedTuple):
    element: str
    x: float
    y: float
    z: float
    atom_type: str  # the Dissolve atom type, e.g. "C1"


def read_dissolve_species(path):
    """Atoms (by Dissolve index), bonds and box length of the species in the input file."""
    atoms, bonds, box_length = {}, [], None
    for line in Path(path).read_text().splitlines():
        words = line.split()
        if len(words) >= 7 and words[0] == "Atom" and words[1].isdigit():
            atoms[int(words[1])] = Atom(
                element=words[2],
                x=float(words[3]),
                y=float(words[4]),
                z=float(words[5]),
                atom_type=words[6].strip("'"),
            )
        elif len(words) >= 3 and words[0] == "Bond" and words[1].isdigit():
            bonds.append((int(words[1]), int(words[2])))
        elif words[:1] == ["BoxLengths"]:
            box_length = float(words[1])
    return atoms, bonds, box_length


def linker_groups(atoms, bonds):
    """Connected groups of atoms when the Cu-O bonds are removed."""
    parent = {index: index for index in atoms}

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in bonds:
        if "Cu" in (atoms[first].atom_type, atoms[second].atom_type):
            continue
        parent[root(first)] = root(second)
    groups = defaultdict(list)
    for index in sorted(atoms):
        groups[root(index)].append(index)
    return [group for group in groups.values() if len(group) > 1]


def paddlewheel_pairs(atoms, box_length):
    """Pairs of Cu atoms that are each other's nearest Cu neighbour (minimum image)."""

    def distance(first, second):
        components = [
            atoms[first].x - atoms[second].x,
            atoms[first].y - atoms[second].y,
            atoms[first].z - atoms[second].z,
        ]
        components = [
            value - box_length * round(value / box_length) for value in components
        ]
        return sum(value * value for value in components) ** 0.5

    copper_atoms = [index for index in sorted(atoms) if atoms[index].atom_type == "Cu"]
    nearest = {}
    for index in copper_atoms:
        other_copper_atoms = [other for other in copper_atoms if other != index]
        nearest[index] = min(
            other_copper_atoms, key=lambda other: distance(index, other)
        )
    pairs = []
    for index in copper_atoms:
        partner = nearest[index]
        if nearest[partner] != index:
            raise SystemExit(f"Cu {index} has no mutual nearest Cu neighbour")
        if index < partner:
            pairs.append((index, partner))
    return pairs


def pdb_atom_line(serial, atom, residue_name, residue_number):
    """One ATOM record (readpdb reads only ATOM records, not HETATM).

    readpdb takes the element from columns 13-14 and the label from columns 13-16,
    so names of one-letter elements start in column 14 (" C1 ") and names of
    two-letter elements in column 13 ("CU  "), upper case as in standard PDB files.
    Copper therefore gets the EPSR label CU.
    """
    # readpdb compares columns 13-14 with upper-case element symbols.
    name = atom.atom_type.upper()
    name_field = f" {name:<3s}" if len(atom.element) == 1 else f"{name:<4s}"
    return (
        f"ATOM  {serial:5d} {name_field} {residue_name:3s} A{residue_number:4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {atom.element:>2s}"
    )


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    atoms, bonds, box_length = read_dissolve_species(sys.argv[1])
    copper_pairs = paddlewheel_pairs(atoms, box_length)
    linkers = linker_groups(atoms, bonds)
    if 2 * len(copper_pairs) + sum(len(linker) for linker in linkers) != len(atoms):
        raise SystemExit("the framework does not split into Cu pairs and linkers")

    lines = [
        "REMARK   Cu-BTC (HKUST-1) unit cell from the Dissolve example "
        "experimental_data/cu-btc/input.txt",
        f"CRYST1{box_length:9.4f}{box_length:9.4f}{box_length:9.4f}"
        f"  90.00  90.00  90.00 P 1           1",
    ]
    serial, residue_number = 0, 0
    for pair in copper_pairs:
        residue_number += 1
        for index in pair:
            serial += 1
            lines.append(pdb_atom_line(serial, atoms[index], "CU2", residue_number))
    for linker in linkers:
        residue_number += 1
        ordered = sorted(
            linker, key=lambda index: (TYPE_ORDER.index(atoms[index].atom_type), index)
        )
        for index in ordered:
            serial += 1
            lines.append(pdb_atom_line(serial, atoms[index], "BTC", residue_number))
    lines.append("END")
    Path(sys.argv[2]).write_text("\n".join(lines) + "\n")
    print(
        f"wrote {serial} atoms in {residue_number} residues "
        f"({len(copper_pairs)} Cu2 paddlewheels, {len(linkers)} linkers) to {sys.argv[2]}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
