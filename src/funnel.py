"""
Testing the funnel hypothesis
"""
import util
import chemcoord as cc
from amino_data import get_amino_bonds, get_amino_atom_count
structure, residue_list, polypeptide = util.main_load('ubiq', '1ubq.pdb')

atom_counter = 0
for res in polypeptide:
    name = res.get_name()
    get_amino_bonds(name, atom_counter)
    atom_counter += get_amino_atom_count(name)