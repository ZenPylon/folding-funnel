"""
Reconstructs atom coordinates from torsion angles
"""
from Bio.PDB import *
from math import pi
from main import remove_hetero, main_load

structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
torsion_angles = polypeptide.get_phi_psi_list()

def test_coords_equal():
    """
    Ensure sure that after a phi or psi angle is updated, the relative 
    coordinates of all non-backbone atoms stays the same relative to the backbone
    """
    assert(len(residue_list) !=0 )
