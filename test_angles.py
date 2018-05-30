"""
Reconstructs atom coordinates from torsion angles
"""
from Bio.PDB import *
from math import pi, isclose
from main import remove_hetero, main_load

structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
torsion_angles = polypeptide.get_phi_psi_list()

def test_coords_equal():
    """
    Ensure sure that after a phi or psi angle is updated, the relative 
    coordinates of all non-backbone atoms stays the same relative to the backbone
    """
    res = residue_list[0]
    nitro_pos = res['N'].get_vector()
    carbon_pos = res['C'].get_vector()
    carbon_a_pos = res['CA'].get_vector()
    
    # Get all atoms in the side-chain (i.e. those not in the backbone)
    backbone = [res['N'], res['C'], res['CA'], res['O']]
    sidechain = [atom for atom in res.get_atoms() if atom not in backbone]
    print(sidechain)

    # center at origin
    nitro_pos = nitro_pos - carbon_a_pos
    carbon_pos = carbon_pos - carbon_a_pos

    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis(-pi * 120.0/180.0, carbon_pos)

    # apply rotation to ca-n vector
    cb_at_origin = nitro_pos.left_multiply(rot)

    # isclose(a, b)
