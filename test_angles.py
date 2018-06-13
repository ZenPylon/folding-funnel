"""
Reconstructs atom coordinates from torsion angles
"""
import numpy as np
from Bio.PDB import *
from math import pi, isclose
from main import remove_hetero, main_load
from angles import rot_atom, rot_backbone

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


def test_rot_atom():
    res0 = polypeptide[0]
    res1 = polypeptide[1]
    res1_n = res1['N'].get_vector()
    new_coord = rot_atom(0, (res0['N'], res0['CA'], res0['C'], res1['N']))

    print(f'rotated vector {new_coord.get_array()}')
    print(f'original vector {res1_n.get_array()}')
    no_rotation = np.allclose(new_coord.get_array(), res1_n.get_array())
    assert(no_rotation)

    # TODO - do gradual rotation on coordinates to explore space

    # res0['N'].coord = Vector(0, 0, 0)
    # res0['CA'].coord = Vector(3, 0, 0)
    # res0['C'].coord = Vector(5, 0, 0)
    # res1['N'].coord = Vector(5, 1, 0)
    # res1_n = res1['N'].get_vector()

    # offset_angle = 3 * pi / 2
    # expected_rotation = Vector(5, -1, 0) 
    # new_coord = rot_atom(offset_angle, (res0['N'], res0['CA'], res0['C'], res1['N']))

    # print(f'rotated vector {new_coord.get_array()}')
    # print(f'expected vector {expected_rotation.get_array()}')
    # rotation_match = np.allclose(new_coord.get_array(), expected_rotation.get_array())
    # assert(rotation_match)

