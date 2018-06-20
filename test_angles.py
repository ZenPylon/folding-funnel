"""
Reconstructs atom coordinates from torsion angles
"""
import numpy as np
from Bio.PDB import *
from math import pi, isclose
from util import remove_hetero, main_load
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
    first_dist = res1['N'] - res0['C']
    res1_n = res1['N'].get_vector()
    torsion = calc_dihedral(
        res0['N'].get_vector(),
        res0['CA'].get_vector(),
        res0['C'].get_vector(),
        res1['N'].get_vector(),
    )
    new_coord = rot_atom(
        torsion, (res0['N'], res0['CA'], res0['C'], res1['N']))

    # print(f'rotated vector {new_coord.get_array()}')
    # print(f'original vector {res1_n.get_array()}')
    no_rotation = np.allclose(new_coord.get_array(), res1_n.get_array())
    assert(no_rotation)

    torsion = 0
    new_coord = rot_atom(
        torsion, (res0['N'], res0['CA'], res0['C'], res1['N']))
    new_dist =  res1['N'] - res0['C']

    print(new_dist)
    print(first_dist)
    assert(isclose(new_dist, first_dist))
    


def test_rot_backbone():
    new_polypeptide = rot_backbone(torsion_angles, polypeptide)
    matches = []

    for index, res in enumerate(polypeptide):
        new_res = new_polypeptide[index]
        matches.append(np.allclose(res['N'].coord, new_res['N'].coord))
        matches.append(np.allclose(res['CA'].coord, new_res['CA'].coord))
        matches.append(np.allclose(res['C'].coord, new_res['C'].coord))

    assert(all(matches))

    # Modify the structure, make sure the change is registered
    torsion_angles[5] = (0.5, 1.4)
    new_polypeptide = rot_backbone(torsion_angles, polypeptide)

    matches = []
    for index, res in enumerate(polypeptide):
        new_res = new_polypeptide[index]
        matches.append(np.allclose(res['N'].coord, new_res['N'].coord))
        matches.append(np.allclose(res['CA'].coord, new_res['CA'].coord))
        matches.append(np.allclose(res['C'].coord, new_res['C'].coord))

    assert(not all(matches))

    # Make more changes, and ensure that all changes match what we expect
    torsion_angles[10] = (1.123, 1.677)
    torsion_angles[35] = (-2, 2)
    torsion_angles[55] = (3.14, .001)
    torsion_angles[70] = (-1.01, -.001)

    new_polypeptide = rot_backbone(torsion_angles, polypeptide)
    new_torsion_angles = new_polypeptide.get_phi_psi_list()

    matches = []
    for index, torsion_pair in enumerate(torsion_angles):
        if torsion_pair[0] is not None and new_torsion_angles[index][0] is not None:
            matches.append(
                isclose(torsion_pair[0], new_torsion_angles[index][0])
            )
        if torsion_pair[1] is not None and new_torsion_angles[index][1] is not None:
            matches.append(
                isclose(torsion_pair[1], new_torsion_angles[index][1])
            )

    assert(all(matches))


def test_rot_backbone_test_preserve():
    """
    Test to make sure that bond length and bond angle are preserved across
    torsion angle modifications.
    """
    structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
    torsion_angles = polypeptide.get_phi_psi_list()

    # Make additional changes - make sure that bond length is preserved
    torsion_angles[5] = (-.68, -2.1)
    torsion_angles[13] = (1.123, 1.677)
    torsion_angles[38] = (-2, 2)
    torsion_angles[50] = (3.14, .001)
    torsion_angles[65] = (-1.01, -.001)

    new_polypeptide = rot_backbone(torsion_angles, polypeptide)
    new_torsion_angles = new_polypeptide.get_phi_psi_list()

    matches = []
    num_res = len(polypeptide)

    for i in range(num_res):
        res = polypeptide[i]
        new_res = new_polypeptide[i]

        # Get bond lengths
        n_to_ca = res['N'] - res['CA']
        new_n_to_ca = new_res['N'] - new_res['CA']
        ca_to_c = res['CA'] - res['C']
        new_ca_to_c = new_res['CA'] - new_res['C']

        # print(n_to_ca)
        # print(new_n_to_ca)
        matches.append(isclose(n_to_ca, new_n_to_ca))
        matches.append(isclose(ca_to_c, new_ca_to_c))

        if i < num_res - 1:
            res_next = polypeptide[i + 1]
            new_res_next = new_polypeptide[i + 1]
            c_to_n = res_next['N'] - res['C']
            new_c_to_n = new_res_next['N'] - new_res['C']

            matches.append(isclose(c_to_n, new_c_to_n))

    assert(all(matches))
