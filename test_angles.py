"""
Reconstructs atom coordinates from torsion angles
"""
import numpy as np
from Bio.PDB import *
from math import pi, isclose
from util import remove_hetero, main_load
from angles import rot_atom, rot_backbone, get_all_backbone_torsions

structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
# torsion_angles = polypeptide.get_phi_psi_list()
torsion_angles = get_all_backbone_torsions(polypeptide)

print([torsion[2] / pi * 180 for torsion in torsion_angles if torsion[2] is not None])
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

    no_rotation = np.allclose(new_coord.get_array(), res1_n.get_array())
    assert(no_rotation)

    torsion = 0
    new_coord = rot_atom(
        torsion, (res0['N'], res0['CA'], res0['C'], res1['N']))
    new_dist = res1['N'] - res0['C']

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
    torsion_angles[5] = (0.5, 1.4, .2)
    new_polypeptide = rot_backbone(torsion_angles, polypeptide)

    matches = []
    for index, res in enumerate(polypeptide):
        new_res = new_polypeptide[index]
        matches.append(np.allclose(res['N'].coord, new_res['N'].coord))
        matches.append(np.allclose(res['CA'].coord, new_res['CA'].coord))
        matches.append(np.allclose(res['C'].coord, new_res['C'].coord))

    assert(not all(matches))

    # Make more changes, and ensure that all changes match what we expect
    torsion_angles[10] = (1.123, 1.677, 1.4)
    torsion_angles[35] = (-2, 2, .05, 2)
    torsion_angles[55] = (3.14, .001, -2)
    torsion_angles[70] = (-1.01, -.001, 3)

    new_polypeptide = rot_backbone(torsion_angles, polypeptide)
    new_torsion_angles = get_all_backbone_torsions(new_polypeptide)

    matches = []
    for index, torsion_triple in enumerate(torsion_angles):
        if torsion_triple[0] is not None and new_torsion_angles[index][0] is not None:
            matches.append(
                isclose(torsion_triple[0], new_torsion_angles[index][0])
            )
        if torsion_triple[1] is not None and new_torsion_angles[index][1] is not None:
            matches.append(
                isclose(torsion_triple[1], new_torsion_angles[index][1])
            )
        if torsion_triple[2] is not None and new_torsion_angles[index][2] is not None:
            matches.append(
                isclose(torsion_triple[2], new_torsion_angles[index][2])
            )

    assert(all(matches))

# After messing around with Avogadro, it's easy to see that bond lengths
# should (in general) NOT be preserved when adjust torsion angles

# def test_rot_backbone_test_preserve():
#     """
#     Test to make sure that bond length and bond angle are preserved across
#     torsion angle modifications.
#     """
#     structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
#     torsion_angles = get_all_backbone_torsions(polypeptide)
    
#     # Introduce rounding error so isclose can work ()
#     polypeptide = rot_backbone(torsion_angles, polypeptide)
    
#     # Make additional changes - make sure that bond length is preserved
#     torsion_angles[9] = (-.68, -2.1, 1)
#     torsion_angles[10] = (1.123, 1.677, -.5)
#     torsion_angles[11] = (-2, 2, -2)
#     torsion_angles[12] = (3.14, .001, -2.5)
#     torsion_angles[13] = (-1.01, -.001, -3)

#     new_polypeptide = rot_backbone(torsion_angles, polypeptide)
#     new_torsion_angles = get_all_backbone_torsions(new_polypeptide)

#     matches = []
#     num_res = len(polypeptide)

#     for i in range(num_res):

#         res = polypeptide[i]
#         new_res = new_polypeptide[i]

#         # Get bond lengths
#         n_to_ca = res['N'] - res['CA']
#         new_n_to_ca = new_res['N'] - new_res['CA']
#         ca_to_c = res['CA'] - res['C']
#         new_ca_to_c = new_res['CA'] - new_res['C']

#         print(f'residue {i}')
#         print('N TO CA')
#         print(n_to_ca)
#         print(new_n_to_ca)
#         print('CA to C')
#         print(ca_to_c)
#         print(new_ca_to_c)
#         matches.append(isclose(n_to_ca, new_n_to_ca))
#         matches.append(isclose(ca_to_c, new_ca_to_c))

#         if i < num_res - 1:
#             res_next = polypeptide[i + 1]
#             new_res_next = new_polypeptide[i + 1]
#             c_to_n = res_next['N'] - res['C']
#             new_c_to_n = new_res_next['N'] - new_res['C']
#             print('C to N')
#             print(c_to_n)
#             print(new_c_to_n)
            

#             matches.append(isclose(c_to_n, new_c_to_n))

#     print(matches)
#     assert(all(matches))
