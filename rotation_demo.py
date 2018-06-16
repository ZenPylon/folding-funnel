"""
Interactive exploration of torsion angle rotation.
"""
from util import main_load
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from matplotlib import animation

from angles import rot_atom

# TODO - make 3d plot with different angles
# Generate coordinates in script, import into folding_funnel notebook
num_res = 4
num_atoms = num_res * 3


def get_backbone(polypeptide):

    coords = np.zeros((3, num_atoms))
    for index, res in enumerate(polypeptide[0:num_res]):
        coords[:, index * 3] = res['N'].coord
        coords[:, index * 3 + 1] = res['CA'].coord
        coords[:, index * 3 + 2] = res['C'].coord

    return coords

# Return first four atoms of polypeptide
def first_four(polypeptide):
    res0 = polypeptide[0]
    res1 = polypeptide[1]
    coords = np.zeros((3, 4))
    coords[:, 0] = res0['N'].coord
    coords[:, 1] = res0['CA'].coord
    coords[:, 2] = res0['C'].coord
    coords[:, 3] = res1['N'].coord
    return coords


structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
first_four = first_four(polypeptide)
res0 = polypeptide[0]
res1 = polypeptide[1]

backbone = get_backbone(polypeptide)
fig = plt.figure()
ax = p3.Axes3D(fig)
# ax.scatter3D(backbone[0, :], backbone[1, :],
#              backbone[2, :], s=100, c=(0, 0, 0))


# lines = [ax.plot(backbone[0, :],
#                  backbone[1, :],
#                  backbone[2, :],
#                  c='b')[0] for index in range(num_atoms - 1)]
ax.scatter3D(first_four[0, :], first_four[1, :],
             first_four[2, :], s=100, c=(0, 0, 0))

lines = [ax.plot(first_four[0, :],
                 first_four[1, :],
                 first_four[2, :],
                 c='b')[0] for index in range(3)]
new_pos = rot_atom(1.0, (res0['N'], res0['CA'], res0['C'], res1['N']))
plt.show()
# animation.FuncAnimation(fig, )
