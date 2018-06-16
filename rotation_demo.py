"""
Interactive exploration of torsion angle rotation.
"""
from util import main_load
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from matplotlib import animation

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


structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
backbone = get_backbone(polypeptide)
fig = plt.figure()
ax = p3.Axes3D(fig)
ax.scatter3D(backbone[0, :], backbone[1, :], backbone[2, :], s=100, c=(0, 0, 0))
lines = [ax.plot(backbone[0, :],
                 backbone[1, :],
                 backbone[2, :],
                 c='b')[0] for index in range(num_atoms - 1)]

plt.show()
# animation.FuncAnimation(fig, )
