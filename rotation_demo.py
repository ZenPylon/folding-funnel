"""
Interactive exploration of torsion angle rotation.
"""
from util import main_load
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from math import pi
from matplotlib import animation

from angles import rot_atom

# TODO - make 3d plot with different angles
# Generate coordinates in script, import into folding_funnel notebook
num_res = 4
num_atoms = num_res * 3
num_frames = 200
interval = 50

# Tuple of atoms, array of positions
def update_anim(frame, atoms, positions, scatter, lines):
    offset = pi / 100
    new_coord = rot_atom(offset, atoms).get_array()

    atoms[3].coord = new_coord
    positions[:, 3] = new_coord
    updated_points = np.array([atom.coord.T for atom in atoms])

    lines[2].set_data(new_coord[0:2])
    lines[2].set_3d_properties(new_coord[2])
    scatter.set_data(positions[0:2, :])
    scatter.set_3d_properties(positions[2, :])

    print(atoms[3] - atoms[2])
    return atoms, positions, scatter, lines

def get_backbone(polypeptide):

    coords = np.zeros((3, num_atoms))
    for index, res in enumerate(polypeptide[0:num_res]):
        coords[:, index * 3] = res['N'].coord
        coords[:, index * 3 + 1] = res['CA'].coord
        coords[:, index * 3 + 2] = res['C'].coord

    return coords

# Return first four atom positions of polypeptide
def first_positions(polypeptide):
    res0 = polypeptide[0]
    res1 = polypeptide[1]
    coords = np.zeros((3, 4))
    coords[:, 0] = res0['N'].coord
    coords[:, 1] = res0['CA'].coord
    coords[:, 2] = res0['C'].coord
    coords[:, 3] = res1['N'].coord
    return coords


structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
first_positions = first_positions(polypeptide)
res0 = polypeptide[0]
res1 = polypeptide[1]

backbone = get_backbone(polypeptide)
fig = plt.figure()
ax = p3.Axes3D(fig)
ax.set_xlim(left=24, right=28)
ax.set_ylim(bottom=24, top=28)
ax.set_zlim(bottom=2, top=6)
# ax.scatter3D(backbone[0, :], backbone[1, :],
#              backbone[2, :], s=100, c=(0, 0, 0))


# lines = [ax.plot(backbone[0, :],
#                  backbone[1, :],
#                  backbone[2, :],
#                  c='b')[0] for index in range(num_atoms - 1)]
scatter = ax.plot(first_positions[0, :], first_positions[1, :],
             first_positions[2, :], linestyle='', marker='o', 
             markersize="10", c=(0, 0, 0))[0]
print(first_positions)
lines = [ax.plot(first_positions[0, :],
                 first_positions[1, :],
                 first_positions[2, :],
                 c='b')[0] for index in range(3)]

anim = animation.FuncAnimation(fig, update_anim, frames=num_frames,
                fargs=((res0['N'], res0['CA'], res0['C'], res1['N']), 
                first_positions, scatter, lines), interval=interval, blit=False)
# anim.save('amino_spin.gif', writer='imagemagick')
plt.show()
