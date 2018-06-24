"""
Interactive exploration of torsion angle rotation.
"""
from copy import deepcopy
from util import main_load
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from math import pi
from matplotlib import animation
import random

from angles import rot_atom, rot_backbone, get_all_backbone_torsions

# TODO - make 3d plot with different angles
# Generate coordinates in script, import into folding_funnel notebook
num_frames = 500
interval = 100
structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
starting_angles = get_all_backbone_torsions(polypeptide)
torsion_angles = deepcopy(starting_angles)
backbone = []
num_atoms = 3 * len(polypeptide)

# Tuple of atoms, array of positions


def update_anim(frame, positions, scatter, lines):
    if frame < 100:
        torsion_index = 12
    elif frame < 200:
        torsion_index = 24
    elif frame < 300:
        torsion_index = 36
    elif frame < 400:
        torsion_index = 48
    else:
        torsion_index = 60

    starting_angle = starting_angles[torsion_index][0]
    print(f'frame {frame}')
    current_angle = starting_angle + (frame / 100) * 2*pi
    print(f'current_angle {current_angle}')
    torsion_angles[torsion_index] = (
        current_angle, torsion_angles[torsion_index][1], torsion_angles[torsion_index][2])

    # Update points
    new_polypeptide = rot_backbone(torsion_angles, polypeptide)
    for index in range(len(new_polypeptide)):
        res = new_polypeptide[index]
        polypeptide[index]['N'].coord = res['N'].coord
        polypeptide[index]['CA'].coord = res['CA'].coord
        polypeptide[index]['C'].coord = res['C'].coord

    positions = get_backbone(new_polypeptide)

    for line in lines:
        line.set_data(positions[0:2, :])
        line.set_3d_properties(positions[2, :])

    scatter.set_data(positions[0:2, :])
    scatter.set_3d_properties(positions[2, :])

    return positions, scatter, lines


def get_backbone(polypeptide):
    coords = np.zeros((3, num_atoms))
    for index, res in enumerate(polypeptide):
        coords[:, index * 3] = res['N'].coord
        coords[:, index * 3 + 1] = res['CA'].coord
        coords[:, index * 3 + 2] = res['C'].coord

    return coords


backbone = get_backbone(polypeptide)
fig = plt.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=40, azim=None)
ax.set_xlim(left=np.min(backbone[0, :]), right=np.max(backbone[0, :]))
ax.set_ylim(bottom=np.min(backbone[1, :]), top=np.max(backbone[1, :]))
ax.set_zlim(bottom=np.min(backbone[2, :]), top=np.max(backbone[2, :]))

scatter = ax.plot(backbone[0, :], backbone[1, :],
                  backbone[2, :], linestyle='', marker='o',
                  markersize="4", c='black')[0]

lines = [ax.plot(backbone[0, :],
                 backbone[1, :],
                 backbone[2, :],
                 c='b')[0] for i in range(num_atoms)]

# anim = animation.FuncAnimation(fig, update_anim, frames=num_frames,
                # fargs=(backbone, scatter, lines), interval=interval,
                # repeat_delay=1000, blit=False)
# anim.save('backbone.gif', writer='imagemagick')
# plt.show()


# Second one with modified torsion angles
# for plot in range(3):
#     for index, torsion in enumerate(torsion_angles):
#         offset = .3 
#         torsion0 = torsion[0] + random.uniform(-offset, offset) if torsion[0] is not None else None
#         torsion1 = torsion[1] + random.uniform(-offset, offset) if torsion[1] is not None else None
#         torsion2 = torsion[2] + random.uniform(-offset, offset) if torsion[2] is not None else None

#         torsion_angles[index] = (torsion0, torsion1, torsion2)

#     new_polypeptide = rot_backbone(torsion_angles, polypeptide)
#     backbone = get_backbone(new_polypeptide)
#     fig = plt.figure()
#     ax = p3.Axes3D(fig)
#     ax.view_init(elev=40, azim=None)
#     ax.set_xlim(left=np.min(backbone[0, :]), right=np.max(backbone[0, :]))
#     ax.set_ylim(bottom=np.min(backbone[1, :]), top=np.max(backbone[1, :]))
#     ax.set_zlim(bottom=np.min(backbone[2, :]), top=np.max(backbone[2, :]))

#     scatter = ax.plot(backbone[0, :], backbone[1, :],
#                     backbone[2, :], linestyle='', marker='o',
#                     markersize="4", c='black')[0]

#     lines = [ax.plot(backbone[0, :],
#                     backbone[1, :],
#                     backbone[2, :],
#                     c='b')[0] for i in range(num_atoms)]

#         # anim = animation.FuncAnimation(fig, update_anim, frames=num_frames,
#                         # fargs=(backbone, scatter, lines), interval=interval,
#                         # repeat_delay=1000, blit=False)
#         # anim.save('backbone.gif', writer='imagemagick')
# plt.show()
