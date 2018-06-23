from util import main_load
import chemcoord as cc
import pandas as pd
from math import pi

structure, residue_list, polypeptide = main_load('1ubq', '1ubq.pdb')

rows = []
for res in polypeptide:
    for atom in res:
        name = atom.get_name() 
        if name == 'CA' or name == 'C' or name == 'N':
            rows.append([
                atom.get_name(), atom.coord[0], atom.coord[1], atom.coord[2]
            ])

torsion_angles = polypeptide.get_phi_psi_list()

for torsion in torsion_angles:
    angle0 = None
    angle1 = None
    if torsion[0] is not None:
        angle0 = torsion[0] * 180 / pi
    if torsion[1] is not None:
        angle1 = torsion[1] * 180 / pi
    print((angle0, angle1))

bonds = {}
num_atoms = len(rows)
bonds[0] = {1}
bonds[num_atoms - 1] = {num_atoms - 2}
for i in range(1, num_atoms - 1):
    bonds[i] = {i - 1, i + 1}
print(bonds)

df = pd.DataFrame(rows, columns=['atom', 'x', 'y', 'z'])
molecule = cc.Cartesian(df)
molecule.set_bonds(bonds)
zmat = molecule.get_zmat(use_lookup=True)
xyz = zmat.to_xyz()

with pd.option_context('display.max_rows', None):
    print(bonds)
    print(zmat)



