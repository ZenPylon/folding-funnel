from util import main_load
import chemcoord as cc
import pandas as pd

structure, residue_list, polypeptide = main_load('1ubq', '1ubq.pdb')

rows = []
for res in polypeptide:
    for atom in res:
        rows.append([
            atom.get_name(), atom.coord[0], atom.coord[1], atom.coord[2]
        ])

df = pd.DataFrame(rows, columns=['atom', 'x', 'y', 'z'])
molecule = cc.Cartesian(df)
# molecule.set_bonds({})
bonds = molecule.get_bonds(use_lookup=True)
zmat = molecule.get_zmat()
xyz = zmat.to_xyz()
with pd.option_context('display.max_rows', None):
    print(xyz.columns)
    print(molecule.columns)
    print(cc.xyz_functions.isclose(xyz, molecule))



