"""
Testing the funnel hypothesis
"""
from math import pi
from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
from simtk.openmm import LangevinIntegrator
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
import chemcoord as cc
import numpy as np
import pandas as pd
import util

pdb = PDBFile('1ubq.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
pdb_bonds = modeller.topology.bonds()
atoms = modeller.topology.atoms()
positions = modeller.getPositions()

cc_bonds = {}
cc_positions = np.zeros((3, modeller.topology.getNumAtoms()))
atom_names = []

# Construct bond dictionary and positions chemcoord 
for index, atom in enumerate(atoms):
    cc_bonds[index] = set()
    pos = positions[index] / nanometer
    atom_names.append(atom.name)
    cc_positions[:, index] = pos

nitro_list = []

for bond in pdb_bonds:
    if bond[0].name == 'N' or bond[1].name == 'N':
        nitro_list.append(bond)
    
    cc_bonds[bond[0].index].add(bond[1].index)
    cc_bonds[bond[1].index].add(bond[0].index)

print(len(nitro_list))
for nitro in nitro_list:
    print(nitro)

with pd.option_context('display.max_rows', None):
    cc_df = pd.DataFrame({
            'atom': atom_names, 
            'x': cc_positions[0, :], 
            'y': cc_positions[1, :], 
            'z': cc_positions[2, :]
    })

    molecule = cc.Cartesian(cc_df)
    molecule.set_bonds(cc_bonds)
    molecule._give_val_sorted_bond_dict(use_lookup=True)
    zmat = molecule.get_zmat(use_lookup=True)
    zmat.to_zmat(buf='zmat.xyz', implicit_index=False)
    molecule2 = zmat.get_cartesian()

    phi_indices = []
    psi_indices = []

    print('getting angles')
    print(zmat.shape)
    for i in range(len(zmat.index)):
        b_index = zmat.loc[i, 'b'] 
        a_index = zmat.loc[i, 'a'] 
        d_index = zmat.loc[i, 'd']

        if isinstance(b_index, str) or isinstance(a_index, str) or isinstance(d_index, str):
            print('skipping')
            print(zmat.loc[i])
            continue

        # Psi angles
        if (zmat.loc[i, 'atom'] == 'N') & \
                (zmat.loc[b_index, 'atom'] == 'CA') & \
                (zmat.loc[a_index, 'atom'] == 'C') & \
                (zmat.loc[d_index, 'atom'] == 'N'):
            psi_indices.append(i)

        elif (zmat.loc[i, 'atom'] == 'N') & \
                (zmat.loc[b_index, 'atom'] == 'C') & \
                (zmat.loc[a_index, 'atom'] == 'CA') & \
                (zmat.loc[d_index, 'atom'] == 'N'):
            psi_indices.append(i)
        
        elif (zmat.loc[i, 'atom'] == 'C') & \
                (zmat.loc[b_index, 'atom'] == 'N') & \
                (zmat.loc[a_index, 'atom'] == 'CA') & \
                (zmat.loc[d_index, 'atom'] == 'C'):
            phi_indices.append(i)

        elif (zmat.loc[i, 'atom'] == 'C') & \
                (zmat.loc[b_index, 'atom'] == 'CA') & \
                (zmat.loc[a_index, 'atom'] == 'N') & \
                (zmat.loc[d_index, 'atom'] == 'C'):
            phi_indices.append(i)

    for psi in psi_indices:
        print(psi)
    
    print(f'phi indices {len(phi_indices)}')
    print(phi_indices)
    print(zmat)
    print(zmat.safe_loc[psi_indices, 'dihedral'])
    zmat.safe_loc[psi_indices, 'dihedral'] = np.zeros(len(psi_indices))
    print(zmat.loc[psi_indices])

    structure, residue_list, polypeptide = util.main_load('ubiq', '1ubq.pdb')
    
    psi_list = []
    for torsions in polypeptide.get_phi_psi_list():
        if torsions[1] is not None:
            psi_list.append(torsions[1] / pi * 180)
    
    print(psi_list)
    # print(zmat.loc[psi_indices, 'dihedral'] - np.array(psi_list))
        



