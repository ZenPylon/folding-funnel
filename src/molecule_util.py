from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
import numpy as np
import pandas as pd
import chemcoord as cc

def get_modeller(pdb_file: str):
    pdb = PDBFile('data/1ubq.pdb')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    return modeller


def calculate_zmat(modeller: Modeller):
    """
    Calculates the zmat from an OpenMM modeller
    """
    # Create new document with field pdb_name, doc_id as random string
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

    for bond in pdb_bonds:
        cc_bonds[bond[0].index].add(bond[1].index)
        cc_bonds[bond[1].index].add(bond[0].index)

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
    return zmat


def get_torsion_indices(zmat):
    """
    Calculates indices into the zmatrix which correspond to phi
    and psi angles.

    Args:
        zmat: the zmatrix specifying the molecule
    Returns:
        a numpy.array, with first column as phi_indices, second column 
        as psi_indices
    """
    phi_indices = []
    psi_indices = []

    for i in range(len(zmat.index)):
        b_index = zmat.loc[i, 'b'] 
        a_index = zmat.loc[i, 'a'] 
        d_index = zmat.loc[i, 'd']

        # If this molecule references a magic string (origin, e_x, e_y, e_z, etc)
        if isinstance(b_index, str) or isinstance(a_index, str) or isinstance(d_index, str):
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
    
    return np.array([phi_indices, psi_indices]).T
