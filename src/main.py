from math import pi
from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
import chemcoord as cc
import numpy as np
import pandas as pd
import util

np.random.seed(20)

# Don't worry about hydrogen atoms - we'll add them to the PDB model later
ubiq_zmat = '1ubq.xyz'
ubiq_pdb = '1ubq.pdb'
# structure, residue_list, polypeptide = util.main_load('ubiq', ubiq_pdb)

# for torsions in polypeptide.get_phi_psi_list():
#     angle0 = None
#     angle1 = None
#     if torsions[0] is not None:
#         angle0 = torsions[0] / pi * 180
#     if torsions[1] is not None:
#         angle1 = torsions[1] / pi * 180

#     print((angle0, angle1))


def get_modeller(pdb_file: str):
    pdb = PDBFile('1ubq.pdb')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    return modeller


def get_zmat(zmat_filename: str):
    """
    Loads a zmat file from firebase.  Returns None if it doesn't exist

    Args:
        zmat_filename - the file name in firebase
    """
    # zmat_file = load_from_firebase(zmat_filename)
    # zmat = cc.ZMat.read_zmat(zmat_file)
    # return zmat
    return None


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


def save_starting_zmat(zmat_name: str, zmat: cc.Zmat):
    """
    Saves the initial zmat data (unmodified internal coordinates 
    of the PDB file) to firebase

    Args:
        zmat_name - name of file to write
        zmat - the zmat data to save
    """
    print('TODO - save zmat')
    # save_to_firebase(zmat)


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

def init_model():
    modeller = get_modeller(ubiq_pdb)
    zmat = get_zmat(ubiq_zmat)
    if zmat is None:
        zmat = calculate_zmat(modeller)
        save_starting_zmat(ubiq_zmat, zmat)

    torsion_indices = get_torsion_indices(zmat)
    print(torsion_indices)
    starting_torsions = np.array([zmat.loc[torsion_indices[:, 0], 'dihedral'],
                                  zmat.loc[torsion_indices[:, 1], 'dihedral']]).T

    
    offsets = np.random.choice([0, 0, -1, 1], starting_torsions.shape)
    run_sequence(zmat, starting_torsions, torsion_indices, offsets)


def run_sequence(zmat: cc.Zmat, starting_torsions, torsion_indices, offsets, offset_size=4, num_configs=3):
    """
    Args:
        pdb_name: name field in firebase
        pdb_file: path to pdb file
        offset_size: amount to deviate on each configuration
        num_configs: number of configurations to try
    """

    print(f'\n\nstarting torsions {starting_torsions.shape}\n')
    print(starting_torsions)
    print(f'\n\noffsets {offsets.shape}\n')
    print(offsets)
    # Modify (on average) half the torsion angles
    for i in range(num_configs):
        zmat.safe_loc[torsion_indices[:, 0], 'dihedral'] = \
            starting_torsions[:, 0] + (offsets[:, 0] * i * offset_size)
        zmat.safe_loc[torsion_indices[:, 1], 'dihedral'] = \
            starting_torsions[:, 1] + (offsets[:, 1] * offset_size * i)
        print(zmat.loc[torsion_indices[:, 0], 'dihedral'])
        print(zmat.loc[torsion_indices[:, 1], 'dihedral'])


init_model()

