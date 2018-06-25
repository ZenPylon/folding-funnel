from math import pi
import numpy as np
import util

np.random.seed(20)

# Don't worry about hydrogen atoms - we'll add them to the PDB model later
structure, residue_list, polypeptide = util.main_load('ubiq', '1ubq.pdb')
for torsions in polypeptide.get_phi_psi_list():
    angle0 = None
    angle1 = None
    if torsions[0] is not None:
        angle0 = torsions[0] / pi * 180
    if torsions[1] is not None:
        angle1 = torsions[1] / pi * 180

    print((angle0, angle1))



def save_starting_zmat(pdb_name: str, zmat):
    """
    Saves the initial zmat data (unmodified internal coordinates 
    of the PDB file) to firebase
    """
    # Create new document with field pdb_name, doc_id as random string
    print('TODO')

def load_starting_zmat(pdb_name: str):
    """
    Loads the initial zmat data (unmodified internal coordinates 
    of the PDB file) from firebase
    """
    print('TODO')

def get_torsion_indices(zmat):
    """
    Return a numpy.array, with first column as phi_indices, second column as psi_indices

    Args:
        zmat: the zmatrix specifying the molecule
    """

def begin_seqeunce(pdb_name: str, pdb_file: str, offset_size=4 num_configs=30):
    """
    Args:
        pdb_name: name field in firebase
        pdb_file: path to pdb file
        num_configs: number of configurations to try
    """
    modeller = get_modeller(pdb_file)
    zmat = get_zmat(modeller)
    save_starting_zmat(pdb_name, zmat)
    torsions = get_torsion_indices(zmat)
    offsets = np.random.choice([0, 0, -1, 1], torsions.shape)

    print('\n\ntorsions\n\n')
    print(torsions)
    print('offsets')
    print(offsets)

    # Modify (on average) half the torsion angles
    for i in range(num_configs):
        zmat.safe_loc[torsions[:, 0], 'dihedral'] += offsets[:, 0] + (i * offset_size)







# TODO - test by changing angles more and more, and plot vs. RMSD distance
#        Expectation is that it RMSD should increase as angle distance increases

# Construct backbone from angles and distances

# 1. Get torsion angles and bond distances for backbone molecules
# 1. Get vectors for sidechain atoms relative to backbone atom
# 1. Iterate through each residue and construct each backbone atom 
#    relative to transformed bond angle (use spherical coordinates)


# THEN
# 1. Save PDB with minor alterations
# 1. Calculate RMSD.  Should be small
# 1. Plot superimposed or side-by-side structures.  Get visual confirmation
# of similarity

# THEN:
# 1. Make increasingly large deviations from loaded structure
# 1. Calculate RMSD
# 1. Plot Angle distance vs RMSD.  Should be positive correlation


# TODO: The Algorithm:
# 
# 1. Specify all phi / psi angle deltas (magnitude M) relative to native structure
#    and get absolute phi / psi coordinates (calculate native angles + deltas).
# 2. Convert from angle-space to cartesian coordinates.
# 3. Pass coordinates to OpenMM (how?).
    # Use topology.addChain, topology.addResidue(), etc based on constructed structure
# 4. Add solvent to model with OpenMM.
# 5. Run energy minimization on model.
# 6. Pass final energy value and coordinates to BioPython, and reset atom coords.
# 7. Recalculate phi / psi angles from coordinates.
# 8. Calculate distance in angle-space from starting point (in this case,
#    native structure) and calculate delta in energy level.
# 9. Output distance and energy-delta to file (i.e track progress).
# 10. Repeat steps 1-9, increasing delta magnitude M each time.
