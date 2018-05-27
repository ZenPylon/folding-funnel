from Bio.PDB import *
from math import pi

def remove_hetero(model, chain):
    """Remove heteroatoms from a chain and return the new residue_list"""
    residue_list = [res for res in chain.get_residues()]
    hetero_residue_ids = []
    for res in residue_list:
        if res.id[0] != ' ':
            hetero_residue_ids.append(res.id)

    print(f'Number of residues in chain: {len(residue_list)}')
    print(f'Removing {len(hetero_residue_ids)} hetero residues')
    [chain.detach_child(hetero_id) for hetero_id in hetero_residue_ids]

    # Update residues after removing hetero residue
    residues = chain.get_residues()
    residue_list = [res for res in residues]
    print(f'Chain now has {len(residue_list)} residues')
    return residue_list

def main_load(model_name, filename):
    """
    Loads a pdb file and ensures it has one model and chain.
    Returns structure, chain, and polypeptide.

    This function will likely change
    """
    parser = PDBParser()
    ppb = PPBuilder()

    # Load PDB file and verify we're working with one model and one chain
    structure = parser.get_structure(model_name, filename)
    models = [model for model in structure.get_models()]
    chains = [model for model in structure.get_chains()]
    if len(models) != 1 or len(chains) != 1:
        print('ERROR: Only PDBs with one model and one chain are currently supported!\nExiting...')
        exit()

    # Remove hetero residues (such as solvent)
    model = models[0]
    chain = chains[0]
    residue_list = remove_hetero(model, chain)
    print(len(residue_list))

    # Build a polypeptide from our structure and get the torsion angles
    polypeptides = ppb.build_peptides(structure)
    if len(polypeptides) != 1:
        print(f'ERROR: expected only one polypeptide from PDB structure')
        exit()

    polypeptide = polypeptides[0]

    return structure, residue_list, polypeptide


structure, residue_list, polypeptide = main_load('ubiq', '1ubq.pdb')
torsion_angles = polypeptide.get_phi_psi_list()

# TODO:
# 1. Specify all phi / psi angle deltas (magnitude M) relative to native structure
#    and get absolute phi / psi coordinates (calculate native angles + deltas).
# 2. Convert from angle-space to cartesian coordinates.
# 3. Pass coordinates to OpenMM (how?).
# 4. Add solvent to model with OpenMM.
# 5. Run energy minimization on model.
# 6. Pass final energy value and coordinates to BioPython, and reset atom coords.
# 7. Recalculate phi / psi angles from coordinates.
# 8. Calculate distance in angle-space from starting point (in this case,
#    native structure) and calculate delta in energy level.
# 9. Output distance and energy-delta to file (i.e track progress).
# 9. Repeat steps 1-8, increasing delta magnitude M each time.
