from Bio.PDB import *
from simtk.openmm.app import PDBFile
from math import pi
import pandas as pd

def remove_hetero(model, chain):
    """Remove heteroatoms from a chain and return the new residue_list"""
    residue_list = [res for res in chain.get_residues()]
    hetero_residue_ids = []
    for res in residue_list:
        if res.id[0] != ' ':
            hetero_residue_ids.append(res.id)
    
    print(f'Number of residues in chain: {len(residue_list)}')
    print(f'Removing {len(hetero_residue_ids)} hetero residues')
    for hetero_id in hetero_residue_ids:
        chain.detach_child(hetero_id)

def remove_hydrogens(model, chain):
    """Remove hydrogens from a chain and return the new residue_list"""
    hydrogen_count = 0
    for res in chain.get_residues():
        hydrogens = [atom.get_name() for atom in res.get_atoms() if atom.element == 'H']
        for hydrogen in hydrogens:
            res.detach_child(hydrogen) 
        hydrogen_count += len(hydrogens)

    print(f'Removed {hydrogen_count} hydrogen atoms')
        
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
        # exit()

    # Remove hetero residues (such as solvent)
    model = models[0]
    chain = chains[0]
    remove_hetero(model, chain)
    remove_hydrogens(model, chain)
    residues = chain.get_residues()
    residue_list = [res for res in residues]
    print(f'Chain now has {len(residue_list)} residues')

    # Build a polypeptide from our structure and get the torsion angles
    polypeptides = ppb.build_peptides(structure)
    if len(polypeptides) != 1:
        print(f'ERROR: expected only one polypeptide from PDB structure')
        exit()

    polypeptide = polypeptides[0]

    return structure, residue_list, polypeptide
