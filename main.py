from Bio.PDB import *
from math import pi

# BioPython setup
parser = PDBParser()
ppb = PPBuilder()

# Load PDB file and verify we're working with one model and one chain
structure = parser.get_structure('ubiq', '1ubq.pdb')
models = [model for model in structure.get_models()]
chains = [model for model in structure.get_chains()]
if len(models) != 1 or len(chains) != 1:
    print('ERROR: Only PDBs with one model and one chain are currently supported!\nExiting...')
    exit()

# Remove hetero residues (such as solvent)
model = models[0]
chain = chains[0]
residue_list = [res for res in chain.get_residues()]
hetero_residue_ids = []

for res in residue_list:
    if res.id[0] != ' ':
        hetero_residue_ids.append(res.id)

print(f'Number of residues in chain: {len(residue_list)}')
print(f'Removing {len(hetero_residue_ids)} hetero residues')
[chain.detach_child(hetero_id) for hetero_id in hetero_residue_ids]

# Update residues after removing hetero residues
residues = structure.get_residues()
residue_list = [res for res in residues]
print(f'Chain now has {len(residue_list)} residues')

# Build a polypeptide from our structure and get the torsion angles
polypeptides = ppb.build_peptides(structure)
if len(polypeptides) != 1:
    print(f'ERROR: expected only one polypeptide from PDB structure')
    exit()

polypeptide = polypeptides[0]
torsion_angles = polypeptide.get_phi_psi_list()

for atom in residue_list[0].get_atoms():
    print(atom.get_coord())
    print(atom.get_vector())
residue = residue_list[0]

# From Biopython FAQ
nitro = residue['N'].get_vector()
carbon = residue['C'].get_vector()
carbon_a = residue['CA'].get_vector()

# center at origin
nitro = nitro - carbon_a
carbon = carbon - carbon_a

# find rotation matrix that rotates n -120 degrees along the ca-c vector
rot = rotaxis(-pi * 120.0/180.0, carbon)

# apply rotation to ca-n vector
cb_at_origin = nitro.left_multiply(rot)

# put on top of ca atom
cb = cb_at_origin + carbon_a

print(f'rotated {cb}')

for atom in residue_list[0].get_atoms():
    print(atom.get_coord())
    print(atom.get_vector())

