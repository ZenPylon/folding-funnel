from Bio.PDB import *

parser = PDBParser()
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

print(residue_list)
for res in residue_list:
    print(res)
    if res.id[0] != ' ':
        hetero_residue_ids.append(res.id)

print(f'Number of residues in chain: {len(residue_list)}')
print(f'Removing {len(hetero_residue_ids)} hetero residues')
[chain.detach_child(hetero_id) for hetero_id in hetero_residue_ids]

# Update residues after removing hetero residues
residues = structure.get_residues()
residue_list = [res for res in residues]
print(f'Chain now has {len(residue_list)} residues')

