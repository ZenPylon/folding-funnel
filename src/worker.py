from cloud_util import CloudUtil
from molecule_util import MoleculeUtil
from settings import AppSettings


cloud_util = CloudUtil(AppSettings.project_id, AppSettings.bucket_name)

try:
    cloud_util.download_pdb(AppSettings.cloud_pdb_path,
                            AppSettings.local_pdb_path)
except Exception as e:
    print(f'Error while initializing PDB file: \n{e}  \nExiting...')
    exit(0)

ubiq_molecule = MoleculeUtil(AppSettings.local_pdb_path)

# TODO - put this in the worker process
# for i in range(num_configs):
#     zmat.safe_loc[torsion_indices[:, 0], 'dihedral'] = \
#         starting_torsions[:, 0] + (offsets[:, 0] * i * offset_size)
#     zmat.safe_loc[torsion_indices[:, 1], 'dihedral'] = \
#         starting_torsions[:, 1] + (offsets[:, 1] * offset_size * i)
