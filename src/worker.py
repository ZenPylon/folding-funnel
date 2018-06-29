import requests
import pickle
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

molecule = MoleculeUtil(AppSettings.local_pdb_path)

for i in range(100):
    job = requests.get(f'{AppSettings.host}/job_request')
    new_torsions = pickle.loads(job.content)
    molecule.set_torsions(new_torsions)
    molecule.run_simulation()

