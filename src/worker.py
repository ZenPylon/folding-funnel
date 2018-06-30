import requests
import pickle
import numpy as np
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

for i in range(2):
    job = requests.get(f'{AppSettings.host}/job_request')
    new_torsions = pickle.loads(job.content)
    molecule.set_torsions(new_torsions)

    p_energy, positions = molecule.run_simulation()
    print(p_energy)
    molecule.set_cc_positions(positions)
    torsions = molecule.get_torsions()
    angle_diff = torsions.flatten() - molecule.starting_torsions.flatten()
    print(angle_diff)
    angle_dist = np.sqrt(np.dot(angle_diff, angle_diff))
    requests.post(f'{AppSettings.host}/complete_job',
                  data=pickle.dumps((p_energy, angle_dist)))
