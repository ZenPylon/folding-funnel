import pickle
from flask import Flask, request
from settings import AppSettings
from cloud_util import CloudUtil
from molecule_util import MoleculeUtil


cloud_util = CloudUtil(AppSettings.project_id, AppSettings.bucket_name)
try:
    cloud_util.download_pdb(AppSettings.cloud_pdb_path,
                            AppSettings.local_pdb_path)
except Exception as e:
    print(f'Error while initializing PDB file: \n{e}  \nExiting...')
    exit(0)

molecule = MoleculeUtil(AppSettings.local_pdb_path)
app = Flask(__name__)
num_jobs_requested = 0
num_jobs_completed = 0


@app.route('/job_request', methods=['GET'])
def handle_job_request():
    """
    Responds to a worker process's request for a job.
    Returns a new set of angles for the worker's simulation to start with.
    """
    global num_jobs_requested

    new_torsions = molecule.get_offset_torsions(num_jobs_requested)
    num_jobs_requested += 1
    return pickle.dumps(new_torsions)


@app.route('/complete_job', methods=['POST'])
def handle_job_finished():
    global num_jobs_completed
    num_jobs_completed += 1
    print(f'Finished job {num_jobs_completed}')
    print(pickle.loads(request.data))



# Init worker processes
