import pickle
from flask import Flask, jsonify
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

ubiq_molecule = MoleculeUtil(AppSettings.local_pdb_path)

# 
# Setup webserver
#

num_jobs_requested = 0
num_jobs_completed = 0
offsets = np.random.choice([0, 0, -1, 1], starting_torsions.shape)
offset_step = 4
app = Flask(__name__)


@app.route('/init_request', methods=['GET'])
def handle_init_request():
    return pickle.dumps(starting_torsions)


@app.route('/job_request', methods=['GET'])
def handle_job_request():
    """
    Responds to a worker process's request for a job.
    Returns a new set of angles for the worker's simulation to start with.
    """
    global num_jobs_requested
    global offset_step

    num_jobs_requested += 1
    new_torsions = np.zeros(shape=starting_torsions.shape)
    total_offset = offset_step * num_jobs_requested
    new_torsions[:, 0] = starting_torsions[:, 0] + (offsets[:, 0] * total_offset)
    new_torsions[:, 1] = starting_torsions[:, 1] + (offsets[:, 1] * total_offset)

    return pickle.dumps(new_torsions)


@app.route('/mark_finished', methods=['POST'])
def handle_job_finished():
    global num_jobs_completed
    num_jobs_completed += 1
    print(f'Finished job #{num_jobs_completed}')


# Init worker processes
