import numpy as np
import pickle
from google.cloud import storage
from molecule_util import get_modeller, calculate_zmat, get_torsion_indices
from flask import Flask, jsonify

np.random.seed(20)

bucket_name = 'funnel-folding.appspot.com'
project_id = 'funnel_folding'
pdb_file = '1ubq.pdb'

storage_client = storage.Client()
bucket = storage_client.bucket(bucket_name)
pdb = None
modeller = None

# Load the PDB file and construct the modeller and zmatrix
try:
    with open(f'data/{pdb_file}', 'wb+') as f:
        print('Downloading PDB file...')
        pdb_blob = bucket.get_blob(f'pdb/{pdb_file}')
        if pdb_blob is None:
            print('ERROR: PDB file not found.  Exiting...')
            exit(0)
        pdb_blob.download_to_file(f)
        f.close()
except Exception as e:
    print(f'Error while initializing PDB file: \n{e}  \nExiting...')
    exit(0)

modeller = get_modeller(f'data/{pdb_file}')
zmat = calculate_zmat(modeller)
torsion_indices = get_torsion_indices(zmat)
starting_torsions = np.array([zmat.loc[torsion_indices[:, 0], 'dihedral'],
                              zmat.loc[torsion_indices[:, 1], 'dihedral']]).T
print(starting_torsions)


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
    print('Starting torsions\n\n')
    print(starting_torsions)
    print('\n offsets \n')
    print(offsets)
    new_torsions = np.zeros(shape=starting_torsions.shape)
    total_offset = offset_step * num_jobs_requested
    print('\nnew torsions\n')
    new_torsions[:, 0] = starting_torsions[:, 0] + (offsets[:, 0] * total_offset)
    new_torsions[:, 1] = starting_torsions[:, 1] + (offsets[:, 1] * total_offset)
    print(new_torsions)
    
    return pickle.dumps(new_torsions)

    # TODO - put this in the worker process
    # for i in range(num_configs):
    #     zmat.safe_loc[torsion_indices[:, 0], 'dihedral'] = \
    #         starting_torsions[:, 0] + (offsets[:, 0] * i * offset_size)
    #     zmat.safe_loc[torsion_indices[:, 1], 'dihedral'] = \
    #         starting_torsions[:, 1] + (offsets[:, 1] * offset_size * i)


@app.route('/mark_finished', methods=['POST'])
def handle_job_finished():
    num_jobs_completed += 1
    print(f'Finished job #{num_jobs_completed}')


# Init worker processes
