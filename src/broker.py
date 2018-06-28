from google.cloud import storage
from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds

bucket_name = 'funnel-folding.appspot.com'
pdb_file = '1ubq.pdb'
client = storage.Client()
bucket = client.bucket(bucket_name)
pdb = None

try:
    with open(f'data/{pdb_file}', 'w+') as f:
        print('Downloading PDB file...')

        pdb_blob = bucket.get_blob(f'pdb/{pdb_file}')
        pdb_blob.download_to_file(f)
        pdb = PDBFile('data/1ubq.pdb')
        if pdb_blob is None:
            print('ERROR: PDB file not found.  Exiting...')
            exit(0)
except:
    print('Error while downloading PDB file.  Exiting...')
    exit(0)



