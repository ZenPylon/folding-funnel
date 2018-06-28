from google.cloud import storage
from molecule_util import get_modeller, calculate_zmat

bucket_name = 'funnel-folding.appspot.com'
pdb_file = '1ubq.pdb'
client = storage.Client()
bucket = client.bucket(bucket_name)
pdb = None
modeller = None

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
print(zmat)
