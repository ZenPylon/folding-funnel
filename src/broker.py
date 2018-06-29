from google.cloud import storage
from google.cloud.pubsub_v1.gapic.publisher_client import PublisherClient
from google.cloud.pubsub_v1.gapic.subscriber_client import SubscriberClient
from molecule_util import get_modeller, calculate_zmat

bucket_name = 'funnel-folding.appspot.com'

publisher = PublisherClient()
subscriber = SubscriberClient()
project_id = 'funnel_folding'

pdb_file = '1ubq.pdb'
client = storage.Client()
bucket = client.bucket(bucket_name)
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
print(zmat)
# Listen for worker events
subscriber.create_subscription(
    f'projects/{project_id}/topics/init_requests',
    f'projects/{project_id}/subscriptions/broker_init_request',
)

# Init worker processes



