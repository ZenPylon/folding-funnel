from google.cloud import storage


class CloudUtil(object):
    """
    A utility class that downloads pdb files from cloud soorage

    Attributes:
        bucket_name (str):
        project_id (str): the google 
    """
    def __init__(self, project_id, bucket_name):
        self.bucket_name = bucket_name
        self.project_id = project_id
        self.storage_client = storage.Client()

    def download_pdb(self, cloud_path: str, local_path: str):
        """
        Downloads a pdb file from cloud storage and saves it locally

        Args:
            cloud_path (str): path to pdb file in storage bucket
            local_path(str): file to be written locally
        """
        bucket = self.storage_client.bucket(self.bucket_name)
        pdb = None

        # Load the PDB file and construct the modeller and zmatrix
        with open(local_path, 'wb+') as f:
            print('Downloading PDB file...')
            pdb_blob = bucket.get_blob(cloud_path)
            if pdb_blob is None:
                print('ERROR: PDB file not found.  Exiting...')
                exit(0)
            pdb_blob.download_to_file(f)
            f.close()