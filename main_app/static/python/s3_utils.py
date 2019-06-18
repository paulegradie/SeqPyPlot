import boto3
import os

def upload_to_s3(file_to_upload):
    s3 = boto3.client('s3', region='us-east-2')
    key = os.path.join('SPLOT_USER_EXPERIMENTS', os.path.basename(file_to_upload))
    s3.upload_file(
        file_to_upload,
        SPPLOT_BUCKET,
        key,
        ExtraArgs={
            'ACL': 'public-read',
            'ContentType': "application/zip"
        }
    )