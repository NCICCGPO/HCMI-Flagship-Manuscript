import pandas as pd
import numpy as np
import os
from google.cloud import storage
import make_auth

class Bicseq2Prep():
    ''''Prep the non-sample-specific and non-filesystem-specific 
    portion of the config files.
    Also returns the Bicseq2 WDL input variables that will be required for the run.'''
    def __init__(self, list_of_chroms_full,
                 inputs,
                 uniq_coords,
                 read_length,
                 upload_bucket,
                 create_config=True
                 ):
        self.list_of_chroms_full = list_of_chroms_full['object']
        self.read_length = self.match_read_length(uniq_coords, read_length)
        self.uniq_coords = uniq_coords['object'][self.read_length]
        self.inputs = inputs
        self.inputs['coordReadLength']['object'] = self.read_length
        if create_config:
            self.upload_bucket_uri = upload_bucket
            self.credentials, self.gcp_project = make_auth.login()
            self.write_sample_configs()
            self.write_pair_config()
        
    def write_sample_configs(self):
        '''create and upload configs for the normalization steps'''
        data = self.prep()
        file = 'sampleId.bicseq2.config'
        config = data.to_csv(sep='\t', index=False)
        uri = self.upload(file, config)
        self.inputs['bicseq2ConfigFile'] = uri
        
    def write_pair_config(self):
        file = 'pairId.bicseq2.seg.config'
        data = self.prep_pair()
        config = data.to_csv(sep='\t', index=False)
        uri = self.upload(file, config)
        self.inputs['bicseq2SegConfigFile'] = uri
            
    def parse_url(self, url):
        '''divide gcp bucket location into parts'''
        bucket_id = url.split('gs://')[-1].split('/')[0]
        project_id = '-'.join(bucket_id.split('-')[0:-1])
        name = '/'.join(url.split('gs://')[-1].split('/')[1:])
        return bucket_id, project_id, name
            
    def upload(self, file, config):
        '''upload file content (as a string) with a given filename'''
        bucket_id, project_id, name = self.parse_url(self.upload_bucket_uri)
        client = storage.Client(credentials=self.credentials,
                                project=self.gcp_project)
        bucket = client.get_bucket(bucket_id)
        blob = bucket.blob(file)
        blob.upload_from_string(config)
        return os.path.join(self.upload_bucket_uri, file)
        
    def match_read_length(self, uniq_coords, read_length):
        '''
        Make config file that includes the path
        to mappability file.
        '''
        length_key = False
        for length in uniq_coords['object']:
            if np.isclose(int(read_length), int(length), atol=2):
                length_key = length
                return length_key
            
    def prep_pair(self):
        ''' file: pairId.bicseq2.seg.config
        prep fasta-specific but sample independent portion of config file'''
        data = pd.DataFrame({'chr' : [chrom for chrom in self.list_of_chroms_full]})
        return data
        
    def prep(self):
        ''' sampleId.bicseq2.config
        prep fasta-specific but sample independent portion of config file'''
        data = pd.DataFrame({'chrom_name' : [chrom for chrom in self.list_of_chroms_full]})
        return data
