'''
Created on May 3, 2022

@author: jshelton
'''
import os
import sys
import json
import argparse
import pandas as pd
import jsonschema
from jsonschema import validate
import logging as log
import pprint
import Colorer
import parse

log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)

def load_pairs(file):
    pairs = pd.read_csv(file)
    assert 'tumorId' in pairs.columns and 'normalId' in pairs.columns, 'Error: can not find tumorId and normalId columns in pairs file: ' + ' '.join(pairs.columns)
    pairs['pairId'] = pairs.apply(lambda row: row.tumorId + '--' + row.normalId, axis=1)
    return pairs[['tumorId', 'normalId', 'pairId']], pairs

def load_sample_ids(file):
    '''Return a deduplicated list of sample ids'''
    sample_ids = pd.read_csv(file)
    assert 'sampleId' in sample_ids.columns , 'Error: can not find "sampleId" columns in samples file: ' + ' '.join(sample_ids.columns)
    return sample_ids.sampleId.unique().tolist()
    
    
def fill_pair_relationship(row):
    '''Add pair and sample level info to object
    PairRelationship '''
    pair_relationship = {'pairId' : row.pairId}
    pair_relationship['normalId'] = row['normalId']
    pair_relationship['tumorId'] = row['tumorId']
    pair_relationship['normalPrefix'] = row['normalId']
    pair_relationship['tumorPrefix'] = row['tumorId']
    return pair_relationship


def populate(args):
    run_info = {}
    run_info['project_name'] = args['project_name']
    run_info['workflow_uuid'] = args['uuid']
    run_info['project_data'] = {'listOfPairRelationships' : [],
                                        'sampleIds' : [],
                                        'pairIds' : [],
                                        'tumorIds' : [],
                                        'normalIds' : [],
                                }
    if args['pairs_file'] in args:
        pairs, full_pairs = load_pairs(args['pairs_file'])
        assert 'pairId' in pairs.columns , 'Error: can not find "pairId" columns in pairs file: ' + ' '.join(pairs.columns)
        assert 'tumorId' in pairs.columns , 'Error: can not find "tumorId" columns in pairs file: ' + ' '.join(pairs.columns)
        assert 'normalId' in pairs.columns , 'Error: can not find "normalId" columns in pairs file: ' + ' '.join(pairs.columns)
        pair_info_relationships = []
        # update basic pairing information with tumorId, normalId and pairId from pairs file
        pairs['listOfPairRelationships'] = pairs.apply(lambda row: fill_pair_relationship(row), axis=1)
        run_info['listOfPairRelationships'] = pairs['listOfPairRelationships'].tolist()
        run_info['pairIds'] = pairs.pairId.tolist()
        run_info['normalIds'] = pairs.normalId.tolist()
        run_info['tumorIds'] = pairs.tumorId.tolist()
    if args['samples_file']:
        run_info['sampleIds'] = load_sample_ids(args['samples_file'])
    return run_info
        


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser(description='Program uses the UUID and sample or pairs files to create a json')
    parser.add_argument('--uuid',
                        help='Cromwell workflow uuid',
                        required=True
                        )
    parser.add_argument('--output',
                        help='Output file name',
                        required=True
                        )
    parser.add_argument('--project-name',
                        help='Project name associated with account. If not supplied '
                        'define genome using --project-data',
                        required=True
                        )
    parser.add_argument('--pairs-file',
                        help='CSV file with items that are required to have '
                        '"tumorId", "normalId", "pairId" as columns '
                        ,
                        required=False
                        )
    parser.add_argument('--samples-file',
                        help='Not generally required. '
                        'If steps run only require sampleId and do not use pairing '
                        'information sample info can be populated with a CSV file. '
                        'The CSV file requires a columns named ["sampleId"]. If not supplied '
                        'define samples using --project-data',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__

def write_json(dictionary, file_out):
    with open(file_out, 'w') as info_file:
        json.dump(dictionary, info_file, indent=4)

def main():
    args = get_args()
    # project, run and sample related metadata (run_info)
    run_info = populate(args)
    file_out = args['output']
    # print run json
    write_json(dictionary=run_info,
               file_out=file_out)

if __name__ == "__main__":
    main()
