
import json
import logging as log
import pandas as pd
import argparse


class Merge():
    def __init__(self, run_info_files, file_out):
        self.file_out = file_out
        self.run_infos = [self.read(file) for file in run_info_files]
        self.write()
        
    @staticmethod
    def read(file):
        with open(file) as output_info_file:
            output_info = json.load(output_info_file)
            return output_info
        
    def write(self):
        with open(self.file_out, 'w') as project_info_file:
                json.dump(self.run_infos, project_info_file, indent=4)
        

def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-info-file',
                        help='JSON file with RunInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library',
                        nargs='*',
                        required=True
                        )
    parser.add_argument('--output',
                        help='Output file name',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    Merge(run_info_files=args['run_info_file'], 
          file_out=args['output'])

if __name__ == "__main__":
    main()      