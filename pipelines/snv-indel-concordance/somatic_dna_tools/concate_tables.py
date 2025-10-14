#!/usr/bin/env python 

import argparse
import pandas as pd
import csv


class Concate():
    
    def __init__(self, 
                 files,
                 output):
        '''Quickly concate any number of tables as long as all files share identical headers'''
        simple = self.test(files)
        if simple:
            self.concate(files, output)
        else:
            self.slow_concate(files, output)
        
    def test(self, files):
        '''Confirm that headers all match'''
        headers = []
        for file in files:
            with open(file) as f:
                header = f.readline()
                headers.append(header)
        if len(set(headers)) == 1:
            return True
        return False
        
    def slow_concate(self, files, output):
        with open(files[0]) as file:
            first_line = file.readline()
        s = csv.Sniffer()
        delimiter = s.sniff(first_line).delimiter
        dfs = []
        for file in files:
            df = pd.read_csv(file, sep=delimiter)
            dfs.append(df)
        data = pd.concat(dfs, axis=0, ignore_index=True)
        data.to_csv(output, sep=delimiter, index=False)
    
    def concate(self, files, output):
        '''Quickly concate any number of tables as long as all files share identical headers'''
        with open(output,"wb") as f:
            # first file:
            with open(files[0], "rb") as fin:
                f.write(fin.read())
            # now the rest:    
            for file in files[1:]:
                with open(file, "rb") as fin:
                    next(fin) # skip the header
                    f.write(fin.read())
                
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--tables',
                        help='All table files',
                        required=True,
                        nargs='+'
                       )
    parser.add_argument('--output',
                        help='Output CSV file',
                        required=True
                       )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    Concate(files=args['tables'],
            output=args['output'])


if __name__ == "__main__":
    main() 