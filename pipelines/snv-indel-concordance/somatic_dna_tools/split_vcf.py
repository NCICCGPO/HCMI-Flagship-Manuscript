#!/usr/bin/env python

import argparse
import pandas as pd
import csv
import math
import numpy as np
import os


class Split():
    
    def __init__(self,
                 file,
                 out_prefix,
                 out_suffix='.vcf',
                 max_rows=1000,
                 min_splits=40,
                 max_splits=50):
        '''Quickly split to any number of tables'''
        self.min_splits = min_splits
        self.max_splits = max_splits
        self.max_rows = max_rows
        self.out_prefix = out_prefix
        self.out_suffix = out_suffix
        self.load_header(file)
        self.variant_count = self.count_variants(file)
        self.chunk_method = self.chunk_method()
        print(self.chunk_method)
        if self.chunk_method == 'max_rows':
            self.intervals = self.get_intervals(start=0,
                                            end=self.variant_count,
                                            step=self.max_rows)
        elif self.chunk_method == 'min_splits':
            step = max(math.floor(self.variant_count / self.min_splits), 1)
            self.intervals = self.get_intervals(start=0,
                                            end=self.variant_count,
                                            step=step)
        self.write_vcf(file)
    
    def write_vcf(self, file):
        feed = self.feed_vcf(file)
        i = 0
        mnv_count = 0
        mnv_file = self.out_prefix + ".mnv" + self.out_suffix
        with open(mnv_file, 'w') as mnv:
            mnv.write(self.header)
            for start, end in self.intervals:
                i += 1
                out_file = self.out_prefix + '.' + str(i) + self.out_suffix
                with open(out_file, 'w') as out:
                    out.write(self.header)
                    for k in list(range(start, end)):
                        line = next(feed)
                        if 'MNV_ID' in line:
                            mnv.write(line)
                            mnv_count += 1
                        else:
                            out.write(line)
        if mnv_count == 0:
            os.remove(mnv_file)
    
    def feed_vcf(self, file):
        with open(file) as input:
            for line in input:
                if not line.startswith('#'):
                    yield line
        
    def get_intervals(self, start, end, step):
        intervals = [(s, (s + step if s + step < end else end))
                    for s in range(start, end, step)]
        return intervals
    
    def chunk_method(self):
        # splits given min_splits
        # divide into n files (plus partial file)
        file_count = self.min_splits
        if self.variant_count < self.min_splits:
            file_count = self.variant_count
        # splits given max_rows
        # split into files with i number of lines in most files
        max_row_file_count = math.ceil((self.variant_count) / self.max_rows)
        # if max rows gives more splits than min splits
        # but not too many splits for max rows...
        if (file_count > max_row_file_count) or (self.max_splits < max_row_file_count):
            self.min_splits = self.max_splits
            return 'min_splits'
        return 'max_rows'
    
    def count_variants(self, file):
        count = 0
        with open(file) as input:
            for line in input:
                if not line.startswith('#'):
                    count += 1
        return count
    
    def load_header(self, file):
        '''load vcf header'''
        self.header = ''
        with open(file) as input:
            for line in input:
                if line.startswith('#'):
                    self.header += line
                else:
                    return True


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf',
                        help='Input vcf file',
                        required=True
                       )
    parser.add_argument('--output-prefix',
                        help='Output file prefix',
                        required=True
                       )
    parser.add_argument('--output-suffix',
                        help='Output file suffix',
                        default='.vcf'
                       )
    parser.add_argument('--max-rows',
                        help='number of rows per output files',
                        type=int,
                        default='30'
                       )
    parser.add_argument('--min-splits',
                    help='Min number of output files to create',
                    type=int,
                    default='30'
                   )
    parser.add_argument('--max-splits',
                help='Max number of output files to create',
                type=int,
                default='30'
               )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    Split(file=args['vcf'],
         out_prefix=args['output_prefix'],
         out_suffix=args['output_suffix'],
         max_rows=args['max_rows'],
         min_splits=args['min_splits'],
         max_splits=args['max_splits'])


if __name__ == "__main__":
    main()
