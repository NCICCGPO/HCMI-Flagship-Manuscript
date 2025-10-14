import subprocess
import sys
import logging as log
import pandas as pd
import json
import os
import numpy as np
import copy
from collections.abc import Iterable
import pprint
import argparse

class JsonModify():
    
    def __init__(self, inputs,
                 dir,
                 convert_path=True):
        '''
        deep copy required to modify complex outputs structure while
        keeping an original
        tool options:

        convert_path: swaps path betore the task call- directory out with
        the final output directory
        '''
        self.convert_path_var = convert_path
        self.files = []
        self.json_obj = copy.deepcopy(inputs)
        self.final_json_obj = copy.deepcopy(inputs)
        self.dir = dir
        # gather all uris in the json_obj
        # or process the json object with choosen tool
        self.hlpr_fnc(self.json_obj, self.final_json_obj)
        
    def convert_path(self, string):
        ''' Old version of 
        Convert from intermediate file URI to final URI of the final output.
        This is '''
        if not string.startswith('gs://'):
            return string
        else:
            if self.convert_path_var: 
                uri = string
                basename = '/'.join(uri.replace('/cacheCopy', '').split('/call-')[-1].split('/')[1:])
                if self.dir.endswith('/'):
                    new_uri = self.dir + basename
                else:
                    new_uri = self.dir + '/' + basename
                new_uri = '/'.join([section for section in new_uri.split('/') 
                                   if not section.startswith('attempt-') ])
                self.files.append(new_uri)
                return new_uri
            else:
                uri = string
                return uri.replace('gs://', 'gs:/').replace('//', '/').replace('gs:/', 'gs://')
    
    def run_tool(self, subnode):
        '''run the function for subnode'''
        result = self.convert_path(subnode)
        return result
        
    def hlpr_fnc_list(self, list1, description):
        '''replace for list
            hlpr_fnc(list_obj, final_json_obj)
        '''
        for i, subnode in enumerate(list1):
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[i])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[i])
            else:
                if isinstance(subnode, str):
                    result = self.run_tool(subnode)
                    if result:
                        description[i] = result

    def hlpr_fnc(self, dict1, description):
        '''replace for dictionary
            hlpr_fnc(json_obj, final_json_obj)
        '''
        nodes = dict1.keys()
        for node in nodes:
            subnode = dict1[node]
            if isinstance(subnode, dict):
                self.hlpr_fnc(subnode, description[node])
            elif not isinstance(subnode, str) and isinstance(subnode, Iterable):
                self.hlpr_fnc_list(subnode, description[node])
            else:
                if isinstance(subnode, str):
                    result = self.run_tool(subnode)
                    if result:
                        description[node] = result


class Lookup():
    def __init__(self,
                 output_info_file,
                 final_workflow_outputs_dir=False):
        '''gather files from the outputInfo 
            sample_association and pair_association section'''
        self.output_info = self.load_input(output_info_file)
        self.pair_ids = False
        self.sample_ids = False
        self.analysis_groups = False # add support for analysis_groups
         
        if 'pairIds' in self.output_info['project_data']:
            self.pair_ids = self.output_info['project_data']['pairIds']
        if 'sampleIds' in self.output_info['project_data']:
            self.sample_ids = self.output_info['project_data']['sampleIds']
        if 'analysisGroups' in self.output_info['project_data'] \
                and len(self.output_info['project_data']['analysisGroups']) > 0:
            self.analysis_groups = self.output_info['project_data']['analysisGroups']
        if not final_workflow_outputs_dir:
            try:
                self.options = self.output_info['project_data']['options']
                self.final_workflow_outputs_dir = self.options['final_workflow_outputs_dir']
            except KeyError:
                log.error('"options" must be specified with a flag or in the outputInfo file')
        else:
            self.final_workflow_outputs_dir = final_workflow_outputs_dir
            self.options = pd.DataFrame({})
        self.sample_association = False
        self.pair_association = False
        self.analysis_group_association = False
        self.outputs = self.remove_workflow(self.output_info['outputs'], outputs=True)
        if 'sample_association' in self.output_info:
            self.sample_association = self.remove_workflow(self.output_info['sample_association'])
        if 'pair_association' in self.output_info:
            self.pair_association = self.remove_workflow(self.output_info['pair_association'])
        if 'analysis_group_association' in self.output_info:
            self.analysis_group_association = self.remove_workflow(self.output_info['analysis_group_association'])
        self.searchable = {'outputs' : self.outputs,
                           'sample_association' : self.sample_association,
                           'pair_association' : self.pair_association,
                           'analysis_group_association' : self.analysis_group_association
                           }
        
    def finalize_uris(self, raw_outputs):
        '''replace workspace with final dir'''
        outputs = JsonModify(inputs=raw_outputs,
                             dir=self.final_workflow_outputs_dir,
                             convert_path=self.convert_path)
        self.renamed_outputs = outputs.final_json_obj
        self.renamed_files = outputs.files
        
    def lookup(self, searchable, variable, convert_path, id=False):
        '''return file object'''
        self.convert_path = convert_path
        if not searchable in self.searchable:
            log.error('possible objects to search are outputs, sample_association, pair_association, analysis_group_association')
            log.error('in this outputInfo file the following are defined: ' + ' '.join([s for s in self.searchable if s]))
            sys.exit(1)
        if self.searchable[searchable]:
            if searchable == 'outputs':
                if variable in self.searchable[searchable]:
#                     if self.convert_path:
                    self.finalize_uris(self.searchable[searchable], self.convert_path)
                    return self.renamed_outputs[variable]
#                     else:
#                         return self.searchable[searchable][variable]
            else:
                try:
                    if variable in self.searchable[searchable][id]:
                        self.finalize_uris(self.searchable[searchable][id])
                        return self.renamed_outputs[variable]
                except KeyError:
                    return False
        return False
    
    def remove_workflow(self, dictionary, outputs=False):
        ''' remove workflow prefix because multiple workflows 
        may create the same output (e.g. start from bam vs start from fastq
        but produce the same vcfs)'''
        new = {}
        if outputs:
            for var_name in dictionary:
                new[var_name.split('.')[-1]] = dictionary[var_name]
        else:
            for id in dictionary:
                new[id] = {}
                for var_name in dictionary[id]:
                    new[id][var_name.split('.')[-1]] = dictionary[id][var_name]
        return new
        
    def load_input(self, output_info_file):
        with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
        return output_info


class Make_GDC_inputs():
    def __init__(self,
                 vcf_map,
                 pairs_file,
                 output_info_file,
                 file_out,
                 convert_path=False,
                 final_workflow_outputs_dir=False,
                 gdc_gcp_bucket="gs://nygc-comp-s-fd4e-input/gdc_vcf/"):
        ''' nygcVcf
            gdcSnvVcf
            gdcIndelVcf
            gdcSvVcf
            nygcFinalBedPe
        '''
        self.convert_path = convert_path
        self.gdc_gcp_bucket = gdc_gcp_bucket
        self.pairs = pd.read_csv(pairs_file)
        if not 'pair_id' in self.pairs.columns:
            self.pairs['pair_id'] = self.pairs.apply(lambda row: row.tumor + '--' + row.normal, axis=1)
        self.looked_up = Lookup(output_info_file,
                                final_workflow_outputs_dir=final_workflow_outputs_dir)
        self.parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        self.vcf_map = pd.read_csv(vcf_map)
        if not 'pair_id' in self.vcf_map.columns:
            self.vcf_map['pair_id'] = self.vcf_map.apply(lambda row: row.tumor + '--' + row.normal, axis=1)
        self.vcf_map.fillna({'caveman': False,
                             'pindel': False}, inplace=True)
        self.file_out = file_out
        self.gdc_relationships = []
        for index, row in self.pairs.iterrows():
            self.fill_in(row)
        self.inputs = { "gdcRelationships" : self.gdc_relationships}
        self.write_wdl_json()
        
    def fill_in(self, row):
        ''' nygcVcf
            gdcSnvVcf
            gdcIndelVcf
            gdcSvVcf
            nygcFinalBedPe
        '''
        finalVcfPairInfos = self.looked_up.lookup(searchable='pair_association',
                              variable='finalVcfPairInfos',
                              id=row.pair_id,
                              convert_path=self.convert_path)
        if finalVcfPairInfos:
            for finalVcfPairInfo in finalVcfPairInfos:
                if 'mainVcf' in finalVcfPairInfo:
                    nygc_vcf = finalVcfPairInfo['mainVcf']
                    nygcFinalBedPe = finalVcfPairInfo['svHighConfidenceFinalBedPe']
            match = self.vcf_map[self.vcf_map.pair_id == row.pair_id]
            caveman_vcf = match.caveman_vcf.tolist()[0]
            pindle_vcf = match.pindle_vcf.tolist()[0]
            brass_vcf = match.brass_vcf.tolist()[0]
            if caveman_vcf \
                    and pindle_vcf \
                    and brass_vcf:
                gdcSnvVcf = self.gdc_gcp_bucket + caveman_vcf.split('/')[-1]
                gdcIndelVcf = self.gdc_gcp_bucket + pindle_vcf.split('/')[-1]
                gdcSvVcf = self.gdc_gcp_bucket + brass_vcf.split('/')[-1]
                gdc_relationship = {
                          "gdcSnvVcf" : gdcSnvVcf,
                          "tumor" : row.tumor,
                          "normal" : row.normal,
                          "gdcIndelVcf" : gdcIndelVcf,
                          "pairId" : row.pair_id,
                          "nygcVcf" : nygc_vcf,
                          "nygcFinalBedPe" : nygcFinalBedPe,
                          "gdcSvVcf" : gdcSvVcf
                          }
                self.gdc_relationships.append(gdc_relationship)
    
    def write_wdl_json(self):
        with open(self.file_out, 'w') as input_info_file:
            json.dump(self.inputs, input_info_file, indent=4)


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf-map',
                        help='CSV file with tumor, normal, '
                        'caveman_vcf, pindle_vcf and brass_vcf '
                        '(fix pindel typo later)',
                        required=True
                        )
    parser.add_argument('--file-out',
                        help='Out customInput  JSON file name.',
                        required=True
                        )
    parser.add_argument('--pairs-file',
                        help='CSV file with tumor, normal',
                        required=True
                        )
    parser.add_argument('--output-info',
                        help='outputInfo JSON file from post nygc_pipeline metrics',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    Make_GDC_inputs(args['vcf_map'],
                    args['pairs_file'],
                    args['output_info'],
                    args['file_out'],
                    gdc_gcp_bucket="gs://nygc-comp-s-fd4e-input/gdc_vcf/")


main()



