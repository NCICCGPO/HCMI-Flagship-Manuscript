import json
import pandas as pd
import sys
import argparse


class Lookup():
    def __init__(self,
                 file_out,
                 output_info_file,
                 goal_name='pairRawVcfInfo',
                 goal_struct='PairRawVcfInfo'):
        possible = ['PairRawVcfInfo', 'MergedPairVcfInfo']
        assert goal_struct in possible, 'goal_struct must be included in: ' + ' '.join(possible)
        self.goal_struct = goal_struct
        self.goal_name = goal_name
        self.file_out = file_out
        self.output_info = self.load_input(output_info_file)
        self.pair_ids = self.output_info['project_data']['pairId']
        self.sample_ids = self.output_info['project_data']['sampleIds']
        self.available = {}
        self.available_by_pair = {}
        self.available_by_sample = {}
        self.note_available_inputs()
        self.note_available_outputs()
        self.inputs = {}
        self.run_lookup()
        self.write_wdl_json()
        
    def run_lookup(self):
        if self.goal_struct == 'PairRawVcfInfo':
            self.lookup_pair_raw_vcf_info()
        elif self.goal_struct == 'MergedPairVcfInfo':
            self.lookup_merged_vcf_info()
            
    def lookup_merged_vcf_info(self):
        self.inputs['mergedPairVcfInfos'] = []
        for pair_id in self.pair_ids:
            pair_info = [pair_info for pair_info in self.output_info['project_data']['pairInfos'] 
                                 if pair_info['pairId'] == pair_id][0]
            try:
                pair_merged_info = {'pairId' : pair_id,
                                     'unannotatedVcf' : self.available_by_pair[pair_id]['mergedVcf'],
                                     'tumorId' : pair_info['tumorId'],
                                     'normalId' : pair_info['normalId']
                                     }
            except KeyError:
                pair_merged_info = {'pairId' : pair_id,
                                 'unannotatedVcf' : self.available_by_pair[pair_id]['mergedVcfs'],
                                 'tumorId' : pair_info['tumorId'],
                                 'normalId' : pair_info['normalId']
                                 }
            self.inputs['mergedPairVcfInfos'].append(pair_merged_info)
            
    def lookup_pair_raw_vcf_info(self):
#             pairId : pairRelationship.pairId,
#             filteredMantaSV : Calling.filteredMantaSV,
#             strelka2Snv : Calling.strelka2Snv,
#             strelka2Indel : Calling.strelka2Indel,
#             mutect2 : Calling.mutect2,
#             lancet : Calling.lancet,
#             svabaSv : Calling.svabaSv,
#             svabaIndel : Calling.svabaIndel,
#             tumorId : pairRelationship.tumorId,
#             normalId : pairRelationship.normalId,
#             tumorFinalBam : Preprocess.finalBam[tumorGetIndex.index],
#             normalFinalBam :

        self.inputs['pairRawVcfInfos'] = []
        for pair_id in self.pair_ids:
            pair_info = [pair_info for pair_info in self.output_info['project_data']['pairInfos'] 
                                 if pair_info['pairId'] == pair_id][0]
            if 'pairRawVcfInfos' in self.available_by_pair[pair_id]:
                self.inputs['pairRawVcfInfos'].append(self.available_by_pair[pair_id]['pairRawVcfInfos'])
            else:
                pair_raw_vcf_info = {'pairId' : pair_id,
                                     'filteredMantaSV' : self.available_by_pair[pair_id]['filteredMantaSV'],
                                     'strelka2Snv' : self.available_by_pair[pair_id]['strelka2Snv'],
                                     'strelka2Indel' : self.available_by_pair[pair_id]['strelka2Indel'],
                                     'mutect2' : self.available_by_pair[pair_id]['mutect2'],
                                     'lancet' : self.available_by_pair[pair_id]['lancet'],
                                     'svabaSv' : self.available_by_pair[pair_id]['svabaSv'],
                                     'svabaIndel' : self.available_by_pair[pair_id]['svabaIndel'],
                                     'tumorId' : pair_info['tumorId'],
                                     'normalId' : pair_info['normalId'],
                                     'tumorFinalBam' : pair_info['tumorFinalBam'],
                                     'normalFinalBam' : pair_info['normalFinalBam']
                                     }
                self.inputs['pairRawVcfInfos'].append(pair_raw_vcf_info)
            
    def write_wdl_json(self):
        with open(self.file_out, 'w') as input_info_file:
            json.dump(self.inputs, input_info_file, indent=4)

    def load_input(self, output_info_file):
        with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
        return output_info
    
    def test_duplicate(self, key, 
                       sample_id=False,
                       pair_id=False):
        if not sample_id:
            sample_id = self.sample_ids[0]
        if not pair_id:
            try:
                pair_id = self.pair_ids[0]
            except IndexError:
                pass
        if pair_id:
            assert ((key not in self.available) 
                    and (key not in self.available_by_sample[sample_id]) 
                    and (key not in self.available_by_pair[pair_id])), key + ' is in both inputs and outputs'
        else:
            assert ((key not in self.available) 
                    and (key not in self.available_by_sample[sample_id])), key + ' is in both inputs and outputs'
    
    def note_available_inputs(self):
        '''Load available variables from prior input'''
        for key in self.output_info['inputs']:
            assert key not in self.available, key + ' is in both inputs and outputs'
            self.available[key.split('.')[-1]] = self.output_info['inputs'][key]
            
    def note_available_outputs(self):
        '''Load available variables from prior output'''
        for pair_id in self.pair_ids:
            for key in self.output_info['pair_association'][pair_id]:
                assert key not in self.available, key + ' is in both inputs and outputs'
                if pair_id in self.available_by_pair:
                    assert key not in self.available_by_pair[pair_id], key + ' is in both inputs and outputs'
                else:
                    self.available_by_pair[pair_id] = {}
                self.available_by_pair[pair_id][key.split('.')[-1]] = self.output_info['pair_association'][pair_id][key][0]
        for sample_id in self.sample_ids:
            for key in self.output_info['sample_association'][sample_id]:
                assert key not in self.available , key + ' is in both inputs and outputs'
                if sample_id in self.available_by_sample:
                    assert key not in self.available_by_sample[sample_id], key + ' is in both inputs and outputs'
                else:
                    self.available_by_sample[sample_id] = {}
                self.available_by_sample[sample_id][key.split('.')[-1]] = self.output_info['sample_association'][sample_id][key][0]
       
       
    
def get_args():
    '''Parse input flags
        Need to add optional task-specific input json
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-file',
                        help='Output customInputs json file',
                        required=True
                       )
    parser.add_argument('--goal-name',
                        help='name of customInputs variable to populate',
                        required=False,
                        default='pairRawVcfInfo'
                       )
    parser.add_argument('--goal-struct',
                        help='Struct of customInputs variable to populate',
                        required=False,
                        default='PairRawVcfInfo'
                       )
    parser.add_argument('--output-info',
                        help='Prior outputInfo.json file. '
                        'Create by running collect.py',
                        required=True
                       )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__
         
                
def main():
    args = get_args()
    results = Lookup(file_out=args['out_file'],
                     output_info_file=args['output_info'],
                     goal_name=args['goal_name'],
                     goal_struct=args['goal_struct'])

    
if __name__ == "__main__":
    main()
