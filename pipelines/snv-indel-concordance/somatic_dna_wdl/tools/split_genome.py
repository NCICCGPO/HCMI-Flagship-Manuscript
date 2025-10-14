import pandas as pd
import json
import sys
import subprocess
from subprocess import Popen, PIPE
import logging as log
import time

gap_file = "GRCh38.gaps.bed"
chrom_file = "wdl_port/config/fasta_references.json"
genome = "Human_GRCh38_full_analysis_set_plus_decoy_hla"
chrom_bed = "/gpfs/commons/groups/nygcfaculty/kancero/data/deconstructSigs/GRCh38_full_analysis_set_plus_decoy_hla.chrom_lengths.bed"
chrom_lengths_file = "/gpfs/commons/datasets/old-nygc-resources/GRCh38_full_analysis_set_plus_decoy_hla/internal/chrom_lengths.txt"
bins = 100


def run_bedtools(chrom_bed, gap_bed, nongap_bed):
    '''invert gap bed
    '''
    print(chrom_bed, gap_bed, nongap_bed)
    args = ['bedtools', 
                'subtract', 
                '-a', chrom_bed,
                '-b', gap_bed,
                '>', nongap_bed
           ]
    print(' '.join(args))
    run = Popen(' '.join(args), stdout=PIPE, shell=True) # split
    while True:
        line = run.stdout.readline().decode('ascii') # convert for python3
        if line != '':
            program = line.rstrip()
        if str(line) == '' and run.poll() != None:
            break
        run.wait()
    if run.returncode == 0:
        log.info('bedtools run complete')
    else:
        return False
    return True
    

# get chroms
with open(chrom_file) as fasta_refs:
    fasta_references = json.load(fasta_refs)

chroms = fasta_references[genome]['listOfChroms']['object']

chrom_lengths = pd.read_csv(chrom_lengths_file, sep='\t', header=None, 
                            names=['chrom', 'length_bp'])



# prep intervals
data = pd.read_csv(gap_file, sep='\t')
# filter for chroms
data = data[data.chrom.isin(chroms)].copy()
# filter for length
data = data[data['size'] > 100].copy()
# invert bed
bed = data[['chrom', 'chromStart', 'chromEnd', 'type', 'size']].copy()
out = 'GRCh38_large.gaps.bed'
with open(out, 'wb') as outfile:
    bed.to_csv(outfile, header=False, index=False, sep='\t')
nongap_bed='GRCh38_min_gap_100.bed'
time.sleep(3)
ran = run_bedtools(chrom_bed=chrom_bed, 
                   gap_bed=out, 
                   nongap_bed=nongap_bed)
# sort nongap_bed
def add_bin(intervals, i, bin):
    intervals[bin].append(i)
    if bin == bins:
        bin = 1
    else:
        bin += 1
    return intervals, bin


if ran:
    print('ran')
    data = pd.read_csv(nongap_bed, sep='\t', header=None, 
                       names=['chrom', 'chromStart', 'chromEnd'])
    data['type'] = 'minGap100Seg'
    data['size'] = data.apply(lambda row : row.chromEnd - row.chromStart, axis=1)
    data = data.sort_values(['size'], ascending=False).reset_index()
    intervals = {i : [] for i in range(1, bins + 1)}
    mean_length = data['size'].sum() / float(bins)
    print('mean_length', mean_length)
    # mean_length 27563642.67
    bin = 1
    for i, row in data.iterrows():
        current_bed = data.iloc[intervals[bin], :][['chrom', 'chromStart', 'chromEnd', 'type', 'size']].copy()
        if not current_bed.empty:
            while current_bed['size'].sum() + row['size'] >= mean_length:
                if bin == bins:
                    bin = 1
                else:
                    bin += 1
                current_bed = data.iloc[intervals[bin], :][['chrom', 'chromStart', 'chromEnd', 'type', 'size']].copy()
        intervals, bin = add_bin(intervals, i, bin)
    
    for bin in range(1, bins + 1):
        out = 'intervals/GRCh38_interval_' + str(bins) +  '_bin_' + str(bin) + '.bed'
        bed = data.iloc[intervals[bin], :][['chrom', 'chromStart', 'chromEnd', 'type', 'size']].copy()
        print(bin, bed['size'].sum())
        bed['chrom'] = pd.Categorical(bed['chrom'], 
                                      categories=chroms,
                                      ordered=True)
        bed = bed.sort_values(['chrom', 'chromStart'])
        bed.to_csv(out, header=False, index=False, sep='\t')
        
    
    