import numpy as np
import pysam
import pandas as pd
import itertools


class CGC():
    '''
        Get GCG annotation from VCF or from lookup table made from orthologs
        to human CGC genes
    '''
    def __init__(self, bcf_in=False,
                 csq_columns=False,
                 cosmic_census_file=False,
                 vcf_backup=False):
        '''
            Set sample-wide values.
        '''
        self.vcf_backup = vcf_backup
        if 'CancerGeneCensus' in bcf_in.header.info:
            source='CancerGeneCensus'
            cgc_columns = self.bcf_in.header.info[source].description.split('Format: ')[-1].lstrip().split('|') # grab the definitions
        else:
            cgc_columns = False
        if cosmic_census_file:
            self.cosmic_census = pd.read_csv(cosmic_census_file, dtype={'Tier': str})
        self.cgc_columns = cgc_columns
        self.csq_columns = csq_columns


    def set_annotation(self, i, record):
        '''
            Get annotations from VCF INFO for human.
            Get CGC annotation from lookup in CGC ortholog list for
            mouse.
        '''
        if self.cgc_columns:
            self.cgc_dicts = self.get_csqs(record, 
                                           columns=self.cgc_columns,
                                           source='CancerGeneCensus')
            self.tier = self.cgc_dicts[i]['Tier']
        else:
            csq_dicts = self.get_csqs(record, 
                                      columns=self.csq_columns,
                                      source='CSQ')
            gene = csq_dicts[i]['SYMBOL']
            self.tier = ''
            if gene in self.cosmic_census['Gene Symbol'].values:
                self.tier = self.cosmic_census.loc[(self.cosmic_census['Gene Symbol'] == gene),
                                                   'Tier'].values[0]

    def set_cgc_dicts(self, record):
        if self.cgc_columns:
            self.cgc_dicts = self.get_csqs(record, 
                                           columns=self.cgc_columns, 
                                           source='CancerGeneCensus')
        else:
            self.cgc_dicts = {i : {'Tier' : ''} for i, alt in enumerate(record.alts)}
            for i, alt in enumerate(record.alts):
                self.set_annotation(i, record)
                self.cgc_dicts[i]['Tier'] = self.tier

    def lookup_role_in_cancer(self, gene):
        '''
            Check if gene is listed as a oncogene or supressor (TSG)
        '''
        self.tsg = 0
        self.oncogene = 0
        self.fusion = 0
        if gene in self.cosmic_census['Gene Symbol'].values:
            self.tier = self.cosmic_census.loc[(self.cosmic_census['Gene Symbol'] == gene),
                                               'Tier'].values[0]
            role = self.cosmic_census.loc[(self.cosmic_census['Gene Symbol'] == gene),
                                               'Role in Cancer'].values[0]
            try:
                roles = [role.strip() for role in role.split(',')]
            except AttributeError:
                roles = []
            if 'TSG' in roles:
                self.tsg = 1
            if 'fusion' in roles:
                self.fusion = 1
            if 'oncogene' in roles:
                self.oncogene = 1


class NYGC(CGC):
    def __init__(self, vcf,
                 source='CSQ',
                 cosmic_census_file=False):
        '''
            Load VCF and backup
            cosmic_census_file provided only for mouse
        '''
        self.cosmic_census_file = cosmic_census_file
        # backout runs when VEP annotation include special characters that break PYSAM
        self.vcf_backup = self.load_utf_backup(vcf)
        self.bcf_in = self.read_vcf(vcf)
        self.set_csq_columns()
        self.source = 'CSQ'
        super().__init__(bcf_in=self.bcf_in,
                         csq_columns=self.csq_columns,
                         cosmic_census_file=self.cosmic_census_file,
                         vcf_backup=self.vcf_backup)
        
    def read_vcf(self, vcf_file):
        '''
            Read in annotated VCF file.
            '''
        bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
        return bcf_in

    
    def load_utf_backup(self, vcf_file, germ=False):
        '''
            Load searchable VCF for lines that are not ASCII.
        '''
        last = self.find_header(vcf_file)
        if germ:
            header = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'INFO', 'format', 'normal']
        else:
            header = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'INFO', 'format', 'normal', 'tumor']
        try:
            data = pd.read_csv(vcf_file, skiprows=last, sep='\t',
                               names=header, encoding = 'utf8', dtype = {'chrom' : str})
        except UnicodeDecodeError:
            data = pd.read_csv(vcf_file, skiprows=last, sep='\t',
                               names=header, encoding = 'latin-1', dtype = {'chrom' : str})
        return data
    
    def find_header(self, vcf_file):
        '''
            Get VCF header line numbers (because pandas can't skip
            based on > 1 character words. ## is a comment in VCF but
            # can occur in the VCF INFO fields.
            '''
        try:
            with open(vcf_file) as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        except UnicodeDecodeError:
            with open(vcf_file, encoding='latin-1') as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        return last + 1

    def set_csq_columns(self, source='CSQ'):
        '''
            get column names from the bar
            separated CSQ VEP annotation
            results. CSQ are Consequence
            annotations from Ensembl VEP.
            '''
        if source == 'CSQ':
            self.csq_columns = self.bcf_in.header.info[source].description.split()[-1].split('|') # grab the definitions
        else:
            self.csq_columns = self.bcf_in.header.info[source].description.split('Format: ')[-1].lstrip().split('|') # grab the definitions

    def get_csqs(self, record, columns=False, source=False):
        '''
            Get new INFO field results. Works with (python3). Skip lines starting with '#' for vcf backup.
            Because # can occur in info lines from VEP it can't be used as a comment character
        '''
        if not columns:
            columns=self.csq_columns
        if not source:
            source=self.source
        alt_count = len(record.alts)
        csq_values = []
        csq_dicts = {}
        for i in range(alt_count):
            try:
                csq_line = record.info[source][i]
            except UnicodeDecodeError: # for CSQ results with accents and other unexpected non-ascii characters (rare)
                line_data = self.vcf_backup[(self.vcf_backup.chrom == record.chrom)
                                            & (self.vcf_backup.pos == record.pos)
                                            & (self.vcf_backup.ref == record.ref)
                                            & (self.vcf_backup.alt == ','.join(record.alts))]
                csq_line = line_data.INFO.values.tolist()[0].split(source + '=')[1]
                csq_line = csq_line.split(';')[0]
            csq_line = csq_line.split(',')[i]
            csq_values = csq_line.split('|')
            try:
                # Works with (python2)
                csq_dict = dict(itertools.izip(columns, csq_values))
            except AttributeError:
                # Works with (python3)
                csq_dict = dict(zip(columns, csq_values))
            csq_dicts[i] = csq_dict
        return csq_dicts
