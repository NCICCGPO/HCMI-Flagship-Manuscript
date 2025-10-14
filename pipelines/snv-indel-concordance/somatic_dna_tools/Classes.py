import argparse
import sys
import os
import re

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

class Variant:
    def __init__(self,chrom,pos,identity,ref,alt,qual,filter_value):
        self.chrom=chrom
        self.pos=pos
        self.identity=identity
        self.ref=ref
        self.alt=alt
        self.qual=qual
        self.filter=filter_value
        self.info_dict=dict()
        self.info_text="."
        self.tumor=dict()
        self.normal=dict()

    def print_variant(self):
        sorted_format_keys=sorted(self.tumor.keys())
        tumor_format=[]
        normal_format=[]
        for field in sorted_format_keys:
            if field in self.tumor:
                tumor_format.append(str(self.tumor[field]))
            else:
                tumor_format.append(".")
            if field in self.normal:
                normal_format.append(str(self.normal[field]))
            else:
                normal_format.append(".")
        self.tumor_format_text=":".join(tumor_format)
        self.normal_format_text=":".join(normal_format)
        self.format_keys_text = ":".join(sorted_format_keys)
        if len(self.info_dict.keys())>0:
            info_text_list=[]
            for key in sorted(self.info_dict.keys()):
                info_text_list.append('{0}={1}'.format(key,self.info_dict[key]))
            self.info_text=";".join(info_text_list)
        return "\t".join([self.chrom, self.pos, self.identity, self.ref, self.alt, self.qual, self.filter, self.info_text, self.format_keys_text, self.normal_format_text, self.tumor_format_text])

