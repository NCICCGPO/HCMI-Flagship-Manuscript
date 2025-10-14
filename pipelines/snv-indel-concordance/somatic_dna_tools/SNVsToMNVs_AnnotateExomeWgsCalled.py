#!/usr/bin/env python

from Classes import *

def assemble_header(header, header_keys=["fileformat","FILTER","FORMAT","INFO","contig","cmdline","col_headers"]):
    header_text=""
    for k in header_keys:
        if k in header:
            if k in ["FILTER","FORMAT","INFO"]:
                header_text=header_text+"\n".join(sorted(set(header[k])))+"\n"
            else:
                header_text=header_text+"\n".join(header[k])+"\n"
    return header_text

def __main__():
    parser = ArgumentParser(prog=os.path.basename(sys.argv[0]), description='Parses file and converts adjacent SNVs to MNVs if they have they match the MNV_ID and called_by fields.', epilog='', formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=100, width=150))
    parser.add_argument('-i', '--input_vcf', help = 'Input VCF file.', required=True)
    parser.add_argument('-o', '--output_vcf', help = 'Output VCF file.', required=True)
    args=parser.parse_args()
    INPUT=args.input_vcf
    OUTPUT=args.output_vcf
    if not os.path.isfile(INPUT):
        print("ERROR: Required file {0} does not exist. Cannot run.".format(INPUT))
        sys.exit(1)
    f=open(INPUT)
    content=f.readlines()
    o=open(OUTPUT,"w")
    info_pattern=re.compile('(\w+)\=([^;]+)')
    MNV=dict()
    SNV=dict()
    mnv_filter=dict()
    header=dict()

    ## Parse input VCF and identify SNVs to merge
    ## Expects sorted VCF
    
    seen_filter_header=0
    seen_info_header=0
    for line in content:
        line=line.strip()
        # parse metadata
        if line.startswith("##"):
            meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
            m=meta_pattern.match(line)
            keyid=m.group('key')
            if keyid not in header:
                header[keyid]=[]
            header[keyid].append(line)
        # parse chrom header
        elif line.startswith("#CHROM"):
            header['col_headers']=[line]
        # parse records
        else:
            (CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,NORMAL,TUMOR)=line.strip().split("\t")
            info_pattern=re.compile('(\w+)\=([^;]+)')
            info_dict=dict(info_pattern.findall(INFO))
            normal_format_dict=dict(zip(FORMAT.split(":"),NORMAL.split(":")))
            tumor_format_dict=dict(zip(FORMAT.split(":"),TUMOR.split(":")))
            varID=":".join([CHROM,POS,REF,ALT])

            ### check if variant is an MNV or part of MNV ##
            if "MNV_ID" in info_dict \
                    and "TYPE" in info_dict:
                #varID=":".join([CHROM,POS,REF,ALT])
                # type may not be MNV if an SNV also existed in the merge
                if info_dict["TYPE"]=="MNV":
                    MNV[varID]=info_dict
                    MNV[varID]["SNVs"]=[]
                    for i in range(len(REF)):
                        ref=REF[i]
                        alt=ALT[i]
                        if not ref == alt:
                            pos=int(POS)+i
                            MNV[varID]["SNVs"].append(":".join([CHROM,str(pos),ref,alt]))
                elif info_dict["TYPE"]=="SNV":
                    SNV[varID]=info_dict
    # for every varID that goes to type MNV
    for varID in MNV.keys():
        called_by=[] # list callers for SNV supporting this MNV
        MNV_ID=[] # list MNV_IDs for this group?
#        keep_mnv=1
        # for every varID that goes to a SNV supporting this MNV
        for snvID in MNV[varID]["SNVs"]:
            if "called_by" in SNV[snvID]:
                called_by.append(SNV[snvID]["called_by"])
            MNV_ID.append(SNV[snvID]["MNV_ID"])
        # Split into SNVs if different callers are calling some of the site
        # note each *ed_by is a string which can become a comma separated list
        # order is conserved in the list
        if len(set(called_by)) > 1 \
                or len(set(MNV_ID)) > 1 \
                or set(MNV_ID) != set([MNV[varID]["MNV_ID"]]):
            mnv_filter[varID]="SplitToSNVs"
        else:
            # keep MNV together because no lone evidence supporting
            # any one individual SNVs was reported
            for snvID in MNV[varID]["SNVs"]:
                mnv_filter[snvID]="PartOfMNV"
            # if the MNV has a called_by update it with support from anything else
            if "called_by" in MNV[varID]:
                snv_callers = called_by[0].split(",")
                mnv_callers=MNV[varID]["called_by"].split(",") # check who called the MNV
                if len(set(mnv_callers)) < len(set(snv_callers)): # check for new callers
                    additional_snv_callers = list(set(mnv_callers)^set(snv_callers))

    header["FILTER"].append('##FILTER=<ID=PartOfMNV,Description="The candidate SNV is part of a MNV">')
    header["FILTER"].append('##FILTER=<ID=SplitToSNVs,Description="The candidate MNV is now split into multiple SNVs">')
    header["INFO"].append('##INFO=<ID=WgsExomeCalled,Number=0,Type=Flag,Description="High Confidence variant">')
    o.write(assemble_header(header))

    ## Read VCF for the second time ##
    for line in content:
        if line.startswith("#"):
            continue
        else:
            (CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,NORMAL,TUMOR) = line.strip().split("\t")
            info_pattern=re.compile('(\w+)\=([^;]+)')
            info_dict=dict(info_pattern.findall(INFO))
            varID=":".join([CHROM,POS,REF,ALT])
            if varID in mnv_filter:
                if FILTER != "PASS":
                    FILTER=FILTER+";"+mnv_filter[varID]
                else:
                    FILTER=mnv_filter[varID]
            info_text_list=[]
            for key in sorted(info_dict.keys()):
                info_text_list.append('{0}={1}'.format(key,info_dict[key]))
            INFO=";".join(info_text_list)
            if "num_callers" in info_dict \
                    and FILTER=="PASS":
                num_callers=int(info_dict["num_callers"])
                if num_callers>1:
                    INFO=INFO+";WgsExomeCalled"
            o.write("\t".join([CHROM,POS,ID,REF,
                               ALT,QUAL,FILTER,
                               INFO,FORMAT,NORMAL,
                               TUMOR])+"\n")
    o.close()
    f.close()

if __name__ == "__main__":
    __main__()
