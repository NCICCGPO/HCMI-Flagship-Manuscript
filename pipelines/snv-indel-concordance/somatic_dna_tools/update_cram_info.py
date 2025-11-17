
import argparse
import json


# takes a bam info files and array of cram sample info and reformats them
# output is normalSampleCramInfos and pairCramInfos objects

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair-infos',
                        help='pair info json file.',
                        required=True
                        )
    parser.add_argument('--normal-infos',
                        help='normal bam info json file.',
                        required=True
                        )
    parser.add_argument('--cram-infos',
                        help='cram info json file.',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    with open(args['pair_infos'], 'r') as fh:
        pair_infos = json.load(fh)
    with open(args['normal_infos'], 'r') as fh:
        normal_infos = json.load(fh)
    with open(args['cram_infos'], 'r') as fh:
        cram_infos = json.load(fh)

    # flatten cram info from list to dict with sampleId keys
    # keys not necessarily unique (ie. single sampleId with 2 different bams) but in practice should be
    cram_dict = {}
    for c in cram_infos:
        cram_dict[c['sampleId']] = c['finalCram']

    # make the pair info file
    pair_cram_info = []
    for p in pair_infos:
        try:
            res = {
                "tumor": p['tumorId'],
                "normal": p['normalId'],
                "pairId": p['pairId'],
                "normalFinalCram": cram_dict[p['normalId']],
                "tumorFinalCram": cram_dict[p['tumorId']],
            }
        except KeyError:
            res = p
    pair_cram_info.append(res)

    # make the normal cram info file
    normal_cram_info = []
    for n in normal_infos:
        try:
            cram = {
                "sampleId": n['sampleId'],
                "finalCram": cram_dict[n['sampleId']]
            }
        except KeyError:
            cram = n
        normal_cram_info.append(cram)

    with open('pair_cram_infos.json', 'w') as fh:
        json.dump(pair_cram_info, fh)

    with open('normal_cram_infos.json', 'w') as fh:
        json.dump(normal_cram_info, fh)


if __name__ == "__main__":
    main()

