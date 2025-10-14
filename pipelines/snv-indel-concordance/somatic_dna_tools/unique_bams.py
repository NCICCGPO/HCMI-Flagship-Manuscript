import argparse
import json


# takes a pairInfos file and returns list of UNIQUE bam objects with their sample_id

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair-infos',
                        help='pair info json file.',
                        required=True
                        )

    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    with open(args['pair_infos'], 'r') as fh:
        infos = json.load(fh)

    # loop through T/N pairs and keep unique bams
    unique_bams = []
    for i in infos:
        normal = {'sampleId': i['normalId'],
                  'finalBam': i['normalFinalBam']}
        # need to check if it is a duplicate before adding
        if normal not in unique_bams:
            unique_bams.append(normal)
        tumor = {'sampleId': i['tumorId'],
                 'finalBam': i['tumorFinalBam']}
        if tumor not in unique_bams:
            unique_bams.append(tumor)
    with open('unique_bams.json', 'w') as fh:
        json.dump(unique_bams, fh)


if __name__ == "__main__":
    main()
