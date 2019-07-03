#! /usr/bin/env python3

import json
import csv
import argparse
from collections import defaultdict	


parser = argparse.ArgumentParser()
parser.add_argument("meta", help="Required. the FULL path to the tab delimited meta file")
args = parser.parse_args()

assert args.meta is not None, "please provide the path to the meta file"

FILES = defaultdict(list)

with open(args.meta) as ifile:
    csv_reader = csv.reader(ifile, delimiter='\t')
    header = next(csv_reader)
    for row in csv_reader:
        sample = row[0]
        bam_path = row[1]
        cluster = row[2]
        white_list = row[3]
        FILES[sample].append(bam_path)
        FILES[sample].append(cluster)
        FILES[sample].append(white_list)

sample_num = len(FILES.keys())
print ("total {} unique samples will be processed".format(sample_num))
print ("------------------------------------------")
print("check the samples.json file for metadata for each sample")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)

