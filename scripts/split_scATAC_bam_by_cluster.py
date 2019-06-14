#! /usr/bin/env python3

import pysam
import csv
import argparse
import os.path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("csv", help="Required. the FULL path to the cluster csv file with header, \
    first column is the cell barcode, second column is the cluster id")
parser.add_argument("bam", help="Required. the FULL path to the 10x scATAC bam file generated \
    by cellranger-atac count")
parser.add_argument("-prefix", help="Optional, the prefix of the output bam, default is cluster_id.bam")
parser.add_argument("-outdir", help="Optional, the output directory for the splitted bams, default is current dir")
args = parser.parse_args()


if os.path.exists(args.csv):
    pass
else:
    print("csv file is not found")
    sys.exit(1)

if os.path.exists(args.bam):
    pass
else:
    print("10x scATAC bam not found")
    sys.exit(1)

if args.outdir:
    if os.path.isdir(args.outdir):
        pass
    else:
        try:
            os.mkdir(args.outdir)
        except OSError:
            print("can not create directory {}".format(args.outdir))

cluster_dict = {}
with open(args.csv) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())


fin = pysam.AlignmentFile(args.bam, "rb")

# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb
fouts_dict = {}
for cluster in clusters:
    if args.prefix:
        fout_name = args.prefix + "_cluster_" + cluster + ".bam"
    else:
        fout_name = "cluster_" + cluster + ".bam"
    if args.outdir:
        fout = pysam.AlignmentFile(os.path.join(args.outdir,fout_name), "wb", template = fin)
    else:
        fout = pysam.AlignmentFile(fout_name, "wb", template = fin)
    fouts_dict[cluster] = fout

for read in fin:
    tags = read.tags
    # the 8th item is the CB tag
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    else: 
        continue
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)

## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
