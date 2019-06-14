shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

localrules: all


######################################################################

######## split bam file by user provided cluster id ##################

######################################################################


## avoid dynamic rules, read the cluster id file and determine the number of splitted bam
## cluster_id.csv should have no headers
import csv
cluster_ids = []
with open(config["cluster_id.csv"]) as ifile:
	csv_reader = csv.reader(ifile, delimiter=',')
	# if there is a header, skip it by
	# header = next(csv_reader)
	#skip header
	header = next(csv_reader)
	for row in csv_reader:
		cluster_ids.append(row[1])

cluster_ids = sorted(set(cluster_ids))


TARGET= []
CLUSTERS_BAMS = expand("01split_bam/{prefix}_cluster_{{cluster_id}}.bam".format(prefix = config["prefix"]), cluster_id = cluster_ids)
CLUSTERS_BIGWIGS = expand("02bigwigs/{prefix}_cluster_{{cluster_id}}.bw".format(prefix = config["prefix"]), cluster_id = cluster_ids)

TARGET.extend(CLUSTERS_BAMS)
TARGET.extend(CLUSTERS_BIGWIGS)

rule all:
	input: TARGET

rule split_scATAC_bam:
	input: config["merged_scATACseq_bam"]
	output:  CLUSTERS_BAMS
	threads: 1
	shell:
		"""
		scripts/split_scATAC_bam_by_cluster.py -prefix {config[prefix]} -outdir 01split_bam {input} {config[cluster_id.csv]}
		"""

rule index_bam:
	input: "01split_bam/{cluster}.bam"
	output: "01split_bam/{cluster}.bam.bai"
	shell:
		"""
		samtools index {input}
		"""


rule make_bigwig:
	input: "01split_bam/{cluster}.bam"
	output: "02bigwigs/{cluster}.bw"
	log: "00log/{cluster}_make_bigwig.log"
	threads: 5
	params:
		custom= config.get("bamCoverage_args", "")
	message: "making bigwig files for {input} with {threads} threads"
	shell:
		"""
		bamCoverage -b {input} -p {threads} {params.custom} -o {output}
		"""


#rule call_peaks:
#	input:
#	output:
#	shell:
#		"""
#		Generich 
#		"""


#rule rgt_hint_motif:
