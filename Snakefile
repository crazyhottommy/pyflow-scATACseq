shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
FILES = json.load(open(config['SAMPLES_JSON']))
localrules: all


######################################################################

######## split bam file by user provided cluster id ##################

######################################################################


## avoid dynamic rules, read the cluster id file and determine the number of splitted bam
## cluster_id.csv should have no headers
import csv
SAMPLES = FILES.keys()

def get_cluster_id(csv_file):
    cluster_ids = []
    with open(csv_file) as ifile:
        csv_reader = csv.reader(ifile, delimiter=',')
        # if there is a header, skip it by
        # header = next(csv_reader)
        #skip header
        header = next(csv_reader)
        for row in csv_reader:
            cluster_ids.append(row[1])
        cluster_ids = sorted(set(cluster_ids))
    return(cluster_ids)


TARGET= []
CLUSTERS_BAMS = []
CLUSTERS_BIGWIGS = []

for sample in SAMPLES:
    CLUSTERS_BAMS.append(expand("01split_bam/{sample}_{{cluster_id}}.bam". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    CLUSTERS_BAMS.append(expand("02bigwigs/{sample}_{{cluster_id}}.bw". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))

TARGET.extend(CLUSTERS_BAMS)
TARGET.extend(CLUSTERS_BIGWIGS)

rule all:
    input: TARGET


# see https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output

for sample in SAMPLES:
    rule:
        input: FILES[sample][0], FILES[sample][1]
        output: expand("01split_bam/{sample}_{{cluster_id}}.bam".format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1]))
        threads: 1
        shell:
            """
            scripts/split_scATAC_bam_by_cluster.py -prefix {sample} -outdir 01split_bam {input[0]} {input[1]}
            """


rule index_bam:
    input: "01split_bam/{sample}_{cluster_id}.bam"
    output: "01split_bam/{sample}_{cluster_id}.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule make_bigwig:
    input: "01split_bam/{sample}_{cluster_id}.bam", "01split_bam/{sample}_{cluster_id}.bam.bai"
    output: "02bigwigs/{sample}_{cluster_id}.bw"
    log: "00log/{sample}_{cluster_id}_make_bigwig.log"
    threads: 5
    params:
        custom= config.get("bamCoverage_args", "")
    message: "making bigwig files for {input} with {threads} threads"
    shell:
        """
        bamCoverage -b {input} -p {threads} {params.custom} -o {output}
        """


#rule call_peaks:
    
#   input:
#   output:
#   shell:
#       """
#       Generich 
#       """


#rule rgt_hint_motif:
