shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
FILES = json.load(open(config['SAMPLES_JSON']))
localrules: all remove_backlist


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
CLUSTERS_BAMS = expand("01split_bam/{sample}.split.touch", sample = SAMPLES)
CLUSTERS_BAIS = []
CLUSTERS_BIGWIGS = []
PEAKS = []

for sample in SAMPLES:
    CLUSTERS_BAIS.append(expand("01split_bam/{sample}/{sample}_{{cluster_id}}.bam.bai". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    CLUSTERS_BIGWIGS.extend(expand("02bigwigs/{sample}/{sample}_{{cluster_id}}.bw". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    PEAKS.exend(expand("05peak_filter/{sample}/{sample}_{cluster_id}_blacklist_removed.bed", \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1]))

TARGET.extend(CLUSTERS_BAMS)
TARGET.extend(CLUSTERS_BAIS)
TARGET.extend(CLUSTERS_BIGWIGS)

rule all:
    input: TARGET


######################################################################

######## split bam file by user provided cluster id ##################

######################################################################

# see https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output

rule split_scATAC_bam_by_cluster:
    input: lambda wildcards: FILES[wildcards.sample]
    output: "01split_bam/{sample}.split.touch"
    threads: 1
    shell:
        """
        scripts/split_scATAC_bam_by_cluster.py -prefix {wildcards.sample} -outdir 01split_bam/{wildcards.sample} {input[1]} {input[0]}
        touch {output}
        """


rule index_bam:
    input: "01split_bam/{sample}.split.touch"
    output: "01split_bam/{sample}/{sample}_{cluster_id}.bam.bai"
    shell:
        """
        samtools index 01split_bam/{wildcards.sample}/{wildcards.sample}_{wildcards.cluster_id}.bam
        """

### think about how to normalize 

rule make_bigwig:
    input: "01split_bam/{sample}.split.touch", "01split_bam/{sample}/{sample}_{cluster_id}.bam.bai"
    output: "02bigwigs/{sample}/{sample}_{cluster_id}.bw"
    log: "00log/{sample}_{cluster_id}_make_bigwig.log"
    threads: 5
    params:
        custom= config.get("bamCoverage_args", "")
    message: "making bigwig files for {input} with {threads} threads"
    shell:
        """
        bamCoverage -b 01split_bam/{wildcards.sample}/{wildcards.sample}_{wildcards.cluster_id}.bam -p {threads} {params.custom} -o {output}
        """

######################################################################

######## call peaks for the splitted bam files using Genrich #########

######################################################################

rule sort_bam_by_name:
    input: "01split_bam/{sample}.split.touch", "01split_bam/{sample}/{sample}_{cluster_id}.bam.bai"
    output:"03name_sorted_bam/{sample}/{sample}_{cluster_id}.name.sorted.bam"
    log: "00log/{sample}_{cluster_id}_name_sort.log"
    threads: 5
    params:
        custom = config.get("name_sort_bam_agrs", "")
    message: "sorting {input} by name"
    shell:
        """
        samtools sort -n -m 2G -@ {threads} -T {wildcards.sample}_{wildcards.cluster_id} \
        -o {output} out-cluster-10.bam \
        01split_bam/{wildcards.sample}/{wildcards.sample}_{wildcards.cluster_id}.bam 2> {log}
        """


rule call_peaks:
    input: "03name_sorted_bam/{sample}/{sample}_{cluster_id}.name.sorted.bam"
    output: "04peak/{sample}/{sample}_{cluster_id}_Genrich.bed"
    log: log1 = "00log/{sample}_{cluster_id}_Genrich.log1",
        log2 ="00log/{sample}_{cluster_id}_Genrich.log2"
    threads: 1
    params:
        custom = config.get("Generich_args", "")
    shell:
       """
       Genrich -t {input} -o {output} {params.custom} -f {log1} 2> {log2}
       """

### remove black listed regions
rule remove_backlist:
    input: "04peak/{sample}/{sample}_{cluster_id}_Genrich.bed", config["black_list"]
    output: "05peak_filter/{sample}/{sample}_{cluster_id}_blacklist_removed.bed"
    log: "00log/{sample}_{cluster_id}_remove_blacklist.log"
    threads: 1
    shell:
        """
        bedtools intersect -a {input[0]} -b <(zcat {input[1]}) -v > {output}

        """
######################################################################

########         motif footprint using rgt_hint                #########

######################################################################

#rule rgt_hint_motif:
