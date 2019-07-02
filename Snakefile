shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
FILES = json.load(open(config['SAMPLES_JSON']))
localrules: all, remove_backlist


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
EXTEND_SUMMIT = []

for sample in SAMPLES:
    CLUSTERS_BAIS.append(expand("01split_bam/{sample}/{sample}_{{cluster_id}}.bam.bai". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    CLUSTERS_BIGWIGS.extend(expand("02bigwigs/{sample}/{sample}_{{cluster_id}}.bw". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    PEAKS.extend(expand("05peak_filter/{sample}/{sample}_{{cluster_id}}_blacklist_removed.bed". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))
    EXTEND_SUMMIT.extend(expand("06extend_summit/{sample}/{sample}_{{cluster_id}}_extend_summit.bed". \
        format(sample = sample), cluster_id = get_cluster_id(FILES[sample][1])))

TARGET.extend(CLUSTERS_BAMS)
TARGET.extend(CLUSTERS_BAIS)
TARGET.extend(CLUSTERS_BIGWIGS)
TARGET.extend(PEAKS)
TARGET.extend(EXTEND_SUMMIT)

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

## see https://twitter.com/lh3lh3/status/1129900399626981382
## use samtools collate not samtools sort -n 
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
        -o {output} \
        01split_bam/{wildcards.sample}/{wildcards.sample}_{wildcards.cluster_id}.bam 2> {log}
        """


rule call_peaks:
    input: "03name_sorted_bam/{sample}/{sample}_{cluster_id}.name.sorted.bam"
    output: "04peak/{sample}/{sample}_{cluster_id}_Genrich.bed"
    log: "00log/{sample}_{cluster_id}_Genrich.log"
    threads: 1
    params:
        custom = config.get("Generich_args", "")
    shell:
       """
       Genrich -t {input} -o {output} {params.custom} 2> {log}
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

rule extend_summit:
    input: "05peak_filter/{sample}/{sample}_{cluster_id}_blacklist_removed.bed"
    output: summit = "06extend_summit/{sample}/{sample}_{cluster_id}_summit.bed",
            extend_summit = "06extend_summit/{sample}/{sample}_{cluster_id}_extend_summit.bed"
    log: "00log/{sample}_{cluster_id}_extend_summit.log"
    threads: 1
    shell:
        """
        # the 10th column is the peak summit
        cat {input} | awk -v OFS="\t" '$2=$2+$10,$3=$2+1' > {output.summit}
        bedtools slop -i {output.summit} -g {config[genome_size]}  -b {config[extend_size]} > {output.extend_summit}
        """

######################################################################

########   merge all peaks and recount reads in peaks per cell #######

######################################################################

def get_extend_summit(wildcards):
    sample = wildcards.sample
    cluster_id = get_cluster_id(FILES[wildcards.sample][1])
    return expand("06extend_summit/{sample}/{sample}_{cluster_id}_extend_summit.bed", \
        sample = sample, cluster_id = cluster_id)


rule merge_extend_summit:
    input: get_extend_summit
    output: "07recount/{sample}/{sample}_merged_peaks.bed"
    shell:
        """
        
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}

        """

rule recount:
    input: bed = "07recount/{sample}/{sample}_merged_peaks.bed",
           bam = lambda wildcards: FILES[wildcards.sample][0]
    output: "07recount/"
    log: "00log/recount.log"
    threads: 1
    shell:
        """
        scripts/rcbbc_v02 -w {config["white_list"]} -p {wildcards.sample} {input.bed} {input.bam} 
        """





######################################################################

########         motif footprint using rgt_hint                #######

######################################################################
