shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"
FILES = json.load(open(config['SAMPLES_JSON']))
localrules: all, remove_backlist, remove_backlist_merge, extend_summit, extend_summit_merge


## avoid dynamic rules, read the cluster id file and determine the number of splitted bam
## cluster_id.csv should have no headers
import csv
from collections import defaultdict 
SAMPLES = list(FILES.keys())

## for samples are not clustered together, different sample may have different number of clusters
def get_cluster_id(csv_file):
    cluster_ids = []
    with open(csv_file) as ifile:
        csv_reader = csv.reader(ifile, delimiter=',')
        # if there is a header, skip it by
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

### for merged bam files
map_sample_to_cluster = defaultdict(list)

for sample in SAMPLES:
    map_sample_to_cluster[sample].extend(get_cluster_id(FILES[sample][1]))

CLUSTERS = []
for sample in map_sample_to_cluster:
    CLUSTERS.extend(map_sample_to_cluster[sample])
CLUSTERS = sorted(set(CLUSTERS))


## some clusters only present in certain samples, when merge the same cluster from different 
## samples, this needs to be taken care of.

## e.g.
## {'A':['1','2','3'], 'B':['1','2'], 'C':['1']}
## {'1':['A', 'B', 'C'], '2':['A','B'], '3':[1]}

## for cluster 1, merge bam from A B C samples
## for cluster 2 merge bam from A B samples
## for cluster 3 merge bam only from A sample.

map_cluster_to_sample = defaultdict(list)
for k,v in map_sample_to_cluster.items():
    for x in v:
        map_cluster_to_sample[x].append(k)

 

MERGED_BIGWIGS = expand("02bigwigs_merge/{cluster_id}.bw", cluster_id = CLUSTERS)
MERGED_EXTEND_SUMMIT = expand("06extend_summit_merge/{cluster_id}_extend_summit.bed", cluster_id = CLUSTERS )
RECOUNT_ALL = expand("07recount_all/{sample}/{sample}.mtx", sample = SAMPLES)
RECOUNT = expand("07recount/{sample}/{sample}.mtx", sample = SAMPLES)
DIFFPEAKS = expand("08diff_peaks/{sample}/{sample}_differential_accessible_peaks.txt", sample =SAMPLES)


TARGET.extend(CLUSTERS_BAMS)
TARGET.extend(CLUSTERS_BAIS)
TARGET.extend(CLUSTERS_BIGWIGS)
TARGET.extend(PEAKS)
TARGET.extend(EXTEND_SUMMIT)
TARGET.extend(MERGED_BIGWIGS)
TARGET.extend(MERGED_EXTEND_SUMMIT)
TARGET.extend(DIFFPEAKS)

if config['recount_all']:
    TARGET.extend(RECOUNT_ALL)
    TARGET.append("08diff_peaks_all/all_sample_differential_accessible_peaks.txt")
TARGET.extend(RECOUNT)

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

###  merge bam files for the same cluster from different samples ####
### different samples are merged together for cluster identification ##

######################################################################

## note 01split_bam/{sample}/{sample}_{{cluster_id}}.bam is not any output from 
## previous rules because I used the touch file to split the bam.
## however, the bam.bai index is produced by previous rules, so I use it as input
## and remove the .bai in params directive for shell to make it work.


def get_merge_bam_input(wildcards):
    samples = map_cluster_to_sample[wildcards.cluster_id]
    return expand("01split_bam/{sample}/{sample}_{{cluster_id}}.bam.bai", sample = samples)

rule merge_bam_per_cluster:
    input: get_merge_bam_input
    output: "01merged_bam/{cluster_id}.bam"
    params:
            bam = lambda wildcards, input: " ".join(input).replace(".bai", "")
    log: "00log/{cluster_id}.merge_bam.log"
    threads: 2
    shell:
        """
        samtools merge -@ 2 -r {output} {params.bam}
        """

rule index_merged_bam:
    input: "01merged_bam/{cluster_id}.bam"
    output: "01merged_bam/{cluster_id}.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule make_merged_bigwig:
    input: "01merged_bam/{cluster_id}.bam", "01merged_bam/{cluster_id}.bam.bai"
    output: "02bigwigs_merge/{cluster_id}.bw"
    log: "00log/{cluster_id}_make_merge_bigwig.log"
    threads: 5
    params:
        custom= config.get("bamCoverage_args", "")
    message: "making bigwig files for {input} with {threads} threads"
    shell:
        """
        bamCoverage -b {input[0]} -p {threads} {params.custom} -o {output}
        """


rule sort_merge_bam_by_name:
    input: "01merged_bam/{cluster_id}.bam", "01merged_bam/{cluster_id}.bam.bai"
    output:"03name_sorted_bam_merge/{cluster_id}.name.sorted.bam"
    log: "00log/{cluster_id}_name_sort_merge.log"
    threads: 5
    params:
        custom = config.get("name_sort_bam_agrs", "")
    message: "sorting {input} by name"
    shell:
        """
        samtools sort -n -m 2G -@ {threads} -T {wildcards.cluster_id} \
        -o {output} \
        {input[0]} 2> {log}
        """

rule call_peaks_merge:
    input: "03name_sorted_bam_merge/{cluster_id}.name.sorted.bam"
    output: "04peak_merge/{cluster_id}_Genrich.bed"
    log: "00log/{cluster_id}_Genrich.log"
    threads: 1
    params:
        custom = config.get("Generich_args", "")
    shell:
       """
       Genrich -t {input} -o {output} {params.custom} 2> {log}
       """

### remove black listed regions
rule remove_backlist_merge:
    input: "04peak_merge/{cluster_id}_Genrich.bed", config["black_list"]
    output: "05peak_filter_merge/{cluster_id}_blacklist_removed.bed"
    log: "00log/{cluster_id}_remove_blacklist.log"
    threads: 1
    shell:
        """
        bedtools intersect -a {input[0]} -b <(zcat {input[1]}) -v > {output}

        """

rule extend_summit_merge:
    input: "05peak_filter_merge/{cluster_id}_blacklist_removed.bed"
    output: summit = "06extend_summit_merge/{cluster_id}_summit.bed",
            extend_summit = "06extend_summit_merge/{cluster_id}_extend_summit.bed"
    log: "00log/{cluster_id}_extend_summit.log"
    threads: 1
    shell:
        """
        # the 10th column is the peak summit
        cat {input} | awk -v OFS="\t" '$2=$2+$10,$3=$2+1' > {output.summit}
        bedtools slop -i {output.summit} -g {config[genome_size]}  -b {config[extend_size]} > {output.extend_summit}
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
    message: "calling peaks for {input}"
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

##   merge all peaks per sample and recount reads in peaks per cell #######

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
        # get rid of unconventional chromosomes and natural sort the chr (chr10 after chr9)
        cat {input} | sort -k1,1V -k2,2n | bedtools merge | grep -v "_"  > {output}

        """

rule recount:
    input: bed = "07recount/{sample}/{sample}_merged_peaks.bed",
           bam = lambda wildcards: FILES[wildcards.sample][0]
    output: "07recount/{sample}/{sample}.mtx"
    log: "00log/recount_{sample}.log"
    params: white_list = lambda wildcards: FILES[wildcards.sample][2]
    threads: 1
    message: "recouting reads per peak for {input}"
    shell:
        """
        scripts/rcbbc -w {params.white_list} -p 07recount/{wildcards.sample}/{wildcards.sample} {input.bed} {input.bam} 
        """

######################################################################

########  merge all peaks for all samples and recount reads in peaks per cell #######

######################################################################

rule merge_extend_summit_all:
    input: expand("06extend_summit_merge/{cluster_id}_extend_summit.bed", cluster_id = CLUSTERS)
    output: "07recount_all/merged_peaks.bed"
    threads: 1
    shell:
        """
        # remove unconventional chromosomes and natural sort by -V
        cat {input} | sort -k1,1V -k2,2n | bedtools merge | grep -v "_"  > {output}
        """


rule recount_all:
    input: bed = "07recount_all/merged_peaks.bed", 
           bam = lambda wildcards: FILES[wildcards.sample][0]
    output: "07recount_all/{sample}/{sample}.mtx"
    log: "00log/recount_all_{sample}.log"
    params: white_list = lambda wildcards: FILES[wildcards.sample][2]
    threads: 1
    message: "recouting reads per peak for {input}"
    shell:
        """
        scripts/rcbbc -w {params.white_list} -p 07recount_all/{wildcards.sample}/{wildcards.sample} {input.bed} {input.bam} 
        """

######################################################################

########         Find differential peaks using presto            #####

######################################################################

rule presto:
    input: mtx = "07recount/{sample}/{sample}.mtx",
           csv = lambda wildcards: FILES[wildcards.sample][1]
    output: "08diff_peaks/{sample}/{sample}_differential_accessible_peaks.txt"
    singularity: "docker://crazyhottommy/seuratv3_presto"
    log: "00log/diff_peaks_{sample}.log"
    message: "Finding differential peaks using {input} "
    shell:
        """
        Rscript scripts/presto.R --mtx {input.mtx} --cluster {input.csv} {output}
        """


rule presto_all:
    input: mtxs = expand("07recount_all/{sample}/{sample}.mtx", sample = SAMPLES),
           csvs = [FILES[sample][1] for sample in SAMPLES]
    output: "08diff_peaks_all/all_sample_differential_accessible_peaks.txt"
    singularity: "docker://crazyhottommy/seuratv3_presto"
    log: "00log/diff_peaks_all_sample.log"
    message: "Finding differential peaks using {input} "
    shell:
        """
        mtxs=""
        for mtx in {input.mtxs}; do
            mtxs+=" --mtx=$mtx"
        done

        csvs=""
        for csv in {input.csvs}; do
            csvs+=" --cluster=$csv"
        done

        Rscript scripts/presto.R $mtxs $csvs {output}> {log} 2>&1 
        """

######################################################################

########         motif footprint using rgt_hint                #######

######################################################################
