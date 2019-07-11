# pyflow-scATACseq
snakemake workflow for post-processing scATACseq data

### What does it do?

for single cell ATACseq experiment, one get a merged bam file for all cells. After going through clustering, one group similar cells into cell types (cell states). This workflow will split the bam by clusters to create a pseudo bulk bam for each cluster, create bigwig tracks for visulization, call peaks for each cluster and merge the peaks across the clusters. Finally it will count reads per peak per cell from the original bam file on the merged peaks.

Thanks [Brent Pedersen](https://github.com/brentp) for the `rcbbc` (Read Count By BarCode) in the scripts folder.

### workflow of the piepline

![](./rulegraph.png)


### Dependencies

### How to use this pipeline

1. clone the repo 

```bash
git clone https://github.com/crazyhottommy/pyflow-scATACseq
cd pyflow-scATACseq
```

2. Prepare a tab delimited `tsv` file with header: `meta.tsv` inside the `pyflow-scATACseq` folder:

| sample  | bam_path         | cluster_id.csv               | white_list                      | 
|---------|------------------|------------------------------|---------------------------------| 
| sample1 | /path/to/sample1 | /path/to/cluster_sample1.csv | /path/to/white_list_sample1.txt | 
| sample2 | /path/to/sample2 | /path/to/cluster_sample2.csv | /path/to/white_list_sample2.txt | 


  * The first column is the sample name.  
  * The second column is the path to the `10x cellranger_atac` produced `possorted_bam.bam`  
  * The `cluster_id.csv` is a two column csv file **with** header (does not matter what name you give). The first column is the cell barcode, the second column is the cluster id (can be numeric or strings).

```
cat example_cluster.csv
name,value
AAACGAAAGCACCATT-1,32
AAACGAAAGTAGTCGG-1,4
AAACGAAGTAACGGCA-1,31
AAACTCGAGAACAGGA-1,17
AAACTCGAGTCCAGAG-1,1
AAACTCGCAAAGGAAG-1,41
AAACTCGCAAGGCGTA-1,5
AAACTCGCAGAAAGAG-1,4
AAACTCGCAGTAGGCA-1,24
```
  * The white_list is a one column txt file with valid barcode you want to include.

you can create it by:

```bash
cat example_cluster.csv | cut -f, -d1 | sed '1d' > example_white_list.txt
```

2. Create a `samples.json` file:

```bash
python sample2json.py meta.tsv

# a samples.json file will be created
cat samples.json
{
    "sample1": [
        "/path/to/sample1",
        "/path/to/cluster_sample1.csv",
        "/path/to/white_list_sample1.txt"
    ],
    "sample2": [
        "/path/to/sample2",
        "/path/to/cluster_sample2.csv",
        "/path/to/white_list_sample2.txt"
    ]
}
```

3. Running the workflow

```bash
#dry run
snakemake -np 

## real-run on local machine
snakemake -j 10

## submit to slurm (you can change this script for your own HPC)
./pyflow-scATACseq.sh 
```


### Other Notes 

more thoughts on pileup and shifting reads
https://github.com/deeptools/deepTools/issues/453
maxFragmentLength

https://github.com/deeptools/deepTools/issues/370

https://groups.google.com/forum/#!topic/deeptools/JU9itiT5rYk

You probably want “–Offset 1”

The above is not shifting exactly, for shifting use
https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html

--ATACshift

filter out fragments?