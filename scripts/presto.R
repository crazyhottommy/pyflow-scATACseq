#! /usr/bin/env Rscript
'Finding Differential accessible regions from one or more sparse matrix
The number of mtx sparse matrix should be the same as the cluster information csv
files and have the same basename. e.g. sample1.mtx and sample1.csv
Usage:
    presto.R --mtx=<mtx>... --cluster=<cluster>... <output>
    presto.R --mtx=sample1.mtx --cluster=sample1.csv sample1.txt
    presto.R --mtx=sample1.mtx --mtx=sample2.mtx --cluster=sample1.csv --cluster=sample2.csv sample.txt
    
Options:
    --mtx=<mtx>...  input file path of the sparse matrix
    --cluster=<cluster>...  input file path of the cluster information. a csv file with a header. column 1 is the
             cell barcode, column 2 is the cluster id.
    -h --help  Show this screen.
    -v --version  Show version.
Arguments:
    output   output file path of the differential peaks
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(presto))
suppressMessages(library(tidyverse))

arguments <- docopt(doc, version = 'presto.R v1.0\n\n')
num_of_matrix<- length(arguments$mtx)
num_of_csv<- length(arguments$cluster)

if (num_of_matrix != num_of_csv){
        stop("the number of matrix must be the same with the number
             of cluster information csv file") 
}


check_basename<- function(mtx, csv){
        return(gsub(".mtx", "", basename(mtx)) == gsub(".csv", "", basename(csv)))
}

if (!all(purrr::map2_lgl(arguments$mtx, arguments$cluster, check_basename))){
        stop("the basename of the mtx and the cluster info csv file do not match")
} 

read_mtx<- function(mtx){
        recount<- readMM(mtx)
        peaks<- read_tsv(paste0(mtx, ".sites"), col_names = FALSE)
        cells<- read_tsv(paste0(mtx, ".cells"), col_names = FALSE)
        prefix<- gsub(".mtx", "", basename(mtx))
        cells<- paste0(prefix, "_", cells$X1)
        rownames(recount)<- peaks$X1
        colnames(recount)<- cells
        return(recount)
}

message(paste("reading in the sparse matrix:", arguments$mtx, "\n"))
mtxs<- purrr::map(arguments$mtx, read_mtx)

if (!length(unique((purrr::map(mtxs, rownames))))==1){
        stop("all matrix should have same row/peak name")
}


recount.mat<- purrr::reduce(mtxs, cbind)

makeSeuratObject<- function(recount.mat){
        recount.TFIDF.mat<- TF.IDF(recount.mat)
        recount.TFIDF.mat <- Seurat:::LogNorm(data = recount.TFIDF.mat, scale_factor = 1e4)
        rownames(recount.TFIDF.mat)<- rownames(recount.mat)
        colnames(recount.TFIDF.mat)<- colnames(recount.mat)
        recount.atac<- CreateSeuratObject(counts = recount.TFIDF.mat, assay = 'ATAC', project = "recount")
        return(recount.atac)
}



### transferred label information
message(paste("reading in the cluster information csv files:", arguments$cluster, "\n"))
read_id<- function(csv){
        id<- read_csv(csv)
        prefix<- gsub(".csv", "", basename(csv))
        id<- id %>% 
                mutate(cell = paste0(prefix, "_", name)) %>%
                dplyr::rename(cluster_id = value)
        return(id)
}

cluster_ids<- purrr::map(arguments$cluster, read_id)
cluster_ids<- bind_rows(cluster_ids)

message("making Seurat ATAC objects")
recount.atac<- makeSeuratObject(recount.mat)
recount.atac[["cell"]]<- rownames(recount.atac@meta.data)
old_meta<- recount.atac@meta.data
new_meta<- left_join(recount.atac@meta.data, cluster_ids) 

if(!identical(rownames(old_meta), new_meta$cell)){
       stop("The name of the cells are different") 
}

rownames(new_meta)<- rownames(old_meta)
recount.atac@meta.data<- new_meta

#library(devtools)
#install_github('immunogenomics/presto')
# this counts is the TF-IDF normalized and log transformed counts
# see how long it takes for 1.3 million features

message("Finding differential accessible regions by presto")
start_time <- Sys.time()
res<- wilcoxauc(recount.atac, 'cluster_id', seurat_assay = 'ATAC', assay = "counts")
end_time <- Sys.time()

message("presto takes")
end_time - start_time

#Time difference of 2.933884 mins !! impressive

message(paste("writing differential peaks to", arguments$output, "\n"))
write_tsv(res %>% arrange(padj, desc(auc)), arguments$output)


