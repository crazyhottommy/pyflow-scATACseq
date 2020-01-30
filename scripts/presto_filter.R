
#! /usr/bin/env Rscript
'Filtering differential accessible regions 
Usage:
    presto_filter.R --presto=<presto_out> --cellnumber=<cell_number> --workers=<number> --n=<number> --pct.in.cutoff=<percent> --pct.each.out.cutoff=<percent> --max.num.off.target.cells=<number> --pct.out.cutoff=<percent> --padj.cutoff=<padj> <output>

Options:
    --presto=<presto_out>  The output of presto for all the peaks and all the clusters without filtering
    --cellnumber=<cell_number> input file path of the cell number information. a tsv file with a header. column 1 is the cluster id, column 2 is number of cells in that cluster.
    --workers=<number> number of workers for furrr to parallize. [default: 1]
    --n=<number>  peak in maximum how many clusters [default: 2]
    --pct.in.cutoff=<percent> percentage of cells have this peak in in group. [default: 15]
    --pct.each.out.cutoff=<percent> percentage of cells have this peak in each out group. [default: 10]
    --max.num.off.target.cells=<number> maximum number of cells have this peak. [default: 10]
    --pct.out.cutoff=<percent> percentage of cells have this peak in all out group combined. [default: 10]
    --padj.cutoff=<padj> adjusted p value for this peak calculated by presto. [default: 0.05]
    -h --help  Show this screen.
    -v --version  Show version.

Arguments:
    output   output file path of the filtered differential peaks
' -> doc

suppressWarnings(library(docopt))
suppressWarnings(library(furrr))
suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(purrr))
suppressWarnings(library(tidyr))

arguments <- docopt(doc, version = 'presto_filter.R v1.0\n\n')

# max 1G, see https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
options(future.globals.maxSize = 1000 * 1024^2)

filter_atac_per_peak<- function(df, cell_number, n = 2, pct.each.out.cutoff = 10, max.num.off.target.cells = 10){
        df<- inner_join(df, cell_number, by = c("group" = "cluster_id"))
        df %>%
                mutate(num_cells = pct_in * num/100)
        num_cluster_off_cells<- sum(df$num_cells >= max.num.off.target.cells)
        num_clusters<- sum(df$pct_in >= pct.each.out.cutoff)
        if (num_clusters <= n && num_clusters > 0){
                return (TRUE)
        } else {
                if (num_cluster_off_cells > 0 && num_cluster_off_cells <= n ) {
                        return (TRUE)
                } else {
                        return(FALSE)
                }
        }
} 


filter_atac_peaks<- function(res_nest, cell_number, n = 2, pct.in.cutoff = 15, 
                             pct.each.out.cutoff = 10,
                             max.num.off.target.cells = 10,
                             pct.out.cutoff = 10,
                             padj.cutoff = 0.05){
        
        indx<- furrr::future_map_lgl(res_nest$data, function(x) filter_atac_per_peak(df = x, cell_number = cell_number, n = n, pct.each.out.cutoff = pct.each.out.cutoff, max.num.off.target.cells = max.num.off.target.cells), .progress = TRUE)
        res_filter<- res_nest[indx,] %>%
                unnest()
        res_filter<- res_filter %>%
                filter(pct_in > pct.in.cutoff, padj < padj.cutoff, logFC >0) %>%
                group_by(group) %>%
                arrange(desc(logFC), padj) %>%
                tidyr::separate(feature, into=c("chr", "start", "end"), sep = "-") %>%
                nest() %>% 
                arrange(group)
        
}

plan(multiprocess, workers = as.numeric(arguments$workers))
message(paste0("using ", arguments$workers, " workers for purrr"))


## read in the txt file with two columns, cluster_id and cell number
cell_number<- read_tsv(arguments$cellnumber)

# this takes ~5mins 
res<- read_tsv(arguments$presto)

## this takes 4 mins
res_nest<- res %>%
        group_by(feature) %>%
        nest()



## takes 11 mins with 24 cpus. but make sure you have 
## 7.6G (size of the dataframe) * 24 = 182G  memory for your node.


filtered_peaks<- filter_atac_peaks(res_nest,  cell_number= cell_number, n = arguments$n, pct.in.cutoff = arguments$pct.in.cutoff, 
                                   pct.each.out.cutoff = arguments$pct.each.out.cutoff,
                                   max.num.off.target.cells = arguments$max.num.off.target.cells,
                                   pct.out.cutoff = arguments$pct.out.cutoff,
                                   padj.cutoff = arguments$padj.cutoff)

saveRDS(filtered_peaks, arguments$output)
walk2(filtered_peaks$data, filtered_peaks$group, function(x,y) write_tsv(x = x, path = paste0("08diff_peaks_all/",y,"_filtered_peaks.tsv") ))


