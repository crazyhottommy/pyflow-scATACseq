The ENCODE black list files were downloaded on 04/24/2019 
at https://sites.google.com/site/anshulkundaje/projects/blacklists

```bash
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
``` 


chromosome size files were gotten by `fetchChromSizes` from [UCSC utilites](http://hgdownload.soe.ucsc.edu/admin/exe/).

```bash
fetchChromSizes mm10 > mm10_chrom_size.txt
fetchChromSizes hg38 > hg38_chrom_size.txt
```