# pyflow-scATACseq
snakemake workflow for post-processing scATACseq data



Notes:

more thoughts on pileup and shifting reads
https://github.com/deeptools/deepTools/issues/453
maxFragmentLength

https://github.com/deeptools/deepTools/issues/370

https://groups.google.com/forum/#!topic/deeptools/JU9itiT5rYk

You probably want “–Offset 1”

The above is not shifting exactly, for shifting use
https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html

--ATACshift

filter out fragments