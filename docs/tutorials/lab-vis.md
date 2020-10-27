---
layout: default
title:  'Signal visualisation with deepTools'
---

### Learning outcomes
Using `deepTools`
- to visualise ChIP signal in relation to annotated TSS

## Signal visualisation with deepTools
One more thing that may come useful when analysing ChIP-seq data is visualising ChIP signal in relation to annotated transcription start sites (TSS), here on chromosomes 1 and 2. To do so we will:
* convert bedgraph to bigWig using UCSC utilities
* calculate scores per genome regions using among others the bigWig file
* plot a heatmap of scores associated with genomic regions


In case you have logged out Uppmax:
```bash

ssh -Y <username>@rackham.uppmax.uu.se
interactive -A g2018030 -p core -n 4 --reservation=g2018030_WED
source ~/chipseq_env.sh

```

Assuming the same files structure as in the main data processing tutorial, create a separate directory in `~/chipseq/analysis` and navigate to it. Copy the files needed for this exercise.

```bash
cd ~/chipseq/analysis/
mkdir ~/chipseq/analysis/vis
cd ~/chipseq/analysis/vis

cp ../../hg19/chrom.sizes.hg19 chrom.sizes.hg19
cp ../bam_preproc/ENCFF000PED.chr12.cov.norm1x.bedgraph ./
```

To calculate scores per genome with deepTools [computeMatrix](http://deeptools.readthedocs.org/en/latest/content/tools/computeMatrix.html) we need [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) file that we can obtain by converting bedgraph using UCSC utilities:

```bash
module load ucsc-utilities/v287

bedGraphToBigWig ENCFF000PED.chr12.cov.norm1x.bedgraph chrom.sizes.hg19 hela_1.bw

module unload ucsc-utilities/v287
```

We can now compute the matrix of scores for visualisation using [computeMatrix](http://deeptools.readthedocs.org/en/latest/content/tools/computeMatrix.html). This tool calculates scores per genome regions and prepares an intermediate file that can be used with `plotHeatmap` and `plotProfiles`. Typically, the genome regions are genes, but any other regions defined in a BED file can be used. `computeMatrix` accepts multiple score files (bigWig format) and multiple regions files (BED format). This tool can also be used to filter and sort regions according to their score.

We will need a BED file with positions of TSS that we can copy to the working directory before running computeMatrix e.g.
```bash
module load deepTools/2.5.1

cp /sw/share/compstore/courses/ngsintro/chipseq/hg19/refGene_hg19_TSS_chr12_sorted_corr.bed ./

computeMatrix reference-point -S hela_1.bw \
-R refGene_hg19_TSS_chr12_sorted_corr.bed -b 5000 -a 5000 \
--outFileName matrix.tss.dat --outFileNameMatrix matrix.tss.txt \
--referencePoint=TSS --numberOfProcessors=max
```

We can now create a heatmap for scores associated with genomic regions, i.e. plot the binding profile around TSS
```bash

plotHeatmap --matrixFile matrix.tss.dat \
--outFileName tss.hela_1.pdf \
--sortRegions descend --sortUsing mean

```

Have a look at the `tss.hela_1.pdf`. What do you think?
