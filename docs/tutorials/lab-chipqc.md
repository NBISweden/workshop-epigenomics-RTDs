---
layout: default
title:  'ChIPQC'
---



## Introduction
Here, we will explore the alternative quality control workflow, using Bioconductor ChIPQC package. ChIPQC computes quality metrics for aligned data from ChIP-seq experiments. It also provides simple ways to generate a ChIP-seq experiment quality report which can be examined to asses the absolute and relative quality of individual ChIP-seq samples (and their associated controls, as well as overall quality of the experimental data.)

### Learning outcomes
Using `ChIPQC` package
- to generate a summary QC report for experimental sample groups
- to be able to understand and assess QC metrics and plots

## Setting-up
In principle one can run `ChIPQC` both on Uppmax or locally. However, today we will test the package locally.

<!-- **Uppmax**

To run on Uppmax, assuming the same files structure as for the [ChIP-seq data processing tutorial](processing) set pathway to R libraries installed on Uppmax, navigate to R directory and open R:
```bash

export R_LIBS="/sw/courses/ngsintro/chipseq/software/zzz_R_lib"

cd ~/chipseq/analysis/R

R

```

**Locally**-->

Follow set-up instructions from [Downstream analysis tutorial](lab-diffBinding-local), differential binding part. We will need the same files and we can work in the same directory.

Install `ChIPQC` library and any required dependencies
```bash

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ChIPQC", version = "3.8")

```

## Running ChIPQC
While running commands feel free to have a look at [ChIPQC package documentation](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf) to learn more about different steps and/or build upon them. Here we will just show you the very basics.

```bash

library(DiffBind)
library(ChIPQC)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#	reading in the sample information (metadata)
samples = read.csv("samples_REST.txt", sep="\t")

#	inspecting the metadata
samples

#	creating an object containing data
res=dba(sampleSheet=samples, config=data.frame(RunParallel=FALSE))

# inspecting the object
res

#	performing quality control
resqc = ChIPQC(res,annotation="hg19", config=data.frame(RunParallel=TRUE))

#	creating the quality control report in html format
ChIPQCreport(resqc)

```

Examine the html report.

What do you think?

Are results in line with the previous quality control workflow?

----------

The report can be also downloaded from Box [here](https://stockholmuniversity.box.com/s/c1lbrr1s1khw4ctiqfq0f9j2m1b6vp90)
