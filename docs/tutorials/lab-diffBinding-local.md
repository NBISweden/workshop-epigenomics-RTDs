---
layout: default
title:  'ChIP-seq down-stream analysis'
---

# ChIP-seq down-stream analysis

### Learning outcomes
- obtain differentially bound sites with `DiffBind`
- annotate differentially bound sites with nearest genes and genomic features with `ChIPpeakAnno`
- perform functional enrichment analysis to identify predominant biological themes with `ChIPpeakAnno` and `reactome.db`

## Content
[Introduction](#Introduction)

[Data & Methods](#DataMethods)

[Setting-up](#Setting-up)

[Differential binding](#DB)
* [Installing DiffBind](#DB_install)
* [Running DiffBind](#DB_run)

[Functional analysis](#FA_local)
* [Installing ChIPpeakAnno](#FA_install)
* [Running ChIPpeakAnno](#FA_run)
* [Loading GO & REACTOME](#FA_go_install)
* [GO and REACTOME](#FA_go_run)


[Concluding remarks and next steps](#Next)

[Appendix: figures](#Next)

## Introduction <a name="Introduction"></a>
Welcome back to the second part of the tutorial. In the first part we have learnt how to access the quality of ChIP-seq data and we how to derive a consensus peakset for downstream analyses.

In this part we will learn how to place our peaks in a biological context, by identifying differentially bound sites between two sample groups and annotating these sites to find out predominant biological themes separating the groups.

## Data & Methods <a name="DataMethods">
We will continue using the same data as in the first part of the tutorial. Please note that usually **three biological replicates** are the **minimum requirement** for statistical analysis such as in factor occupancy.

_The ENCODE data we are using have only two replicates and we are using them to demonstrate the tools and methodologies. No biological conclusions should be drawn from them, or as a matter of fact, from any other dataset with duplicates only. Just because the tool computes does not make it right!_

## Setting-up  <a name="Setting-up">
If you have not done it already, install R and R-Studio. Refer back to [pre-course](../precourse) preparations for instructions.


## Differential binding <a name="DB">
### Intro

We will usage `Bioconductor` package [DiffBind](http://bioconductor.org/packages/release/bioc/html/DiffBind.html) to **identify sites that are differentially bound** between two sample groups.

The package includes _"functions to support the processing of peak sets, including overlapping and merging peak sets, counting sequencing reads overlapping intervals in peak sets, and identifying statistically significantly differentially bound sites based on evidence of binding affinity (measured by differences in read densities). To this end it uses statistical routines developed in an RNA-Seq context (primarily the Bioconductor packages edgeR and DESeq2 ). Additionally, the package builds on Rgraphics routines to provide a set of standardized plots to aid in binding analysis."_

 This means that we will **repeat deriving a consensus peakset** in a more powerful way before identifying differentially bound sites. Actually, defying the consensus peaks is an important step that takes up entire chapter in the `DiffBind` manual. We recommend reading entire section:  [6.2 Deriving consensus peaksets](http://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf).

So how does the differential binding affinity analysis work?

_"The core functionality of DiffBind is the differential binding affinity analysis, which enables binding sites to be identified that are statistically significantly differentially bound between sample groups. To accomplish this, first a contrast (or contrasts) is established, dividing the samples into groups to be compared. Next the core analysis routines are executed, by default using DESeq2. This will assign a p-value and FDR to each candidate binding site indicating confidence that they are differentially bound."_

### Setting-up `DiffBind`  <a name="DB_install">
Create a separate directory on your local computer where you would like to work and name it e.g. `diffBind`.

To be able to run `DiffBind` locally, we will need access to a set of files
- BAM files
- BED files with called peaks regions
- sample sheet information `.txt` file

Download these from Box (link under Alternative Files locations on the main site) or download these from Uppmax with `_scp_`command:

```bash

scp -r <username>@rackham.uppmax.uu.se:~/chipseq/data/bam/* .
scp -r <username>@rackham.uppmax.uu.se:~/chipseq/results/peaks_bed/* .
scp -r <username>@rackham.uppmax.uu.se:~/chipseq/analysis/R/samples_REST.txt .

```

You may want to place the downloaded files in the `diffBind` directory or at least keep a track of their location.

Also we need to modify `samples_REST.txt` so the pathways are pointing to the BAM and BED files on your local computer. Adjust the pathways in any editor your like.

Now, we can open R-Studio and set working directory to working folder e.g. `diffBind` folder by `Session -> Set Working Directory -> Choose Directory`. Now, all R commands will be in respect to this directory.

You can type commands directly in the Console window. A bit smarter way is to open a new R script under `File -> New File -> R Script` and type commands there, saving it from time to time. This way if you want to go back and repeat commands you can. To execute commands written in script, copy and paste commands to Console window and press Enter, press `Run` button in R-Studio or ask for a demo.

To use `DiffBind` package we need to install it first. To do so:
```bash

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffBind", version = "3.8")

```

If the above worked, we should be able to load DiffBind library:
```bash

library(DiffBind)

```

### Running `DiffBind`  <a name="DB_run">

We will now follow `DiffBind` example to obtain differentially bound sites, given our samples. You may want to open `DiffBind` tutorial and read section [3 Example: Obtaining differentially bound sites](http://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) while typing the command to get more information about each step.

```bash

# reading in the sample information (metadata)
samples = read.csv("samples_REST.txt", sep="\t")

#	inspecting the metadata
samples

#	creating an object containing data
res=dba(sampleSheet=samples, config=data.frame(RunParallel=FALSE))

# inspecting the object: how many peaks are identified given the default settings?
res

# counting reads mapping to intervals (peaks)
# at this step a normalisation is applied by the default set to: score=DBA_SCORE_TMM_MINUS_FULL
res.cnt = dba.count(res, minOverlap=2, score=DBA_SCORE_TMM_MINUS_FULL, fragmentSize=130)

# inspecting the object: notice the FRiP values!
res.cnt

# plotting the correlation of libraries based on normalised counts of reads in peaks
pdf("correlation_libraries_normalised.pdf")
plot(res.cnt)
dev.off()

# PCA scores plot: data overview
pdf("PCA_normalised_libraries.pdf")
dba.plotPCA(res.cnt,DBA_TISSUE,label=DBA_TISSUE)
dev.off()

# setting the contrast
res.cnt2 = dba.contrast(res.cnt, categories=DBA_TISSUE, minMembers=2)

# inspecting the object: how many contrasts were set in the previous step
res.cnt2

# performing analysis of differential binding
res.cnt3 = dba.analyze(res.cnt2)

# inspecting the object: which condition are most alike, which are most different, is this in line with part one of the tutorial?
dba.show(res.cnt3, bContrasts = T)

# correlation heatmap  using only significantly differentially bound sites
# choose the contrast of interest e.g. HeLa vs. neuronal (#1)
pdf("correlation_HeLa_vs_neuronal.pdf")
plot(res.cnt3, contrast=1)
dev.off()

# boxplots to view how read distributions differ between classes of binding sites
# are reads distributed evenly between those that increase binding affinity HeLa vs. in neuronal?
pdf("Boxplot_HeLa_vs_neuronal.pdf")
pvals <- dba.plotBox(res.cnt3, contrast=1)
dev.off()

# extracting differentially binding sites in GRanges
res.db1 = dba.report(res.cnt3, contrast=1)
head(res.db1)

# plotting overlaps of sites bound by REST in different cell types
pdf("binding_site_overlap.pdf")
dba.plotVenn(res.cnt3, 1:4, label1="HeLa",label2="neuron",label3="HepG2",label4="sknsh")
dev.off()

# finally, let's save our R session including the generated data. We will need everything in the next section
save.image("diffBind.RData")
```

## Functional analysis <a name="FA">

So now we have list of differentially bound sites for comparisons of interest but we do not know much about them besides the genomic location. It is time to them in a biological context. To do so, we will use another `Bioconductor` package [ChIPpeakAnno](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html).

ChIPpeakAnno _"is for facilitating the downstream analysis for ChIP-seq experiments. It includes functions to find the nearest gene, exon, miRNA or custom features such as the most conserved elements and other transcription factor binding sites supplied by users, retrieve the sequences around the peak, obtain enriched Gene Ontology (GO) terms or pathways. Starting 2.0.5, new functions have been added for finding the peaks with bi-directional promoters with summary statistics (peaksNearBDP), for summarizing the occurrence of motifs in peaks (summarizePatternInPeaks) and for adding other IDs to annotated peaks or enrichedGO (addGeneIDs). Starting 3.4, permutation test has been added to determine whether there is a significant overlap between two sets of peaks. In addition, binding patterns of multiple transcription factors (TFs) or distributions of multiple epigenetic markers around genomic features could be visualized and compared easily using a side-by-side heatmap and density plot._

Here, we will annotate deferentially bound sites, summarise them in a genomic feature context and obtain enriched GO terms and pathways.


### Setting-up `ChIPpeakAnno`  <a name="FA_install">

We will continue our R-Studio session. If you have logged-out or lost connection or simply want to start fresh follow setting up instructions for running DiffBind locally.

To install ChIPpeakAnno
```bash

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno", version = "3.8")

```

We will also need to load DiffBind results saved in the differential binding session. We will build on them.
```bash

load("diffBind.RData")

```

### Running `ChIPpeakAnno`  <a name="FA_run">

Like with DiffBind package there is a nice [ChIPpeakAnno tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html#annotate-peaks) that you can view along this exercise to read more about the various steps.

```bash

# Loading DiffBind library
# we will need it to extract interesting peaks for down-stream analysis
library(DiffBind)

# Loading ChIPpeakAnno library
library(ChIPpeakAnno)

# Loading TSS Annotation For Human Sapiens (GRCh37) Obtained From BiomaRt
data(TSS.human.GRCh37)

# Choosing the peaks for the interesting comparison, e.g.
data.peaks = dba.report(res.cnt3, contrast=1)
head(data.peaks)

# Annotate peaks with information on closest TSS using precompiled annotation data
data.peaksAnno=annotatePeakInBatch(data.peaks, AnnotationData=TSS.human.GRCh37)

# View annotated peaks: can you see the added information in comparsition to data.peaks?
head(data.peaksAnno)

# Saving results
write.table(data.peaksAnno, file="peaks_HeLa_vs_neuronal.txt", sep="\t", row.names=F)
```


### Loading GO and REACTOME database <a name="FA_go_run">
Locally, we can install few more R libraries and annotation data to inspect our peaks a bit more. We will need libraries `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg19.knownGene` and `reactome.db`. To install:


```bash

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("reactome.db", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")

```

### Enriched GO / REACTOME terms  <a name="FA_go_run">

```bash

library(org.Hs.eg.db)
library(reactome.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Peak distribution over genomic features
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peaks.featuresDist<-assignChromosomeRegion(data.peaksAnno, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf("peaks_featuresDistr_HeLa_vs_neuronal.pdf")
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(peaks.featuresDist$percentage, las=1, horiz=T)
dev.off()

# GO ontologies
peaks.go <- getEnrichedGO(data.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)

# Preview GO ontologies results
head(peaks.go$bp[, 1:2])
head(peaks.go$mf[, 1:2])
head(peaks.go$cc[, 1:2])

# REACTOME pathways
peaks.pathways <- getEnrichedPATH(data.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

# REACTOME pathways: preview data
head(peaks.pathways)

# REACTOME pathways: list all pathways
print(unique(peaks.pathways$PATH))

```

Feel free to build more on the exercises. Follow the [ChIPpeakAnno tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html#annotate-peaks) for ideas.


## Concluding remarks and next steps <a name="Next">

The workflow presented in the tutorials is quite common and it includes recommended steps for analysis of ChIP-seq data. Naturally, there may be different tools or ways to preform similar tasks. New tools are being developed all the time and no single tool can do it all.

In the extra labs we have prepared you can find for instance an alternative way of quality control of ChIP-seq data with R package called `ChIPQC` as well as alternative differential binding workflow with a packaged called `csaw`. Note, these labs were not extensively tested so you may need to experiment and draw from the knowledge gained in the main labs.

Also, there are more types of analyses one can do beyond the one presented here. A common further analysis, for instance, includes identification of short sequence motifs enriched in regions bound by the assayed factor (peaks). There are several tools available here and we recommend you test one or two with on the tutorial data: [Homer](http://homer.salk.edu/homer/), [GEM](http://groups.csail.mit.edu/cgs/gem/), [RSAT](http://floresta.eead.csic.es/rsat/peak-motifs_form.cgi)m [MEME](http://meme-suite.org/)

Above all, we recommend that you keep trying to analyze your own data. Practice makes perfect :)

----

## Appendix: figures <a name="Appendix">

![correlation_librarires_normalised](../figures/lab-diffBinding/correlation_libraries_normalised.pdf)

Fig: Correlation of libraries based on normalised counts of reads in peaks

----

![PCA](../figures/lab-diffBinding/PCA_normalised_libraries.pdf)

Fig: PCA scores plot: data overview using normalised counts of reads in peaks

----

![Heatmap](../figures/lab-diffBinding/correlation_HeLa_vs_neuronal.pdf)

Fig: Correlation heatmap  using only significantly differentially bound sites for HeLa and neuronal


----

![Boxplot](../figures/lab-diffBinding/Boxplot_HeLa_vs_neuronal.pdf)

Fig: Boxplots of reads distributions between HeLa and neuronal

----

![Venn](../figures/lab-diffBinding/binding_site_overlap.pdf)

Fig: Venn diagram of overlapping sites bound by REST in different cell types

----

![Features](../figures/lab-diffBinding/peaks_featuresDistr_HeLa_vs_neuronal.pdf)

Fig: Boxplots of reads distributions between HeLa and neuronal
