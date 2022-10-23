.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


=================
Peak Annotation
=================



**Learning outcomes**


- to perform GC-aware normalisation of ATAC-seq data using ``EDASeq``

- to detect differentially accessible regions using ``edgeR``



.. contents:: Contents
    :depth: 1
    :local:



Introduction
=============

In this tutorial we use an R / Bioconductor packages ``EDAseq`` and ``edgeR`` to perform normalisation and analysis of differential accessibility in ATAC-seq data.



Data & Methods
===============

We will build upon the main lab :doc:`ATACseq data analysis <../ATACseq/lab-atacseq-bulk>`:

* first we will summarise reads to detected peaks using **subset data** as an example to obtain a counts table; 

* we will use the counts table encompassing **complete data** for differential accessibility analysis; 




Setting-up
===========

We need access to file ``nk_merged_peaks.saf`` which we created in the :doc:`ATACseq data analysis <../ATACseq/lab-atacseq-bulk>`, part "Merged Peaks" and bam files with alignments.

Assuming we start at ``analysis``:


.. code-block:: bash

	mkdir counts_table
	cd counts_table

	ln -s ../peaks/consensus/nk_merged_peaks.saf
	ln -s ../../data_proc/* .

Data Summarisation 
=======================

We can now summarise the reads:

.. code-block:: bash

	module load subread/2.0.0
	featureCounts -p -F SAF -a nk_merged_peaks.saf --fracOverlap 0.2 -o nk_merged_peaks_macs_broad.counts ENCFF363HBZ.chr14.proc.bam ENCFF398QLV.chr14.proc.bam ENCFF828ZPN.chr14.proc.bam ENCFF045OAB.chr14.proc.bam

Let's take a look inside the counts table using ``head nk_merged_peaks_macs_broad.counts``.

.. admonition:: nk_merged_peaks_macs_broad.counts

   .. code-block:: bash

	# Program:featureCounts v2.0.0; Command:"featureCounts" "-p" "-F" "SAF" "-a" "nk_merged_peaks.saf" "--fracOverlap" "0.2" "-o" "nk_merged_peaks_macs_broad.counts" "ENCFF363HBZ.chr14.proc.bam" "ENCFF398QLV.chr14.proc.bam" "ENCFF828ZPN.chr14.proc.bam" "ENCFF045OAB.chr14.proc.bam" 
	Geneid	Chr	Start	End	Strand	Length	ENCFF363HBZ.chr14.proc.bam	ENCFF398QLV.chr14.proc.bam	ENCFF828ZPN.chr14.proc.bam	ENCFF045OAB.chr14.proc.bam
	nk_merged_macsBroadPeak_1	chr1	10004	10442	.	439	5	10	15	20
	nk_merged_macsBroadPeak_2	chr1	28945	29419	.	475	9	7	4	2
	nk_merged_macsBroadPeak_3	chr1	180755	181858	.	1104	7	8	5	2
	nk_merged_macsBroadPeak_4	chr1	191246	191984	.	739	0	3	1	0
	nk_merged_macsBroadPeak_5	chr1	778381	779290	.	910	78	88	34	30
	nk_merged_macsBroadPeak_6	chr1	817270	817490	.	221	0	0	0	0
	nk_merged_macsBroadPeak_7	chr1	826976	827974	.	999	63	40	28	23
	nk_merged_macsBroadPeak_8	chr1	838077	838576	.	500	0	1	0	0


We should remove the first line, as it can interfere with the way R reads in data:

.. code-block:: bash

	awk '(NR>1)' nk_merged_peaks_macs_broad.counts > nk_merged_peaks_macs_broad.counts.tsv


Differential Accessibility
============================

We now load R and packages:


.. code-block:: bash

	module load R_packages/4.1.1


We activate R console upon typing ``R`` in the terminal.


We begin by loading necessary libraries:

.. code-block:: R

	library(edgeR)
	library(EDASeq)

	library(GenomicAlignments)
	library(GenomicFeatures)

	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	library(wesanderson)

	library(Hmisc)
	library(dplyr)

	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

	ff = FaFile("/proj/epi2022/atacseq_proc/hg38ucsc/hg38.fa")

:raw-html:`<br />`



We can read in the data, format it and define experimental groups:

.. code-block:: R

	cnt_table = read.table("../counts/nk_merged.macs_broad.counts", sep="\t", header=TRUE, blank.lines.skip=TRUE)
	rownames(cnt_table)=cnt_table$Geneid

	groups = factor(c(rep("NK",2),rep("NKstim",2)))

	#this data frame contains only read counts to peaks on assembled chromosomes
	reads.peak = cnt_table[,c(7:10)]


We now prepare data with GC content of the peak regions for GC-aware normalisation.

.. code-block:: R

	gr = GRanges(seqnames=cnt_table$Chr, ranges=IRanges(cnt_table$Start, cnt_table$End), strand="*", mcols=data.frame(peakID=cnt_table$Geneid))

	peakSeqs = getSeq(x=ff, gr)
	
	gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
	
	#divide into 20 bins by GC content
	gcGroups = Hmisc::cut2(gcContentPeaks, g=20)
	mcols(gr)$gc = gcContentPeaks


Figure below shows that the accessibility measure of a particular genomic region is associated with its GC content. However, the slope and shape of the curves may differ between samples, which indicates that GC content effects are sample–specific and can therefore bias between–sample comparisons. 

To visualise GC bias in peaks:

.. code-block:: R


	lowListGC = list()
	for(kk in 1:ncol(reads.peak)){
	  set.seed(kk)
	  lowListGC[[kk]] = lowess(x=gcContentPeaks, y=log1p(reads.peak[,kk]), f=1/10)
	}

	names(lowListGC)=colnames(reads.peak)

	dfList = list()
	for(ss in 1:length(lowListGC)){
	  oox = order(lowListGC[[ss]]$x)
	  dfList[[ss]] = data.frame(x=lowListGC[[ss]]$x[oox], y=lowListGC[[ss]]$y[oox], sample=names(lowListGC)[[ss]])
	}
	dfAll = do.call(rbind, dfList)
	dfAll$sample = factor(dfAll$sample)

	p1.1 = ggplot(dfAll, aes(x=x, y=y, group=sample, color=sample)) +
	  geom_line(size = 1) +
	  xlab("GC-content") +
	  ylab("log(count + 1)") +
	  theme_classic()

	pdf("GCcontent_peaks.pdf")
	## plot just the average GC content
	p1.1
	dev.off()


.. admonition:: Counts vs GC contents in ATAC-seq peaks.
   :class: dropdown, warning

   .. image:: figures/GCcontent_peaks.png
          :width: 300px


We can see that GC content has an effect on counts within the peaks.



We have seen from analyses presented on lecture slides (https://www.biorxiv.org/content/10.1101/2021.01.26.428252v2)
that full quantile normalisation (FQ-FQ) implemented in package ``EDASeq`` is one of the methods which can mitigate the GC bias in detection of DA regions.

We'll detect differentially accessible regions using ``edgeR``. We will input the normalised GC content as offset to ``edgeR``.

To calculate the offsets, which correct for library size as well as GC content (full quantile normalisation in both cases):

.. code-block:: R
	
	reads.peak=as.matrix(reads.peak)

	dataOffset = withinLaneNormalization(reads.peak,y=gcContentPeaks,num.bins=20,which="full",offset=TRUE)
	dataOffset = betweenLaneNormalization(reads.peak,which="full",offset=TRUE)

We now use the statistical framework of ``edgeR``. We do not perform the internal normalisation (TMM) as usually, and instead we provide the offsets calculated by EDASeq.

.. code-block:: R

	design = model.matrix(~groups)

	d = DGEList(counts=reads.peak, group=groups)

	keep = filterByExpr(d)

	> summary(keep)
   		Mode   FALSE    TRUE 
	logical      21   54743 


	d=d[keep,,keep.lib.sizes=FALSE]

	d$offset = -dataOffset[keep,]
	d.eda = estimateGLMCommonDisp(d, design = design)
	d.eda = estimateGLMCommonDisp(d, design = design)
	fit = glmFit(d.eda, design = design)
	lrt.EDASeq = glmLRT(fit, coef = 2)

	DA_res=as.data.frame(topTags(lrt.EDASeq, nrow(lrt.EDASeq$table)))

The top DA peaks in stimulated vs non-stimulated NK cells::

	> head(DA_res)

	                                 logFC   logCPM       LR PValue FDR
	nk_merged_macsBroadPeak_29593 7.743577 4.014714 1648.016      0   0
	nk_merged_macsBroadPeak_9796  6.501470 4.527986 2801.485      0   0
	nk_merged_macsBroadPeak_20351 6.490681 4.934009 3551.762      0   0
	nk_merged_macsBroadPeak_12067 6.260759 4.441109 2593.194      0   0
	nk_merged_macsBroadPeak_11203 6.165875 4.511952 2684.820      0   0
	nk_merged_macsBroadPeak_53036 6.153595 4.023089 1922.240      0   0


Let's add more peak information:

.. code-block:: R

	DA_res$Geneid = rownames(DA_res)
	DA.res.coords = left_join(DA_res,cnt_table[1:4],by="Geneid")


We can check how well the GC correction worked:


.. code-block:: R

	gcGroups.sub=gcGroups[keep]
	dfEdgeR = data.frame(logFC=log(2^lrt.EDASeq$table$logFC), gc=gcGroups.sub)

	pedgeR = ggplot(dfEdgeR) +
	  aes(x=gc, y=logFC, color=gc) +
	  geom_violin() +
	  geom_boxplot(width=0.1) +
	  scale_color_manual(values=wesanderson::wes_palette("Zissou1", nlevels(gcGroups), "continuous")) +
	  geom_abline(intercept = 0, slope = 0, col="black", lty=2) +
	  ylim(c(-1,1)) +
	  ggtitle("log2FCs in bins by GC content") +
	  xlab("GC-content bin") +
	  theme_bw()+ 
	  theme(aspect.ratio = 1)+
	  theme(axis.text.x = element_text(angle = 45, vjust = .5),
	        legend.position = "none",
	        axis.title = element_text(size=16))

	ggsave(filename="log2FC_vs_GCcontent.pdf",plot=pedgeR ,path=".",device="pdf")


:raw-html:`<br />`

.. admonition:: Dependence of log2FC on GC content in ATAC-seq.
   :class: dropdown, warning

   .. image:: figures/log2FC_vs_GCcontent.png
          :width: 300px


It seems that FQ-FQ normalisation did not completely remove the effect of GC content on log2FC in thie dataset. However, these effects are somewhat mitigated, you can compare this plot to one obtained by using the standard TMM normalisation.

The reason why the GC effects are not completely removed in this case may be that the DA analysis is not performed on properly replicated data; we should have at least 3 replicates per condition, and we only have two.

