.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


========================================
Transcription Factor Footprinting
========================================



**Learning outcomes**


- detect transcription factor binding signatures in ATAC-seq data


.. list-table:: Figure 1. Overview of TF footprinting in ATAC-seq data.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/TFfootprinting_tobias.png
   			:width: 600px




:raw-html:`<br />`




.. contents:: Contents
    :depth: 1
    :local:




:raw-html:`<br />`


Introduction
=============


While ATAC-seq can uncover accessible regions where transcription factors (TFs) *might* bind, reliable identification of specific TF binding sites (TFBS) still relies on chromatin immunoprecipitation methods such as ChIP-seq.
ChIP-seq methods require high input cell numbers, are limited to one TF per assay, and are further restricted to TFs for which antibodies are available. Therefore, it remains costly, or even impossible, to study the binding of multiple TFs in one experiment.

Similarly to nucleosomes, **bound TFs hinder cleavage of DNA**, resulting in **footprints**: defined regions of decreased signal strength within larger regions of high signal (Figure 1.).

Despite its compelling potential, a number of issues have rendered footprinting of ATAC-seq data cumbersome.
It has been described that enzymes used in chromatin accessibility assays (e.g., DNase-I, Tn5) are biased towards certain sequence compositions. If unaccounted for, this impairs the discovery of true TF footprints.

In this tutorial we use an R / Bioconductor packages ``ATACseqQC`` and ``MotifDb`` to detect TF binding signatures in ATAC-seq data. Please note this tutorial is merely an early attempt to determine whether TF bindng sites can be identified in ATAC-seq data, and *not a statistical framework for TF footrpinting*. It can be used as a QC step to detect expected TF sites in the data rather than for discovery of novel binding patterns.


:raw-html:`<br />`



Data & Methods
===============

We will build upon the lab :doc:`ATACseq specifc QC <../data-preproc/data-qc-atac>` and use the same data as for other ATAC-seq labs.


:raw-html:`<br />`


Setting-up
===========

We need access to bam file with shifted alignments ``shifted.bam`` which we created in the :doc:`ATACseq specifc QC <../data-preproc/data-qc-atac>`. This file contains alignments shifted +4 bps on the + strand and -5 bps on the - strand, to account for Tn5 transposition.


Assuming we start at ``analysis``:


.. code-block:: bash

	mkdir TF_footprnt
	cd TF_footprnt

	ln -s ../QC/splitBam/shifted.bam .
	ln -s ../QC/splitBam/shifted.bam.bai .

	module load R_packages/4.1.1


.. Hint::

	Please check first that file ``shifted.bam`` exists in this location: ``ls ../QC/splitBam/shifted.bam``. If it does, the output of this command is the path; if it does not you get "file does not exist". Depending on the directory structure, you may need to link the file like this:

	.. code-block:: bash

		ln -s ../../QC/splitBam/shifted.bam .
		ln -s ../../QC/splitBam/shifted.bam.bai .


We activate R console upon typing ``R`` in the terminal.

:raw-html:`<br />`

Detection of TF Binding Signatures
======================================

We begin by loading necessary libraries:


.. code-block:: R

	library(ATACseqQC)
	library(MotifDb)

	library(BSgenome.Hsapiens.UCSC.hg38)
	genome <- Hsapiens


We load data from shifted bam file:

.. code-block:: R

	shifted.bamFile="shifted.bam"

	# we will limit the analysis to chr14
	seqlev <- "chr14"


Let's first check signatures for a general TF **CTCF**. This is its motif as position weight matrix (PWM):

.. code-block:: R


	CTCF <- query(MotifDb, c("CTCF"))
	CTCF <- as.list(CTCF)
	print(CTCF[[1]], digits=2)

We now summarise the signal in the vicinity of CTCF motifs (100 bps up- and down-stream):

.. code-block:: R

	ctcf <- factorFootprints(shifted.bamFile, pfm=CTCF[[1]], 
	                         genome=genome,
	                         min.score="90%", seqlev=seqlev,
	                         upstream=100, downstream=100)



This function outputs signal mean values of coverage for positive strand and negative strand in feature regions, and other information which you can inspect using ``str(ctcf)``:

* ``spearman.correlation`` spearman correlations of cleavage counts in the highest 10-nucleotide-window

* ``bindingSites`` - GRanges object with detected bindng sites

* ``Profile.segmentation``


:raw-html:`<br />`


Let's inspect the statistics ``ctcf$spearman.correlation``::


	> ctcf$spearman.correlation
	$`+`

		Spearman's rank correlation rho

	data:  predictedBindingSiteScore and highest.sig.windows
	S = 3761558200, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	      rho 
	0.2653951 


	$`-`

		Spearman's rank correlation rho

	data:  predictedBindingSiteScore and highest.sig.windows
	S = 3801340300, p-value < 2.2e-16
	alternative hypothesis: true rho is not equal to 0
	sample estimates:
	      rho 
	0.2576259 


.. Strength of the site can be summarised by profile segmentation. Distal and proximal abundance (shown as red dash lines in the plot produce dby this function) are calculated by averaging the signal at the center of binding sites to the end of distal sites, and then calculated in the proximal and distal locations with respect to the motif.

.. Profile segmentation::

.. 	 ctcf$Profile.segmentation
.. 	          pos   distal_abun proximal_abun       binding 
.. 	   56.0000000     0.1481738     0.2627657     0.1392823 




The plot produced by this function is of signal mean values of coverage for positive strand and negative strand in feature regions.


.. code-block:: R

	pdf("ctcf_footprnt.pdf")
	sigs <- factorFootprints(shifted.bamFile, pfm=CTCF[[1]], 
	                         genome=genome,
	                         min.score="90%", seqlev=seqlev,
	                         upstream=100, downstream=100)
	dev.off()



We can generate similar plots for other TFs.

RFX5

.. code-block:: R

	RFX5 <- query(MotifDb, c("RFX5"))
	RFX5 <- as.list(RFX5)


	rfx5 <- factorFootprints(shifted.bamFile, pfm=RFX5[[1]], 
	                         genome=genome,
	                         min.score="90%", seqlev=seqlev,
	                         upstream=100, downstream=100)

	 rfx5$spearman.correlation$`+`$estimate
	 rfx5$spearman.correlation$`+`$p.value


	pdf("rfx5_footprnt.pdf")
	rfx5 <- factorFootprints(shifted.bamFile, pfm=RFX5[[1]], 
	                         genome=genome,
	                         min.score="90%", seqlev=seqlev,
	                         upstream=100, downstream=100)
	dev.off()



STAT3

.. code-block:: R


	STAT3 <- query(MotifDb, c("STAT3"))
	STAT3 <- as.list(STAT3)


	stat3 <- factorFootprints(shifted.bamFile, pfm=STAT3[[1]], 
		                         genome=genome,
		                         min.score="90%", seqlev=seqlev,
		                         upstream=100, downstream=100)

	stat3$spearman.correlation$`+`$estimate
	stat3$spearman.correlation$`+`$p.value

	stat3$Profile.segmentation
	
	pdf("stat3_footprnt.pdf")
	stat3 <- factorFootprints(shifted.bamFile, pfm=STAT3[[1]], 
		                         genome=genome,
		                         min.score="90%", seqlev=seqlev,
		                         upstream=100, downstream=100)
	dev.off()


.. profile segmentation for stat3

	.. ``stat3$Profile.segmentation``::

	..  stat3$Profile.segmentation
	..           pos   distal_abun proximal_abun       binding 
	..   88.00000000    0.04103757    0.05341195    0.04394870 



Which factors show evidence of binding enrichment in this data set?


.. list-table:: Figure 2. Examples of TF footprints.
   :widths: 40 40 40 
   :header-rows: 1

   * - CTCF
     - RFX5
     - STAT3
   * - .. image:: figures/ctcf_footprnt.png
   			:width: 300px
     - .. image:: figures/rfx5_footprnt.png
   			:width: 300px
     - .. image:: figures/stat3_footprnt.png
   			:width: 300px



image source: *https://doi.org/10.1038/s41467-020-18035-1* (Figure 1.)


.. spearman correlation of binding score (scoring by pwm) and signal in sliding windows



.. For transcription factor footprint analysis, pwmScore, plotFootprints and factorFootprints are implemented in ATACseqQC. It makes use of genomic sequences as BSgenome objects, available for various reference genomes, which can be efficiently accessed by methods in the BSgenome package [29], and of the position frequency matrices (PFMs) of binding motifs of transcription factors from the Jaspar database in the MotifDb package [30]. The footprint analysis also leverages the matchPWM function in the BSgenome package [29, 31] to search potential binding sites for a given DNA-binding protein, represent the matched genomic coordinates as GenomicRanges objects, and plot the motif as a sequence logo using the motifStack package [32]. The factorFootprints function first uses the matchPWM function to predict the binding sites with an input position weight matrix (PWM) for a DNA-binding protein. Next, it calculates and plots the average cutting rates for those binding sites and 100-bp flanking sequences.

