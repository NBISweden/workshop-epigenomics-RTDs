.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html




==========================================
Quality Control for Functional Genomics
==========================================


This tutorial is a continuatin of :doc:`Data preprocessing <data-preproc>`.



**Learning outcomes**

- understand and apply signal structure agnostic quality control methods used in functional genomics on an example of ATAC-seq data

- become accustomed to work on Rackham cluster


:raw-html:`<br />`


.. contents:: Table of Contents
   :depth: 1
   :local:
   :backlinks: none



:raw-html:`<br />`




.. image:: figures/workflow-proc.png
   			:width: 600px


The aim of this part of the data analysis workflow is to perform general signal structure agnostic (i.e. peak - independent) quality control (lower-right part of the concept map). These include:

* assessment of read coverage along the genome;

* replicate congruency.



:raw-html:`<br />`





Cumulative Enrichment
========================


`Cumulative enrichment <http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html>`_, aka BAM **fingerprint**, is a way of assesing the quality of signal concentrated predeminantly at small fraction of a reference (such as peaks detected in ATAC-seq and ChIP-seq). It determines how well the signal in the sample can be differentiated from the background.

Cumulative enrichment is obtained by sampling indexed BAM files and plotting a profile of cumulative read coverages for each sample. All reads overlapping a window (bin) of the specified length are counted; these counts are sorted and the cumulative sum is finally plotted.


To compute cumulative enrichment for processed bam files in our ATAC-seq data set (assuming we are in drectory ``analysis``):


.. code-block:: bash

	mkdir deepTools
	cd deepTools

	#link necessary files to avoind long paths in commands
	ln -s ../../data/proc/* .

	module load deepTools/3.3.2

	plotFingerprint --bamfiles ENCFF363HBZ.chr14.bam ENCFF398QLV.chr14.bam ENCFF045OAB.chr14.bam ENCFF828ZPN.chr14.bam \
	 --binSize=1000 --plotFile NKcellsATAC_chr14.fingerprint.pdf \
	 --labels ENCFF363HBZ ENCFF398QLV ENCFF045OAB ENCFF828ZPN -p 8 &> fingerprint.log


You can copy the resulting file to your local system to view it.


Have a look at ``NKcellsATAC_chr14.fingerprint.pdf``, read ``deepTools`` `What the plots tell you <http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#what-the-plots-tell-you>`_ and answer

- does it indicate a good sample quality, i.e. signal present in narrow regions?


.. admonition:: Fingerprint for ATAC-seq signal in NK cells.
   :class: dropdown, warning

   .. image:: figures/NKcellsATAC_chr14.fingerprint.png
          :width: 300px



Replicate Clustering
========================

**To assess overall similarity between libraries from different samples** one can compute sample clustering heatmaps using
`multiBamSummary <http://deeptools.readthedocs.io/en/latest/content/tools/multiBamSummary.html>`_ and `plotCorrelation <http://deeptools.readthedocs.io/en/latest/content/tools/plotCorrelation.html>`_ in bins mode from ``deepTools``.

In this method the genome is divided into bins of specified size (``--binSize`` parameter) and reads mapped to each bin are counted. The resulting signal profiles are used to cluster libraries to identify groups of similar signal profile.

We chose to compute pairwise Spearman correlation coefficients for this step, as they are based on ranks of each bin rather than signal values.



.. code-block:: bash

	multiBamSummary bins --bamfiles ENCFF363HBZ.chr14.bam ENCFF398QLV.chr14.bam ENCFF045OAB.chr14.bam ENCFF828ZPN.chr14.bam \
	 --labels ENCFF363HBZ ENCFF398QLV ENCFF045OAB ENCFF828ZPN \
	 --outFileName multiBamArray_NKcellsATAC_chr14.npz --binSize 5000 -p 8 &> multiBamSummary.log


	plotCorrelation --corData multiBamArray_NKcellsATAC_chr14.npz \
	 --plotFile NKcellsATAC_chr14_correlation_bin.pdf --outFileCorMatrix NKcellsATAC_chr14_correlation_bin.txt \
	 --whatToPlot heatmap --corMethod spearman


You can copy the resulting file to your local system to view it.

What do you think?

- which samples are similar?

- are the clustering results as you would have expected them to be?


.. admonition:: Correlation of binned ATAC-seq signal in NK cells.
   :class: dropdown, warning

   .. image:: figures/NKcellsATAC_chr14_correlation_bin.png
          :width: 300px




In addition to these general procedures, several specialised assay - specific quality metrics exist, which probe signal characteristics related to each method. These are **key QC metrics** to evaluate the experiment and should always be colleced during the QC step. The method specific tutorials are: :doc:`ATACseq <data-qc-atac>` and :doc:`ChIPseq <data-qc-chip.rst>`. 

We can now follow with :doc:`ATACseq specifc <data-qc-atac>` QC methods.



