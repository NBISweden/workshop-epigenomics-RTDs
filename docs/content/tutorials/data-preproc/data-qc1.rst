.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html




==========================================
Quality Control for Functional Genomics
==========================================


This tutorial is a continuation of :doc:`Data preprocessing <data-preproc>`.



**Learning outcomes**

- apply quality control methods agnostic to signal structure, which are used in functional genomics, on an example of ATAC-seq data

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


Basic read count statistics were already collected in :doc:`Data preprocessing <data-preproc>`.


:raw-html:`<br />`


.. Important::

	We assume that the environment and directory structure has been already set in :doc:`Data preprocessing <data-preproc>`.




Cumulative Enrichment
========================


`Cumulative enrichment <https://deeptools.readthedocs.io/en/3.5.2/content/tools/plotFingerprint.html>`_, aka BAM **fingerprint**, is a way of assesing the quality of signal concentrated predeminantly in a small fraction of a genome (such as peaks detected in ATAC-seq and ChIP-seq). It determines how well the signal in the sample can be differentiated from the background.

Cumulative enrichment is obtained by sampling indexed BAM files and plotting a profile of cumulative read coverages for each sample. All reads overlapping a window (bin) of the specified length are counted; these counts are sorted and the cumulative sum is plotted.


To compute cumulative enrichment for processed bam files in our ATAC-seq data set (assuming we are in drectory ``analysis``, so if you have followed the previous tutorial, you should move one directory level up ``cd ..``). Here we use files preprocessed earlier:


.. code-block:: bash

	mkdir deepTools
	cd deepTools

	#link necessary files to avoid long paths in commands
	ln -s ../../data_proc/* .

	module load deepTools/3.3.2

	plotFingerprint --bamfiles ENCFF363HBZ.chr14.proc.bam ENCFF398QLV.chr14.proc.bam ENCFF045OAB.chr14.proc.bam ENCFF828ZPN.chr14.proc.bam \
	 --binSize=1000 --plotFile NKcellsATAC_chr14.fingerprint.pdf \
	 --labels ENCFF363HBZ ENCFF398QLV ENCFF045OAB ENCFF828ZPN -p 5 &> fingerprint.log


You can copy the resulting file to your local system to view it.


.. admonition:: Copying files from Rackham
   :class: dropdown, warning

   To copy files from Rackham you need to know the path to the file on Rackham  (i.e. the remote side). Type in the terminal::

   	pwd

   This gives you the path to the working directory, e.g::

   	$pwd

   	/proj/epi2023/nobackup/private/agata/tests/atacseq/analysis/deepTools

   To copy file ``NKcellsATAC_chr14.fingerprint.pdf`` to *current directory*, type in the **local** terminal::

   	scp <username>@rackham.uppmax.uu.se:/path/to/file .

   E.g.::

   scp agata@rackham.uppmax.uu.se:/proj/epi2023/nobackup/private/agata/tests/atacseq/analysis/deepTools/NKcellsATAC_chr14.fingerprint.pdf .

   **when connecting from abroad**

   You need to login in another session to be able to copy files, as 2FA does not work with ``scp``. This mock-login serves only to refresh your credentials and results in a few minutes of grace time, during which each session from the same host is accepted without the need to provide 2FA. This time window is sufficient for copying files.




Have a look at ``NKcellsATAC_chr14.fingerprint.pdf``, read ``deepTools`` `What the plots tell you <http://deeptools.readthedocs.io/en/3.5.2/content/tools/plotFingerprint.html#what-the-plots-tell-you>`_ and answer

- does it indicate a good sample quality, i.e. signal present in narrow regions?


.. admonition:: Fingerprint for ATAC-seq signal in NK cells.
   :class: dropdown, warning

   .. image:: figures/NKcellsATAC_chr14.fingerprint.png
          :width: 300px



Replicate Clustering
========================

**To assess overall similarity between libraries from different samples** one can compute sample clustering heatmaps using
`multiBamSummary <http://deeptools.readthedocs.io/en/3.5.2/content/tools/multiBamSummary.html>`_ and `plotCorrelation <http://deeptools.readthedocs.io/en/3.5.2/content/tools/plotCorrelation.html>`_ in bins mode from ``deepTools``.

In this method the genome is divided into bins of specified size (``--binSize`` parameter) and reads mapped to each bin are counted. The resulting signal profiles are used to cluster libraries to identify groups of similar signal profile.

We chose to compute pairwise Spearman correlation coefficients for this step, as they are based on ranks of each bin rather than signal values.

In this part we use bam files filtered previously, to save time.


.. code-block:: bash

	multiBamSummary bins --bamfiles ENCFF363HBZ.chr14.proc.bam ENCFF398QLV.chr14.proc.bam ENCFF045OAB.chr14.proc.bam ENCFF828ZPN.chr14.proc.bam \
	 --labels ENCFF363HBZ ENCFF398QLV ENCFF045OAB ENCFF828ZPN \
	 --outFileName multiBamArray_NKcellsATAC_chr14.npz --binSize 5000 -p 5 &> multiBamSummary.log


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



