.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

.. modified from
.. https://training.galaxyproject.org/topics/epigenetics/tutorials/atac-seq/tutorial.html

============
ATAC-seq
============

ATAC-seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) is a method for determining chromatin accessibility across the genome. It utilizes a hyperactive Tn5 transposase to insert sequencing adapters into open chromatin regions. High-throughput sequencing then yields reads that indicate these regions of increased accessibility.



.. list-table:: Figure 1. Overview of ATAC-seq (Buenrostro et al., 2015).
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/atac-seq-figure1.png
   			:width: 400px




The sample processing workflow is similar to that of ChIP-seq. There are some small changes, however, and we'll discuss them here.

First of all, the sample should be sequenced using the paired end (PE) protocol. Single end sequencing for ATAC-seq is strongly discouraged, for reasons mentioned below.

Prior to aligning reads to the reference genome, the reads must be properly trimmed off adapters - because we will allow "dovetailing" (with the mates seemingly extending "past" each other) of read pairs during alignment:

.. code-block:: bash
	
	<--------------------Mate 1-----------------------
	AGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCATAGC
	  CTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCATAGCCT
	  ----------------------Mate 2--------------------->




Data
======

In this tutorial we will use data from the study of Buenrostro et al. 2013, the first paper on the ATAC-Seq method. The data is from a human cell line of purified CD4+ T cells. The original dataset had 2 x 200 million reads and would be too big to process in a training session, so the original dataset was downsampled to 200,000 randomly selected reads. Additionally, about 200,000 reads pairs that will map to chromosome 22 were included to have a good profile on this chromosome, similar to what you might get with a typical ATAC-Seq sample.


