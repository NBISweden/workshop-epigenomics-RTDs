.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


=========
ATAC-seq
=========


ATAC-seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) is a method for determining chromatin accessibility across the genome. It utilizes a hyperactive Tn5 transposase to insert sequencing adapters into open chromatin regions (Figure 1). High-throughput sequencing then yields reads that cluster within these regions of increased accessibility.

ATAC-seq can been used to

* profile open chromatin regions;

* identify differentially accessible regions (for multiple-sample experiments);

* identify transcription factor (TF) footprints;

* infer nucleosome positioning (in the TSS region).

:raw-html:`<br />`

.. list-table:: Figure 1. Overview of ATAC-seq (Buenrostro et al., 2015).
   :widths: 60
   :header-rows: 0

   * - .. image:: ATACseq/figures/atac-seq-figure1.png
   			:width: 400px


:raw-html:`<br />`

The tutorials in this section cover the steps oulined on Figure 2.


.. list-table:: Figure 2. Overview of ATAC-seq analysis.
   :widths: 60
   :header-rows: 0

   * - .. image:: ATACseq/figures/workflow-analysis.png
   			:width: 400px


:raw-html:`<br />`

These are tutorials for analysis and peak-dependent QC of ATAC-seq data:

.. toctree::
   :maxdepth: 1

   ATACseq: peak detection <ATACseq/lab-atacseq-bulk.rst>

