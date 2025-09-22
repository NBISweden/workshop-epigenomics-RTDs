.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html



====================
Downstream Analyses
====================


Ater the initial data analysis when the occupancy sites are determined, it is time to annotate the results with nearest genes, identify sequence motifs enriched in the peaks, produce visualisations, perform statistical analyses of differential occupancy and functionally annotate the results with gene sets. This part of the course will guide you through some of these steps.


.. image:: atac-chip-downstream/figures/downstream-proc.png
   			:width: 400px

*source for images for the collage: vignettes for ChIPseeker, deepTools, enrichplot*


:raw-html:`<br />`


Transcription factors (TFs) regulate gene expression by binding sites in the genome that often harbor a specific DNA motif. Identifying which TFs regulate defined set of genes is important to understanding the physiological consequences of altered chromatin status.

Tutorials pertaining to probing various aspects of TF activity are linked under:

.. toctree::
   :maxdepth: 1

   Transcription Factors <TF_tutorials.rst>


:raw-html:`<br />`

Other tutorials for downstream processing of ChIP-seq and ATAC-seq data are linked below. These analyses can be applied to other data types as well.


.. toctree::
   :maxdepth: 1

   Peak Annotation <atac-chip-downstream/PeakAnnot_tsao2022.fulldata_rtds.12ix2025.rst>
   Differential Accessibility in ATAC-seq data <atac-chip-downstream/PeakDA_tsao2022.fulldata_rtds.12ix2025.rst>
   Detection of differential occupancy in ChIP-seq data using DiffBind <diffBind/lab-diffBinding.rst>
   Signal Visualisation in Functional Genomics <visualisation/lab-visualisation.rst>
   

.. Genomic overlaps <genomeOverlap/lab-genomicGoverlaps.md>


:raw-html:`<br />`


Workflow managers have been used to control execution of complex analytical pipelines; in this tutorial we offer an introduction to standard data processing pipelines from the nf-core collection.


.. toctree::
   :maxdepth: 1

   Nextflow and nf-core pipelines <nextflow.rst>
