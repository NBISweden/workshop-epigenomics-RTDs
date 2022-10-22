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

These are tutorials for downstream processing of ChIP-seq and ATAC-seq data.


.. toctree::
   :maxdepth: 1

   Peak Annotation <atac-chip-downstream/lab-PeakAnnot.rst>
   Identification of sequence motifs <motifs/lab-motifs.rst>
   Differential Accessibility in ATAC-seq data <atac-chip-downstream/lab-DifAcces.rst>



.. ADD VISULALISATION

.. ADD DIFFBIND
