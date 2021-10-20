.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


======================
Dependencies
======================

This page lists software dependencies for the tutorials. This is useful for people without access to HPC Rackham on Uppmax.


.. contents:: 
   :depth: 1
   :local:
   :backlinks: none


ChIP-seq
============


.. list-table:: Requirements for ChIP-seq labs.
   :widths: 25 25 25 25
   :header-rows: 1

   * - Software
     - Version
     - Lab
     - Comment
   * - samtools
     - 0.1.19 & 1.18
     - processing and QC
     - both versions
   * - phantompeakqualtools
     - ?
     - processing and QC
     - https://github.com/kundajelab/phantompeakqualtools
   * - Picard
     - 2.23.4
     - processing and QC
     - 
   * - NGSUtils
     - 0.5.9
     - processing and QC
     -
   * - deepTools
     - 3.3.2
     - processing and QC; visualisation
     - 
   * - MACS
     - 3?
     - peak calling
     - not tested yet
   * - BEDTools
     - 2.29.2
     - peak calling
     - 
   * - BEDOPS
     - 2.4.3
     - peak calling
     - 



ATAC-seq
============

.. list-table:: Requirements for ATAC-seq labs.
   :widths: 25 25 25 25
   :header-rows: 1

   * - Software
     - Version
     - Lab
     - Comment
   * - samtools
     - 1.18
     - 
     - 
   * - 
     - 
     - 
     - 



.. NNN
.. ============

.. .. list-table:: Requirements for NNN labs.
..    :widths: 25 25 25 25
..    :header-rows: 1

..    * - Software
..      - Version
..      - Lab
..      - Comment
..    * - 
..      - 
..      - 
..      - 




Methylation
============

.. list-table:: Requirements for both methylation labs (all of these are R packages).
   :widths: 25 25 25 25
   :header-rows:    
   
   * - Software
     - Version
   * - R
     - 4.0.0
   * - limma
     - 3.44.3
   * - minfi
     - 1.34.0
   * - RColorBrewer
     - 1.1-2
   * - missMethyl
     - 1.22.0
   * - minfiData
     - 0.34.0
   * - Gviz
     - 1.32.0
   * - DMRcate
     - 2.2.3
   * - DMRcatedata
     - 2.10.0
   * - stringr
     - 1.4.0
   * - mCSEA
     - 1.12.0
   * - methylKit
     - 1.14.2
   * - genomation
     - 1.20.0
   * - GenomicRanges
     - 1.40.0