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
     - 1.2 (*)
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
     - 2.2.6
     - peak calling
     - 
   * - BEDTools
     - 2.29.2
     - peak calling
     - 
   * - BEDOPS
     - 2.4.3
     - peak calling
     - 
   * - epic2
     - 0.0.52 (**)
     - broad peaks
     - https://github.com/biocore-ntnu/epic2
   * - csaw
     - 1.24.3 (***)
     - broad peaks
     - 


In addition, the following Bioconductor libraries are used in the tutorials:

* annotation: ``org.Hs.eg.db``, ``TxDb.Hsapiens.UCSC.hg19.knownGene``, ``TxDb.Hsapiens.UCSC.hg38.knownGene``,``reactome.db``,``biomaRt``

* analysis:  ``DiffBind``, ``edgeR``, ``ChIPQC ``, ``ChIPpeakAnno``,``ChIPseeker``,``ReactomePA``,``clusterProfiler``


(*) 
You will notice that we use a conda environment for ``phantompeakqualtools``, along with two versions of of ``samtools`` loaded as modules. This is because we have encountered incompatibilities with how ``phantompeakqualtools`` saves temporary data and recent versions of ``R``, and also syntax differences between ``samtools`` versions and libraries required for processing *bam* files produced by newer ``samtools`` versions. Long sotry short, while a bit clunky, this solution works on our system.

Conda environment was produced using::

  module load conda/latest

  conda create --prefix /path/to/your/environments/xcor -c bioconda phantompeakqualtools=1.2.1.1


(**)
Newer versions of ``Pysam`` seem to throw errors when used with ``epic2``. We have tested the configuration below and it works on our system.

Conda environment was produced using::

  module load conda/latest

  conda create --prefix /path/to/your/environments/epic_2b -c conda-forge python=3.7

  pip install Cython
  pip install pysam==0.15.2
  pip3 install epic2


(***)
This was tested under R run from conda environment with this package installed:

``r-v8                      3.4.0             r40h4312b11_0    conda-forge``


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


Downstream Processing
=======================

.. list-table:: Requirements for Downstream Processing labs.
   :widths: 25 25 25 25
   :header-rows: 1

   * - Software
     - Version
     - Lab
     - Comment
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
   :widths: 25 25
   :header-rows: 1   
   
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