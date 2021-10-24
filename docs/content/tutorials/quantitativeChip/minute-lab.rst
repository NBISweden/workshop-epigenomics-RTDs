.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

MINUTE-ChIP
===========

.. Contents
.. ========

.. contents:: 
    :local:


Background
----------

Quantitative ChIP-seq methods allow us to compare the effect of DNA-protein interactions across samples. In our lab we developed **MINUTE-ChIP**, a multiplexed quantitative ChIP-seq method. 

To process the data with produce with this method, we developed Minute. Minute takes multiplexed MINUTE-ChIP FASTQ files and outputs a set of normalized, scaled bigWig files that can be directly compared across conditions.

We are going to analyse H3K27m3 data from human ESCs, where we have shown that Naive cells H3K27m3 landscape looks very different when compared to Primed cells in a quantitative manner.

This tutorial has two parts:

1. **Primary analysis**: Run minute workflow on a set of raw FASTQ files.
2. **Downstream analysis**: Perform analyses on the resulting bigWig files using usual bioinformatics tools.

In order to make this tutorial last a reasonable time, the **primary analysis** part is independent from the second part, as it runs only on a subset of the data you can see in the second part. During the **primary analysis**
you will see H3K27m3 data output for replicate 1 of our experiments. In the **downstream analysis** you can look at the three replicates together.


What you will learn:

- Apply the Minute data-processing pipeline to a set of MINUTE-ChIP fastq files.
- Visualize differences between unscaled and scaled bigWig files, using standard tools like IGV, seqtools. 
- Understand how quantitative ChIP results are different from traditional ChIP, using as example H3K27m3 data on human ESCs.


Setup
-----

Minute runs on `Snakemake <https://snakemake.readthedocs.io/en/stable>`_, a workflow management system that is based on python and inspired on makefile rules. 


Conda environment
^^^^^^^^^^^^^^^^^

To make sure you have all the necessary tools to run Minute, you can set up a conda environment using a :code:`environment.yml` specification file:


.. code-block:: bash

  module load bioinfo-tools
  module load conda

  conda env create -f /proj/epi2021/nobackup/minute_chip/minute-main/environment.yml -n minute_lab 


.. warning::
  There seems to be a problem with recent versions of :code:`tbb`, a :code:`bowtie2` dependency. We haven’t fixed yet this in the main repo, but you can downgrade tbb to a compatible version:

  .. code-block:: bash

    conda activate minute_lab
    conda install tbb=2020.2

  This will ask you if you want to downgrade :code:`tbb`. Say yes.



Primary analysis
----------------

Files
^^^^^ 

We will be using H3K27m3 data from our `publication <https://www.biorxiv.org/content/10.1101/2021.08.21.457215v1>`_:

Kumar, B., Navarro, C., Winblad, N., Schell, J. P., Zhao, C., Lanner, F., & Elsässer, S. J. (2021). Polycomb Repressive Complex 2 shields naïve human pluripotent cells from trophectoderm differentiation. bioRxiv.

In particular, we are going to look at Naïve vs Primed cells, and as control we have EZH2-inhibitor treatment, which removes H3K27m3 from the cells, creating a baseline for technical background.

There are 3 replicates for each condition. In the first part of the tutorial you will run replicate 1, and results for all replicates are available in the second part.


.. code-block:: bash
  
  # Create your primary analysis directory
  mkdir -p my_primary/fastq
  cd my_primary/fastq

  # Create symlinks to our fastq files
  for i in /proj/epi2021/nobackup/minute_chip/primary/*.fastq.gz; do ln -s ${i}; done
  cd ..
  cp /proj/epi2021/nobackup/minute_chip/primary/*.tsv
  cp /proj/epi2021/nobackup/minute_chip/primary/*.yaml


Now, this is how your file structure should look like:

- fastq/ - Contains all the fastq.gz files in the table below.
- libraries.tsv
- groups.tsv
- config.yaml


Minute needs three configuration files to run:

:code:`config.yaml`: Contains information about reference mapping: where the fasta files and bowtie2 indexes are, and a blocklist to remove artifact-prone regions before scaling:

.. code-block:: yaml

  references:
    hg38:  # Arbitrary name for this reference. This is also used in output file names.
      # Path to a reference FASTA file (may be gzip-compressed).
      # A matching Bowtie2 index must exist in the same location.
      fasta: "/proj/epi2021/nobackup/minute_chip/reference/hg38.fa"

      # Path to a BED file with regions to exclude
      exclude: "/proj/epi2021/nobackup/minute_chip/reference/hg38.blocklist.bed"

  # Length of the 5' UMI
  umi_length: 6

  # Fragment length (insert size)
  fragment_size: 150

:code:`libraries.tsv`: Contains information about the demultiplexing. In our case, the barcodes are skipped because we have the already demultiplexed FASTQ files. The raw FASTQ
mate 1 contains a 6nt UMI followed by a 8nt barcode that identifies the sample.

.. code-block::

  H3K27m3_Naive   1       .       H3K27m3-ChIP_H9_naive_rep1
  H3K27m3_Primed  1       .       H3K27m3-ChIP_H9_primed_rep1
  H3K27m3_Naive_EZH2i     1       .       H3K27m3-ChIP_H9_naive_EZH2i_rep1
  H3K27m3_Primed_EZH2i    1       .       H3K27m3-ChIP_H9_primed_EZH2i_rep1
  Input_Naive     1       .       IN-ChIP_H9_naive_rep1
  Input_Primed    1       .       IN-ChIP_H9_primed_rep1
  Input_Naive_EZH2i       1       .       IN-ChIP_H9_naive_EZH2i_rep1
  Input_Primed_EZH2i      1       .       IN-ChIP_H9_primed_EZH2i_rep1


:code:`groups.tsv`: Contains *scaling* information. Reads are normalized to matching sample input read counts, and in each scaling group, the first sample is used as reference. This has two implications:

1. Reference sample is normalized to 1x genome coverage.
2. Rest of samples values are directly comparable to the reference and across themselves.

Additionally, we may have some spike-in data from another reference, so Minute allows to map to different references in the same run. So :code:`groups.tsv` has also 

.. code-block::

  H3K27m3_Naive   1       Input_Naive     H3K27m3 hg38
  H3K27m3_Naive_EZH2i     1       Input_Naive_EZH2i       H3K27m3 hg38
  H3K27m3_Primed  1       Input_Primed    H3K27m3 hg38
  H3K27m3_Primed_EZH2i    1       Input_Primed_EZH2i      H3K27m3 hg38


.. list-table:: Table 1. Files used in this tutorial.
   :widths: 25 25 20 10 40
   :header-rows: 1

   * - IP
     - Cell type
     - Treatment
     - Rep
     - File
   * - H3K27m3
     - Naive
     - Untreated
     - 1
     - H3K27m3-ChIP_H9_naive_rep1_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - Untreated
     - 2
     - H3K27m3-ChIP_H9_naive_rep2_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - Untreated
     - 3
     - H3K27m3-ChIP_H9_naive_rep3_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - EZH2i
     - 1
     - H3K27m3-ChIP_H9_naive_EZH2i_rep1_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - EZH2i
     - 2
     - H3K27m3-ChIP_H9_naive_EZH2i_rep2_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - EZH2i
     - 3
     - H3K27m3-ChIP_H9_naive_EZH2i_rep3_R{1,2}.fastq.gz
   * - H3K27m3
     - Naive
     - Untreated
     - 1
     - H3K27m3-ChIP_H9_primed_rep1_R{1,2}.fastq.gz
   * - H3K27m3
     - Primed
     - Untreated
     - 2
     - H3K27m3-ChIP_H9_primed_rep2_R{1,2}.fastq.gz
   * - H3K27m3
     - Primed
     - Untreated
     - 3
     - H3K27m3-ChIP_H9_primed_rep3_R{1,2}.fastq.gz
   * - H3K27m3
     - Primed
     - EZH2i
     - 1
     - H3K27m3-ChIP_H9_primed_EZH2i_rep1_R{1,2}.fastq.gz
   * - H3K27m3
     - Primed
     - EZH2i
     - 2
     - H3K27m3-ChIP_H9_primed_EZH2i_rep2_R{1,2}.fastq.gz
   * - H3K27m3
     - Primed
     - EZH2i
     - 3
     - H3K27m3-ChIP_H9_primed_EZH2i_rep3_R{1,2}.fastq.gz

   * - Input
     - Naive
     - Untreated
     - 1
     - IN-ChIP_H9_naive_rep1_R{1,2}.fastq.gz
   * - Input
     - Naive
     - Untreated
     - 2
     - IN-ChIP_H9_naive_rep2_R{1,2}.fastq.gz
   * - Input
     - Naive
     - Untreated
     - 3
     - IN-ChIP_H9_naive_rep3_R{1,2}.fastq.gz
   * - Input
     - Naive
     - EZH2i
     - 1
     - IN-ChIP_H9_naive_EZH2i_rep1_R{1,2}.fastq.gz
   * - Input
     - Naive
     - EZH2i
     - 2
     - IN-ChIP_H9_naive_EZH2i_rep2_R{1,2}.fastq.gz
   * - Input
     - Naive
     - EZH2i
     - 3
     - IN-ChIP_H9_naive_EZH2i_rep3_R{1,2}.fastq.gz
   * - Input
     - Primed
     - Untreated
     - 1
     - IN-ChIP_H9_primed_rep1_R{1,2}.fastq.gz
   * - Input
     - Primed
     - Untreated
     - 2
     - IN-ChIP_H9_primed_rep2_R{1,2}.fastq.gz
   * - Input
     - Primed
     - Untreated
     - 3
     - IN-ChIP_H9_primed_rep3_R{1,2}.fastq.gz
   * - Input
     - Primed
     - EZH2i
     - 1
     - IN-ChIP_H9_primed_EZH2i_rep1_R{1,2}.fastq.gz
   * - Input
     - Primed
     - EZH2i
     - 2
     - IN-ChIP_H9_primed_EZH2i_rep2_R{1,2}.fastq.gz
   * - Input
     - Primed
     - EZH2i
     - 3
     - IN-ChIP_H9_primed_EZH2i_rep3_R{1,2}.fastq.gz


Normally, we run Minute on the multiplexed data. However, our pipeline can skip that step and go directly to the mapping and scaling.


Running Minute
^^^^^^^^^^^^^^ 

So if you already got your files, you need to run:

.. code-block:: bash

  conda activate minute_lab

  # Move to the directory where you copied the files
  cd my_primary

  # Run snakemake
  snakemake -p /proj/epi2021/nobackup/minute_chip/minute-main/Snakefile -j 4


:code:`-j` is the number of jobs or cores used by Snakemake. Depending on how many there are available on your node, you can raise this value.
The amount of files in this part of the tutorial is small enough to be possible to run in a local computer. For 4 out of 8 cores running on my laptop
(intel i7), this took around 4 hours to run. If you run this locally, consider not to use all the available cores you have, since you still need
to run other things on the side and it may eat up your RAM memory as well (more tasks means usually more memory use).


Essentially, the steps performed by Minute are:

- Demultiplex the reads and remove contaminated sequences (this is skipped in this execution).
- Map each condition to a reference genome.
- Deduplicate the reads.
- Remove excluded regions (such as artifact-prone regions, repeats, etc).
- Calculate scaling factors based on number of reads mapped and matching input conditions.
- Generate 1x coverage and scaled bigWig files from alignment using the calculated scaling factors.
- QC at every step (fastQC, Picard insert size metrics, duplication rates, etc) are gathered and output in the form of MultiQC report.


.. note::
  If the pipeline crashes at some point and you want to resume where it ran:

  .. code-block:: bash

  snakemake -p /proj/epi2021/nobackup/minute_chip/minute-main/Snakefile -j 4 --rerun-incomplete


After the pipeline is run, you will have the following folders:

- :code:`final/`: Contains final files: bigWig files, BAM files and demultiplexed FASTQ files (in this case, the same as your input).
- :code:`reports/`: Some reports on QC and scaling.
- :code:`log/`: Log output from each step.
- :code:`stats/`: Some stats files generated at each step.


Scaling info
^^^^^^^^^^^^ 

Scaling info is very relevant output, you will see the following figure under reports:


.. image:: Figures/minute_00_scaling.png
  :target: Figures/minute_00_scaling.png
  :alt:

*Fig. 1: Global scaling for H3K27m3 replicate 1*

What you see here is that Naive has around 3 times as much H3K27m3 than Primed cells, and that EZH2i treatment removes the majority of H3K27m3.

IGV tracks
^^^^^^^^^^ 

You can take the final/bigwig files and look at them on IGV. Here you can see IGF2 gene, where once scaled, H3K27m3 decoration seems around the same values
Primed vs Naïve, information that is lost in unscaled files.


.. image:: Figures/minute_01_igv.png
  :target: Figures/minute_01_igv.png
  :alt:

*Fig. 2: IGV screenshot of bigWig tracks at IGF2 gene. Gray tracks are unscaled, blue tracks are scaled. Here, Primed looks higher than Naïve, but upon scaling, values are similar.*



.. image:: Figures/minute_02_igv.png
  :target: Figures/minute_02_igv.png
  :alt:

*Fig. 3: Overview of tracks. Gray tracks are unscaled, blue tracks are scaled.*


**Q: How is the global distribution of primed H3K27m3 changing upon scaling? Why do Naïve samples look the same both scaled and unscaled?**



.. note::
  Make sure you select all your tracks and click on *group autoscale* so all the scales match.



Downstream analysis
-------------------

Files
^^^^^ 

Now you will get a copy of all the bigWig files.


.. code-block:: bash
  
  # Create your primary analysis directory
  mkdir my_downstream
  cd my_downstream

  cp /proj/epi2021/nobackup/minute_chip/downstream/*.bw
  cp /proj/epi2021/nobackup/minute_chip/downstream/*.bed


There should be :code:`unscaled` and :code:`scaled` bigWig files, plus a set of genes marked as Bivalent: :code:`Bivalent_Court2017.hg38.bed`. This annotation
comes from:

Court, F., & Arnaud, P. (2017). An annotated list of bivalent chromatin regions in human ES cells: a new tool for cancer epigenetic research. Oncotarget, 8(3), 4110.


Additionally, some bigWig tracks are pooled. These ones are all the replicates pooled together.


Looking at bivalent genes
^^^^^^^^^^^^^^^^^^^^^^^^^

You can look at these tracks over a BED file using ``seqplots``. ``seqplots`` is an ``R`` package that can be installed
from ``Bioconductor``. It can be run as a shiny app on a browser or from ``Rstudio``. You can check
how to use it `in this link <https://bioconductor.org/packages/release/bioc/vignettes/seqplots/inst/doc/SeqPlotsGUI.html>`_.

.. note::

  You can look at the bigWig files locally in your computer. These files are not so large and the processing will not be
  so computationally demanding.


Basically you need to open ``Rstudio``, install it:

.. code-block:: R

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("seqplots")

And run it:

.. code-block:: R

  library(seqplots)
  seqplots()

This will open a shiny app where you can load bigWig files and BED files, and plot profile plots and heatmaps.


**Q: How do the scaled vs unscaled plots differ? What do you think that means?**

**Q: How do replicates look?**

**Q: Is this a general effect, or is it dominated by a few loci? Hint: You can use a heatmap for this**

