.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


=================================================
In-situ methods: CUT&RUN and ChIP data for CTCF
=================================================

.. Contents
.. ========

.. contents:: 
    :local:

Background
=============

CTCF (CCCTC-binding factor) is a very general architectural protein that mediates interchromosomal or intrachromosomal interactions. By recruiting additional structural proteins like condensin, they create chromatin loops and contact points. CTCF is a sequence-specific binding factor with a well-known sequence motif. There are 55,000-65,000 theoretical binding sites in the human genome (predicted by the consensus sequence). Of these, 5000 appear to be occupied by CTCF in every tissue/cell type, whereas up to 60% are used cell-type specifically. DNA methylation modulates CTCF binding (lowering CTCF affinity). Given the plethora of knowledge and datasets available one would think that there is a universal agreement on where CTCF binds in which cell type. However, different ChIP-Seq data and alternative profiling methods actually give different qualitative and certainly quantitative answers. 

So we will look at a comparison of CUT&RUN and different ChIP datasets to understand why and how different datasets report different binding sites.

.. image:: Figures/00_CTCF.png
	:target: Figures/00_CTCF.png
	:alt:

*Fig. 1*: Features of `CTCF binding sites in the genome <https://doi.org/10.1038/nrg3663>`_ [1]_.


Datasets
==========

The datasets used correspond to two publications, Skene and Henikoff [2]_ and Pugacheva et al. [3]_. They can be found on GEO (GSE84474 and GSE137216). See below a table with the samples included in this example:

.. list-table:: Table 1. Files in Skene et al. 2017 used in this exercise.
   :widths: 25 40
   :header-rows: 1

   * - GEO Accession
     - Sample
   * - GSM2247150
     - NBIS_Skene2017_K562_CTCF_ChIP_loMN
   * - GSM2247152 
     - NBIS_Skene2017_K562_CTCF_ChIP_medMN    
   * - GSM2247154 
     - NBIS_Skene2017_K562_CTCF_ChIP_hiMN     
   * - GSM2433137 
     - NBIS_Skene2017_K562_CTCF_CnR_15s       
   * - GSM2433139 
     - NBIS_Skene2017_K562_CTCF_CnR_3m        
   * - GSM2433140 
     - NBIS_Skene2017_K562_CTCF_CnR_45s       
   * - GSM2433141 
     - NBIS_Skene2017_K562_CTCF_CnR_5s        
   * - GSM2433142 
     - NBIS_Skene2017_K562_CTCF_CnR_9m        



In Skene and Henikoff, they compare CUT&Run to their own ChIP (special optimized native ChIP for CTCF they say). In Pugacheva et. all, they perform a crosslinking ChIP for CTCF. They’ve been using different antibodies just to be sure that they are ChIPing the right thing.

.. list-table:: Table 2. Files in Pugacheva et al. 2020 used in this exercise.
   :widths: 25 40
   :header-rows: 1

   * - GEO Accession
     - Sample
   * - GSM4073237
     - NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoC_Abs
   * - GSM4073240 
     - NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoN_Abs_Rep1    
   * - GSM4073238 
     - NBIS_Pugacheva2020_K562_ChIP_CTCF_RabbitC_Abs     
   * - GSM4073242 
     - NBIS_Pugacheva2020_K562_ChIP_Mix_of_CTCF_PolyC_Abs       
   * - GSM4073239 
     - NBIS_Pugacheva2020_K562_ChIP_IGG_Abs_Control        


You can find many more CTCF datasets In the `Cistrome database <http://cistrome.org/db>`_!


Data preprocessing
====================

Primary analysis of the initial FASTQ files was performed beforehand. Reads were mapped with ``bowtie2`` with default parameters. Resulting BAM files were deduplicated using Picard and blacklisted regions were removed. Peak files were generated using ``MACS2`` using very standard parameters: ``--qval cutoff 0.01``. In this case no replicates were used to call the peaks. BigWig files were generated with ``deepTools`` to 1x coverage. Resulting files are available under ``/sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR``

.. attention::
   You can download bigWig and peak annotations. Most of what we are going to do can be done locally on a regular laptop. When this is not the case, Uppmax-specific instructions will be given. In case something does not work properly, the output of most of these is available also in the workshop folder and in this documentation.

**Uppmax**: When running things on Uppmax, copy the files to your home directory:

.. code-block:: bash

    cd
    mkdir -p cnr_chip/bw
    cd cnr_chip/bw
    
    cp /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/bw/*.bw .
    
    # Or you can create symlinks for the bigWig files instead:
    for i in /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/bw/*.bw; do ln -s ${i}; done



Manual inspection of bigwig profiles with IGV
==============================================

Besides quality control tools and exploratory analyses, it is usually a good idea to have a look at the data without any further processing to have an idea about how it looks like. It can help spot technical issues or processing problems.

.. attention::
    When running these commands, the same relative path structure is expected to be kept. So if you
    made a directory for this tutorial, be sure you are in that directory when starting.

**Local**: If you have not done it yet, download a copy of the bigWig and peaks to your local laptop.

.. code-block:: bash

    # Wherever you have this directory, it's from now on your main working directory
    mkdir cnr_chip
    cd cnr_chip

    mkdir bw
    mkdir peaks

    scp <youruser>@rackham.uppmax.uu.se:/sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/bw/*.bw ./bw/

    scp <youruser>@rackham.uppmax.uu.se:/sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/peaks/*.narrowPeak ./peaks/
    

Open IGV and import the bigWig files (file > load from file...). You should see something like this:


.. image:: Figures/01_IGV.png
	:target: Figures/01_IGV.png
	:alt:

**Q: Note down some observations concerning following questions:**

- How are the samples similar and how are they different?
- Can you make a prediction about the peak calling result?


Go to this site ``chr1:27,389,733-27,973,162``
(corresponding to `Figure 8 <https://elifesciences.org/articles/21856/figures#fig8>`_ ).

**Q: Are peaks shared, or unique to a specific technique?**


Quality Control
================

Quality Control metrics are very important to understand whether the experiment worked or not and spot possible caveats that may come up before further analysis. As mentioned in other tutorials, some metrics that can be useful are:

- Library complexity estimates and duplication rates.
- Sample clustering (PCA, correlation).
- Insert size distribution.
- Cumulative enrichment.


Correlation plots
-----------------

One way to look at similarity between ChIP datasets is to partition the signal in bins of a fixed size and compute a pairwise correlation value between distributions. This can be done using ``deepTools``.

**Uppmax**: This is a time consuming step that would need to be done on Uppmax. Move to the directory where you copied bigWig files before. If you followed the same names, it would be:

.. code-block:: bash
    
    cd ~/cnr_chip

Load ``deepTools`` (and ``bioinfo-tools``) module beforehand:

.. code-block:: bash
    
    module load bioinfo-tools
    module load deepTools

``deepTools`` needs that you first compute a bin matrix. From this, many other things can be done, such as correlation plots and PCA.

You can calculate a bin matrix using the command below. Since it is a very long one-line command, it is split into several lines
so it is easier to read. You can split any shell command in lines using the back slash ``\``, which tells the shell parser that the
command continues on the next line. Just be sure that there are no white spaces after the back slash. 

When working on a data analysis with deepTools, if you are not sure if a command is correct, it is a good
idea to use the ``--region`` parameter, which will do the analysis only on a given genomic region,
and will run faster (failing fast is a good philosophy that will save you loads of time). 

For instance, this command has a ``--region chr1:300000:900000`` parameter, which will only calculate 
bins on that genomic region and will run immediately. Then you can check the output and verify that it is what you expected.

You can copy and paste the command below as it is on the terminal:

.. code-block:: bash

    # Labels basically follow alphabetical order (which is the way blob *.bw will expand).
    # But if you want to be totally sure, it’s better to spell out the order, as there may be 
    # some surprises: https://unix.stackexchange.com/questions/368318/does-the-bash-star-wildcard-always-produce-an-ascending-sorted-list
    
    multiBigwigSummary bins -b \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoC_Abs.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoN_Abs_Rep1.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_RabbitC_Abs.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_IGG_Abs_Control.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_Mix_of_CTCF_PolyC_Abs.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_hiMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_loMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_medMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_15s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_3m.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_45s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_5s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_9m.GRCh38.bw \
        -bs 5000 -p 4 -o ./bins_table.npz --outRawCounts ./bins_table.tab \
        --labels \
        Pugacheva2020_ChIP_MonoC \
        Pugacheva2020_ChIP_MonoN \
        Pugacheva2020_ChIP_RabbitC \
        Pugacheva2020_ChIP_Control \
        Pugacheva2020_ChIP_Mix \
        Skene2017_ChIP_hiMN \
        Skene2017_ChIP_loMN \
        Skene2017_ChIP_medMN \
        Skene2017_CnR_15s \
        Skene2017_CnR_3m \
        Skene2017_CnR_452 \
        Skene2017_CnR_5s \
        Skene2017_CnR_9m \
        --region chr1:300000:900000


.. note::
    The parameter ``--outRawCounts`` is not necessary and usually not generated, as the same values are saved in ``bins_table.npz`` in a way they occupy less space. But raw counts are text, so you can basically peek at the values directly using ``head`` or ``more``.


If the small test ran successfully, you can run ``multiBigWigSummary`` on the whole genome by running:

.. code-block:: bash
        
    multiBigwigSummary bins -b \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoC_Abs.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoN_Abs_Rep1.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_CTCF_RabbitC_Abs.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_IGG_Abs_Control.GRCh38.bw \
        ./bw/NBIS_Pugacheva2020_K562_ChIP_Mix_of_CTCF_PolyC_Abs.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_hiMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_loMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_ChIP_medMN.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_15s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_3m.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_45s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_5s.GRCh38.bw \
        ./bw/NBIS_Skene2017_K562_CTCF_CnR_9m.GRCh38.bw \
        -bs 5000 -p 4 -o ./bins_table.npz --outRawCounts ./bins_table.tab \
        --labels \
        Pugacheva2020_ChIP_MonoC \
        Pugacheva2020_ChIP_MonoN \
        Pugacheva2020_ChIP_RabbitC \
        Pugacheva2020_ChIP_Control \
        Pugacheva2020_ChIP_Mix \
        Skene2017_ChIP_hiMN \
        Skene2017_ChIP_loMN \
        Skene2017_ChIP_medMN \
        Skene2017_CnR_15s \
        Skene2017_CnR_3m \
        Skene2017_CnR_452 \
        Skene2017_CnR_5s \
        Skene2017_CnR_9m


The ``.npz`` matrix is then used by ``deepTools`` to produce other plots. For our correlation plot:

.. code-block:: bash

    plotCorrelation --corData ./bins_table.npz \
        --plotFile ./correlation_spearman_all.png \
        --whatToPlot heatmap \
        --corMethod spearman \
        --plotNumbers

This will generate a correlation plot based on genome-wide 5kb bins.


.. image:: Figures/02_corrplot.png
	:target: Figures/02_corrplot.png
	:alt:

**Q: Check out how the datasets cluster - does it make sense? Is the overall clustering following the biological target/control or underlying batch effect?**

.. note:: 
    If some step did not work, you can get the generated genome-wide bins files from: ``cp /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/tmp/bins* .``.



Cumulative enrichment
---------------------

Also known as fingeprint plots, these give a feeling about the signal to noise ratio of each signal. You
can understand more about what they exactly mean in `deepTools` `documentation <https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html#id6>`_.

**Uppmax**: You can plot this with ``deepTools`` as well. This requires the BAM files and takes quite a bit to compute. You can symlink the bam files from: ``/sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/bam/`` the same way as before:

.. code-block:: bash

    # Or your preferred folder
    cd ~/cnr_chip
    mkdir bam
    cd bam

    # In this case symlink is strongly preferred, as these files are bigger
    for i in /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/bam/*; do ln -s ${i}; done

.. attention:: 
    This is a time consuming step that is recommended to do on Uppmax. You can also leave this for later and look at the results.

You can plot them.

.. code-block:: bash
    
    cd ~/cnr_chip

    plotFingerprint -b \
      ./bam/NBIS_Skene2017_K562_CTCF_ChIP_hiMN.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_ChIP_loMN.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_ChIP_medMN.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_CnR_15s.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_CnR_3m.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_CnR_45s.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_CnR_5s.GRCh38.bam \
      ./bam/NBIS_Skene2017_K562_CTCF_CnR_9m.GRCh38.bam \
      -o ./fingerprint_Skene2017.pdf \
      --labels \
      Skene2017_ChIP_hiMN \
      Skene2017_ChIP_loMN \
      Skene2017_ChIP_medMN \
      Skene2017_CnR_15s \
      Skene2017_CnR_3m \
      Skene2017_CnR_452 \
      Skene2017_CnR_5s \
      Skene2017_CnR_9m \
      -p 4

The resulting plot should look like:

.. image:: Figures/03_fingerprint1.png
	:target: Figures/03_fingerprint1.png
	:alt:

**Q: What do you think it means in terms of quality of experiments?**
**Do you see different groups of samples? Go to IGV and browse through the tracks.**
**How does each group look? Fingerprint supports the difference CUT&RUN authors argue - That it has less background.**

Clearly, the CUT&Run data scores better by QC compared to the CTCF ChIP presented in this paper. But how about comparing to the Pugacheva CTCF ChIP datasets?

.. code-block:: bash

    plotFingerprint -b \
      ./bam/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoC_Abs.GRCh38.bam \
      ./bam/NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoN_Abs_Rep1.GRCh38.bam \
      ./bam/NBIS_Pugacheva2020_K562_ChIP_CTCF_RabbitC_Abs.GRCh38.bam \
      ./bam/NBIS_Pugacheva2020_K562_ChIP_IGG_Abs_Control.GRCh38.bam \
      ./bam/NBIS_Pugacheva2020_K562_ChIP_Mix_of_CTCF_PolyC_Abs.GRCh38.bam \
      -o ./fingerprint_Pugacheva2020.pdf \
      --labels Pugacheva2020_ChIP_MonoC Pugacheva2020_ChIP_MonoN Pugacheva2020_ChIP_RabbitC Pugacheva2020_ChIP_Control Pugacheva2020_ChIP_Mix \
      -p 4

.. image:: Figures/04_fingerprint2.png
	:target: Figures/04_fingerprint2.png
	:alt:

**Q: How does Pugacheva CTCF ChIP measure up with the CUT&Run? Given that the Skene CTCF ChIP was done under native conditions (no crosslinker) and Pugacheva CTCF ChIP was formaldehyde crosslinked, what do you think could be the problem with the native ChIP?**


Peak calling
============

Peaks were called with MACS2 using standard parameters.

.. attention::
    This is a step that could be better fine-tuned to specific experimental settings. 

Again, it is usually a good idea to visually inspect the tracks, so you can have a feeling on whether the peaks were correctly called and how the samples look like.

**Local**: You can visualize the peaks and bigwig files you downloaded before to your local computer.

.. image:: Figures/05_IGV_peaks.png
	:target: Figures/05_IGV_peaks.png
	:alt:

**Q: Given the QC you did above, does the peak calling confirm the quality differences amongst the samples? Does a higher signal/noise ratio allow to identify more peaks? Are peaks more confidently called?**

Number of peaks per sample
--------------------------

A simple ``wc -l peaks/*.narrowPeak`` count allows you to quickly check how many peaks you got:

.. code-block:: bash

    94019 NBIS_Pugacheva2020_K562_ChIP_Mix_of_CTCF_PolyC_Abs_peaks.narrowPeak
    65837 NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoC_Abs_peaks.narrowPeak
    50420 NBIS_Pugacheva2020_K562_ChIP_CTCF_RabbitC_Abs_peaks.narrowPeak
    40077 NBIS_Pugacheva2020_K562_ChIP_CTCF_MonoN_Abs_Rep1_peaks.narrowPeak
    828   NBIS_Pugacheva2020_K562_ChIP_IGG_Abs_Control_peaks.narrowPeak
    
    67164 NBIS_Skene2017_K562_CTCF_CnR_9m_peaks.narrowPeak
    52647 NBIS_Skene2017_K562_CTCF_CnR_3m_peaks.narrowPeak
    34678 NBIS_Skene2017_K562_CTCF_CnR_452_peaks.narrowPeak
    3782  NBIS_Skene2017_K562_CTCF_ChIP_loMN_peaks.narrowPeak
    3348  NBIS_Skene2017_K562_CTCF_ChIP_medMN_peaks.narrowPeak
    2897  NBIS_Skene2017_K562_CTCF_CnR_15s_peaks.narrowPeak
    2413  NBIS_Skene2017_K562_CTCF_CnR_5s_peaks.narrowPeak
    1611  NBIS_Skene2017_K562_CTCF_ChIP_hiMN_peaks.narrowPeak


Peaks overlap using intervene
------------------------------

``intervene`` is an easy-to-use tool to look for overlaps between BED files. It relies on ``bedtools``, but it saves some work when looking at different sets of files. You can install it using ``pip`` as they `explain <https://intervene.readthedocs.io/en/latest/install.html>`_.

.. attention::
    There is no ``intervene`` module on Uppmax. If you want to run it there, you can activate a conda environment that is precomputed: ``conda activate /sw/courses/epigenomics/quantitative_chip_simon/condaenv/intervene``. Otherwise you can download the peaks files to your local computer and install intervene there, if you prefer.

**Uppmax**: You can generate venn diagrams (pairwise or more). For example, we may want to look at how much two of the CTCF ChIP peaks from Pugacheva 2020 agree:

.. code-block::
    
    cd cnr_chip
    mkdir peaks
    cp /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/peaks/*.narrowPeak peaks/

    conda activate /sw/courses/epigenomics/quantitative_chip_simon/condaenv/intervene

    intervene venn --in ./peaks/NBIS_Pugacheva2020_K562_ChIP_CTCF_Mono*.narrowPeak

This will output a ``intervene_results`` folder with a pdf file:

.. image:: Figures/06_venn_1.png
	:target: Figures/06_venn_1.png
	:alt:

As you can see there is a lot of overlap but MonoC dataset called much more peaks than MonoN.

Note that this can also be refined by setting overlap thresholds. See intervene documentation for the possibilities. The default behavior is to count any overlap as overlap.

A way of looking broadly at a set of BED files an the overlap between them is to do pairwise comparison:

.. code-block:: bash

    intervene pairwise --in ./peaks/*.narrowPeak -o pairwise_results

will generate a plot like this in ``./pairwise_results/``.

.. image:: Figures/07_venn_pairwise.png
	:target: Figures/07_venn_pairwise.png
	:alt:

Very different numbers of peaks! 

**Q: How do you think the number of peaks relates to the fingerprint? (the ones with the most accentuated fingerprint are the ones showing larger amounts of peaks, however the difference is not such for some of the CnR samples. Why can that be?**

This has actually been noted in the CnR paper. They say that too scarce background read coverage can throw off traditional peak callers, thus they develop their own peak caller for CnR data. 


**Q: Now that you have peaks, think about what you could do with the peak information. How to make sense of the peaks? Would you use the dataset with the most or the least peaks for downstream analysis?**

.. note::
    ``MACS`` doesn’t just give peaks, it also assigns a score. High confidence peaks have higher score. Small peaks have low score. It is not apparent from the analysis above, but it is quite likely that if you would pick the top 5000 scored peaks from each dataset, the overlap would be better.

Comparison between methods
==========================

As an example, we are going to look at one peak set representative for each method: CnR_45s, CTCF_ChIP_medMN from Skene 2017 and ChIP_MonoC from Pugacheva 2020.

We can look at the overlap between them in a venn diagram:

.. image:: Figures/08_three_venn.png
	:target: Figures/08_three_venn.png
	:alt:


There are quite some differences between approaches, and one of the peak sets is very small compared to the others.

.. note::
    These venn diagrams are not size proportional. Automatically drawing size-proportional set intersections is a complex problem. `eulerr` is an `R package <https://cran.r-project.org/web/packages/eulerr/vignettes/introduction.html>`_ that does a really nice job at approximating this. However it’s only for drawing, not for computing the actual intersections.


Heatmaps
=========

The fact that there are loci marked as peak in one dataset that do not appear in another does not mean that there is no signal there. It could be that a peak is not called due to lack of statistical power, the signal being weaker or a combination of other factors. Remember the example in the slides - some peaks don’t make the threshold, but still show enrichment:

.. image:: Figures/09_peak_example.png
	:target: Figures/09_peak_example.png
	:alt:

So another way to look at this is to plot the signal of every dataset in a given peakset, for instance. This was run on a subset of 20000 peaks from Skane2017_CnR_45s. 

.. image:: Figures/10_heatmaps_subsample.png
	:target: Figures/10_heatmaps_subsample.png
	:alt:


.. warning::
    These plots were generated in the past using ``seqplots``. Unfortunately, this package `is now deprecated <https://bioconductor.org/packages/release/bioc/html/seqplots.html>`_. You can still
    use older versions of R and this package, or just generate plots like this using another tool, like :code:`deepTools`. You can see more in their `documentation <https://deeptools.readthedocs.io/en/develop/>`_. And there are also examples in other tutorials, such as the MINUTE-ChIP data tutorial in this section.


The subset file this was created on can be copied from the `tmp` folder:

.. code-block:: bash
    
    scp <youruser>@rackham.uppmax.uu.se:/sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/tmp/Skene2017_CTCF_CnR_45s_small.bed .


Or from your Uppmax node:

.. code-block:: bash
    
    cp /sw/courses/epigenomics/quantitative_chip_simon/K562_CTCF_CnR/tmp/Skene2017_CTCF_CnR_45s_small.bed . 


**Q: how does the data compare? Note the difference in scale amongst the samples. What other features are different? Which dataset would be better to identify the CTCF binding motif?**

.. note::
    In Skene and Henikoff paper, the authors make a biological interpretation regarding the difference in peaks seen in native ChIP and CUT&Run. They call the peaks uniquely present in CUT&Run “indirect binding sites” because they infer that those peaks are not directly bound by CTCF. The tagging of sites potentially in intact nuclei works through space, thus tagging not only the loci that are bound by CTCF, but also those regions that are nearby in space (Hi-C contacts). So it is important to note that CUT&Run may report ‘indirect’ or ‘shadow’ peaks that do not represent bona fide binding sites. For CTCF in principle the distinction should be easy since the ‘shadow’ peaks should **not have CTCF binding motif**.

.. image:: Figures/11_cut_run.png
	:target: Figures/11_cut_run.png
	:alt:


References
===============

.. [1] Ong, CT., Corces, V. CTCF: an architectural protein bridging genome topology and function. Nat Rev Genet 15, 234–246 (2014).
.. [2] Skene PJ, Henikoff S. An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites. Elife 2017.
.. [3] Pugacheva, Elena M., et al. CTCF mediates chromatin looping via N-terminal domain-dependent cohesin retention. Proceedings of the National Academy of Sciences 117.4 (2020).