.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


============================================
Detection of broad peaks from ChIP-seq data
============================================


**Learning outcomes**

- be able to assess the quality of ChIP-seq data for factors with broad occupancy pattern
- be able to detect regions of enrichment for factors with broad occupancy pattern



.. Contents
.. =========

.. contents:: 
    :local:




Requirements
==============

* ``MACS 3.0.0a6``
* ``epic2 0.0.52``
* ``R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"``
* ``csaw 1.24.3`` and its dependencies


Bioconductor packages required for annotation:

* ``org.Hs.eg.db``
* ``TxDb.Hsapiens.UCSC.hg38.knownGene``


Please note that this lab consists of three parts: 

(i) calling broad peaks using ``MACS`` (Uppmax);

(ii) calling broad peaks using ``epic2`` (Uppmax); and 

(iii) finding enriched genomic windows using  ``csaw`` in ``R``  (Uppmax / local).


Notes on software
-------------------

We provide a conda environment to run ``epic2``. This package proved a bit tricky to install because of dependency incompatibilities. To find how this environment was constructed, please visit :doc:`Dependencies <../../dependencies>`.


Please note that this workflow has been tested using ``R 4.0.3`` and ``csaw 1.24.3`` on Uppmax only.


Instructions how to install **R and Bioconductor packages** (including dependencies for ``csaw``) can be found in instructions to previous labs, for example :doc:`csaw tutorial <../csaw/lab-csaw>`.



Data
=====

We will use ChIP-seq of **H3K79me2** from ENCODE. H3K79me2 is enriched at **active promoters** and linked to transcriptional activation as it tends to accumulate in **transcribed regions of active genes**. 

The samples are from GM23338 cell line (human induced pluripotent stem cell from skin fibroblasts) before and after stimulation with doxycycline to form bipolar neurons.

These data are paired-end with read length 76 bp (2x76) and consist of duplicates of H3K79me2 ChIP and matching input samples. Raw data was processed using ENCODE3 pipeline and mapped to ``GRCh38`` genome assembly. We have further processed the data to remove unwanted signal and subset to chromosomes 1 and 2, as in the tutorial :doc:`ChIPseq data processing <../chipseqProc/lab-chipseq-processing>`. Accession details are summarised in Table 1.



.. list-table:: Table 1. Samples and accessions in the data set used in this exercise.
   :widths: 40 25 25
   :header-rows: 1

   * - Sample
     - ENCODE Accession
     - Replicate
   * - GM23338, H3K79me2 ChIP
     - ENCFF263XBT
     - rep 1
   * - GM23338, H3K79me2 ChIP
     - ENCFF531WSU
     - rep 2
   * - GM23338, input
     - ENCFF992ZBS
     - rep 1
   * - GM23338, input
     - ENCFF237PRF
     - rep 2
   * - neuro GM23338, H3K79me2 ChIP
     - ENCFF395DAJ
     - rep 1
   * - neuro GM23338, H3K79me2 ChIP
     - ENCFF806KRA
     - rep 2
   * - neuro GM23338, input
     - ENCFF956GLJ
     - rep 1
   * - neuro GM23338, input
     - ENCFF687LIL
     - rep 2



We will detect occupancy regions from *one sample only*, and compare the results to other samples processed earlier.



:raw-html:`<br />`


Quality control
================

As always, one should start the analysis from assesment of data quality. This is already performed to save time, and the plots and metrics are included below.


Cross-correlation and related metrics
----------------------------------------

The files discussed in this section can be accessed at 
``/sw/courses/epigenomics/broad_peaks2021/results/qc``

These metrics have been developed with application to TF ChIP-seq in mind, and you can see that the results for broad domains are not as easy to interpret as for point-source factors. Below are cross correlation plots for the IP and input you are going to use for the exercise. 

Already from these plots alone it is evident that the data has some quality issues. At this point you should be able to identify them.


.. list-table:: Figure 1. Cross correlations in H3K79me2 ChIP  and input samples in GM23338-derived neuron cells (ENCODE).
   :widths: 40 40
   :header-rows: 1

   * - K79me2 ChIP (ENCFF395DAJ)
     - input (ENCFF956GLJ)
   * - .. image:: figures/ENCFF395DAJ-xcor.png
   			:width: 300px
     - .. image:: figures/ENCFF956GLJ-xcor.png
   			:width: 300px



The cross correlation profile of factors with broad occupancy patterns is not going to be as sharp as for TFs, and the values of NSC and RSC tend to be lower, which does not mean that the ChIP failed. In fact, the developers of the tool do not recommend using the same NSC / RSC values as quality cutoffs for broad marks. However, input samples should not display signs of enrichment, as is the case here.

:raw-html:`<br />`

Cumulative enrichment
----------------------

Another plot worth examining is cumulative enrichment (aka fingerprint from deepTools):

.. list-table:: Figure 2. Cumulative enrichment (bamFingerprint) in H3K79me2 ChIP and input samples in GM23338-derived neuron cells (ENCODE).
   :widths: 60
   :header-rows: 1

   * - all samples
   * - .. image:: figures/nGM23338_fingerprint.png
   			:width: 600px





You can see that even though the cross correlation metrics don't look great, a significant enrichment can be observed for the ChIP samples (ENCFF395DAJ, ENCFF806KRA), and not for the input samples.


:raw-html:`<br />`
:raw-html:`<br />`


Broad peak calling using MACS
===============================

MACS: Model-based Analysis for ChIP-Seq is one of the leading peak calling algorithms. It has been excellent for detection of point-source peaks. However, until the recent version 3, it somewhat underperformed when used for detection of broad signal. Fortunatley, version 3, which is still under active development and hasn't been officially released, seems to fix issues with calling broad peaks. We will use this new version in this tutorial.



You will call peaks using sample GM23338 neuro - H3K79me2 ChIP  ``ENCFF395DAJ`` and its matching input ``ENCFF956GLJ``.

Effective genome size for chr 1 and 2 in ``hg38`` is ``4.9e8``.



.. code-block:: bash
	
  mkdir -p results/macs3
  cd results/macs3

  ln -s /sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338/ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam
  ln -s /sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338/ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam

  module load bioinfo-tools #if needed
  module load MACS/3.0.0a6

  macs3 callpeak --broad \
  -t ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
  -c ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
  -f BAMPE  -g 04.9e8 --broad-cutoff 0.1 -n GM23338_rep1



The main difference here, in comparison to detecting narrow peaks, is using the options ``--broad --broad-cutoff 0.1``. With the option ``--broad`` on, MACS will try to composite broad regions in BED12 (gene-model-like format) by putting nearby highly enriched regions into a broad region with loose cutoff. The broad region is controlled by another cutoff through ``--broad-cutoff``. If ``-p`` is set, this is a p-value cutoff, otherwise, it's a q-value (FDR) cutoff.

Because we use PE data, there is no need to build a model to estimate fragment length (similar to cross correlation) necessary for extending the SE reads. We know precisely how long each fragment is because its both ends are sequenced and mapped to the reference.



.. If you would like to compare the results of two different methods of finding broad peaks, repeat this with another data set:

.. .. code-block:: bash

.. 	ln -s /sw/courses/epigenomics/broad_peaks/bam/SRR1536561.bwt.hg38_dm6.sorted.hg38.BLfilt.bam
.. 	ln -s /sw/courses/epigenomics/broad_peaks/bam/SRR1584493.bwt.hg38_dm6.sorted.hg38.BLfilt.bam

.. 	macs2 callpeak -t SRR1536561.bwt.hg38_dm6.sorted.hg38.BLfilt.bam -c SRR1584493.bwt.hg38_dm6.sorted.hg38.BLfilt.bam -n 100_R1 --outdir 100_R1 -f BAM --gsize 3.0e9 -q 0.1 --nomodel --extsize 180 --broad --broad-cutoff 0.1




You can now inspect the results in the output folder ``macs3``. The structure is alike the output for calling narrow peaks. The file ``*.broadPeak`` is in ``BED6+3`` format which is similar to ``narrowPeak`` file used for point-source factors, except for missing the 10th column for annotating peak summits. Look at `MACS repository homepage <https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md>`_ for details.

The meaning of columns in ``NAME_peaks.xls`` files:

chr
  chromosome name

start
  start position of peak

end
   end position of peak

length
  length of peak region

pileup
  pileup height at peak summit

-log10(pvalue)
  -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)

fold_enrichment
  fold enrichment for this peak summit against random Poisson distribution with local lambda

-log10(qvalue)
  -log10(qvalue) at peak summit

name
  peak id
 
    

Let's take a look at another format of the output ``broadPeak``. It is compatible with major genome browsers (IGV, UCSC Genome Browser) and easier to work with because it does not contain a long header.

This is an example::

  head neuroGM23338_macs3_rep1_peaks.broadPeak

  chr1  777491  778262  neuroGM23338_macs3_rep1_peak_1  34  . 3.542 4.93525 3.48401
  chr1  779812  780867  neuroGM23338_macs3_rep1_peak_2  10  . 2.28884 2.27839 1.03252
  chr1  782000  784521  neuroGM23338_macs3_rep1_peak_3  17  . 2.6654  3.01765 1.70342
  chr1  820548  826643  neuroGM23338_macs3_rep1_peak_4  36  . 3.5486  5.10182 3.65624
  chr1  828271  830128  neuroGM23338_macs3_rep1_peak_5  34  . 3.4958  4.87316 3.42798
  chr1  831350  833671  neuroGM23338_macs3_rep1_peak_6  22  . 2.7518  3.55309 2.20097
  chr1  882552  890194  neuroGM23338_macs3_rep1_peak_7  34  . 3.21783 4.86863 3.43262
  chr1  925794  926897  neuroGM23338_macs3_rep1_peak_8  18  . 2.71963 3.12803 1.80546
  chr1  957085  959246  neuroGM23338_macs3_rep1_peak_9  60  . 4.54986 7.61848 6.03443
  chr1  999291  999914  neuroGM23338_macs3_rep1_peak_10 16  . 2.63811 2.95948 1.65064



The meaning of columns in ``NAME.broadPeak`` files:

    
chrom
  Name of the chromosome (or contig, scaffold, etc.).
chromStart
  The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd
  The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores were "0" when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
name
  Name given to a region (preferably unique). Use "." if no name is assigned.
score
  Indicates how dark the peak will be displayed in the browser (0-1000).
strand
  +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
signalValue
  Measurement of overall (usually, average) enrichment for the region.
pvalue
  Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
qValue
  Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.



How many peaks were identified in replicate 1?


.. code-block:: bash

  wc -l neuroGM23338_macs3_rep1_peaks.broadPeak
  6826 neuroGM23338_macs3_rep1_peaks.broadPeak



.. HINT::

	You can also copy the results from
	``/sw/courses/epigenomics/broad_peaks2021/results/macs3/neuroGM23338``

This is a preliminary peak list, and in case of broad domains, it often needs some processing or filtering.


Firstly, let's select the peaks reproducible in both replicates. 


.. admonition:: Select reproducible peaks (MACS).
   :class: dropdown, warning


   .. code-block:: bash

      #link the files if you are in a different directory
      ln -s /sw/courses/epigenomics/broad_peaks2021/results/macs3/neuroGM23338/neuroGM23338_macs3_rep1_peaks.broadPeak
      ln -s /sw/courses/epigenomics/broad_peaks2021/results/macs3/neuroGM23338/neuroGM23338_macs3_rep2_peaks.broadPeak

      #make bed
      cut -f 1-6 neuroGM23338_macs3_rep1_peaks.broadPeak >neuroGM23338_macs3_rep1_peaks.bed
      cut -f 1-6 neuroGM23338_macs3_rep2_peaks.broadPeak >neuroGM23338_macs3_rep2_peaks.bed


      #intersect bed files
      module load BEDTools/2.29.2
      bedtools intersect -a neuroGM23338_macs3_rep1_peaks.bed -b neuroGM23338_macs3_rep2_peaks.bed -f 0.50 -r > peaks_macs3_neuroGM23338.chr12.bed

      #how many peaks which overlap?
      wc -l peaks_macs3_neuroGM23338.chr12.bed 
      2679 peaks_macs3_neuroGM23338.chr12.bed






Visual inspection of the peaks (MACS)
======================================

You will use ``IGV`` for this step, and it is recommended that you run it locally on your own computer. Please load ``hg38`` reference genome.

Required files are:

* ChIP ``ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam`` and matching ``bai``
* input ``ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam`` and matching ``bai``
* signal domains ``neuroGM23338_macs3_rep1_peaks.broadPeak``
* reproducible signal domains ``peaks_macs3_neuroGM23338.chr12.bed``


.. HINT::

	You can access the bam and bai files from
	``/sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338``



You can look at some locations of interest. Peaks with low FDR (q value) or high fold enrichment may be worth checking out. Or check your favourite gene.

.. admonition:: Potentially interesting locations to view (MACS peaks).
   :class: dropdown, warning

   Let's sort the broadPeak file using the score column to find the peaks with the strongest signal

   .. code-block:: bash

    sort -k5,5rn neuroGM23338_macs3_rep1_peaks.broadPeak | head

    chr1 226062814 226073870 neuroGM23338_macs3_rep1_peak_3341 518 . 13.0292 54.709  51.8498
    chr1  234598782 234610959 neuroGM23338_macs3_rep1_peak_3515 513 . 12.139  54.2276 51.3993
    chr2  101698297 101748719 neuroGM23338_macs3_rep1_peak_5147 462 . 13.375  49.0116 46.2392
    chr2  47158830  47176361  neuroGM23338_macs3_rep1_peak_4479 449 . 12.5555 47.6423 44.96
    chr1  204403186 204412701 neuroGM23338_macs3_rep1_peak_2999 431 . 10.9654 45.8227 43.1784
    chr1  220527779 220538029 neuroGM23338_macs3_rep1_peak_3239 414 . 11.5724 44.1237 41.4922
    chr1  244049535 244060201 neuroGM23338_macs3_rep1_peak_3688 403 . 10.2749 42.9757 40.397
    chr2  54970096  55050304  neuroGM23338_macs3_rep1_peak_4577 401 . 11.3781 42.7092 40.1167
    chr2  5692693 5703228 neuroGM23338_macs3_rep1_peak_3837 399 . 10.6746 42.502  39.9488
    chr1  150568640 150579241 neuroGM23338_macs3_rep1_peak_2121 394 . 11.0116 42.0139 39.4717




Below you see IGV visualisations of the following regions (top two peaks and one of the bottom rank):

.. code-block:: bash

  chr1:226,055,687-226,080,997
  chr1:234,592,216-234,617,526
  chr1:777,176-783,503

IGV settings for this visualiation: Group alignments (by read strand); Colour alignments (by read strand); Squished.

Regions detected by ``MACS3`` are the topmost purple track, two bam files are ChIP and input (with their pileups calculated by IGV), and the bottom panel are gene models and, finally the regions reproducible between both replicates in green.

Please note the length of these detected domains. 



.. list-table:: Figure 3. Results of peak calling in H3K79me2 ChIP-seq in GM23338-derived neuron cells (ENCODE). Tracks from the top: peaks in rep1, ChIP, input, gene models, reproducible peaks
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/macs3-rep1-peak1.png
   			:width: 600px
   * - .. image:: figures/macs3-rep1-peak2.png
        :width: 600px
   * - .. image:: figures/macs3-rep1-peak3.png
        :width: 600px




.. Postprocessing of peak candidates
.. ====================================

.. Please note that this step is only an example, as **any postprocessing of peak calling results is highly project specific**.

.. Normally, you would work with replicated data. As in the case of TFs earlier, it is recommended to continue working with peaks reproducible between replicates.

.. The peak candidate lists can and should be further filtered, based on fold enrichment and pileup value, to remove peaks which could have a high fold enrichment but low signal, as these are likely non-informative. Any filtering, however has to be performed having in mind the biological characteristics of the signal.

.. You can merge peaks which are close to one another using `bedtools <https://bedtools.readthedocs.io/en/latest/>`_. You will control the distance of features to be merged using option ``-d``. Here we arbitrarily choose 1 kb.


.. .. code-block:: bash

.. 	cp 50_R1_peaks.broadPeak 50_r1.bed

.. 	module load bioinfo-tools
.. 	module load BEDTools/2.27.1

.. 	bedtools merge -d 1000 -i 50_r1.bed > 50_r1.merged.bed

.. 	#how many peaks?
.. 	wc -l *bed
.. 	27699 50_r1.bed
..   	11732 50_r1.merged.bed

:raw-html:`<br />`


Broad peak calling using epic2
===============================

epic2 is an ultraperformant reimplementation of SICER, an algorithm developed especially for detection of broad marks. It focuses on speed, low memory overhead and ease of use. It also contains a reimplementation of the SICER-df scripts for differential enrichment and a script to create many kinds of bigwigs from your data. In this tutorial we will use it to detect domains in the same data as we used earlier for MACS. At the end we will compare the results.

Here again we use a prepared conda environment. Newer versions of ``Pysam`` seem to throw errors when used with ``epic2``. For details please consult :doc:`Dependencies <../../dependencies>`.



.. code-block:: bash
  
  mkdir ../epic2
  cd ../epic2

  ln -s /sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338/ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam
  ln -s /sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338/ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam

  conda activate /sw/courses/epigenomics/software/conda/epic_2b

  epic2 --treatment ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
   --control ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
    -fdr 0.05 --effective-genome-fraction 0.95 \
    --chromsizes /sw/courses/epigenomics/broad_peaks2/annot/hg38_chr12.chromsizes \
    --guess-bampe --output neuroGM23338.rep1.epic2


The result looks like this::

  head neuroGM23338.rep1.epic2
  #Chromosome Start End PValue  Score Strand  ChIPCount InputCount  FDR log2FoldChange
  chr1  777400  778199  1.525521715195486e-17 302.351431362403  . 28  4 3.6469084588658344e-17  3.02351431362403
  chr1  821600  822599  9.82375064925635e-15  309.0628509482567 . 22  3 2.149649196967778e-14 3.090628509482567
  chr1  823400  826599  1.912606813336636e-50 258.96177870938703  . 114 22  7.16207410279741e-50  2.58961778709387
  chr1  828200  830399  5.824989404621774e-19 201.13395996779278  . 59  17  1.4498521203360967e-18  2.0113395996779277
  chr1  831400  833599  2.388427396703667e-14 151.73289262869903  . 69  28  5.17334134771682e-14  1.5173289262869905
  chr1  880600  885799  2.569723734705377e-37 153.80874864537878  . 195 78  8.520333236155976e-37 1.538087486453788
  chr1  886600  890199  1.481460204837024e-22 165.05622157122005  . 100 37  3.9835062883707675e-22  1.6505622157122006
  chr1  925800  926999  8.491192455372113e-17 280.11218922875815  . 30  5 1.9845609170824598e-16  2.801121892287582
  chr1  957000  959199  4.1356938759907084e-64  337.1437617044337 . 98  11  1.688095302270476e-63 3.3714376170443368


The meaning of the columns:

PValue
  Poisson-computed PValue based on the number of ChIP count vs. library-size normalized Input count in the region
Score
  Log2FC * 100 (capped at 1000). Regions with a larger relative ChIP vs. Input count will show as darker in the UCSC genome browser
ChIPCount
  The number of ChIP counts in the region (also including counts from windows with a count below the cutoff)
InputCount
  The number of Input counts in the region
FDR
  Benjamini-Hochberg correction of the p-values
log2FoldChange
  Log2 of the region ChIP count vs. the library-size corrected region Input count



How many domains were found? (the first line is a header)

.. code-block:: bash

  wc -l neuroGM23338.rep1.epic2
  5242 neuroGM23338.rep1.epic2




How many domains reproducible between replicates?


.. admonition:: Select reproducible peaks (epic2).
   :class: dropdown, warning


   .. code-block:: bash

      #link the files if you are in a different directory
      ln -s /sw/courses/epigenomics/broad_peaks2021/results/epic2/neuroGM23338/neuroGM23338.rep1.epic2
      ln -s /sw/courses/epigenomics/broad_peaks2021/results/epic2/neuroGM23338/neuroGM23338.rep2.epic2

      #make bed
      cut -f 1-3 neuroGM23338.rep1.epic2 >neuroGM23338_epic2_rep1_peaks.bed
      cut -f 1-3 neuroGM23338.rep2.epic2 >neuroGM23338_epic2_rep2_peaks.bed

      #intersect bed files
      module load bioinfo-tools #if necessary
      module load BEDTools/2.29.2
      bedtools intersect -a neuroGM23338_epic2_rep1_peaks.bed -b neuroGM23338_epic2_rep2_peaks.bed -f 0.50 -r > peaks_epic2_neuroGM23338.chr12.bed

      #how many peaks which overlap?
      wc -l peaks_epic2_neuroGM23338.chr12.bed
      2692 peaks_epic2_neuroGM23338.chr12.bed


How about the overlap between different methods?


.. admonition:: Compare MACS3 and epic2.
   :class: dropdown, warning

   (please make sure the relative path to macs3 results is correct in the command below)

   .. code-block:: bash

      #intersect bed files
      module load bioinfo-tools #if necessary
      module load BEDTools/2.29.2
      bedtools intersect -a peaks_epic2_neuroGM23338.chr12.bed -b ../macs3/peaks_macs3_neuroGM23338.chr12.bed \
      -f 0.50 -r > peaks_epic2macs3_neuroGM23338.chr12.bed
      
      #how many peaks which overlap?
      wc -l peaks_epic2macs3_neuroGM23338.chr12.bed
      1629 peaks_epic2macs3_neuroGM23338.chr12.bed


You can visualise the peaks as for MACS. Below are some of the locations as before, with peaks detected by both epic2 and MACS marked in orange.


.. list-table:: Figure 4. Results of peak calling in H3K79me2 ChIP-seq in GM23338-derived neuron cells (ENCODE). Comparison of MACS3 and epic2. Tracks from the top: peaks in rep1, ChIP, input, gene models, reproducible peaks (MACS3), peaks detected by epic2 and MACS3 (orange)
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/epic2-macs3-peak1.png
        :width: 600px
   * - .. image:: figures/epic2-macs3-peak2.png
        :width: 600px
   * - .. image:: figures/epic2-macs3-peak3.png
        :width: 600px



.. admonition:: Locations plotted.
   :class: dropdown, warning

   .. code-block:: bash

    chr1:226,055,687-226,080,997
    chr1:237,310,726-237,323,381
    chr2:202,204,546-202,229,168


:raw-html:`<br />`


Alternative approach: window-based enrichment analysis (csaw)
===============================================================

This workflow is similar to the one using ``csaw`` designed for TF peaks. The differences pertain to analysis of signal from diffuse marks. Please check the :doc:`csaw tutorial <../csaw/lab-csaw>` for setup and more detailed comments on each step.

.. You will use data from the same dataset, however, the files were processed in a different manner: the alignments were not filtered to remove duplictate reads nor the reads mapping to the ENCODE blacklisted regions. To reduce the computational burden, the bam files were subset to contain alignments to ``chr1``.

.. NOTE::
  
  This exercise was tested on Rackham using pre-installed R libraries. Local installation of recommended R packages may require additional software dependecies.


.. Requirements Local
.. ----------------------

.. * ``csaw``
.. * ``edgeR``

.. R packages required for annotation:

.. * ``org.Hs.eg.db``
.. * ``TxDb.Hsapiens.UCSC.hg38.knownGene``

.. Recommended:

.. * R-Studio to work in



.. **Getting the data**


.. First, you need to copy the necessary files to your laptop:

.. .. code-block:: bash

.. 	cd /desired/location

.. 	scp <USERNAME>@rackham.uppmax.uu.se:/sw/courses/epigenomics/broad_peaks/broad_peaks_bam.tar.gz .

.. 	#type your password at the prompt

.. 	tar zdvf broad_peaks_bam.tar.gz




Requirements Remote (Uppmax)
--------------------------------

The software is configured.

To prepare the files, assuming you are in ``~/broad_peaks/results``:

.. code-block:: bash
  
   mkdir csaw
   cd csaw

   mkdir bam
   ln -s  /sw/courses/epigenomics/broad_peaks2021/data_sub_preproc/neuron_GM23338/* bam



.. Remote:

.. code-block:: bash

    module load conda/latest
    conda activate /sw/courses/epigenomics/software/conda/v8
    R


The remaining part of the exercise is performed in ``R``.


.. Remote:


.. code-block:: R

  # provide the tutorial specific path to R libraries
  assign(".lib.loc", "/sw/courses/epigenomics/software/R", envir = environment(.libPaths))

  # verify that the tutorial-specific R library path is added
  .libPaths()
  [1] "/sw/courses/epigenomics/software/R"


Sort out the working directory and file paths:

.. code-block:: R

	setwd("/path/to/workdir")

	dir.data = "/path/to/desired/location/bam"

	#for example when in broad_peaks/csaw
	dir.data = "./bam"	

	k79_1=file.path(dir.data,"ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam")
	input_1=file.path(dir.data,"ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam")
	k79_2=file.path(dir.data,"ENCFF806KRA.chr12.MAPQ30.blcklst.rh.sorted.bam")
	input_2=file.path(dir.data,"ENCFF687LIL.chr12.MAPQ30.blcklst.rh.sorted.bam")

	bam.files <- c(k79_1,k79_2,input_1,input_2)



.. HINT:: Setting the paths in R

  To find the path to your current location type ``pwd`` in the terminal. You can use this path in R like this:

  .. code-block:: bash

    setwd("/path/to/where_you_are")

  All the paths will be then relative to ``/path/to/where_you_are``.

  You can also find it directly from R using ``getwd``::

    > getwd()
    [1] "/crex/course_data/epigenomics/broad_peaks2021/results/csaw"



Read in the data:

.. code-block:: R

	frag.len=200

	library(csaw)

	data <- windowCounts(bam.files, ext=frag.len, width=100) 



.. admonition:: data
   :class: dropdown, warning


   .. code-block:: R

     > data
      class: RangedSummarizedExperiment 
      dim: 8808414 4 
      metadata(6): spacing width ... param final.ext
      assays(1): counts
      rownames: NULL
      rowData names(0):
      colnames: NULL
      colData names(4): bam.files totals ext rlen





You will identify the enrichment windows by performing a differential occupancy analysis between ChIP and input samples.

Information on the contrast to test:

.. code-block:: R

	grouping <- factor(c('chip', 'chip', 'input', 'input'))
	design <- model.matrix(~0 + grouping)
	colnames(design) <- levels(grouping)
	library(edgeR)
	contrast <- makeContrasts(chip - input, levels=design)


.. admonition:: contrast
   :class: dropdown, warning


   .. code-block:: R

      > contrast
       Contrasts
      Levels  chip - input
        chip             1
        input           -1




Next, you need to filter out uninformative windows with low signal prior to further analysis. Selection of appropriate filtering strategy and cutoff is key to a successful detection of differential occupancy events, and is data dependent. Filtering is valid so long as it is independent of the test statistic under the null hypothesis.
One possible approach involves choosing a filter threshold based on the fold change over
the level of non-specific enrichment (background). The degree of background enrichment is estimated
by counting reads in large bins across the genome.

The function ``filterWindowsGlobal`` returns the increase in the abundance of
each window over the global background. 
Windows are filtered by setting some minimum threshold on this increase. Here, a **fold change of 3** is necessary for a window to be considered as containing a binding site. 

In this example, you estimate the global background using ChIP samples only. You can do it using the entire dataset including inputs of course.

.. code-block:: R

	bam.files_chip <- c(k79_1,k79_2)

	bin.size <- 2000L
	binned.ip <- windowCounts(bam.files_chip, bin=TRUE, width=bin.size, ext=frag.len)
	data.ip=data[,1:2]
	filter.stat <- filterWindowsGlobal(data.ip, background=binned.ip)

	keep <- filter.stat$filter > log2(3)
	data.filt <- data[keep,]


To examine how many windows passed the filtering:

.. code-block:: R

	summary(keep)
  	 Mode   FALSE    TRUE 
  logical 7892070  916344 

To normalise the data for different library sizes you need to calculate normalisation factors based on large bins:

.. code-block:: R

	binned <- windowCounts(bam.files, bin=TRUE, width=10000)
	data.filt <- normFactors(binned, se.out=data.filt)

	data.filt$norm.factors
	## [1] 0.6100099 0.6649707 1.5649880 1.5752503




Detection of DB windows (in our case, the occupancy sites, as we test for differences in ChIP vs. input):

.. code-block:: R

	data.filt.calc <- asDGEList(data.filt)
	data.filt.calc <- estimateDisp(data.filt.calc, design)
	fit <- glmQLFit(data.filt.calc, design, robust=TRUE)
	results <- glmQLFTest(fit, contrast=contrast)



You can inspect the raw results:

.. code-block:: R

	head(results$table)
       logFC     logCPM        F       PValue
  1 2.814162 -0.2316357 15.13205 1.056751e-03
  2 3.221807 -0.3432297 14.69448 1.199631e-03
  3 3.588308 -0.4335388 17.64303 5.279200e-04
  4 3.790889 -0.1623742 27.22414 5.636337e-05
  5 4.372525  0.5046133 54.82505 6.826942e-07
  6 4.075533  0.5386263 49.87994 1.305226e-06


The following steps will calculate the FDR for each peak, merge peaks within 1 kb and calculate the FDR for resulting composite peaks.

.. code-block:: R

	merged <- mergeWindows(rowRanges(data.filt), tol=1000L)
	table.combined <- combineTests(merged$id, results$table)


Short inspection of the results:

.. code-block:: R

	head(table.combined)

  DataFrame with 6 rows and 8 columns
    num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
    <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
  1        17           17              0 1.20248e-06 1.20918e-05          up
  2        26           26              0 1.19007e-04 6.57639e-04          up
  3         3            3              0 1.56689e-03 3.68082e-03          up
  4         3            3              0 1.36304e-03 3.39817e-03          up
  5        78           78              0 5.10443e-07 5.61850e-06          up
  6        85           85              0 2.42116e-05 1.71896e-04          up
     rep.test rep.logFC
    <integer> <numeric>
  1        13   4.14334
  2        29   3.91717
  3        44   2.32107
  4        49   2.94867
  5       114   4.43809
  6       139   3.48985


How many regions are up (i.e. enriched in chip compared to input)?

.. code-block:: R

	is.sig.region <- table.combined$FDR <= 0.1
	table(table.combined$direction[is.sig.region])

     up 
  10794 


Does this make sense? How does it compare to results obtained from MACS and epic2 runs?

You can now annotate the results as in the csaw TF exercise:

.. code-block:: R

	library(org.Hs.eg.db)
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)

	anno <- detailRanges(merged$region, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene,
	orgdb=org.Hs.eg.db, promoter=c(3000, 1000), dist=5000)

	merged$region$overlap <- anno$overlap
	merged$region$left <- anno$left
	merged$region$right <- anno$right

	all.results <- data.frame(as.data.frame(merged$region)[,1:3], table.combined, anno)

	sig=all.results[all.results$FDR<0.05,]
	all.results <- all.results[order(all.results$PValue),]

	head(all.results)

	filename="k79me2_neuroGM_csaw.txt"
	write.table(all.results,filename,sep="\t",quote=FALSE,row.names=FALSE)

Let's inspect the results:

.. code-block:: R

  head(all.results)
      seqnames  start    end num.tests num.up.logFC num.down.logFC       PValue
    1     chr1 777251 778200        17           17              0 1.202476e-06
    2     chr1 779251 784600        26           26              0 1.190072e-04
    3     chr1 811251 811450         3            3              0 1.566886e-03
    4     chr1 816101 816300         3            3              0 1.363039e-03
    5     chr1 820551 826950        78           78              0 5.104427e-07
    6     chr1 828201 833800        85           85              0 2.421162e-05
               FDR direction rep.test rep.logFC                         overlap
    1 1.209181e-05        up       13  4.143336 100133331:-:PI,LOC100288069:-:P
    2 6.576393e-04        up       29  3.917169                   100133331:-:P
    3 3.680820e-03        up       44  2.321065                                
    4 3.398175e-03        up       49  2.948666                      FAM87B:+:P
    5 5.618504e-06        up      114  4.438095   LINC01128:+:PE,LINC00115:-:PE
    6 1.718961e-04        up      139  3.489853    LINC01128:+:PE,LINC00115:-:P
                                      left           right
    1 100133331:-:2971,LOC100288069:-:2971  100133331:-:84
    2  100133331:-:625,LOC100288069:-:4971                
    3                                                     
    4                                        FAM87B:+:1071
    5                         FAM87B:+:714 LINC01128:+:648
    6       LINC01128:+:74,LINC00115:-:679                


To compare with peaks detected by MACS it is convenient to save the results in ``BED`` format:

.. code-block:: R

	sig.up=sig[sig$direction=="up",]

	starts=sig.up[,2]-1

	sig.up[,2]=starts

	sig_bed=sig.up[,c(1,2,3)]

	filename="k79me2_neuroGM_csaw_peaks.bed"
	write.table(sig_bed,filename,sep="\t",col.names=FALSE,quote=FALSE,row.names=FALSE)

.. nrow(sig_bed)
.. 10711

You can now load the ``bed`` file to ``IGV`` along with the appropriate ``broad.Peak`` file and zoom in to your favourite location on chromosomes 1 and 2.


Below is the IGV snapshot of top peak, this time with csaw peaks added in light blue.


.. list-table:: Figure 5. Results of broad peak calling in H3K79me2 ChIP-seq in GM23338-derived neuron cells (ENCODE). Comparison of MACS3, epic2 and csaw. Tracks from the top: peaks in rep1, ChIP, input, gene models, reproducible peaks (MACS3), peaks detected by epic2 and MACS3 (orange), peaks deteced by csaw (light blue).
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/csaw_peak1.png
        :width: 600px


As you can see the regions with strong signal (high enrichment in ChIP over input) are detected by all methods tested. What about the sites with weak signal?

In this tutorial we have worked with good quality data which was sequenced to a recommended depth. All three methods tested in this tutorial perform well is such scenario. However, their preformace deteriorates with decreasing sequening depth (less data to rely on) and decreasing quality of the sample preparation (more noise).


.. ----

.. Written by: Agata Smialowska
