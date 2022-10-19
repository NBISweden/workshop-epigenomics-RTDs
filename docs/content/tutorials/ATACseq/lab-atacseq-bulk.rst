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


:raw-html:`<br />`


.. contents:: 
    :local:


This tutorial is a continuation of :doc:`Data preprocessing <../data-preproc/data-preproc>`, :doc:`General QC <../data-preproc/data-qc1>`, and :doc:`ATACseq specifc QC <../data-preproc/data-qc-atac>`. 

:raw-html:`<br />`


Data
======


We will use the same data as before: **ATAC-seq** libraries prepared to analyse chromatin accessibility status in natural killer (NK) cells without and with stimulation (in duplicates) from the `ENCODE <www.encodeproject.org>`_ project.

Natural killer (NK) cells are innate immune cells that show strong cytolytic function against physiologically stressed cells such as tumor cells and virus-infected cells. NK cells express several activating and inhibitory receptors that recognize the altered expression of proteins on target cells and control the cytolytic function. To read more about NK cells please refer to `Paul and Lal <https://doi.org/10.3389/fimmu.2017.01124>`_ . The interleukin cocktail used to stimulate NK cells induces proliferation and activation (`Lauwerys et al <https://doi.org/10.1006/cyto.1999.0501>`_ ).

ENCODE sample accession numbers are listed in Table 1.


.. list-table:: Table 1. ENCODE accession numbers for data set used in this tutorial.
   :widths: 10 25 25 50
   :header-rows: 1

   * - No
     - Accession
     - Cell type
     - Description
   * - 1
     - ENCFF398QLV
     - Homo sapiens natural killer cell female adult
     - untreated
   * - 2
     - ENCFF363HBZ
     - Homo sapiens natural killer cell female adult
     - untreated
   * - 3
     - ENCFF045OAB
     - Homo sapiens natural killer cell female adult
     - Interleukin-18, Interleukin-12-alpha, Interleukin-12-beta, Interleukin-15
   * - 4
     - ENCFF828ZPN
     - Homo sapiens natural killer cell female adult
     - Interleukin-18, Interleukin-12-alpha, Interleukin-12-beta, Interleukin-15


We have processed the data, starting from reads aligned to **hg38** reference assembly using **bowtie2**. The alignments were obtained from ENCODE in *bam* format and further processed:

* alignments were subset to include chromosome 14 and 1% of reads mapped to chromosomes 1 to 6 and chrM.

This allows you to see a realistic coverage of one selected chromosome and perform data analysis while allowing shorter computing times. Non-subset ATAC-seq data contains 100 - 200 M PE reads, too many to conveniently process during a workshop.

In this workshop, we have filtered and quality-controlled the data (parts :doc:`Data preprocessing <../data-preproc/data-preproc>`, :doc:`General QC <../data-preproc/data-qc1>`, and :doc:`ATACseq specifc QC <../data-preproc/data-qc-atac>`).




Peak Calling
=================

To find regions corresponding to potential open chromatin, we want to identify ATAC-seq "peaks" where reads have piled up to a greater extent than the background read coverage.

The tools which are currently used are `Genrich <https://github.com/jsh58/Genrich>`_ and `MACS2 <https://github.com/taoliu/MACS>`_.

* **Genrich** has a mode dedicated to ATAC-Seq; however, Generich is still not published;

* **MACS2** is widely used so lots of help is available online; however it is designed for ChIP-seq rather than ATAC-seq (MACS3 has more ATAC-seq oriented features, but still lacks a stable release);

The differences between these two peak callers in the context of ATAC-seq data are discussed `here <https://informatics.fas.harvard.edu/atac-seq-guidelines.html#peak>`_.



Shifting Alignments
-----------------------

We have already discussed (and performed) this step in the :doc:`ATACseq specifc QC <../data-preproc/data-qc-atac>` tutorial. Briefly, the alignments are shifted to account for the duplication created as a result of DNA repair after Tn5-introduced DNA nicks.


When Tn5 cuts an accessible chromatin locus it inserts adapters separated by 9bp, see Figure 2. This means that to have the read start site reflect the centre of where Tn5 bound, the reads on the **positive strand** should be **shifted 4 bp to the right** and reads on the **negative strand** should be **shifted 5 bp to the left** as in Buenrostro et al. 2013. 


.. list-table:: Figure 2. Nextera Library Construction.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/NexteraLibraryConstruction.jpg
   			:width: 400px



**To shift or not to shift?** It, as always, depends on the downstream application of the alignments.

If we use the ATAC-seq peaks for **differential accessibility**, and detect the peaks in the "broad" mode, then shifting does not play any role: the peaks are hundreds of bps long, reads are summarised to these peaks allowing a partial overlap, so 9 basepairs of difference has a neglibile effect. 

However, when we plan to use the data for any **nucleosome-centric analysis** (positioning at TSS or TF footprinting), shifting the reads allows to center the signal in peaks around the nucleosomes and not directly on the nucleosome.

If we only assess the coverage of the start sites of the reads, the data would be too sparse and it would be impossible to call peaks. Thus, we will extend the start sites of the reads by 100bp (50 bp in each direction) to assess coverage. This is performed automatically by Genrich, and using command line options ``extsize`` and ``shift`` in MACS2.



Genrich
---------

We start this tutorial in directory ``atacseq/analysis/``:

.. code-block:: bash
	
	mkdir peaks
	cd peaks


We need to link necessary files first.

.. code-block:: bash

	mkdir genrich
	cd genrich

	# we link the pre-processed bam file
	ln -s ../../../processedData/ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam ENCFF045OAB.chr14.proc.bam
	ln -s ../../../processedData/ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam.bai ENCFF045OAB.chr14.proc.bam.bai


.. Hint ::

	If you got lost on the way, you can link files preprocessed by us:

	.. code-block:: bash

		ln -s ../../../data_proc/ENCFF045OAB.chr14.proc.bam ENCFF045OAB.chr14.proc.bam
		ln -s ../../../data_proc/ENCFF045OAB.chr14.proc.bam.bai ENCFF045OAB.chr14.proc.bam.bai


Genrich requires bam files to be name-sorted rather than the default coordinate-sorted. Also, we remove all reference sequences other than chr14 from the header, as this is where our data is subset to. Genrich uses the reference sequence length from bam header in its statistical model, so retaining the original bam header would impair peak calling statistics.


.. code-block:: bash

	# in case not already loaded
	module load bioinfo-tools
	module load samtools/1.8


	#subset bam and change header
	samtools view -h ENCFF045OAB.chr14.proc.bam chr14 | grep -P "@HD|@PG|chr14" | samtools view -Shbo ENCFF045OAB.chr14.proc_rh.bam


	# sort the bam file by read name (required by generich)
	samtools sort -n -o ENCFF045OAB.chr14.proc_rh.nsort.bam -T sort.tmp  ENCFF045OAB.chr14.proc_rh.bam



Genrich can apply the read shifts when ATAC-seq mode ``-j`` is selected. We detect peaks by:

.. code-block:: bash

	/sw/courses/epigenomics/ATACseq_bulk/software/Genrich/Genrich -j -t ENCFF045OAB.chr14.proc_rh.nsort.bam  -o ENCFF045OAB.chr14.genrich.narrowPeak


The output file produced by Genrich is in `ENCODE narrowPeak format <https://genome.ucsc.edu/FAQ/FAQformat.html#format12>`_, listing the genomic coordinates of each peak called and various statistics.


.. code-block:: bash
	
	chr start end name score strand signalValue pValue qValue peak

	signalValue - Measurement of overall (usually, average) enrichment for the region.
	pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
	qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.

How many peaks were detected?

.. code-block:: bash
	
	wc -l ENCFF045OAB.chr14.genrich.narrowPeak
	1027 ENCFF045OAB.chr14.genrich.narrowPeak


.. admonition:: ENCFF045OAB.chr14.genrich.narrowPeak
   :class: dropdown, warning


   .. code-block:: bash

		head ENCFF045OAB.chr14.genrich.narrowPeak



MACS2
-----

MACS2 is widely used for peak calling in ATAC-seq, as evidenced by literature and many data processing pipelines. Several different peak calling protocols / commands have been encountered in various sources:

1. macs narrow peak (default for ``callpeak``), ``--nomodel``, shifted reads using PE reads as SE (bed file);

2. macs narrow peak, ``--nomodel``, shifted reads, using PE reads (bedpe file);

3. macs narrow peak, ``--nomodel``, using PE reads (bedpe file); ``shift`` and ``extsize`` similar to Genrich;

4. macs broad peak, using PE reads (bampe file) as in nf-core;

5. macs narrow peak, ``--nomodel``,  using PE reads as SE (bam file) as in Encode, ``shift`` and ``extsize`` similar to Genrich;

6. macs narrow peak, unshifted reads in bampe file

7. macs broad peak, BAM (PE reads used as SE reads)


The peaks obtained by these above commands have been compared to peaks detected by Genrich, and examples are shown on Figures 3 - 6.


Figure 3 depicts large genomic region. In general Genrich detects less peaks (shown in green) than MACS2 (navy). MACS2 commands 1, 2 and 3 result in many peaks in regions where Genrich detects none. MACS commands 4, 5, 6 and 7 produce less peaks which are somewhat similar to the result of Genrich. If we zoom in, we can see that commands 1, 2 and 3 detect spurious peaks which do not have strong evidence in alignment pipeups (Figures 4 to 6).


.. list-table:: Figure 3. Comparison of peaks detected by different algorithms. Overview of a large genomic region.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/igv_5.png
   			:width: 400px


.. list-table:: Figure 4. Comparison of peaks detected by different algorithms.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/igv1.png
   			:width: 400px


.. list-table:: Figure 5. Comparison of peaks detected by different algorithms.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/igv3.png
   			:width: 400px


.. list-table:: Figure 6. Comparison of peaks detected by different algorithms.
   :widths: 60
   :header-rows: 0

   * - .. image:: figures/igv4.png
   			:width: 400px



So which method to choose? You can test them all, we selected 4 (which produced some spurious peaks) and 5 (which we found had most support in read coverage in regions we inspected) for this tutorial.


.. code-block:: bash

	mkdir ../macs
	cd ../macs

	module load MACS/2.2.6

	macs2 callpeak -t ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam -n ENCFF045OAB.chr14.macs.broad_peaks --broad -f BAMPE -g 107043718 -q 0.1 --nomodel  --keep-dup all

	macs2 callpeak -t ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam -n ENCFF045OAB.chr14.macs.shift-75.extsize150_peaks -f BAM -g 107043718 -q 0.05 --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all



We chose genome size ``-g 107043718`` - because it is the length of chromosome 14, which is the only one included in the bam file.


.. Please note that we selected ``--extsize 100``  to match the behaviour of Genrich. Normally ``--extsize 200`` would be selected. ``--shift`` needs to be minus half of the size of ``--extsize`` to be centered on the 5’, so normally -100. ``--shift -100 --extsize 200`` will amplify the cutting sites' enrichment from ATAC-seq data. So in the end, the peak is where Tn5 transposase likes to attack.



.. Hint:: How to shift reads in BED files

	If you would like to test the effect of shifting reads, this how you do it on bed and bedpe files:

	.. code-block:: bash

		bedtools bamtobed -bedpe -i ENCFF045OAB.chr14.bam >ENCFF045OAB.chr14_pe.bed
		bedtools bamtobed -i ENCFF045OAB.chr14.bam >ENCFF045OAB.chr14.bed

		cat ENCFF045OAB.chr14.bed | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' >| ENCFF045OAB.chr14.proc.tn5.bed
	
		cat ENCFF045OAB.chr14_pe.bed | awk -F $'\t' 'BEGIN {OFS = FS} {$2 = $2 + 4; $6 = $6 - 5; print $0}' >| ENCFF045OAB.chr14.proc.tn5.pe.bed



How many peaks were detected?

.. code-block:: bash
	
	wc -l *Peak
	   2428 ENCFF045OAB.chr14.macs.shift-75.extsize150_peaks.narrowPeak
	   2011 ENCFF045OAB.chr14.macs.broad_peaks.broadPeak



Quite similar number of peaks, and double than what Genrich has detected.


.. admonition:: ENCFF045OAB.chr14.macs.broad_peaks.broadPeak
   :class: dropdown, warning

   .. code-block:: bash

		head ENCFF045OAB.chr14.macs.broad_peaks.broadPeak




Comparing results of MACS and Genrich
----------------------------------------

How many peaks actually overlap?

.. code-block:: bash
	
	cd ..

	bedtools intersect -a macs/ENCFF045OAB.chr14.macs.broad_peaks.broadPeak  -b genrich/ENCFF045OAB.chr14.genrich.narrowPeak  -f 0.50 -r >peaks_overlap.macs_broad.genrich.bed


	bedtools intersect -a macs/ENCFF045OAB.chr14.macs.shift-75.extsize150_peaks.narrowPeak -b genrich/ENCFF045OAB.chr14.genrich.narrowPeak  -f 0.50 -r >peaks_overlap.macs_narrow.genrich.bed

	wc -l peaks_common.bed 
	747 peaks_overlap.macs_broad.genrich.bed
	2327 peaks_overlap.macs_narrow.genrich.bed



:raw-html:`<br />`

:raw-html:`<br />`

:raw-html:`<br />`


