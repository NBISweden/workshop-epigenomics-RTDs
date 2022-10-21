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

.. Important::

	We assume that the environment and directory structure has been already set in :doc:`Data preprocessing <../data-preproc/data-preproc>`.



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
	ln -s ../../processedData/ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam ENCFF045OAB.chr14.proc.bam
	ln -s ../../processedData/ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam.bai ENCFF045OAB.chr14.proc.bam.bai


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


	# sort the bam file by read name (required by Genrich)
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

	You can inspect file contents:

   .. code-block:: bash

		head ENCFF045OAB.chr14.genrich.narrowPeak
		chr14	18332340	18333050	peak_0	490	.	347.820770	3.316712	-1	372
		chr14	18583390	18584153	peak_1	1000	.	1267.254150	6.908389	-1	474
		chr14	19061839	19062647	peak_2	732	.	591.112671	4.484559	-1	472
		chr14	19337360	19337831	peak_3	1000	.	517.304626	4.484559	-1	373
		chr14	19488517	19489231	peak_4	393	.	280.354828	2.916375	-1	210
		chr14	20216750	20217291	peak_5	1000	.	625.121826	4.537151	-1	441




MACS2
-----

MACS2 is widely used for peak calling in ATAC-seq, as evidenced by literature and many data processing pipelines. Several different peak calling protocols / commands have been encountered in various sources (and more combinations of parameters exist):

1. macs narrow peak (default for ``callpeak``), ``--nomodel``, shifted reads using PE reads as SE (BED file);

2. macs narrow peak, ``--nomodel``, shifted reads, using PE reads (BEDPE file);

3. macs narrow peak, ``--nomodel``, using PE reads (BEDPE file); ``shift`` and ``extsize`` similar to Genrich;

4. macs broad peak, using PE reads (BAMPE file) as in nf-core;

5. macs narrow peak, ``--nomodel``,  using PE reads as SE (BAM file) as in Encode, ``shift`` and ``extsize`` similar to Genrich;

6. macs narrow peak, unshifted reads in BAMPE file

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



So which method to choose? You can test them all, for this tutorial we selected 4 and 5, which we found had most support in read coverage in regions we inspected yet still produced some spurious peaks. We'll call the results ``broad`` and ``narrow``, respectively. We use ``-g 107043718`` for peak calling because this is the length of chr14, which is the only one included in the bam file.

.. code-block:: bash

	mkdir ../macs
	cd ../macs

	module load MACS/2.2.6

	macs2 callpeak -t ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam \
	-n ENCFF045OAB.chr14.macs.broad --broad -f BAMPE \
	-g 107043718 -q 0.1 --nomodel  --keep-dup all

	macs2 callpeak -t ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam \
	-n ENCFF045OAB.chr14.macs.narrow -f BAM \
	-g 107043718 -q 0.05 --nomodel --shift -75 --extsize 150 \
	--call-summits --keep-dup all



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
	 2011 ENCFF045OAB.chr14.macs.broad_peaks.broadPeak
  	 2428 ENCFF045OAB.chr14.macs.narrow_peaks.narrowPeak
 


Quite similar number of peaks for both methods, and double than what Genrich has detected.



.. admonition:: ENCFF045OAB.chr14.macs.broad_peaks.broadPeak
   :class: dropdown, warning

   .. code-block:: bash

		head ENCFF045OAB.chr14.macs.broad_peaks.broadPeak
		chr14	18674026	18674550	ENCFF045OAB.chr14.macs.broad_peak_1	21	.	3.15137	4.25082	2.12270
		chr14	19096643	19097148	ENCFF045OAB.chr14.macs.broad_peak_2	31	.	3.61047	5.39085	3.19165
		chr14	19098499	19098851	ENCFF045OAB.chr14.macs.broad_peak_3	16	.	2.95556	3.77788	1.68306
		chr14	19105556	19105809	ENCFF045OAB.chr14.macs.broad_peak_4	27	.	3.42266	4.93515	2.76765
		chr14	19161075	19162012	ENCFF045OAB.chr14.macs.broad_peak_5	24	.	3.29408	4.58492	2.42759
		chr14	19172872	19173211	ENCFF045OAB.chr14.macs.broad_peak_6	22	.	3.19527	4.34825	2.20892


.. admonition:: ENCFF045OAB.chr14.macs.narrow_peaks.narrowPeak
   :class: dropdown, warning

   .. code-block:: bash

		head ENCFF045OAB.chr14.macs.narrow_peaks.narrowPeak
		chr14	19372856	19373058	ENCFF045OAB.chr14.macs.narrow_peak_1	17	.	3.18021	3.90652	1.78714	160
		chr14	19374426	19374806	ENCFF045OAB.chr14.macs.narrow_peak_2	65	.	5.27241	8.87504	6.51449	120
		chr14	19388860	19389063	ENCFF045OAB.chr14.macs.narrow_peak_3	62	.	4.76858	8.64439	6.29014	89
		chr14	19889924	19890074	ENCFF045OAB.chr14.macs.narrow_peak_4	49	.	4.67007	7.24319	4.93076	131
		chr14	20093651	20093822	ENCFF045OAB.chr14.macs.narrow_peak_5	25	.	3.59236	4.74566	2.54154	91




Comparing results of MACS and Genrich
----------------------------------------

How many peaks actually overlap?

.. code-block:: bash
	
	cd ..

	module load BEDTools/2.25.0

	bedtools intersect -a macs/ENCFF045OAB.chr14.macs.broad_peaks.broadPeak  -b genrich/ENCFF045OAB.chr14.genrich.narrowPeak  -f 0.50 -r >peaks_overlap.macs_broad.genrich.bed


	bedtools intersect -a macs/ENCFF045OAB.chr14.macs.narrow_peaks.narrowPeak -b genrich/ENCFF045OAB.chr14.genrich.narrowPeak  -f 0.50 -r >peaks_overlap.macs_narrow.genrich.bed

	wc -l *bed
	   747 peaks_overlap.macs_broad.genrich.bed
 	 1613 peaks_overlap.macs_narrow.genrich.bed


Fraction of Reads in Peaks
-----------------------------




:raw-html:`<br />`

Consensus Peaks
===================


Merged Peaks
===================


:raw-html:`<br />`

:raw-html:`<br />`


