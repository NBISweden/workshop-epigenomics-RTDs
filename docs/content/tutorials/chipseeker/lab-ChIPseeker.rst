.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


===========================================
ChIP-seq downstream analysis: ChIPseeker
===========================================



**Learning outcomes**

Using ``ChIPseeker`` package

- to profile ChIP signal by genomics location and by ChIP binding to TSS regions

- to annotate peaks, visualise and compare annotations

- to run and compare functional enrichment


.. contents::
    :local:


Introduction
=============


In this tutorial we use another package, ``ChIPseeker``, to have a look at the ChIP profiles, annotate peaks and visualise annotations as well as to run functional enrichment. In a way, ``ChIPseeker`` can be seen as an alternative and newer workflow to ``ChIPpeakAnno`` (introduced in :doc:`differential binding <../diffBind/lab-diffBinding-remote>`). It also offers additional functionality, e.g. especially when it comes to visualising ChIP profiles and comparing functional annotations.


*It supports annotating ChIP peaks and provides functions to visualize ChIP peaks coverage over chromosomes and profiles of peaks binding to TSS regions. Comparison of ChIP peak profiles and annotation are also supported. Moreover, it supports evaluating significant overlap among ChIP-seq datasets. Currently, ChIPseeker contains 17,000 bed file information from GEO database. These datasets can be downloaded and compare with user’s own data to explore significant overlap datasets for inferring co-regulation or transcription factor complex for further investigation.*


Data & Methods
===============

We will build upon the main labs, using the same dataset and results from ``DiffBind`` analyses that we have saved under ``DiffBind.RData``. The tutorial is based on the `ChIPseeker package tutorial <https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html>`_ so feel free to have this open alongside to read and experiment more.

:raw-html:`<br />`

Setting-up
===========


You can continue working in the ``diffBind`` directory. We need access to ``diffBind.RData`` object, and some libraries, whcih are preinstalled. We access them via:


.. code-block:: bash

	   module load R_packages/4.0.4


In an ``R`` session:


.. code-block:: R

	# Load libraries (install if needed)
	library(DiffBind)
	library(ChIPseeker)
	library(ReactomePA)
	library(clusterProfiler)
	library(biomaRt)

	library(org.Hs.eg.db)
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



:raw-html:`<br />`

ChIP profile
==============


ChIP peaks coverage plot
------------------------

After peak calling one may want to visualise distribution of peaks locations over the whole genome. Function ``covplot`` calculates coverage of peaks regions over chromosomes.

Let's use data saved in ``DiffBind.RData`` objects. From this object we can easily extract peaks called for all our libraries as well as consensus peakset. In principle, we could also use ``ChIPseeker`` on raw ``.BED`` files.

.. code-block:: R

	# Let's start fresh removing all objects from R environment
	rm(list = ls())

	# loading diffBind.RData
	load("diffBind.RData")

	# Do you remember what objects we have saved in the diffBind.RData
	ls()

	# res.cnt3 object was the final one containing consensus peaks and differential binding results

	# viewing all samples
	dba.show(res.cnt3)

	# this should show you our 8 libraries
	> dba.show(res.cnt3)
	          ID Tissue Factor Replicate Caller Intervals   Reads FRiP
	1 REST_chip1   HeLa   REST         1 counts      5343 1637778 0.09
	2 REST_chip2   HeLa   REST         2 counts      5343 1991560 0.06
	3 REST_chip3 neural   REST         1 counts      5343 3197782 0.04
	4 REST_chip4 neural   REST         2 counts      5343 4924672 0.05
	5 REST_chip5  HepG2   REST         1 counts      5343 2988915 0.03
	6 REST_chip6  HepG2   REST         2 counts      5343 4812034 0.04
	7 REST_chip7  sknsh   REST         1 counts      5343 2714033 0.07
	8 REST_chip8  sknsh   REST         2 counts      5343 4180463 0.04


Please note the number of intervals (i.e. peaks) is **5343**. This is different from the original consensus peakset which had 6389 peaks. This original data is present in object ``cnt.res2``. This is because 1046 peaks fall in the internal blacklisted regions. At the same time, the object which holds the results ``res.cnt3`` contains information which peaks are detected in which sample, and this matrix is still in the original peakset format (i.e. has 6389 rows).

Information on the consensus peakset in ``res.cnt3``::

	> head(res.cnt3$called, n=3)
	     REST_chip1 REST_chip2 REST_chip3 REST_chip4 REST_chip5 REST_chip6
	[1,]          0          0          0          0          1          1
	[2,]          0          0          0          0          1          1
	[3,]          0          0          0          0          1          1
	     REST_chip7 REST_chip8
	[1,]          0          1
	[2,]          0          1
	[3,]          1          1

	> nrow(res.cnt3$called)
	[1] 6389


We have to do some data wrangling on this matrix to extract the rows of interest to us, i.e. rows corresponding to peaks which were NOT blacklisted. The code to do this is presented below. You may copy - paste it and you'll arrive at the correct object to continue working. If you would like to understand what's happening, you can inspect the objects created in each step using commands ``head``, ``nrow`` etc.

.. code-block:: R

	#all peaks including blacklisted, this corresponds to our object of interest res.cnt3$called
	peaks.all=res.cnt2$peaks[[1]]

	#peaks in blacklists
	peaks.blck=as.data.frame(res.cnt3$peaks.blacklisted[1])

	#clean up colnames
	library(janitor)
	peaks.all=peaks.all %>% clean_names()

	#merge between the two objects
	library(dplyr)
	peaks.all.blck=left_join(peaks.all,peaks.blck, by=c("start", "chr"="seqnames" ))

	#indices of peaks NOT in blacklists > to keep in the data
	peaks.all.no_blck.ind=is.na(peaks.all.blck$group)

	#examine if the numbers add up
	table(peaks.all.no_blck.ind)
	##	peaks.all.no_blck.ind
	##	FALSE  TRUE
	##	 1046  5343

	#subset res.cnt3$called
	called.peaks=res.cnt3$called[peaks.all.no_blck.ind,]

	nrow(called.peaks)
	## [1] 5343


To plot peaks over genomic locations we need to extract from ``res.cnt3`` peaks of interest, e.g. consensus peaks or present in a single replicate etc. Here, we will focus on peaks present in HeLa replicates.

.. code-block:: R

	# extracting consensus peak set with 5343 peaks
	peaks.consensus <- dba.peakset(res.cnt3, bRetrieve = T)

``peaks.consensus`` is a ``GRangers`` object::

	> peaks.consensus
	GRanges object with 5343 ranges and 8 metadata columns:
	       seqnames              ranges strand | REST_chip1 REST_chip2 REST_chip3
	          <Rle>           <IRanges>  <Rle> |  <numeric>  <numeric>  <numeric>
	     1     chr1         29190-29590      * |          0          0          0
	     2     chr1       100300-100700      * |          0          0          0
	     3     chr1       151013-151413      * |          0          0          0
	     4     chr1       246634-247034      * |          0          0          0
	     5     chr1       408268-408668      * |          0          0          0
	   ...      ...                 ...    ... .        ...        ...        ...
	  5339     chr2 242910501-242910901      * |    0.00000    0.00000   48.35128
	  5340     chr2 243011991-243012391      * |    0.00000    0.00000    0.00000
	  5341     chr2 243030594-243030994      * |    0.00000    9.46756    4.83513
	  5342     chr2 243093019-243093419      * |    2.56875    0.00000    0.00000
	  5343     chr2 243184803-243185203      * |   77.06258  156.21468    0.00000
	       REST_chip4 REST_chip5 REST_chip6 REST_chip7 REST_chip8
	        <numeric>  <numeric>  <numeric>  <numeric>  <numeric>
	     1          0   58.99879    59.8714    36.0184    45.3907
	     2          0    1.07271    71.7064    57.9120    14.8551
	     3          0    2.14541   110.6924   120.0614    44.5654
	     4          0    1.07271    75.1873    58.6182    14.0298
	     5          0    4.29082   130.1854   101.6991    48.6918
	   ...        ...        ...        ...        ...        ...
	  5339   28.73111    1.07271     0.0000    0.00000     0.0000
	  5340    0.00000    0.00000    12.5312    7.76868     0.0000
	  5341    6.91675    2.14541    18.7968    9.88741    12.3793
	  5342    0.00000    6.43623    67.5293   48.73080    23.9333
	  5343    0.00000    2.14541   188.6644  160.31726   120.4916
	  -------
	  seqinfo: 2 sequences from an unspecified genome; no seqlengths


We select interesting peaks and work on them. First let's check the peak locations and scores along the chromosomes.


.. code-block:: R


	# extracting HeLA peaks
	peaks.HeLa_rep1 <- peaks.consensus[called.peaks[,1]==1] # peaks called in rep 1
	peaks.HeLa_rep2 <- peaks.consensus[called.peaks[,2]==1] # peaks called in rep 2

	# adding an unified affinity scores column (re-formatting data)
	peaks.HeLa_rep1$Score <- peaks.HeLa_rep1$REST_chip1
	peaks.HeLa_rep2$Score <- peaks.HeLa_rep2$REST_chip2

	# plotting coverage for replicate 1, using affinity scores as a weight for peaks height
	covplot(peaks.HeLa_rep1, weightCol = "Score")

	# zooming in to a selected region is also possible
	covplot(peaks.HeLa_rep1, weightCol = "Score", xlim=c(0, 1e07))

	#save the plots
	pdf("chipseeker-coverage-plots-HeLa-r1.pdf")
	covplot(peaks.HeLa_rep1, weightCol = "Score")
	covplot(peaks.HeLa_rep1, weightCol = "Score", xlim=c(0, 1e07))
	dev.off()

.. admonition:: chipseeker-coverage-plots-HeLa-r1.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-coverage-plots-HeLa-r1-0.png
            :width: 600px


We can also compare peaks across replicates. This should give us visual assessment of variability between replicates: peaks locations and strength should match in an ideal scenario.

.. code-block:: R

	# creating genomicRangesList object holding replicates 1 and 2
	grL.HeLa = GRangesList(HeLa_rep1=peaks.HeLa_rep1, HeLa_rep2=peaks.HeLa_rep2, compress=FALSE)


	# plotting using affinity scores as a weight for peaks height
	covplot(grL.HeLa, weightCol = "Score")

	# zooming in
	covplot(grL.HeLa, weightCol = "Score", xlim=c(0, 1e07))

	#save the plots
	pdf("chipseeker-coverage-plots-HeLa-r1-r2.pdf")
	covplot(grL.HeLa, weightCol = "Score")
	covplot(grL.HeLa, weightCol = "Score", xlim=c(0, 1e07))
	dev.off()


.. admonition:: chipseeker-coverage-plots-HeLa-r1-r2.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-coverage-plots-HeLa-r1-r2-0.png
            :width: 600px


What do you think?

- are these peaks reproducible?

- which pair of replicates is most consistent, HeLa, neural, HepG2 or sknsh? (hint: you may need to generate more plots to answer this)

- why is it good to always look at the data instead of simply trusting the output of the summary statistics, after all, we do rely on ``diffBind`` to call peaks being consistent?


Profile of ChIP peaks binding to TSS regions
--------------------------------------------


For calculating the profile of ChIP peaks binding to TSS regions, we need to prepare the TSS regions, which are defined as the flanking sequence of the TSS sites. Then we can align the peaks that are mapping to these regions, and generate the tagMatrix used for plotting.

Here, we will select peaks present per cell type, i.e. found in two replicates. We will also create tagMatrix list to enable group comparisons across cell lines.

.. code-block:: R

	# extracting peaks for each cell line present across replicates
	peaks.HeLa <- peaks.consensus[called.peaks[,1]==1 & called.peaks[,2]==1]
	peaks.neural <- peaks.consensus[called.peaks[,3]==1 & called.peaks[,4]==1]
	peaks.HepG2 <- peaks.consensus[called.peaks[,5]==1 & called.peaks[,6]==1]
	peaks.sknsh <- peaks.consensus[called.peaks[,7]==1 & called.peaks[,8]==1]

	# getting TSS regions
	promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

	# calculating tagMatrix
	tagMatrix.1 <- getTagMatrix(peaks.HeLa, windows=promoter)
	tagMatrix.2 <- getTagMatrix(peaks.neural, windows=promoter)
	tagMatrix.3 <- getTagMatrix(peaks.HepG2, windows=promoter)
	tagMatrix.4 <- getTagMatrix(peaks.sknsh, windows=promoter)

	# preparing tagMatrix list to enable cell lines comparisions
	tagMatrixList <- list(HeLa=tagMatrix.1, neural=tagMatrix.2, HepG2=tagMatrix.3, sknsh=tagMatrix.4)

	# plotting tagMatrix heatmaps for each cell line
	tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

	# plotting average profile of ChIP peaks among different cell lines
	plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

	#save the plot
	pdf("chipseeker-average-peak-profile.pdf")
	plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
	dev.off()


.. admonition:: chipseeker-average-peak-profile.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-average-peak-profile.png
            :width: 600px



:raw-html:`<br />`

Peaks Annotation
===================

**Peak annotations is performed by annotatePeak() function**

Here, we can define TSS region, by default set to -3kb to 3kb. The output of ``annotatePeak`` is ``csAnno`` object than we can convert to ``GRanges`` with ``as.GRanges()`` function or to data frame with ``as.data.frame()`` function.

Similar to annotations with ``ChIPpeakAnno`` we will need ``TxDB`` object containing annotations, transcript-related features of a particular genome. We can use Bioconductor packages providing annotations for various model organisms. It may be however **good to know that one can also prepare their own TxDb object** by retrieving information from UCSC or BioMart using ``GenomicFeature`` package. Here, we will use ``TxDb.Hsapiens.UCSC.hg19.knownGene`` annotations provided by Bioconductor.

Some **annotations may overlap** and by default ChIPseeker annotates peaks with the priority: promoter, 5' UTR, 3' UTR, exon, intron, downstreamn, intergenic, where downstream is defined as the downstream of gene end. This priority can be changed with ``genomicAnnotationPriority`` parameter.

While annotating peaks we can include optional parameter ``annoDb`` containig further genome wide annotation data. If added, this will add SYMBOL, GENENAME, ENSEMBL/ENTREZID to the peaks annotations. Again, we will use Bioconductor ``org.Hs.eg.db`` for human genome wide annotation data.


.. code-block:: R

	# extracting all consensus peaks (repeating commands for clarity)
	peaks.consensus <- dba.peakset(res.cnt3, bRetrieve = T)

	# extracting peaks for each cell line present across replicates (repeating commands for clarity)
	peaks.HeLa <- peaks.consensus[res.cnt3$called[,1]==1 & res.cnt3$called[,2]==1]
	peaks.neural <- peaks.consensus[res.cnt3$called[,3]==1 & res.cnt3$called[,4]==1]
	peaks.HepG2 <- peaks.consensus[res.cnt3$called[,5]==1 & res.cnt3$called[,6]==1]
	peaks.sknsh <- peaks.consensus[res.cnt3$called[,7]==1 & res.cnt3$called[,8]==1]

	# annotating peaks
	peaks.HeLa_ann <- annotatePeak(peaks.HeLa, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
	peaks.neural_ann <- annotatePeak(peaks.neural, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
	peaks.HepG2_ann <- annotatePeak(peaks.HepG2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
	peaks.sknsh_ann <- annotatePeak(peaks.sknsh, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

	# previewing annotations summary for HeLa peaks
	peaks.HeLa_ann

	> peaks.HeLa_ann
	Annotated peaks generated by ChIPseeker
	996/996  peaks were annotated
	Genomic Annotation Summary:
	              Feature  Frequency
	9    Promoter (<=1kb) 12.6506024
	10   Promoter (1-2kb)  5.2208835
	11   Promoter (2-3kb)  3.9156627
	4              5' UTR  0.2008032
	3              3' UTR  1.5060241
	1            1st Exon  0.1004016
	7          Other Exon  3.3132530
	2          1st Intron 11.4457831
	8        Other Intron 21.2851406
	6  Downstream (<=300)  1.2048193
	5   Distal Intergenic 39.1566265


	# previewing peaks annotations for HeLa peaks
	head(as.data.frame(peaks.HeLa_ann))


.. admonition:: peaks.HeLa_ann
   :class: dropdown, warning


   .. code-block:: R

			  seqnames   start     end width strand REST_chip1 REST_chip2 REST_chip3
		1     chr1 1234433 1234833   401      *  652.46315  347.93269   6.043910
		2     chr1 1270240 1270640   401      *  783.46953  639.06004   8.461474
		3     chr1 1408222 1408622   401      *   51.37505   68.63978   0.000000
		4     chr1 1563883 1564283   401      *   97.61260   23.66889   0.000000
		5     chr1 1714063 1714463   401      *   43.66879   30.76956   4.835128
		6     chr1 1875369 1875769   401      *  431.55043  142.01334   6.043910
		  REST_chip4 REST_chip5 REST_chip6 REST_chip7 REST_chip8        annotation
		1   7.980863   0.000000  66.833141  125.00509   88.30549  Promoter (<=1kb)
		2   9.044979   4.290821 228.346566  194.92319  127.91917 Downstream (<1kb)
		3   0.000000   5.363527  18.100642   10.59365   26.40912  Promoter (1-2kb)
		4   0.000000   0.000000   2.088536    0.00000    4.95171  Promoter (<=1kb)
		5   2.128230   0.000000   3.480893    0.00000    3.30114  Promoter (2-3kb)
		6   2.660288   7.508937  58.478998  120.76763  100.68476 Distal Intergenic
		  geneChr geneStart geneEnd geneLength geneStrand geneId transcriptId
		1       1   1227764 1234335       6572          2 116983   uc001ady.2
		2       1   1266726 1269844       3119          1  83756   uc010nyk.2
		3       1   1407164 1431582      24419          1  83858   uc001afv.3
		4       1   1564486 1565990       1505          1 142678   uc001ago.3
		5       1   1682671 1711508      28838          2  65220   uc001aie.3
		6       1   1849029 1850740       1712          2 339456   uc001aij.2
		  distanceToTSS         ENSEMBL SYMBOL
		1           -98 ENSG00000131584  ACAP3
		2          3514 ENSG00000169962 TAS1R3
		3          1058 ENSG00000160072 ATAD3B
		4          -203 ENSG00000197530   MIB2
		5         -2555 ENSG00000008130   NADK
		6        -24629 ENSG00000178821 TMEM52
		                                                  GENENAME
		1 ArfGAP with coiled-coil, ankyrin repeat and PH domains 3
		2                                taste 1 receptor member 3
		3                   ATPase family AAA domain containing 3B
		4                   mindbomb E3 ubiquitin protein ligase 2
		5                                               NAD kinase
		6                                 transmembrane protein 52




We find our genomic annotations in _annotation_ column. Plots, pie and barplot, are supported to visualise these annotations.

.. code-block:: R

	# creating barplot for HeLa peaks genomics annotations
	plotAnnoBar(peaks.HeLa_ann)

	# creating vennpie plot
	vennpie(peaks.HeLa_ann)

	# creating upsetplot showing overlapping annotations
	upsetplot(peaks.HeLa_ann)




Let's save these plots:

.. code-block:: R

	pdf("chipseeker-annoplots-hela.pdf")
	plotAnnoBar(peaks.HeLa_ann)
	vennpie(peaks.HeLa_ann)
	upsetplot(peaks.HeLa_ann)
	dev.off()

.. admonition:: chipseeker-annoplots-hela.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-annoplots-hela-0.png
            :width: 600px

   .. image:: figures/chipseeker-annoplots-hela-1.png
            :width: 600px

   .. image:: figures/chipseeker-annoplots-hela-2.png
            :width: 600px



We can also use ``plotAnnoBar`` to compare annotations between different datasets, here cell lines. For that, we just need to create a list containing peaks annotations of datasets to compare.

.. code-block:: R

	# creating list holding annotations for different cell lines
	list.annotations <- list(HeLa=peaks.HeLa_ann, neural=peaks.neural_ann, HepG2=peaks.HepG2_ann, sknskh=peaks.sknsh_ann)

	# creating barplot for HeLa, neural, HepG2 and sknsh peaks genomic annotations
	plotAnnoBar(list.annotations)



Finally, we can also visualise distribution of TF-binding loci relative to TSS, for single annotation set or using annotations list for comparisons.

.. code-block:: R

	# plotting distance to TSS for HeLa peaks
	plotDistToTSS(peaks.HeLa_ann)

	# plotting distance to TSS for all cell lines in our annotation list
	plotDistToTSS(list.annotations)


	pdf("chipseeker-DistToTSS-hela.pdf")
	plotDistToTSS(peaks.HeLa_ann)
	plotDistToTSS(list.annotations)
	dev.off()


.. admonition:: chipseeker-DistToTSS-hela.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-DistToTSS-hela-1.png
            :width: 600px


What do you think?

- would you expect such distribution of features?

- do these distributions differ between cell-lines?

:raw-html:`<br />`


Functional analysis
===================

Having obtained annotations to nearest genes, we can perform **functional enrichment analysis to identify predominant biological themes** among these genes by incorporating biological knowledge provided by biological ontologies, incl. GO (Gene Ontology, Ashburner et al. 2000), KEGG (Kyoto Encyclopedia of Genes and Genomes, Kanehisa et al. 2004), DO (Disease Ontology, Schriml et al. 2011) or Reactome (Croft et al. 2013).

Here, we can also use ``seq2gene`` function for linking genomic regions to genes in a **many-to-many mapping**. This function consider host gene (exon/intron), promoter region and flanking gene from intergenic region that may undergo control via cis-regulation.

One can **build on** using ChIPseeker for functional enrichment and annotation as there are several packages by the same author to identify biological themes, i.e. ``ReactomePA`` for reactome pathways enrichment, ``DOSE`` for Disease Ontology, ``clusterProfiler`` for Gene Ontology and KEGG enrichment analysis. Especially `clustserProfiler <http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html>`_
comes handy when **visualising and comparing** biological themes, also when comparing functions derived from other omics technologies for integrative analyses.

Here, we will experiment with few functions only. We will search for enriched reactome pathways using genes annotated to peaks by nearest location and allowing for many-to-many mapping. We will also learn how to compare functional annotations between peak sets using GO terms as an example.

We will start by defying our genes background, i.e. genes on chromosome 1 and 2. For this we can use functions from ``biomaRt``


.. code-block:: R

	# defining chromosomes
	chrom=c(1,2)

	# defining source
	ensembl=useMart("ensembl")
	ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

	# running query: extracting ENTREZID for genes on chromosome 1 and 2
	genes.chr1chr2 <- getBM(attributes= "entrezgene_id",
	        filters=c("chromosome_name"),
	        values=list(chrom), mart=ensembl)

	# reformatting output to character string (as required later on by clusterProfiler functions)
	genes.universe <- as.character(as.numeric(as.matrix(genes.chr1chr2)))



Reactome pathway enrichment of genes defined as a) nearest feature to the peaks and b) allowing for many-to-many mapping

.. code-block:: R

	# a: selecting annotated peaks for functional enrichment in object
	data.peaks_ann <- peaks.neural_ann

	# a: finding enriched Reactome pathways using chromosome 1 and 2 genes as a background
	pathway.reac1 <- enrichPathway(as.data.frame(data.peaks_ann)$geneId, universe = genes.universe)

	# a: previewing enriched Reactome pathways
	head(pathway.reac1)


This is the overrepresented pathway::

	> head(pathway.reac1)
	                       ID     Description GeneRatio BgRatio      pvalue
	R-HSA-112316 R-HSA-112316 Neuronal System    29/435 61/1797 4.56747e-05
	               p.adjust     qvalue
	R-HSA-112316 0.02101036 0.02101036
	                                                                                                                                                                 geneID
	R-HSA-112316 2782/8514/57576/58512/2899/55970/5567/3737/3752/3782/777/2752/127833/3756/3776/3775/3754/3790/170850/9378/80059/347730/60482/785/3760/90134/2571/2744/1385
	             Count
	R-HSA-112316    29


.. code-block:: R


	# b: selecting peaks
	data.peaks <- peaks.HeLa

	# b: running seq2gene function for many-to-many mapping based on sequence regions (note: no prior peaks annotations here, many-to-many mapping is done from the sequence)
	genes.m2m <- seq2gene(data.peaks, tssRegion = c(-3000, 3000), flankDistance = 3000, TxDb=txdb)

	# b: finding enriched Reactome pathways given many to many mapping and chromosome 1 and 2 genes as a background
	pathway.reac2 <- enrichPathway(genes.m2m, universe = genes.universe)

	# b: creating dotplot to visualise enrichment results
	dotplot(pathway.reac2)

	#save the plot
	pdf("chipseeker-dotplot-reactome-HeLa.pdf")
	dotplot(pathway.reac2)
	dev.off()

.. admonition:: chipseeker-dotplot-reactome-HeLa.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-dotplot-reactome-HeLa.png
            :width: 600px




Let's search for enriched GO terms, and let's see how we can do it for all the peak sets together so we can easily compare the results on a ``dotplot``. Also, let's learn how to simplify the output of GO terms using ``simplify`` function, useful in cases where lots of GO terms turn-up to be significant and it becomes difficult to interpret results. ``simply`` function removes redundant GO terms obtained from ``encrichGO`` calling internally ``GoSemSim`` function to calculate similarities among GO terms and removes those highly similar terms by keeping one representative term.

.. code-block:: R

	# creating a gene list with ENTREZID ideas extracted from our annotation list, containing annotated peaks for all four cell lines
	list.genes = lapply(list.annotations, function(i) as.data.frame(i)$geneId)
	names(list.genes) <- sub("_", "\n", names(list.genes))

	# running enrichedGO function to find enriched MF correlation_libraries_normalised on the gene list

	compMF <- compareCluster(geneCluster = list.genes,
	                       fun           = "enrichGO",
	                       pvalueCutoff  = 0.05,
	                       pAdjustMethod = "BH",
	                       OrgDb='org.Hs.eg.db',
	                       ont="MF")

	# comparing results on a dotplot
	dotplot(compMF)

	# simplifying results although here we do not have problems with too many GO terms
	compMF.flr <- simplify(compMF, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)

	# creating a dotplot on reduced GO terms
	dotplot(compMF.flr)


And let's save the plots::

	pdf("chipseeker-GO-MF.pdf")
	 dotplot(compMF)
	 dotplot(compMF.flr)
	dev.off()


.. admonition:: chipseeker-GO-MF.pdf
   :class: dropdown, warning

   .. image:: figures/chipseeker-GO-MF-0.png
            :width: 600px

   .. image:: figures/chipseeker-GO-MF-1.png
            :width: 600px




Concluding remarks and next steps
======================================

There are different flavours to functional annotations, and what and how functional annotations should be done is context dependent, i.e. they should be adjusted given available data and biological question being asked. There are many methods out there, all relying on the available annotations and databases, being constantly improved and developed. As a rule of thumb to understand the results and be able to draw biological conclusions, it may be good to think about i) the statistical test behind the method, ii) what is compared against what (i.e. genes vs. background) and which databases are being used (i.e. Reactome, GO, DO, KEGG).

For more examples on what can be done in terms on functional annotations, we recommend reading tutorials on `clusterProfiler <http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#reduce-redundancy-of-enriched-go-terms>`_ and `DOSE <https://bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/DOSE.html>`_, where you can further learn about semantic similarity analysis, disease enrichment analysis, GSEA analysis and much more.


:raw-html:`<br />`



.. ----

.. Written by: Olga Dethelefsen
.. Modified by: Agata Smialowska
