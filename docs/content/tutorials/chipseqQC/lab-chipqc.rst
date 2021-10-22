.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

============================
Alternative QC using ChIPQC
============================

**Learning outcomes**

Using ``ChIPQC`` package

- to generate a summary QC report for experimental sample groups

- to be able to understand and assess QC metrics and plots

.. contents:: 
    :local:


Introduction
==============

Here, we will explore the alternative quality control workflow, using Bioconductor ``ChIPQC`` package. ``ChIPQC`` computes quality metrics for aligned data from ChIP-seq experiments. It also provides simple ways to generate a ChIP-seq experiment quality report which can be examined to asses the absolute and relative quality of individual ChIP-seq samples (and their associated controls, as well as overall quality of the experimental data.)



Setting-up
==============

In principle one can run ``ChIPQC`` both on Uppmax or locally. 
The instructions below are 
to use the package on Uppmax. For local usage, please make sure the paths to files are correct, and that all dependecies are properly installed for the version of R you are running.

.. We provide a conda environment to run the version of ``R`` for which the tutorial was set up and tested. To find how this environment was constructed, please visit :doc:`Dependencies <../../dependencies>`.


.. Follow set-up instructions from :doc:`Downstream analysis tutorial <../diffBind/lab-diffBinding-remote>`, differential binding part. We will need the same files and we can work in the same directory.



.. HINT::
	
	If you start from a different location, you should ``cd chipseq/analysis/R``


You can now load the version of R for which we tested this class along with other dependencies:


.. .. code-block:: bash
	
.. 	cd analysis/R

.. 	module load conda/latest
	
.. 	conda activate /sw/courses/epigenomics/software/conda/v8

.. 	R


.. code-block:: bash

   module load R_packages/4.0.4

The remaining part of the exercise is performed in ``R``.



.. HINT::

	We are running 
	``R version 4.0.4 (2021-02-15) -- "Lost Library Book"``



Running ChIPQC
================

While running commands (in ``R``) feel free to have a look at `ChIPQC package documentation <http://bioconductor.org/packages/devel/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf>`_ to learn more about different steps and/or build upon them. Here we will just show you the very basics.


.. code-block:: R

	library(DiffBind)
	library(ChIPQC)
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)


	#	reading in the sample information (metadata)
	samples = read.csv("samples_REST.txt", sep="\t")

	#	inspecting the metadata
	samples

	#	creating an object containing data
	res=dba(sampleSheet=samples, config=data.frame(RunParallel=FALSE))

	# inspecting the object
	res


.. admonition:: res
   :class: dropdown, warning

   .. code-block:: R

	   > res
		8 Samples, 6518 sites in matrix (17056 total):
		          ID Tissue Factor Replicate Intervals
		1 REST_chip1   HeLa   REST         1      2252
		2 REST_chip2   HeLa   REST         2      2344
		3 REST_chip3 neural   REST         1      5948
		4 REST_chip4 neural   REST         2      3003
		5 REST_chip5  HepG2   REST         1      2663
		6 REST_chip6  HepG2   REST         2      4326
		7 REST_chip7  sknsh   REST         1      8700
		8 REST_chip8  sknsh   REST         2      3524


We will now perform QC using a wrapper function which does all the heavy lifting (well, typing...) for us. This function only uses the data mapped the first chromosome in the reference. It is enough to get representative results.


.. code-block:: R

	#	performing quality control
	resqc = ChIPQC(res,annotation="hg19", config=data.frame(RunParallel=TRUE))


Finally, we save the results for later viewing:

.. code-block:: R

	#	creating the quality control report in html format
	ChIPQCreport(resqc)


.. WARNING::
	
	If you run this tutorial on Rackham, you may see an error

		``ChIPQCreport(resqc)``

		``Error in browseURL...`` : ``'browser' must be a non-empty character string``

  	This is because the html report cannot be open in a browser directly from Rackham. You can download it to your computer and view it locally.


You need to copy the report to your local computer (copy the entire ``ChIPQCreport`` folder):

.. code-block:: bash
	
	scp -r <USER>@rackham.uppmax.uu.se:/path/to/ChIPQCreport .

	#if you follow the paths used in this tutorial
	scp -r <USER>@rackham.uppmax.uu.se:~/chipseq/analysis/R/ChIPQCreport .


Examine the html report.

What do you think?

Are these results in line with the previous quality control workflow?




.. admonition:: relevant information from sessionInfo()
   :class: dropdown, warning


   .. code-block:: R

	   other attached packages:
	 [1] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
	 [2] GenomicFeatures_1.42.3                 
	 [3] AnnotationDbi_1.52.0                   
	 [4] ChIPQC_1.26.0                          
	 [5] ggplot2_3.3.3                          
	 [6] DiffBind_3.0.15                        
	 [7] SummarizedExperiment_1.20.0            
	 [8] Biobase_2.50.0                         
	 [9] MatrixGenerics_1.2.1                   
	[10] matrixStats_0.58.0                     
	[11] GenomicRanges_1.42.0                   
	[12] GenomeInfoDb_1.26.7                    
	[13] IRanges_2.24.1                         
	[14] S4Vectors_0.28.1                       
	[15] BiocGenerics_0.36.0 



----------

.. The report can be also downloaded from Box [here](https://stockholmuniversity.box.com/s/c1lbrr1s1khw4ctiqfq0f9j2m1b6vp90)


.. ----

.. Written by: Agata Smialowska
