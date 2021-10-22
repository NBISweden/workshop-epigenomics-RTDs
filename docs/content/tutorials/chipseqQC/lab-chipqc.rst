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

We provide a conda environment to run the version of ``R`` for which the tutorial was set up and tested. To find how this environment was constructed, please visit :doc:`Dependencies <../../dependencies>`.


Follow set-up instructions from :doc:`Downstream analysis tutorial <../diffBind/lab-diffBinding-remote>`, differential binding part. We will need the same files and we can work in the same directory.


.. Install ``ChIPQC`` library and any required dependencies

.. ```bash

.. if (!requireNamespace("BiocManager", quietly = TRUE))
..     install.packages("BiocManager")
.. BiocManager::install("ChIPQC")

.. ```

.. HINT::
	
	If you start from a different location, you should ``cd chipseq/analysis/R``


You can now load the version of R for which we tested this class along with other dependencies:


.. code-block:: bash
	
	cd analysis/R

	module load conda/latest
	
	conda activate /sw/courses/epigenomics/software/conda/v8

	R


.. HINT::

	We are running 
	``R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"``



Running ChIPQC
================

While running commands (in ``R``) feel free to have a look at `ChIPQC package documentation <http://bioconductor.org/packages/devel/bioc/vignettes/ChIPQC/inst/doc/ChIPQC.pdf>`_ to learn more about different steps and/or build upon them. Here we will just show you the very basics.


.. code-block:: R

	#	add the path to the tutorial-specific R libraries
	assign(".lib.loc", "/sw/courses/epigenomics/software/R", envir = environment(.libPaths))

	#	you can see that the tutorial-specific R library path is added
	.libPaths()
	[1] "/sw/courses/epigenomics/software/R"

	#
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

	#	performing quality control
	resqc = ChIPQC(res,annotation="hg19", config=data.frame(RunParallel=TRUE))

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


----------

.. The report can be also downloaded from Box [here](https://stockholmuniversity.box.com/s/c1lbrr1s1khw4ctiqfq0f9j2m1b6vp90)


.. ----

.. Written by: Agata Smialowska
