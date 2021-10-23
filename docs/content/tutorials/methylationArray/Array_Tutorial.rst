DNA Methylation: Array Workflow
===============================

**Learning Outcomes**

In this tutorial, we will provide examples of the steps involved in analyzing 450K methylation array data using R and Bioconductor. The different steps include: importing the raw data, quality control checks, data filtering, different normalization methods and probe-wise differential methylation analysis. Additional approaches such as differential methylation analysis of regions, gene ontology analysis and estimating cell type composition will also be presented. 

.. Contents
.. ========

.. contents:: 
    :local:

Introduction
------------

Despite the increasing popularity of sequencing based methods, methylation arrays remain the platform of choice for many epigenome-wide association studies. Their user-friendly and more streamlined data analysis workflow in combination with a lower price per sample make them the preferred tool for - especially larger scale - studies. In this tutorial, an overview of a typical analysis of a Illumina HumanMethylation450 array will be presented. 

But first; a bit of history. Measurement of DNA methylation by Infinium technology (Infinium I) was first employed by Illumina on the HumanMethylation27 (27k) array, which measured methylation at approximately 27,000 CpGs, primarily in gene promoters. Like bisulfite sequencing, the Infinium assay detected methylation status at single base resolution. However, due to its relatively limited coverage the array platform was not truly considered “genome-wide” until the arrival of the 450k array. Introduced in 2011, the 450k array increased the genomic coverage of the platform to over 450,000 gene-centric sites by combining the original Infinium I probes with the novel Infinium II probes. Both probe types employ 50bp probes that query a [C/T] polymorphism created by bisulfite conversion of unmethylated cytosines in the genome. However, the Infinium I and II assays differ in the number of beads required to detect methylation at a single locus. Infinium I assays use two beads per CpG, one for each of the methylated and unmethylated states. If bisulfite converted DNA matches the probe, the probe is extended with a nucleotide  attached to a red or green dye. In contrast, the Infinium II design uses one bead type and the methylated state is determined at the single base extension step after hybridization (the methylated signal is measured in the green channel and the unmethylated signal in the red channel) [See Figure 1]. In 2016, the 850k array (also called EPIC array) was introduced. This array also uses a combination of the Infinium I and II assays but builds upon the 450k slide with >90% of the original 450K CpGs plus an additional 350,000 CpGs in mainly enhancer regions. As a result of this increase coverage a 450k slide can contain 12 arrays for 12 samples whilst the EPIC has only 8 spaces for 8 samples per array. The EPIC array is replacing the 450K array as the *de facto* standard for methylation analyses; the data processing for both is however fairly similar.


.. image:: Figures/Infinium.png
   :target: Figures/Infinium.png
   :alt: 
 
*Fig. 1: Infinium I and II design.*

Regardless of array type, both the 450k and EPIC record two measurements for each CpG: a methylated intensity (M) and an unmethylated intensity (U). Using these values, the proportion of methylation at each site CpG locus can be determined. The level of methylation at a locus is commonly reported as the Beta-value, *i.e.* the ratio of the methylated probe intensity and the overall intensity:

.. math::
   \beta = M/(M + U)

Illumina recommends adding a constant offset α (by default, α = 100) to the denominator to regularize Beta value when both methylated and unmethylated probe intensities are low. The Beta-value statistic results in a number between 0 and 1, or 0 and 100%. Under ideal conditions, a value of zero indicates that all copies of the CpG site in the sample were completely unmethylated (no methylated molecules were measured) and a value of one indicates that every copy of the site was methylated.

A second common metric to describe the methylation level is the M-value, *i.e* the log2 ratio of the intensities of methylated probe versus unmethylated probe:

.. math::
   Mvalue = log2(M/U)

A M-value close to 0 indicates a similar intensity between the methylated and unmethylated probes, which means the CpG site is about half-methylated, assuming that the intensity data has been properly normalized. Positive M-values mean that more molecules are methylated than unmethylated, while negative M-values mean the opposite. 

Beta and M-values are related to each other but Beta-values are generally preferable for the graphical representation of methylation levels as *percentage methylation* has a more intuitive biological interpretation. Due to their distributional properties, M-values are more statistically valid for the differential analysis of methylation levels. A thorough comparison of both metrics, can be found `here <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587>`_.


.. image:: Figures/Beta_M.png
   :target: Figures/Beta_M.png
   :alt: 

*Fig. 2: Relationship between Beta and M-values.*


Datasets
--------

To demonstrate the various aspects of analysing methylation data, we will be using a small, publicly available 450k methylation dataset (\ `GSE49667 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49667>`_). The dataset contains 10 samples in total: there are 4 different sorted T-cell types (naive, rTreg, act_naive, act_rTreg, collected from 3 different individuals: M28, M29, M30). Not all individuals contributed all 4 cell types, so there are 10 samples in total. An additional birth sample (individual VICS-72098-18-B) is included from another study (`GSE51180 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51180>`_) to illustrate approaches for identifying and excluding poor quality samples.

Set Up R environment
--------------------

This exercise has been set up to run on Uppmax, so connect to the server as described in :doc:`../setup/lab-setup`. Now, there are two options to set up the R environment. 

**Option A**

The easiest option makes use of the module system on Uppmax. This is the best way to avoid problems with dependencies between packages and avoids the issue of missing system libraries. Sometimes, this option suffers from slow response times when using Rstudio or has issues rendering figures. Becuse of the easy setup it might still be worth trying out this option first.

On Uppmax, most packages are already installed, and can be loaded into R after the *R/4.0.0* and  *R_packages/4.0.0* modules have been loaded. If you are running on Uppmax, start by loading the following modules:

.. code-block:: bash

   module load R/4.0.0
   module load R_packages/4.0.0
   module load RStudio

Start the analysis by initiating *RStudio*... This might take a few seconds and a :code:`libGL error` can be shown before loading the RStudio graphical interface.

.. code-block:: bash

   rstudio

.. note::

   If rstudio runs too slow, you could also decide to run the whole tutorial in the normal R terminal. Instead of ``rstudio`` do

   .. code-block:: bash

      R

   If you do this, you might have to take a few extra steps to show the graphics produced by R. First, check if plotting works by trying ``plot(1:10)`` in the R terminal. If you see the plot, you are good to start the tutorial. If not and you are on Mac; install and open `Xquartz  <https://www.xquartz.org>`_ on your Mac before ssh-ing to rackham. If you are on a PC, follow the instructions on `this website <https://uisapp2.iu.edu/confluence-prd/pages/viewpage.action?pageId=280461906>`_ (under the section "Alternate methods using OS-soecific tools").

Next, run the R commands by copying them from this website into the Rstudio terminal or R terminal and pressing *Enter*. 

**Option B**

Alternatively, we provide a containerized environment consisting of R, Rstudio and the necessary packages for this session. Containers are a relatively new method to package software together with all its dependencies and an operating system. This means the software can easily run within the container on almost any computer or server, greatly simplifying software installation and management. Containers will be discussed in a bit more detail on Thursday. A benefit of using it here is that Rstudio runs a whole lot faster using the container approach. However, to access it from Uppmax, a few more steps are necessary. First, make sure you are connected to your alloted node (described in :doc:`../setup/lab-setup`) and then perform following steps.

.. code-block:: bash

   # Run the startup script; this will start the container and run Rstudio
   sh /sw/courses/epigenomics/DNAmethylation/startup_script.sh

You should see something like this:

.. code-block:: bash

   1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:r37.uppmax.uu.se:35616 vincent@rackham.uppmax.uu.se
   
   and point your web browser to http://localhost:8787

   2. log in to RStudio Server using the following credentials:

   user: vincent
   password: epi2021

Now, open a second terminal and run **your** ssh command from 1. Then open your web browser (Safari, Chrome, ...) and go to http://localhost:8787. Here, fill in **your** user and password as in 2. and Rstudio will start.

**Load Libraries**

After setting up Rstudio with by either option A or B start by loading the set of R packages that will be needed during the analysis: *limma* provides the statistical framework for testing differential methylation. *minfi*\ , *missMethyl*\ , *minfiData* and *DMRcate* are packages developed to work with methylation data. *Gviz* and *RColorBrewer* provide functions for the visualization of the data.

.. code-block:: r

   # Set the correct library path
   .libPaths("/sw/apps/R_packages/4.0.0/rackham")
   # load packages required for analysis
   library("limma")
   library("minfi")
   library("RColorBrewer")
   library("missMethyl") # Can take a short time...
   library("minfiData")
   library("Gviz")
   library("DMRcate")
   library("DMRcatedata")
   library("stringr")
   library("mCSEA")

Included with *minfi* is the *IlluminaHumanMethylation450kanno.ilmn12.hg19* package; it contains all the annotation information for each of the CpG probes on the 450k array. This will be useful later to to determine where the differentially methylated probes (hereafter referred to as DMP) are located in a genomic context and to link the Red and Green raw data to methylated and unmethylated status.

.. code-block:: r

   ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   # Use the head command to get a quick overview of the data and see what types of annotations are available
   head(ann450k)

.. note::

   These packages are of course also available for the later array versions. The EPIC array annotation package is called *IlluminaHumanMethylationEPICanno.ilm10b2.hg19* and also included in *minfi*.

Load Datasets
-------------

The datasets have been uploaded to Uppmax prior to the workshop, so you just need to point R towards the directory they are saved. The ``list.files`` command will return the list of files in the specified directory.

.. code-block:: r

   dataDirectory <- "/sw/courses/epigenomics/DNAmethylation/array_data/"
   # list the files
   list.files(dataDirectory, recursive = TRUE)

Illumina methylation data is usually obtained in the form of Intensity Data (IDAT) Files. This is a proprietary format that is output by the slide scanner and stores the intensities for each probe on the array. Typically, each IDAT file is approximately 8MB in size. The simplest way to import the raw methylation data into R is using the minfi function ``read.metharray.sheet``\ , along with the path to the IDAT files and a sample sheet. The sample sheet is a CSV (comma-separated) file containing one line per sample, with a number of columns describing each sample. The format expected by the ``read.metharray.sheet`` function is based on the sample sheet file that usually accompanies Illumina methylation array data. It is also very similar to the targets file described by the limma package. Importing the sample sheet into R creates a dataframe with one row for each sample and several columns. The ``read.metharray.sheet`` function uses the specified path and other information from the sample sheet to create a column called Basename which specifies the location of each individual IDAT file in the experiment. Import the metadata and have a look at the different samples.

.. code-block:: r

   # read in the sample sheet for the experiment
   targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
   targets

Now we know where the data is located and we have essential information on each samples identity, we can read in the raw intensity data into R using the ``read.metharray.exp`` function. This creates an *RGChannelSet* object that contains all the raw intensity data, from both the red and green colour channels, for each of the samples. This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. Note that this object contains the intensities of the internal control probes as well. Because we read the data from a data sheet experiment, the phenotype data is also stored in the *RGChannelSet* and can be accessed via the accessor command ``pData``. Also the probed design can be summarized by querying this object. Before starting the actual analysis it is good practice to get a feel of the structure and content of the *RGChannelSet* object in this way.

.. code-block:: r

   # read in the raw data from the IDAT files; warnings can be ignored.
   rgSet <- read.metharray.exp(targets=targets)

   # Get an overview of the data
   rgSet
   pData(rgSet)
   getManifest(rgSet)

It might be useful to change the names of the samples into something a little more descriptive.

.. code-block:: r

   # give the samples descriptive names
   targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
   sampleNames(rgSet) <- targets$ID

   # Check the names have been updated by looking at the rownames of the phenoData
   pData(rgSet)

.. note::

   If you prefer to run this tutorial locally, you can also download the dataset to your personal computer. To do this, navigate to the folder on your own conputer where you want to deposit the data and execute :code:`scp -r <username>@rackham.uppmax.uu.se:/sw/courses/epigenomics/DNAmethylation/array_data .`. Then you can point the :code:`dataDirectory` to this local directory. Of course, you will also have to install all packages locally!

A Note on Class Structure
-------------------------

minfi generates a number of classes corresponding to various transformations of the raw data. It is important to understand how these classes relate to each other. Figure 2 provides a useful overview. In a first step, IDAT files are collected in a *RGChannelSet* object, transformed in a *MethylSet* through a preprocess function and via two functions *ratioConvert* and *mapToGenome* (order does not matter) converted into an analysis-ready *GenomicRatioSet*.


.. image:: Figures/Classes.png
   :target: Figures/Classes.png
   :alt: 
   
*Fig. 2: Flowchart of the different *minfi* class objects.*

As of now, our dataset is an *RGChannelSet* object containing the raw green and red intensity data. To proceed, this needs to be transformed into a *MethylSet* object containing the methylated and unmethylated signals. The most basic way to construct a *MethylSet* is to use the function *preprocessRaw* which uses the array design to match up the different probes and color channels to construct the methylated and unmethylated signals. This function does not do any normalization (in a later step we will add normalization, but this step is useful for initial quality control). Do this now for your object and have a look at the changes in the metadata. Notice that the red and green assays have been transformed in Meth and Unmeth signals.

.. code-block:: r

   MSet <- preprocessRaw(rgSet)
   MSet
   # Compare to previous object
   rgSet

The accessors *getMeth* and *getUnmeth* can now be used on the *MethylSet* to get the methylated and unmethylated intensities matrices, if necessary.

.. code-block:: r

   head(getMeth(MSet)[,1:3])
   head(getUnmeth(MSet)[,1:3])

A *RatioSet* object is class designed to store Beta and/or M-values instead of the (un)methylated signals. An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. Mapping a *MethylSet* to a *RatioSet* is irreversible, i.e. one cannot technically retrieve the methylated and unmethylated signals from a *RatioSet*. A *RatioSet* can be created with the function ratioConvert. The function *mapToGenome* applied to a *RatioSet* object will add genomic coordinates to each probe together with some additional annotation information. The output object is a *GenomicRatioSet* 

.. code-block:: r

   ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
   # Observe the change of the assays
   ratioSet

   gset <- mapToGenome(ratioSet)
   gset

The functions *getBeta*\ , *getM* and *getCN* work on the *GenomicRatioSet* return respectively the Beta value matrix, M value matrix and a the Copy Number matrix.

.. code-block:: r

   beta <- getBeta(gset)
   head(beta)
   m <- getM(gset)
   head(m)
   cn <- getCN(gset)
   head(cn)

Much more annotation data can be extracted from this object (see the *minfi* `documentation <http://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html>`_\ ). Now we have a analysis ready object, albeit unnormalized. As we will see in a later section, there are several normalization options that automatically take care of the preprocessing and conversion of a *RGChannelSet* to a *GenomicRatioSet*. But before doing this, an important step is Quality Control

Quality control
---------------

*minfi* provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels. When plotting these two medians against each other, good samples tend to cluster together, while failed samples tend to separate and have lower median intensities. In general, users should make the plot and make a judgement. The line separating ”bad” from ”good” samples represent a useful cutoff, which is not always very clear and may have to be adapted to a specific dataset. The functions *getQC* and *plotQC)* are designed to extract and plot the quality control information from the *MethylSet*. 

.. code-block:: r

   qc <- getQC(MSet)
   plotQC(qc)

Here, the cutoff line suggests 3 "bad" samples. Can you determine which samples these are? 
   
.. hint:: 
   the *pData* function might be of help here. 
   
In general, a decision of good versus bad quality should be based on multiple metrics, not just one. Therefore, we can additionally look at the detection p-values for every CpG in every sample, which is indicative of the quality of the signal. The method used by *minfi* to calculate detection p-values compares the total signal (M+U) for each probe to the background signal level, which is estimated from the negative control probes. Very small p-values are indicative of a reliable signal whilst large p-values, for example >0.01, generally indicate a poor quality signal.

Plotting the mean detection p-value for each sample allows us to gauge the general quality of the samples in terms of the overall signal reliability. Samples that have many failed probes will have relatively large mean detection p-values.

.. code-block:: r

   # calculate the detection p-values
   detP <- detectionP(rgSet)
   head(detP)

These p-values can be summarized in a single plot to simplify the comparison between samples

.. code-block:: r

   # examine mean detection p-values across all samples to identify any failed samples
   barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
   abline(h=0.05,col="red")

Poor quality samples can be easily excluded from the analysis using a detection p-value cutoff, for example >0.05. For this particular dataset, the *birth* sample shows a very high mean detection p-value.

The overall density distribution of Beta values for each sample is another useful metric to determine sample quality. Usually, one would expect to see most Beta values to be either close to 0 or 1, indicating most of the CpG sites in the sample are unmethylated or methylated. The *densityPlot* function plots these distribution for each sample.

.. code-block:: r

   phenoData <- pData(MSet)
   densityPlot(MSet, sampGroups = phenoData$Sample_Group)

The 450k array contains several internal control probes that can be used to assess the quality control of different sample preparation steps (bisulfite conversion, hybridization, etc.). The values of these control probes are stored in the initial *RGChannelSet* and can be plotted by using the function *controlStripPlot* and by specifying the control probe type. We will not go into the details of each control probe type, but these might be useful to determine the exact reason a sample failed QC.

.. code-block:: r

   controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")
   # The plots of the different control probes can be exported into a pdf file in one step using the function qcReport
   #qcReport(rgSet, pdf= "qcReport.pdf")

Taking these different metrics into account, it seems clear that the *birth* sample is of lower quality than the other samples. Therefore, we can decide to exclude it from the initial *rgSet* prior to further analysis.

.. code-block:: r

   # select the samples to keep for further analysis
   keep <- !colnames(rgSet) == "birth.11"
   # subset rgSet
   rgSet <- rgSet[,keep]
   # Check the sample has been removed by looking at the number of colnames
   rgSet
   # subset target as well
   targets <- targets[keep,]

Normalization
-------------

So far, we did not perform any normalization to process the data. Due to the intrinsic chip design of 2 types of chemistry probes, data normalization is an important step to think about before data analysis. Given the higher dynamic range of type I probes, one expects that  - when left uncorrected - there would be a relative overenrichment of type I over type II probes in a top ranked list of probes correlating with a phenotype. So, if you are comparing probes within an array or are ranking them for example, a within-array normalization might be necessary. 

Additionally, there is often systematic bias between arrays due to a variety of variable experimental conditions such as concentrations of reagents or temperature, especially when the experiments are carried out in several batches. Relevant biological signals may be masked by technical differences, also called batch effects and there are two fundamental ways to deal with them. One possibility is to consider batch effects in the statistical analysis, for instance by introducing a dummy variable for the batch in a linear model. However, batch effects may alter the data in complicated ways for which the statistical model in mind may not be adequate. It might therefore be preferable to remove these technical differences in a preprocessing step. 
 
Several distinct preprocessing and normalization procedures are therefore available in *minfi* (see below). A choice of different options raise of course the question which one is best or most optimal for your particular dataset. This is a difficult question to answer beforehand and selecting the best option is in practice often an iterative procedure while looking at the distribution of the Beta values (see example of different methods in Figure 4). Nevertheless, there are some general guidelines and the authors of *minfi* have the following to say about this:

.. note::

    "Many people have asked us which normalization they should apply to their dataset. Our rule of thumb is the following. If there exist global biological methylation differences between your samples, as for instance a dataset with cancer and normal samples, or a dataset with different tissues/cell types, use the preprocessFunnorm function as it is aimed for such datasets. On the other hand, if you do not expect global differences between your samples, for instance a blood dataset, or one-tissue dataset, use the preprocessQuantile function. In our experience, these two normalization procedures perform always better than the functions preprocessRaw, preprocessIllumina and preprocessSWAN discussed below. For convenience, these functions are still implemented in the minfi package."

So, try different methods and compare the normalized data. Do the Beta values of the different probes or different samples look more comparable after normalization?

.. image:: Figures/norms.jpg
   :target: Figures/norms.jpg
   :alt: 
   
*Fig. 4: (A) No normalization. (B) Lumi-based classical quantile normalization. (C) Peak-based correction followed by quantile normalization. (D) Subset quantile normalization with a unique set of reference quantiles computed from Infinium I signals. (E) Subset quantile normalization with a reference quantiles set computed from Infinium I signals for each kind of probe category according to the ‘relation to CpG’ annotations provided by Illumina (CA, USA). (F) Subset quantile normalization with a reference quantiles set computed from Infinium I signals for each kind of probe category. NT: Density plot of the median β-value profile for nontumoral samples; T: Density plot of the median β-value profile for tumoral samples.*

Below a short overview of the normalization methods included in *minfi*.

preprocessRaw
^^^^^^^^^^^^^

As seen before, this function converts a *RGChannelSet* to a *MethylSet* by converting the Red and Green channels into a matrix of methylated signals and a matrix of unmethylated signals. No normalization is performed.

.. attention::
   | Input: *RGChannelSet* 
   | Output: *MethylSet*

preprocessIllumina
^^^^^^^^^^^^^^^^^^

Convert a *RGChannelSet* to a *MethylSet* by implementing the preprocessing choices as available in Genome Studio: background subtraction and control normalization. Both of them are optional and turning them off is equivalent to raw preprocessing (\ *preprocessRaw*\ ):

.. attention::
   | Input: *RGChannelSet* 
   | Output: *MethylSet*

preprocessSWAN
^^^^^^^^^^^^^^

Perform Subset-quantile within array normalization (SWAN), a within-array normalization correction for the technical differences between the Type I and Type II array designs. The algorithm matches the Beta-value distributions of the Type I and Type II probes by applying a within-array quantile normalization separately for different subsets of probes (divided by CpG content). The input of SWAN is a *MethylSet*\ , and the function returns a *MethylSet* as well. If an *RGChannelSet* is provided instead, the function will first call *preprocessRaw* on the *RGChannelSet*\ , and then apply the SWAN normalization. 

.. attention::
   | Input: *RGChannelSet* or *MethylSet* 
   | Output: *MethylSet*

preprocessFunnorm
^^^^^^^^^^^^^^^^^

The function *preprocessFunnorm* uses the internal control probes present on the array to infer between-array technical variation. It is particularly useful for studies comparing conditions with known large-scale differences, such as cancer/normal studies, or between-tissue studies. It has been shown that for such studies, functional normalization outperforms other existing approaches. By default, is uses the first two principal components of the control probes to infer the unwanted variation.

.. attention::
   | Input: *RGChannelSet*
   | Output: *GenomicRatioSet*

preprocessQuantile
^^^^^^^^^^^^^^^^^^

This function implements stratified `quantile normalization <https://en.wikipedia.org/wiki/Quantile_normalization>`_ preprocessing. The normalization procedure is applied to the Meth and Unmeth intensities separately. The distribution of type I and type II signals is forced to be the same by first quantile normalizing the type II probes across samples and then interpolating a reference distribution to which we normalize the type I probes. Since probe types and probe regions are confounded and we know that DNA methylation varies across regions we stratify the probes by region before applying this interpolation. Note that this algorithm relies on the assumptions necessary for quantile normalization to be applicable and thus is not recommended for cases where global changes are expected such as in cancer-normal comparisons as these would be removed by the normalization. 

.. attention::
   | Input: *RGChannelSet* 
   | Output: *GenomicRatioSet*

As we are comparing different blood cell types, which are globally relatively similar, we will apply the preprocessQuantile method to our data. 

.. warning::
   This assumption might not be true; in an actual analysis it would be advised to try and compare different normalization methods. 

Note that after normalisation, the data is housed in a GenomicRatioSet object; automatically running the steps we did manually to do an initial quality control. 

.. code-block:: r

   # normalize the data; this results in a GenomicRatioSet object
   mSetSq <- preprocessQuantile(rgSet)


Compare with the unnormalized data to visualize the effect of the normalization. First a comparison of the Beta distributions for the different probe designs. This will give an indication of the effectiveness of the within-array normalization.

.. code-block:: r

   par(mfrow=c(1,2))
   # Plot distributions prior to normalization for sample 1
   plotBetasByType(MSet[,1],main="Raw")
   # The normalized object is a GenomicRatioSet which does not contain
   # the necessary probe info, we need to extract this from the MethylSet first.
   typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
   typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
   probeTypes <- rbind(typeI, typeII)
   probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
   # Now plot the distributions of the normalized data for sample 1
   plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized",)
   

Does it look like the normalization brought the distributions closer to each other? Now let's see how the between-array normalization worked...

.. code-block:: r

   # visualise what the data looks like before and after normalisation
   par(mfrow=c(1,2))
   densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
   legend("top", legend = levels(factor(targets$Sample_Group)), 
          text.col=brewer.pal(8,"Dark2"))
   densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
               main="Normalized", legend=FALSE)
   legend("top", legend = levels(factor(targets$Sample_Group)), 
          text.col=brewer.pal(8,"Dark2"))
   

.. hint::
   Click on Zoom above the RStudio plot panel to watch a larger version of the plotted figure.

Data exploration
----------------

After normalization of your data is a good time to look at the similarities and differences between the various samples. One way to do this is by creating a MDS or Multi-Dimenional Scaling plot. This is a method to graphically represent relationships between objects (here the different samples) in multidimensional space onto 2 or 3 dimensional space. Dimension one (or principal component one) captures the greatest source of variation in the data, dimension two captures the second greatest source of variation in the data and so on. Colouring the data points or labels by known factors of interest can often highlight exactly what the greatest sources of variation are in the data. In a good quality dataset, one would hope that biological differences would show up as one of the greatest sources of variation. It is also possible to use MDS plots to decipher sample mix-ups. The following code creates the MDS plot twice but the samples in the left plot are colored according to celltype, while the plot on the right is colored according to "individual". Before you proceed think a moment about what this figure tells you about the sources in variation in the data. Try changing the ``dim=c(1,2)`` parameter to for example ``dim=c(1,3)`` or other values to get an even deeper understanding of the variation in the data. 

.. code-block:: r

   # MDS plots to look at largest sources of variation
   # Create color panel
   pal <- brewer.pal(8,"Dark2")
   # Plot figures
   par(mfrow=c(1,2))
   plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
           col=pal[factor(targets$Sample_Group)], dim=c(1,2))
   legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
          bg="white", cex=0.7)

   plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
           col=pal[factor(targets$Sample_Source)], dim=c(1,2))
   legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
          bg="white", cex=0.7)

Examining the MDS plots for this dataset demonstrates that the largest source of variation is the difference between individuals. The higher dimensions reveal that the differences between cell types are largely captured by the third and fourth principal components. This type of information is useful in that it can inform downstream analysis. If obvious sources of unwanted variation are revealed by the MDS plots, we can include them in our statistical model to account for them. In the case of this particular dataset, we will include individual to individual variation in our statistical model.

Filtering
---------

Poor performing probes can obscure the biological signals in the data and are generally filtered out prior to differential methylation analysis. As the signal from these probes is unreliable, by removing them we perform fewer statistical tests and thus lower the multiple testing penalty. We filter out probes that have failed in one or more samples based on detection p-value.

.. code-block:: r

   # ensure probes are in the same order in the mSetSq and detP objects
   detP <- detectionP(rgSet)
   detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

   # remove any probes that have failed in one or more samples; this next line
   # checks for each row of detP whether the number of values < 0.01 is equal 
   # to the number of samples (TRUE) or not (FALSE)
   keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
   table(keep)
   # Subset the GenomicRatioSet
   mSetSqFlt <- mSetSq[keep,]
   mSetSqFlt

Because the presence of short nucleotide polymorphisms (or SNPs) inside the probe body or at the nucleotide extension can have important consequences on the downstream analysis, *minfi* offers the possibility to remove such probes. 

.. note::
   Can you see why SNP can be a problem in methylation data analysis (Hint: C to T conversions are the most common type of SNP in the human genome)? 

There is a function in *minfi* that provides a simple interface for the removal of probes where common SNPs may affect the CpG. You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.

.. code-block:: r

   mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
   mSetSqFlt

Once the data has been filtered and normalised, it is often useful to re-examine the MDS plots to see if the relationship between the samples has changed. From the new MDS plots it is apparent that much of the inter-individual variation has been removed as this is no longer the first principal component, likely due to the removal of the SNP-affected CpG probes. However, the samples do still cluster by individual in the second dimension and thus a factor for individual should still be included in the model.

.. code-block:: r

   par(mfrow=c(1,2))
   plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
           col=pal[factor(targets$Sample_Group)], cex=0.8)
   legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
          cex=0.65, bg="white")

   plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
           col=pal[factor(targets$Sample_Source)])
   legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
          cex=0.7, bg="white")
   # Close double plotting window
   dev.off()

Probe-Wise Differential Methylation
-----------------------------------

After all this preprocessing and filtering, the time has come to address the actual biological question of interest! Namely, which CpG sites are differentially differentially methylated between the different cell types? To do this, we will design a linear model in *limma*.

As was apparent from the MDS plots, there is an additional factor that we need to take into account when performing the statistical analysis needed to solve this question. In the targets file, there is a column called Sample_Source, which refers to the individuals that the samples were collected from. Hence, when we specify our design matrix, we need to include two factors: individual and cell type. This style of analysis is called a paired analysis; differences between cell types are calculated within each individual, and then these differences are averaged across individuals to determine whether there is an overall significant difference in the mean methylation level for each CpG site. 

.. warning::
   This design is fit for this dataset, and this dataset only. For future analyses, you will have to adapt the analysis style and design to your particular dataset. The `limma User’s Guide <https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf>`_ extensively covers the different types of designs that are commonly used for microarray experiments and how to analyse them in R.

.. code-block:: r

   # calculate M-values for statistical analysis: as previously mentioned, M-values have nicer statistical properties and are thus better for use in statistical analysis of methylation data 
   mVals <- getM(mSetSqFlt)

   # Set up the design matrix for the Differential Methylation analysis
   # Define the factor of interest
   cellType <- factor(targets$Sample_Group)
   # Define is the individual effect that we need to account for
   individual <- factor(targets$Sample_Source) 
   # use the above to create a design matrix
   design <- model.matrix(~0+cellType+individual, data=targets)
   colnames(design) <- c(levels(cellType),levels(individual)[-1])

   # fit the actual linear model to the data
   fit <- lmFit(mVals, design)

We are interested in pairwise comparisons between the four cell types, taking into account variation between individuals. We perform this analysis on the matrix of M-values in *limma*\ , obtaining t-statistics and associated p-values for each CpG site. A convenient way to set up the model when the user has many comparisons of interest that they would like to test is to use a contrasts matrix in conjunction with the design matrix. A contrasts matrix will take linear combinations of the columns of the design matrix corresponding to the comparisons of interest, essentially subsetting the data to these comparisons.

.. code-block:: r

   # create a contrast matrix for specific comparisons
   contMatrix <- makeContrasts(naive-rTreg,
                              naive-act_naive,
                              rTreg-act_rTreg,
                              act_naive-act_rTreg,
                              levels=design)
   contMatrix

Next, these contrasts are fitted to the model and the statistics and p-values of differential expression are calculated by the function *eBayes*. this function is used to rank genes in order of evidence for differential methylation. We will not delve too deep into the background of this statistical testing framework; if you are interested in this more info can be found `here <Linear models and empirical bayes methods for assessing differential expr…>`_. 

.. code-block:: r

   # fit the contrasts
   fit2 <- contrasts.fit(fit, contMatrix)
   # Rank genes
   fit2 <- eBayes(fit2)

Using the *topTable* function in *limma*\ , the differentially methylated genes per comparison/contrast can be extracted. To order these by p-value, the user can specify sort.by="p". The results of the analysis for the first comparison, naive vs. rTreg, can be saved as a data.frame by setting *coef=1*. The *coef* parameter explicitly refers to the column in the contrasts matrix which corresponds to the comparison of interest.

.. code-block:: r

   # get the table of results for the first contrast (naive - rTreg)
   DMPs <- topTable(fit2, num=Inf, coef=1)
   head(DMPs)

We can add a bit more annotation to this list of CpGs, by adding a *genelist* parameter to the *topTable* function. This can be useful to retrieve the location of the CpG, the nearest gene or CpG island and other information.

.. code-block:: r

   # Retrieve data from the array annotation package; this is array-specific
   ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                         c(1:4,12:19,24:ncol(ann450k))]
   DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
   head(DMPs)

   # The resulting data.frame can easily be written to a CSV file, which can be opened in Excel.
   # write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

It is always a good idea to plot the most differentially methylated sites as a quick sanity check; if the plot does not make sense there might have been an issue with the model design or setup of the contrast matrix. To do this, we first extract the Beta-values (remember these are the preferential values to visualize).

.. code-block:: r

   # eXtract Beta-values
   bVals <- getBeta(mSetSqFlt)

   # Plot most significant differentially methylated CpG
   plotCpg(bVals, cpg="cg07499259", pheno=targets$Sample_Group, ylab = "Beta values")

Does this plot makes sense? Are the differences in methylation percentage as expected? 

Regional Differential Methylation (DMR)
---------------------------------------

Location-based Regions
^^^^^^^^^^^^^^^^^^^^^^

Often, differential methylation of a single CpG is not so informative or can be hard to detect. Therefore, knowing whether several CpGs near to each other (or *regions*\ ) are concordantly differentially methylated can be of greater interest.

There are several Bioconductor packages that have functions for identifying differentially methylated regions from 450k data. Some of the most popular are the *dmrFind* function in the *charm* package, which has been somewhat superseded for 450k arrays by the *bumphunter* function in *minfi*\ , and, the *dmrcate* in the *DMRcate* package. They are each based on different statistical methods, but we will be using *dmrcate* here, as it is based on *limma* and thus we can use the design and contrast matrix we defined earlier.

We will again start from our matrix of M-values. For this kind of analysis, this matrix has to be annotated with the chromosomal position of the CpGs and their gene annotations. Because in a first step the *limma* differential methylation analysis for single CpGs will be run again, we need to specify the design matrix, contrast matrix and contrast of interest. 

.. note::
   More info on the different options can always be found in the manual; *i.e* by using *?cpg.annotate* in R.

.. code-block:: r

   myAnnotation <- cpg.annotate(object = mVals, 
                                datatype = "array", 
                                what = "M", 
                                analysis.type = "differential", 
                                design = design, 
                                contrasts = TRUE, 
                                cont.matrix = contMatrix, 
                                coef = "naive - rTreg", 
                                arraytype = "450K")
   myAnnotation

Once we have the relevant statistics for the individual CpGs, we can then use the *dmrcate* function to combine them to identify differentially methylated regions. Of particular interest here is the *lambda* parameter; this value is the number of nucleotides that is allowed between significant CpGs before splitting them up in different regions. So a smaller *lambda* will result in more but smaller regions. For array data, the authors of the *dmrcate* package currently recommend a lambda of 1000. The main output table DMRs contains all of the regions found, along with their genomic annotations and p-values. To inspect this object and further visualization, you can best use the *extractRanges* function to create a *GRanges* object.

.. code-block:: r

   DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
   DMRs
   # Create GRanges object; create directory when prompted
   results.ranges <- extractRanges(DMRs)
   results.ranges

Just as for the single CpG analysis, it is a good idea to visually inspect the results to make sure they make sense. For this, use the *DMR.plot* function. By default, this plot draws the location of the DMR in the genome, the position of nearby genes, the positions of the CpG probes, the Beta value levels of each sample as a heatmap and the mean methylation levels for the various sample groups in the experiment.

.. code-block:: r

   # set up the grouping variables and colours
   pal <- brewer.pal(8,"Dark2")
   groups <- pal[1:length(unique(targets$Sample_Group))]
   names(groups) <- levels(factor(targets$Sample_Group))
   cols <- groups[as.character(factor(targets$Sample_Group))]
   # draw the plot for the second DMR - first gives error for some reason...
   DMR.plot(ranges = results.ranges, 
            dmr = 2, 
            CpGs = mSetSqFlt,
            phen.col = cols, 
            genome = "hg19")

Interestingly, the hypomethylation of the second DMR, near TIGIT, in Treg was  one of the main conclusions of the paper base don this dataset:  

.. note::
   ...In support of the view that methylation limits access of FOXP3 to its DNA targets, we showed that increased expression of the immune suppressive receptor T-cell immunoglobulin and immunoreceptor tyrosine-based inhibitory motif domain (TIGIT), which delineated Treg from activated effector T cells, was associated with hypomethylation and FOXP3 binding at the TIGIT locus... 

Functional Regions
^^^^^^^^^^^^^^^^^^

An alternative approach to detect DMRs is to predefine the regions to be tested; so, as opposed to the previous approach where the regions are defined according to heuristic distance rules we can define regions based on a shared function. For this, we will used the package *mCSEA* which contains three types of regions for 450K and EPIC arrays: promoter regions, gene body and CpG Islands. *mCSEA* is based on Gene Set Enrichment analysis (GSEA), a popular methodology for functional analysis that was specifically designed to avoid some drawbacks in the field of gene expression. Briefly, CpG sites are ranked according to a metric (logFC, t-statistic, ...) and an enrichment score (ES) is calculated for each region. This is done by running through the entire ranked CpG list, increasing the score when a CpG in the region is encountered and decreasing the score when the gene encountered is not in the region. A high ES indicates these probes are found high up in the ranked list. In other words, a high (N)ES value means that for the CpG sites in this region there is - on average - a shift towards a higher methylation level. This approach has been `shown <https://academic.oup.com/bioinformatics/article/35/18/3257/5316232>`_ to be more effective to detect smaller but consistent methylation differences.

Here, we will apply this method to the output of the "naive-rTreg" comparison, ranking the CpGs by logFC differences. We specify "promoters" as the type of regions to be considered, but other options such as CpG Islands or gene bodies are possible. 

.. note::
   "Promoters" are not really restricted to pure promoters, but also include UTR, 1st Exon and a region upstream of the TSS.

.. code-block:: r

   # Create a named vector containing the rank metric (here: logFC)
   myRank <- DMPs$logFC
   names(myRank) <- rownames(DMPs)

   # Reshape the phenotype data to a format suitable for mCSEA
   pheno <- as.data.frame(pData(mSetSqFlt))
   pheno <- pheno[,"Sample_Group", drop=FALSE]

   # Run the mCSEA 
   myResults <- mCSEATest(myRank, 
                          bVals, 
                          pheno,
                          regionsTypes = "promoters", 
                          platform = "450k")
   head(myResults$promoters)

The main results are found in *myResults$promoters*. This data frame contains the (normalized) enrichment score, p-values, total number of associated CpGs and the leading edge CpGs. The leading edge CpGs are the real drivers of the ES; these can be considered the most important CpGs with the largest logFC.
The results of selected results can be visualized using *mCSEAPlot*\ , by specifying the *regionType* and the *dmrName*. Here an example of the second hit of the DMRs based on location; the promoter of TIGIT. Note that the gene name indicates the promoter of said gene, since we specified we only consider promoter regions in this analysis. The result of this visualization are the chromosomal location, Beta levels per CpG per sample, leading edge status (green if in leading edge set) and gene annotation.

.. code-block:: r

    mCSEAPlot(myResults, 
              regionType = "promoters", 
              dmrName = "TIGIT",
              transcriptAnnotation = "symbol", 
              makePDF = FALSE)

Gene Ontology Testing
---------------------

After obtaining a - potentially long - list of significantly differentially methylated CpG sites, one might wonder whether there is a (or multiple) specific biological pathway(s) over-represented in this list. In some cases it is relatively straightforward to link the top differentially methylated CpGs to genes that make biological sense in terms of the cell types or samples being studied, but there may be many thousands of CpGs significantly differentially methylated. Gene-set analysis (GSA) is frequently used to discover meaningful biological patterns from lists of genes generated from high-throughput experiments, including genome-wide DNA methylation studies. The objective is typically to identify similarities between the genes, with respect to annotations available from sources such as the Gene Ontology (GO) or Kyoto Encyclopedia of Genes and Genomes (KEGG).

We can perform this type of analysis using the *gometh* function in the *missMethyl* package. This function takes as input a character vector of the names (e.g. cg20832020) of the significant CpG sites, and optionally, a character vector of all CpGs tested. This is recommended particularly if extensive filtering of the CpGs has been performed prior to analysis as it constitutes the "background" out of which any significant CpG could be chosen. For gene ontology testing, the user can specify collection="GO” (which is the default option). For testing KEGG pathways, specify collection="KEGG”. In this tutorial, we will continue with the results from the single-probe "naive vs rTreg" comparison and select all CpG sites that have an adjusted p-value of less than 0.05.

.. code-block:: r

   # Get the significant CpG sites at less than 5% FDR
   sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
   # First 10 significant CpGs
   sigCpGs[1:10]
   # Total number of significant CpGs at 5% FDR
   length(sigCpGs)
   # Get all the CpG sites used in the analysis to form the background
   all <- DMPs$Name
   # Total number of CpG sites tested
   length(all)

..

.. warning::
   A key assumption of GSA methods is that all genes have, *a priori*\ , the same probability of appearing in the list of significant genes. If this is not true, that is, if certain genes are more likely to appear in the list, regardless of the treatments or conditions being investigated, this has the potential to cause misleading results from GSA. This has been `shown <https://academic.oup.com/bioinformatics/article/29/15/1851/265573>`_ to be a major source of bias in genome-wide methylation gene set analysis. Essentially it comes down to this: genes that have more CpGs associated with them will have a much higher probability of being identified as differentially methylated compared to genes with fewer CpGs. As a result gene sets containing many "highly covered" genes will be found to be significantly enriched much easier than other gene sets, regardless of the treatment or condition. For the 450k array, the numbers of CpGs mapping to genes can vary from as few as 1 to as many as 1200. The *gometh* function takes into account the varying numbers of CpGs associated with each gene on the Illumina methylation arrays. If you want to try alternative methods, keep in mind to check how they handle this source of bias. 


After having defined the significant and background sites, it is time to run the enrichment analysis itself.

.. code-block:: r

   # Run enrichment - Can take a bit of time...
   gst <- gometh(sig.cpg=sigCpGs, all.cpg=all)
   # Top 10 GO categories
   topGSA(gst, number=10)

Can you find the top 10 KEGG pathways? Do they make sense biologically?

While gene set testing is useful for providing some biological insight in terms of what pathways might be affected by abberant methylation, care should be taken not to over-interpret the results. Gene set testing should be used for the purpose of providing some biological insight that ideally would be tested and validated in further laboratory experiments. It is important to keep in mind that we are not observing gene level activity such as in RNA-Seq experiments, and that we have had to take an extra step to associate CpGs with genes.

Cell Type Composition
---------------------

As methylation is cell type specific and methylation arrays provide CpG methylation values for a population of cells, biological findings from samples that are comprised of a mixture of cell types, such as blood, can be confounded with cell type composition. In order to estimate the confounding levels between phenotype and cell type composition, the function *estimateCellCounts* (depending on the package *FlowSorted.Blood.450k*\ ) can be used to estimate the cell type composition of blood samples. The function takes as input a *RGChannelSet* and returns a cell counts vector for each samples. If there seems to be a large difference in cell type composition in the different levels of the phenotype, it might be needed to include the celltype proportions in the *limma* model to account for this confounding. Since we have been working with sorted populations of cells, this was not necessary for our data.

Alternative Workflows
---------------------

`RnBeads <https://rnbeads.org>`_ 
   R-based and user-friendly; includes modules for data import, quality control, filtering and normalization (“preprocessing”), export of processed data (“tracks and tables”), covariate inference (e.g., predicting epigenetic age and cell type heterogeneity from DNA methylation data), exploratory analysis (e.g., dimension reduction, global distribution of DNA methylation levels, hierarchical clustering), and differential DNA methylation analysis. Each analysis module generates an HTML report that combines method descriptions, results tables, and publication-grade plots. These reports provide the user with a comprehensive and readily sharable summary of the dataset.

`COHCAP <https://www.bioconductor.org/packages/release/bioc/html/COHCAP.html>`_ 
   R-based; provides a pipeline to analyze single-nucleotide resolution methylation data (Illumina 450k/EPIC methylation array, targeted BS-Seq, etc.). It provides differential methylation for CpG Sites, differential methylation for CpG Islands, integration with gene expression data, with visualizaton options. 
