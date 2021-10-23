.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

====================================================
Detection of differential binding sites using csaw
====================================================

This is an alternative workflow for detection of differential binding / occupancy in ChIP-seq data. In contrast to working with reads counted within peaks detected in a peak calling step (as in the earlier example with ``DiffBind``), this approach uses a **sliding window** to count reads across the genome. Each window is then tested for significant differences between libraries from different conditions, using the methods in the ``edgeR`` package. This package also offers an FDR control strategy more appropriate for ChIP-seq experiments than simple BH adjustment.

It can be used for point-source binding (TFs) as well as for broad signal (histones). However, it can only be used for cases where **global occupancy levels are unchanged**.

As this method is agnostic to signal structure, it requires careful choice of strategies for filtering and normalisation. Here, we show a very simple workflow. More details can be found in the `Csaw User Guide <https://bioconductor.org/packages/3.12/workflows/vignettes/csawUsersGuide/inst/doc/csaw.pdf>`_.

Requirements Local
====================



* ``R version 4.0.4 (2021-02-15) -- "Lost Library Book"``

* ``csaw``

* ``edgeR``

R packages required for annotation:

* ``org.Hs.eg.db``

* ``TxDb.Hsapiens.UCSC.hg19.knownGene``

Recommended (local):

* R-Studio to work in




.. NOTE::
  
  This exercise was tested on Rackham using pre-installed R libraries. Local installation of recommended R packages may require additional software dependecies.



.. Getting the data
.. ------------------

.. We will examine differences in REST binding in two cell types: SKNSH and HeLa, subset to chromosome 1. To run the tutorial locally need to download required files. Let's use the Box links for simplicity. 

.. HeLa:

.. * [zip](https://stockholmuniversity.box.com/s/2o3lchp61kzxpil1y1snn4onjk4e0sjo)
.. * [tar.gz](https://stockholmuniversity.box.com/s/wmx4uhgo3esuessr4f8g9kmvarm42bxz)

.. SKNSH:

.. * [zip](https://stockholmuniversity.box.com/s/dkurmi5suwh3qnxnx0ysfhh5c0g2d1ti)
.. * [tar.gz](https://stockholmuniversity.box.com/s/8m3rgtakx8h8rmhnwltitqccbx7h0wy3)


.. .. HINT::
  
..   You can also ``scp -r`` the files from rackham at ``/proj/g2020022/chipseq_proc/data/bam/``



.. To extract ``tar.gz`` files 


.. .. code-block:: bash

..   tar -zxvf archive_name.tar.gz




Requirements Remote (Uppmax)
--------------------------------

The software is configured.

To prepare the files, assuming you are in ``~/chipseq/analysis``:

.. code-block:: bash
  
   mkdir csaw
   cd csaw

   ln -s  ../../data/bam/hela/
   ln -s  ../../data/bam/sknsh/


Loading data and preparing the contrast to test
=================================================


You can now load the version of R for which we tested this class along with other dependencies:



.. code-block:: bash

   module load R_packages/4.0.4

The remaining part of the exercise is performed in ``R``.


.. .. HINT::

..   Modify the paths to folders with respective data to match your local setup:

..   Local:

..   .. code-block:: R

..     dir.sknsh = "/path/to/data/sknsh"
..     dir.hela = "/path/to/data/hela"


.. Remote:

.. code-block:: R

  dir.sknsh = "./sknsh"
  dir.hela = "./hela"
  hela.1=file.path(dir.hela,"ENCFF000PED.chr12.rmdup.sort.bam")
  hela.2=file.path(dir.hela,"ENCFF000PEE.chr12.rmdup.sort.bam")
  sknsh.1=file.path(dir.sknsh,"ENCFF000RAG.chr12.rmdup.sort.bam")
  sknsh.2=file.path(dir.sknsh,"ENCFF000RAH.chr12.rmdup.sort.bam")

  bam.files <- c(hela.1,hela.2,sknsh.1,sknsh.2)


We need to provide the information about the design of the experiment using ``model.matrix`` function:

.. code-block:: R

  grouping <- factor(c('hela', 'hela', 'sknsh', 'sknsh'))
  design <- model.matrix(~0 + grouping)
  colnames(design) <- levels(grouping)


The design should look like this:

.. code-block:: R

  > design
    hela sknsh
  1    1     0
  2    1     0
  3    0     1
  4    0     1
  attr(,"assign")
  [1] 1 1
  attr(,"contrasts")
  attr(,"contrasts")$grouping
  [1] "contr.treatment"


Let's test which peaks are differentially occupied in HeLa cells vs in sknsh cells.
We prepare the information on contrast to be tested using ``makeContrasts`` function from package ``limma``. This is not the only way to do so, and examples are given in ``csaw`` and ``edgeR`` manuals. In this case we want to test for the differences in REST binding in HeLa vs. SKNSH cell lines:

.. code-block:: R

  library(edgeR)
  contrast <- makeContrasts(hela - sknsh, levels=design)


The contrast should look like this

.. code-block:: R

  > contrast
         Contrasts
  Levels  hela - sknsh
    hela             1
    sknsh           -1


Now we are ready to load data and create an object with counted reads:

.. code-block:: R

  library(csaw)
  data <- windowCounts(bam.files, ext=100, width=10) 


Parameters for file loading can be modified (examples in the ``csaw`` User Guide), depending on how the data was processed. Here we explicitely input the value for fragment length as we have this information from the cross correlation analysis performed earlier during :doc:`ChIP-seq data processing tutorial <../chipseqProc/lab-chipseq-processing>`. It is 100 for Hela and 95 & 115 for sknsh. We've used 100 as it seems a reasonable averge value.


We can inspect the resulting ``data`` object, e.g.:

.. code-block:: R
  
  > data
  class: RangedSummarizedExperiment 
  dim: 155408 4 
  metadata(6): spacing width ... param final.ext
  assays(1): counts
  rownames: NULL
  rowData names(0):
  colnames: NULL
  colData names(4): bam.files totals ext rlen

  > data$totals
  [1] 1637778 2009932 2714033 4180463



Filtering out regions with very low coverage
===============================================

The next step is to filter out uninformative regions, i.e. windows with low read count, which represent background. There are many strategies to do it, depending on the biology of the experiment, IP efficiency and data processing. Here, we filter out lowest 99.9% of the windows, retaining the 0.1% windows with highest signal. The rationale is that for TF experiments only 0.1% of the genome can be bound, hence the remaining must represent background.

.. code-block:: R

  keep <- filterWindowsProportion(data)$filter > 0.999
  data.filt <- data[keep,]



To investigate the effectiveness of our filtering strategy:

.. code-block:: R

  > summary(keep)
     Mode   FALSE    TRUE 
  logical  145558    9850 


Normalisation
===============

Assigning reads into larger bins for normalisation:

.. code-block:: R
  
  binned <- windowCounts(bam.files, bin=TRUE, width=10000)

The TMM method trims away putative DB bins (i.e., those with extreme M-values) and computes normalization factors from the remainder to use in edgeR. The size of each library is scaled by the corresponding factor to obtain an effective library size for modelling.

Calculating the normalisation factors using a modified TMM method:

.. code-block:: R
  
  data.filt <- normFactors(binned, se.out=data.filt)


Inspecting the normalisation factors:

.. code-block:: R

  > data.filt$norm.factors
  [1] 0.9727458 1.0718693 0.9279702 1.0335341



Detecting differentially binding (DB) sites
============================================

This part of the procedure follows the logic developed for transcriptomics data in ``edgeR`` package. The steps are described in great detail in `csawBook <http://bioconductor.org/books/3.13/csawBook/chap-stats.html#setting-up-for-edger>`_.

Detecting DB windows:


.. code-block:: R

  data.filt.calc <- asDGEList(data.filt)
  data.filt.calc <- estimateDisp(data.filt.calc, design)

  >summary(data.filt.calc$trended.dispersion)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.2136  0.3565  0.4022  0.3984  0.4742  0.4892 

  fit <- glmQLFit(data.filt.calc, design, robust=TRUE)
  results <- glmQLFTest(fit, contrast=contrast)

  rowData(data.filt) <- cbind(rowData(data.filt), results$table)


Inspecting the results table:


.. code-block:: R

  > head(results$table)
      logFC   logCPM         F       PValue
  1 -7.239404 2.165639 17.229173 3.327018e-05
  2 -5.244217 2.783211  9.484909 2.074540e-03
  3 -3.023888 2.755437  4.721852 2.979352e-02
  4 -2.050617 2.612401  2.684560 1.013412e-01
  5 -1.827703 2.459979  2.459072 1.168638e-01
  6 -4.336717 2.052296 14.330442 1.538194e-04


Correcting for multiple testing
================================


First we merge adjacent DB windows into longer clusters. Windows that are less than ``tol`` apart are considered to be adjacent and are grouped into the same cluster. The chosen ``tol``
represents the minimum distance at which two binding events are treated as separate sites.
Large values (500 - 1000 bp) reduce redundancy and favor a region-based interpretation of
the results, while smaller values (< 200 bp) allow resolution of individual binding sites.

.. code-block:: R

  merged <- mergeWindows(rowRanges(data.filt), tol=1000L)


Next, we apply the multiple testing correction to obtain FDR. We combine p-values across clustered tests using Simes method to control the cluster FDR.

.. code-block:: R

  table.combined <- combineTests(merged$id, results$table)


The resulting ``table.combined`` object contains FDR for each cluster:

.. code-block:: R
  
  > head(table.combined)
  DataFrame with 6 rows and 8 columns
    num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
    <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
  1         7            0              5 2.32891e-04 0.004039711        down
  2         3            0              3 6.98933e-06 0.000482289        down
  3         3            0              3 1.94804e-04 0.003679925        down
  4         5            0              5 4.10817e-05 0.001168075        down
  5         3            0              3 6.67458e-05 0.001720419        down
  6         5            0              5 1.88055e-04 0.003620735        down
     rep.test rep.logFC
    <integer> <numeric>
  1         1  -7.23940
  2         8  -7.00091
  3        13  -7.50151
  4        14  -7.12133
  5        19  -7.20842
  6        23  -8.90988


* ``num.tests`` - the total number of windows in each cluster;
* fields ``num.up.logFC`` and ``num.down.logFC`` - for each log-FC column in ``results$table``; contain the number of windows with log-FCs above 0.5 or below -0.5, respectively;
* ``PValue`` - the combined p value;
* ``FDR`` - the q-value corresponding to the combined p value;
* ``direction`` - the dominant direction of change for windows in each cluster.


Each combined p value represents evidence against the global null hypothesis,
i.e., all individual nulls are true in each cluster. This may be more relevant than examining each
test individually when multiple tests in a cluster represent parts of the same underlying event, i.e.,
genomic regions consisting of clusters of windows. The BH method is then applied to control the
FDR across all clusters.


Inspecting the results
=========================

We select statistically significant DB events at FDR 0.05:

.. code-block:: R
  
  is.sig.region <- table.combined$FDR <= 0.05
  table(table.combined$direction[is.sig.region])


How many regions were detected as differentialy bound?

.. code-block:: R

    down   up 
     231  201 

out of

.. code-block:: R

  > length(table.combined$FDR)
  [1] 2758



We can also obtain information on the best window in each cluster:

.. code-block:: R

    tab.best <- getBestTest(merged$id, results$table)

  > head(tab.best)
  DataFrame with 6 rows and 8 columns
  num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
  <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
  1         7            0              4 2.32891e-04 0.004310830        down
  2         3            0              3 6.98933e-06 0.000529342        down
  3         3            0              3 2.56536e-04 0.004535416        down
  4         5            0              5 4.10817e-05 0.001168075        down
  5         3            0              3 6.67458e-05 0.001720419        down
  6         5            0              5 6.41895e-04 0.008158277        down
     rep.test rep.logFC
    <integer> <numeric>
  1         1  -7.23940
  2         8  -7.00091
  3        11  -7.33950
  4        14  -7.12133
  5        19  -7.20842
  6        22  -7.47709



We can inspect congruency of the replicates on multi-dimensional scaling (MDS) plots. 
The distance between each pair of libraries is computed as the square root of the mean squared log-fold change across the top set of bins with the highest absolute log-fold changes.
A small top set visualizes the most extreme differences whereas a large set visualizes overall differences.

.. code-block:: R

  par(mfrow=c(2,2))
  adj.counts <- cpm(data.filt.calc, log=TRUE)
  for (top in c(100, 500, 1000, 5000)) {
  plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"),labels=c("hela", "hela", "sknsh", "sknsh"), top=top)
  }


Let's save this plot::

  pdf("csaw-MDS.pdf")
    par(mfrow=c(2,2))
  adj.counts <- cpm(data.filt.calc, log=TRUE)
  for (top in c(100, 500, 1000, 5000)) {
  plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"),labels=c("hela", "hela", "sknsh", "sknsh"), top=top)
  }
  dev.off()


.. admonition:: csaw-MDS.pdf
   :class: dropdown, warning

   .. image:: figures/csaw-MDS.png
          :width: 600px




Annotation of the results
===========================

.. code-block:: R

  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  anno <- detailRanges(merged$region, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,
  orgdb=org.Hs.eg.db, promoter=c(3000, 1000), dist=5000)

  merged$region$overlap <- anno$overlap
  merged$region$left <- anno$left
  merged$region$right <- anno$right


Creating the final object with results and annotation
======================================================

Now we bring it all together:

.. code-block:: R

  all.results <- data.frame(as.data.frame(merged$region)[,1:3], table.combined, anno)


.. admonition::  all.results, unsorted
   :class: dropdown, warning


   .. code-block:: R

        > head(all.results)
      seqnames  start    end num.tests num.up.logFC num.down.logFC       PValue
    1     chr1  10051  10360         7            0              5 2.328912e-04
    2     chr1  29301  29410         3            0              3 6.989334e-06
    3     chr1 100451 100560         3            0              3 1.948039e-04
    4     chr1 151101 151310         5            0              5 4.108169e-05
    5     chr1 246751 246860         3            0              3 6.674578e-05
    6     chr1 408351 408560         5            0              5 1.880546e-04
               FDR direction rep.test rep.logFC     overlap          left
    1 0.0040397108      down        1 -7.239404 DDX11L1:+:P              
    2 0.0004822892      down        8 -7.000913 WASH7P:-:PE WASH7P:-:4410
    3 0.0036799249      down       13 -7.501513                          
    4 0.0011680754      down       14 -7.121331                          
    5 0.0017204192      down       19 -7.208420                          
    6 0.0036207355      down       23 -8.909876                          
                             right
    1 DDX11L1:+:1514,WASH7P:-:4002
    2                 WASH7P:-:414
    3                             
    4                             
    5                             
    6    




All significant regions are in:

.. code-block:: R
  
  sig=all.results[all.results$FDR<0.05,]


How many significantly different pekas?::

  > nrow(sig)
  [1] 432



To view the top of the ``all.results`` table:

.. code-block:: R
  
  all.results <- all.results[order(all.results$PValue),]

  > head(all.results)
     seqnames     start       end num.tests num.up.logFC num.down.logFC
  1726     chr2  25642751  25642760         1            1              0
  822      chr1 143647051 143647060         1            0              1
  876      chr1 149785201 149785210         1            0              1
  386      chr1  40530701  40530710         1            1              0
  2519     chr2 199778551 199778560         1            1              0
  1613     chr2   8683951   8683960         1            0              1
             PValue          FDR direction rep.test rep.logFC
  1726 7.875683e-07 0.0004407602        up     6126  7.360348
  822  1.197351e-06 0.0004407602      down     3065 -6.938285
  876  1.197351e-06 0.0004407602      down     3263 -6.938285
  386  1.574877e-06 0.0004407602        up     1616  7.389491
  2519 1.574877e-06 0.0004407602        up     8988  7.389491
  1613 2.198223e-06 0.0004407602      down     5724 -6.997646
                     overlap          left           right
  1726              DTNB:-:I    DTNB:-:347                
  822                                      100286793:-:579
  876  H2BC18:-:P,H3C13:-:PE H2BC18:-:1273                
  386               CAP1:+:I    CAP1:+:470     CAP1:+:1177
  2519                                                    
  1613                                                    


We of course discourage ranking the results by p value only ;-).

Now you are ready to save the results as a table, inspect further and generate a compelling scientific hypothesis.
You can also compare the outcome with results obtained from peak-based couting approach.

One final note: In this example we have used preprocessed bam files, i.e. reads mapped to the regions of spurious high signal in ChIP-seq (i.e. the ENCODE "blacklisted regions") were removed, as were the so called **duplicated reads** - reads mapped to the same genomic positions. While filtering out the blacklisted regions is always recommended, **removal of duplicated reads is not recommended** for DB analysis, as they may represent true signal. As always, your mileage may vary, depending on the project (library sequencing depth and complexity being the factors to watch out for in this context), so exploring several options is essential for obtaining meaningful results.



.. admonition:: relevant information from sessionInfo()
   :class: dropdown, warning

   .. code-block:: R

     other attached packages:
     [1] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
     [2] GenomicFeatures_1.42.3                 
     [3] org.Hs.eg.db_3.12.0                    
     [4] AnnotationDbi_1.52.0                   
     [5] csaw_1.24.3                            
     [6] SummarizedExperiment_1.20.0            
     [7] Biobase_2.50.0                         
     [8] MatrixGenerics_1.2.1                   
     [9] matrixStats_0.58.0                     
     [10] GenomicRanges_1.42.0                   
     [11] GenomeInfoDb_1.26.7                    
     [12] IRanges_2.24.1                         
     [13] S4Vectors_0.28.1                       
     [14] BiocGenerics_0.36.0                    
     [15] edgeR_3.32.1                           
     [16] limma_3.46.0         





.. <!-- #### for sanity reasons: (need to dig deeper to find a better example)
.. check against macs2 peaks

.. bedtools intersect -a hela_1_peaks.chr12.bed -b hela_2_peaks.chr12.bed -f 0.50 -r > peaks_hela.chr12.bed
.. bedtools intersect -a sknsh_1_peaks.chr12.bed -b sknsh_2_peaks.chr12.bed -f 0.50 -r > peaks_sknsh.chr12.bed

.. bedtools intersect -a peaks_sknsh.chr12.bed -b peaks_hela.chr12.bed -f 0.50 -r > peaks_sknsh_hela.chr12.bed

..     1088 peaks_hela.chr12.bed
..     2031 peaks_sknsh.chr12.bed
..  473 peaks_sknsh_hela.chr12.bed



.. all.results <- all.results[order(all.results$start),]

.. macs2 in sknsh 1:
.. chr1 1270265 1270622 sknsh_1_REST.enc.macs2_peak_25  2714  . 80.09766  275.94952 271.41241 304

.. csaw DB:
.. 11       chr1   1270251   1270610        8        7
..  -->



.. ----

.. Written by: Agata Smialowska


