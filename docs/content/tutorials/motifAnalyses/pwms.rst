====================================================================
Motifs: representing them and searching for hits against a reference
====================================================================

:Author: Dania Machlab
:Authors:
   Dania Machlab
:Date: 2025-07-15

.. contents::
   :depth: 3
..

Overview
========

In this lab, we will examine how transcription factor binding sites
(TFBSs) or motifs can be represented in the form of Position Weight
Matrices (PWMs). We will see how these matrices can be created,
downloaded from public databases like JASPAR, and use them to scan
against a reference sequence in search for motif hits.

Learning outcomes
=================

-  understand how motifs can be represented in matrices.
-  understand the differences between the PFM, PPM, PWM, and ICM as well
   as how these matrices are calculated.
-  know how to load PWMs in ``R`` and scan for motif hits against a
   reference.

Libraries
=========

We start by loading the needed packages. If necessary, use
``BiocManager::install()`` to install missing packages.

.. container:: cell

   .. code:: r

      suppressPackageStartupMessages({
        #library(TFBSTools)
        library(Biostrings)
        #library(monaLisa) # use findMotifHits which uses matchPWM (optimised for parallel processing of multiple PWMs)
        library(JASPAR2024)
        library(RSQLite)
      })

Background
==========

TFs and their activity

explain PFMs, PPMs, PWMS, ICMs

Here we will explore PFMs, PPMs, PWMs and ICMs. Note that representing
motifs in these matrix forms assumes independence of positions.

We will explore these matrices and how they are created using the CTCFL
motif from `JASPAR2024 <https://jaspar.elixir.no/matrix/MA1102.3/>`__ as
an example.

Motif databases
---------------

jaspar and hocomoco

PFMs or PCMs?
=============

The position frequency matrix (PFM), also known as the position count
matrix (PCM), depicts the frequencies of the nucleotides i.e. the number
of times each nucleotide occurs at each position of the motif.

Using CTCFL as an example, we will read in the sequences here and
produce the corresponding PFM. We will use the ``Biostrings`` package to
read the sequences and represnt them as a ``DNAStringSet`` object, which
is a useful way to represent and manipulate DNA sequences in ``R``.

.. container:: cell

   .. code:: r

      # read in the TFBS sequences
      CTCFLsequencesFile <- "data/MA1102.3.sites"
      CTCFLsequences <- readDNAStringSet(CTCFLsequencesFile)
      CTCFLsequences

   .. container:: cell-output cell-output-stdout

      ::

         DNAStringSet object of length 18037:
                 width seq                                           names               
             [1]     8 CAGGGGGC                                      hg38_chr1:869925-...
             [2]     8 CAGGGGGC                                      hg38_chr1:904775-...
             [3]     8 GAGGGGGC                                      hg38_chr1:925040-...
             [4]     8 CAGGGGGC                                      hg38_chr1:945418-...
             [5]     8 GAGGGGGC                                      hg38_chr1:951563-...
             ...   ... ...
         [18033]     8 AAGGGGGC                                      hg38_chrX:1550575...
         [18034]     8 GAGGGGGC                                      hg38_chrX:1552166...
         [18035]     8 CAGGGGGA                                      hg38_chrX:1552168...
         [18036]     8 CAGGGGGC                                      hg38_chrX:1554355...
         [18037]     8 GAGGGGGC                                      hg38_chrX:1556122...

   .. code:: r

      # create a PFM by counting nucleotide occurrences per position
      pfm <- consensusMatrix(CTCFLsequences)
      pfm <- pfm[c("A", "C", "G", "T"), ]
      pfm

   .. container:: cell-output cell-output-stdout

      ::

            [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]
         A  1301 17270    86   684   349   193   351   181
         C 12867   413   481   600  1141   767   270 17593
         G  2033   329 17285 16511 15425 16819 17001    68
         T  1836    25   185   242  1122   258   415   195

PPMs
====

We can now calculate the probability of observing each nucleotide at a
particular position by dividing count of each nucleotide by the total
count per position.

.. math::


   PPM_{ij} = \frac{count_{ij}}{\sum_{i}{count_{ij}}}

\ Where :math:`i` is the nucleotide and :math:`i \in \{A, C, G, T\}`,
and :math:`j` is the positoin along the motif.

.. container:: cell

   .. code:: r

      # calculate ppm
      ppm <- sweep(x = pfm, MARGIN = 2, STATS = colSums(pfm), FUN = "/")
      ppm

   .. container:: cell-output cell-output-stdout

      ::

                 [,1]       [,2]        [,3]       [,4]       [,5]       [,6]       [,7]
         A 0.07212951 0.95747630 0.004767977 0.03792205 0.01934912 0.01070023 0.01946000
         C 0.71336697 0.02289738 0.026667406 0.03326496 0.06325886 0.04252370 0.01496923
         G 0.11271276 0.01824028 0.958307923 0.91539613 0.85518656 0.93247214 0.94256251
         T 0.10179076 0.00138604 0.010256695 0.01341687 0.06220547 0.01430393 0.02300826
                  [,8]
         A 0.010034928
         C 0.975383933
         G 0.003770028
         T 0.010811110

   .. code:: r

      # all positions now sum to 1
      colSums(ppm)

   .. container:: cell-output cell-output-stdout

      ::

         [1] 1 1 1 1 1 1 1 1

We can now calculate the probability of observing a certain motif
sequence by multiplying the probabilities of each nucleotide per
position.

.. container:: cell

   .. code:: r

      p_AATTGGTT <- ppm["A", 1] * ppm["A", 2] * ppm["T", 3] * ppm["T", 4] * 
        ppm["G", 5] * ppm["G", 6] * ppm["T", 7] * ppm["T", 8] 
      p_AATTGGTT

   .. container:: cell-output cell-output-stdout

      ::

                    A 
         1.885169e-09 

In this example, we do not have any zero counts for a particular
nculeotide and position. What would happen if we did? Let us suppose
that ``pfm["T", 7]`` had a count of zero and thus ``ppm["T", 7]`` is
also zero. Multiplying by zero would result in ``p_AATTGGTT = 0``. This
can particularly be a problem when starting from a low number of
sequences to calculate the PFM. To avoid low count issues, we will add a
peudo-count when calculating the PPM. We will re-calculate the PPM but
adding a pseudo-count :math:`p` of 1 per position. This corresponds to a
pseudo-count of :math:`p/N` where :math:`N` is the number of
nucleotides, and :math:`N=4` in our case.

.. math::


   PPM_{ij} = \frac{count_{ij}+\frac{p}{N}}{\sum_{i}{count_{ij}+p}}

.. container:: cell

   .. code:: r

      pseudooCount <- 1
      N <- nrow(pfm)

      pfmWithPseudo <- pfm + pseudooCount/N
      ppm <- sweep(x = pfmWithPseudo, MARGIN = 2, STATS = colSums(pfm), FUN = "/")
      ppm

   .. container:: cell-output cell-output-stdout

      ::

                 [,1]       [,2]        [,3]       [,4]       [,5]       [,6]       [,7]
         A 0.07214337 0.95749016 0.004781837 0.03793591 0.01936298 0.01071409 0.01947386
         C 0.71338083 0.02291124 0.026681266 0.03327882 0.06327272 0.04253756 0.01498309
         G 0.11272662 0.01825414 0.958321783 0.91540999 0.85520042 0.93248600 0.94257637
         T 0.10180462 0.00139990 0.010270555 0.01343073 0.06221933 0.01431779 0.02302212
                  [,8]
         A 0.010048789
         C 0.975397793
         G 0.003783889
         T 0.010824971

PWMs
====

The position weight matrices (PWM) is also known as the
position-specific scoring matrix or the logodds scoring matrix. Here,
log-odds scores are calculated by comparing the probabilities we have in
the PPM to the probabilities of observing each nucleotide outside of a
binding site (background nucleotide probabilities). Assuming a uniform
background, in which each nucleotide has an equal probability would give
us the following background probabilities for each nucleotide:
:math:`p(A) = p(C) = p(G) = p(T) = 0.25`. The log-odds scores can be
obtained as follows:

.. math::


   PWM_{ij}=log2\Bigl(\frac{PPM_{ij}}{B_i}\Bigr) 

Where :math:`i` is the specific nucleotide, :math:`j` is the position
along the motif, and :math:`B_i` is the background probability for
nucleotide :math:`i`. Thanks to the pseudo-count we have added, we will
avoid situations where we are taking the :math:`log2(0)` which is
:math:`-Inf`.

.. container:: cell

   .. code:: r

      (B <- c("A" = 0.25, "C" = 0.25, "G" = 0.25, "T" = 0.25))

   .. container:: cell-output cell-output-stdout

      ::

            A    C    G    T 
         0.25 0.25 0.25 0.25 

   .. code:: r

      pwm <- log2(sweep(x = ppm, MARGIN = 2, STATS = B, FUN = "/"))
      pwm

   .. container:: cell-output cell-output-stdout

      ::

                [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
         A -1.792989  1.937330 -5.708219 -2.720292 -3.690555 -4.544347 -3.682317
         C  1.512744 -3.447801 -3.228029 -2.909252 -1.982273 -2.555119 -4.060521
         G -1.149100 -3.775632  1.938582  1.872490  1.774334  1.899154  1.914681
         T -1.296125 -7.480460 -4.605342 -4.218319 -2.006493 -4.126047 -3.440835
                [,8]
         A -4.636835
         C  1.964063
         G -6.045915
         T -4.529493

The score for a particular sequence can be calculated by combining the
PWM scores at each position. For example the score for ``AATTGGTT`` is:

.. container:: cell

   .. code:: r

      pwm["A", 1] + pwm["A", 2] + pwm["T", 3] + pwm["T", 4] + pwm["G", 5] + pwm["G", 6]+ pwm["T", 7] + pwm["T", 8]

   .. container:: cell-output cell-output-stdout

      ::

                 A 
         -12.97616 

EXERCISE: calculate the score for sequence XXX

ICMs
====

JASPAR1014 and useful functions
===============================

We could have directly gotten the PFM matrix from the ``JASPAR2024``
package as well. Let us get it here and compare it to the one we have
generated.

.. container:: cell

   .. code:: r

      JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

      # extracting the PFM
      TFBSTools::getMatrixSet(JASPARConnect, opts = list(tax_group = "vertebrates", 
                                                         collection="CORE", 
                                                         matrixtype = "PFM", 
                                                         name = "CTCFL"))




      # converting to a PWM (comment on the default pseudo-count)


      # directly extracting the PWM


      # 
      # getMatrixByName(JASPARConnect, "CTCFL")

Scanninf for motif hits
=======================

exercise: scan hits for a list of PWMs (get them from jaspar)

Scanning for motif hits
=======================

all throughout, if they come up, give a little info on some useful bioC
objects like DNAstringset and link to Bioconductor and the package

refer to section X on PWMs and how the scores are summed for a
particular sequence, and that this is used when scanning for motif hits.

Session
=======

.. container:: cell

   .. code:: r

      date()

   .. container:: cell-output cell-output-stdout

      ::

         [1] "Tue Jul 15 11:08:28 2025"

   .. code:: r

      sessionInfo()

   .. container:: cell-output cell-output-stdout

      ::

         R version 4.5.1 (2025-06-13)
         Platform: aarch64-apple-darwin20
         Running under: macOS Sequoia 15.5

         Matrix products: default
         BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
         LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

         locale:
         [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

         time zone: Europe/Stockholm
         tzcode source: internal

         attached base packages:
         [1] stats4    stats     graphics  grDevices utils     datasets  methods  
         [8] base     

         other attached packages:
          [1] RSQLite_2.4.1        JASPAR2024_0.99.7    BiocFileCache_2.16.0
          [4] dbplyr_2.5.0         Biostrings_2.76.0    GenomeInfoDb_1.44.0 
          [7] XVector_0.48.0       IRanges_2.42.0       S4Vectors_0.46.0    
         [10] BiocGenerics_0.54.0  generics_0.1.4      

         loaded via a namespace (and not attached):
          [1] bit_4.6.0               jsonlite_2.0.0          dplyr_1.1.4            
          [4] compiler_4.5.1          crayon_1.5.3            filelock_1.0.3         
          [7] tidyselect_1.2.1        blob_1.2.4              yaml_2.3.10            
         [10] fastmap_1.2.0           R6_2.6.1                curl_6.4.0             
         [13] knitr_1.50              tibble_3.3.0            GenomeInfoDbData_1.2.14
         [16] DBI_1.2.3               pillar_1.11.0           rlang_1.1.6            
         [19] cachem_1.1.0            xfun_0.52               bit64_4.6.0-1          
         [22] memoise_2.0.1           cli_3.6.5               magrittr_2.0.3         
         [25] digest_0.6.37           lifecycle_1.0.4         vctrs_0.6.5            
         [28] evaluate_1.0.4          glue_1.8.0              rmarkdown_2.29         
         [31] httr_1.4.7              tools_4.5.1             pkgconfig_2.0.3        
         [34] htmltools_0.5.8.1       UCSC.utils_1.4.0       

References
==========

https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf
stormo paper (in intro on PWMs) biostrings package (when mentioning
DNAStringset) TFBStools package (when mentioning the functions and
package) JASPAr original
publication(https://watermark.silverchair.com/gkh012.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA3QwggNwBgkqhkiG9w0BBwagggNhMIIDXQIBADCCA1YGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMd1c0A9OaUNYKpGaVAgEQgIIDJwwUmaHjFou4o64xzqwnzE_Os5--EEm_ueATvyN6JBYtLUFgdM2LrmZwHC_RB2zWy16XOuZOhGFlJ8MW0ccOEek118cSYHMiMjkT1257Ge5Hec3CLYhoBDvjiC84wx8xsr8LsI7mHuizgl_yMZDj16gKnJASjb7XLZHEdhpe0uwvcA41t4B_UWRAEPUJGBrbxiHHs93QAD2NL_bsCG7NWkdaFpgTktsoXN2XXwN7SwTNDX4IkvqEhqFByY9OpIndE3Nwd2LgWzRq1lnmqfajfNoO7mT4RJAmN_tUGDxTR5jNjWro9yt0JvdrpEIivXROKmImpjy0fyKCCYMyZYDHRo6qihSKGLkJa7kQipXrg4tNvIJ3yWz1tLc3c3t5BrcrSawKDuxnEYc1oZWGwULFqztkm2g4M0qY-GT_NQSpjeu83L4iWTSP_NAiLVAlwGqb3oTXDRbE0QBE69WbG5AqcHk7GK8Oec3AjadMXszBJp32M42F2qpxBDRK-a7O4rnIgXHrSN1wJZyNbQjgFh9KS4dCXGvIZWzFyRnG7SA_gEJDGKDXxMAS8mqJBPwRtLNYMauMIB4Z8Iqabr6q64vyysktCngYFovyorKNfoSBnlf8GBLpgPHdKWkIFIlallWAFvt5qcaIG_-3rgSqUBfM9BVVbid7zxIWL0hEQRO-bPVRKBJdoittpFRLXbywlXb-CTU51nQDKuo8va4rHQ5c0lPAvR_pHI9JR4mFS7icZEESnqkVJ4YGG0kXE3JpwFrYgOhCoJydAU8R8AZEtFjVZ_YpUOlgn-sBgUXxtl3Afz-FGIQsZOQBxskbmWSueE4wdVvCxCcuU_yldq043V8kaLWBFX8AXWaBaFnCATln1T1WVJJANXMof5OwjMgLaFlzOMip-jNxQybq56InT1ttMQADzmtauaycr6fJ47_wacJRV3BiOursLhZLOTBi3X39JnzypFtncJmAZGCAw-ohWepL0MO59shETxg6F6dmcCVNOc7SQ96fmyke1saOTuauhVgCOAehnwnzzO8GIUp0b8AYuckz2N26Dabb9UqkuHyFpicqe_K7dQ)

when mentioning the database in the intro
