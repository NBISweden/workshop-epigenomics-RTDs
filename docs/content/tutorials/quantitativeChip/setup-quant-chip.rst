.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

Setup
-----

In these tutorials we will be using the following software:

- **deepTools**. Analysis suite for sequencing data. You can read about the installation `in their documentation <https://deeptools.readthedocs.io/en/develop/content/installation.html>`_.
- **IGV**. `Integrative Genomics Viewer <http://software.broadinstitute.org/software/igv/>`_.
- **intervene**. Easy to use R tool that allows to plot overlap between bed files. It can be installed through bioconda and pip (for me it worked using pip, but not conda though). Bedtools is a dependency. See more detailed info `in the intervene documentation site <https://intervene.readthedocs.io/en/latest/install.html>`_.
- **R**. :code:`ggplot2` and a couple Bioconductor R packages will be needed, mostly as dependencies: :code:`GenomicRanges`, :code:`rtracklayer`. Furthermore, the Spike-in tutorial requires some other packages: :code:`ChIPSeqSpike` and :code:`BSgenome.Hsapiens.UCSC.hg38`.

You can run R in any preferred way you have been running during this workshop. The dependencies are small so you probably have these packages installed. Most definitely these are available through Uppmax module system.

Bioconductor packages can be installed in the usual way:

.. code-block:: r

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GenomicRanges")
  BiocManager::install("rtracklayer")


You can run things from Uppmax as well as on your local laptop. The relevant module load calls on Uppmax are:

.. code-block:: bash

  module load deepTools
  module load IGV
  module load R_packages

:code:`R_packages` contains all the mentioned R packages and their dependencies: :code:`GenomicRanges`, :code:`rtracklayer`. Be sure to have the corresponding modules loaded (and :code:`module load bioinfo-tools` before anything) before running them.

Across the tutorials, the environment considered best for running each step will be specified as **Uppmax** or **local** (your computer).

In general, Uppmax is preferrable for computationally demanding steps, however working on a local environment tends to be more straightforward.

Things that can easily be done locally: IGV visualization (you can download the bigWig files, peaks or genomic annotations, since these are not large files), R figures.

Things that will run faster on Uppmax: The first part of the minute tutorial is more demanding. It is still small enough that can be run on a laptop, but it will take a few hours. Some deepTools calls, especially fingerprint plots, are also slow, so it may be better to run on Uppmax.

.. note:: 
    Computationally demanding steps have been precalculated and the resulting plots and data are available throughout the tutorials, so if there is some issue
    in any of these steps, you can pick up from the precalculated data. In those cases it will be noted within the corresponding section. 
