.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

Setup
-----

In these tutorials we will be using the following software:

- **deepTools**. Analysis suite for sequencing data. You can read about the installation `in their documentation <https://deeptools.readthedocs.io/en/develop/content/installation.html>`_.
- **IGV**. `Integrative Genomics Viewer <http://software.broadinstitute.org/software/igv/>`_.
- **seqplots**. R package for the analysis of bigWig files. It can be run as a shiny app and can be installed from Bioconductor. Check `seqplots Bioconductor web <https://bioconductor.org/packages/release/bioc/html/seqplots.html>`_ for more details.
- **intervene**. Easy to use R tool that allows to plot overlap between bed files. It can be installed through bioconda and pip (for me it worked using pip, but not conda though). Bedtools is a dependency. See more detailed info `in the intervene documentation site <https://intervene.readthedocs.io/en/latest/install.html>`_.

Additionally, other Bioconductor R packages will be needed, mostly as dependencies: :code:`GenomicRanges`, :code:`rtracklayer`.

Bioconductor packages can be installed in the usual way:

.. code-block:: r

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("seqplots")
  BiocManager::install("GenomicRanges")
  BiocManager::install("rtracklayer")


You can run everything from Uppmax as well as on your local laptop. If you choose to run examples on Uppmax, the relevant module load calls are:

.. code-block:: bash

  module load deepTools
  module load IGV
  module load R_packages

:code:`R_packages` contains all the mentioned R packages and their dependencies: :code:`GenomicRanges`, :code:`rtracklayer`, :code:`seqplots`. Be sure to have the corresponding modules loaded (and :code:`bioinfo-tools` before anything) before running them.

.. hint:: 
    Some of the results here can be computationally demanding. In those cases Uppmax will be recommended. In
    any case, the resulting files of performing these operations will also be available in the course directory, to avoid technical issues delaying the whole tutorial.