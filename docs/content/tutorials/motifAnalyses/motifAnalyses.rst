.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html
    
======================================================
Predicting Relevant Transcription Factors
======================================================

:raw-html:`<br />`

DNA-binding proteins known as transcription factors (TFs) can regulate gene
expression by binding to gene promoters and/or enhancers. The binding site is 
typically 6-12 base pairs long and the mechanisms by which TFs
regulate and influence gene expression are quite varied and still not fully
understood (`Spitz & Furlong <https://www.nature.com/articles/nrg3207>`_).

The effect of TF binding can be indirectly observed via associated changes in 
transcription, chromatin accessibility, DNA methylation and histone 
modifications. Understanding which TFs are likely to explain these observed 
experimental differences could help shed light on finding TFs that could be
key players in our system of interest.

This tutorial covers two main parts. First, we examine how TF binding sites 
(TFBSs), also known as motifs, can be represnted. We will then use this along 
with other tools to look for TFs that are likely to explain the experimental 
read outs we see. In that part, we will focus on ATAC-seq data, and look for 
motifs that could potentially explain the differences in accessibility we see 
between two conditions.



.. toctree::
   :maxdepth: 1

   Part 1: Representing Motifs <representingMotifs.rst>
   Part 2: Finding relevant motifs <motifAnalysesWithMonalisa.rst>


