.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html
    
======================================================
Identifying Relevant Transcription Factors
======================================================

:raw-html:`<br />`

DNA-binding proteins known as transcription factors (TFs) add another facet from 
which gene expression can be controlled. Gene expression can be regulated by the
interaction between enhacer elements and promoters, and TFs play key roles in 
these interactions. TF binding sites (TFBSs) tend to be in the range of 6-12 
bases (`Spitz & Furlong <https://www.nature.com/articles/nrg3207>`_). 

The state of the epigenome also influences TF binding. The binding of certain
TFs can be methylations sensistive, as showcased by 
`Domcke, Bardet et al. <https://www.nature.com/articles/nature16462>`_

So far, we have learnt how to process various assays like DNA methylation, 
ATAC-seq, and ChIP-seq. Now we want to pose the question: 
Which TFs could be driving these changes in methylation,
chromatin accessibility or histone modification that we observe? 

This tutorial covers two main parts that tackle these kinds of analyses. First,
we examine how TFBSs, also known as motifs, can be represnted. We will
then use this along with other tools to look for TFs that are likely to explain
the experimental read outs we see. In that part, we will focus on ATAC-seq
data, and look for motifs that could potentially explain the differences in 
accessibility we see between two conditions.



delete below: 


TFs come with their own flavors in terms of 
functionality (REF here nature rev). Some may for example act as remodellers, 
opening up previously inaccessible parts of the genome which may include new
binding sites for other TFs. Others may bind enhancers or promoters and be
incolved in the gene regulation 

On top of this, TF binding to open chromatin can also be governed by epigentic
modifications like DNA methylation (ref)







Genome regulation and control can be observed from sevral aspects. Chemical 
modifications to the nucleotides like DNA methylation or the various flavors of
histone modifications have an influence on gene expression as well as 
affecting binding capabilites of certain 

These can control cell fate or diference between normal and disease state


The genome can be regulated at several different levels. So far, we have seen 
examples of 

In this course our 
focus has been on level beyond the raw DNA sequence, the epigenetics level. 
While this can include chemical modifications to the DNA like methylation 
and histone modifications in various falvors, gene expression 

describe some intro/background (general here) on this aspect of downstream analyses



Gene regulation and control can be achieved by the interaction of distal regions
in the genome, called enhancers, with the gene promoters. These enhancers 
typically have binding sites for DNA-binding proteins called transcription
factors (TFs). The binding of TFs to enhancers in turn can result in the 
recruitment of additional TFs and modify nearby 
chromatin [Shlyueva et al](https://doi.org/10.1038/nrg3682).

The mechanisms by which gene expression is regulated by enhancers and TF binding
is still an active field of research. furlong review and suggestions of some mechanisms? 

One way in which the activity of genes is known to be regulated is through 
interactions with distal regulatory regions called enhancers, facilitated by the 
binding of regulatory DNA-binding proteins called transcription factors. 
The binding of TFs to enhancers in turn can result in the recruitment of
additional TFs and modify nearby 
chromatin [Shlyueva et al](https://doi.org/10.1038/nrg3682). In looping to 
promoters, distal enhancers can activate target genes.



tutorials:

.. toctree::
   :maxdepth: 1

   Representing Motifs <representingMotifs.rst>
   PWMs2 <CopyOfpwms/pwms.rst>


