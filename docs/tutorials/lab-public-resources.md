---
layout: default
title:  'Public data'
---

# Public ChIP-seq resources

Aim of this exercise

* To see what public datasets are out there

* To be able to find and download public ChIP-seq data in useful formats

This exercise is mostly run in a web browser, so it’s easiest to run it on you local computer.



## 1. ENOCDE

[https://www.encodeproject.org](https://www.encodeproject.org)

Besides ENCODE data, the ENCODE data portal also contains data from the Roadmap Epigenome project, as well as the modENCODE & modERN projects (for fly and worm).

To explore this data repostioty, go to the [encode website](https://www.encodeproject.org) and select *Data* and then *Search*. Say that we want to see all data sets for the histone mark H3K27me3. Start by by typing "H3K27me3" in the search box in the top right corner. **How many results do you see?** This refers to everything in the encode data base: experiments, series of experiments, publications etc.

You can select subsets of the results from the panel on the left. Select *Experiment* to only see experiments. **How many results do you see now? Are all these ChIP-seq experiments?**

Now, let’s make a finer selection: Select only released experiments (from *Experiment Status*), and then only experiments using the GRCh38 genome. **How many experiments do we have now?**

Perhaps we are only interested in data from the brain. Under *Organ*, click on *See more*, and the select *brain*. **How many experiments do we have now?**

You can see a list of all experiments to the right. Click on the first one, and open the page in a new browser window or tab. This will take you to a page describing this experiment, and what protocols and analysis pipelines were used. If you scroll down and select the tab *File details* you will se a list of all files that are available for download.

**Do you know what these files are?** Try downloading some if you want to. But since some of these files are large, remember to remove them when you are done looking at them.

Now have seen how to find and download a single ChIP-seq data set. If we instead want to download all H3K27me3 data from the human brain mapped to the GRCh38 genome, go back to the ENCODE search page, where you had selected the relevant experiments. Then click the button *Download*. This will download a file to your computer.

Open this file. It contains URLs to all data files for the selected experiments. If we want to download e.g. only the bed files with peaks that are stable in both replicates of each experiment, we need to do some extra steps.

First, we download the meta data table. The URL is the first line in the file you just downloaded:

```
wget "https://www.encodeproject.org/metadata/searchTerm=H3K27me3&type=Experiment&status=released&assembly=GRCh38&organ_slims=brain/metadata.tsv"
```

You can open this file, e.g. in excel to have a look. We want to select all lines corresponding to bed files with stable/replicated peaks, using the GRCh38 genome, and save these in a new file `metadata_peak_files.tsv`:

```
grep bed.gz metadata.tsv | grep "stable\|replicated" | grep GRCh38 > metadata_peak_files.tsv
```

From this file we can now get the URLs, which are in column 37:

```
cut -d$'\t' -f37 metadata_peak_files.tsv > metadata_peak_files_urls.txt
```

Finally, we can now download all these files:

```
wget -i metadata_peak_files_urls.txt
```

This will download load all peak files. These will still have non-informative names, e.g. `ENCFF591RMN.bed.gz`. To see which experiment each file corresponds to, look in `metadata_peak_files.tsv`


### (Roadmap epigenomics)
There are several ways to download data from the Roadmap Epigenomics web site. But since these data sets are also available through ENCODE, it's probably easiest to use the ENCODE web site. In case you want to have a look, these pages also host the Roadmap Epigenomics data:

[https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/](https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/)

[http://genboree.org/EdaccData/Release-9/](http://genboree.org/EdaccData/Release-9/)



## 2. Cistrome

[http://cistrome.org/db/#/](http://cistrome.org/db/#/)

Cistrome is another database where ChIP-seq data has been collected and processed uniformly. This is a good complement to the big projects (ENCODE & Roadmap Epigenomics) since it collects data from many smaller studies.

As an example, we will look for data on the three human transcription factors *Grhl1*, *Grhl2* and *Grhl3*. Grhl stands for “Grainy head like”, which means that these proteins are similar to the *Grainy head* protein first found in fruit fly. In human there are 3 Grhl homologs: *Grhl1*, *Grhl2* and *Grhl3*. They are involved in development and would healing, and have been implicated in hearing loss and cancer.

To see which ChIP-seq data sets are available for the Grhl proteins, go to [http://cistrome.org/db/#/](http://cistrome.org/db/#/), type "Grhl" in the search box and click on *Search*. We can then refine the search further by selecting *Homo sapiens* under *Species*. **Which Grhl proteins do we find data for?**

There is also an option to filter data on quality measures. To try this, click on *Options* and then *Samples passing peak quality controls*. **Which Grhl proteins do we still have data for after this filtering?**

Now, select the first data set in the list. It will be highlighted in blue. Scroll down to see details about this sample. Select the tab *QC reports*. **Can you make sense of this information? Does this look like an experiment that worked?**

In the cistrome database, motif finding programs were run on all transcription factor data sets. Select the tab *QC motifs*, and have a look at some of the top motifs. **Do they look similar the known Grhl site? Can you find what the Grhl site is?** Hint: Go to JASPAR and search for "Grhl".

It’s also possible to download files from cistrome. For each experiment, bigwig files with the read coverage signal and bed files with peaks are available. Note that batch download is not available for manually selected batches, just for e.g. all human transcription factor data or all mouse chromatin data etc.
