### Monday
09.15 - 09.30 Introduction to the course (*online session*) (Olga, Agata)

09.30 - 09.45 Introduction to methylation (*online session*) (Vincent)

09.45 - 10.00 Break-out session (*online session*) (Vincent)

10.00 - 11.30 Methylation methods & technologies (*online session*) (Jessica Nordlund)

11.30 - 12.00 Computer exercises set-up (*online support*)

12.00 - 13.00 lunch (*offline*)

13.00 - 13.30 Methylation Exercises Overview I: Array workflow (Vincent) [internal notes: load data, QC, Normalisation, probe and region DM, GO enrichment]

13.30 - 15.30 Exercises (*online support*)

15.30 - 16.00 Methylation Exercises Overview II: Methylation Sequencing (Vincent)

16.00 - 17.00 Exercises (*online support*)

17.00 - 17.30 Test yourself (*offline*)

### Tuesday
09.15 - 10.00 Recap of the previous day and activation exercise (*online session*) (Vincent)

10.00 - 11.30 Introduction to ChIP-seq & ATAC-seq (2 x 30 min, 1 break-out session) (Agata) (*online session*) [internal notes: intro, experimental design, data processing workflow incl. correlation, ATAC-seq]

11.30 - 13.00 lunch (*offline*)

13.00 - 13.30 Calling broad peaks, differential expression, motifs (*online session*) (Agata, Olga, Jakub)

13.30 - 15.30 Exercises (*online support*)

15.30 - 17.00 Exercises (*offline*)

17.00 - 17.30 Test yourself (*offline*)

### Wednesday (Simon's day: needs revisit)
09.15 - 10.30 Recap of the previous day and activation exercise (group session) (*online session*)

10.30 - 11.30 Quantitative ChIP-seq methods (Simon)

11.30 - 13.00 lunch (*offline*)

13.00 - 13.30 TBD

13.00 - 15.30 Exercises (*online support*)

15.30 - 17.00 Exercises (*offline*)

17.00 - 17.30 Test yourself (*offline*)

### Thursday

Phil / Harshil

09.15 - 10.00 Recap of the previous day (*online session*)

10:00 - 10:40 - Intro to Nextflow (Phil)
* Intro to me and NGI
* Get Nextflow, get Singularity
* Config files
* -resume, other core options
* `-bg` to run in the background (`screen`, `tmux`)

10:40 - 11.30 - Exercises
* Run atacseq with `-profile test,uppmax`
* Rerun with `-resume`
* Use nf-core launch to run

11.30 - 13.00 lunch (*offline*)

13:00 - 13:45 - Intro to the ChIP-seq and ATAC-seq pipelines (Harshil)
* What they do
* Where to find documentation
* Run a test dataset (-profile test)
* nf-core launch (needs testing - installed as a module?)

13.45 - 15.30 Exercises (*online support*)

* Run chipseq with chr1 subset data from previous days
  * Write a sample sheet, use `-r` for pipeline version
* Rerun with `-resume` and BroadPeaks
* Email results with MultiQC results
* Use nf-core launch to run with above data / sample sheet
* Run either pipeline with nf-core AWS megatest dataset
* Run a pipeline on their own data
* Visual inspection of results, basic quality control and any other results based interrogation
* Ask users to compare results to see if they are the same! Reproducible results. Woooohoooo!

> **TODO:** Try running the subsampled test datasets from FastQ and time, using job submission with the project ID
> * Agata: Generate FastQ files from subsampled BAMs
> * Someone: Try to run exercises

> **TODO:** Download AWS-megatest FastQ files to UPPMAX
> * Harshil: Give paths to agata from samplesheet: https://raw.githubusercontent.com/nf-core/test-datasets/chipseq/design_full.csv
> * Agata: Download and save to UPPMAX

> **TODO:** Harshil / Phil - write up to 10 questions
> * At least one theory based
> * At least one practical (eg. listing pipelines, find latest release)

> **TODO:** Harshil / Phil - write up the exercises in markdown on the website
> * https://github.com/NBISweden/workshop-epigenomics-RTDs/

15.30 - 17.00 Exercises (*offline*)

17.00 - 17.30 Test yourself (*offline*)

### Friday
9.15 - 09.30 Recap of the previous day (*online session*)

09.30 - 10.20 scATAC-seq / scChIP-seq (Jakub) (*online session, including 10 minute break*)

10.20 - 10.30 Break

10.30 - 11.30 Data integration (Jakub) (*online session, including 10 minute break*)

11.30 - 13.00 lunch (*offline*)

13.00 - 13.20 Intro to exercises (*online session*)

13.20 - 15.00 Exercises (*online support*)

15.00 - 15.30 Test yourself (*offline*)

15.30 - 16.30 Course wrap-up (*online session*)
