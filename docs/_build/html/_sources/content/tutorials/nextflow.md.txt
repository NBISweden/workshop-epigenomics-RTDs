# Introduction to Nextflow and nf-core

## Learning outcomes

* To understand the basic components in a Nextflow run
* To understand the difference between Nextflow and nf-core
* To be able to find relevant nf-core pipelines
* To be able to download and update nf-core pipelines
* To be able to run a nf-core pipeline

## Setup: Nextflow

These steps typically only need to be done once, to get you set up with Nextflow.

### Installation

To use Nextflow, either install it yourself in your home directory or use the UPPMAX environment module.

To install yourself:

```bash
cd ~/bin    # Your home directory bin folder - full of binary executable files, already on your PATH
curl -s https://get.nextflow.io | bash
```

Alternatively, to use the environment modules:

```bash
module purge                        # Clears all existing loaded modules, to start fresh
module load uppmax bioinfo-tools    # Base UPPMAX environment modules, needed for everything else
module load Nextflow                # Note: Capital N!
```

### Environmet variable setup

Nextflow has a large list of bash _environment variables_ that can be set to configure how it runs.
There are some that you want to always be set the same, and the most convenient way to set these
is to add them to your `.bashrc` file in your home directory. Once here, they will be applied
every time you log in automatically.

```bash
# Nextflow singularity image cachedir
export NXF_SINGULARITY_CACHEDIR=/proj/MYPROJ/nobackup/YOURNAME/singularity-images

# Don't let Java get carried away and use huge amounts of memory
export NXF_OPTS='-Xms1g -Xmx4g'

# Don't fill up your home directory with cache files
export NXF_HOME=/proj/MYPROJ/nobackup/YOURNAME/nxf-home
export NXF_TEMP=${SNIC_TMP:-$HOME/glob/nxftmp}
```

> NB: Replace `MYPROJ` and `YOURNAME` so that you have a sensible file path.

#### Singularity cache dir and this course

To avoid everyone sitting and downloading the same Singularity images for the pipelines,
we have done this in advance. You can find them at `/sw/courses/epigenomics/nextflow/singularity-images/`.

Please copy these two files to the directory that you have set for `NXF_SINGULARITY_CACHEDIR`.

If you try to run a different nf-core pipeline or a different version, Nextflow will try to download
the missing Singularity image to this directory for you. This can often fail. If so, try running
an interactive UPPMAX job and pulling the image to the cache directory manually:

```bash
cd $NXF_SINGULARITY_CACHEDIR
interactive -A MYPROJ    # Your UPPMAX project ID. Starts an interactive job
singularity pull --name nfcore-PIPELINE-VERSION.img docker://nfcore/PIPELINE:VERSION
```

Two things to note:

* `docker://nfcore` does _not_ have a hyphen (`-`)
* The image filename must exactly match this pattern, otherwise Nextflow will ignore it and try to download again

### Check that Nextflow works

It's always good to have a mini test to check that everything works.

Make a new directory and run the Nextflow test command as follows:

```bash
mkdir nextflow-hello-test
cd nextflow-hello-test
nextflow run hello
```

You should see something like this:

```
N E X T F L O W  ~  version 20.10.0
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [sharp_sammet] - revision: 96eb04d6a4 [master]
executor >  local (4)
[7d/f88508] process > sayHello (4) [100%] 4 of 4 ✔
Bonjour world!

Ciao world!

Hello world!

Hola world!
```

## Setup: nf-core

The nf-core package is totally separate to Nextflow and is not required to run Nextflow pipelines.
It does however offer some conveinece functions to make your life a little easier.

The nf-core helper package (often referred to as `nf-core/tools`) is written in Python.

### Setting up Python on UPPMAX

It's highly recommended that you set up a custom Python environment on UPPMAX so that you
can install whatever Python packages you want (such as nf-core).

If in doubt, use [_miniconda_](https://docs.conda.io/en/latest/miniconda.html).
First, download and install the correct version for UPPMAX:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```

Follow the instructions and say yes to everything.

Nearly all common bioinformatics tools are available through [Bioconda](https://bioconda.github.io/),
including nf-core/tools. To use bioconda, you must first set up the conda channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Conda uses different _environments_, so you now want to create a new environment that you can
use to install packages. To make a new environment with Python v3.8, called `work`, do the following:

```bash
conda create --name work python=3.8 --yes
```

You may be tidy and use different environments for different tasks, or you can be lazy and just
use one environment for everything - it's up to you! If you fall into the latter category,
you may want to add the following line to your `~/.bashrc` file so that you always activate this
conda environment every time you log in:

```bash
conda activate work
```

### Install nf-core

If you're using conda, the easiest way to install nf-core is via bioconda:

```bash
conda install nf-core
```

Alternatively, you can install from the Python Package Index with `pip`:

```bash
pip install nf-core
```

You can check to see if the installation has worked by trying the following command:

```bash
nf-core --version
```

You should see something like this:

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 1.12

nf-core, version 1.12
```

## Fetching pipelines

Nextflow is pretty clever about working with GitHub and can manage a local cache of pipeline repositories.
Unless you plan on writing your own pipelines, this is nearly always sufficient.

### Listing nf-core pipelines

nf-core can help you find pipelines and pipeline updates, and tells you whether your Nextflow cache is up to date.

First, list all nf-core pipelines:

```bash
nf-core list
```

There are quite a lot, so now filter just those for `chipseq` (try also `rna`, `proteomics` etc)

```bash
nf-core list chipseq
```

You should see something like this:

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 1.12

┏━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name ┃ Stars ┃ Latest Release ┃     Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ chipseq       │    62 │          1.2.1 │ 4 months ago │           - │ -                    │
└───────────────┴───────┴────────────────┴──────────────┴─────────────┴──────────────────────┘
```

### Downloading a Nextflow pipeline

Just running a pipeline with `nextflow run <name>` will fetch it if you don't already have it in the cache,
however you can also explicitly pull a pipeline if you wish. Try by pulling the chipseq pipeline:

```
nextflow pull nf-core/chipseq
```

Now if you re-run the `nf-core list chipseq` command you should see that you have the latest release locally:

```
┏━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name ┃ Stars ┃ Latest Release ┃     Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ chipseq       │    62 │          1.2.1 │ 4 months ago │   yesterday │ Yes (v1.2.1)         │
└───────────────┴───────┴────────────────┴──────────────┴─────────────┴──────────────────────┘
```

### Downloading a specific version of a Nextflow pipeline

It's always a good idea to be specific about the release (or branch) of the Nextflow pipeline that you are running.
That way, if you reuse the command you can be sure that you will be using the same version of the pipeline.

Try downloading an older release of the chipseq pipeline:

```bash
nextflow pull nf-core/chipseq -r 1.2.0
```

You should now see that `nf-core list chipseq` shows that your local cache is not the latest release:

```
┏━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name ┃ Stars ┃ Latest Release ┃     Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ chipseq       │    62 │          1.2.1 │ 4 months ago │   yesterday │ No (v1.2.0)          │
└───────────────┴───────┴────────────────┴──────────────┴─────────────┴──────────────────────┘
```

**Important!** If you don't specify `-r` when using the `nextflow run` command, Nextflow will
use whatever version of the pipeline you have in your local cache. This is another good reason
to always specify it when running a pipeline.

## Running a test workflow

It's always a good idea to start working with a tiny test workflow when using a new Nextflow pipeline.
This confirms that everything is set up and working properly, before you start moving around massive data files.

All nf-core pipelines come with a configuration profile called `test` which will run a minimal test dataset
through the pipeline without needing any other pipeline parameters.

### Trying out atacseq

To try out the nf-core/atacseq pipeline and see if everything is working, let's try the test dataset.

Remember the key points:

* Start with a fresh new empty directory
* Specify the pipeline version number with `-r`
* Run in the background with `-bg`
* Use the `uppmax` configuration profile to run on UPPMAX
  * If using this, also specify an UPPMAX project with `--project` (two hyphens!)

```
mkdir atacseq-test
cd atacseq-test
nextflow run nf-core/atacseq -r 1.2.1 -profile test,uppmax --project XXX -bg
```

Now, I'll be honest, there's a pretty good chance that something will go wrong at this point.
But that's ok, that's why we run a small test dataset!
This is where you ask for help on Slack instead of suffering in silence.

If all goes well, you should start seeing some log output from Nextflow appearing on your console.
If you run `jobinfo -u USERNAME` you will hopefully see lots of jobs entering the queue.

### Generated files

The pipeline will create a bunch of files in your working directory as it goes:

```console
$ ls -a1
./
../
.nextflow/
.nextflow.log
.nextflow.pid
results/
work/
```

The hidden `.nextflow` files and folders contain information for the cache and detailed logs.

Each task of the pipeline runs in its own isolated directory, these can be found under `work`.
The name of each `work` directory corresponds to the task hash which is listed in the Nextflow log.

As the pipeline runs, it saves the final files it generates to `results` (customise this location with `--outdir`).
Once you are happy that the pipeline has finished properly, you can delete the temporary files in `work`:

```bash
rm -rf work/
```

### Re-running a pipeline with `-resume`

Nextflow is very clever about using cached copies of pipeline steps if you re-run a pipeline.

Once the test workflow has finished, try running the same command again with the `-resume` flag.
Hopefully almost all steps will use the previous cached copies of results and the pipeline
will finish extremely quickly.

This option is very useful if a pipeline fails unexpectedly, as it allows you to start again
and pick up where you left off.

### Read the docs

The documentation for nf-core pipelines is a big part of the community ethos.

Whilst the test dataset is running (it's small, but the UPPMAX job queue can be slow), check out
the nf-core website. Every pipeline has its own page with extensive documentation.
For example, the atacseq docs are at <https://nf-co.re/atacseq>

nf-core pipelines also have some documentation on the command line. You can run this
as you would a real pipeline run, but with the `--help` option.

In a new fresh directory(!), try this out:

```bash
mkdir atacseq-help
cd atacseq-help
nextflow run nf-core/atacseq -r 1.2.1 --help
```

## Running a real workflow

Now we get to the real deal! Once you've gotten this far, you start to leave behind the
generalisations that apply to all nf-core pipelines. Now you have to rely on your wits
and the nf-core documentation.

### Example data

We have prepared some example data for you that comes from the exercises you've worked on
earlier in the week. The files have been subsampled to make them small and quick to run,
and are supplied as gzipped (compressed) FastQ files here: `/sw/courses/epigenomics/nextflow/fastq_sub12_gz/`

Copy this directory to your own working space: (it's 4.2GB so may take a minute or two)

```bash
mkdir epigenomics_nextflow_data
cd epigenomics_nextflow_data
cp -R /sw/courses/epigenomics/nextflow/fastq_sub12_gz/ .
```

### Preparing the sample sheet

The nf-core/chipseq pipeline uses a comma-separated sample sheet file to list all of the input
files and which replicate / condition they belong to.

Take a moment to [read the documentation](https://nf-co.re/chipseq/1.2.1/usage#input) and
make sure that you understand the fields and structure of the file.

We have made a sample sheet for you which describes the different conditions, which can be found
in the folder you copied: `fastq_sub12_gz/exp_design_HeLa_neuron.txt`

_HOWEVER_ - it's not yet complete. The file paths are incorrect and need to be fixed.
The file paths are **relative to where you launch Nextflow**, so will depend a little
on where you copied the data and where you made your new empty directory to run Nextflow.

Remember that you can softlink files in Linux (kind of like making a shortcut link).
This is often a good way of keeping a nicely organised work directory without having to
move the data itself. For example:

```
cd my_analysis
mkdir input_files
cd input_files
ln -s /path/to/my/epigenomics_nextflow_data/*/*.fastq.gz .
```

You can check whether this has worked by using `ls -l`, which should show the full link.
For example:

```console
$ ls -l

total 0
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000OWM.chr12.rmdup.fastq.gz -> ../../fastq_sub12/neural/ENCFF000OWM.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000OWQ.chr12.rmdup.fastq.gz -> ../../fastq_sub12/neural/ENCFF000OWQ.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000OXB.chr12.rmdup.fastq.gz -> ../../fastq_sub12/neural/ENCFF000OXB.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000OXE.chr12.rmdup.fastq.gz -> ../../fastq_sub12/neural/ENCFF000OXE.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PED.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hela/ENCFF000PED.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PEE.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hela/ENCFF000PEE.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PET.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hela/ENCFF000PET.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PMG.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hepg2/ENCFF000PMG.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PMJ.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hepg2/ENCFF000PMJ.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000POM.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hepg2/ENCFF000POM.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000PON.chr12.rmdup.fastq.gz -> ../../fastq_sub12/hepg2/ENCFF000PON.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000RAG.chr12.rmdup.fastq.gz -> ../../fastq_sub12/sknsh/ENCFF000RAG.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000RAH.chr12.rmdup.fastq.gz -> ../../fastq_sub12/sknsh/ENCFF000RAH.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000RBT.chr12.rmdup.fastq.gz -> ../../fastq_sub12/sknsh/ENCFF000RBT.chr12.rmdup.fastq.gz
lrwxrwxrwx 1 phil snic-xxx Nov 25 07:57 ENCFF000RBU.chr12.rmdup.fastq.gz -> ../../fastq_sub12/sknsh/ENCFF000RBU.chr12.rmdup.fastq.gz
```

#### Things to look out for

The following things are easy mistakes when working with chipseq sample sheets - be careful!

* File paths are relative to where you launch Nextflow, not relative to the sample sheet
* Do not have any blank newlines at the end of the file
* Use Linux line endings (`\n`), not windows (`\r\n`)
* If using single end data, keep the empty column for the second FastQ file

### Running the pipeline

Once you've got your sample sheet ready, you can launch the analysis!

Remember the core Nextflow flags that you will need (one hyphen!)

* `-r 1.2.1`
* `-bg`
* `-profile uppmax`

Remember the pipeline paprameter flags that you will need (two hyphens!)

* `--project xxx`
* `--genome GRCh37`
* `--input ./my_samplesheet.csv`
* `--single_end`

If all goes well, your pipeline will run and kick off lots of jobs and merrily process the data!
Once it's finished, take a look in the `results` folder and see what it generated.

## Getting help

Remember that you're not on your own! If you're still struggling after checking the
documentation, jump on to the nf-core Slack and ask for help.

Every pipeline has it's own Slack channel (eg. `#atacseq`, `#chipseq` etc) where people will be happy to help.
