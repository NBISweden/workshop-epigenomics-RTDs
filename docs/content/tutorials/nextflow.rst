Nextflow and nf-core
====================

**Learning outcomes**

* To understand the basic components in a Nextflow run
* To understand the difference between Nextflow and nf-core
* To be able to find relevant nf-core pipelines

.. Contents
.. ========

.. contents:: 
    :local:


.. attention::

    For this tutorial **DO NOT** run the ``interactive`` command! These nextflow exercises are run from the login node (so just ``ssh -Y`` to rackham). If you later today want to work on some earlier exercises you can request cores with ``interactive`` after this tutorial.

tmux
-----

Before starting the tutorial it might be useful to introduce you to **tmux**. **Tmux** is a very useful little tool that will allow you to start a process and run it in the background, allowing you to use the terminal for other processes during long calculations. As an added bonus, it will keep your processes going if you leave the server or your connection is unstable and crashes. It first needs to be loaded using UPPMAX's module system, after which you can initiate a new terminal in tmux.

.. attention::

    Do these next steps **after** connecting to rackham!

.. code-block:: bash

    module load tmux
    tmux new -s nf_tutorial # or any other name you like
    tmux set mouse on # enable mouse support for things like scrolling and selecting text

Now, anything you do in this new tmux terminal session is "protected". When the connection to the server crashes mid-session, just reconnect to UPPMAX and do

.. code-block:: bash

    module load tmux
    tmux attach -t nf_tutorial

To put **tmux** in background and keep the processes inside running, press ``Ctrl+b``, release, press ``d``. With ``tmux ls`` you can see which sessions are ongoing (can be multiple ones!) and you could connect to. To reattach to your earlier session type ``tmux attach -t nf_tutorial`` as shown above.

Each session is persistent and will survive accidental disconnection (such as ssh(1) connection timeout) or intentional detaching (with the `C-b d` key strokes). To kill a **tmux** session and stop any process running in it, press ``Ctrl+b``, release, press ``x`` followed by ``y``.

All of this might seem to add unnecessary hassle but tmux is extremely valuable when working on a server. Instead of having to redo a long list of computational step when the connection to a server inevitably crashes, just reconnect to the ongoing tmux session and you are back exactly where you were when the crash happened. Tmux actually can do even more useful things, so if you want to know more, have a look at this `quick and easy guide to tmux <https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/>`_.


Setup: Nextflow
----------------

These steps typically only need to be done once, to get you set up with Nextflow.


Installation
#############

To use Nextflow, either install it yourself in your home directory or use the UPPMAX environment module.

To use the environment modules (preferable when on UPPMAX - do all of the steps in this tutorial in the earlier started tmux session):

.. code-block:: bash

    module purge                        # Clears all existing loaded modules, to start fresh
    module load uppmax bioinfo-tools    # Base UPPMAX environment modules, needed for everything else
    module load Nextflow/22.10.1                # Note: Capital N!


Alternatively, to install yourself (when not on UPPMAX for example):

.. code-block:: bash

    cd ~/bin    # Your home directory bin folder - full of binary executable files, already on your PATH
    curl -s https://get.nextflow.io | bash

This code won't work on Uppmax or another HPC! This is strictly for local installations; you may need admin priviledges.


Environment variable setup
###########################

Nextflow has a large list of bash environment variables that can be set to configure how it runs.

.. note::

    If you don't want to enter these commands every time you log in, the most convenient way to set these is to add them to the end of your ``.bashrc`` file in your home directory. Once here, they will be applied every time you log in automatically. [You don't need to do that for this exercise session]

.. code-block:: bash

    # Don't let Java get carried away and use huge amounts of memory
    export NXF_OPTS='-Xms1g -Xmx4g'

    # Don't fill up your home directory with cache files
    export NXF_HOME=$HOME/nxf-home
    export NXF_TEMP=${SNIC_TMP:-$HOME/glob/nxftmp}

Upon execution of the command, ``$USER`` will be replaced with your login user name.



Check that Nextflow works
#########################

It’s always good to have a mini test to check that everything works.

These pipelines can create large temporary files and large result files, so we will do these exercises in the course project directory. Make a new directory there and run the Nextflow test command as follows:

.. code-block::

    mkdir /proj/epi2025/nobackup/private/$USER # create personal folder in project directory
    cd /proj/epi2025/nobackup/private/$USER
    mkdir nextflow-hello-test
    cd nextflow-hello-test
    nextflow run hello

You should see something like this:

.. code-block:: bash

    N E X T F L O W  ~  version 22.10.2
    Pulling nextflow-io/hello ...
    downloaded from https://github.com/nextflow-io/hello.git
    Launching `nextflow-io/hello` [sharp_sammet] - revision: 96eb04d6a4 [master]
    executor >  local (4)
    [7d/f88508] process > sayHello (4) [100%] 4 of 4 ✔
    Bonjour world!

    Ciao world!

    Hello world!

    Hola world!

Succes!

Setup: nf-core
---------------

**nf-core** is a collection of best practice pipelines for use in genomics, open-source and with extensive `documentation <https://nf-co.re/>`_ .


Recently, all nf-core pipelines have been made available on UPPMAX (Rackham and Bianca) so they can be run on these servers without any additional setup besides loading the ``nf-core-pipelines`` module.


.. code-block:: bash

    module load nf-core-pipelines/latest

Loading this module exposes the variable ``$NF_CORE_PIPELINES``. This is the location on the server where all pipelines are stored. Have a look at all pipelines and versions that are available

.. code-block:: bash

    tree -L 2 $NF_CORE_PIPELINES -I 'singularity_cache_dir'

This directory also contains all necessary software for all pipelines in a folder called ``singularity_cache_dir``. This means you do not have to install any tools at all; they all are here packaged in singularity containers!

.. note::

    nf-core also comes as a Python package that is totally separate to Nextflow and is not required to run Nextflow pipelines. It does however offer some convenience functions to make your life a little easier. A description on how to install this package can be found `here <https://github.com/nf-core/tools>`_. This is useful if you want to run nf-core pipelines outside of UPPMAX or want to use some of the convenience functions included in the nf-core package. [not necessary for running the current exercises on UPPMAX; but the students not on UPPMAX might give this a try]


Running a test workflow
------------------------

It’s always a good idea to start working with a tiny test workflow when using a new Nextflow pipeline. This confirms that everything (environmental variables, nextflow config, metadata, ...) is set up and working properly, before you start moving around massive data files. To accommodate this, all nf-core pipelines come with a configuration profile called test which will run a minimal test dataset through the pipeline without needing any other pipeline parameters.



Trying out atacseq
####################

To try out for example the nf-core `atacseq <https://nf-co.re/atacseq/2.1.2/>`_ pipeline and see if everything is working, let’s try the test dataset.

Remember the key points:

* Start with a fresh new empty directory
* ``$NF_CORE_PIPELINES`` specifies the path where all pipelines are stored
* Specify the pipeline with ``$NF_CORE_PIPELINES/[name]/[version]/workflow``
* Use the ``uppmax`` configuration profile to run on UPPMAX from a login node
    - If using this, also specify an UPPMAX compute allocation with ``--project`` (two hyphens!)
* Use the test configuration profile to run a small test


.. code-block:: bash

    cd /proj/epi2025/nobackup/private/$USER
    mkdir atacseq-test
    cd atacseq-test
    nextflow run $NF_CORE_PIPELINES/atacseq/1.2.1/workflow -profile test,uppmax --project uppmax2025-2-309

There’s a pretty good chance that something will go wrong at this point. But that’s ok, that’s why we run a small test dataset. This is where you check the configuration, environment, and ask for help on Slack.

If all goes well, you should start seeing some log output from Nextflow appearing on your console (in standard output). Nextflow informs you which step of the pipeline it is doing and the percentage completed.


After the completion you'll see message similar to:

.. code-block:: bash
    
    [nf-core/atacseq] Pipeline completed successfully-
    Completed at: 18-Sep-2024 13:21:44
    Duration    : 5m 17s
    CPU hours   : 0.5
    Succeeded   : 88


Even though the datasets in a test run are small, this pipeline can take a while because it submits jobs to the UPPMAX server via the resource manager SLURM. Depending on how busy the server is at the moment (and it might be quite busy if you all run this at the same time!), it may take a while before your jobs are executed. It might therefore be necessary to cancel the pipeline once Nextflow seems to progress though the different steps slowly but steadily. If you want to cancel the pipeline execution to progress with the tutorial, press ``CTRL-C`` which cancels the currently running process. Or alternatively, put it in the background using tmux, do some other things and reattach later to check in on the progress.



Generated files
################

The pipeline will create some files in your directory:

.. code-block:: bashrc

    $ ls -a1
    ./
    ../
    .nextflow/
    .nextflow.log
    .nextflow.pid
    results/
    work/

The hidden ``.nextflow`` files and folders contain information for the cache and detailed logs.

Each task of the pipeline runs in its own isolated directory, these can be found under ``work``. The name of each ``work`` directory corresponds to the task hash which is listed in the Nextflow log.

As the pipeline runs, it saves the final files it generates to ``results`` (customise this location with ``--outdir``). Once you are happy that the pipeline has finished properly, you can delete the temporary files in ``work``:

.. code-block:: bash

    rm -rf work/



Re-running a pipeline with ``-resume``
#######################################

Nextflow is capable of using cached copies of pipeline steps if you re-run a pipeline (i.e. resume using the same starting data).

Once the test workflow has finished or you have canceled it the middle of its execution, try running the same command again with the ``-resume`` flag. Hopefully almost all steps will use the previous cached copies of results and the pipeline will finish extremely quickly.

This option is very useful if a pipeline fails unexpectedly, as it allows you to start again and pick up where you left off.



Read the docs
##############

The documentation for nf-core pipelines is a big part of the community ethos.

Whilst the test dataset is running (it’s small, but the UPPMAX job queue can be slow), check out the nf-core website. Every pipeline has its own page with extensive documentation. For example, the atacseq docs are at https://nf-co.re/atacseq

nf-core pipelines also have some documentation on the command line. You can run this as you would a real pipeline run, but with the ``--help`` option.

In a new fresh directory(!), try this out:

.. code-block:: bash

    cd /proj/epi2025/nobackup/private/$USER
    mkdir atacseq-help
    cd atacseq-help
    nextflow run $NF_CORE_PIPELINES/atacseq/1.2.1/workflow --help


For technical reasons, there cannot be two *simultaneous* nextflow runs in the same working directory. This is to prevent conflicts between different runs.


Running a real workflow
-----------------------------

Now we get to the real deal! Once you’ve gotten this far, you start to leave behind the generalisations that apply to all nf-core pipelines. Now you have to rely on your wits and the nf-core documentation. We have prepared small datasets for a ChIP-seq analysis and a BS-seq analysis. You can choose to do the one that interests you most or if you have time you can try both!


ChIP-seq
---------

Example data
##############

We have prepared some example data for you that comes from the exercises you’ve worked on earlier in the week. The files have been subsampled to make them small and quick to run, and are supplied as gzipped (compressed) FastQ files here: ``/sw/courses/epigenomics/nextflow/fastq_sub12_gz/``

Make a new directory for this ChIP seq analysis and link the data files to a data folder in this directory. We link to these files in this tutorial instead of copying them (which would also be an option) so as not to fill up the filesystem.

.. code-block:: bash

    cd /proj/epi2025/nobackup/private/$USER
    mkdir chip_seq_analysis
    cd chip_seq_analysis
    mkdir input_files
    cd input_files
    ln -s /sw/courses/epigenomics/nextflow/fastq_sub12_gz/neural/*.fastq.gz .
    ls

The last command should show you the 4 neural fastq.gz files in this folder.


Preparing the sample sheet
###########################

The nf-core/chipseq pipeline uses a comma-separated sample sheet file to list all of the input files and which replicate / condition they belong to.

Take a moment to read the documentation and make sure that you understand the fields and structure of the file.

We have made a sample sheet for you which describes the different condition: ``samplesheet.csv``. Copy it to you ``chip_seq_analysis`` folder.

.. code-block:: bash

    cd .. # move up one directory
    cp /sw/courses/epigenomics/nextflow/samplesheet.csv .
    cat samplesheet.csv


The cat command shows you the contents of the sample sheet.


Things to look out for
#######################

The following things are easy mistakes when working with chipseq sample sheets - be careful!

* File paths of the ``fastq.gz`` files are relative to where you launch Nextflow (i.e. the ``chip_seq_analysis`` folder), not relative to the sample sheet
* Do not have any blank newlines at the end of the file
* Use Linux line endings (``\n``), not Windows (``\r\n``) - in practice this means be mindful of this when creating this file using spreadsheet applications e.g. MS Word, MacOS Numbers etc.
* If using single end data, keep the empty column for the second FastQ file


Running the pipeline
#####################

Once you’ve got your sample sheet ready, you can launch the analysis. For this, try to figure out the command you should run standing at ``chip_seq_analysis`` directory. Try to execute the chipseq pipeline version 1.2.2 using the FastQ files you just linked to.

Remember the core Nextflow flags that you will need (one hyphen!)

* ``-profile uppmax``

Remember the pipeline specific parameter flags that you will need (two hyphens!)

* ``--project uppmax2025-2-309``
* ``--genome GRCh38``
* ``--input samplesheet.csv``
* ``--single_end``

If all goes well, your pipeline will run and kick off lots of jobs and merrily process the data! Once it’s finished, take a look in the ``results`` folder and see what it generated. Again, this might take a while due to the job queue (1 hour +), so feel free to detach from the tmux session and return later.

.. admonition:: ChIP command
   :class: dropdown, note

    ``nextflow run $NF_CORE_PIPELINES/chipseq/1.2.2/workflow -profile uppmax --project naiss2024-22-842 --genome GRCh38 --input samplesheet.csv --single_end``


Standard output of this run::

    
    N E X T F L O W  ~  version 22.10.2
    Launching `/sw/bioinfo/nf-core-pipelines/latest/rackham/chipseq/1.2.2/workflow/main.nf` [jolly_hoover] DSL1 - revision: 39c073c68c
    ----------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/chipseq v1.2.2
    ----------------------------------------------------

    Run Name            : jolly_hoover
    Data Type           : Single-End
    Design File         : samplesheet.csv
    Genome              : GRCh38
    Fasta File          : /sw/data/igenomes//Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
    GTF File            : /sw/data/igenomes//Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
    Gene BED File       : /sw/data/igenomes//Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.bed
    BWA Index           : /sw/data/igenomes//Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa
    Blacklist BED       : /sw/bioinfo/nf-core-pipelines/latest/rackham/chipseq/1.2.2/workflow/assets/blacklists/hg38-blacklist.bed
    MACS2 Genome Size   : 2.7e9
    Min Consensus Reps  : 1
    MACS2 Narrow Peaks  : No
    MACS2 Broad Cutoff  : 0.1
    Trim R1             : 0 bp
    Trim R2             : 0 bp
    Trim 3' R1          : 0 bp
    Trim 3' R2          : 0 bp
    NextSeq Trim        : 0 bp
    Fragment Size       : 200 bp
    Fingerprint Bins    : 500000
    Save Genome Index   : Yes
    Max Resources       : 500 GB memory, 16 cpus, 10d time per job
    Container           : singularity - nfcore/chipseq:1.2.2
    Output Dir          : ./results
    Launch Dir          : /crex/proj/epi2023/nobackup/private/agata/chipseq
    Working Dir         : /crex/proj/epi2023/nobackup/private/agata/chipseq/work
    Script Dir          : /sw/bioinfo/nf-core-pipelines/latest/rackham/chipseq/1.2.2/workflow
    User                : agata
    Config Profile      : uppmax
    Config Description  : UPPMAX (Bianca) cluster profile provided by nf-core/configs.
    Config Contact      : Phil Ewels (@ewels)
    Config URL          : https://www.uppmax.uu.se/



Methyl-seq
-----------

nf-core/methylseq is an analysis pipeline used for methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and extract methylation calls and performs extensive quality-control on the results. The default workflow uses Bismark with Bowtie2 as alignment tool: unless specified otherwise, nf-core/methylseq will run this pipeline.

Example data
##############

We have prepared some example data that has been subsampled to make them small and quick to run, and are supplied as gzipped (compressed) FastQ files here: ``/sw/courses/epigenomics/DNAmethylation/pipeline_bsseq_data/Sample1_PE_R[1,2].fastq.gz``. This is mouse data so remember to use the correct genome to map to.

Running the pipeline
#####################

Begin with making a fresh analysis directory in your home directory

.. code-block:: bash

    cd /proj/epi2025/nobackup/private/$USER
    mkdir methylseq_analysis
    cd methylseq_analysis

In this folder you can launch the analysis! For this, try to figure out the command you should run. Try to execute the methylseq pipeline with version 1.6.1 using the FastQ files mentioned above.

Remember the core Nextflow flags that you will need (one hyphen!)

* ``-profile uppmax``

Figure out the pipeline specific parameter flags that you will need (two hyphens!). Have a look at the `list of parameters <https://nf-co.re/methylseq/1.6.1/parameters>`_ to get an idea which options are possible and make sure to use the essential parameters.

* ``--input '/sw/courses/epigenomics/DNAmethylation/pipeline_bsseq_data/Sample1_PE_R{1,2}.fastq.gz'``
* ``--aligner bismark``
* ``--project uppmax2025-2-309``
* ``--genome mm10``

If all goes well, your pipeline will run and kick off lots of jobs and merrily process the data! Once it’s finished, take a look in the ``results`` folder and see what it generated. A description of the outputs can be seen `here <https://nf-co.re/methylseq/1.6.1/docs/output>`_.  Again, this might take a while due to the job queue (1 hour +), so feel free to detach from the tmux session and return later.

.. admonition:: minimal methylseq command
   :class: dropdown, note

   ``nextflow run $NF_CORE_PIPELINES/methylseq/1.6.1/workflow -profile uppmax --input '/sw/courses/epigenomics/DNAmethylation/pipeline_bsseq_data/Sample1_PE_R{1,2}.fastq.gz' --aligner bismark --project naiss2024-22-842 --genome mm10``


Standard output of this run::


    N E X T F L O W  ~  version 22.10.2
    Launching `/sw/bioinfo/nf-core-pipelines/latest/rackham/methylseq/1.6.1/workflow/main.nf` [elated_shannon] DSL1 - revision: 18179884b6


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/methylseq v1.6.1
    ------------------------------------------------------

    Core Nextflow options
      runName                   : elated_shannon
      containerEngine           : singularity
      container                 : nfcore/methylseq:1.6.1
      launchDir                 : /crex/proj/epi2023/nobackup/private/agata/methylseq
      workDir                   : /crex/proj/epi2023/nobackup/private/agata/methylseq/work
      projectDir                : /sw/bioinfo/nf-core-pipelines/latest/rackham/methylseq/1.6.1/workflow
      userName                  : agata
      profile                   : uppmax
      configFiles               : /sw/bioinfo/nf-core-pipelines/latest/rackham/methylseq/1.6.1/workflow/nextflow.config

    Input/output options
      input                     : /sw/courses/epigenomics/DNAmethylation/pipeline_bsseq_data/Sample1_PE_R{1,2}.fastq.gz

    Reference genome options
      genome                    : mm10
      fasta                     : /sw/data/igenomes//Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
      bismark_index             : /sw/data/igenomes//Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/
      save_reference            : true
      igenomes_base             : /sw/data/igenomes/

    Generic options
      max_multiqc_email_size    : 25 MB

    Max job request options
      max_memory                : 500 GB
      max_time                  : 10d

    Institutional config options
      custom_config_base        : /sw/bioinfo/nf-core-pipelines/latest/rackham/methylseq/1.6.1/workflow/../configs/
      config_profile_description: UPPMAX (Bianca) cluster profile provided by nf-core/configs.
      config_profile_contact    : Phil Ewels (@ewels)
      config_profile_url        : https://www.uppmax.uu.se/

    ------------------------------------------------------
     Only displaying parameters that differ from defaults.
    ------------------------------------------------------
    WARN: Access to undefined parameter `bwa_meth_index` -- Initialise it to a default value eg. `params.bwa_meth_index = some_value`
    WARN: Access to undefined parameter `fasta_index` -- Initialise it to a default value eg. `params.fasta_index = some_value`



Getting help
-------------

Please have a look at the nf-core website to see which pipelines are available (66 as of today!) and browse their thorough documentation.

Remember that you’re not on your own! If you’re still struggling after checking the documentation, jump on to the nf-core Slack and ask for help.

Every pipeline has it’s own Slack channel (eg. ``#atacseq``, ``#chipseq`` etc) where people will be happy to help.


How to modify pipeline parameters
-----------------------------------------

**not a tutorial**

Some parameters can be specified from the command line, as described in pipeline documentation, e.g. for `atacseq pipeline <https://nf-co.re/atacseq/2.1.2/docs/usage/#running-the-pipeline>`_. You can also provide a project-specific config file ``params.yaml``  with the parameters and their values.

What happens though if you would like to make changes to parameters not available via command line? This usually means introducing changes to actual command calls for tools of interest. If this is what you want to do, read on - only if you know what you're doing (you should have some idea after this workshop).

This is not a tutorial, merely hints on how to go about customizing nf-core pipelines beyond documented command line changes.

Generally, the code for nf-core pipelines can look quite complicated, with lots of options, parameters and control steps, several config files and profiles, modules, subworkflows... This structure is practical for developers, but also makes them less readable for a casual user. Fortunately though, tool running parameters are exposed in a dedicated config file ``conf/modules.config``, so there is no need to read all the code. 

First off, all edits of this nature should be done on **your own copy of the pipeline code**. Ideally, this would be a github repository with the fork of the original pipeline (version control!). If this seems daunting, you can simply clone the repository to your directory on HPC, make necessary modifications, and run the pipeline from this new location rather than from server-wide source directory.

For example, if we want to change read mapping parameters in the atacseq pipeline (https://github.com/nf-core/atacseq/), the file we need to edit is its ``conf/modules.config`` (https://github.com/nf-core/atacseq/blob/master/conf/modules.config). 

First, we obtain our local copy of this pipeline::

    git clone https://github.com/nf-core/atacseq.git

We then locate the file to be edited, open it in text editor (e.g. ``nano``, ``vi``) and introduce our edits.

Let's say we want to change read alignment command for ``bowtie2``: we want to add option ``--dovetail`` to the existing command. Details on `dovetail <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other>`_ can be found in bowtie2 manual.

We search for arguments defined for process ``BOWTIE2_ALIGN`` (line 258), and we add ``--dovetail`` to existing ``ext.args``. 

After this, the pipeline can be run like this::

   nextflow run /path/to/new-atacseq/main.nf -profile uppmax --project naiss2024-22-842 --genome GRCh38 --input samplesheet.csv --aligner bowtie2  <other command line parameters you want to specify>

Naturally, it is advisable to first **test the performance and verify that the output is what we intended**, on a small test data.


