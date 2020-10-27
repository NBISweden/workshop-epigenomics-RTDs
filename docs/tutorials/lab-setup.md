## Setting-up  <a name="Setting-up">

Before we start tutorial, we need to set up our work environment. In particular, we need to:


1. [be able to login to Uppmax and use the node allocation](#UppmaxNode)
<br />
2. [access files prepared for the tutorial](#FilesStructure) and keep folder structures organised
<br />
3. [learn how to read commands and use module system on Uppmax](#CommandsAndModules)
<br />


### Using computational resources <a name="UppmaxNode"></a>

We have booked half a node on Rackham per course participant. To run the tutorial in the interactive mode log to Rackham and run `salloc` command:

```bash
ssh -Y <username>@rackham.uppmax.uu.se
salloc -A g2020XXX -t 04:00:00 -p core -n 8 --no-shell --reservation=g2020XXX_X
```
<font color='red'>where X should be 1 for day one, 2 for day, 3 for day 3 and so on. This gives you access for four hours, so you will repeat this in the afternoon.<font color='black'>

Please make sure you do it only **ONCE**, as by repeating this command you will use up resources reserved for others.


Check which node you were assigned
```bash
$ squeue -u <username>
```


And connect to your node with
```bash
ssh -Y <nodename>
```


### Directory structure <a name="FilesStructure"></a>

There are many files which are part of the data set as well as there are additional files with annotations that are required to run various steps in this tutorial. Therefore saving files in a structured manner is essential to keep track of the analysis steps (and always a good practice). We have preset data access and environment for you. To use these settings run:

* `chiseq_data.sh` that sets up directory structure and creates symbolic links to data as well as copies smaller files **[RUN ONLY ONCE]**

* `chipseq_env.sh` that sets several environmental variables you will use in the exercise: **[RUN EVERY TIME when the connection to Uppmax has been broken, i.e. via logging out]**

Copy the scripts to your home directory and execute them:


```bash
cp /sw/share/compstore/courses/ngsintro/chipseq/scripts/setup/chipseq_data.sh ./
cp /sw/share/compstore/courses/ngsintro/chipseq/scripts/setup/chipseq_env.sh ./

source chipseq_data.sh
source chipseq_env.sh
```

You should see a directory named `chipseq`:

```bash
ls ~
cd ~/chipseq/analysis
```


### Commands and modules <a name="CommandsAndModules"></a>

Many commands are quite long as there are many input files as well as several parameters to set up. Consequently a command may span over several lines of text. The **backslash character 
<!-- ("\\")**  -->
** `\`
indicates to the interpreter that the command input will continue in the following line, not executing the command prematurely.

To see all options for applications used throughout the class type `command -h` to view usage help.

You will notice that you will load and unload modules practically before and after each command. This is done because there are often dependency conflicts between the modules used in this exercise. If not unloaded, some modules will cause error messages from the module system on Rackham. More on module system is here: [https://www.uppmax.uu.se/resources/software/module-system/](https://www.uppmax.uu.se/resources/software/module-system/)

<font color='red'> The commands in this tutorial contain pathways as we set them up in the above steps. If you change file locations, you will need to adjust pathways to match these changes when running commands. </font>

