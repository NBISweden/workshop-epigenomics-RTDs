.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html

===============
Setting-up
===============



Before we start tutorial, we need to set up our work environment. In particular, we need to:


1. `be able to login to Uppmax and use the node allocation <Using computational resources>`_
2. `access files prepared for the tutorial <Directory structure>`_ and keep folder structures organised
3. `learn how to read commands and use module system on Uppmax <Commands and modules>`_
4. if you work without access to Uppmax: `check the dependencies <Software Dependencies>`_

:raw-html:`<br />`


Using computational resources
==========================================================

We have booked half a node on Rackham per course participant. To run the tutorial in the interactive mode log to Rackham and run ``salloc`` command:

.. code-block:: bash

	ssh -Y <username>@rackham.uppmax.uu.se
	salloc -A g2020XXX -t 04:00:00 -p core -n 8 --no-shell --reservation=g2020XXX_X


.. HINT::

	where X should be 1 for day one, 2 for day, 3 for day 3 and so on. This gives you access for four hours, so you will repeat this in the afternoon.


Please make sure you do it only **ONCE**, as by repeating this command you will use up resources reserved for others.


Check which node you were assigned

.. code-block:: bash

	squeue -u <username>



And connect to your node with

.. code-block:: bash

	ssh -Y <nodename>



Directory structure
======================

There are many files which are part of the data set as well as there are additional files with annotations that are required to run various steps in this tutorial. Therefore saving files in a structured manner is essential to keep track of the analysis steps (and always a good practice). We have preset data access and environment for you.


Commands and modules
======================

Many commands are quite long as there are many input files as well as several parameters to set up. Consequently a command may span over several lines of text. The (backslash character ``\`` indicates to the interpreter that the command input will continue in the following line, not executing the command prematurely.

.. HINT::

	To see all options for applications used throughout the class type ``command -h`` to view usage help.

You will notice that you will load and unload modules practically before and after each command. This is done because there are often dependency conflicts between the modules used in this exercise. If not unloaded, some modules will cause error messages from the module system on Rackham. More on `module system used on Uppmax <https://www.uppmax.uu.se/resources/software/module-system/](https://www.uppmax.uu.se/resources/software/module-system/>`_.

.. HINT::

	The commands in this tutorial contain paths as we set them up in the above steps. If you change file locations, you will need to adjust pathways to match these changes when running commands.



Software Dependencies
=====================

If you have no access to Uppmax, where all software is preinstalled, you can configure your own system to follow the exercises.

The dependencies are listed in :doc:`Dependencies <./dependencies>` .

