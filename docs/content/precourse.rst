.. below role allows to use the html syntax, for example :raw-html:`<br />`
.. role:: raw-html(raw)
    :format: html


======================
Pre-course Information
======================


.. .. contents:: 
..     :local:


.. contents:: 
   :depth: 1
   :local:
   :backlinks: none


:raw-html:`<br />`


There are few things to do ***before*** the course starts. Please read carefully and follow the instructions so we can have a good start to the course. Contact us in case anything in unclear.


:raw-html:`<br />`
:raw-html:`<br />`

Computational Resources
=======================

During the course we will be using Uppsala University’s high performance computer cluster (Uppmax) as well as run scripts locally on laptops using ``R`` and ``RStudio``.

To be able to follow exercises we ask you to

- configure access to `Uppmax <https://uppmax.uu.se/>`_;

- install `R <https://cran.r-project.org/>`_ and `RStudio <https://rstudio.com/>`_ on your laptop;

- install `Integrative Genomics Viewer <https://software.broadinstitute.org/software/igv/>`_ on your laptop.

- install the `ThinLinc <https://www.cendio.com/thinlinc/download/>`_ client to easily access uppmax.

Uppmax
------

Computational resources are provided by SNIC / `Uppmax <https://uppmax.uu.se/>`_ . To be able to use them

* create a user account (only if you do not already have one);

* associate your user account with the course project.


Create user account & request membership in the course project (for those new to Uppmax)
*****************************************************************************************

*This section in only related to accepted workshop participants.*


.. This involves few steps. Briefly,

.. * registering at `SUPR <https://supr.naiss.se/>`_;

.. * accepting the User Agreement;

.. * becoming a member in both the storage project ( **NAISS 2024/23-388**) and the compute project ( **NAISS 2024/22-842**);

.. * applying for an account at Uppmax.


.. .. note::

..	To go through the steps follow the instructions on Uppmax website. While at this keep this information handy:

..	* Cluster name: rackham

..	* Project ID:  **NAISS 2024/23-388** and **NAISS 2024/22-842**




Request membership in the course project (for those already having Uppmax account)
***********************************************************************************

*This section in only related to accepted workshop participants.*

.. * log in to SUPR;

.. * under Projects: Requesting Membership in Projects, request membership in **NAISS 2024/23-388** and **NAISS 2024/22-842**


Check configuration (everyone)
*******************************

*This section in only related to accepted workshop participants.*


After you complete setting-up and you receive a notification from SUPR that **your account have been added to the course allocation**

* log in to ``rackham.uppmax.uu.se``

* type ``id`` in the command line

* copy the output of the command and email back (to the course organisers at edu.epigenomics@nbis.se)

A guide on how to `log in for the first time <https://www.uu.se/en/centre/uppmax/get-started/first-login>`_


R and RStudio
---------------

We will also be using the latest version of ``R`` and ``RStudio`` locally. Both of these work on computers running Linux, Windows and Macintosh operating systems. ``RStudio`` is a set of tools as well as an editor that facilitates the use of ``R`` (R ICE). Over the last years it has become a very popular tool and in many ways become a *de-facto* standard for working with ``R``.

Note that on same operative systems it will be easier to install and run ``R`` and ``RStudio`` if you are administrator of your own computer and hence are allowed to install software on your machine. If you do not have these privileges please ask your system administrator to install the latest version of ``R`` and ``RStudio``.


IGV
----

To install, follow the instructions on the `IGV website <https://software.broadinstitute.org/software/igv/>`_.


:raw-html:`<br />`
:raw-html:`<br />`

Further optional preparations
==============================

It might be useful to set up two factor authorization for UPPMAX. This will be necessary to acces the webversion of ThinLinc at https://rackham-gui.uppmax.uu.se and makes it easier to connect from abroad. The instructions to set up 2FA can be found here: https://www.uu.se/en/centre/uppmax/get-started/2-factor.    

For those of you wanting to start ahead and/or brush up on various skills before the course.

:raw-html:`<br />`

Computer skills
------------------

* `Unix <http://www.ee.surrey.ac.uk/Teaching/Unix/>`_: especially the first three chapters.

* DataCamp free `Introduction to R <https://www.datacamp.com/blog/all-about-r>`_. 

* `A short introduction to R <https://cran.r-project.org/doc/contrib/Torfs+Brauer-Short-R-Intro.pdf>`_. A very short introduction to using ``R``.

* How to install and use RStudio from `Data Camp RStudio Tutorial <https://www.datacamp.com/tutorial/r-studio-tutorial>`_.

.. A nice self learn tutorial to ``R``, introducing many central concepts to ``R``.


ChIP-seq
----------

* Introduction to ChIP-seq data analysis `video <https://www.youtube.com/watch?v=zwuUveGgmS0>`_ by Dr. Carl Hermann, University of Heidelberg.

* ChIP-seq and beyond: new and improved methodologies to detect and characterize protein-DNA interactions: `article <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3591838/>`_.

* Q&A: ChIP-seq technologies and the study of gene regulation `article <https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-56>`_.


Software Dependencies
=====================

If you have no access to Uppmax, where all software is preinstalled, you can configure your own system to follow the exercises.

The dependencies are listed in :doc:`Dependencies <./dependencies>` .

