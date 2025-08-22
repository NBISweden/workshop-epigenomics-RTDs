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


There are few things to do ***before*** the course starts. Please read carefully and follow the instructions so we can have a good start to the course. Contact us in case anything is unclear.

:raw-html:`<br />`
:raw-html:`<br />`

Computational Resources
=======================

During the course we will be using Uppsala University’s high performance computer cluster (`UPPMAX <https://docs.uppmax.uu.se/cluster_guides/uppmax_cluster/>`_)  
as well as run scripts locally on laptops using ``R`` and ``RStudio``.  


Register at SUPR
-----------------

A SUPR/NAISS account is needed to use UPPMAX resources. 
If you do not already have one, please register for an account at `SUPR/NAISS <https://supr.naiss.se/person/register/>`_ 
and then you must accept the `user agreement <https://supr.naiss.se/person/user_agreement/naiss/>`_ either online or in paper form.  
For those within Swedish academia, it is recommended to register with your SWAMID. 
If you don’t have a SWAMID connected account, you will have to send in your signed user agreement in paper form together with a copy of your passport. 
Note that this manual process can take a week or more.
You can follow this detailed `instruction <https://docs.uppmax.uu.se/getting_started/supr_register/>`_.  Once you have a SUPR/NAISS account, please 
`request <https://docs.uppmax.uu.se/getting_started/join_existing_project/>`_ membership to the project **UPPMAX 2025/2-309**.



Apply for an UPPMAX user account
--------------------------------

After you have been added to a SUPR project, you can now apply for a user account in UPPMAX if you do not have an existing account yet.
This will be the account you use to log in to UPPMAX so it is not the same as your SUPR account.  
Note that it might take up to 2 working days for your account to be created.  You will then receive 2 emails with information on how to login to UPPMAX.
To apply, please follow this `instruction <https://docs.uppmax.uu.se/getting_started/user_account/#apply-for-an-account-at-uppmax>`_.

:raw-html:`<br />`
:raw-html:`<br />`


Connect to UPPMAX cluster Rackham
=================================

There are three ways to connect to Rackham: a text-based SSH connection from a terminal, a graphical remote desktop using a Thinlinc client on your laptop and
a remote desktop environment via a web browser.

`SSH connection using a terminal <https://docs.uppmax.uu.se/getting_started/login_rackham_console_password/>`_
----------------------------------------------------------------------------------------------------------------------------

A straightforward way to connect to  Rackham is through a terminal using ``ssh`` connection.

- For Linux users, use Terminal (included by default).
- For Mac users, use Terminal (included by default).  You need to also install `XQuartz  <https://www.xquartz.org/>`_  to enable X11 forwarding from a terminal, 
  i.e., to run graphical applications on a remote server while displaying them on your laptop.

- For Windows system, we recommend using `MobaXterm <http://mobaxterm.mobatek.net/>`_.  Please install the program and not use the portable version.
  MobaXterm also has a built-in X11 and an integrated SFTP file browser to transfer data between the remote cluster and your computer.

Other SSH clients and X-servers for MacOS and Windows are described `here <https://docs.uppmax.uu.se/software/ssh_client/#ssh-client>`_.

`GUI desktop using ThinLinc locally <https://docs.uppmax.uu.se/getting_started/login_rackham_remote_desktop_local_thinlinc_client/>`_
----------------------------------------------------------------------------------------------------------------------------------------

You can also connect to Rackham via a graphical remote desktop using a local ThinLinc client on your laptop.  
Thinlinc is useful if you need to view images or documents in GUI programs without having to first download the image/document to your own computer first. 
Since it is using graphics, it will require you to have an internet connection that is good and stable.  

- Install `ThinLinc <https://www.cendio.com/thinlinc/download/>`_ .

`GUI desktop using ThinLinc in a web browser <https://docs.uppmax.uu.se/getting_started/login_rackham_remote_desktop_website/>`_
---------------------------------------------------------------------------------------------------------------------------------

Lastly,  you can connect to Rackham using a remote desktop environment via a web browser at `https://rackham-gui.uppmax.uu.se/ <https://rackham-gui.uppmax.uu.se/>`_. 
This is the easiest option as no software installation is needed but it will give you the slowest connection.  It will require you to have an internet connection that is good and stable.
For this option,  you have to set up your UPPMAX two factor authorization (2FA). 
See `how to get an UPPMAX 2FA <https://docs.uppmax.uu.se/getting_started/get_uppmax_2fa/>`_.

:raw-html:`<br />`

This `tutorial <https://docs.uppmax.uu.se/getting_started/login_rackham/>`_ will guide you to connect to Rackham using the three methods as described above.  
For Windows users with MobaXterm, please follow this `instruction <https://docs.uppmax.uu.se/software/ssh_client/#using-ssh-with-different-terminals-that-allow-for-graphics>`_.


:raw-html:`<br />`

Some useful tutorials:

- `Troubleshoot your MobaXterm X11 connection to UPPMAX <https://hackmd.io/@pmitev/UPPMAX-MobaXterm-X11>`_
- `File transfer to/from Rackham <https://docs.uppmax.uu.se/cluster_guides/transfer_rackham/>`_
- `Change your UPPMAX password <https://docs.uppmax.uu.se/getting_started/change_uppmax_password/>`_
- `How to get an UPPMAX 2FA <https://docs.uppmax.uu.se/getting_started/get_uppmax_2fa/>`_.


   



:raw-html:`<br />`



Check configuration
--------------------

After you complete setting-up and you receive a notification from SUPR that **your account have been added to the course allocation**

* log in to ``rackham.uppmax.uu.se``

* type ``id`` in the command line

* copy the output of the command and email back (to the course organisers at edu.epigenomics@nbis.se)


:raw-html:`<br />`
:raw-html:`<br />`



Install Tools
=================

To be able to follow exercises we ask you to

- install `R <https://cran.r-project.org/>`_ and `RStudio <https://rstudio.com/>`_ on your laptop.

- install `Integrative Genomics Viewer <https://software.broadinstitute.org/software/igv/>`_ on your laptop.

We will also be using the latest version of ``R`` and ``RStudio`` locally. Both of these work on computers running Linux, Windows and Macintosh operating systems. ``RStudio`` is a set of tools as well as an editor that facilitates the use of ``R`` (R ICE). Over the last years it has become a very popular tool and in many ways become a *de-facto* standard for working with ``R``.

Note that on same operative systems it will be easier to install and run ``R`` and ``RStudio`` if you are administrator of your own computer and hence are allowed to install software on your machine. If you do not have these privileges please ask your system administrator to install the latest version of ``R`` and ``RStudio``.




:raw-html:`<br />`
:raw-html:`<br />`

Further optional preparations
==============================


For those of you wanting to start ahead and/or brush up on various skills before the course.


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

