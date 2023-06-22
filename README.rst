KATObitPipe
===========

Obit based calibration and imaging scripts for MeerKAT.

Quick Install Guide
*******************

Download from github and pip install the package on your own machine with an installation and setup of `Obit < https://www.cv.nrao.edu/~bcotton/Obit.html>`_ . For help setting up Obit for yourself from a binary distribution see the section :ref:`ComInstall` below.

.. code-block::
	$ git clone https://github.com/mauch/KATObitPipe
	$ pip install ./KATObitPipe

Run the Calibration Pipeline
****************************

Once the Obit environment is setup and KATObitPipe is installed as above you can run the calibration script with

.. code-block::
	$ KATCalibPipe.py RDB

The `RDB` can be an https URL to an observation in the MeerKAT archive (with associated token) or simply a `.rdb` filename downloaded from the archive.

There are numerous options to the script, doing

.. code-block::
	$ KATCalPipe.py --help

will list them and describe what they do.

A typical example of a run on the observation with CBID 

.. _ComInstall:
Installing Obit binary distribution on the SARAO comXX machines
***************************************************************