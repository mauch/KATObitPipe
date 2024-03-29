KATObitPipe
===========

Obit based calibration and imaging scripts for MeerKAT.

Quick Install Guide
-------------------

Create a working install directory and install from there eg.:

```
	$ mkdir ~/Obit
	$ cd ~Obit
```

Download from github and pip install the package on your own machine with an installation and setup of [Obit](https://www.cv.nrao.edu/~bcotton/Obit.html). For help setting up Obit for yourself from a binary distribution see the section [Installing Obit binary distribution on the SARAO comXX machines](https://github.com/mauch/KATObitPipe/edit/master/README.rst#installing-obit-binary-distribution-on-the-sarao-comxx-machines) below.

```
	$ git clone https://github.com/mauch/KATObitPipe
	$ pip install ./KATObitPipe
```

Run the Calibration Pipeline
----------------------------

Once the Obit environment is setup, usually using `source setup.sh` with `setup.sh` containing environment variables pointing to the appropriate location of the Obit installation (see [setup.sh](/setup.sh) for an example). KATObitPipe should be installed as above. After this you can run the calibration script from anywhere.

You will need disk space in the current working directory when you run the script in order to store the data as it is copied over from the archive. A typical 8 hour, 4096 channel, 8s dump observation will require approx 1.5 TB of disk space for processing.

In your desired run location you run the calibration script with:

```
	$ KATCalPipe.py RDB
```

The `RDB` can be an https URL to an observation in the MeerKAT archive (with associated token) or simply a `.rdb` filename downloaded from the archive.

There are numerous options to the script, doing

```
	$ KATCalPipe.py --help
```

will list them and describe what they do.

A typical example of a run on the observation with CBID=1592263862 stored in a file downloaded from the archive and renamed to `1592263862_sdp_l0.full.rdb` is:

```
	$ KATCalPipe.py --flag --gzip --polcal 1592263862_sdp_l0.full.rdb
```
This will run the calibration pipeline, re-flag the data (ie. throw away all the cal_rfi flags and recompute them on the fly), gzip the output UV table file and run in `polcal` mode, which means downloading the associated delaycal observation and performing XY-Phase calibration on it before doing the full calibration on the average of the H&V polarisations.

Using the option `--reuse` you can bypass the downloading of the data into the aipsdisk from the archive, and use the data already downloaded from a previous run. This takes the data from *after* the Hanning step.

Typical Output of the Pipeline
------------------------------

A directory listing of what you should see after the pipeline has run with the example above is:
```
1592263862.CalTab.uvtab  1592263862.refAnt.pickle  1592263862_APCal2.ps     1592263862_Spec.ps          da00
1592263862.Parms.pickle  1592263862.uvtab.gz       1592263862_DelayCal.ps   1592263862_sdp_l0.full.rdb  aipsdisk
1592263862.log           1592263862_APCal.ps       1592263862_DelayCal2.ps  31DEC22
```

`1592263862.CalTab.uvtab` : Contains the AIPS calibration tables that were derived and applied to the data.

`1592263862.uvtab.gz` : Countains the calibrated UV data in UVTab format (it is gzipped because of the `--gzip` option used in the example above.

`*.pickle` : Serialised outputs of pipeline parameters `Parms.pickle` and derived reference antenna `refAnt.pickle`.

`*.ps` : Plots of delays, bandpasses, Gains derived during the pipeline run.

`1592263862.log` : logfile output by the run

`aipsdisk` : This cotains the UV data in AIPS UV format, and can be used to run the pipeline in `--reuse` mode.

`31DEC22 da00` : These are the residual AIPS run files used by the AIPS tasks in the pipeline.

Easy Installation and Running using Docker
------------------------------------------

You can build and run a Docker image with the appropriate Obit and `KATObitPipe` scrips installed by:

1. Get the latest `KATObitPipe` from github:
```
	$ git clone https://www.github.com/mauch/KATObitPipe
```

2. Change directory into the repo
```
	$ cd KATObitPipe
```
4. Build the docker image (you will find out if your docker installation is working at this point)
```
	$ docker build -t katobitpipe .
```

5. To RUN the script in the docker image cd to where you want to run (where your .rdb file and parameter file is) - you must have write permission here (obviously). Then the command is:
```   
	$ docker run -t --rm -v ${PWD}:/scratch -e LOCAL_USER_ID=$(id -u) -e LOCAL_GROUP_ID=$(id -g) katobitpipe KATCalPipe.py <RDB URL> <OPTIONS>
```

Just a quick explanation of the options to docker run:

`-t`: Means emulate a pseudo-TTY in the container. Obit needs a TTY for log messages.

`--rm`: Means delete the container after the end of the run (otherwise you end up with loads of dangling containers in your system)

`-v ${PWD}:/scratch` : Means mount the current directory inside the container as /scratch. This is where the script will run in the container.

`-e LOCAL_USER_ID=$(id -u) -e LOCAL_GROUP_ID=$(id -g)` : Means make your current user id and user group inside the container so that output files are written as your user rather than 'root'.

Installing Obit binary distribution on the SARAO comXX machines
---------------------------------------------------------------

Please note that binary Obit requires a Python > 3.6 and <= 3.8 to run. If you require a newer Python, you will have to install Obit from Source.

1. Make a working directory for your Obit installation and work from there eg:
```
	$ mkdir ~/Obit
	$ cd ~/Obit
```

2. Download and untar the desired binary package (r648 in the example below). Available Obit binary distributions can be found at https://www.cv.nrao.edu/~bcotton/ObitBin/linux_distro/:
```
	$ export OBIT_URL=https://www.cv.nrao.edu/~bcotton/ObitBin/linux_distro/
	$ curl ${OBIT_URL}/Obit.AVX-1.1.648.tar.gz | tar xzf -
```

3. Download KATObitPipe:
```
	$ git clone https://github.com/mauch/KATObitPipe
```

4. Copy the setup.sh Obit script to the base Obit install dir (`~/Obit` in this example):
```
	$ cp KATObitPipe/setup.sh ~/Obit
```

5. Modify the first line of setup.sh so that OBIT_ROOT points to the dir you have put the Obit distro - In this example:
```
	$ export OBIT_ROOT=/home/tmauch/Obit/obit-distro-1.1.648
```

6. Source the Obit setup script:
```
	$ source setup.sh
```

7. Copy the static metadata from KATObitPipe distro into Obit, unzip the models, and chmod them to 777 (otherwise Obit can’t see them):
```
	$ cp ./KATObitPipe/FITS/* ${OBIT_ROOT}/share/obit/data
	$ gunzip ${OBIT_ROOT}/share/obit/data/*.fits.gz
	$ chmod -R 777 ${OBIT_ROOT}/share/obit/data/*.fits
```

8. Copy the Obit static metadata dir to the place KATObitPipe sees it by default:
```
	$ mkdir ${OBIT}/share
	$ mv ${OBIT_ROOT}/share/obit/data ${OBIT}/share
```

9. Remove the `libgfortran` from the Obit binaries and use the system installed one instead.

   The version that ships with binary Obit distributions causes the AIPS tasks downloaded on-the-fly by KATObitPipe to crash. It is best to use the system `libgfortran` instead. This assumes `libgfortran` is installed on your system. In Debian based Linux use `apt install libgfortran` or something similar.
```
	$ rm -rf ${OBIT_ROOT}/lib/libgfortran*
```

10. Install KATObitPipe:
```
	$ pip install ./KATObitPipe
```
