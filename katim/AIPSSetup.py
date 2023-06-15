# Define AIPS and FITS disks

import os,shutil
from . import AIPSLite
import configparser
import OErr, OSystem, AIPS, ObitTalkUtil

# Setup FITS, AIPS_VERSION and AIPS_DISK from a configuration file if it exists,
# If a specific config doesn't exist then use aips if it is installed, otherwise
# default AIPS_VERSION to cwd. Default AIPS_DISK is set to cwd and default FITS to
# $OBIT_EXEC/share/data.

def AIPSSetup(err,configfile=None,scratchdir=None,overwrite=True,aipsdisk='aipsdisk'):
    cwd          = './'
    home         = os.path.expanduser("~")
    # Is Obit installed and set up correctly? If not make data cwd/FITS.
    try:
        OBIT_EXEC    = os.environ['OBIT']
        OBIT_DATA    = OBIT_EXEC+'/share/data'
    except:
        OBIT_EXEC    = './'
        OBIT_DATA    = OBIT_EXEC+'/FITS'

    #AIPS user number (arbitrary)
    user = 100

    # Check if AIPS is defined and use it as the default if it is.
    try:
        aips_dir = os.environ['AIPS_ROOT']
        aips_version = os.environ['AIPS_VERSION'][-7:]
    except:
        aips_dir = './'
        aips_version = '31DEC22'    # Should sort out where to change this if necessary!!

    configdefaults   = {'aips_dir': aips_dir, 'obit_dir': OBIT_EXEC, 'aips_version': aips_version, 'scratch_area': cwd, 'metadata_dir': OBIT_DATA}
    config = configparser.ConfigParser(configdefaults)
    config.add_section('KATPIPE')

    if configfile:
        #Get the specified config file
        config.read(configfile)
    else:
        # Get the config file in the users home directory if it exists and overwrite the defaults
        config.read(home + '/.katimrc')

    if scratchdir:
        config.set('KATPIPE' ,'scratch_area', scratchdir)

    ############################# Initialize AIPS ##########################################

    # Sort out the AIPS environment variables for the defined configuration.
    AIPSLite.get_aips(basedir=config.get('KATPIPE','aips_dir'),version=config.get('KATPIPE','aips_version'))

    # Make the aips scratch disk
    DA00         = config.get('KATPIPE','scratch_area') + '/da00'
    AIPS_DISK    = config.get('KATPIPE','scratch_area') + aipsdisk

    if overwrite:
        #Overwrite any previous AIPS disks.
        if os.path.exists(DA00): shutil.rmtree(DA00)
        if os.path.exists(AIPS_DISK): shutil.rmtree(AIPS_DISK)

        # Create the aips disk and the AIPS environment for the disks.
        AIPSLite.make_disk(disk_path=AIPS_DISK)
        AIPSLite.filaip(force=True,data_dir=AIPS_DISK)
        AIPSLite.make_da00(da00_path=DA00)

    # Get the set up AIPS environment.
    AIPS_ROOT    = os.environ['AIPS_ROOT']
    AIPS_VERSION = os.environ['AIPS_VERSION']

    adirs = [(None, AIPS_DISK)]
    fdirs = [(None, config.get('KATPIPE','metadata_dir')), (None, './')]

    ############################# Initialize OBIT ##########################################
    ffdirs = []
    for f in fdirs: ffdirs.append(f[1])
    aadirs = []
    for a in adirs: aadirs.append(a[1])

    ObitSys = OSystem.OSystem ("Pipeline", 1, user, len(aadirs), aadirs, \
                                   len(ffdirs), ffdirs, True, False, err)
    OErr.printErrMsg(err, "Error with Obit startup")

    # setup environment
    ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                                OBIT_EXEC=config.get('KATPIPE','obit_dir'), DA00=DA00, ARCH=AIPSLite.arch(), \
                                aipsdirs=adirs, fitsdirs=fdirs)

    # List directories - removed for scripting
    ObitTalkUtil.ListAIPSDirs()
    ObitTalkUtil.ListFITSDirs()

    # Disks to avoid
    noScrat     = [0]          # AIPS disks to avoid

    return ObitSys
