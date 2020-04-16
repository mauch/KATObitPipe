import sys, pydoc
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil
from AIPS import AIPSDisk
from FITS import FITSDisk
from PipeUtil import *
from .KATCal import *
import h5py
import dask.array as da
from . import MakeIFs
from . import FITS2jpeg
import katpoint
import katdal as katfile
import subprocess
from . import KATH5toAIPS
import os
from . import AIPSSetup
import shutil
import numpy as np
from .KATImExceptions import KATUnimageableError

#possible kwargs: scratchdir
def MKContPipeline(files, outputdir, **kwargs):
    """MeerKAT Continuum pipeline.

    Parameters
    ----------
    files : list
        h5 filenames (note: support for multiple h5 files 
        i.e. ConcatenatedDataSet is not currently supported)
    outputdir : string
        Directory location to write output data, 
    scratchdir : string, optional
        The directory location of the aips disk
    parmFile : string, optional
        Overwrite the default imaging parameters using this parameter file.
    """
    if len(files) == 1:
        h5file = files[0]
    else:
        h5file = files
    ############### Initialize katfile object #########################
    OK = False
    # Open the h5 file as a katfile object
    try:
        #open katfile and perform selection according to kwargs
        katdal_ref_ant = kwargs.get('katdal_refant', '')
        katdata = katfile.open(h5file, ref_ant=katdal_ref_ant)
        OK = True
    except Exception as exception:
        print(exception)
    if not OK:
        raise KATUnimageableError("Unable to read MVF data in " + str(h5file))

    # If we are doing polcal- search for the most recent delaycal observation
    if kwargs.get('polcal'):
        delay_katdata = KATGetDelayCal(katdata, h5file)
        kwargs["delay_katdata"] = delay_katdata   

    # Die gracefully if we cannot write to the output area...
    if not os.path.exists(outputdir):
        print('Specified output directory: '+ outputdir + 'does not exist.')
        exit(-1)

    # Obit error logging
    err = OErr.OErr()

    #################### Initialize filenames #######################################################
    fileRoot      = os.path.join(outputdir,
                                 katdata.obs_params.get('capture_block_id', 
                                                        katdata.experiment_id)) # root of file name
    logFile       = fileRoot + ".log"   # Processing log file
    avgClass      = ("UVAv")[0:6]  # Averaged data AIPS class
    manifestfile  = outputdir + '/manifest.pickle'

    ############################# Initialize OBIT and AIPS ##########################################
    noScrat     = []
    # Logging directly to logFile
    OErr.PInit(err, 2, logFile)
    EVLAAddOutFile(os.path.basename(logFile), 'project', 'Pipeline log file')
    if kwargs.get('reuse'):
        ObitSys = AIPSSetup.AIPSSetup(err,configfile=kwargs.get('configFile'),scratchdir=kwargs.get('scratchdir'),aipsdisk=kwargs.get('aipsdisk'),overwrite=False)
    else:
        ObitSys = AIPSSetup.AIPSSetup(err,configfile=kwargs.get('configFile'),scratchdir=kwargs.get('scratchdir'),aipsdisk=kwargs.get('aipsdisk'))

    # Get the set up AIPS environment.
    AIPS_ROOT    = os.environ['AIPS_ROOT']
    AIPS_VERSION = os.environ['AIPS_VERSION']

    nThreads = 72
    user = OSystem.PGetAIPSuser()
    AIPS.userno = user
    disk = 1
    fitsdisk = 1
    nam = katdata.obs_params.get('capture_block_id', katdata.experiment_id)[:10]
    clss = "Raw"
    seq = 1

    ############### Condition data #########################

    #Get calibrator models
    fluxcals = katpoint.Catalogue(open(FITSDir.FITSdisks[1]+"/PERLEY_BUTLER_2013.csv"))
    #Condition data (get bpcals, update names for aips conventions etc)
    KATh5Condition(katdata, fluxcals, err)

    ############################# Initialise Parameters ##########################################
    ####### Initialize parameters dictionary ##### 
    parms = KATInitContParms()
    parms['PolCal'] = kwargs.get('polcal')
    parms['XYtarg'] = kwargs.get('XYtarg')
    # Get default XYtarg if it is not set
    targs = [targ.name for targ in katdata.catalogue.targets]
    if parms['XYtarg'] is None:
        GOTTARG = False
        for targ in ['1934-638', '0408-65']:
            if targ in targs:
                parms['XYtarg'] = targ
                GOTTARG = True
                break
        if not GOTTARG:
            raise RuntimeError('No default targets (1934-638, 0408-65) for XYFix. Cannot run in PolCal mode.')
    else:
        if parms['XYtarg'] not in targs:
            raise RuntimeError('XYtarg target %s not in observation. Cannot run in PolCal mode.' % (parms['XYtarg']))
    ####### User defined parameters ######
    if kwargs.get('parmFile'):
        print("parmFile",kwargs.get('parmFile'))
        exec(open(kwargs.get('parmFile')).read())
        EVLAAddOutFile(os.path.basename(kwargs.get('parmFile')), 'project', 'Pipeline input parameters')

    ###################### Data selection and static edits ############################################
    # Select data based on static imageable parameters
    KATh5Select(katdata, parms, err, **kwargs)

    # General AIPS data parameters at script level
    dataClass = ("UVDa")[0:6]      # AIPS class of raw uv data
    delayClass = "DELA"
    band      = katdata.spectral_windows[0].product #Correlator product
    project   = os.path.basename(os.path.splitext(files[0])[0])[0:10]  # Project name (12 char or less, used as AIPS Name)
    outIClass = parms["outIClass"] # image AIPS class
    debug     = parms["debug"]
    check     = parms["check"]

    ####################### Import data into AIPS #####################################################
    # Reuse or nay?
    sw = katdata.spectral_windows[katdata.spw]
    # Pick up static flags
    if sw.band == 'L':
        sflags = FetchObject(ObitTalkUtil.FITSDir.FITSdisks[fitsdisk] + 'maskred.pickle')
        if kwargs.get('flag', None):
            mess = 'Using static RFI mask in file %s for L-band' % (ObitTalkUtil.FITSDir.FITSdisks[fitsdisk] + 'maskred.pickle',)
            printMess(mess, logFile)
    elif sw.band == 'UHF':
        sflags = FetchObject(ObitTalkUtil.FITSDir.FITSdisks[fitsdisk] + 'maskredUHF.pickle')
        if kwargs.get('flag', None):
            mess = 'Using static RFI mask in file %s for UHF-band' % (ObitTalkUtil.FITSDir.FITSdisks[fitsdisk] + 'maskredUHF.pickle',)
            printMess(mess, logFile)
    else:
        sflags = np.zeros(sw.num_chans, dtype=np.bool)
    sflags = sflags[katdata.channels]
    # Construct a template uvfits file from master template
    mastertemplate = ObitTalkUtil.FITSDir.FITSdisks[fitsdisk] + 'MKATTemplate.uvtab.gz'
    outtemplate = nam + '.uvtemp'
    if kwargs.get('reuse'):
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project), dataClass, disk, seq, True, err)
        obsdata = KATH5toAIPS.GetKATMeta(katdata, err)
        # Extract AIPS parameters of the uv data to the metadata
        obsdata["Aproject"] = uv.Aname
        obsdata["Aclass"] = uv.Aclass
        obsdata["Aseq"] = uv.Aseq
        obsdata["Adisk"] = disk
        obsdata["calInt"] = katdata.dump_period
        obsdata["fitsdisk"] = fitsdisk
        # TODO: Check if the input data has been Hanned.
        doneHann = True
    else:
        mess = '\nLoading UV data with CBID: %s' % (katdata.obs_params['capture_block_id'],)
        printMess(mess, logFile)
        KATH5toAIPS.MakeTemplate(mastertemplate, outtemplate, katdata)
        uv = OTObit.uvlod(outtemplate, 0, EVLAAIPSName(project), clss, disk, seq, err)
        obsdata = KATH5toAIPS.KAT2AIPS(katdata, uv, disk, fitsdisk, err, calInt=katdata.dump_period, static=sflags, **kwargs)
        MakeIFs.UVMakeIF(uv,8,err,solInt=katdata.dump_period)
        os.remove(outtemplate)

    if parms["PolCal"]:
        mess = '\nLoading delay calibration with CBID: %s' % (delay_katdata.obs_params['capture_block_id'],)
        printMess(mess, logFile)
        # Load the delay cal observation
        KATH5toAIPS.MakeTemplate(mastertemplate, outtemplate, katdata)
        delay_uv = OTObit.uvlod(outtemplate, 0, EVLAAIPSName(project), delayClass, disk, seq, err)
        KATH5toAIPS.KAT2AIPS(delay_katdata, delay_uv, disk, fitsdisk, err, calInt=katdata.dump_period, static=sflags, flag=False)
        MakeIFs.UVMakeIF(delay_uv, 8, err, solInt=katdata.dump_period)     
        os.remove(outtemplate)
    
    # Print the uv data header to screen.
    uv.Header(err)
    ############################# Set Project Processing parameters ###################################
    # Parameters derived from obsdata and katdata
    KATGetObsParms(obsdata, katdata, parms, logFile)

    ###### Initialise target parameters #####
    KATInitTargParms(katdata, parms, err)

    # Load the outputs pickle jar
    EVLAFetchOutFiles()

    OSystem.PAllowThreads(nThreads)   # Allow threads in Obit/oython
    retCode = 0
    doBand = -1
    BPVer = 0
    maxgap = max(parms["CalAvgTime"], 160.*katdata.dump_period)/60.
    ################### Start processing ###############################################################

    mess = "Start project "+parms["project"]+" AIPS user no. "+str(AIPS.userno)+\
           ", KAT7 configuration "+parms["KAT7Cfg"]
    printMess(mess, logFile)
    if debug:
        pydoc.ttypager = pydoc.plainpager # don't page task input displays
        mess = "Using Debug mode "
        printMess(mess, logFile)
    if check:
        mess = "Only checking script"
        printMess(mess, logFile)

    # Log parameters
    printMess("Parameter settings", logFile)
    for p in parms:
        mess = "  "+p+": "+str(parms[p])
        printMess(mess, logFile)
    clist = []
    for DCal in parms["DCals"]:
        if DCal["Source"] not in clist:
            clist.append(DCal["Source"])
    for PCal in parms["PCals"]:
        if PCal["Source"] not in clist:
            clist.append(PCal["Source"])
    for ACal in parms["ACals"]:
        if ACal["Source"] not in clist:
            clist.append(ACal["Source"])
    if kwargs.get('targets') is not None:
        targets = [targ.name for targ in katdata.catalogue if (targ.name not in clist) and (targ.name in kwargs.get('targets').split(','))]
    else:
        targets = [targ.name for targ in katdata.catalogue if (targ.name not in clist)]
    
    refAnt = kwargs.get('refant')
    if refAnt is not None:
        try:
            SaveObject(obsdata['antLookup'][refAnt], fileRoot+".refAnt.pickle", True)
        except:
            mess = "Select reference antenna " + refAnt + " not in antenna table."
            printMess(mess, logFile)
            print(mess)
    refAnt = FetchObject(fileRoot+".refAnt.pickle")

    # Save parameters to pickle jar, manifest
    ParmsPicklefile = fileRoot+".Parms.pickle"   # Where results saved
    SaveObject(parms, ParmsPicklefile, True)
    EVLAAddOutFile(os.path.basename(ParmsPicklefile), 'project', 'Processing parameters used' )
    loadClass = dataClass

    # Hanning - only if not reusing
    doneHann = False
    if not kwargs.get('reuse'):
        if parms["doHann"]:
            uv = KATHann(uv, EVLAAIPSName(project), dataClass, disk, seq, err, \
                      doDescm=parms["doDescm"], flagVer=-1, logfile=logFile, zapin=True, check=check, debug=debug)
            doneHann = True
    
    if parms["PolCal"] and parms["doHann"]:
        mess = "Hanning delay calibration scan"
        printMess(mess, logFile)
        delay_uv = KATHann(delay_uv, EVLAAIPSName(project), delayClass, disk, seq + 1, err, \
                        doDescm=parms["doDescm"], flagVer=-1, logfile=logFile, zapin=True, check=check, debug=debug)

    if doneHann:
        # Halve channels after hanning.
        parms["selChan"]=int(parms["selChan"]/2)
        parms["BChDrop"]=int(parms["BChDrop"]/2)
        parms["EChDrop"]=int(parms["EChDrop"]/2)
        if uv==None and not check:
            raise RuntimeError("Cannot Hann data ")

    # Clear any old calibration/editing 
    if parms["doClearTab"] or kwargs.get('reuse'):
        mess =  "Clear previous calibration"
        printMess(mess, logFile)
        EVLAClearCal(uv, err, doGain=parms["doClearGain"], doFlag=parms["doClearFlag"], doBP=parms["doClearBP"], check=check)
        OErr.printErrMsg(err, "Error resetting calibration")

    # Copy FG 1 to FG 2
    if parms["doCopyFG"]:
        mess =  "Copy FG 1 to FG 2"
        printMess(mess, logFile)
        retCode = KATCopyFG(uv, err, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error Copying FG table")
  
    # Flag antennas shadowed by others?
    if parms["doShad"]:
        retCode = EVLAShadow (uv, err, shadBl=parms["shadBl"], \
                              logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error Shadow flagging data")
    
    # Median window time editing, for RFI impulsive in time
    if parms["doMednTD1"]:
        mess =  "Median window time editing, for RFI impulsive in time:"
        printMess(mess, logFile)
        retCode = EVLAMedianFlag (uv, clist, err, noScrat=noScrat, nThreads=nThreads, \
                                  avgTime=parms["mednAvgTime"], avgFreq=parms["mednAvgFreq"],  chAvg= parms["mednChAvg"], \
                                  timeWind=parms["mednTimeWind"],flagVer=2, flagTab=2,flagSig=parms["mednSigma"], \
                                  logfile=logFile, check=check, debug=False)
        if retCode!=0:
            raise RuntimeError("Error in MednFlag")
    
    # Median window frequency editing, for RFI impulsive in frequency
    if parms["doFD1"]:
        mess =  "Median window frequency editing, for RFI impulsive in frequency:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, clist, err, flagVer=2, flagTab=2, doCalib=-1, doBand=-1,   \
                                timeAvg=parms["FD1TimeAvg"], \
                                doFD=True, FDmaxAmp=1.0e20, FDmaxV=1.0e20, FDwidMW=parms["FD1widMW"],  \
                                FDmaxRMS=[1.0e20,0.1], FDmaxRes=parms["FD1maxRes"],  \
                                FDmaxResBL= parms["FD1maxRes"],  FDbaseSel=parms["FD1baseSel"],\
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError("Error in AutoFlag")
    
    # Parallactic angle correction?
    if parms["doPACor"]:
        retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error in Parallactic angle correction")
    
    # Need to find a reference antenna?  See if we have saved it?
    if (parms["refAnt"]<=0):
        refAnt = FetchObject(fileRoot+".refAnt.pickle")
        if refAnt:
            parms["refAnt"] = refAnt

    # Use bandpass calibrator and center half of each spectrum
    if parms["refAnt"]<=0:
        mess = "Find best reference antenna: run Calib on BP Cal(s) "
        printMess(mess, logFile)
        parms["refAnt"] = EVLAGetRefAnt(uv, parms["BPCals"], err, flagVer=0, \
                                        solInt=parms["bpsolint1"], nThreads=nThreads, \
                                        logfile=logFile, check=check, debug=debug)
        if err.isErr:
                raise  RuntimeError("Error finding reference antenna")
        if parms["refAnts"][0]<=0:
            parms["refAnts"][0] = parms["refAnt"]
        mess = "Picked reference antenna "+str(parms["refAnt"])
        printMess(mess, logFile)
        # Save it
        ParmsPicklefile = fileRoot+".Parms.pickle"   # Where results saved
        SaveObject(parms, ParmsPicklefile, True)
        refAntPicklefile = fileRoot+".refAnt.pickle"   # Where results saved
        SaveObject(parms["refAnt"], refAntPicklefile, True)


    # Plot Raw, edited data?
    if parms["doRawSpecPlot"] and parms["plotSource"]:
        mess =  "Raw Spectral plot for: "+' '.join(parms["BPCal"])
        printMess(mess, logFile)
        plotFile = fileRoot+"_RawSpec.ps"
        retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                               Stokes=["RR","LL"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError("Error in Plotting spectrum")
        EVLAAddOutFile(plotFile, 'project', 'Pipeline log file' )
    
    if parms["PolCal"]:
        mess = "XYphase bandpass calibration"
        printMess(mess, logFile)
        retCode = KATXPhase(delay_uv, uv, err, logfile=logFile, check=check, debug=debug,
                            doCalib=-1, flagVer=0, doBand=-1, refAnt=parms['refAnt'])
        doBand = 1
        BPVer += 1
        if retCode!=0:
            raise RuntimeError("Error in Xphase calibration")

    # delay calibration
    if parms["doDelayCal"] and parms["DCals"] and not check:
        plotFile = fileRoot+"_DelayCal.ps"
        retCode = EVLADelayCal(uv, parms["DCals"], err,  \
                               BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                               doCalib=-1, flagVer=0, doBand=doBand, BPVer=BPVer, \
                               solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                               refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], 
                               doZeroPhs=parms["delayZeroPhs"], \
                               doAvgIF=parms["delayAvgIF"], doAvgPol=parms["delayAvgPol"], \
                               doPlot=parms["doSNPlot"], plotFile=plotFile, \
                               nThreads=nThreads, noScrat=noScrat, \
                               logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error in delay calibration")

        # Plot corrected data?
        if parms["doSpecPlot"] and parms["plotSource"]:
            plotFile = fileRoot+"_DelaySpec.ps"
            retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, \
                                   plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=doBand,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError("Error in Plotting spectrum")

    # Bandpass calibration
    if parms["doBPCal"] and parms["BPCals"]:
        retCode = KATBPCal(uv, parms["BPCals"], err, doBand=doBand, BPVer=BPVer, newBPVer=0,
                            noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error in Bandpass calibration")

        # Plot corrected data?
        if parms["doSpecPlot"] and  parms["plotSource"]:
            plotFile = fileRoot+"_BPSpec.ps"
            retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, \
                                   parms["refAnt"], err, Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError("Error in Plotting spectrum")
    

    # Amp & phase Calibrate
    if parms["doAmpPhaseCal"]:
        plotFile = fileRoot+"_APCal.ps"
        retCode = KATCalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], 
                             doCalib=2, doBand=1, BPVer=0, flagVer=0, \
                             BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                             solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                             doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                             ampEditFG=parms["ampEditFG"], avgPol=parms["PolCal"], \
                             doPlot=parms["doSNPlot"], plotFile=plotFile,  refAnt=parms["refAnt"], \
                             nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)

        if retCode!=0:
            raise RuntimeError("Error calibrating")

    # More editing
    if parms["doAutoFlag"]:
        mess =  "Post calibration editing:"
        printMess(mess, logFile)
        # if going to redo then only calibrators
        if parms["doRecal"]:
            # Only calibrators
            clist = []
            for DCal in parms["DCals"]:
                if DCal["Source"] not in clist:
                    clist.append(DCal["Source"])
            for PCal in parms["PCals"]:
                if PCal["Source"] not in clist:
                    clist.append(PCal["Source"])
            for ACal in parms["ACals"]:
                if ACal["Source"] not in clist:
                    clist.append(ACal["Source"])
        else:
            clist = []

        retCode = EVLAAutoFlag (uv, clist, err, flagVer=0, flagTab =2, \
                                doCalib=2, gainUse=0, doBand=1, BPVer=BPVer,  \
                                IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                doFD=parms["doFirstAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                FDmaxRes=parms["FDmaxRes"],  FDmaxResBL=parms["FDmaxResBL"], \
                                FDbaseSel=parms["FDbaseSel"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError("Error in AutoFlag")

    # Redo the calibration using new flagging?
    if parms["doBPCal2"]==None:
        parms["doBPCal2"] = parms["doBPCal"]
    if parms["doDelayCal2"]==None:
        parms["doDelayCal2"] = parms["doDelayCal2"]
    if parms["doAmpPhaseCal2"]==None:
        parms["doAmpPhaseCal2"] = parms["doAmpPhaseCal"]
    if parms["doAutoFlag2"]==None:
        parms["doAutoFlagCal2"] = parms["doAutoFlag"]
    if parms["doRecal"]:
        mess =  "Redo calibration:"
        printMess(mess, logFile)
        EVLAClearCal(uv, err, doGain=True, doFlag=False, doBP=True, check=check, logfile=logFile)
        OErr.printErrMsg(err, "Error resetting calibration")
        BPVer = 0
        # Parallactic angle correction?
        if parms["doPACor"]:
            retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError("Error in Parallactic angle correction")


        # Run MKXPhase on delaycal data and attach BP table to UV data
        if parms["PolCal"]:
            mess = "XYphase bandpass calibration"
            printMess(mess, logFile)
            retCode = KATXPhase(delay_uv, uv, err, logfile=logFile, check=check, debug=debug,
                            doCalib=-1, flagVer=0, doBand=-1, refAnt=parms['refAnt'])
            BPVer += 1
        if retCode!=0:
            raise RuntimeError("Error in Xphase calibration")


        # Delay recalibration
        if parms["doDelayCal2"] and parms["DCals"] and not check:
            plotFile = fileRoot+"_DelayCal2.ps"
            retCode = EVLADelayCal(uv, parms["DCals"], err, \
                                   BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                                   doCalib=-1, flagVer=0, doBand=doBand, BPVer=BPVer, \
                                   solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                                   refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], \
                                   doZeroPhs=parms["delayZeroPhs"], \
                                   doAvgIF=parms["delayAvgIF"], doAvgPol=parms["delayAvgPol"], \
                                   doPlot=parms["doSNPlot"], plotFile=plotFile, \
                                   nThreads=nThreads, noScrat=noScrat, \
                                   logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError("Error in delay calibration")

            # Plot corrected data?
            if parms["doSpecPlot"] and parms["plotSource"]:
                plotFile = fileRoot+"_DelaySpec2.ps"
                retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                       Stokes=["RR","LL"], doband=doband,          \
                                       check=check, debug=debug, logfile=logFile )
                if retCode!=0:
                    raise  RuntimeError("Error in Plotting spectrum")

        # Bandpass calibration
        if parms["doBPCal2"] and parms["BPCals"]:
            retCode = KATBPCal(uv, parms["BPCals"], err, doBand=doBand, BPVer=BPVer, newBPVer=0, \
                            noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError("Error in Bandpass calibration")
        
            # Plot corrected data?
            if parms["doSpecPlot"] and parms["plotSource"]:
                plotFile = fileRoot+"_BPSpec2.ps"
                retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError("Error in Plotting spectrum")


        # Amp & phase Recalibrate
        if parms["doAmpPhaseCal2"]:
            plotFile = fileRoot+"_APCal2.ps"
            retCode = KATCalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], \
                                 doCalib=2, doBand=1, BPVer=0, flagVer=0, \
                                 BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                                 solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                                 doAmpEdit=True, ampSigma=parms["ampSigma"], \
                                 ampEditFG=parms["ampEditFG"], avgPol=parms["PolCal"], \
                                 doPlot=parms["doSNPlot"], plotFile=plotFile, refAnt=parms["refAnt"], \
                                 noScrat=noScrat, nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError("Error calibrating")

        # More editing
        if parms["doAutoFlag2"]:
            mess =  "Post recalibration editing:"
            printMess(mess, logFile)
            retCode = EVLAAutoFlag (uv, [], err, flagVer=0, flagTab=2, \
                                    doCalib=2, gainUse=0, doBand=1, BPVer=0,  \
                                    IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                    doFD=parms["doSecAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                    FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                    FDmaxRes=parms["FDmaxRes"],  FDmaxResBL= parms["FDmaxResBL"], \
                                    FDbaseSel=parms["FDbaseSel"], \
                                    nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise  RuntimeError("Error in AutoFlag")
    # end recal
    # Calibrate and average data
    # Overwrite avgStokes from command line
    if kwargs.get('halfstokes'):
        parms["avgStokes"] = 'HALF'
    if parms["doCalAvg"] == 'Splat':
        retCode = KATCalAvg (uv, avgClass, parms["seq"], parms["CalAvgTime"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=0, doPol=False, \
                              avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], Stokes=parms["avgStokes"], \
                              BChan=1, EChan=parms["selChan"] - 1, doAuto=parms["doAuto"], \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError("Error in CalAvg")
    elif parms["doCalAvg"] == 'BL':
        retCode = KATBLCalAvg (uv, avgClass, parms["seq"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=0, doPol=False, \
                              avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], FOV=parms['FOV'], \
                              maxInt=min(parms["solPInt"],parms["solAInt"]), Stokes=parms["avgStokes"], \
                              BChan=1, EChan=parms["selChan"] - 1, timeAvg=parms["CalAvgTime"], \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError("Error in BLCalAvg")

    if parms["doSaveTab"]:
        filename = project+".CalTab.uvtab"
        _ = EVLAUVFITSTab (uv, filename, 0, err, logfile=logFile)

    #Zap unaveraged data if requested
    if kwargs.get('zapraw'):
        uv.Zap(err)
    
    # Get calibrated/averaged data
    if not check:
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project), avgClass[0:6], \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")

    plotFile = fileRoot+"_Spec.ps"
    retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, \
                               plotFile, parms["refAnt"], err, \
                               Stokes=["I"], doband=-1, docalib=-1,      \
                               check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError("Error in Plotting spectrum")

    # KATUVFITS(uv, 'preimage.uvfits', 0, err, exclude=["AIPS HI", "AIPS SL", "AIPS PL"], 
    # include=["AIPS AN", "AIPS FQ"], compress=parms["Compress"], logfile=logFile)
    KATUVFITab(uv, project+'.uvtab', 0, err)
    #Gzip the data?
    if kwargs.get('gzip'):
        os.system('pigz -p %d %s'%(nThreads, project+'.uvtab'))
        os.system('rm -f %s'%(project+'.uvtab'))

class DataProductError(Exception):
    """ Exception for data product (output file) errors. """
    pass

class TooManyKatfilesException(Exception):
    """ Exception in KATPipe. """
    pass


