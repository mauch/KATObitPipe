import sys, pydoc
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil
from AIPS import AIPSDisk
from FITS import FITSDisk
from PipeUtil import *
from KATCal import *
import h5py
import dask.array as da
import MakeIFs
import FITS2jpeg
import katpoint
import katdal as katfile
import subprocess
import KATH5toAIPS
import os
import AIPSSetup
import shutil
from KATImExceptions import KATUnimageableError

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
    #if len(files) > 1:
    #    raise TooManyKatfilesException('Processing multiple katfiles are not currently supported')
    h5file = files

    # Die gracefully if we cannot write to the output area...
    if not os.path.exists(outputdir):
        print 'Specified output directory: '+ outputdir + 'does not exist.'
        exit(-1)

    # Obit error logging
    err = OErr.OErr()

    #################### Initialize filenames #######################################################
    fileRoot      = os.path.join(outputdir, os.path.basename(os.path.splitext(files[0])[0])) # root of file name
    logFile       = fileRoot+".log"   # Processing log file
    avgClass      = ("UVAv")[0:6]  # Averaged data AIPS class
    manifestfile  = outputdir + '/manifest.pickle'

    ############################# Initialize OBIT and AIPS ##########################################
    noScrat     = []
    # Logging directly to logFile
    OErr.PInit(err, 2, logFile)
    EVLAAddOutFile(os.path.basename(logFile), 'project', 'Pipeline log file')

    ObitSys = AIPSSetup.AIPSSetup(err,configfile=kwargs.get('configFile'),scratchdir=kwargs.get('scratchdir'))

    # Get the set up AIPS environment.
    AIPS_ROOT    = os.environ['AIPS_ROOT']
    AIPS_VERSION = os.environ['AIPS_VERSION']

    nThreads = 24
    user = OSystem.PGetAIPSuser()
    AIPS.userno = user
    disk = 1
    fitsdisk = 0
    nam = os.path.basename(os.path.splitext(files[0])[0])[0:10]
    cls = "Raw"
    seq = 1

    ############################# Initialise Parameters ##########################################
    ####### Initialize parameters dictionary ##### 
    parms = KATInitContParms()
    ####### User defined parameters ######
    if kwargs.get('parmFile'):
        print "parmFile",kwargs.get('parmFile')
        exec(open(kwargs.get('parmFile')).read())
        EVLAAddOutFile(os.path.basename(kwargs.get('parmFile')), 'project', 'Pipeline input parameters' )

    ############### Initialize katfile object, uvfits object and condition data #########################
    OK = False
    # Open the h5 file as a katfile object
    try:
        #open katfile and perform selection according to kwargs
        katdata = katfile.open(h5file)
        OK = True
    except Exception, exception:
        print exception
    if not OK:
        OErr.PSet(err)
        OErr.PLog(err, OErr.Fatal, "Unable to read KAT HDF5 data in " + str(h5file))
        raise KATUnimageableError("Unable to read KAT HDF5 data in " + str(h5file))

    #We have a katdal object- read some flags and ad them in in available
    if kwargs.get('flags') is not None:
        flags=kwargs.get('flags')
        fa = flags.split(',')
        for fn,ff in enumerate(fa):
	        ex_flags_file = h5py.File(ff)
	        ex_flags = da.from_array(ex_flags_file['flags'], chunks=(1,342,katdata.shape[2]))
	        #Sum the new flags 
        	katdata.datasets[fn].source.data.flags = ex_flags

    #Are we MeerKAT or KAT-7
    telescope = katdata.ants[0].name[0]
    if telescope=='m':
        sefd=500.
    else:
        sefd=1200.
    #Get calibrator models
    fluxcals = katpoint.Catalogue(file(FITSDir.FITSdisks[0]+"/"+parms["fluxModel"]))
    #Condition data (get bpcals, update names for aips conventions etc)
    KATh5Condition(katdata,fluxcals,err)

    ###################### Data selection and static edits ############################################
    # Select data based on static imageable parameters
    MKATh5Select(katdata, parms, err, **kwargs)

    ####################### Import data into AIPS #####################################################
    # Construct a template uvfits file from master template
    mastertemplate=ObitTalkUtil.FITSDir.FITSdisks[fitsdisk]+'MKATTemplate.uvtab.gz'
    outtemplate=nam+'.uvtemp'
    KATH5toAIPS.MakeTemplate(mastertemplate,outtemplate,len(katdata.channel_freqs))
    uv=OTObit.uvlod(outtemplate,0,nam,cls,disk,seq,err)

    obsdata = KATH5toAIPS.KAT2AIPS(katdata, uv, disk, fitsdisk, err, calInt=1.0, **kwargs)
    MakeIFs.UVMakeIF(uv,8,err)

    # Print the uv data header to screen.
    uv.Header(err)
    os.remove(outtemplate)
    ############################# Set Project Processing parameters ###################################
    # Parameters derived from obsdata and katdata
    MKATGetObsParms(obsdata, katdata, parms, logFile)

    ###### Initialise target parameters #####
    KATInitTargParms(katdata,parms,err)

    # General AIPS data parameters at script level
    dataClass = ("UVDa")[0:6]      # AIPS class of raw uv data
    band      = katdata.spectral_windows[0].product #Correlator product
    project   = parms["project"][0:12]  # Project name (12 char or less, used as AIPS Name)
    outIClass = parms["outIClass"] # image AIPS class
    debug     = parms["debug"]
    check     = parms["check"]

    # Load the outputs pickle jar
    EVLAFetchOutFiles()

    OSystem.PAllowThreads(nThreads)   # Allow threads in Obit/oython
    retCode = 0

    maxgap = max(parms["CalAvgTime"],20*katdata.dump_period)/60.
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
    refAnt = FetchObject(fileRoot+".refAnt.pickle")

    # Save parameters to pickle jar, manifest
    ParmsPicklefile = fileRoot+".Parms.pickle"   # Where results saved
    SaveObject(parms, ParmsPicklefile, True)
    EVLAAddOutFile(os.path.basename(ParmsPicklefile), 'project', 'Processing parameters used' )
    loadClass = dataClass

    retCode = KATCalAvg (uv, "PREAVG", parms["seq"], parms["CalAvgTime"], err, \
                              flagVer=-1, doCalib=-1, gainUse=-1, doBand=-1, BPVer=-1, doPol=False, \
                              avgFreq=0, chAvg=1, BChan=1, EChan=0, doAuto=parms["doAuto"], Stokes=' ',\
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError,"Error in CalAvg"
    uv.Zap(err)
    # Get initially averaged data
    if not check:
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project), "PREAVG", \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating initial avg AIPS data")

    # Hanning
    parms["doHann"]=True
    if parms["doHann"]:
       # Set uv if not done
        if uv==None and not check:
            uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project), loadClass[0:6], disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating AIPS data")

        uv = KATHann(uv, EVLAAIPSName(project), dataClass, disk, parms["seq"], err, \
                      doDescm=parms["doDescm"], flagVer=0, logfile=logFile, check=check, debug=debug)
        #Halve channels after hanning.
        parms["selChan"]=int(parms["selChan"]/2)
        parms["BChDrop"]=int(parms["BChDrop"]/2)
        parms["EChDrop"]=int(parms["EChDrop"]/2)
        if uv==None and not check:
            raise RuntimeError,"Cannot Hann data "
 
    # Clear any old calibration/editing 
    if parms["doClearTab"]:
        mess =  "Clear previous calibration"
        printMess(mess, logFile)
        EVLAClearCal(uv, err, doGain=parms["doClearGain"], doFlag=parms["doClearFlag"], doBP=parms["doClearBP"], check=check)
        OErr.printErrMsg(err, "Error resetting calibration")

        # Quack to remove data from start and end of each scan
    if parms["doQuack"]:
        retCode = EVLAQuack (uv, err, begDrop=parms["quackBegDrop"], endDrop=parms["quackEndDrop"], \
                             Reason=parms["quackReason"], \
                             logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Quacking data"
    
    # Flag antennas shadowed by others?
    if parms["doShad"]:
        retCode = EVLAShadow (uv, err, shadBl=parms["shadBl"], \
                              logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Shadow flagging data"
    

    # Median window time editing, for RFI impulsive in time
    if parms["doMednTD1"]:
        mess =  "Median window time editing, for RFI impulsive in time:"
        printMess(mess, logFile)
        retCode = EVLAMedianFlag (uv, clist, err, noScrat=noScrat, nThreads=nThreads, \
                                  avgTime=parms["mednAvgTime"], avgFreq=parms["mednAvgFreq"],  chAvg= parms["mednChAvg"], \
                                  timeWind=parms["mednTimeWind"],flagVer=2, flagTab=2,flagSig=parms["mednSigma"], \
                                  logfile=logFile, check=check, debug=False)
        if retCode!=0:
            raise RuntimeError,"Error in MednFlag"
    
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
           raise  RuntimeError,"Error in AutoFlag"
    
    # Parallactic angle correction?
    if parms["doPACor"]:
        retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Parallactic angle correction"
    
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
                raise  RuntimeError,"Error finding reference antenna"
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
    parms["doRawSpecPlot"]=False
    parms["doSpecPlot"]=False
    if parms["doRawSpecPlot"] and parms["plotSource"]:
        mess =  "Raw Spectral plot for: "+' '.join(parms["BPCal"])
        printMess(mess, logFile)
        plotFile = fileRoot+"_RawSpec.ps"
        retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                               Stokes=["RR","LL"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
        EVLAAddOutFile(plotFile, 'project', 'Pipeline log file' )

    # delay calibration
    if parms["doDelayCal"] and parms["DCals"] and not check:
        plotFile = fileRoot+"_DelayCal.ps"
        retCode = EVLADelayCal(uv, parms["DCals"], err,  \
                               BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                               doCalib=2, flagVer=0, doBand=-1, \
                               solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                               refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], 
                               doZeroPhs=parms["delayZeroPhs"], \
                               doPlot=parms["doSNPlot"], plotFile=plotFile, \
                               nThreads=nThreads, noScrat=noScrat, \
                               logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in delay calibration"

        # Plot corrected data?
        if parms["doSpecPlot"] and parms["plotSource"]:
            plotFile = fileRoot+"_DelaySpec.ps"
            retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, \
                                   plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=-1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"
	print parms["bpBChan1"],parms["bpEChan1"],parms["bpBChan2"],parms["bpEChan2"],parms["bpChWid2"]
    # Bandpass calibration
    if parms["doBPCal"] and parms["BPCals"]:
        retCode = KATBPCal(uv, parms["BPCals"], err, noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Bandpass calibration"

        # Plot corrected data?
        if parms["doSpecPlot"] and  parms["plotSource"]:
            plotFile = fileRoot+"_BPSpec.ps"
            retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, \
                                   parms["refAnt"], err, Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

    # Amp & phase Calibrate
    if parms["doAmpPhaseCal"]:
        plotFile = fileRoot+"_APCal.ps"
        retCode = KATCalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], 
                             doCalib=2, doBand=1, BPVer=1, flagVer=0, \
                             BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                             solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                             doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                             ampEditFG=parms["ampEditFG"], \
                             doPlot=parms["doSNPlot"], plotFile=plotFile,  refAnt=parms["refAnt"], \
                             nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
        #print parms["ACals"],parms["PCals"]
        if retCode!=0:
            raise RuntimeError,"Error calibrating"

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
                                doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                doFD=parms["doAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                FDmaxRes=parms["FDmaxRes"],  FDmaxResBL=parms["FDmaxResBL"], \
                                FDbaseSel=parms["FDbaseSel"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"

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
        # Parallactic angle correction?
        if parms["doPACor"]:
            retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error in Parallactic angle correction"

        # Delay recalibration
        if parms["doDelayCal2"] and parms["DCals"] and not check:
            plotFile = fileRoot+"_DelayCal2.ps"
            retCode = EVLADelayCal(uv, parms["DCals"], err, \
                                   BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                                   doCalib=2, flagVer=0, doBand=-1, \
                                   solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                                   refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], \
                                   doZeroPhs=parms["delayZeroPhs"], \
                                   doPlot=parms["doSNPlot"], plotFile=plotFile, \
                                   nThreads=nThreads, noScrat=noScrat, \
                                   logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error in delay calibration"

            # Plot corrected data?
            if parms["doSpecPlot"] and parms["plotSource"]:
                plotFile = fileRoot+"_DelaySpec2.ps"
                retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                       Stokes=["RR","LL"], doband=-1,          \
                                       check=check, debug=debug, logfile=logFile )
                if retCode!=0:
                    raise  RuntimeError,"Error in Plotting spectrum"

        # Bandpass calibration
        if parms["doBPCal2"] and parms["BPCals"]:
            retCode = KATBPCal(uv, parms["BPCals"], err, noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error in Bandpass calibration"
        
            # Plot corrected data?
            if parms["doSpecPlot"] and parms["plotSource"]:
                plotFile = fileRoot+"_BPSpec2.ps"
                retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

        # Amp & phase Recalibrate
        if parms["doAmpPhaseCal2"]:
            plotFile = fileRoot+"_APCal2.ps"
            retCode = KATCalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], \
                                 doCalib=2, doBand=2, BPVer=1, flagVer=0, \
                                 BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                                 solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                                 doAmpEdit=True, ampSigma=parms["ampSigma"], \
                                 ampEditFG=parms["ampEditFG"], \
                                 doPlot=parms["doSNPlot"], plotFile=plotFile, refAnt=parms["refAnt"], \
                                 noScrat=noScrat, nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error calibrating"

        # More editing
        if parms["doAutoFlag2"]:
            mess =  "Post recalibration editing:"
            printMess(mess, logFile)
            retCode = EVLAAutoFlag (uv, [], err, flagVer=0, flagTab=2, \
                                    doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                    IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                    doFD=parms["doAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                    FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                    FDmaxRes=parms["FDmaxRes"],  FDmaxResBL= parms["FDmaxResBL"], \
                                    FDbaseSel=parms["FDbaseSel"], \
                                    nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise  RuntimeError,"Error in AutoFlag"
    # end recal
    # Calibrate and average data
    if parms["doCalAvg"]:
        retCode = KATCalAvg (uv, avgClass, parms["seq"], parms["CalAvgTime"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1, doPol=False, \
                              avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], \
                              BChan=1, EChan=parms["selChan"] - 1, doAuto=parms["doAuto"], \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in CalAvg"

    # Get calibrated/averaged data
    if not check:
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project), avgClass[0:6], \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")

    KATUVFITS(uv, 'preimage.uvfits', 0, err, exclude=["AIPS HI", "AIPS SL", "AIPS PL"], include=["AIPS AN", "AIPS FQ"], compress=parms["Compress"], logfile=logFile)
    KATUVFITab(uv, 'preimage.uvtab', 0, err)
    # XClip
    if parms["XClip"] and parms["XClip"]>0.0:
        mess =  "Cross Pol clipping:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, [], err, flagVer=-1, flagTab=1, \
                                doCalib=2, gainUse=0, doBand=-1, maxBad=1.0,  \
                                XClip=parms["XClip"], timeAvg=1./60., \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError,"Error in AutoFlag"
    
    # R-L  delay calibration cal if needed,
    if parms["doRLDelay"] and parms["RLDCal"][0][0]!=None:
        if parms["rlrefAnt"]<=0:
            parms["rlrefAnt"] =  parms["refAnt"]
        # parms["rlDoBand"] if before average, BPVer=parms["rlBPVer"], 
        retCode = EVLARLDelay(uv, err,\
                              RLDCal=parms["RLDCal"], BChan=parms["rlBChan"], \
                              EChan=parms["rlEChan"], UVRange=parms["rlUVRange"], \
                              soucode=parms["rlCalCode"], doCalib=parms["rlDoCal"], gainUse=parms["rlgainUse"], \
                              timerange=parms["rltimerange"], \
                              # NOT HERE doBand=parms["rlDoBand"], BPVer=parms["rlBPVer"],  \
                              flagVer=parms["rlflagVer"], \
                              refAnt=parms["rlrefAnt"], doPol=False,  \
                              nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                              check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in R-L delay calibration"
    
    # Polarization calibration
    if parms["doPolCal"]:
        if parms["PCRefAnt"]<=0:
            parms["PCRefAnt"] =  parms["refAnt"]
        retCode = EVLAPolCal(uv, parms["PCInsCals"], err, \
                             doCalib=2, gainUse=0, doBand=-1, flagVer=0, \
                             fixPoln=parms["PCFixPoln"], pmodel=parms["PCpmodel"], avgIF=parms["PCAvgIF"], \
                             solInt=parms["PCSolInt"], refAnt=parms["PCRefAnt"], solType=parms["PCSolType"], \
                             ChInc=parms["PCChInc"], ChWid=parms["PCChWid"], \
                             nThreads=nThreads, check=check, debug=debug, noScrat=noScrat, logfile=logFile)
        if retCode!=0 and (not check):
           raise  RuntimeError,"Error in polarization calibration: "+str(retCode)
        # end poln cal.
    
    
    # R-L phase calibration cal., creates new BP table
    if parms["doRLCal"] and parms["RLDCal"][0][0]!=None:
        plotFile = fileRoot+"_RLSpec2.ps"
        if parms["rlrefAnt"]<=0:
            parms["rlrefAnt"] =  parms["refAnt"]
        retCode = EVLARLCal(uv, err,\
                            RLDCal=parms["RLDCal"], BChan=parms["rlBChan"],
                            EChan=parms["rlEChan"], UVRange=parms["rlUVRange"], \
                            ChWid2=parms["rlChWid"], solInt1=parms["rlsolint1"], solInt2=parms["rlsolint2"], \
                            RLPCal=parms["RLPCal"], RLPhase=parms["RLPhase"], \
                            RM=parms["RLRM"], CleanRad=parms["rlCleanRad"], \
                            calcode=parms["rlCalCode"], doCalib=parms["rlDoCal"], gainUse=parms["rlgainUse"], \
                            timerange=parms["rltimerange"], FOV=parms["rlFOV"], \
                            doBand=-1, BPVer=1, flagVer=parms["rlflagVer"], \
                            refAnt=parms["rlrefAnt"], doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                            doPlot=parms["doSpecPlot"], plotFile=plotFile, \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                            check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in RL phase spectrum calibration"
    
    # VClip
    if parms["VClip"] and parms["VClip"]>0.0:
        mess =  "VPol clipping:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, [], err, flagVer=-1, flagTab=1, \
                                doCalib=2, gainUse=0, doBand=-1,  \
                                VClip=parms["VClip"], timeAvg=parms["timeAvg"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError,"Error in AutoFlag VClip"
    
    # Plot corrected data?
    parms["doSpecPlot"]=True
    if parms["doSpecPlot"] and parms["plotSource"]:
        plotFile = fileRoot+"_Spec.ps"
        retCode = EVLASpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, \
                               plotFile, parms["refAnt"], err, \
                               Stokes=["I"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
    
    
    # Image targets
    if parms["doImage"]:
        # If targets not specified, image all
        if len(parms["targets"])<=0:
            slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = targets
        slist=targets
        KATImageTargets (uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, OutlierArea=parms["outlierArea"],\
                          doCalib=-1, doBand=-1,  flagVer=-1, doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                          Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], Niter=parms["Niter"], \
                          CleanRad=parms["CleanRad"], minFlux=parms["minFlux"], OutlierSize=parms["OutlierSize"], \
                          xCells=parms["xCells"], yCells=parms["yCells"], Reuse=parms["Reuse"], minPatch=parms["minPatch"], \
                          maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], noNeg=parms["noNeg"], \
                          solPInt=parms["solPInt"], solPMode=parms["solPMode"], solPType=parms["solPType"], \
                          maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], nx=parms["nx"], ny=parms["ny"], \
                          solAInt=parms["solAInt"], solAMode=parms["solAMode"], solAType=parms["solAType"], \
                          avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR = parms["minSNR"], refAnt=parms["refAnt"], \
                          do3D=parms["do3D"], BLFact=parms["BLFact"], BLchAvg=parms["BLchAvg"], \
                          doMB=parms["doMB"], norder=parms["MBnorder"], maxFBW=parms["MBmaxFBW"], \
                          PBCor=parms["PBCor"],antSize=parms["antSize"], autoCen=parms["autoCen"], \
                          nTaper=parms["nTaper"], Tapers=parms["Tapers"], sefd=sefd, \
                          nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=False)
        # End image
    
    # Get report on sources
    if parms["doReport"]:
        # If targets not specified, do all
        if len(parms["targets"])<=0:
            slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = parms["targets"]
        Report = EVLAReportTargets(uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, \
                                       Stokes=parms["Stokes"], logfile=logFile, check=check, debug=debug)
        # Save to pickle jar
        ReportPicklefile = fileRoot+"_Report.pickle"   # Where results saved
        SaveObject(Report, ReportPicklefile, True) 
       
    # Write results, cleanup    
    # Save cal/average UV data? 
    if parms["doSaveUV"] and (not check):
        Aname = EVLAAIPSName(project)
        cno = AIPSDir.PTestCNO(disk, user, Aname, avgClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, avgClass, disk, parms["seq"], True, err)
            filename = fileRoot+"_Cal.uvtab"
            KATUVFITS (uv, filename, 0, err, exclude=["AIPS HI", "AIPS SL", "AIPS PL"], include=["AIPS AN", "AIPS FQ"], compress=parms["Compress"], logfile=logFile)
            EVLAAddOutFile(os.path.basename(filename), 'project', "Calibrated Averaged UV data" )
            # Save list of output files
            EVLASaveOutFiles(manifestfile)
            del uvt
    # Imaging results
    # If targets not specified, save all
    if len(parms["targets"])<=0:
        slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = parms["targets"]
    for target in slist:
        if parms["doSaveImg"] and (not check):
            for s in parms["Stokes"]:
                oclass = s+outIClass[1:]
                outname = target
                # Test if image exists
                cno = AIPSDir.PTestCNO(disk, user, outname, oclass, "MA", parms["seq"], err)
                #print cno
                if cno <= 0 :
                    continue
                x = Image.newPAImage("out", outname, oclass, disk, parms["seq"], True, err)
                outfilefits = fileRoot+'_'+target+"."+oclass+".fits"
                xf = KATImFITS(x, outfilefits, 0, err, logfile=logFile)
                x = Image.newPAImage("out", outname, oclass, disk, parms["seq"], True, err)
                outfile = fileRoot+'_'+target+"."+oclass+".fittab.fits"
                xf = EVLAImFITS (x, outfile, 0, err, logfile=logFile)
                EVLAAddOutFile(outfile, target, 'Image of '+ target)
                # Statistics
                zz=imstat(x, err, logfile=logFile)
                # Make a Jpeg image
                FITS2jpeg.fits2jpeg(outfilefits,chans=1,contrast=0.05,cmap='jet',area=0.7)
                EVLAAddOutFile(outfile.replace('.fits','.jpeg'), target, 'Jpeg image of '+ target)
    # end writing loop
    
    # Save list of output files
    EVLASaveOutFiles(manifestfile)
    OErr.printErrMsg(err, "Writing output")
    
    # Contour plots
    if parms["doKntrPlots"]:
        mess = "INFO --> Contour plots (doKntrPlots)"
        printMess(mess, logFile)
        EVLAKntrPlots( err, imName=parms["targets"], project=fileRoot,
                       disk=disk, debug=debug )
        # Save list of output files
        EVLASaveOutFiles(manifestfile)
    elif debug:
        mess = "Not creating contour plots ( doKntrPlots = "+str(parms["doKntrPlots"])+ " )"
        printMess(mess, logFile)

    # Source uv plane diagnostic plots
    if parms["doDiagPlots"]:
        mess = "INFO --> Diagnostic plots (doDiagPlots)"
        printMess(mess, logFile)
        # Get the highest number avgClass catalog file
        Aname = EVLAAIPSName( project )
        uvc = None
        if not check:
            uvname = project+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
        EVLADiagPlots( uvc, err, cleanUp=parms["doCleanup"], \
                           project=fileRoot, \
                           logfile=logFile, check=check, debug=debug )
        # Save list of output files
        EVLASaveOutFiles(manifestfile)
    elif debug:
        mess = "Not creating diagnostic plots ( doDiagPlots = "+str(parms["doDiagPlots"])+ " )"
        printMess(mess, logFile)
    
    # Save metadata
    srcMetadata = None
    projMetadata = None
    if parms["doMetadata"]:
        mess = "INFO --> Save metadata (doMetadata)"
        printMess(mess, logFile)
        uvc = None
        if not uvc:
            # Get calibrated/averaged data
            Aname = EVLAAIPSName(project)
            uvname = project+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    
        # Get source metadata; save to pickle file
        srcMetadata = EVLASrcMetadata( uvc, err, Sources=parms["targets"], seq=parms["seq"], \
                                       sclass=outIClass, Stokes=parms["Stokes"],\
                                       logfile=logFile, check=check, debug=debug )
        picklefile = fileRoot+".SrcReport.pickle" 
        SaveObject( srcMetadata, picklefile, True ) 
        EVLAAddOutFile(os.path.basename(picklefile), 'project', 'All source metadata' )
    
        # Get project metadata; save to pickle file
        projMetadata = KATProjMetadata( uvc, AIPS_VERSION, err, \
            PCals=parms["PCals"], ACals=parms["ACals"], \
            BPCals=parms["BPCals"], DCals=parms["DCals"], \
            project = project, band = band, \
            dataInUVF = parms["archRoot"], archFileID = fileRoot )
        picklefile = fileRoot+".ProjReport.pickle"
        SaveObject(projMetadata, picklefile, True) 
        EVLAAddOutFile(os.path.basename(picklefile), 'project', 'Project metadata' )
    else:
        # Fetch from pickle jar
         picklefile = fileRoot+".SrcReport.pickle"
         srcMetadata = FetchObject(picklefile)
         picklefile = fileRoot+".ProjReport.pickle"
         projMetadata = FetchObject(picklefile)
   
    # Write report
    if parms["doHTML"]:
        mess = "INFO --> Write HTML report (doHTML)"
        printMess(mess, logFile)
        KATHTMLReport( projMetadata, srcMetadata, \
                            outfile=fileRoot+"_report.html", \
                            logFile=logFile )
    
    # Write VOTable
    if parms["doVOTable"]:
        mess = "INFO --> Write VOTable (doVOTable)"
        printMess(mess, logFile)
        EVLAAddOutFile( 'VOTable.xml', 'project', 'VOTable report' ) 
        EVLAWriteVOTable( projMetadata, srcMetadata, filename=fileRoot+'_VOTable.xml' )
    
    # Save list of output files
    EVLASaveOutFiles(manifestfile)
    
    # Cleanup - delete AIPS files
    if parms["doCleanup"] and (not check):
        mess = "INFO --> Clean up (doCleanup)"
        printMess(mess, logFile)
        # Delete target images
        # How many Stokes images
        nstok = len(parms["Stokes"])
        for istok in range(0,nstok):
            oclass = parms["Stokes"][istok:istok+1]+outIClass[1:]
            AllDest(err, disk=disk,Aseq=parms["seq"],Aclass=oclass)
        
        # Delete initial UV data
        Aname = EVLAAIPSName(project)
        # Test if data exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, dataClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS RAW UV DATA", Aname, dataClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting raw AIPS data")
        # Zap calibrated/averaged data
        # Test if data exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, avgClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, avgClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting cal/avg AIPS data")
        # Zap UnHanned data if present
        loadClass = "Raw"
        # Test if image exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, loadClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, loadClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting cal/avg AIPS data")
        OErr.printErrMsg(err, "Writing output/cleanup")


    # Delete AIPS scratch DA00 and disk
    if os.path.exists(os.environ['DA00']): shutil.rmtree(os.environ['DA00'])
    for disk in ObitTalkUtil.AIPSDir.AIPSdisks:
        if os.path.exists(disk): shutil.rmtree(disk)


    # Shutdown
    mess = "Finished project "+parms["project"]+ \
    " AIPS user no. "+str(AIPS.userno)
    printMess(mess, logFile)
    OErr.printErr(err)
    OSystem.Shutdown(ObitSys)
    # end pipeline

class DataProductError(Exception):
    """ Exception for data product (output file) errors. """
    pass

class TooManyKatfilesException(Exception):
    """ Exception in KATPipe. """
    pass


