# General data parameters
parms["check"]         = False      # Only check script, don't execute tasks
parms["debug"]         = False      # run tasks debug
parms["Compress"]      = False      # Use compressed UV data?
parms["doSwPwr"]       = False      # Make EVLA Switched power corr?
parms["calInt"]        = 1.0        # Calibration table interval (min)
parms["doLoadArchive"] = False       # Load AIPS data from archive?

# User Supplied data parameters
parms["BPCal"] = []
parms["ACal"]  = []
parms["PCal"]  = []
parms["targets"] = []
parms["polcal"] = []

# Archive parameters
parms["archRoot"]      = "NOT VLA DATA" # Root of ASDM/BDF data
parms["selBand"]       = "L"   # Selected band, def = first
parms["selConfig"]     = 1     # Selected frequency config, def = first
parms["selNIF"]        = 1     # Selected number of IFs, def = first

# Observation parameters
parms["selChan"]       = obsdata["numchan"]  # Selected number of channels, def = first
parms["project"]       = obsdata["Aproject"]      # Project name (12 char or less, used as AIPS Name)
parms["dataClass"]     = obsdata["Aclass"]        # AIPS class of raw uv data
parms["fluxModel"]     = "PERLEY_BUTLER_2013.csv" # Filename of flux calibrator model (in FITS)
parms["staticflags"]   = "KAT7_SRFI"              # Filename containing a list of frequencies to be flagged (in FITS)
parms["KAT7Freq"]      = obsdata["centerfreq"]    # Representive frequency
parms["KAT7Cfg"]       = obsdata["corrmode"]      # KAT-7 correlator configuraton

# Hanning
parms["doHann"]       = True        # Hanning needed for RFI?
parms["doDescm"]      = True        # Descimate Hanning?

# Parallactic angle correction
parms["doPACor"] =     False         # Make parallactic angle correction

# Special editing list
parms["doEditList"] =  True        # Edit using editList?
parms["editFG"] =      1            # Table to apply edit list to

# Editing
parms["doClearTab"]   = True        # Clear cal/edit tables
parms["doClearGain"]  = True        # Clear SN and CL tables >1
parms["doClearFlag"]  = True        # Clear FG tables > 1
parms["doClearBP"]    = True        # Clear BP tables?
parms["doCopyFG"]     = True        # Copy FG 1 to FG 2
parms["doQuack"]      = True        # Quack data?
parms["quackBegDrop"] = 0.1         # Time to drop from start of each scan in min
parms["quackEndDrop"] = 0.0         # Time to drop from end of each scan in min
parms["quackReason"]  = "Quack"     # Reason string
parms["doShad"]       = None        # Shadow flagging (config dependent)
parms["shadBl"]       = 25.0        # Minimum shadowing baseline (m)
parms["doElev"]       = False       # Do elevation flagging
parms["minElev"]      = 15.0        # Minimum elevation to keep.
parms["doFD1"]        = True        # Do initial frequency domain flagging
parms["FD1widMW"]     = 40          # Width of the initial FD median window
parms["FD1maxRes"]    = 5.0         # Clipping level in sigma
parms["FD1TimeAvg"]   = 1.0         # time averaging in min. for initial FD flagging
parms["FD1baseSel"]   = None        # Baseline fitting region for FD1 (updates by KAT7CorrParms)
parms["doMednTD1"]    = True        # Median editing in time domain?
parms["mednSigma"]    = 5.0         # Median sigma clipping level
parms["mednTimeWind"] = 1.0         # Median window width in min for median flagging
parms["mednAvgTime"]  = 0.0         # Median Averaging time in min
parms["mednAvgFreq"]  = 1           # Median 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
parms["mednChAvg"]    = 5           # Median flagger number of channels to average
parms["doRMSAvg"]    = True         # Edit calibrators by RMSAvg?
parms["RMSAvg"]      = 3.0          # AutoFlag Max RMS/Avg for time domain RMS filtering
parms["RMSTimeAvg"]  = 0.0          # AutoFlag time averaging in min.
parms["doAutoFlag"]  = True         # Autoflag editing after first pass calibration?
parms["doAutoFlag2"] = True         # Autoflag editing after final (2nd) calibration?
parms["IClip"]       = None         # AutoFlag Stokes I clipping
parms["VClip"]       = None
parms["XClip"]       = None         # AutoFlag cross-pol clipping
parms["timeAvg"]     = 0.33         # AutoFlag time averaging in min.
parms["doAFFD"]      = True         # do AutoFlag frequency domain flag
parms["FDwidMW"]     = 40           # Width of the median window
parms["FDmaxRMS"]    = [6.0,0.1]    # Channel RMS limits (Jy)
parms["FDmaxRes"]    = 5.0          # Max. residual flux in sigma
parms["FDmaxResBL"]  = 5.0          # Max. baseline residual
parms["FDbaseSel"]   = None         # Channels for baseline fit (Updated by KAT7CorrParms)
parms["FDmaxAmp"]    = None         # Maximum average amplitude (Jy)
parms["FDmaxV"]      = 2.0          # Maximum average VPol amp (Jy)
parms["minAmp"]      = 1.0e-5       # Minimum allowable amplitude
parms["BChDrop"]     = None         # number of channels to drop from start of each spectrum
# NB: based on original number of channels, halved for Hanning
parms["EChDrop"]     = None         # number of channels to drop from end of each spectrum
# NB: based on original number of channels, halved for Hanning

# Delay calibration
parms["doDelayCal"]   =  False      # Determine/apply delays from contCals
parms["delaySolInt"]  =  0.5        # delay solution interval (min)
parms["delaySmoo"]    =  0.25       # Delay smoothing time (hr)
parms["doTwo"]        =  True       # Use two baseline combinations in delay cal
parms["delayZeroPhs"] =  False      # Zero phase in Delay solutions?
parms["delayBChan"]   =  None       # first channel to use in delay solutions
parms["delayEChan"]   =  None       # highest channel to use in delay solutions

# Bandpass Calibration?
parms["doBPCal"] =       True       # Determine Bandpass calibration
parms["bpBChan1"] =      1          # Low freq. channel,  initial cal
parms["bpEChan1"] =      0          # Highest freq channel, initial cal, 0=>all
parms["bpDoCenter1"] =   0.1        # Fraction of  channels in 1st, overrides bpBChan1, bpEChan1
parms["bpBChan2"] =      1          # Low freq. channel for BP cal
parms["bpEChan2"] =      0          # Highest freq channel for BP cal,  0=>all
parms["bpChWid2"] =      1          # Number of channels in running mean BP soln
parms["bpdoAuto"] =      False      # Use autocorrelations rather than cross?
parms["bpsolMode"] =     'A&P'      # Band pass type 'A&P', 'P', 'P!A'
parms["bpsolint1"] =     10.0/60.0  # BPass phase correction solution in min
parms["bpsolint2"] =     10.0       # BPass bandpass solution in min
parms["bpUVRange"] =    [0.0,0.0]   # uv range for bandpass cal
parms["specIndex"] =    -0.7        # Spectral index of BP Cal
parms["doSpecPlot"] =    True       # Plot the amp. and phase across the spectrum

# Amp/phase calibration parameters
parms["doAmpPhaseCal"] = True
parms["refAnt"]  =       0          # Reference antenna
parms["refAnts"] =      [0]         # List of Reference antenna for fringe fitting
parms["solInt"]  =      0.5        # solution interval (min)
parms["ampScalar"]=    False        # Ampscalar solutions?
parms["solSmo"]   =    0.0          # Smoothing interval for Amps (min)

# Apply calibration and average?
parms["doCalAvg"] =      True       # calibrate and average cont. calibrator data
parms["avgClass"] =      "UVAvg"    # AIPS class of calibrated/averaged uv data
parms["CalAvgTime"] =    10.0/60.0  # Time for averaging calibrated uv data (min)
parms["CABIF"] =         1          # First IF to copy
parms["CAEIF"] =         0          # Highest IF to copy
parms["CABChan"] =       None       # First Channel to copy
parms["CAEChan"] =       None       # Highest Channel to copy
parms["chAvg"] =         None       # No channel average
parms["avgFreq"] =       None       # No channel average

# Right-Left delay calibration
parms["doRLDelay"] =  False             # Determine/apply R-L delays
parms["RLDCal"]    = [(None,None,None)] # Array of triplets of (name, R-L phase (deg at 1 GHz),
# RM (rad/m**2)) for calibrators
parms["rlBChan"]   = 1                  # First (1-rel) channel number
parms["rlEChan"]   = 0                  # Highest channel number. 0=> high in data.
parms["rlUVRange"] = [0.0,0.0]          # Range of baseline used in kilowavelengths, zeros=all
parms["rlCalCode"] = '  '               # Calibrator code
parms["rlDoCal"]   = 2                  # Apply calibration table? positive=>calibrate
parms["rlgainUse"] = 0                  # CL/SN table to apply, 0=>highest
parms["rltimerange"]= [0.0,1000.0]      # time range of data (days)
parms["rlDoBand"]  = 1                  # If > 0 apply bandpass calibration
parms["rlBPVer"]   = 0                  # BP table to apply, 0=>highest
parms["rlflagVer"] = 2                  # FG table version to apply
parms["rlrefAnt"]  = 0                  # Reference antenna, defaults to refAnt

# Instrumental polarization cal?
parms["doPolCal"]  =  False      # Determine instrumental polarization from PCInsCals?
parms["PCInsCals"] = []          # instrumental poln calibrators, name or list of names
parms["PCFixPoln"] = False       # if True, don't solve for source polarization in ins. cal
parms["PCpmodel"]  = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]  # Instrumental poln cal source poln model.
parms["PCAvgIF"]   = False       # if True, average IFs in ins. cal.
parms["PCSolInt"]  = 2.          # instrumental solution interval (min), 0=> scan average(?)
parms["PCRefAnt"]  = 0           # Reference antenna, defaults to refAnt
parms["PCSolType"] = "    "      # solution type, "    ", "LM  "
parms["doPol"]     = False       # Apply polarization cal in subsequent calibration?
parms["PDVer"]     = 1           # Apply PD table in subsequent polarization cal?
parms["PCChInc"] = 5             # Channel increment in instrumental polarization
parms["PCChWid"] = 5             # Channel averaging in instrumental polarization

# Right-Left phase (EVPA) calibration, uses same  values as Right-Left delay calibration
parms["doRLCal"]    = False    # Set RL phases from RLCal - RLDCal or RLPCal
parms["RLPCal"]     = None     # RL Calibrator source name, in None no IF based cal.
parms["RLPhase"]    = 0.0      # R-L phase of RLPCal (deg) at 1 GHz
parms["RLRM"]       = 0.0      # R-L calibrator RM (NYI)
parms["rlChWid"]    = 3        # Number of channels in running mean RL BP soln
parms["rlsolint1"]  = 10./60   # First solution interval (min), 0=> scan average
parms["rlsolint2"]  = 10.0     # Second solution interval (min)
parms["rlCleanRad"] = None     # CLEAN radius about center or None=autoWin
parms["rlFOV"]      = 0.05     # Field of view radius (deg) needed to image RLPCal

# Recalibration
parms["doRecal"]       = True        # Redo calibration after editing
parms["doDelayCal2"]   = False       # Group Delay calibration of averaged data?, 2nd pass
parms["doBPCal2"]      = True        # Determine Bandpass calibration, 2nd pass
parms["doAmpPhaseCal2"]= True        # Amplitude/phase calibration, 2nd pass
parms["doAutoFlag2"]   = True        # Autoflag editing after final calibration?

# Imaging  targets
parms["doImage"]     = True         # Image targets
parms["targets"]     = []           # List of target sources
parms["outIClass"]   = "IClean"     # Output target final image class
parms["Stokes"]      = "I"          # Stokes to image
parms["Robust"]      = 0.0          # Weighting robust parameter
parms["FOV"]         = None         # Field of view radius in deg.
parms["Niter"]       = 800          # Max number of clean iterations
parms["minFlux"]     = 0.0          # Minimum CLEAN flux density
parms["minSNR"]      = 4.0          # Minimum Allowed SNR
parms["solPMode"]    = "P"          # Phase solution for phase self cal
parms["solPType"]    = "    "       # Solution type for phase self cal
parms["solAMode"]    = "A&P"        # Delay solution for A&P self cal
parms["solAType"]    = "    "       # Solution type for A&P self cal
parms["avgPol"]      = True         # Average poln in self cal?
parms["avgIF"]       = False        # Average IF in self cal?
parms["maxPSCLoop"]  = 2            # Max. number of phase self cal loops
parms["minFluxPSC"]  = 0.025        # Min flux density peak for phase self cal
parms["solPInt"]     = 0.25         # phase self cal solution interval (min)
parms["maxASCLoop"]  = 1            # Max. number of Amp+phase self cal loops
parms["minFluxASC"]  = 0.5          # Min flux density peak for amp+phase self cal
parms["solAInt"]     = 3.0          # amp+phase self cal solution interval (min)
parms["nTaper"]      = 0            # Number of additional imaging multiresolution tapers
parms["Tapers"]      = [20.0,0.0]   # List of tapers in pixels
parms["do3D"]        = False        # Make ref. pixel tangent to celest. sphere for each facet
parms["noNeg"]       = False        # F=Allow negative components in self cal model
parms["BLFact"]      = 1.01         # Baseline dependent time averaging
parms["BLchAvg"]     = True         # Baseline dependent frequency averaging
parms["doMB"]        = None         # Use wideband imaging?
parms["MBnorder"]    = None         # order on wideband imaging
parms["MBmaxFBW"]    = None         # max. MB fractional bandwidth (Set by KAT7InitContFQParms
parms["PBCor"]       = False        # Pri. beam corr on final image
parms["antSize"]     = 12.0         # ant. diameter (m) for PBCor
parms["CleanRad"]    = None         # CLEAN radius (pix?) about center or None=autoWin
parms["xCells"]      = None
parms["yCells"]      = None

# Final
parms["doReport"]  =     True       # Generate source report?
parms["outDisk"]   =     0          # FITS disk number for output (0=cwd)
parms["doSaveUV"]  =     True       # Save uv data
parms["doSaveImg"] =     True       # Save images
parms["doSaveTab"] =     True       # Save Tables
parms["doCleanup"] =     True       # Destroy AIPS files

# diagnostics
parms["plotSource"]    =  'None'      # Name of source for spectral plot
parms["plotTime"]      =  [0.,1000.]  # timerange for spectral plot
parms["doRawSpecPlot"] =  True       # Plot diagnostic raw spectra?
parms["doSpecPlot"]    =  True       # Plot diagnostic spectra?
parms["doSNPlot"]      =  True       # Plot SN tables etc
parms["doDiagPlots"]   =  True       # Plot single source diagnostics
parms["doKntrPlots"]   =  True       # Contour plots
parms["prtLv"]         =  2          # Amount of task print diagnostics
parms["doMetadata"]    =  True       # Save source and project metadata
parms["doHTML"]        =  True       # Output HTML report
parms["doVOTable"]     =  True       # VOTable
