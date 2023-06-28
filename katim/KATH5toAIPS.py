""" Utility to convert KAT MVF format data to AIPS (or FITS)

This module requires katdal and katpoint and their dependencies
"""
# $Id: KATH5toAIPS.py 430 2012-11-02 02:00:09Z bill.cotton $
#-----------------------------------------------------------------------
#  Copyright (C) 2012
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
try:
    import katdal
    from katdal import averager
    from katdal.lazy_indexer import DaskLazyIndexer
    from katdal.chunkstore import StoreUnavailable
    import katpoint
except Exception as exception:
    print(exception)
    print("katdal software not available")
    raise  RuntimeError("katdal software unavailable")
else:
    pass
import socket
from collections import namedtuple
import time,math,os
import UV, UVVis, OErr, UVDesc, Table, History
from OTObit import day2dhms
import numpy
import itertools
from astropy.io import fits as pyfits
import multiprocessing
import concurrent.futures
import dask
# Workaround for warnings from dask when using 'dropants' option.
dask.config.set(**{'array.slicing.split_large_chunks': False})
import dask.array as da
import numba
from katsdpsigproc.rfi.twodflag import SumThresholdFlagger
from textwrap import TextWrapper
from .KATImExceptions import KATUnimageableError

""" TextWrapper for wrapping AIPS history text to 70 chars """
_history_wrapper = TextWrapper(width=70, initial_indent='',
                               subsequent_indent='  ',
                               break_long_words=True)

CORR_ID_MAP = {('h', 'h'): 0,
               ('v', 'v'): 1,
               ('h', 'v'): 2,
               ('v', 'h'): 3}

def getmaxstokes(katdata):
    dpmax = -1
    for idx, d in enumerate(katdata.corr_products):
        p1 = d[0][-1:].lower()
        p2 = d[1][-1:].lower()
        dpmax = max(dpmax, CORR_ID_MAP[(p1, p2)])
    # AIPS Stokes index is Fortran.
    dpmax = dpmax + 1
    return dpmax

def KAT2AIPS (katdata, outUV, disk, fitsdisk, err, \
              calInt=1.0, static=None, antphase_adjust_filename=None, \
              quack=1, **kwargs):
    """Convert MeerKAT MVF data set to an Obit UV.

    This module requires katdat and katpoint and their dependencies
    contact Ludwig Schwardt <schwardt@ska.ac.za> for details.

    Parameters
    ----------
    katdata : string
        input katdal object
    outUV : ??
        Obit UV object, shoud be a KAT template for the
        appropriate number of IFs and poln.
    disk  : int
        AIPS Disk number
    fitsdisk: int
        FITS Disk number
    err : ??
        Obit error/message stack
    calInt : 
        Calibration interval in min.
    antphase_adjust_filename : string or None
        filename to numpy .npz file containing phase adjustment per input per channel 
        to be applied to all raw visibilities before further processing.
        npz file contains, e.g. antphasedict_rad={'freqMHz':np.linspace(856,1712,4096,endpoint=False),'m000h':np.zeros(nchans),'m000v':np.zeros(nchans), etc for other inputs};np.savez('antphase_file.npz',**(antphasedict_rad));antphasedict_rad=dict(np.load('antphase_file.npz'))
    quack : number of dumps to drop from the start of each scan
    """
    ################################################################
    OErr.PLog(err, OErr.Info, "Converting MVF data to AIPS UV format.")
    OErr.printErr(err)
    print("Converting MVF data to AIPS UV format.\n")

    # Extract metadata
    meta = GetKATMeta(katdata, err)

    # Extract AIPS parameters of the uv data to the metadata
    meta["Aproject"] = outUV.Aname
    meta["Aclass"] = outUV.Aclass
    meta["Aseq"] = outUV.Aseq
    meta["Adisk"] = disk
    meta["calInt"] = calInt
    meta["fitsdisk"] = fitsdisk
    # Update descriptor
    UpdateDescriptor (outUV, meta, err)
    # Write AN table
    WriteANTable (outUV, meta, err)
    # Write FQ table
    WriteFQTable (outUV, meta, err)
    # Write SU table
    WriteSUTable (outUV, meta, err)

    # Convert data
    ConvertKATData(outUV, katdata, meta, err, static=static, blmask=kwargs.get('blmask',1.e10),
                   antphase_adjust_filename=antphase_adjust_filename, timeav=kwargs.get('timeav',1),
                   flag=kwargs.get('flag',False), doweight=kwargs.get('doweight',True),
                   doflags=kwargs.get('doflags',True), quack=quack)

    # Index data
    OErr.PLog(err, OErr.Info, "Indexing data")
    OErr.printErr(err)
    UV.PUtilIndex (outUV, err)

    # Open/close UV to update header
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")

    # initial CL table
    OErr.PLog(err, OErr.Info, "Create Initial CL table")
    OErr.printErr(err)
    print("Create Initial CL table\n")
    UV.PTableCLfromNX(outUV, meta["maxant"], err, calInt=calInt)
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")

    # History
    outHistory = History.History("outhistory", outUV.List, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp("Convert MeerKAT MVF data to Obit", err)
    for name in katdata.name.split(','):
        for line in _history_wrapper.wrap("datafile = "+name):
            outHistory.WriteRec(-1, line, err)
    outHistory.WriteRec(-1,"calInt   = "+str(calInt), err)
    outHistory.Close(err)
    outUV.Open(UV.READONLY,err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, message="Update UV header failed")
    # Return the metadata for the pipeline
    return meta
   # end KAT2AIPS

def GetKATMeta(katdata, err):
    """
    Get KAT metadata and return as a dictionary.

    Parameters
    ----------
     * katdata  = input KAT dataset
     * err      = Python Obit Error/message stack to init

     Returns : dictionary
     "spw"     Spectral window array as tuple (nchan, freq0, chinc)
               nchan=no. channels, freq0 = freq of channel 0,
               chinc=signed channel increment, one tuple per SPW
     "targets" Array of target tuples:
               (index, name, ra2000, dec2000, raapp, decapp)
     "bpcal"   List of source indices of Bandpass calibrators
     "gaincal" List of source indices of Gain calibrators
     "source" List of source indices of imaging targets
     "targLookup" dict indexed by source number with source index
     "tinteg"   Integration time in seconds
     "obsdate"  First day as YYYY-MM-DD
     "observer" name of observer
     "ants"     Array of antenna tuples (index, name, X, Y, Z, diameter)
     "nstokes"  Number of stokes parameters
     "products" Tuple per data product (ant1, ant2, offset)
                where offset is the index on the Stokes axis (XX=0...)
    """
    ################################################################
    # Checks
    out = {}
    # Spectral windows
    sw = []
    out["spw"] = [(len(katdata.channels), katdata.channel_freqs[0], katdata.channel_freqs[1]-katdata.channel_freqs[0])]
    # targets
    tl = []
    tb = []
    tg = []
    tt = []
    td = {}
    tn = []
    i = 0

    for ti in katdata.target_indices:
        t = katdata.catalogue.targets[ti]
        #Aips doesn't like spaces in names!!
        name = (t.name+"                ")[0:16]
        ra, dec = t.radec().ra.deg, t.radec().dec.deg
        #dec = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        #ra  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
        # Apparent posn
        raa, deca = t.apparent_radec().ra.deg, t.apparent_radec().dec.deg
        #deca = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        #raa  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
        #Avoid duplicates
        if name in tn:
            continue
        tn.append(name)
        i += 1
        tl.append((i, name, ra, dec, raa, deca))
        if 'bpcal' in t.tags:
            tb.append(t)
        elif 'gaincal' in t.tags:
            tg.append(t)
        else:                   # Assume all other targets are for imaging
            tt.append(t)
        td[name.rstrip()] = i
    out["targets"] = tl
    out["targLookup"] = td
    out["bpcal"] = tb
    out["gaincal"] = tg
    out["source"] = tt
    # Antennas
    al = []
    alook = {}
    i = 0
    katdata.ants.sort()
    out["newants"] = katdata.ants
    out["nants"]   = len(katdata.ants)
    antnums = []
    for a in out["newants"]:
        name  = a.name
        x,y,z = a.position_ecef
        diam  = a.diameter.to_value('m')
        i = int(name[1:]) + 1
        antnums.append(i)
        al.append((i, name, x, y, z, diam))
        alook[name] = i
    out["ants"] = al
    out["antLookup"] = alook
    out["maxant"] = max([alook[i] for i in alook])
    # Data products
    nstokes = getmaxstokes(katdata) 
    if not nstokes in [2, 4]:
        msg = f'Only HH,VV (Half pol) or HH,VV,HV,VH (Full pol) datasets can be processed.'
        OErr.printErrMsg(err, msg)
        raise KATUnimageableError(msg)
    #Set up array linking corr products to indices
    dl = numpy.empty((out["nants"], out["nants"], nstokes), dtype=int)
    bl = []
    for idx, d in enumerate(katdata.corr_products):
        a1 = antnums.index(alook[d[0][:-1]])
        a2 = antnums.index(alook[d[1][:-1]])
        p1 = d[0][-1:].lower()
        p2 = d[1][-1:].lower()
        dp = CORR_ID_MAP[(p1, p2)]
        dl[a1, a2, dp] = idx
        #Fill the matrix with the inverse baselines
        dl[a2, a1, dp] = idx
    out["baselines"] = numpy.array([(b[0],b[1]) for b in itertools.combinations_with_replacement(antnums,2)])
    out["blineind"] = numpy.array([(b[0],b[1]) for b in itertools.combinations_with_replacement(list(range(out["nants"])),2)])
    out["products"] = dl
    out["nstokes"]  = nstokes
    # integration time
    out["tinteg"] = katdata.dump_period
    # observing date
    start=time.gmtime(katdata.timestamps[0])
    out["obsdate"] = time.strftime('%Y-%m-%d', start)
    # Observer's name
    out["observer"] = katdata.observer
    # Number of channels
    numchan = len(katdata.channels)
    out["numchan"] = numchan
    # Correlator mode (assuming 1 spectral window KAT-7)
    out["corrmode"] = katdata.spectral_windows[0].product
    # Central frequency (in Hz)
    out["centerfreq"] = katdata.channel_freqs[numchan // 2]
    # Expose all KAT-METADATA to calling script
    out["katdata"] = katdata
    out["RX"] = katdata.spectral_windows[katdata.spw].band
    return out
    # end GetKATMeta

def UpdateDescriptor (outUV, meta, err):
    """
    Update information in data descriptor

    NB: Cannot change geometry of visibilities
    * outUV    = Obit UV object
    * meta     = dict with data meta data
    * err      = Python Obit Error/message stack to init
    """
    ################################################################
    chinc   =  meta["spw"][0][2]   # Frequency increment
    reffreq =  meta["spw"][0][1]   # reference frequency
    nchan   =  meta["spw"][0][0]   # number of channels
    nif     = len(meta["spw"])     # Number of IFs
    nstok   = meta["nstokes"]      # Number of Stokes products
    desc = outUV.Desc.Dict
    outUV.Desc.Dict = desc
    desc['obsdat']   = meta["obsdate"]
    desc['observer'] = meta["observer"]
    desc['JDObs']    = UVDesc.PDate2JD(meta["obsdate"])
    desc['naxis']    = 6
    desc['inaxes']   = [3,nstok,nchan,nif,1,1,0]
    desc['cdelt']    = [1.0,-1.0,chinc, 1.0, 0.0, 0.0, 0.0]
    desc['crval']    = [1.0, -5.0,reffreq, 1.0, 0.0, 0.0, 0.0]
    desc['crota']    = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    outUV.Desc.Dict = desc
    outUV.UpdateDesc(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error updating UV descriptor")
    # end UpdateDescriptor

def WriteANTable(outUV, meta, err):
    """
    Write data in meta to AN table

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    antab = outUV.NewTable(Table.READWRITE, "AIPS AN",1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with AN table")
    antab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening AN table")
    # Update header
    antab.keys['RefDate'] = meta["obsdate"]
    antab.keys['Freq']    = meta["spw"][0][1]
    JD                    = UVDesc.PDate2JD(meta["obsdate"])
    antab.keys['GSTiat0'] = UVDesc.GST0(JD)*15.0
    antab.keys['DEGPDY']  = UVDesc.ERate(JD)*360.0
    Table.PDirty(antab)
    # Force update
    row = antab.ReadRow(1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading AN table")
    OErr.printErr(err)
    irow = 0
    for ant in meta["ants"]:
        irow += 1
        row['NOSTA']    = [ant[0]]
        row['ANNAME']   = [ant[1]+"    "]
        row['STABXYZ']  = [ant[2],ant[3],ant[4]]
        row['DIAMETER'] = [ant[5]]
        row['POLAA']    = [90.0]
        antab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing AN table")
    antab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing AN table")
    # end WriteANTable

def WriteFQTable(outUV, meta, err):
    """
    Write data in meta to FQ table
    An old FQ table is deleted

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    # If an old table exists, delete it
    if outUV.GetHighVer("AIPS FQ")>0:
        zz = outUV.ZapTable("AIPS FQ", 1, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error zapping old FQ table")
    reffreq =  meta["spw"][0][1]   # reference frequency
    noif    = 1     # Number of IFs (1 always for KAT7)
    fqtab = outUV.NewTable(Table.READWRITE, "AIPS FQ",1,err,numIF=noif)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FQ table")
    fqtab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening FQ table")
    # Update header
    fqtab.keys['NO_IF'] = 1  # Structural so no effect
    Table.PDirty(fqtab)  # Force update
    # Create row
    row = {'FRQSEL': [1], 'CH WIDTH': [0.0], 'TOTAL BANDWIDTH': [0.0], \
           'RXCODE': ['L'], 'SIDEBAND': [-1], 'NumFields': 7, 'Table name': 'AIPS FQ', \
           '_status': [0], 'IF FREQ': [0.0]}
    if err.isErr:
        OErr.printErrMsg(err, "Error reading FQ table")
    OErr.printErr(err)
    irow = 0
    for sw in meta["spw"]:
        irow += 1
        row['FRQSEL']    = [irow]
        row['IF FREQ']   = [sw[1] - reffreq]
        row['CH WIDTH']  = [sw[2]]
        row['TOTAL BANDWIDTH']  = [abs(sw[2])*sw[0]]
        row['RXCODE']  = [meta['RX']]
        if sw[2]>0.0:
            row['SIDEBAND']  = [1]
        else:
            row['SIDEBAND']  = [-1]
        fqtab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing FQ table")
    fqtab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing FQ table")
    # end WriteFQTable

def WriteFGTable(outUV, katdata, meta, err):
    """
    Get the flags from the h5 file and convert to FG table format.
    UNUSED- This implimentation is too slow!
    outUV    = Obit UV object
    meta     =  dict with data meta data
    err      = Python Obit Error/message stack to init
    """
    ###############################################################

    # work out Start time in unix sec
    tm = katdata.timestamps[1:2]
    tx = time.gmtime(tm[0])
    time0   = tm[0] - tx[3]*3600.0 - tx[4]*60.0 - tx[5]
    int_time = katdata.dump_period

    #Loop through scans in h5 file
    row = 0
    for scan, state, target in katdata.scans():
        name=target.name.replace(' ','_')
        if state != 'track':
            continue
        tm = katdata.timestamps[:]
        nint=len(tm)
        el=target.azel(tm[int(nint/2)])[1]*180./math.pi
        if el<15.0:
            continue
        row+=1
        flags = katdata.flags()[:]
        numflag = 0
        for t, chan_corr in enumerate(flags):
            for c, chan in enumerate(chan_corr):
                cpflagged=[]
                for p, cp in enumerate(chan):
                #for cpaverage in meta['pairLookup']:
                    flag=cp #numpy.any(chan[meta['pairLookup'][cpaverage]])
                    product=meta['products'][p]
                    if product[0] == product[1]:
                        continue
                    if flag and (not product[0:2] in cpflagged):
                        cpflagged.append(product[0:2])
                        numflag+=1
                        starttime=float((tm[t]-time0 - (int_time/2))/86400.0)
                        endtime=float((tm[t]-time0 + (int_time/2))/86400.0)
                        UV.PFlag(outUV,err,timeRange=[starttime,endtime], flagVer=1, \
                                     Ants=[product[0],product[1]], \
                                     Chans=[c+1,c+1], IFs=[1,1], Stokes='1111', Reason='Online flag')
        numvis=t*c*(p/meta["nstokes"])
        msg = "Scan %4d %16s   Online flags: %7s of %8s vis"%(row,name,numflag,numvis)
        OErr.PLog(err, OErr.Info, msg);
        OErr.printErr(err)

def WriteSUTable(outUV, meta, err):
    """
    Write data in meta to SU table

     * outUV    = Obit UV object
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
    """
    ################################################################
    sutab = outUV.NewTable(Table.READWRITE, "AIPS SU",1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with SU table")
    sutab.Open(Table.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening SU table")
    # Update header
    sutab.keys['RefDate'] = meta["obsdate"]
    sutab.keys['Freq']    = meta["spw"][0][1]
    Table.PDirty(sutab)  # Force update
    row = sutab.ReadRow(1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading SU table")
    OErr.printErr(err)
    irow = 0
    for tar in meta["targets"]:
        irow += 1
        row['ID. NO.']   = [tar[0]]
        row['SOURCE']    = [tar[1]]
        row['RAEPO']     = [tar[2]]
        row['DECEPO']    = [tar[3]]
        row['RAOBS']     = [tar[2]]
        row['DECOBS']    = [tar[3]]
        row['EPOCH']     = [2000.0]
        row['RAAPP']     = [tar[4]]
        row['DECAPP']    = [tar[5]]
        row['BANDWIDTH'] = [meta["spw"][0][2]]
        sutab.WriteRow(irow, row,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error writing SU table")
    sutab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing SU table")
    # end WriteSUtable

def StopFringes(visData,freqData,cable_delay):
    """
    Stop the fringes for a KAT-7 antenna/polarisation pair.

    * visData  = input array of visibilities
    * freqData = the frequencies corresponding to each channel in visData
    * wData    = the w term for each channel in VisData
    * polProd  = the polarisation product which is being stopped

    Outputs an array of the same size as visData with fringes stopped
    """

    # KAT-7 Antenna Delays (From h5toms.py)
    # Result of delay model in turns of phase. This is now frequency dependent so has shape (tstamps, channels)
    turns = numpy.outer(cable_delay, freqData)
    outVisData = visData*numpy.exp(-2j * numpy.pi * turns)

    return outVisData


def load_phase_correction(antphase_adjust_filename, katdata, err):
    """Blame mattieu@sarao.ska.ac.za"""
    if antphase_adjust_filename is None:
        return None
    msg='Adjusting visibility phases from file: %s\n'%(antphase_adjust_filename)
    OErr.PLog(err, OErr.Info, msg)
    OErr.printErr(err)
    print(msg)
    #dictionary with items e.g. antphasedict_rad={'freqMHz':np.linspace(856,1712,4096,endpoint=False),'m000h':np.zeros(nchans),'m000v':np.zeros(nchans)}
    antphasedict_rad=dict(numpy.load(antphase_adjust_filename))
    #check all is OK (any inputs missing, channels match)
    missinginputs=[]
    for input0 in katdata.inputs:
        if input0 not in antphasedict_rad:
            missinginputs.append(input0)
    if len(missinginputs):
        msg='Warning: The phase adjustment is set to zero for the following inputs that are missing from %s: %s'% \
            (antphase_adjust_filename,','.join(missinginputs))
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
        print(msg)
        for input0 in missinginputs:
            antphasedict_rad[input0]=np.zeros(len(antphasedict_rad['freqMHz']))
    if (katdata.channel_width/1e6!=antphasedict_rad['freqMHz'][1]-antphasedict_rad['freqMHz'][0]):
        msg='Error: channel_width (%g MHz) does not match that of %s (%g MHz)'% \
            (katdata.channel_width/1e6,antphase_adjust_filename,antphasedict_rad['freqMHz'][1]-antphasedict_rad['freqMHz'][0])
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
        print(msg)
        return None
    if (katdata.channels[-1]>len(antphasedict_rad['freqMHz'])):
        msg='Error: incorrect number of channels (%d) in %s'% \
            (len(antphasedict_rad['freqMHz']),antphase_adjust_filename)
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
        print(msg)
        return None
    phase_corr = numpy.zeros(katdata.shape[1:], dtype=katdata.vis.dtype)
    for iprod,(input0,input1) in enumerate(katdata.corr_products):
        phase_corr[:,iprod]=numpy.exp(1j*(antphasedict_rad[input0][katdata.channels] -
                                      antphasedict_rad[input1][katdata.channels]))
    return phase_corr


def ConvertKATData(outUV, katdata, meta, err, static=None, blmask=1.e10, antphase_adjust_filename=None,
                   timeav=1, flag=False, doweight=True, doflags=True, quack=1):
    """
    Read KAT HDF data and write Obit UV

     outUV    = Obit UV object
     katdata  = input KAT dataset
     meta     = dict with meta data
     err      = Python Obit Error/message stack to init
     static   = Numpy array of static flags per channel
     blmask   = Baseline length in metres to apply the mask
     antphase_adjust_filename = Path to npz file contain per antenna phase adjustment
     timeav   = number of dumps to average in time
     flag     = Run the SDP flagger on the data?
     doweight = Read the weights in the katdal file? (else they will all be 1.0)
     doflags  = read the flags from the katdal file?
     quack    = number of dumps to 'quack' before each scan.
    """
    ################################################################
    reffreq =  meta["spw"][0][1]    # reference frequency
    lamb    = 2.997924562e8 / reffreq # wavelength of reference freq
    nchan   =  meta["spw"][0][0]    # number of channels
    nif     = len(meta["spw"])      # Number of IFs
    nstok   = meta["nstokes"]       # Number of Stokes products
    newants = meta["newants"]
    p       = meta["products"]      # baseline stokes indices
    b       = meta["baselines"]
    bi      = meta["blineind"]
    nbase   = b.shape[0]                # number of correlations/baselines
    nprod   = nbase * nstok
    uvw_indices = p[bi[:, 0], bi[:, 1], 0]
    antslookup = meta["antLookup"]
    # work out Start time in unix sec
    tm = katdata.timestamps[0]
    tx = time.gmtime(tm)
    time0   = tm - tx[3]*3600.0 - tx[4]*60.0 - tx[5]

    # Set data to read one timestamp per IO
    outUV.List.set("nVisPIO", nbase)
    d = outUV.Desc.Dict
    d.update(numVisBuff=nbase)
    outUV.Desc.Dict = d
    # Open data
    zz = outUV.Open(UV.READWRITE, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening output UV")
    # visibility record offsets
    idb = {}
    idb['ilocu']   = d['ilocu']
    idb['ilocv']   = d['ilocv']
    idb['ilocw']   = d['ilocw']
    idb['iloct']   = d['iloct']
    idb['ilocb']   = d['ilocb']
    idb['ilocsu']  = d['ilocsu']
    idb['nrparm']  = d['nrparm']

    visshape = nchan * nstok * 3
    count = 0.0
    # Get IO buffers as numpy arrays
    buff =  numpy.frombuffer(outUV.VisBuf, dtype=numpy.float32)
    #Set up a flagger if needs be
    if flag:
        flagger = SumThresholdFlagger(outlier_nsigma=4.5, freq_chunks=7,
                                      spike_width_freq=1.5e6/katdata.channel_width,
                                      spike_width_time=100./katdata.dump_period,
                                      time_extend=3, freq_extend=3, average_freq=1)
    #Set up thestatic flags
    if static is not None:
        static_flags = get_static_flags(katdata, blmask, static)

    # Template vis
    vis = outUV.ReadVis(err, firstVis=1)
    first = True
    visno = 1
    numflags = 0
    numvis = 0
    #load per-antenna,per-channel phase adjustment, if any, from file to apply to visibilities later
    visphase_corr=load_phase_correction(antphase_adjust_filename, katdata, err)

    max_scan = 151
    quack = quack
    # Generate arrays for storage
    scan_vs = numpy.empty((max_scan, nchan, nprod), dtype=katdata.vis.dtype)
    scan_fg = numpy.empty((max_scan, nchan, nprod), dtype=katdata.flags.dtype)
    scan_wt = numpy.empty((max_scan, nchan, nprod), dtype=katdata.weights.dtype)
    for scan, state, target in katdata.scans():
        # Don't read at all if all will be "Quacked"
        if katdata.shape[0] < ((quack + 1) * timeav):
            continue
        # Chunk data into max_scan dumps
        if katdata.shape[0] > max_scan:
            scan_slices = [slice(i, i + max_scan, 1) for i in range(quack * timeav, katdata.shape[0], max_scan)]
            scan_slices[-1] = slice(scan_slices[-1].start, katdata.shape[0], 1)
        else:
            scan_slices = [slice(quack * timeav, katdata.shape[0])]

        # Number of integrations
        num_ints = katdata.timestamps.shape[0] - quack * timeav
        msg = "Scan:%4d Int: %4d %16s Start %s"%(scan, num_ints, target.name,
                                                 day2dhms((katdata.timestamps[0] - time0) / 86400.0)[0:12])
        OErr.PLog(err, OErr.Info, msg);
        OErr.printErr(err)
        print(msg)
        for sl in scan_slices:
            tm = katdata.timestamps[sl]
            nint = tm.shape[0]
            load(katdata, numpy.s_[sl.start:sl.stop, :, :], scan_vs[:nint], scan_wt[:nint], scan_fg[:nint], err)
            # Make sure we've reset the weights
            wt = scan_wt[:nint]
            if doweight==False:
                wt[:] = 1.
            vs = scan_vs[:nint]
            if visphase_corr is not None: # applies if antphasedict if supplied
                vs *= visphase_corr[numpy.newaxis, :]
            fg = scan_fg[:nint]
            if doflags==False:
                fg[:] = False
            if static is not None:
                fg |= static_flags[numpy.newaxis, :]
            if flag:
                fg |= flag_data(vs, fg, flagger)
            if timeav>1:
                vs, wt, fg, tm, _ = averager.average_visibilities(vs, wt, fg, tm, katdata.channel_freqs, timeav=int(timeav), chanav=1)
                # Update number of integrations for averaged data.
                nint = tm.shape[0]
            # Get target suid
            # Only on targets in the input list
            try:
                suid = meta["targLookup"][target.name[0:16]]
            except:
                continue

            numflags += numpy.sum(fg)
            numvis += fg.size
            # uvw calculation
            u = katdata.u[sl, uvw_indices]
            v = katdata.v[sl, uvw_indices]
            w = katdata.w[sl, uvw_indices]
            uvw_coordinates = numpy.stack((u, v, w), axis=-1)
            # Convert to aipsish
            uvw_coordinates /= lamb

            # Convert to AIPS time
            tm = (tm - time0) / 86400.0

            #Get random parameters for this scan
            rp = get_random_parameters(idb, b, uvw_coordinates, tm, suid)
            # Loop over integrations
            for iint in range(0, nint):
                # Fill the buffer for this integration
                buff = fill_buffer(vs[iint], fg[iint], wt[iint], rp[iint], p, bi, buff)
                # Write to disk
                outUV.Write(err, firstVis=visno)
                visno += nbase
                firstVis = None
            # end loop over integrations
            if err.isErr:
                OErr.printErrMsg(err, "Error writing data")
    # end loop over scan
    if numvis>0:
        msg= "Applied %s online flags to %s visibilities (%.3f%%)"%(numflags,numvis,(float(numflags)/float(numvis)*100.))
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
    outUV.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing data")
    # end ConvertKATData

def MakeTemplate(inuv, outuv, katdata):
    """
    Construct a template file with the correct channel range and write it to outuv.
    """
    numchans = len(katdata.channel_freqs)
    numstokes = getmaxstokes(katdata)
    nvispio = len(numpy.unique([(cp[0][:-1] + cp[1][:-1]).upper() for cp in katdata.corr_products]))
    uvfits = pyfits.open(inuv)
    #Resize the visibility table
    vistable = uvfits[1].columns
    vistable.del_col('VISIBILITIES')
    newvis = pyfits.Column(name='VISIBILITIES', format='%dE'%(3*numstokes*numchans),
                           dim='(3,%d,%d,1,1,1)'%(numstokes,numchans,),
                           array=numpy.zeros((1,1,1,numchans,numstokes,3,), dtype=numpy.float32))
    vistable.add_col(newvis)
    vishdu = pyfits.BinTableHDU.from_columns(vistable)
    for key in list(uvfits[1].header.keys()):
        if (key not in list(vishdu.header.keys())) and (key != 'HISTORY'):
            vishdu.header[key]=uvfits[1].header[key]

    newuvfits = pyfits.HDUList([uvfits[0],vishdu,uvfits[2],uvfits[3],uvfits[4],uvfits[5],uvfits[6]])
    #Add nvispio rows
    newuvfits.writeto(outuv, overwrite=True)
    uvfits = pyfits.open(outuv, mode='update')
    uvfits[1].data = uvfits[1].data.repeat(nvispio)
    uvfits.flush()

def flag_data(vs,fg,flagger):
    """
    Flag the data using flagger.
    """
    with concurrent.futures.ThreadPoolExecutor(multiprocessing.cpu_count()) as pool:
        detected_flags = flagger.get_flags(vs, fg, pool)
    return detected_flags

def get_uvw_coordinates(array_centre, baseline_vectors, tm, target, b):
    # uvw calculation
    a1 = b[:,0]
    a2 = b[:,1]
    uvw_basis = target.uvw_basis(tm, array_centre)
    uvw_ant = numpy.tensordot(baseline_vectors, uvw_basis, ([1], [1]))
    uvw_ant = numpy.transpose(uvw_ant, (2, 0, 1))
    uvw_coordinates = (numpy.take(uvw_ant, a1, axis=1)
                            - numpy.take(uvw_ant, a2, axis=1))
    return uvw_coordinates

def get_random_parameters(dbiloc, bl, uvws, tm, suid):
    """
    Construct an array of shape (len(tm), len(bl),len(dbiloc)) containing the
    AIPS random parameters as a function of baseline and timestamp.
    Array returned has order: specified in dbiloc
    """

    rp = numpy.empty((len(tm), len(bl), dbiloc['nrparm'],), dtype=numpy.float32)

    # Get baselines in aips format
    aips_bl = 256.0*(bl[:, 0]) + (bl[:, 1])
    # uvw
    rp[..., dbiloc['ilocu']] = uvws[..., 0]
    rp[..., dbiloc['ilocv']] = uvws[..., 1]
    rp[..., dbiloc['ilocw']] = uvws[..., 2]
    # time
    rp[..., dbiloc['iloct']] = tm[:, numpy.newaxis]
    # baseline
    rp[..., dbiloc['ilocb']] = aips_bl[numpy.newaxis, :]
    # source
    rp[..., dbiloc['ilocsu']] = suid

    return rp


# Number of times to try reading each chunk before giving up
NUM_RETRIES = 3
# Number of dumps to read at at time
CHUNK_SIZE = 2

def load(dataset, indices, vis, weights, flags, err):
    """Load data from lazy indexers into existing storage.
    This is optimised for the MVF v4 case where we can use dask directly
    to eliminate one copy, and also load vis, flags and weights in parallel.
    In older formats it causes an extra copy.
    Parameters
    ----------
    dataset : :class:`katdal.DataSet`
        Input dataset, possibly with an existing selection
    indices : tuple
        Slice expression for subsetting the dataset
    vis, flags : array-like
        Outputs, which must have the correct shape and type
    """

    t_min = indices[0].start
    t_max = indices[0].stop
    in_time_slices = [slice(ts, min(ts+CHUNK_SIZE, t_max)) for ts in range(t_min, t_max, CHUNK_SIZE)]
    for in_ts in in_time_slices:
        out_ts = slice(in_ts.start - t_min, in_ts.stop - t_min)
        out_vis = vis[out_ts]
        out_weights = weights[out_ts]
        out_flags = flags[out_ts]
        for i in range(NUM_RETRIES):
            try:
                if isinstance(dataset.vis, DaskLazyIndexer):
                    DaskLazyIndexer.get([dataset.vis, dataset.weights, dataset.flags], in_ts, out=[out_vis, out_weights, out_flags])
                else:
                    out_vis[:] = dataset.vis[in_ts]
                    out_weights[:] = dataset.weights[in_ts]
                    out_flags[:] = dataset.flags[in_ts]
                break
            except (StoreUnavailable, socket.timeout):
                msg = 'Timeout when reading dumps %d to %d. Try %d/%d....' % (out_ts.start + 1, out_ts.stop, i + 1, NUM_RETRIES)
                OErr.PLog(err, OErr.Warn, msg);
                OErr.printErr(err)
                print(msg)
        # Flag the data and warn if we can't get it
        if i == NUM_RETRIES - 1:
            msg = 'Too many timeouts, flagging dumps %d to %d' % (out_ts.start + 1, out_ts.stop)
            OErr.PLog(err, OErr.Warn, msg);
            OErr.printErr(err)
            print(msg)
            flags[out_ts] = True

@numba.jit(nopython=True, parallel=True)
def fill_buffer(in_vis, in_flags, in_weights, in_rparm, cp_index, bls_index, out_buffer, or_flags_pols=True):
    """Reorganise baselines and axis order.
    The inputs have dimensions (channel, pol-baseline), and the output 
    is a 1d array buffer written to aips, random parameters are written 
    before each visibility. Flags are applied by negating the weights.
    cp_index is a 3D array which is indexed by
    ant1, ant2 and pol to get the input pol-baseline.
    Stolen from katdal, mvftoms
    """
    in_flags_u8 = in_flags.view(numpy.uint8)
    n_bls = bls_index.shape[0]
    n_chans = in_vis.shape[0]
    n_pols = cp_index.shape[2]
    n_rparm = in_rparm.shape[1]
    vis_step = n_rparm + (n_chans * n_pols * 3)
    bstep = 128
    bblocks = (n_bls + bstep - 1) // bstep
    for bblock in numba.prange(bblocks):
        bstart = bblock * bstep
        bstop = min(n_bls, bstart + bstep)
        for b in range(bstart, bstop):
            a1, a2 = bls_index[b]
            thisrparm = in_rparm[b]
            for r in range(n_rparm):
                out_buffer[(b * vis_step) + r] = thisrparm[r]
            for c in range(n_chans):
                if or_flags_pols:
                    p_flag = False
                    # OR the flags over pols
                    for p in range(n_pols):
                        idx = cp_index[a1, a2, p]
                        p_flag |= in_flags_u8[c, idx] > 0
                for p in range(n_pols):
                    idx = cp_index[a1, a2, p]
                    vis = in_vis[c, idx]
                    if or_flags_pols:
                        flg = p_flag
                    else:
                        flg = in_flags_u8[c, idx] > 0
                    if flg:
                        weight = -32767.
                    else:
                        weight = in_weights[c, idx]
                    buff_idx = (b * vis_step) + n_rparm + (c * n_pols * 3) + (p * 3)
                    out_buffer[buff_idx] = vis.real
                    out_buffer[buff_idx + 1] = vis.imag
                    out_buffer[buff_idx + 2] = weight
    return out_buffer

def get_static_flags(katdata, limit, static):
    """
    Compute a mask of the same shape as katdata.corr_products that indicates
    whether the basline length of the given correlation product is
    shorter than limit in meters
    """
    static_flags = numpy.zeros(katdata.shape[1:], dtype = katdata.flags.dtype)
    antlookup = {}
    for ant in katdata.ants:
        antlookup[ant.name] = ant
    for prod, baseline in enumerate(katdata.corr_products):
        bl_vector = antlookup[baseline[0][:4]].baseline_toward(antlookup[baseline[1][:4]])
        bl_length = bl_vector.norm().to_value('m')
        if bl_length < limit:
            static_flags[:, prod] = static
    return static_flags

