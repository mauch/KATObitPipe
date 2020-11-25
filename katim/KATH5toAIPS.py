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
    print("KAT software not available")
    raise  RuntimeError("KAT software unavailable")
else:
    pass
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
import dask.array as da
import numba
from katsdpsigproc.rfi.twodflag import SumThresholdFlagger


def KAT2AIPS (katdata, outUV, disk, fitsdisk, err, \
              calInt=1.0, static=None, **kwargs):
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
    targets : list, optinal
        List of targetnames to extract from the file
    stop_w : bool
        Fring stop data? (Values only for KAT-7)
    """
    ################################################################
    OErr.PLog(err, OErr.Info, "Converting MVF data to AIPS UV format.")
    OErr.printErr(err)
    print("Converting MVF data to AIPS UV format.\n")

    # Extract metadata
    meta = GetKATMeta(katdata, err)

    # TODO: Fix this all up so that the below isn't the case!
    if meta["products"].size != meta["nants"] * meta["nants"] * 4:
        raise ValueError("Only full stokes and all correlation products are supported.")

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
    ConvertKATData(outUV, katdata, meta, err, static=static, blmask=kwargs.get('blmask',1.e10), stop_w=kwargs.get('stop_w',False), timeav=kwargs.get('timeav',1), flag=kwargs.get('flag',False), doweight=kwargs.get('doweight',True), doflags=kwargs.get('doflags',True))

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
        outHistory.WriteRec(-1,"datafile = "+name, err)
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
        ras, decs = t.radec()
        dec = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        ra  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
        # Apparent posn
        ras, decs = t.apparent_radec()
        deca = UVDesc.PDMS2Dec(str(decs).replace(':',' '))
        raa  = UVDesc.PHMS2RA(str(ras).replace(':',' '))
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
        diam  = a.diameter
        i = int(name[1:]) + 1
        antnums.append(i)
        al.append((i, name, x, y, z, diam))
        alook[name] = i
    out["ants"] = al
    out["antLookup"] = alook
    out["maxant"] = max([alook[i] for i in alook])
    # Data products
    nstokes = 4
    #Set up array linking corr products to indices
    dl = numpy.empty((out["nants"], out["nants"], nstokes), dtype=numpy.int)
    bl = []
    for idx, d in enumerate(katdata.corr_products):
        a1 = antnums.index(alook[d[0][:4]])
        a2 = antnums.index(alook[d[1][:4]])
        if d[0][4:]=='h' and d[1][4:]=='h':
            dp = 0
        elif d[0][4:]=='v' and d[1][4:]=='v':
            dp = 1
        elif d[0][4:]=='h' and d[1][4:]=='v':
            dp = 2
        else:
            dp = 3
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
    out["centerfreq"] = katdata.channel_freqs[numchan //2]
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


def ConvertKATData(outUV, katdata, meta, err, static=None, blmask=1.e10, stop_w=False, timeav=1, flag=False, doweight=True, doflags=True):
    """
    Read KAT HDF data and write Obit UV

     * outUV    = Obit UV object
     * katdata  = input KAT dataset
     * meta     = dict with data meta data
     * err      = Python Obit Error/message stack to init
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
    #Set up the baseline mask
    blmask = get_baseline_mask(newants, katdata.corr_products, blmask)

    # Template vis
    vis = outUV.ReadVis(err, firstVis=1)
    first = True
    visno = 1
    numflags = 0
    numvis = 0
    # Do we need to stop Fringes
    if stop_w:
        msg = "W term in UVW coordinates will be used to stop the fringes."
        OErr.PLog(err, OErr.Info, msg)
        OErr.printErr(err)
        print(msg)

    # Set up baseline vectors of uvw calculation
    array_centre = katpoint.Antenna('', *newants[0].ref_position_wgs84)
    baseline_vectors = numpy.array([array_centre.baseline_toward(antenna)
                                for antenna in newants])

    max_scan = 151
    QUACK = 1
    # Generate arrays for storage
    scan_vs = numpy.empty((max_scan, nchan, nprod), dtype=katdata.vis.dtype)
    scan_fg = numpy.empty((max_scan, nchan, nprod), dtype=katdata.flags.dtype)
    scan_wt = numpy.empty((max_scan, nchan, nprod), dtype=katdata.weights.dtype)
    for scan, state, target in katdata.scans():
        # Don't read at all if all will be "Quacked"
        if katdata.shape[0] < ((QUACK + 1) * timeav):
            continue
        # Chunk data into max_scan dumps
        if katdata.shape[0] > max_scan:
            scan_slices = [slice(i, i + max_scan, 1) for i in range(QUACK * timeav, katdata.shape[0], max_scan)]
            scan_slices[-1] = slice(scan_slices[-1].start, katdata.shape[0], 1)
        else:
            scan_slices = [slice(QUACK * timeav, katdata.shape[0])]

        # Number of integrations
        num_ints = katdata.timestamps.shape[0] - QUACK * timeav
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
            fg = scan_fg[:nint]
            if doflags==False:
                fg[:] = False
            if static is not None:
                fg[:, :, blmask] |= static[numpy.newaxis, :, numpy.newaxis]
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
            uvw_coordinates = get_uvw_coordinates(array_centre, baseline_vectors, tm, target, bi)

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
    nvispio = len(numpy.unique([(cp[0][:-1] + cp[1][:-1]).upper() for cp in katdata.corr_products]))
    uvfits = pyfits.open(inuv)
    #Resize the visibility table
    vistable = uvfits[1].columns
    vistable.del_col('VISIBILITIES')
    newvis = pyfits.Column(name='VISIBILITIES',format='%dE'%(3*4*numchans),dim='(3,4,%d,1,1,1)'%(numchans,),array=numpy.zeros((nvispio,1,1,1,numchans,4,3,),dtype=numpy.float32))
    vistable.add_col(newvis)
    vishdu = pyfits.BinTableHDU.from_columns(vistable)
    for key in list(uvfits[1].header.keys()):
        if (key not in list(vishdu.header.keys())) and (key != 'HISTORY'):
            vishdu.header[key]=uvfits[1].header[key]

    newuvfits = pyfits.HDUList([uvfits[0],vishdu,uvfits[2],uvfits[3],uvfits[4],uvfits[5],uvfits[6]])
    #Add nvispio rows
    newuvfits.writeto(outuv, overwrite=True)

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
            except StoreUnavailable:
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

def get_baseline_mask(ants, corr_prods, limit):
    """
    Compute a mask of the same shape as corr_products that indicates
    whether the basline length of the given correlation product is
    shorter than limit in meters
    """
    baseline_mask = numpy.zeros(corr_prods.shape[0], dtype=numpy.bool)
    antlookup = {}
    for ant in ants:
        antlookup[ant.name] = ant
    for prod, baseline in enumerate(corr_prods):
        bl_vector = antlookup[baseline[0][:4]].baseline_toward(antlookup[baseline[1][:4]])
        bl_length = numpy.linalg.norm(bl_vector)
        if bl_length < limit:
            baseline_mask[prod] = True
    return baseline_mask

