#! /usr/bin/env python
import os
import shutil
import numpy as np
import katdal
import OErr, OTObit, ObitTalkUtil, OSystem
import AIPS
from katim import AIPSSetup
from katim import KATH5toAIPS
from katim import KATCal
# Option parser
from optparse import OptionParser
import pyfits

import warnings
warnings.simplefilter('ignore')

usage = "%prog [options] h5file"
description = "Write uvfits file from h5 file"
parser = OptionParser( usage=usage, description=description)
parser.add_option("--write-flags", action="store_true", default=False, help="Write flags into uvfits file (this negates weights and cannot be undone inside AIPS!!)")
parser.add_option("--channel-range", default=None, help="Range of frequency channels to keep (zero-based inclusive 'first_chan,last_chan', default is all channels)")
(options, args) = parser.parse_args()

h5file=args[0]

filebase=os.path.basename(os.path.splitext(h5file)[0])
# Obit error logging
err = OErr.OErr()
OErr.PInit(err,2,'/dev/null')
ObitSys = AIPSSetup.AIPSSetup(err)
# Get the set up AIPS environment.
AIPS_ROOT    = os.environ['AIPS_ROOT']
AIPS_VERSION = os.environ['AIPS_VERSION']

user = OSystem.PGetAIPSuser()
AIPS.userno = user
disk = 1
fitsdisk = 0
nam = filebase
cls = "Raw"
seq = 1

katdata = katdal.open(h5file)

if not options.write_flags:
    #Reset flags in data file
    katdata._flags=np.zeros(katdata.shape,dtype=np.uint8)

if options.channel_range:
    channel_range = [int(chan_str) for chan_str in options.channel_range.split(',')]
    first_chan, last_chan = channel_range[0], channel_range[1]

    if (first_chan < 0) or (last_chan >= katdata.shape[1]):
        raise RuntimeError("Requested channel range outside data set boundaries. Set channels in the range [0,%s]" % (katdata.shape[1]-1,))

    chan_range = slice(first_chan, last_chan + 1)
    print "\nChannel range %s through %s." % (first_chan, last_chan)
    katdata.select(channels=chan_range)

katdata.select(scans='track',targets='1934-638')
numchans = len(katdata.channel_freqs)

#Condition the uvfits template
templatefile=ObitTalkUtil.FITSDir.FITSdisks[fitsdisk]+'MKATTemplate.uvtab.gz'

uvfits = pyfits.open(templatefile)
#Resize the visibility table
vistable = uvfits[1].columns
vistable.del_col('VISIBILITIES')
newvis = pyfits.Column(name='VISIBILITIES',format='%dE'%(3*4*numchans),dim='(3,4,%d,1,1,1)'%(numchans,),array=np.zeros((1,1,1,numchans,4,3,),dtype=np.float32))
vistable.add_col(newvis)
vishdu = pyfits.BinTableHDU.from_columns(vistable)
for key in uvfits[1].header.keys():
    if (key not in vishdu.header.keys()) and (key != 'HISTORY'):
        vishdu.header[key]=uvfits[1].header[key]

newuvfits = pyfits.HDUList([uvfits[0],vishdu,uvfits[2],uvfits[3],uvfits[4],uvfits[5],uvfits[6]])

newuvfits.writeto(filebase+'.uvfits',clobber=True)

uv=OTObit.uvlod(filebase+'.uvfits',0,nam,cls,disk,seq,err)

obsdata = KATH5toAIPS.KAT2AIPS(katdata, uv, disk, fitsdisk, err, calInt=1.0, stop_w=False)

uv.Header(err)

KATCal.KATUVFITab(uv, filebase+'.uv', 0, err)

if os.path.exists(os.environ['DA00']): shutil.rmtree(os.environ['DA00'])
for disk in ObitTalkUtil.AIPSDir.AIPSdisks:
    if os.path.exists(disk): shutil.rmtree(disk)

