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
from astropy.io import fits as pyfits

import warnings
warnings.simplefilter('ignore')

usage = "%prog [options] mvffile"
description = "Write uvfits file from mvfv4 file"
parser = OptionParser( usage=usage, description=description)
parser.add_option("--write-flags", action="store_true", default=False, help="Write flags into uvfits file. Writing flags negates weights and cannot be undone inside AIPS!!")
parser.add_option("--channel-range", default=None, help="Range of frequency channels to keep (zero-based inclusive 'first_chan,last_chan', default is all channels)")
parser.add_option("--leave-aips", action="store_true", default=False, help="Don't write out uvfits, just leave the result in an aipsdisk")
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
fitsdisk = 1
nam = filebase
cls = "Raw"
seq = 1

katdata = katdal.open(h5file)

if options.channel_range:
    channel_range = [int(chan_str) for chan_str in options.channel_range.split(',')]
    first_chan, last_chan = channel_range[0], channel_range[1]

    if (first_chan < 0) or (last_chan >= katdata.shape[1]):
        raise RuntimeError("Requested channel range outside data set boundaries. Set channels in the range [0,%s]" % (katdata.shape[1]-1,))

    chan_range = slice(first_chan, last_chan + 1)
    print("\nChannel range %s through %s." % (first_chan, last_chan))
    katdata.select(channels=chan_range)

katdata.select(scans='track')
nbl = len(np.unique([(cp[0][:-1] + cp[1][:-1]).upper() for cp in katdata.corr_products]))

#Condition the uvfits template
templatefile=ObitTalkUtil.FITSDir.FITSdisks[fitsdisk]+'MKATTemplate.uvtab.gz'

KATH5toAIPS.MakeTemplate(templatefile,filebase+'.uvfits',katdata)

uv=OTObit.uvlod(filebase+'.uvfits',0,nam,cls,disk,seq,err)

os.remove(filebase+'.uvfits')

obsdata = KATH5toAIPS.KAT2AIPS(katdata, uv, disk, fitsdisk, err, calInt=1.0, doflags=options.write_flags)

uv.Header(err)
if not options.leave_aips:
    KATCal.KATUVFITS(uv, filebase+'.uv', 0, err)

    if os.path.exists(os.environ['DA00']): shutil.rmtree(os.environ['DA00'])
    d = ObitTalkUtil.AIPSDir.AIPSdisks[disk]
    if os.path.exists(d): shutil.rmtree(d)

