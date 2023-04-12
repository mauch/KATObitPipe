#! /usr/bin/env python
import sys
from katim import KATPipe
from optparse import OptionParser

usage = "%prog [options] h5fileToImage"
description = "Make an image from an h5 file."
parser = OptionParser( usage=usage, description=description)
parser.add_option("--outputdir", default='./', help="Specify the output data directory.")
parser.add_option("--parms", dest='parmFile', help="Overwrite the default imaging parameters using a parameter file.")
parser.add_option("--scratchdir", help="Specify the scratch directory.")
parser.add_option("--targets", help="List of targets to load (You'll need calibrators in this list!!).")
parser.add_option("--config", dest='configFile', help="Location of .katimrc configuration file.")
parser.add_option("--timeav", default=1, help="Dumber of dumps to average when making uvfits file")
parser.add_option("--flags", default=None, help="External flags to add")
parser.add_option("--antphase-adjust-filename", default=None, help="string or None\nfilename to numpy .npz file containing phase adjustment per input per channel\nto be applied to all raw visibilities before further processing.")
(options, katfilenames) = parser.parse_args()

if len(katfilenames) == 0:
    parser.print_help()
    sys.exit()

kwargs = {}
for k in ['parmFile', 'scratchdir', 'targets', 'configFile', 'timeav', 'flags', 'antphase_adjust_filename']:
	if getattr(options,k) != None:
		kwargs[k] =  getattr(options,k)
try:
    KATPipe.MKContPipeline(katfilenames, options.outputdir, **kwargs)
finally:
    pass
