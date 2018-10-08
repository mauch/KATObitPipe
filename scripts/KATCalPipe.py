#! /usr/bin/env python
import sys
from katim import KATCalibPipe
from optparse import OptionParser

S3_URL = 'http://stgr2.sdp.mkat.chpc.kat.ac.za:8081/link/s3'
S3_HEAD = {'Content-Type': 'application/json;charset=UTF-8', 'Accept': 'application/json, text/plain, */*'}

def get_archive(katfilenames):
	"""If the names listed in katfilenames don't exist
	try and find them in the archive and download them locally.
	Then return references to existent files on the local disk.
	"""
	import requests

	file_refs = []
	for filename in katfilenames:
		if filename.startswith('s3'):
			res = requests.post(S3_URL, headers=S3_HEAD, data='{"s3_ref":"%s","ref_key":"Nope"}'%(filename,))
			url = res.json()['url']
			res1 = requests.get(url)
			outfile = filename.split('/')[-1]
			open(outfile, 'wb').write(res1.content)
			file_refs.append(outfile)
		else:
			file_refs.append(filename)
	return file_refs

usage = "%prog [options] mvffileToImage"
description = "Calibrate an MVF file."
parser = OptionParser( usage=usage, description=description)
parser.add_option("--outputdir", default='./', help="Specify the output data directory.")
parser.add_option("--parms", dest='parmFile', help="Overwrite the default imaging parameters using a parameter file.")
parser.add_option("--scratchdir", help="Specify the scratch directory.")
parser.add_option("--targets", help="List of targets to load (You'll need calibrators in this list!!).")
parser.add_option("--config", dest='configFile', help="Location of .katimrc configuration file.")
parser.add_option("--timeav", default=1, type='int', help="Number of dumps to average when making uvfits file")
parser.add_option("--flag", action='store_true', default=False, help="Flag data on the fly during conversion")
parser.add_option("--reuse", action='store_true', default=False, help="Reuse already loaded data from aipsdisk")
parser.add_option("--zapraw", action='store_true', default=False, help="zap raw and intermediate uv files")
parser.add_option("--aipsdisk", default='aipsdisk', help='Name of aipsdisk (in "scratchdir" - or cwd) to use - default is "aipsdisk"')
parser.add_option("--halfstokes", default=False, action='store_true', help='Only write out HH,VV when saving uv data')
parser.add_option("--gzip", default=False, action='store_true', help='Gzip the output UV data')
parser.add_option("--dropants", help='List of antennas to remove from pbservation.')
parser.add_option("--blmask", type='float', default=1.e10, help='Baseline length cutoff for the static mask (default apply to all baselines)')
parser.add_option("--refant", type='str', default=None, help='Reference antenna to use for calibration')
parser.add_option("--katdal_refant", type='str', default='', help='Reference antenna to use for activity in katdal')
parser.add_option("--XYfix", default=False, action='store_true', help='Fix the X & Y gains. This will cause the selected target (via --XTYarg) alone'
 									'to be used for delay and bandpass then subsequent gains to be computed with avgPol=True')
parser.add_option("--XYtarg", default='1934-638', help='Name of target to use to fix the X & Y gains (default is 1934-638')

(options, katfilenames) = parser.parse_args()

if len(katfilenames) == 0:
    parser.print_help()
    sys.exit()

katfilenames = get_archive(katfilenames)

kwargs = {}
for k in ['parmFile', 'scratchdir', 'targets', 'configFile', 'timeav', 'flag', 'reuse', 'zapraw', 'aipsdisk', 'halfstokes', 'gzip', 'dropants', 'blmask', 'refant', 'katdal_refant', 'XYfix', 'XYtarg']:
	if getattr(options,k) != None:
		kwargs[k] =  getattr(options,k)
try:
    KATCalibPipe.MKContPipeline(katfilenames, options.outputdir, **kwargs)
finally:
    pass
