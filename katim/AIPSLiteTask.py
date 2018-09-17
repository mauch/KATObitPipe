# Module to assist in running AIPS/ParselTongue on machines without
# an AIPS installation.
#
# Stephen Bourke, JIVE
# Version: 20090204
# Modified by tmauch to work with Obit instead of parseltongue. 20/3/2013

from AIPSTask import *
from AIPSTask import AIPSTask as obAIPSTask
import AIPSLite
import os

class AIPSTask(obAIPSTask):
    """AIPSTask with that uses AIPSLite to fetch Tasks that
    are not available locally."""
    # Package.
    _package = 'AIPS'

    # List of adverbs referring to data.
    _data_adverbs = ['indata', 'outdata',
                     'in2data', 'in3data', 'in4data', 'out2data']

    # List of adverbs referring to disks.
    _disk_adverbs = ['indisk', 'outdisk',
                     'in2disk', 'in3disk', 'in4disk', 'out2disk']

    # List of adverbs referring to file names.
    _file_adverbs = ['infile', 'infile2', 'outfile', 'outprint',
                     'ofmfile', 'boxfile', 'oboxfile']

    # Default version.
    version = os.environ.get('VERSION', 'NEW')

    # Default user number.
    userno = 0

    # Default verbosity level.
    msgkill = 0

    # Default to batch mode.
    isbatch = 32000

    # Run synchronous?
    doWait = False

    # Logging file?
    logFile = ""

    def __init__(self, name, **kwds):
        AIPSLite.get_task(name)
        obAIPSTask.__init__(self, name, **kwds)
