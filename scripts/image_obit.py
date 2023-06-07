#! /usr/bin/env python
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil
import ObitTask
from PipeUtil import setname
import katim.AIPSLiteTask as AIPSTask
from AIPS import AIPSDisk
from FITS import FITSDisk
from PipeUtil import *
from katim.KATCal import *
from katim import KATH5toAIPS
from katim import AIPSSetup
from optparse import OptionParser

usage = "%prog [options] uvfitsfile"
description = "Image a uvfits dataset"
parser = OptionParser( usage=usage, description=description)
parser.add_option("--targets", default=None, help="List of targets to image")
parser.add_option("--refant",type='int',default=0,help="Ref ant to use")
parser.add_option("--blavg",action="store_true",help="BL dep. avg. and Stokes I?")
parser.add_option("--scratch",help="Location of aips disk to create")
(options, args) = parser.parse_args()

filebase=os.path.basename(os.path.splitext(args[0])[0])
# Obit error logging
err = OErr.OErr()
ObitSys = AIPSSetup.AIPSSetup(err,scratchdir=options.scratch)
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
logFile       = "Image.log"   # Processing log file

#astertemplate=ObitTalkUtil.FITSDir.FITSdisks[fitsdisk]+'MKATTemplate.uvtab.gz'
#outtemplate=nam+'.uvtemp'
#uv=OTObit.uvlod(mastertemplate,0,'TEMP',cls,disk,seq,err)
#uv2=OTObit.uvlod(args[0],0,nam,cls,disk,seq,err)
# Add output directory to the environment so AIPS can see it
fn=args[0]
pth,fnn=os.path.split(fn)
if not pth:
    os.environ['FTD']=os.environ['PWD']
else:
    os.environ['FTD']=pth

#fitld=AIPSTask.AIPSTask("uvlod")
print(OSystem.PGetAIPSuser())
#try:
#    fitld.userno = OSystem.PGetAIPSuser()   # This sometimes gets lost
#except Exception, exception:
#    pass
#fitld.doall=1
#fitld.datain='FTD:'+fnn
#fitld.outname='IMAGE'
#fitld.outclass='UV'
#fitld.outseq=1
#fitld.outdisk=1
#fitld.optype='UV'
#fitld.g
uv=OTObit.uvlod(args[0],0,nam,cls,disk,seq,err)
#uv=UV.newPAUV('IMAGE DATA','IMAGE','UV',1,1, True, err)

# Initialize parameters dictionary
####### Initialize parameters dictionary ##### 
parms = KATInitContParms()
targets=options.targets or EVLAAllSource(uv,err,check=False,debug=False)

#BLavg?
if options.blavg:
    blavg=ObitTask.ObitTask('UVBlAvg')
    try:
        blavg.userno = OSystem.PGetAIPSuser()   # This sometimes gets lost
    except Exception as exception:
        pass
    setname(uv,blavg)
    blavg.Sources = targets.split(',')
    #blavg.Stokes=' '
    blavg.FOV=1.3
    blavg.maxInt=1.0
    blavg.maxFact=1.01
    blavg.outDType='AIPS'
    blavg.outName='IMAGE'
    blavg.outClass='BLAVG'
    blavg.outSeq=1
    blavg.outDisk=1
    blavg.logFile=""
    blavg.avgFreq=1
    blavg.chAvg=2
    blavg.g

    uv=UV.newPAUV('BL IMAGE DATA','IMAGE','BLAVG',1,1,True,err)

#uv=OTObit.uvlod(args[0],0,nam,cls,disk,seq,err)

####### Initialize parameters dictionary ##### 
#parms = KATInitContParms()


#Now we should have a uv file in our aips disk
KATImageTargets (uv, err, Sources=targets, seq=1, sclass="IClean", OutlierArea=parms["outlierArea"],\
                          doCalib=-1, doBand=-1,  flagVer=-1, doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                          Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], Niter=parms["Niter"], \
                          CleanRad=parms["CleanRad"], minFlux=parms["minFlux"], OutlierSize=parms["OutlierSize"], \
                          xCells=parms["xCells"], yCells=parms["yCells"], Reuse=parms["Reuse"], minPatch=parms["minPatch"], \
                          maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], noNeg=parms["noNeg"], \
                          solPInt=parms["solPInt"], solPMode=parms["solPMode"], solPType=parms["solPType"], \
                          maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], nx=parms["nx"], ny=parms["ny"], \
                          solAInt=parms["solAInt"], solAMode=parms["solAMode"], solAType=parms["solAType"], \
                          avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR = parms["minSNR"], refAnt=options.refant, \
                          do3D=parms["do3D"], BLFact=1., BLchAvg=parms["BLchAvg"], \
                          doMB=True, norder=parms["MBnorder"], maxFBW=parms["MBmaxFBW"], \
                          PBCor=parms["PBCor"],antSize=parms["antSize"], autoCen=parms["autoCen"], \
                          nTaper=parms["nTaper"], Tapers=parms["Tapers"], sefd=1500,\
                          nThreads=72, noScrat=[], logfile="", check=False, debug=False)

x = Image.newPAImage("out", targets, "IClean", disk, 1, True, err)
xf = KATImFITS(x,targets+".fits", 0 , err, logfile="IMAGE.log")
x = Image.newPAImage("out", targets, "IClean", disk, 1, True, err)
xf = EVLAImFITS (x, targets+".fitab.fits", 0, err, logfile="IMAGE.log")
