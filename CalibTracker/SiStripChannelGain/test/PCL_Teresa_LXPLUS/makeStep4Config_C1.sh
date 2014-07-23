#!/bin/bash/

file=C1

export filelistName=filelist${file}.txt
export configfileName=step4_ALCAHARVEST_${file}.py
export outputRootFile=Gains_Tree_${file}_Mean.root
MYVARS='$filelistName:$configfileName:$outputRootFile'

cd /afs/cern.ch/user/t/telenz/work/PixelCalibration/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa 

# put this list into the config file
envsubst "$MYVARS" <step4_ALCAHARVEST_withVar.py > ${configfileName}

# to submit this job to the bash system
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`;
cmsRun ${configfileName}

rm ${filelistName}~

# Delete all DQM files besides the last one 
lastline=$(tac ${filelistName} |egrep -m 1 .)
echo $lastline
aux=${#lastline}
echo $aux
runNumber=${lastline:aux-39:6}
echo $runNumber
#find DQM*.root ! -name "*${runNumber}*" -delete