#!/bin/sh

file=${step}

export filelistName=filelist${file}.txt
export configfileName=step4_ALCAHARVEST_${file}.py
export outputRootFile=Gains_Tree_${file}.root
MYVARS='$filelistName:$configfileName:$outputRootFile'

# put this list into the config file
envsubst "$MYVARS" <step4_ALCAHARVEST_withVar.py > ${configfileName}

# to submit this job to the bash system
# Settings to run cmsRUn 
shopt -s expand_aliases
module(){ eval `/usr/bin/modulecmd bash $*`;}
module use -a /afs/desy.de/group/cms/modulefiles/
module load cmssw
cd /nfs/dust/cms/user/tlenz/PixelCalibrationCode/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa 
eval `scramv1 runtime -sh`
cmsRun ${configfileName}

rm ${filelistName}~

# Delete all DQM files besides the last one
lastline=$(tac ${filelistName} |egrep -m 1 .)
aux=${#lastline}
lastRunNumber=${lastline:aux-39:6}
i=0
while read line
  do
  runNumber=${line:aux-39:6}
  if [ "$runNumber" == "$lastRunNumber" ]
      then
      echo "break"
      break
  fi
  rm DQM*${runNumber}*
  i=$(( $i + 1 ))
  done<$filelistName
