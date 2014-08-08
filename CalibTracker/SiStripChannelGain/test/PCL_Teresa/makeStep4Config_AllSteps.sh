#!/bin/sh

file=${step}

export filelistName=filelist${file}.txt
export configfileName=step4_ALCAHARVEST_${file}.py
export outputRootFile=Gains_Tree_${file}_CL2_New.root
MYVARS='$filelistName:$configfileName:$outputRootFile'

# put this list into the config file
envsubst "$MYVARS" <step4_ALCAHARVEST_withVar.py > ${configfileName}

# to submit this job to the bash system
# Settings to run cmsRUn 
shopt -s expand_aliases
module(){ eval `/usr/bin/modulecmd bash $*`;}
module use -a /afs/desy.de/group/cms/modulefiles/
module load cmssw
SCRIPTPATH=`pwd`
cd $SCRIPTPATH 
eval `scramv1 runtime -sh`
cmsRun ${configfileName}

rm ${filelistName}~

# Delete all DQM files besides the last one
root -b -l -q addDQMFiles.C'("'${step}'")'
rm addDQMFiles_C.* 
while read line
  do
  rm DQM*${line}*
  done<runlist${step}.txt
