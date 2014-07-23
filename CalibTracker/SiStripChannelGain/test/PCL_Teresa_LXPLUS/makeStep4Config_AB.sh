#!/bin/bash/

file=AB

export filelistName=filelist${file}.txt
export configfileName=step4_ALCAHARVEST_${file}.py
export outputRootFile=Gains_Tree_${file}_Mean.root
MYVARS='$filelistName:$configfileName:$outputRootFile'

cd /afs/cern.ch/user/t/telenz/work/PixelCalibration/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa 

# find the path to the runs
find Run2012A/ -name "PromptCalibProdSiStripGains.root" > filelistA.txt
find Run2012B/ -name "PromptCalibProdSiStripGains.root" > filelistB.txt

# put both together
cat filelistA.txt > ${filelistName}
cat filelistB.txt >> ${filelistName}
sed -i -e 's/^/file:/' ${filelistName}


# put this list into the config file

envsubst "$MYVARS" <step4_ALCAHARVEST_withVar.py > ${configfileName}
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh`;

# to submit this job to the bash system
cmsRun ${configfileName}

#rm filelistA.txt
#rm filelistB.txt
rm ${filelistName}~


lastline=$(tac ${filelistName} |egrep -m 1 .)
runNumber=${lastline:14:6}
#find DQM*.root ! -name "*${runNumber}*" -delete