#!/bin/bash/

cmsenv

minimumsize=70

fileCheck=true
sizeCheck=true
zombieCheck=true

submit=false 

echo zombieCheck=$zombieCheck
echo sizeCheck=$sizeCheck
echo fileCheck=$fileCheck
echo

find /afs/cern.ch/user/t/telenz/work/PixelCalibration/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa/Run2012CValidation* -name "cmsDriver.sh" -type f > allJobs.txt

while read line 
  do

  aux=${#line}
  path=${line:0:aux-13}
  runNumber=${line:aux-19:6}

  #echo
  #echo "line=" $line
  #echo "path=" $path
  #echo "run number=" $runNumber

# check if file exists
  if [ "$fileCheck" == true ]; then
      if ! [ -a ${path}/PromptCalibProdSiStripGains.root  ] 
	  then
	  echo file does not exist
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      bsub -q 2nd -J gainPCLrun${runNumber} "sh $path/cmsDriver.sh"
	  fi
	  continue
      fi
  fi

# check if file is smaller than 70 kb
  if [ "$sizeCheck" == true ]; then 

      actualsize=$(du -k "${path}/PromptCalibProdSiStripGains.root" | cut -f 1)
      if [ $actualsize -le $minimumsize ]
	  then
	  echo size is below $minimumsize kilobytes
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      #rm $path/PromptCalibProdSiStripGains.root
	      bsub -q 2nd -J gainPCLrun${runNumber} "sh $path/cmsDriver.sh"
	  fi
	  continue
      fi
  fi

# Check if file is a zombie
  if [ "$zombieCheck" == true ]; then
      root -b -l -q isZombieFile.C'("'${path}/PromptCalibProdSiStripGains.root'")' 1> output.txt 2>error_f.txt 
      var=$( tail -1 output.txt )
      #echo $var
      if [ "$var" == "(int)1" ]
	  then 
	  echo "is a Zombie"
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      #rm $path/PromptCalibProdSiStripGains.root
	      bsub -q 2nd -J gainPCLrun${runNumber} sh "${path}/cmsDriver.sh"
	  fi
      fi
  fi

done<allJobs.txt



rm output.txt
rm error_f.txt
rm isZombieFile_C.* 
rm allJobs.txt