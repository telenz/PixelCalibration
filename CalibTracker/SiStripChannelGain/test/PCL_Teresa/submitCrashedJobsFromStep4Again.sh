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
SCRIPTPATH=`pwd`

find $SCRIPTPATH/MC2012* -name "cmsDriver.sh" -type f > allJobs.txt

while read line 
  do

  aux=${#line}
  path=${line:0:aux-13}
  runNumber=${line:aux-19:6}
  # only for MC 
  runNumber=$(echo $runNumber | cut -f2 -d"/")
   

# echo
# echo "line=" $line
# echo "path=" $path
# echo "run number=" $runNumber

# check if file exists
  if [ "$fileCheck" == true ]; then
      if ! [ -a ${path}/DQM_V0001_R000000001__Express__PCLTest__ALCAPROMPT.root ] 
	  then
	  echo file does not exist
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      qsub -cwd -l h_vmem=10G -m ae -q long.q -N gain${runNumber} $path/cmsDriver.sh
	  fi
	  continue
      fi
  fi

# check if file is smaller than 70 kb
  if [ "$sizeCheck" == true ]; then 

      actualsize=$(du -k "${path}/DQM_V0001_R000000001__Express__PCLTest__ALCAPROMPT.root" | cut -f 1)
      if [ $actualsize -le $minimumsize ]
	  then
	  echo size is below $minimumsize kilobytes
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      qsub -cwd -l h_vmem=10G -m ae -q long.q -N gain${runNumber} $path/cmsDriver.sh
	  fi
	  continue
      fi
  fi

# Check if file is a zombie
  if [ "$zombieCheck" == true ]; then
      root -b -l -q isZombieFile.C'("'${path}/DQM_V0001_R000000001__Express__PCLTest__ALCAPROMPT.root'")' 1> output.txt 2>error_f.txt 
      var=$( tail -1 output.txt )
      #echo $var
      if [ "$var" == "(int)1" ]
	  then 
	  echo "is a Zombie"
	  echo $path
	  if [ "$submit" == true ]; then
	      echo submit
	      qsub -cwd -l h_vmem=10G -m ae -q long.q -N gain${runNumber} $path/cmsDriver.sh
	  fi
      fi
  fi

done<allJobs.txt



rm output.txt
rm error_f.txt
rm isZombieFile_C.* 
rm allJobs.txt