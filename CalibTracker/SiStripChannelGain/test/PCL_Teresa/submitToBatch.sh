#!/bin/sh

cmsenv

SCRIPTPATH=`pwd -P`
echo ${SCRIPTPATH}

steps=("AB" "C1" "C2" "D1" "D2")


for i in "${steps[@]}"
  do
  echo $i
  cat makeStep4Config_AllSteps.sh | sed -e "s/\${step}/${i}/" > makeStep4Config_${i}.sh
  qsub -cwd -l h_vmem=10G -m ae -q long.q -N step${i} ${SCRIPTPATH}/makeStep4Config_${i}.sh
done

