#!/bin/bash/

cmsenv

SCRIPTPATH=`pwd -P`

echo ${SCRIPTPATH}

bsub -q 2nd -J step4AB "sh ${SCRIPTPATH}/makeStep4Config_AB.sh"
bsub -q 1nd -J step4C1 "sh ${SCRIPTPATH}/makeStep4Config_C1.sh"
bsub -q 2nd -J step4C2 "sh ${SCRIPTPATH}/makeStep4Config_C2.sh"
bsub -q 2nd -J step4D1 "sh ${SCRIPTPATH}/makeStep4Config_D1.sh"
bsub -q 2nd -J step4D2 "sh ${SCRIPTPATH}/makeStep4Config_D2.sh"