#!bin/bash

steps=("AB" "C1" "C2" "D1" "D2" "mc")


root -b -l <<EOF
.L makeFitsForCalibration.C+
.q
EOF

for i in "${steps[@]}"
do
  echo ${i}
  root -l -b -q makeFitsForCalibration.C+'("'${i}'")' 1> output_${i}.txt 2>error_f_${i}.txt &
done
rm makFitsForCalibration.C~
rm makeFitsForCalibration_C.*



