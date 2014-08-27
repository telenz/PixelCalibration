#!/bin/sh

cmsenv

root -b -l -q addDQMFiles_MC.C+ 1> output_1.txt 2>error_1.txt 

#rm output_1.txt
#rm error_1.txt
rm addDQMFiles_MC_*

