#!/usr/bin/env python

import urllib
import string
import os
import sys
import commands

DATASET = []
DATASETMC="/MinBias_TuneZ2star_8TeV-pythia6/Summer12_DR53X-DEBUG_PU_S10_START53_V7A-v1/GEN-SIM-RECODEBUG"
DATASET.append(DATASETMC)

print DATASET

for sample in DATASET:
   f = open('filesMC.txt', 'w')
   results = commands.getstatusoutput('python2.6 das_client.py --query="file dataset='+sample+'" --limit=100000')[1].splitlines()
   print results
   
   runs = []
   results.sort()
   for line in results:
      run=line
      runs.append(run)
      
   for r in runs:
      f.write(str(r) + '\n')
  
