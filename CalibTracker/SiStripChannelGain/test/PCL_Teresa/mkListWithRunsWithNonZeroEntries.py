#!/usr/bin/env python

import urllib
import string
import os
import sys
import commands

DATASET = []
DATASETA="/MinimumBias/Run2012A-SiStripCalMinBias-v1/ALCARECO"
DATASET.append(DATASETA)
DATASETB="/MinimumBias/Run2012B-SiStripCalMinBias-v1/ALCARECO"
DATASET.append(DATASETB)
DATASETCv1="/MinimumBias/Run2012C-SiStripCalMinBias-v1/ALCARECO"
DATASET.append(DATASETCv1)
DATASETCv2="/MinimumBias/Run2012C-SiStripCalMinBias-v2/ALCARECO"
DATASET.append(DATASETCv2)
DATASETD="/MinimumBias/Run2012D-SiStripCalMinBias-v1/ALCARECO"
DATASET.append(DATASETD)
print DATASET

f = open('runlistPromptReco.txt', 'w')

for sample in DATASET:
   results = commands.getstatusoutput('python2.6 das_client.py --query="run dataset='+sample+' | grep run.run_number" --limit=100000')[1].splitlines()
   print results
   folder = sample[13:21]
   
   runs = []
   results.sort()
   for line in results:
      run=int(line)
      runs.append(run)

   commands.getstatusoutput('mkdir -p '+folder+'')
   for r in runs:
      print r
      out=commands.getstatusoutput('python2.6 das_client.py --query="summary dataset='+sample+' run='+str(r)+' | grep summary.nevents" --limit=100000')[1].splitlines()
      if int(out[0]) == 0:
         print "continued"
         continue;
      else:
         f.write(str(r) + '\n')
         print out[0]
         print r
