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

for sample in DATASET:
   f = open('files' + str(sample[13:21]) + '.txt', 'w')
   results = commands.getstatusoutput('python2.6 das_client.py --query="file dataset='+sample+'" --limit=100000')[1].splitlines()
   print results
   
   runs = []
   results.sort()
   for line in results:
      run=line
      runs.append(run)
      
   for r in runs:
      f.write(str(r) + '\n')
  
