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

   results = commands.getstatusoutput('python2.6 das_client.py --query="run dataset='+sample+' | grep run.run_number" --limit=100000')[1].splitlines()

   folder = sample[13:21]
   print folder
   commands.getstatusoutput('mkdir -p '+folder+'')

   runs = []
   results.sort()
   for line in results:
      run=int(line)
      runs.append(run)
      
   for r in runs:
      print str(r)
      if str(r) in open('runlistPromptReco.txt').read():
         print str(r)
      else:
         continue
      print r
      
      initEnv='#!/bin/sh\n'
      initEnv+='shopt -s expand_aliases;'
      initEnv+='module(){ eval `/usr/bin/modulecmd bash $*`;};'
      initEnv+='module use -a /afs/desy.de/group/cms/modulefiles/;'
      initEnv+='module load cmssw;'
      initEnv+='cd ' + os.getcwd() + '/'+folder+'/'+str(r)+'/;'
      initEnv+='eval `scramv1 runtime -sh`;'
      listFiles='python2.6 /nfs/dust/cms/user/tlenz/PixelCalibrationCode/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa/das_client.py --query="file dataset='+sample+' run='+str(r)+' " --limit=100000 > query.txt;'
      listFiles+='awk "{print $1}" < query.txt | paste -s -d,>files.txt;'
      listFiles+="var=$(cat files.txt);"
      
      commands.getstatusoutput('mkdir -p '+folder+'/'+str(r))
      print "submitting jobs for run " + str(r)
      config_file=open(''+folder+'/'+str(r)+'/cmsDriver.sh','w')
      config_file.write( initEnv + listFiles + 'cmsDriver.py '+folder+'' +str(r)+ ' --datatier ALCARECO --conditions auto:com10 -s ALCA:PromptCalibProdSiStripGains --eventcontent ALCARECO -n -1 --filein=${var} --fileout file:'+folder+''+str(r)+'_out.root' + '; rm core.*;')
      config_file.close()
      print('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + str(r) + ' ' + os.getcwd() + '/'+folder+'/'+str(r)+'/cmsDriver.sh')
      #out = commands.getstatusoutput('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + str(r) + ' ' + os.getcwd() + '/'+folder+'/'+str(r)+'/cmsDriver.sh')
      
