#!/usr/bin/env python

import urllib
import string
import os
import sys
import commands

DATASET = []
DATASETA="/MinimumBias/Run2012A-SiStripCalMinBias-v1/ALCARECO"
#DATASET.append(DATASETA)
DATASETB="/MinimumBias/Run2012B-SiStripCalMinBias-v1/ALCARECO"
#DATASET.append(DATASETB)
DATASETC="/MinimumBias/Run2012C-SiStripCalMinBias-v2/ALCARECO"
DATASET.append(DATASETC)
DATASETD="/MinimumBias/Run2012D-SiStripCalMinBias-v1/ALCARECO"
#DATASET.append(DATASETD)

print DATASET

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
      if int(r)<203648:
         continue
      print r
      out=commands.getstatusoutput('python2.6 das_client.py --query="summary dataset='+sample+' run='+str(r)+' | grep summary.nevents" --limit=100000')[1].splitlines()
      if int(out[0]) == 0:
         print "continued"
         continue;
      else: 
         print out[0]
         print r
      initEnv=''
      initEnv+='cd ' + os.getcwd() + '/'+folder+'/'+str(r)+'/;'
      initEnv+='source /afs/cern.ch/cms/cmsset_default.sh' + ';'
      initEnv+='eval `scramv1 runtime -sh`;'
      listFiles='python2.6 /afs/cern.ch/work/t/telenz/PixelCalibration/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa/das_client.py --query="file dataset='+sample+' run='+str(r)+' " --limit=100000 > query.txt;'
      listFiles+='awk "{print $1}" < query.txt | paste -s -d,>files.txt;'
      listFiles+="var=$(cat files.txt);"
      
      commands.getstatusoutput('mkdir -p '+folder+'/'+str(r))
      print "submitting jobs for run " + str(r)
      config_file=open(''+folder+'/'+str(r)+'/cmsDriver.sh','w')
      config_file.write( initEnv + listFiles + 'cmsDriver.py '+folder+'' +str(r)+ ' --datatier ALCARECO --conditions auto:com10 -s ALCA:PromptCalibProdSiStripGains --eventcontent ALCARECO -n -1 --filein=${var} --fileout file:'+folder+''+str(r)+'_out.root' + '; rm core.*; rm *.txt')
      config_file.close()
      print('bsub -q 2nd -J gainPCLrun' + str(r) +' "sh ' +  os.getcwd() + '/'+folder+'/'+str(r)+'/cmsDriver.sh"')
      out = commands.getstatusoutput('bsub -q 2nd -J gainPCLrun' + str(r) +' "sh ' +  os.getcwd() + '/'+folder+'/'+str(r)+'/cmsDriver.sh"')


