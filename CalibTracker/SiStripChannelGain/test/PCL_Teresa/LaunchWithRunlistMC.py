#!/usr/bin/env python

import urllib
import string
import os
import sys
import commands

folder='MC2012'
f=open('filesMC.txt','r').read().split('\n')
commands.getstatusoutput('mkdir -p '+folder+'')

filelist = []
a=0
for i in f:
      
   filelist.append(i)
   
   if len(filelist)<4 :
      continue
   
   myString = ",".join(filelist )
      
            
   initEnv='#!/bin/sh \n'
   initEnv+='shopt -s expand_aliases;'
   initEnv+='module(){ eval `/usr/bin/modulecmd bash $*`;};'
   initEnv+='module use -a /afs/desy.de/group/cms/modulefiles/;'
   initEnv+='module load cmssw;'
   initEnv+='cd ' + os.getcwd() + '/'+folder+'/'+str(a)+'/;'
   initEnv+='eval `scramv1 runtime -sh`;'
   listFiles='var="' + myString + '";'
   
   commands.getstatusoutput('mkdir -p '+folder+'/'+str(a))
   config_file=open(''+folder+'/'+str(a)+'/cmsDriver.sh','w')
   config_file.write( initEnv + listFiles + 'cmsDriver.py '+folder+'' +str(a)+ ' --no_exec --mc --datatier ALCARECO --conditions START53_V7G::All -s ALCA:PromptCalibProdSiStripGains --eventcontent ALCARECO -n -1 --filein=${var} --fileout file:'+folder+''+str(a)+'_out.root' + ';echo -e "process.pathALCARECOPromptCalibProdSiStripGains.remove(process.ALCARECOCalMinBiasFilterForSiStripGains)\nprocess.ALCARECOCalibrationTracks.src = \'generalTracks\'" >> '+folder+''+str(a)+'_ALCA.py;cmsRun '+folder+''+str(a)+'_ALCA.py;rm core.*;')
   config_file.close()
   commands.getstatusoutput('cd '+folder+'/'+str(a))
   print('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + folder + '_' + str(a) + ' ' + os.getcwd() + '/'+folder+'/'+str(a)+'/cmsDriver.sh')
   command_file=open(''+folder+'/'+str(a)+'/command.txt','w')
   command_file.write('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + folder + '_' + str(a) + ' ' + os.getcwd() + '/'+folder+'/'+str(a)+'/cmsDriver.sh')
   command_file.close()
   out = commands.getstatusoutput('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + folder + '_' + str(a) + ' ' + os.getcwd() + '/'+folder+'/'+str(a)+'/cmsDriver.sh')
   qsubCommand=open('qsubCommand.txt','w')
   qsubCommand.write('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + folder + '_' + str(a) + ' ' + os.getcwd() + '/'+folder+'/'+str(a)+'/cmsDriver.sh')
   qsubCommand.close()
   commands.getstatusoutput('cd ../..')
   filelist = []
   a+=1
