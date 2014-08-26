#!/usr/bin/env python

import os
import shutil

for dirname in os.walk('MC2012'):
    print dirname[0]
    #print dirname[0]+'/step4_ALCAHARVEST_MC.py'
    shutil.copy2('step4_ALCAHARVEST_MC.py',dirname[0]+'/step4_ALCAHARVEST_MC.py')

    initEnv='#!/bin/sh \n'
    initEnv+='shopt -s expand_aliases;'
    initEnv+='module(){ eval `/usr/bin/modulecmd bash $*`;};'
    initEnv+='module use -a /afs/desy.de/group/cms/modulefiles/;'
    initEnv+='module load cmssw;'
    initEnv+='cd ' + os.getcwd() + '/'+dirname[0]+'/;'
    initEnv+='eval `scramv1 runtime -sh`;'
    config_file=open(''+dirname[0]+'/step4Driver.sh','w')
    config_file.write( initEnv + 'cmsRun step4_ALCAHARVEST_MC.py;')
    config_file.close()

    command_file=open(''+dirname[0]+'/'+'command.txt','w')
    command_file.write('qsub -cwd -l h_vmem=8G -m ae -q long.q -N step4' + dirname[0] + ' ' + os.getcwd() + '/'+dirname[0]+'/step4Driver.sh')
    command_file.close()
    print('qsub -cwd -l h_vmem=8G -m ae -q long.q -N step4' + dirname[0] + ' ' + os.getcwd() + '/'+dirname[0]+'/step4Driver.sh')
    #out = commands.getstatusoutput('qsub -cwd -l h_vmem=8G -m ae -q long.q -N gain' + folder + '_' + str(a) + ' ' + os.getcwd() + '/'+folder+'/'+str(a)+'/cmsDriver.sh')

