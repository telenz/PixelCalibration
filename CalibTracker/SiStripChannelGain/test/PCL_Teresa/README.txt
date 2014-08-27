# steps to get calibration factors:


###### For Data ######
# for step3:
python mkListWithRunsWithNonZeroEntries.py # makes a list with all runs during 2012 with nonzero entries
python LaunchWithRunlistPromptReco.py      # Run step3 on batch
source submitCrashedJobsAgain.sh           # change to right folder

# for step4:
submitToBatch.sh #(calls: source makeStep4Config\_AllSteps.sh)  # Runs step4 on batch
######

###### For MC ######
# for step3
python mkFileListMC.py                    # makes MC file list output: filesMC.txt
python LaunchWithRunlistMC.py             # Run step3 on batch 
source submitCrashedJobsAgain.sh          # change to right folder

# for step4:
python submitToBatchMC.py                 # Run step4 on batch
source submitCrashedJobsFromStep4Again.sh # checks whether all DQM files are ok      
source addDQMFiles_MC.sh                  # add all DQM histograms
######

##### Fits #####
source makeAllFits.sh                     # make all Fits
makeCalibrationTreeFromFitsMC.C           # make Calibration Tree for MC
makeCalibrationTreeFromFitsData.C         # make Calibration Tree for Data
#####