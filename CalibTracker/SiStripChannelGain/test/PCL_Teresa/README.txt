# steps to get calibration factors:


######
# for step3:
python mkListWithRunsWithNonZeroEntries.py
python LaunchWithRunlistPromptReco.py
source submitCrashedJobsAgain.sh   # change to right folder

# for step4:
submitToBatch.sh #(calls: source makeStep4Config\_AllSteps.sh)
source makeAllFits.sh

#####
# For MC
# for step3
python mkFileListMC.py
python LaunchWithRunlistMC.py
source submitCrashedJobsAgain.sh  # change to right folder

# for step4:
python submitToBatchMC.py
#add all DQM histograms