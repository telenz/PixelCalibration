import FWCore.ParameterSet.Config as cms

# FIXME: the safest option would be to import the basic cfi from a place mantained by the developer
alcaSiStripGainsHarvester =  cms.EDAnalyzer("SiStripGainFromCalibTree",

                                  OutputGains         = cms.string('Gains_ASCII.txt'),
                                  Tracks              = cms.untracked.InputTag('ALCARECOCalibrationTracksRefit'),
                                  AlgoMode            = cms.untracked.string('PCL'),

                                  #Gain quality cuts
                                  minNrEntries        = cms.untracked.double(25),
                                  maxChi2OverNDF      = cms.untracked.double(9999999.0),
                                  maxMPVError         = cms.untracked.double(25.0),

                                  #track/cluster quality cuts
                                  minTrackMomentum    = cms.untracked.double(2),
                                  maxNrStrips         = cms.untracked.uint32(8),

                                  Validation          = cms.untracked.bool(False),
                                  OldGainRemoving     = cms.untracked.bool(False),
                                  FirstSetOfConstants = cms.untracked.bool(True),

                                  CalibrationLevel    = cms.untracked.int32(2), # 0==APV, 1==Laser, 2==module

                                  InputFiles          = cms.vstring(),

                                  UseCalibration     = cms.untracked.bool(False),
                                  calibrationPath    = cms.untracked.string(""),

                                  SinceAppendMode     = cms.bool(True),
                                  IOVMode             = cms.string('Job'),
                                  Record              = cms.string('SiStripApvGainRcd'),
                                  doStoreOnDB         = cms.bool(True),
                                  )

alcaSiStripGainsHarvester.FirstSetOfConstants = cms.untracked.bool(False)
alcaSiStripGainsHarvester.CalibrationLevel    = cms.untracked.int32(2) # 0==APV, 1==Laser, 2==module

# THIS is the crucial parameter to use this in harvesting mode                                  
alcaSiStripGainsHarvester.harvestingMode      = cms.untracked.bool(True)
#alcaSiStripGainsHarvester.UseCalibration      = cms.untracked.bool(False)
#alcaSiStripGainsHarvester.calibrationPath     = cms.untracked.string("/afs/cern.ch/user/t/telenz/work/PixelCalibration/CMSSW_5_3_15/src/CalibTracker/SiStripChannelGain/test/PCL_Teresa/Gains_Tree_C1_CL2.root")

