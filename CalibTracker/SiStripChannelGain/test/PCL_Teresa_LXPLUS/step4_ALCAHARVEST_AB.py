# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step4 --data --conditions auto:com10 --scenario pp -s ALCAHARVEST:SiStripGains --filein file:PromptCalibProdSiStripGains.root -n -1 --no_exec
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('filelistAB.txt')
readFiles = cms.untracked.vstring( *mylist)

process = cms.Process('ALCAHARVEST')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.AlCaHarvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = readFiles,
    processingMode = cms.untracked.string('RunsAndLumis')
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    fileMode = cms.untracked.string('FULLMERGE')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('step4 nevts:-1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

# Additional output definition

# Other statements
process.PoolDBOutputService.toPut.append(process.ALCAHARVESTSiStripGains_dbOutput)
process.pclMetadataWriter.recordsToMap.append(process.ALCAHARVESTSiStripGains_metadata)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10', '')

# Path and EndPath definitions
process.BeamSpotByRun = cms.Path(process.ALCAHARVESTBeamSpotByRun)
process.ALCAHARVESTDQMSaveAndMetadataWriter = cms.Path(process.dqmSaver+process.pclMetadataWriter)
process.SiStripGains = cms.Path(process.ALCAHARVESTSiStripGains)
process.BeamSpotByLumi = cms.Path(process.ALCAHARVESTBeamSpotByLumi)
process.SiStripQuality = cms.Path(process.ALCAHARVESTSiStripQuality)

# Schedule definition
process.schedule = cms.Schedule(process.SiStripGains,process.ALCAHARVESTDQMSaveAndMetadataWriter)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Gains_Tree_AB_Mean.root')
)
