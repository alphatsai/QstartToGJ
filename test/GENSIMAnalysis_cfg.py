import FWCore.ParameterSet.Config as cms

process = cms.Process("GENSIMAnalysis")

#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'MCRUN2_74_V1::All'  #https://twiki.cern.ch/twiki/bin/viewauth/CMS/RelValGT


# Default Parameter options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('MaxEvents', -1,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.int,
	"Run events max"
	)
options.register('OutFilename', 'GENSIMAnalysis.root',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Output File name"
	)
options.register('InputGENSIM', True,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.bool,
	"Input with GEN-SIM"
	)
options.register('genInfoLabel', 'genParticles',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Generator level infomation"
	)
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.MaxEvents) )

### input
from inputFiles_cfi import * #FileNames 
process.source = cms.Source("PoolSource",
    #skipEvents = cms.untracked.uint32(0),
    #firstEvent = cms.untracked.uint32(1),
    #fileNames = cms.untracked.vstring(FileNames_QstarM500Test)
    #fileNames = cms.untracked.vstring(FileNames_QstarM500_100K)
    fileNames = cms.untracked.vstring(FileNames_QstarM1000_75K)
)

### output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.OutFilename) 
)

### Input parameters
process.GENSIMAnalysis = cms.EDAnalyzer('GENSIMAnalysis')
process.p = cms.Path(process.GENSIMAnalysis)
