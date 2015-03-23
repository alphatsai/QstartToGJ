import FWCore.ParameterSet.Config as cms

process = cms.Process("GENSIMAnalysis")

#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'MCRUN2_74_V1::All'  #https://twiki.cern.ch/twiki/bin/viewauth/CMS/RelValGT

# Default Parameter options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('maxEvts', -1,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.int,
	"Run events max"
	)
options.register('inputFileName', 'FileNames',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"input File name"
	)
options.register('outputFileName', 'GENSIMAnalysis',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Output File name"
	)
options.register('inputGENSIM', True,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.bool,
	"Input with GEN-SIM"
	)
options.register('genInfoLabel', 'genParticles',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Generator level infomation"
	)
options.register('selectQstarStatus', -1,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.int,
	"Select Qstar status"
	)
options.register('numEventListsPrint', 2,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.int,
	"Number of event list be printed out"
	)
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )

### input
from inputFiles_cfi import * #FileNames 

print 'NOTE: Using filename from inputFiles.py: '+options.inputFileName
if options.inputFileName == 'FileNames':
	InputFiles = FileNames
elif options.inputFileName == 'FileNames_QstarM500_100K':
	InputFiles = FileNames_QstarM500_100K
elif options.inputFileName == 'FileNames_QstarM500_TSchanel_100K':
	InputFiles = FileNames_QstarM500_TSchanel_100K

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(InputFiles)
    #skipEvents = cms.untracked.uint32(0),
    #firstEvent = cms.untracked.uint32(1),
)

### output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string((options.outputFileName+'.root')) 
)

### Input parameters
process.GENSIMAnalysis = cms.EDAnalyzer('GENSIMAnalysis',
	genInfoLabel       = cms.InputTag(options.genInfoLabel),
	outputFileName	   = cms.string(options.outputFileName),
	selectQstarStatus  = cms.int32(options.selectQstarStatus),
	numEventListsPrint = cms.int32(options.numEventListsPrint),
)

process.p = cms.Path(process.GENSIMAnalysis)
