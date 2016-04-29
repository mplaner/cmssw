import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            '/store/relval/CMSSW_7_1_0_pre1/RelValProdTTbar/GEN-SIM-RECO/START70_V5-v1/00000/14842A6B-2086-E311-B5CB-02163E00E8DA.root'
                )
                            )
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('testTFileService.root')
                                   )


process.testAna = cms.EDAnalyzer("TestTFileServiceAnalyzer",
                                 dir1 = cms.string("mydir1"),
                                 dir2 = cms.string("mydir2")
                                 )

process.p = cms.Path( process.testAna )
