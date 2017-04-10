import FWCore.ParameterSet.Config as cms

process = cms.Process("TUPLE")
process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

from Configuration.AlCa.autoCond import autoCond 
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.GlobalTag.globaltag = cms.string( "MCRUN2_74_V9A::All" )
print process.GlobalTag.globaltag

process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring("/store/mc/Fall14DR73/DoubleElectron_FlatPt-1To300/GEN-SIM-RECO/EGM_Flat0to50bx50_raw_MCRUN2_73_V15-v1/00000/00371A42-CEF3-E411-B666-002590D0B0D4.root")
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/mplaner/CMSSW_9_0_0_pre6/src/HLT_RECO_pat.root')
    )


# Event output
#process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("regression_ntuple_746.root")
    )

### ntuple makers
genParticleID = cms.int32(11)
matchMinERatio = cms.double(0.25)
matchMaxDR = cms.double(0.3)
matchingType = cms.int32(2)

process.mustacheSCTree = cms.EDAnalyzer(
    "PFSuperClusterRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    genParticleID = genParticleID,
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALBarrel'),
    superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALEndcapWithPreshower'),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )

process.gedGsfElectronTree = cms.EDAnalyzer(
    "GedGsfElectronRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    gsfElectronSrc = cms.InputTag('gedGsfElectrons', ''),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )


process.egSCTree = cms.EDAnalyzer(
    "EGSuperClusterRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    genParticleID = genParticleID,
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    superClusterSrcEB = cms.InputTag('correctedHybridSuperClusters', ''),
    superClusterSrcEE = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower', ''),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )

process.regSeq = cms.Sequence( #process.egSCTree +
                               process.gedGsfElectronTree +
                               process.mustacheSCTree )

process.regPath = cms.Path(process.regSeq)

process.schedule = cms.Schedule(process.regPath)
