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
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/user/lgray/SinglePhotonGunWithPU_700pre12reRECO-step4_RECO_EI/step4_RECO_EI-000A6123-D369-E211-B8EF-003048F9EB46.root')
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/mplaner/CMSSW_9_0_0_pre6/src/relval.root')
     fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/mplaner/CMSSW_9_0_0_pre6/src/HLT_RECO_pat.root')
    )


# Event output
#process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('test_test.root')
    )

### ntuple makers
genParticleID = cms.int32(22)
matchMinERatio = cms.double(0.25)
matchMaxDR = cms.double(0.3)
matchingType = cms.int32(2)

process.mustacheSCTree = cms.EDAnalyzer(
    "PFSuperClusterRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    rhoSrc = cms.InputTag('hltKT6CaloJets','rho'),
    genParticleID = genParticleID,
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    #superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECALMustache', 'particleFlowSuperClusterECALBarrel'),
    #superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECALMustache', 'particleFlowSuperClusterECALEndcapWithPreshower'),
    #superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALBarrel','reRECO'),
    #superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALEndcapWithPreshower','reRECO'),
    superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALBarrel'),
    superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALEndcapWithPreshower'),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )

process.gedPhotonTree = cms.EDAnalyzer(
    "GedPhotonRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    photonSrc = cms.InputTag('gedPhotons', ''),
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
    #superClusterSrcEB = cms.InputTag('correctedHybridSuperClusters', '','reRECO'),
    #superClusterSrcEE = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower', '','reRECO'),
    superClusterSrcEB = cms.InputTag('correctedHybridSuperClusters', ''),
    superClusterSrcEE = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower', ''),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )

process.regSeq = cms.Sequence( process.egSCTree +
                               process.gedPhotonTree +
                               process.mustacheSCTree )

process.regPath = cms.Path(process.regSeq)

process.schedule = cms.Schedule(process.regPath)
