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
    fileNames = cms.untracked.vstring([#"/store/mc/RunIISpring15DR74/DoubleElectron_FlatPt-1To300/GEN-SIM-RECO/AsymptNoPURawReco_MCRUN2_74_V9A-v1/00000/00B9171B-3511-E511-B4AC-002590D9D956.root",
                                       "/store/mc/RunIISpring15DR74/DoubleElectron_FlatPt-1To300/GEN-SIM-RECO/AsymptNoPURawReco_MCRUN2_74_V9A-v1/00000/066C6E62-3811-E511-BF46-008CFA197C58.root",
                                       #"/store/mc/RunIISpring15DR74/DoubleElectron_FlatPt-1To300/GEN-SIM-RECO/AsymptNoPURawReco_MCRUN2_74_V9A-v1/00000/10DF6151-3311-E511-B266-0025905A60E4.root",
                                       #"/store/mc/RunIISpring15DR74/DoubleElectron_FlatPt-1To300/GEN-SIM-RECO/AsymptNoPURawReco_MCRUN2_74_V9A-v1/00000/124BD3E8-3E11-E511-9B53-002590D9D8A4.root"
                                       ])
    )


# Event output
#process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("test.root")
    )

### ntuple makers
genParticleID = cms.int32(11)
matchMinERatio = cms.double(0.25)
matchMaxDR = cms.double(0.3)
matchingType = cms.int32(2)
#
#matchMinERatio = cms.double(0.)
#matchMaxDR = cms.double(0.3)
#matchingType = cms.int32(2)

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
    #superClusterSrcEB = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALBarrel','reRECO'),
    #superClusterSrcEE = cms.InputTag('particleFlowSuperClusterECAL', 'particleFlowSuperClusterECALEndcapWithPreshower','reRECO'),
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
    superClusterSrcEB = cms.InputTag('correctedHybridSuperClusters', '','reRECO'),
    superClusterSrcEE = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower', '','reRECO'),
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    )

process.regSeq = cms.Sequence( #process.egSCTree +
                               process.gedGsfElectronTree +
                               process.mustacheSCTree )

process.regPath = cms.Path(process.regSeq)

process.schedule = cms.Schedule(process.regPath)
