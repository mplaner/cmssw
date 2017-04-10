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


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)

process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000000

process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring(
          $inputFileNames ),
)

#process.source = cms.Source(
#    "PoolSource",
#    fileNames = cms.untracked.vstring([
##'file:/afs/hep.wisc.edu/cms/tperry/lastJanv_slc6_493_CMSSW_7_4_16_patch2/src/HLTrigger/Configuration/test/HLT_RECO_pat.root',
#'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_pat_74x_25ns/160203_155009/0000/HLT_RECO_pat_10.root',
###'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_pat_74x_25ns/160203_155009/0000/HLT_RECO_pat_11.root',
###'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_pat_74x_25ns/160203_155009/0000/HLT_RECO_pat_12.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2001.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2002.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2003.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2004.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2005.root',
##'/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_74x_25ns_rh/160131_163316/0002/RelVal_HLT2_GRun_MC_2006.root',
#                                       ])
#    )

process.options= cms.untracked.PSet(
SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Event output
#process.load("Configuration.EventContent.EventContent_cff")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("HLT_RECO_SCUns_ntuple.root")
    #fileName = cms.string("HLT_RECO_SCL1s_ntuple.root")
    )
process.TFileService.fileName=cms.string("$outputFileName")


### ntuple makers
genParticleID = cms.int32(11)
matchMinERatio = cms.double(0.25)
matchMaxDR = cms.double(0.3)
matchingType = cms.int32(2)
#
#matchMinERatio = cms.double(0.)
#matchMaxDR = cms.double(0.3)
#matchingType = cms.int32(2)

process.RECO_mustacheSCTree = cms.EDAnalyzer(
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
    rhoSrc = cms.InputTag('fixedGridRhoFastjetAll'),
    )

process.HLT_mustacheSCTree = cms.EDAnalyzer(
    "PFSuperClusterRegTreeMaker",
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    genParticleID = genParticleID,
    matchMinERatio = matchMinERatio,
    matchMaxDR = matchMaxDR,
    matchingType = matchingType,
    #superClusterSrcEB = cms.InputTag('hltParticleFlowSuperClusterECALUnseeded', 'hltParticleFlowSuperClusterECALBarrel'),
    #superClusterSrcEE = cms.InputTag('hltParticleFlowSuperClusterECALUnseeded', 'hltParticleFlowSuperClusterECALEndcapWithPreshower'),
    superClusterSrcEB = cms.InputTag('hltParticleFlowSuperClusterECALL1Seeded', 'hltParticleFlowSuperClusterECALBarrel'),
    superClusterSrcEE = cms.InputTag('hltParticleFlowSuperClusterECALL1Seeded', 'hltParticleFlowSuperClusterECALEndcapWithPreshower'),
    ebReducedRecHitCollection = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
    eeReducedRecHitCollection = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
    VtxLabel = cms.InputTag("offlinePrimaryVertices"),
    rhoSrc = cms.InputTag('fixedGridRhoFastjetAll'),
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

process.regSeq = cms.Sequence( process.RECO_mustacheSCTree + process.HLT_mustacheSCTree )

process.regPath = cms.Path(process.regSeq)

process.schedule = cms.Schedule(process.regPath)
