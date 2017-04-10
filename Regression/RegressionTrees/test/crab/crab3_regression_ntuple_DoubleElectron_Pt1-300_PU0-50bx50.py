from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'regression_ntuple_DoubleElectron_Pt1-300_PU0-50bx50'
config.General.workArea = 'crab3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runElectronRegressionTrees_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/DoubleElectron_FlatPt-1To300/RunIISpring15DR74-AsymptFlat0to50bx50RawReco_MCRUN2_74_V9A-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/jsauvan/Regression/Ntuples/'
config.Data.publication = False
config.Data.publishDataName = 'regression_ntuple_RunIISpring15DR74-AsymptFlat0to50bx50RawReco_MCRUN2_74_V9A-v3_2015_10_12'

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'
