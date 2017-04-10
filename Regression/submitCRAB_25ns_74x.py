from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'HLT_RECO_pat_74x_25ns'
config.General.workArea = 'CRABLOGS'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HLT_RECO_pat_74x_25ns.py'
#config.JobType.psetName = 'fullconfig_HLT_RECO.py'

config.Data.inputDataset = '/DoubleElectron_FlatPt-1To300/RunIISpring15DR74-AsymptFlat0to50bx25RawReco_MCRUN2_74_V9-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.totalUnits = 10
config.Data.useParent = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HLT_RECO_pat_74x_25ns'

#config.Site.blacklist = 'T2_EE_'
config.Site.whitelist = ['T2_US_*']
config.Site.storageSite = 'T2_US_Wisconsin'
