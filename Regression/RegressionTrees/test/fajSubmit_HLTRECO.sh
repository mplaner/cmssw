

#farmoutAnalysisJobs $1 --input-files-per-job=100 --use-osg --resubmit-failed-jobs \
farmoutAnalysisJobs $1 --input-files-per-job=100 --use-osg  \
  --input-dir=root://cmsxrootd.hep.wisc.edu/store/user/tmperry/DoubleElectron_FlatPt-1To300/HLT_RECO_pat_74x_25ns/160203_155009/ \
  FullerMars_74x_25ns_L1seed_pat_ntuple $CMSSW_BASE $CMSSW_BASE/src/Regression/RegressionTrees/test/runElectronRegressionTrees_HLTRECO_cfg.py


  #FullerMars_74x_25ns_Unseed_pat_ntuple $CMSSW_BASE $CMSSW_BASE/src/Regression/RegressionTrees/test/runElectronRegressionTrees_HLTRECO_cfg.py
