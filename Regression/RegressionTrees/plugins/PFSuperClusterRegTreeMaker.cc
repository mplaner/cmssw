//
// Class: PFSuperClusterRegTreeMaker.cc
//
// Info: Processes a track into histograms of delta-phis and such
//
// Author: L. Gray (FNAL)
//


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"

#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

#include "Regression/RegressionTrees/plugins/PFSuperClusterRegTreeMaker.h"

#include "TTree.h"
#include "TVector2.h"

#include <map>


namespace MK = reco::MustacheKernel;

PFSuperClusterRegTreeMaker::PFSuperClusterRegTreeMaker(const PSet& p) 
{
    _calib.reset(new PFEnergyCalibration());
    _N_ECALClusters   = 1;
    _N_PSClusters     = 1;
    _tree = _fs->make<TTree>("SuperClusterTree","Dump of all available SC info");
    _tree->Branch("eventNumber"                 , &_eventNumber                , "eventNumber/I");
    _tree->Branch("N_ECALClusters"              , &_N_ECALClusters             , "N_ECALClusters/I");
    _tree->Branch("N_RecHits_ee"                , &_N_RecHits_ee               , "N_RecHits_ee/I");
    _tree->Branch("N_RecHits_eb"                , &_N_RecHits_eb               , "N_RecHits_eb/I");
    _tree->Branch("N_RecHits"                   , &_N_RecHits                  , "N_RecHits/I");
    _tree->Branch("N_RecHits_ee_p5"             , &_N_RecHits_ee_p5            , "N_RecHits_ee_p5/I");
    _tree->Branch("N_RecHits_eb_p5"             , &_N_RecHits_eb_p5            , "N_RecHits_eb_p5/I");
    _tree->Branch("N_RecHits_p5"                , &_N_RecHits_p5               , "N_RecHits_p5/I");
    _tree->Branch("N_RecHits_ee_1"              , &_N_RecHits_ee_1             , "N_RecHits_ee_1/I");
    _tree->Branch("N_RecHits_eb_1"              , &_N_RecHits_eb_1             , "N_RecHits_eb_1/I");
    _tree->Branch("N_RecHits_1"                 , &_N_RecHits_1                , "N_RecHits_1/I");
    _tree->Branch("N_RecHits_ee_1p5"            , &_N_RecHits_ee_1p5           , "N_RecHits_ee_1p5/I");
    _tree->Branch("N_RecHits_eb_1p5"            , &_N_RecHits_eb_1p5           , "N_RecHits_eb_1p5/I");
    _tree->Branch("N_RecHits_1p5"               , &_N_RecHits_1p5              , "N_RecHits_1p5/I");
    _tree->Branch("N_RecHits_ee_2"              , &_N_RecHits_ee_2             , "N_RecHits_ee_2/I");
    _tree->Branch("N_RecHits_eb_2"              , &_N_RecHits_eb_2             , "N_RecHits_eb_2/I");
    _tree->Branch("N_RecHits_2"                 , &_N_RecHits_2                , "N_RecHits_2/I");
    _tree->Branch("N_RecHits_ee_2p5"            , &_N_RecHits_ee_2p5           , "N_RecHits_ee_2p5/I");
    _tree->Branch("N_RecHits_eb_2p5"            , &_N_RecHits_eb_2p5           , "N_RecHits_eb_2p5/I");
    _tree->Branch("N_RecHits_2p5"               , &_N_RecHits_2p5              , "N_RecHits_2p5/I");
    _tree->Branch("N_RecHits_ee_3"              , &_N_RecHits_ee_3             , "N_RecHits_ee_3/I");
    _tree->Branch("N_RecHits_eb_3"              , &_N_RecHits_eb_3             , "N_RecHits_eb_3/I");
    _tree->Branch("N_RecHits_3"                 , &_N_RecHits_3                , "N_RecHits_3/I");
    _tree->Branch("N_RecHits_ee_3p5"            , &_N_RecHits_ee_3p5           , "N_RecHits_ee_3p5/I");
    _tree->Branch("N_RecHits_eb_3p5"            , &_N_RecHits_eb_3p5           , "N_RecHits_eb_3p5/I");
    _tree->Branch("N_RecHits_3p5"               , &_N_RecHits_3p5              , "N_RecHits_3p5/I");
    _tree->Branch("nVtx"                        , &_nVtx                       , "nVtx/I");
    _tree->Branch("rho"                         , &_rho                        , "rho/F");
    _tree->Branch("N_PSClusters"                , &_N_PSClusters               , "N_PSClusters/I");
    _tree->Branch("scIndex"                     , &_scIndex                    , "scIndex/I");
    _tree->Branch("scRawEnergy"                 , &_scRawEnergy                , "scRawEnergy/F");
    _tree->Branch("scCalibratedEnergy"          , &_scCalibratedEnergy         , "scCalibratedEnergy/F");
    _tree->Branch("scPreshowerEnergy"           , &_scPreshowerEnergy          , "scPreshowerEnergy/F");
    _tree->Branch("scPreshowerEnergyPlane1"     , &_scPreshowerEnergyPlane1    , "scPreshowerEnergyPlane1/F");
    _tree->Branch("scPreshowerEnergyPlane2"     , &_scPreshowerEnergyPlane2    , "scPreshowerEnergyPlane2/F");
    _tree->Branch("scIsEB"                      , &_scIsEB                     , "scIsEB/I");
    _tree->Branch("scEta"                       , &_scEta                      , "scEta/F");
    _tree->Branch("scPhi"                       , &_scPhi                      , "scPhi/F");
    _tree->Branch("scR"                         , &_scR                        , "scR/F");
    _tree->Branch("scPhiWidth"                  , &_scPhiWidth                 , "scPhiWidth/F");
    _tree->Branch("scEtaWidth"                  , &_scEtaWidth                 , "scEtaWidth/F");
    _tree->Branch("scSeedRawEnergy"             , &_scSeedRawEnergy            , "scSeedRawEnergy/F");
    _tree->Branch("scSeedCalibratedEnergy"      , &_scSeedCalibratedEnergy     , "scSeedCalibratedEnergy/F");
    _tree->Branch("scSeedEta"                   , &_scSeedEta                  , "scSeedEta/F");
    _tree->Branch("scSeedPhi"                   , &_scSeedPhi                  , "scSeedPhi/F");
    _tree->Branch("scSeedSize"                  , &_scSeedSize                 , "scSeedSize/I");
    _tree->Branch("scSeedR9"                    , &_scSeedR9                   , "scSeedR9/F");
    _tree->Branch("scSeedEmax"                  , &_scSeedEmax                 , "scSeedEmax/F");
    _tree->Branch("scSeedE2nd"                  , &_scSeedE2nd                 , "scSeedE2nd/F");
    _tree->Branch("scSeedLeftRightAsym"         , &_scSeedLeftRightAsym        , "scSeedLeftRightAsym/F");
    _tree->Branch("scSeedTopBottomAsym"         , &_scSeedTopBottomAsym        , "scSeedTopBottomAsym/F");
    _tree->Branch("scSeedE2x5max"               , &_scSeedE2x5max              , "scSeedE2x5max/F");
    _tree->Branch("scSeed2x5LeftRightAsym"      , &_scSeed2x5LeftRightAsym     , "scSeed2x5LeftRightAsym/F");
    _tree->Branch("scSeed2x5TopBottomAsym"      , &_scSeed2x5TopBottomAsym     , "scSeed2x5TopBottomAsym/F");
    _tree->Branch("scSeedSigmaIetaIeta"         , &_scSeedSigmaIetaIeta        , "scSeedSigmaIetaIeta/F");
    _tree->Branch("scSeedSigmaIetaIphi"         , &_scSeedSigmaIetaIphi        , "scSeedSigmaIetaIphi/F");
    _tree->Branch("scSeedSigmaIphiIphi"         , &_scSeedSigmaIphiIphi        , "scSeedSigmaIphiIphi/F");
    _tree->Branch("scSeedR9NoZS"                , &_scSeedR9NoZS               , "scSeedR9NoZS/F");
    _tree->Branch("scSeedEmaxNoZS"              , &_scSeedEmaxNoZS             , "scSeedEmaxNoZS/F");
    _tree->Branch("scSeedE2ndNoZS"              , &_scSeedE2ndNoZS             , "scSeedE2ndNoZS/F");
    _tree->Branch("scSeedLeftRightAsymNoZS"     , &_scSeedLeftRightAsymNoZS    , "scSeedLeftRightAsymNoZS/F");
    _tree->Branch("scSeedTopBottomAsymNoZS"     , &_scSeedTopBottomAsymNoZS    , "scSeedTopBottomAsymNoZS/F");
    _tree->Branch("scSeedE2x5maxNoZS"           , &_scSeedE2x5maxNoZS          , "scSeedE2x5maxNoZS/F");
    _tree->Branch("scSeed2x5LeftRightAsymNoZS"  , &_scSeed2x5LeftRightAsymNoZS , "scSeed2x5LeftRightAsymNoZS/F");
    _tree->Branch("scSeed2x5TopBottomAsymNoZS"  , &_scSeed2x5TopBottomAsymNoZS , "scSeed2x5TopBottomAsymNoZS/F");
    _tree->Branch("scSeedSigmaIetaIetaNoZS"     , &_scSeedSigmaIetaIetaNoZS    , "scSeedSigmaIetaIetaNoZS/F");
    _tree->Branch("scSeedSigmaIetaIphiNoZS"     , &_scSeedSigmaIetaIphiNoZS    , "scSeedSigmaIetaIphiNoZS/F");
    _tree->Branch("scSeedSigmaIphiIphiNoZS"     , &_scSeedSigmaIphiIphiNoZS    , "scSeedSigmaIphiIphiNoZS/F");
    _tree->Branch("scSeedCryEta"                , &_scSeedCryEta               , "scSeedCryEta/F");
    _tree->Branch("scSeedCryPhi"                , &_scSeedCryPhi               , "scSeedCryPhi/F");
    _tree->Branch("scSeedCryIeta"               , &_scSeedCryIeta              , "scSeedCryIeta/I");
    _tree->Branch("scSeedCryIphi"               , &_scSeedCryIphi              , "scSeedCryIphi/I");
    _tree->Branch("scSeedCryIetaV2"             , &_scSeedCryIetaV2            , "scSeedCryIetaV2/I");
    _tree->Branch("scSeedCryIphiV2"             , &_scSeedCryIphiV2            , "scSeedCryIphiV2/I");
    _tree->Branch("scSeedCryX"                  , &_scSeedCryX                 , "scSeedCryX/F");
    _tree->Branch("scSeedCryY"                  , &_scSeedCryY                 , "scSeedCryY/F");
    _tree->Branch("scSeedCryIx"                 , &_scSeedCryIx                , "scSeedCryIx/I");
    _tree->Branch("scSeedCryIy"                 , &_scSeedCryIy                , "scSeedCryIy/I");
    _tree->Branch("scSeedCryIxV2"               , &_scSeedCryIxV2              , "scSeedCryIxV2/I");
    _tree->Branch("scSeedCryIyV2"               , &_scSeedCryIyV2              , "scSeedCryIyV2/I");

    // ecal cluster information
    _clusterRawEnergy     .reset(new float[1], array_deleter<float>());
    _clusterCalibEnergy   .reset(new float[1], array_deleter<float>());
    _clusterEta           .reset(new float[1], array_deleter<float>());
    _clusterPhi           .reset(new float[1], array_deleter<float>());
    _clusterDEtaToSeed    .reset(new float[1], array_deleter<float>());  
    _clusterDPhiToSeed    .reset(new float[1], array_deleter<float>());
    _clusterDPhiToCentroid.reset(new float[1], array_deleter<float>());
    _clusterDEtaToCentroid.reset(new float[1], array_deleter<float>());
    _clusterInMustache    .reset(new int[1],   array_deleter<int>());
    _clusterInDynDPhi     .reset(new int[1],   array_deleter<int>());
    _clusterLeakage       .reset(new int[1],   array_deleter<int>());
    _tree->Branch("clusterRawEnergy"      , _clusterRawEnergy.get()     , "clusterRawEnergy[N_ECALClusters]/F");
    _tree->Branch("clusterRawEnergy_0"    , &_clusterRawEnergy_0        , "clusterRawEnergy_0/F");
    _tree->Branch("clusterRawEnergy_1"    , &_clusterRawEnergy_1        , "clusterRawEnergy_1/F");
    _tree->Branch("clusterRawEnergy_2"    , &_clusterRawEnergy_2        , "clusterRawEnergy_2/F");
    _tree->Branch("clusterCalibEnergy"    , _clusterCalibEnergy.get()   , "clusterCalibEnergy[N_ECALClusters]/F");
    _tree->Branch("clusterEta"            , _clusterEta.get()           , "clusterEta[N_ECALClusters]/F");
    _tree->Branch("clusterPhi"            , _clusterPhi.get()           , "clusterPhi[N_ECALClusters]/F");
    _tree->Branch("clusterDPhiToSeed"     , _clusterDPhiToSeed.get()    , "clusterDPhiToSeed[N_ECALClusters]/F");
    _tree->Branch("clusterDPhiToSeed_0"   , &_clusterDPhiToSeed_0       , "clusterDPhiToSeed_0/F");
    _tree->Branch("clusterDPhiToSeed_1"   , &_clusterDPhiToSeed_1       , "clusterDPhiToSeed_1/F");
    _tree->Branch("clusterDPhiToSeed_2"   , &_clusterDPhiToSeed_2       , "clusterDPhiToSeed_2/F");
    _tree->Branch("clusterDEtaToSeed"     , _clusterDEtaToSeed.get()    , "clusterDEtaToSeed[N_ECALClusters]/F");  
    _tree->Branch("clusterDEtaToSeed_0"   , &_clusterDEtaToSeed_0       , "clusterDEtaToSeed_0/F");
    _tree->Branch("clusterDEtaToSeed_1"   , &_clusterDEtaToSeed_1       , "clusterDEtaToSeed_1/F");
    _tree->Branch("clusterDEtaToSeed_2"   , &_clusterDEtaToSeed_2       , "clusterDEtaToSeed_2/F");
    _tree->Branch("clusterDPhiToCentroid" , _clusterDPhiToCentroid.get(), "clusterDPhiToCentroid[N_ECALClusters]/F");
    _tree->Branch("clusterDEtaToCentroid" , _clusterDEtaToCentroid.get(), "clusterDEtaToCentroid[N_ECALClusters]/F");
    _tree->Branch("clusterInMustache"     , _clusterInMustache.get()    , "clusterInMustache[N_ECALClusters]/I");
    _tree->Branch("clusterInDynDPhi"      , _clusterInDynDPhi.get()     , "clusterInDynDPhi[N_ECALClusters]/I");
    _tree->Branch("clusterLeakage"        , _clusterLeakage.get()       , "clusterLeakage[N_ECALClusters]/I");

    _tree->Branch("clusterMaxDR"          , &_clusterMaxDR              , "clusterMaxDR/F");
    _tree->Branch("clusterMaxDRDPhi"      , &_clusterMaxDRDPhi          , "clusterMaxDRDPhi/F");
    _tree->Branch("clusterMaxDRDEta"      , &_clusterMaxDRDEta          , "clusterMaxDRDEta/F");
    _tree->Branch("clusterMaxDRRawEnergy" , &_clusterMaxDRRawEnergy     , "clusterMaxDRRawEnergy/F");

    _tree->Branch("clustersMeanRawEnergy" , &_clustersMeanRawEnergy     , "clustersMeanRawEnergy/F");
    _tree->Branch("clustersRMSRawEnergy"  , &_clustersRMSRawEnergy      , "clustersRMSRawEnergy/F");
    _tree->Branch("clustersMeanDRToSeed"  , &_clustersMeanDRToSeed      , "clustersMeanDRToSeed/F");
    _tree->Branch("clustersMeanDEtaToSeed", &_clustersMeanDEtaToSeed    , "clustersMeanDEtaToSeed/F");
    _tree->Branch("clustersMeanDPhiToSeed", &_clustersMeanDPhiToSeed    , "clustersMeanDPhiToSeed/F");

    // preshower information
    _psClusterRawEnergy.reset(new float[1], array_deleter<float>());
    _psClusterEta      .reset(new float[1], array_deleter<float>());
    _psClusterPhi      .reset(new float[1], array_deleter<float>());
    _tree->Branch("psClusterRawEnergy"    , _psClusterRawEnergy.get()   , "psClusterRawEnergy[N_PSClusters]/F");
    _tree->Branch("psClusterEta"          , _psClusterEta.get()         , "psClusterEta[N_PSClusters]/F");
    _tree->Branch("psClusterPhi"          , _psClusterPhi.get()         , "psClusterPhi[N_PSClusters]/F");


    _tree->Branch("isMatched"      , &_isMatched      , "isMatched/I");
    _tree->Branch("genPt"          , &_genPt          , "genPt/F");
    _tree->Branch("genEta"         , &_genEta         , "genEta/F");
    _tree->Branch("genPhi"         , &_genPhi         , "genPhi/F");
    _tree->Branch("genEnergy"      , &_genEnergy      , "genEnergy/F");
    _tree->Branch("genDEoE"        , &_genDEoE        , "genDEoE/F");
    _tree->Branch("genDRToCentroid", &_genDRToCentroid, "genDRToCentroid/F");
    _tree->Branch("genDRToSeed"    , &_genDRToSeed    , "genDRToSeed/F");

    _clusterDPhiToGen.reset(new float[1],array_deleter<float>());
    _clusterDEtaToGen.reset(new float[1],array_deleter<float>());
    _tree->Branch("clusterDPhiToGen", _clusterDPhiToGen.get(), "clusterDPhiToGen[N_ECALClusters]/F");
    _tree->Branch("clusterDEtaToGen", _clusterDEtaToGen.get(), "clusterDPhiToGen[N_ECALClusters]/F");


//    _scInputEB                 = p.getParameter<edm::InputTag>("superClusterSrcEB");
//    _scInputEE                 = p.getParameter<edm::InputTag>("superClusterSrcEE"); 
//    _ebReducedRecHitToken      = consumes<EcalRecHitCollection>(p.getParameter<edm::InputTag>("ebReducedRecHitCollection"));
//    _eeReducedRecHitToken      = consumes<EcalRecHitCollection>(p.getParameter<edm::InputTag>("eeReducedRecHitCollection"));
//    _PrimaryVertex             = p.getParameter<edm::InputTag>("VtxLabel"); 
//    _geninput                  = p.getParameter<edm::InputTag>("genSrc");    
    _rhoToken                    = consumes<double>(p.getParameter<edm::InputTag>("rhoSrc"));
    _scInputEBToken              = consumes<std::vector<reco::SuperCluster>>(p.getParameter<edm::InputTag>("superClusterSrcEB"));
    _scInputEEToken              = consumes<std::vector<reco::SuperCluster>>(p.getParameter<edm::InputTag>("superClusterSrcEE")); 
    _ebReducedRecHitToken        = consumes<EcalRecHitCollection>(p.getParameter<edm::InputTag>("ebReducedRecHitCollection"));
    _eeReducedRecHitToken        = consumes<EcalRecHitCollection>(p.getParameter<edm::InputTag>("eeReducedRecHitCollection"));
    _PrimaryVertexToken          = consumes<std::vector<reco::Vertex>>(p.getParameter<edm::InputTag>("VtxLabel")); 
    _geninputToken               = consumes<std::vector<reco::GenParticle>>(p.getParameter<edm::InputTag>("genSrc"));   
    _genParticleID             = p.getParameter<int>("genParticleID");
    _matchMinERatio            = p.getParameter<double>("matchMinERatio");
    _matchMaxDR                = p.getParameter<double>("matchMaxDR");
    _matchingType              = p.getParameter<int>("matchingType");// 0=DR; 1=DE/E, 2=DR+DE/E

}


void PFSuperClusterRegTreeMaker::analyze(const edm::Event& e, const edm::EventSetup& es) 
{

    //_lazyTool = new EcalClusterLazyTools(e, es, _ebReducedRecHitCollection, _eeReducedRecHitCollection); 
    _lazyTool = new EcalClusterLazyTools(e, es, _ebReducedRecHitToken, _eeReducedRecHitToken); 
    _lazyToolNoZS = new noZS::EcalClusterLazyTools(e, es, _ebReducedRecHitToken, _eeReducedRecHitToken);

    edm::Handle<reco::SuperClusterCollection> ebSCs, eeSCs;
    edm::Handle<reco::VertexCollection> PV;
    e.getByToken(_scInputEBToken, ebSCs);  
    e.getByToken(_scInputEEToken, eeSCs);  
    e.getByToken(_PrimaryVertexToken, PV);
//    e.getByLabel(_scInputEB, ebSCs);  
//    e.getByLabel(_scInputEE, eeSCs);  
//    e.getByLabel(_PrimaryVertex, PV);
    _eventNumber = e.id().event();
    _nVtx = PV->size();
    _scIndex = 0;
    _matchedParticles.clear();
    _matchedSuperClusters.clear();

    //Get Rho
    edm::Handle<double> hRhoKt6PFJets;
    e.getByToken(_rhoToken, hRhoKt6PFJets);
//    e.getByLabel(edm::InputTag("fixedGridRhoFastjetAll",""), hRhoKt6PFJets);
    _rho = (*hRhoKt6PFJets);

    // nRecHit
    //const int N_ECAL = sc.clustersEnd() - sc.clustersBegin();

  //ev.getByToken( esRHToken_, pESRecHits );
  //esRecHits_ = pESRecHits.product();
  //// make the map of rechits
  //rechits_map_.clear();
  //if (pESRecHits.isValid()) {
  //  EcalRecHitCollection::const_iterator it;
  //  for (it = pESRecHits->begin(); it != pESRecHits->end(); ++it) {

    edm::Handle<EcalRecHitCollection> ebRHs, eeRHs;
    e.getByToken(_ebReducedRecHitToken, ebRHs);  
    e.getByToken(_eeReducedRecHitToken, eeRHs);  
    //ebRecHits = ebRHs.product();
    //eeRecHits = eeRHs.product();


    EcalRecHitCollection::const_iterator itRHee;
    EcalRecHitCollection::const_iterator itRHeb;

    _N_RecHits_ee     = 0;
    _N_RecHits_eb     = 0;
    _N_RecHits        = 0;
    _N_RecHits_ee_p5  = 0;
    _N_RecHits_eb_p5  = 0;
    _N_RecHits_p5     = 0;
    _N_RecHits_ee_1   = 0;
    _N_RecHits_eb_1   = 0;
    _N_RecHits_1      = 0;
    _N_RecHits_ee_1p5 = 0;
    _N_RecHits_eb_1p5 = 0;
    _N_RecHits_1p5    = 0;
    _N_RecHits_ee_2   = 0;
    _N_RecHits_eb_2   = 0;
    _N_RecHits_2      = 0;
    _N_RecHits_ee_2p5 = 0;
    _N_RecHits_eb_2p5 = 0;
    _N_RecHits_2p5    = 0;
    _N_RecHits_ee_3   = 0;
    _N_RecHits_eb_3   = 0;
    _N_RecHits_3      = 0;
    _N_RecHits_ee_3p5 = 0;
    _N_RecHits_eb_3p5 = 0;
    _N_RecHits_3p5    = 0;

    for (itRHee = eeRHs->begin(); itRHee != eeRHs->end(); ++itRHee) {
     if(itRHee->energy() > 0.5 ){ _N_RecHits_ee_p5++; }
     if(itRHee->energy() > 1.0 ){ _N_RecHits_ee_1++;  _N_RecHits_ee++; }
     if(itRHee->energy() > 1.5 ){ _N_RecHits_ee_1p5++; }
     if(itRHee->energy() > 2.0 ){ _N_RecHits_ee_2++; }
     if(itRHee->energy() > 2.5 ){ _N_RecHits_ee_2p5++; }
     if(itRHee->energy() > 3.0 ){ _N_RecHits_ee_3++; }
     if(itRHee->energy() > 3.5 ){ _N_RecHits_ee_3p5++; }
    }
    for (itRHeb = ebRHs->begin(); itRHeb != ebRHs->end(); ++itRHeb) {
     if(itRHeb->energy() > 0.5 ){ _N_RecHits_eb_p5++; }
     if(itRHeb->energy() > 1.0 ){ _N_RecHits_eb_1++;  _N_RecHits_eb++; }
     if(itRHeb->energy() > 1.5 ){ _N_RecHits_eb_1p5++; }
     if(itRHeb->energy() > 2.0 ){ _N_RecHits_eb_2++; }
     if(itRHeb->energy() > 2.5 ){ _N_RecHits_eb_2p5++; }
     if(itRHeb->energy() > 3.0 ){ _N_RecHits_eb_3++; }
     if(itRHeb->energy() > 3.5 ){ _N_RecHits_eb_3p5++; }
    }
    _N_RecHits = _N_RecHits_ee + _N_RecHits_eb;
    _N_RecHits_p5 = _N_RecHits_ee_p5 + _N_RecHits_eb_p5;
    _N_RecHits_1 = _N_RecHits_ee_1 + _N_RecHits_eb_1;
    _N_RecHits_1p5 = _N_RecHits_ee_1p5 + _N_RecHits_eb_1p5;
    _N_RecHits_2 = _N_RecHits_ee_2 + _N_RecHits_eb_2;
    _N_RecHits_2p5 = _N_RecHits_ee_2p5 + _N_RecHits_eb_2p5;
    _N_RecHits_3 = _N_RecHits_ee_3 + _N_RecHits_eb_3;
    _N_RecHits_3p5 = _N_RecHits_ee_3p5 + _N_RecHits_eb_3p5;


    edm::Handle<reco::GenParticleCollection> gens;
    e.getByToken(_geninputToken,gens);
    //std::cout<<"ebsc: "<<ebSCs.isValid()<<" eesc: "<<eeSCs.isValid()<<" gen: "<<gens.isValid()<<std::endl;
//    e.getByLabel(_geninput,gens);
    if(ebSCs.isValid() && eeSCs.isValid() && gens.isValid())
    {
        matchSuperClustersToGen(*ebSCs, *eeSCs, *gens);

        std::vector< std::pair<reco::GenParticle, int> >::const_iterator it = _matchedParticles.begin();
        std::vector< std::pair<reco::GenParticle, int> >::const_iterator itE = _matchedParticles.end();
        for(;it!=itE;++it)
        {
            if(it->second!=-1)
            {
                processSuperClusterFillTree(e, es, _matchedSuperClusters[it->second], it->first);
            }
            else
            {
                processNoMatchFillTree(e, es, it->first);
            }
        }

    }
    //else 
    //{
    //    std::cout<<"hack"<<std::endl;
    //    //throw cms::Exception("PFSuperClusterRegTreeMaker")
    //    //    << "Product ID was invalid! xxx"
    //    //    << std::endl;
    //}
//    std::vector< std::pair<reco::GenParticle, int> >::const_iterator it = _matchedParticles.begin();
//    std::vector< std::pair<reco::GenParticle, int> >::const_iterator itE = _matchedParticles.end();
//    for(;it!=itE;++it)
//    {
//        if(it->second!=-1)
//        {
//            processSuperClusterFillTree(e, es, _matchedSuperClusters[it->second], it->first);
//        }
//        else
//        {
//            processNoMatchFillTree(e, es, it->first);
//        }
//    }

    delete _lazyTool;
    delete _lazyToolNoZS;
}

void PFSuperClusterRegTreeMaker::processSuperClusterFillTree(const edm::Event& e, const edm::EventSetup& es, const reco::SuperCluster& sc, const reco::GenParticle& gen) 
{
    const int N_ECAL = sc.clustersEnd() - sc.clustersBegin();
    const int N_PS   = sc.preshowerClustersEnd() - sc.preshowerClustersBegin();

    _N_ECALClusters = std::max(0,N_ECAL - 1); // minus 1 because of seed
    _N_PSClusters   = N_PS;



    // generated information
    _isMatched       = 1;
    _genEnergy       = gen.energy();
    _genPt           = gen.pt();
    _genEta          = gen.eta();
    _genPhi          = gen.phi();
    _genDEoE         = (gen.energy()-sc.energy())/gen.energy();
    _genDRToCentroid = reco::deltaR(gen,sc);
    _genDRToSeed     = reco::deltaR(gen,*(sc.seed()));

    // supercluster information
    _scIndex                += 1;
    _scRawEnergy             = sc.rawEnergy();
    _scCalibratedEnergy      = sc.energy();
    _scPreshowerEnergy       = sc.preshowerEnergy();
    _scPreshowerEnergyPlane1 = sc.preshowerEnergyPlane1();
    _scPreshowerEnergyPlane2 = sc.preshowerEnergyPlane2();
    _scIsEB                  = (sc.seed()->hitsAndFractions().at(0).first.subdetId()==EcalBarrel);
    _scEta                   = sc.position().Eta();
    _scPhi                   = sc.position().Phi();
    _scR                     = sc.position().R();
    _scPhiWidth              = sc.phiWidth();
    _scEtaWidth              = sc.etaWidth();


    // sc seed information
    edm::Ptr<reco::CaloCluster> theseed = sc.seed();
    _scSeedRawEnergy        = theseed->energy();
    _scSeedCalibratedEnergy = theseed->correctedEnergy();
    _scSeedEta              = theseed->eta();
    _scSeedPhi              = theseed->phi();
    _scSeedSize             = theseed->hitsAndFractions().size();
    // shapes with zero suppression
    float e3x3              = _lazyTool->e3x3(*(theseed));
    _scSeedR9               = e3x3/_scRawEnergy;
    _scSeedEmax             = _lazyTool->eMax(*(theseed));
    _scSeedE2nd             = _lazyTool->e2nd(*(theseed));
    _scSeedE2x5max          = _lazyTool->e2x5Max(*(theseed));
    float eLeft             = _lazyTool->eLeft(*(theseed));
    float eRight            = _lazyTool->eRight(*(theseed));
    float eTop              = _lazyTool->eTop(*(theseed));
    float eBottom           = _lazyTool->eBottom(*(theseed));
    float e2x5Left          = _lazyTool->e2x5Left(*(theseed));
    float e2x5Right         = _lazyTool->e2x5Right(*(theseed));
    float e2x5Top           = _lazyTool->e2x5Top(*(theseed));
    float e2x5Bottom        = _lazyTool->e2x5Bottom(*(theseed));
    _scSeedLeftRightAsym    = (eLeft+eRight!=0. ? (eLeft-eRight)/(eLeft+eRight) : 0.);
    _scSeedTopBottomAsym    = (eTop+eBottom!=0. ? (eTop-eBottom)/(eTop+eBottom) : 0.);
    _scSeed2x5LeftRightAsym = (e2x5Left+e2x5Right!=0. ? (e2x5Left-e2x5Right)/(e2x5Left+e2x5Right) : 0.);
    _scSeed2x5TopBottomAsym = (e2x5Top+e2x5Bottom!=0. ? (e2x5Top-e2x5Bottom)/(e2x5Top+e2x5Bottom) : 0.);
    std::vector<float> vCov = _lazyTool->localCovariances(*(theseed));
    double see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    double spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    double sep = 0.;
    if (see*spp > 0)
        sep = vCov[1] / (see * spp);
    else if (vCov[1] > 0)
        sep = 1.0;
    else
        sep = -1.0;
    _scSeedSigmaIetaIeta    = see;
    _scSeedSigmaIphiIphi    = spp;
    _scSeedSigmaIetaIphi    = sep;
    // shapes without zero suppression
    float e3x3NoZS              = _lazyToolNoZS->e3x3(*(theseed));
    _scSeedR9NoZS               = e3x3NoZS/_scRawEnergy;
    _scSeedEmaxNoZS             = _lazyToolNoZS->eMax(*(theseed));
    _scSeedE2ndNoZS             = _lazyToolNoZS->e2nd(*(theseed));
    _scSeedE2x5maxNoZS          = _lazyToolNoZS->e2x5Max(*(theseed));
    float eLeftNoZS             = _lazyToolNoZS->eLeft(*(theseed));
    float eRightNoZS            = _lazyToolNoZS->eRight(*(theseed));
    float eTopNoZS              = _lazyToolNoZS->eTop(*(theseed));
    float eBottomNoZS           = _lazyToolNoZS->eBottom(*(theseed));
    float e2x5LeftNoZS          = _lazyToolNoZS->e2x5Left(*(theseed));
    float e2x5RightNoZS         = _lazyToolNoZS->e2x5Right(*(theseed));
    float e2x5TopNoZS           = _lazyToolNoZS->e2x5Top(*(theseed));
    float e2x5BottomNoZS        = _lazyToolNoZS->e2x5Bottom(*(theseed));
    _scSeedLeftRightAsymNoZS    = (eLeftNoZS+eRightNoZS!=0. ? (eLeftNoZS-eRightNoZS)/(eLeftNoZS+eRightNoZS) : 0.);
    _scSeedTopBottomAsymNoZS    = (eTopNoZS+eBottomNoZS!=0. ? (eTopNoZS-eBottomNoZS)/(eTopNoZS+eBottomNoZS) : 0.);
    _scSeed2x5LeftRightAsymNoZS = (e2x5LeftNoZS+e2x5RightNoZS!=0. ? (e2x5LeftNoZS-e2x5RightNoZS)/(e2x5LeftNoZS+e2x5RightNoZS) : 0.);
    _scSeed2x5TopBottomAsymNoZS = (e2x5TopNoZS+e2x5BottomNoZS!=0. ? (e2x5TopNoZS-e2x5BottomNoZS)/(e2x5TopNoZS+e2x5BottomNoZS) : 0.);
    std::vector<float> vCovNoZS = _lazyToolNoZS->localCovariances(*(theseed));
    double seeNoZS = (isnan(vCovNoZS[0]) ? 0. : sqrt(vCovNoZS[0]));
    double sppNoZS = (isnan(vCovNoZS[2]) ? 0. : sqrt(vCovNoZS[2]));
    double sepNoZS = 0.;
    if (seeNoZS*sppNoZS > 0)
        sepNoZS = vCovNoZS[1] / (seeNoZS * sppNoZS);
    else if (vCovNoZS[1] > 0)
        sepNoZS = 1.0;
    else
        sepNoZS = -1.0;
    _scSeedSigmaIetaIetaNoZS    = seeNoZS;
    _scSeedSigmaIphiIphiNoZS    = sppNoZS;
    _scSeedSigmaIetaIphiNoZS    = sepNoZS;

    EcalClusterLocal ecalLocal;
    DetId seedid = theseed->seed();
    if(theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel)
    {
        float cryPhi, cryEta, thetatilt, phitilt;
        int ieta, iphi;
        ecalLocal.localCoordsEB(*(theseed), es, cryEta, cryPhi, ieta, iphi, thetatilt, phitilt);
        _scSeedCryEta  = cryEta;
        _scSeedCryPhi  = cryPhi;
        _scSeedCryIeta = ieta;
        _scSeedCryIphi = iphi;
        _scSeedCryX  = 0;
        _scSeedCryY  = 0;
        _scSeedCryIx = 0;
        _scSeedCryIy = 0;
        //
        EBDetId ebseedid(seedid);
        int ietaV2 = ebseedid.ieta();
        int iphiV2 = ebseedid.iphi();
        _scSeedCryIetaV2 = ietaV2;
        _scSeedCryIphiV2 = iphiV2;
        _scSeedCryIxV2 = 0;
        _scSeedCryIyV2 = 0;
    }
    else
    {
        float cryX, cryY, thetatilt, phitilt;
        int ix, iy;
        ecalLocal.localCoordsEE(*(theseed), es, cryX, cryY, ix, iy, thetatilt, phitilt);
        _scSeedCryX  = cryX;
        _scSeedCryY  = cryY;
        _scSeedCryIx = ix;
        _scSeedCryIy = iy;
        _scSeedCryEta  = 0;
        _scSeedCryPhi  = 0;
        _scSeedCryIeta = 0;
        _scSeedCryIphi = 0;
        //
        EEDetId eeseedid(seedid);
        int ixV2 = eeseedid.ix();
        int iyV2 = eeseedid.iy();
        _scSeedCryIetaV2 = 0;
        _scSeedCryIphiV2 = 0;
        _scSeedCryIxV2 = ixV2;
        _scSeedCryIyV2 = iyV2;
    }


    // loop over all clusters that aren't the seed
    setTreeArraysForSize(_N_ECALClusters, _N_PSClusters);
    _clusterMaxDR     = 999.;
    _clusterMaxDRDPhi = 999.;
    _clusterMaxDRDEta = 999.;
    _clustersMeanDRToSeed = 999.;
    _clustersMeanDEtaToSeed = 999.;
    _clustersMeanDPhiToSeed = 999.;
    _clusterMaxDRRawEnergy      = 0.;
    _clustersMeanRawEnergy      = 0.;
    _clustersRMSRawEnergy       = 0.;
    float clustersMeanSquareRawEnergy = 0.;
    float subclustersRawEnergy = 0.;
    auto clusend = sc.clustersEnd();
    size_t iclus = 0;
    float maxDR = 0;
    edm::Ptr<reco::CaloCluster> pclus;
    for( auto clus = sc.clustersBegin(); clus != clusend; ++clus ) 
    {
        pclus = *clus;

        _clustersMeanRawEnergy += pclus->energy();
        clustersMeanSquareRawEnergy += (pclus->energy()*pclus->energy());

        if( theseed == pclus ) continue;
        _clusterRawEnergy.get()[iclus]      = pclus->energy();

        if(iclus==0){
         _clusterRawEnergy_0 = pclus->energy();
         _clusterDPhiToSeed_0 = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi());
         _clusterDEtaToSeed_0 = pclus->eta() - theseed->eta();
        }
        if(iclus==1){
         _clusterRawEnergy_1 = pclus->energy();
         _clusterDPhiToSeed_1 = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi());
         _clusterDEtaToSeed_1 = pclus->eta() - theseed->eta();
        }
        if(iclus==2){
         _clusterRawEnergy_2 = pclus->energy();
         _clusterDPhiToSeed_2 = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi());
         _clusterDEtaToSeed_2 = pclus->eta() - theseed->eta();
        }
        _clusterCalibEnergy.get()[iclus]    = pclus->correctedEnergy();
        _clusterEta.get()[iclus]            = pclus->eta();
        _clusterPhi.get()[iclus]            = pclus->phi();
        _clusterDPhiToSeed.get()[iclus]     = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi());
        _clusterDEtaToSeed.get()[iclus]     = pclus->eta() - theseed->eta();
        _clusterDPhiToCentroid.get()[iclus] = TVector2::Phi_mpi_pi(pclus->phi() - sc.phi());
        _clusterDEtaToCentroid.get()[iclus] = pclus->eta() - sc.eta();
        // find cluster with max dR
        if(reco::deltaR(*pclus, *theseed) > maxDR)
        {
            maxDR = reco::deltaR(*pclus, *theseed);
            _clusterMaxDR = maxDR;
            _clusterMaxDRDPhi = _clusterDPhiToSeed.get()[iclus];
            _clusterMaxDRDEta = _clusterDEtaToSeed.get()[iclus];
            _clusterMaxDRRawEnergy = _clusterRawEnergy.get()[iclus];
        }

        subclustersRawEnergy += pclus->energy();
        _clustersMeanDRToSeed   = reco::deltaR(*pclus, *theseed)*pclus->energy();
        _clustersMeanDEtaToSeed = (pclus->eta() - theseed->eta())*pclus->energy();
        _clustersMeanDPhiToSeed = TVector2::Phi_mpi_pi(pclus->phi() - theseed->phi())*pclus->energy();

        _clusterDPhiToGen.get()[iclus] = TVector2::Phi_mpi_pi(pclus->phi() - gen.phi());
        _clusterDEtaToGen.get()[iclus] = pclus->eta() - gen.eta();
        _clusterInMustache.get()[iclus] = (int) MK::inMustache(theseed->eta(),
                theseed->phi(),
                pclus->energy(),
                pclus->eta(),
                pclus->phi());
        _clusterInDynDPhi.get()[iclus] = (int) MK::inDynamicDPhiWindow(
                    theseed->hitsAndFractions().at(0).first.subdetId()==EcalBarrel,
                    theseed->phi(),
                    pclus->energy(),
                    pclus->eta(),
                    pclus->phi());
        //if(IsLinkedByRecHit(theseed, pclus,50., 0.5))
        //    _clusterLeakage.get()[iclus] = 1;
        //else 
        //    _clusterLeakage.get()[iclus] = 0;

        ++iclus;
    }
    // compute mean and rms values
    _clustersMeanRawEnergy /= (double)(_N_ECALClusters+1);
    clustersMeanSquareRawEnergy /= (double)(_N_ECALClusters+1);
    _clustersRMSRawEnergy = sqrt(clustersMeanSquareRawEnergy - _clustersMeanRawEnergy*_clustersMeanRawEnergy);

    if(subclustersRawEnergy>0.)
    {
        _clustersMeanDRToSeed /= subclustersRawEnergy;
        _clustersMeanDEtaToSeed /= subclustersRawEnergy;
        _clustersMeanDPhiToSeed /= subclustersRawEnergy;
    }

    // loop over all preshower clusters 
    auto psclusend = sc.preshowerClustersEnd();
    size_t ipsclus = 0;
    edm::Ptr<reco::CaloCluster> ppsclus;
    for( auto psclus = sc.preshowerClustersBegin(); psclus != psclusend; ++psclus ) 
    {
        ppsclus = *psclus;
        _psClusterRawEnergy.get()[ipsclus] = ppsclus->energy();    
        _psClusterEta.get()[ipsclus]       = ppsclus->eta();    
        _psClusterPhi.get()[ipsclus]       = ppsclus->phi();
        ++ipsclus;
    }
    _tree->Fill();
}

void PFSuperClusterRegTreeMaker::processNoMatchFillTree(const edm::Event& e, const edm::EventSetup& es, const reco::GenParticle& gen) 
{

    _N_ECALClusters = 0; 
    _N_PSClusters   = 0;

    // generated information
    _isMatched       = 0;
    _genEnergy       = gen.energy();
    _genPt           = gen.pt();
    _genEta          = gen.eta();
    _genPhi          = gen.phi();
    _genDEoE         = 999.;
    _genDRToCentroid = 999.;
    _genDRToSeed     = 999.;

    // supercluster information
    _scIndex                += 1;
    _scRawEnergy             = 0;
    _scCalibratedEnergy      = 0;
    _scPreshowerEnergy       = 0;
    _scPreshowerEnergyPlane1 = 0;
    _scPreshowerEnergyPlane2 = 0;
    _scIsEB                  = 0;
    _scEta                   = 0;
    _scPhi                   = 0;
    _scR                     = 0;
    _scPhiWidth              = 0;
    _scEtaWidth              = 0;


    // sc seed information
    _scSeedRawEnergy        = 0;
    _scSeedCalibratedEnergy = 0;
    _scSeedEta              = 0;
    _scSeedPhi              = 0;
    _scSeedSize             = 0;
    _scSeedR9               = 0;
    _scSeedEmax             = 0;
    _scSeedE2nd             = 0;
    _scSeedE2x5max          = 0;
    _scSeedLeftRightAsym    = 0;
    _scSeedTopBottomAsym    = 0;
    _scSeed2x5LeftRightAsym = 0;
    _scSeed2x5TopBottomAsym = 0;
    _scSeedSigmaIetaIeta    = 0;
    _scSeedSigmaIphiIphi    = 0;
    _scSeedSigmaIetaIphi    = 0;
    _scSeedR9NoZS               = 0;
    _scSeedEmaxNoZS             = 0;
    _scSeedE2ndNoZS             = 0;
    _scSeedE2x5maxNoZS          = 0;
    _scSeedLeftRightAsymNoZS    = 0;
    _scSeedTopBottomAsymNoZS    = 0;
    _scSeed2x5LeftRightAsymNoZS = 0;
    _scSeed2x5TopBottomAsymNoZS = 0;
    _scSeedSigmaIetaIetaNoZS    = 0;
    _scSeedSigmaIphiIphiNoZS    = 0;
    _scSeedSigmaIetaIphiNoZS    = 0;

    _scSeedCryEta    = 0;
    _scSeedCryPhi    = 0;
    _scSeedCryIeta   = 0;
    _scSeedCryIphi   = 0;
    _scSeedCryIetaV2 = 0;
    _scSeedCryIphiV2 = 0;
    _scSeedCryX      = 0;
    _scSeedCryY      = 0;
    _scSeedCryIx     = 0;
    _scSeedCryIy     = 0;
    _scSeedCryIxV2   = 0;
    _scSeedCryIyV2   = 0;



    // loop over all clusters that aren't the seed
    setTreeArraysForSize(0, 0);
    _clusterMaxDR     = 999.;
    _clusterMaxDRDPhi = 999.;
    _clusterMaxDRDEta = 999.;
    _clustersMeanDRToSeed = 999.;
    _clustersMeanDEtaToSeed = 999.;
    _clustersMeanDPhiToSeed = 999.;
    _clusterMaxDRRawEnergy      = 0.;
    _clustersMeanRawEnergy      = 0.;
    _clustersRMSRawEnergy       = 0.;
    

    _tree->Fill();
}




void PFSuperClusterRegTreeMaker::setTreeArraysForSize(const size_t N_ECAL, const size_t N_PS) 
{
    float* cRE_new        = new float[N_ECAL];
    float* cCE_new        = new float[N_ECAL];
    float* cEta_new       = new float[N_ECAL];
    float* cPhi_new       = new float[N_ECAL];
    float* cDPhiSeed_new  = new float[N_ECAL];
    float* cDEtaSeed_new  = new float[N_ECAL];
    float* cDPhiCntr_new  = new float[N_ECAL];
    float* cDEtaCntr_new  = new float[N_ECAL];
    int*   cInMust_new    = new int[N_ECAL];
    int*   cInDynDPhi_new = new int[N_ECAL];
    float* psRE_new       = new float[N_PS];
    float* psEta_new      = new float[N_PS];
    float* psPhi_new      = new float[N_PS];

    _clusterRawEnergy     .reset(cRE_new       , array_deleter<float>());
    _clusterCalibEnergy   .reset(cCE_new       , array_deleter<float>());
    _clusterEta           .reset(cEta_new      , array_deleter<float>());
    _clusterPhi           .reset(cPhi_new      , array_deleter<float>());
    _clusterDPhiToSeed    .reset(cDPhiSeed_new , array_deleter<float>());
    _clusterDEtaToSeed    .reset(cDEtaSeed_new , array_deleter<float>());  
    _clusterDPhiToCentroid.reset(cDPhiCntr_new , array_deleter<float>());
    _clusterDEtaToCentroid.reset(cDEtaCntr_new , array_deleter<float>());
    _clusterInMustache    .reset(cInMust_new   , array_deleter<int>());
    _clusterInDynDPhi     .reset(cInDynDPhi_new, array_deleter<int>());
    _psClusterRawEnergy   .reset(psRE_new      , array_deleter<float>());
    _psClusterEta         .reset(psEta_new     , array_deleter<float>());
    _psClusterPhi         .reset(psPhi_new     , array_deleter<float>());

    _tree->GetBranch("clusterRawEnergy")     ->SetAddress(_clusterRawEnergy.get());
    _tree->GetBranch("clusterCalibEnergy")   ->SetAddress(_clusterCalibEnergy.get());
    _tree->GetBranch("clusterEta")           ->SetAddress(_clusterEta.get());
    _tree->GetBranch("clusterPhi")           ->SetAddress(_clusterPhi.get());
    _tree->GetBranch("clusterDPhiToSeed")    ->SetAddress(_clusterDPhiToSeed.get());
    _tree->GetBranch("clusterDEtaToSeed")    ->SetAddress(_clusterDEtaToSeed.get());
    _tree->GetBranch("clusterDPhiToCentroid")->SetAddress(_clusterDPhiToCentroid.get());
    _tree->GetBranch("clusterDEtaToCentroid")->SetAddress(_clusterDEtaToCentroid.get());
    _tree->GetBranch("clusterInMustache")    ->SetAddress(_clusterInMustache.get());
    _tree->GetBranch("clusterInDynDPhi")     ->SetAddress(_clusterInDynDPhi.get());
    _tree->GetBranch("psClusterRawEnergy")   ->SetAddress(_psClusterRawEnergy.get());
    _tree->GetBranch("psClusterEta")         ->SetAddress(_psClusterEta.get());
    _tree->GetBranch("psClusterPhi")         ->SetAddress(_psClusterPhi.get());

    float* cDPhiGen_new = new float[N_ECAL];
    float* cDEtaGen_new = new float[N_ECAL];
    _clusterDPhiToGen.reset(cDPhiGen_new, array_deleter<float>());
    _clusterDEtaToGen.reset(cDEtaGen_new, array_deleter<float>());
    _tree->GetBranch("clusterDPhiToGen") ->SetAddress(_clusterDPhiToGen.get());
    _tree->GetBranch("clusterDEtaToGen") ->SetAddress(_clusterDEtaToGen.get());
}


bool PFSuperClusterRegTreeMaker::IsLinkedByRecHit(const edm::Ptr<reco::PFCluster>& the_seed, const edm::Ptr<reco::PFCluster>& x, const double threshold, const double majority)
{

    if( the_seed->energy() < threshold ) 
        return false; 
    const auto& seedHitsAndFractions = the_seed->hitsAndFractions();
    const auto& xHitsAndFractions    = x->hitsAndFractions();      
    float x_rechits_tot              = xHitsAndFractions.size();
    float x_rechits_match            = 0.0;      
    for( const std::pair<DetId, float>& seedHit : seedHitsAndFractions ) 
    {
        for( const std::pair<DetId, float>& xHit : xHitsAndFractions ) 
        {
            if( seedHit.first == xHit.first ) 
            {	    
                x_rechits_match += 1.0;
            }
        }	
    }      
    return x_rechits_match/x_rechits_tot > majority;

}

void PFSuperClusterRegTreeMaker::matchSuperClustersToGen(const reco::SuperClusterCollection& ebSCs, const reco::SuperClusterCollection& eeSCs, const reco::GenParticleCollection& gens)
{
    for( const auto& gen : gens ) 
    {
        if( abs(gen.pdgId()) == _genParticleID &&  gen.status()==1) // looking at final state electrons or photons
        {
            double minDr = 1e6;
            double minDe = 1e6;
            double minDeDr = 1e6;
            double this_dr;
            double this_de;
            double this_dedr;
            const reco::SuperCluster* bestmatch = NULL;
            // loop on EB SCs
            for( const auto& sc : ebSCs ) 
            {
                math::XYZVector direction = sc.position() - gen.vertex();
                double energy = sc.energy();
                math::XYZVector momentum = direction.unit() * energy;
                math::XYZTLorentzVector physicalSCP4(momentum.x(), momentum.y(), momentum.z(), energy );
                this_dr = reco::deltaR(gen,physicalSCP4);
                this_de = fabs(gen.energy()-sc.energy())/gen.energy();
                this_dedr = sqrt(this_dr*this_dr + this_de*this_de);
                //if(this_dr < _matchMaxDR && sc.energy()/gen.energy()>_matchMinERatio && (this_dr < minDr || this_de<minDe))
                if(this_dr < _matchMaxDR && sc.energy()/gen.energy()>_matchMinERatio)
                {
                    if( (_matchingType==0 && this_dr<minDr) || 
                            (_matchingType==1 && this_de<minDe) ||
                            (_matchingType==2 && this_dedr<minDeDr))
                    {
                        minDr = this_dr;
                        minDe = this_de;
                        minDeDr = this_dedr;
                        bestmatch = &sc;
                    }
                }
            }
            // loop on EE SCs
            for( const auto& sc : eeSCs ) 
            {
                math::XYZVector direction = sc.position() - gen.vertex();
                double energy = sc.energy();
                math::XYZVector momentum = direction.unit() * energy;
                math::XYZTLorentzVector physicalSCP4(momentum.x(), momentum.y(), momentum.z(), energy );
                this_dr = reco::deltaR(gen,physicalSCP4);
                this_de = fabs(gen.energy()-sc.energy())/gen.energy();
                this_dedr = sqrt(this_dr*this_dr + this_de*this_de);
                //if(this_dr < _matchMaxDR && sc.energy()/gen.energy()>_matchMinERatio && (this_dr < minDr || this_de<minDe))
                if(this_dr < _matchMaxDR && sc.energy()/gen.energy()>_matchMinERatio)
                {
                    if( (_matchingType==0 && this_dr<minDr) ||
                            (_matchingType==1 && this_de<minDe) ||
                            (_matchingType==2 && this_dedr<minDeDr))
                    {
                        minDr = this_dr;
                        minDe = this_de;
                        minDeDr = this_dedr;
                        bestmatch = &sc;
                    }
                }
            }
            if(bestmatch)
            {
                _matchedSuperClusters.push_back(*bestmatch);
                _matchedParticles.push_back( std::make_pair(gen, (int)_matchedSuperClusters.size()-1) );
            }
            else
            {
                _matchedParticles.push_back( std::make_pair(gen, -1) );
            }
        }
    }
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFSuperClusterRegTreeMaker);
