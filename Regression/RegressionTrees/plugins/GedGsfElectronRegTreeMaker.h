

#ifndef GEDGSFELECTRONREGTREEMAKER_H
#define GEDGSFELECTRONREGTREEMAKER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include <memory>
#include <algorithm>

typedef edm::ParameterSet PSet;

namespace 
{  
  template<typename T>
  struct array_deleter
  {
    void operator () (T* arr) { delete [] arr; }
  };
}

class TTree;

class GedGsfElectronRegTreeMaker : public edm::EDAnalyzer 
{
    public:
        GedGsfElectronRegTreeMaker(const PSet&);
        ~GedGsfElectronRegTreeMaker() {}

        void analyze(const edm::Event&, const edm::EventSetup&);

    private:    
        edm::Service<TFileService> _fs;
        bool _dogen;
        edm::InputTag _geninput;
        double _matchMinERatio;
        double _matchMaxDR;
        int _matchingType;
        edm::InputTag _eleInput;
        edm::EDGetTokenT<EcalRecHitCollection> _ebReducedRecHitToken;
        edm::EDGetTokenT<EcalRecHitCollection> _eeReducedRecHitToken;
        //edm::InputTag _ebReducedRecHitCollection;
        //edm::InputTag _eeReducedRecHitCollection;
        edm::InputTag _PrimaryVertex;
        std::shared_ptr<PFEnergyCalibration> _calib;
        void processSuperClusterFillTree(const edm::Event&, const edm::EventSetup&, const reco::GsfElectron&, const reco::SuperCluster&, const reco::GenParticle&);
        void processNoMatchFillTree(const edm::Event&, const edm::EventSetup&, const reco::GenParticle&);
        void matchElectronToGen(const reco::GsfElectronCollection& eles, const reco::GenParticleCollection& gens);

        EcalClusterLazyTools *_lazyTool;
        noZS::EcalClusterLazyTools *_lazyToolNoZS;

        std::vector< reco::GsfElectron > _matchedElectrons;
        std::vector< std::pair<reco::GenParticle, int> > _matchedParticles;

        // the tree  
        void setTreeArraysForSize(const size_t N_ECAL,const size_t N_PS);
        bool IsLinkedByRecHit(const edm::Ptr<reco::PFCluster>& the_seed, const edm::Ptr<reco::PFCluster>& x,const double threshold, const double majority);
        TTree* _tree;

        int _eventNumber;
        int _nVtx;
        int _nElectrons;
        float _rho;

        // supercluster variables
        int   _scIndex;
        float _scRawEnergy;
        float _scCalibratedEnergy;
        float _scPreshowerEnergy;
        float _scPreshowerEnergyPlane1;
        float _scPreshowerEnergyPlane2;
        int   _scIsEB;
        float _scEta;
        float _scPhi;
        float _scR;
        float _scPhiWidth;
        float _scEtaWidth;
        float _scSeedRawEnergy; 
        float _scSeedCalibratedEnergy;
        float _scSeedEta;
        float _scSeedPhi;
        float _scSeedSize;
        float _scSeedR9;
        float _scSeedEmax;
        float _scSeedE2nd;
        float _scSeedLeftRightAsym;
        float _scSeedTopBottomAsym;
        float _scSeedE2x5max;
        float _scSeed2x5LeftRightAsym;
        float _scSeed2x5TopBottomAsym;
        float _scSeedSigmaIetaIeta;
        float _scSeedSigmaIetaIphi;
        float _scSeedSigmaIphiIphi;
        float _scSeedR9NoZS;
        float _scSeedEmaxNoZS;
        float _scSeedE2ndNoZS;
        float _scSeedLeftRightAsymNoZS;
        float _scSeedTopBottomAsymNoZS;
        float _scSeedE2x5maxNoZS;
        float _scSeed2x5LeftRightAsymNoZS;
        float _scSeed2x5TopBottomAsymNoZS;
        float _scSeedSigmaIetaIetaNoZS;
        float _scSeedSigmaIetaIphiNoZS;
        float _scSeedSigmaIphiIphiNoZS;
        float _scSeedCryEta;
        float _scSeedCryPhi;
        int   _scSeedCryIeta;
        int   _scSeedCryIphi;
        int   _scSeedCryIetaV2;
        int   _scSeedCryIphiV2;
        float _scSeedCryX;
        float _scSeedCryY;
        int   _scSeedCryIx;
        int   _scSeedCryIy;
        int   _scSeedCryIxV2;
        int   _scSeedCryIyV2;

        // generator variables
        int   _isMatched;
        float _genEnergy; 
        float _genPt;
        float _genEta;
        float _genPhi;
        float _genDEoE;
        float _genDRToCentroid;
        float _genDRToSeed;

        // clusters variables
        int _N_ECALClusters;
        float _clusterMaxDR;
        float _clusterMaxDRDPhi;
        float _clusterMaxDRDEta;
        float _clusterMaxDRRawEnergy;
        float _clustersMeanRawEnergy;
        float _clustersRMSRawEnergy;
        float _clustersMeanDRToSeed;
        float _clustersMeanDEtaToSeed;
        float _clustersMeanDPhiToSeed;
        std::shared_ptr<float> _clusterRawEnergy;
        std::shared_ptr<float> _clusterCalibEnergy;
        std::shared_ptr<float> _clusterEta;
        std::shared_ptr<float> _clusterPhi;
        std::shared_ptr<float> _clusterDPhiToSeed;
        std::shared_ptr<float> _clusterDEtaToSeed;
        std::shared_ptr<float> _clusterDPhiToCentroid;
        std::shared_ptr<float> _clusterDEtaToCentroid; 
        std::shared_ptr<float> _clusterDPhiToGen;
        std::shared_ptr<float> _clusterDEtaToGen;
        std::shared_ptr<int>   _clusterInMustache;
        std::shared_ptr<int>   _clusterInDynDPhi;
        std::shared_ptr<int>   _clusterLeakage;
        int _N_PSClusters;
        std::shared_ptr<float> _psClusterRawEnergy;
        std::shared_ptr<float> _psClusterEta;
        std::shared_ptr<float> _psClusterPhi;

        // track-combination variables
        float _trkMomentum;
        float _trkMomentumError;
        float _trkMomentumRelError;
        int   _eleClass;
        int   _eleEcalDriven;
        int   _eleTrackerDriven;

};


#endif
