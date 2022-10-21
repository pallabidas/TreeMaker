#ifndef __PHOTON_INFO_H_
#define __PHOTON_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "TreeMaker/TM/interface/electronInfo.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h" 
//#include "Geometry/CaloTopology/interface/CaloTopology.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h" 
using namespace std;

class photonInfo : public baseTree{

 public:
  photonInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~photonInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();

  
 private:
  photonInfo(){};
  edm::InputTag photonLabel_;

  //variables which would become branches
  int Photon_n;
  TClonesArray *Photon_4Momentum, *Photon_4Momentum_corr, *Photon_Vposition; 
 
  std::vector<float> Photon_HoverE;
  std::vector<float> Photon_sigmaIetaIeta;
  std::vector<float> Photon_chargedHadronIso;
  std::vector<float> Photon_neutralHadronIso;
  std::vector<float> Photon_photonIso;
  std::vector<float> Photon_r9;
  std::vector<bool> Photon_hasPixelSeed;
  std::vector<bool> Photon_passElectronVeto;
  std::vector<float> Photon_hOVERe;
  std::vector<float> Photon_full5x5_r9;
  std::vector<float> Photon_full5x5_sigmaIetaIeta;
  std::vector<float> Photon_full5x5_e5x5;
  std::vector<float> Photon_scEta;
  std::vector<float> Photon_scEnergy;
  
};

#endif
