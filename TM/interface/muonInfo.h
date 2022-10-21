#ifndef __MUON_INFO_H_
#define __MUON_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

class muonInfo : public baseTree{

 public:
  muonInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~muonInfo();
  void Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx, float& rhoLepton);
  void SetBranches();
  void Clear();

 private:
  muonInfo(){};
  edm::InputTag muonLabel_;

  int Muon_n;
  TClonesArray *Muon_4Momentum;
  std::vector<bool> Muon_isLooseMuon;
  std::vector<float> Muon_idPreselection;
  double pfIsoCharged, pfIsoNeutral, EffArea, correction, pfIsoPUSubtracted, relIso, dxy, dz;
  float idPreselection; 
};

#endif

