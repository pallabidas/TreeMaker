#ifndef __ELECTRON_INFO_H_
#define __ELECTRON_INFO_H_

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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

class electronInfo : public baseTree{

 public:
  electronInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~electronInfo();
  void Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx, float& rhoLepton);
  void SetBranches();
  void Clear();
  
 private:
  electronInfo(){};
  edm::InputTag electronLabel_;

  int Electron_n;
  TClonesArray *Electron_4Momentum;
  std::vector<float> Electron_idPreselection;
  double pfIsoCharged, pfIsoNeutral, EffArea, correction, pfIsoPUSubtracted, relIso, dxy, dz;
  float scEta, dEtaInSeed, numMissingHits, idPreselection;
};

#endif

