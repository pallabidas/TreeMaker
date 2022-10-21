#ifndef __TAU_INFO_H_
#define __TAU_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"

class tauInfo : public baseTree{

 public:
  tauInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~tauInfo();
  void Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx);
  void SetBranches();
  void Clear();

 private:
  tauInfo(){};
  edm::InputTag tauLabel_;

  //variables which would become branches
  int Tau_n;
  TClonesArray *Tau_4Momentum;
  std::vector<int>   Tau_charge;
  std::vector<bool>  Tau_decayModeNewDM;
  std::vector<float>  Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw, Tau_byIsolationMVArun2v1DBnewDMwLTraw,\
	Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, \
	Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, Tau_againstMuonLoose3, Tau_againstMuonTight3, Tau_againstElectronLooseMVA6, \
	Tau_againstElectronMediumMVA6, Tau_againstElectronTightMVA6;
  std::vector<float> Tau_dxy, Tau_dz, Tau_dxy_approx, Tau_dz_approx, Tau_dxy_new, Tau_dz_new, Tau_pvDiff, Tau_z_impact;
  std::vector<float> Tau_chargedIsoPtSum, Tau_chargedIsoPtSumdR03, Tau_footprintCorrection, Tau_footprintCorrectiondR03, Tau_neutralIsoPtSum, Tau_neutralIsoPtSumWeight, Tau_neutralIsoPtSumWeightdR03, Tau_neutralIsoPtSumdR03, Tau_photonPtSumOutsideSignalCone, Tau_photonPtSumOutsideSignalConedR03, Tau_puCorrPtSum; 
};

#endif

