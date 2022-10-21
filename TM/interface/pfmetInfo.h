#ifndef __PFMET_INFO_H_
#define __PFMET_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"


class pfmetInfo : public baseTree{

 public:
  pfmetInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~pfmetInfo();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  float correct_phi(float phi);
  bool isMC_;

 private:
  pfmetInfo(){};
  edm::InputTag pfmetLabel_, puppimetLabel_, ecalBadCalibLabel_, metfilterspatLabel_, metfiltersrecoLabel_;

  //variables which would become branches
  std::vector<float> PFMetPt, PFMetPx, PFMetPy, PFMetPhi;
  std::vector<float> PFMetPt_JetResUp, PFMetPt_JetResDown, PFMetPt_JetEnUp, PFMetPt_JetEnDown, PFMetPt_UnclusteredEnUp, PFMetPt_UnclusteredEnDown;
  std::vector<float> PFMetPx_JetResUp, PFMetPx_JetResDown, PFMetPx_JetEnUp, PFMetPx_JetEnDown, PFMetPx_UnclusteredEnUp, PFMetPx_UnclusteredEnDown;
  std::vector<float> PFMetPy_JetResUp, PFMetPy_JetResDown, PFMetPy_JetEnUp, PFMetPy_JetEnDown, PFMetPy_UnclusteredEnUp, PFMetPy_UnclusteredEnDown;
  std::vector<float> PFMetPhi_JetResUp, PFMetPhi_JetResDown, PFMetPhi_JetEnUp, PFMetPhi_JetEnDown, PFMetPhi_UnclusteredEnUp, PFMetPhi_UnclusteredEnDown;
  std::vector<float> RawMetPt, RawMetPx, RawMetPy, RawMetPhi, RawChsMetPt, RawChsMetPx, RawChsMetPy, RawChsMetPhi, RawTrkMetPt, RawTrkMetPx, RawTrkMetPy, RawTrkMetPhi;
  std::vector<float> PuppiMetPt, PuppiMetPx, PuppiMetPy, PuppiMetPhi, PuppiMetcorPt, PuppiMetcorPx, PuppiMetcorPy, PuppiMetcorPhi;
  std::vector<float> PuppiMetPt_JetResUp, PuppiMetPt_JetResDown, PuppiMetPt_JetEnUp, PuppiMetPt_JetEnDown, PuppiMetPt_UnclusteredEnUp, PuppiMetPt_UnclusteredEnDown;
  std::vector<float> PuppiMetPx_JetResUp, PuppiMetPx_JetResDown, PuppiMetPx_JetEnUp, PuppiMetPx_JetEnDown, PuppiMetPx_UnclusteredEnUp, PuppiMetPx_UnclusteredEnDown;
  std::vector<float> PuppiMetPy_JetResUp, PuppiMetPy_JetResDown, PuppiMetPy_JetEnUp, PuppiMetPy_JetEnDown, PuppiMetPy_UnclusteredEnUp, PuppiMetPy_UnclusteredEnDown;
  std::vector<float> PuppiMetPhi_JetResUp, PuppiMetPhi_JetResDown, PuppiMetPhi_JetEnUp, PuppiMetPhi_JetEnDown, PuppiMetPhi_UnclusteredEnUp, PuppiMetPhi_UnclusteredEnDown;
  std::vector<bool> filter_goodVertices, filter_globaltighthalo2016, filter_globalsupertighthalo2016, filter_hbher2t, filter_hbheiso, filter_ecaltp, filter_ecalsc, filter_badPFMuon, filter_ecalBadCalib;
};

#endif

