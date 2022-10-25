#ifndef __JET_INFO_H_
#define __JET_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//#include "FWCore/Framework/interface/ConsumesCollector.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Framework/interface/ConsumesCollector.h"
//#include "FWCore/Framework/interface/EventSetup.h"

using namespace btagbtvdeep;

class jetInfo : public baseTree{

 public:
  jetInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~jetInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, TrackInfoBuilder trackinfo);
  void SetBranches();
  void Clear();

 private:
  jetInfo(){};
  edm::InputTag jetLabel_;
  edm::InputTag genpartLabel_;
  edm::InputTag genjetLabel_;
  edm::InputTag vtxLabel_;
  edm::InputTag svLabel_;

  //edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

  float jet_pt;
  float jet_eta;
  float jet_phi;
  float jet_mass;
  //float jet_mass_raw;
  //std::vector<float> jet_pnet_probHbb;
  //std::vector<float> jet_pnet_probHcc;
  //std::vector<float> jet_pnet_probHqq;
  //std::vector<float> jet_pnet_probQCDbb;
  //std::vector<float> jet_pnet_probQCDb;
  //std::vector<float> jet_pnet_probQCDcc;
  //std::vector<float> jet_pnet_probQCDc;
  //std::vector<float> jet_pnet_probQCDothers;
  //std::vector<float> jet_pnet_regmass;
  float jet_ncand;
  float jet_nel;
  float jet_nmu;
  float jet_nbhad;
  float jet_nchad;
  float jet_hflav;
  float jet_pflav;
  //float jet_muflav;
  //float jet_elflav;
  float jet_lepflav;
  //float jet_tauflav;
  //float jet_taudecaymode;
  //float jet_flav_2prong_partonjet_match;
  //float jet_flav_2prong_parton_match;
  float jet_deepcsv_probb;
  float jet_deepcsv_probc;
  float jet_deepcsv_probudsg;
  float jet_deepjet_probb;
  float jet_deepjet_probc;
  float jet_deepjet_probuds;
  float jet_deepjet_probg;
  float jet_deepjet_probbb;
  float jet_deepjet_problepb;
  float jet_genmatch_pt;
  float jet_genmatch_eta;
  float jet_genmatch_phi;
  float jet_genmatch_mass;

  float nsv;
  std::vector<float> jet_sv_pt;
  std::vector<float> jet_sv_pt_log;
  std::vector<float> jet_sv_ptrel;
  std::vector<float> jet_sv_ptrel_log;
  std::vector<float> jet_sv_eta;
  std::vector<float> jet_sv_phi;
  std::vector<float> jet_sv_mass;
  std::vector<float> jet_sv_energy;
  std::vector<float> jet_sv_energy_log;
  std::vector<float> jet_sv_erel;
  std::vector<float> jet_sv_erel_log;
  std::vector<float> jet_sv_deta;
  std::vector<float> jet_sv_dphi;
  std::vector<float> jet_sv_chi2;
  std::vector<float> jet_sv_dxy;
  std::vector<float> jet_sv_dxysig;
  std::vector<float> jet_sv_d3d;
  std::vector<float> jet_sv_d3dsig;
  std::vector<float> jet_sv_ntrack;

  std::vector<float> jet_pfcand_px;
  std::vector<float> jet_pfcand_py;
  std::vector<float> jet_pfcand_pz;
  std::vector<float> jet_pfcand_energy;
  std::vector<float> jet_pfcand_pt;
  std::vector<float> jet_pfcand_pt_log;
  std::vector<float> jet_pfcand_e_log;
  std::vector<float> jet_pfcand_phirel;
  std::vector<float> jet_pfcand_etarel;
  std::vector<float> jet_pfcand_abseta;
  std::vector<float> jet_pfcand_puppiw;
  std::vector<float> jet_pfcand_charge;
  std::vector<float> jet_pfcand_isEl;
  std::vector<float> jet_pfcand_isMu;
  std::vector<float> jet_pfcand_isChargedHad;
  std::vector<float> jet_pfcand_isGamma;
  std::vector<float> jet_pfcand_isNeutralHad;
  std::vector<float> jet_pfcand_hcalFrac;
  std::vector<float> jet_pfcand_hcalFracCalib;
  std::vector<float> jet_pfcand_VTX_ass;
  std::vector<float> jet_pfcand_fromPV;
  std::vector<float> jet_pfcand_lostInnerHits;
  std::vector<float> jet_pfcand_trackHighPurity;
  std::vector<float> jet_pfcand_dz;
  std::vector<float> jet_pfcand_dzsig;
  std::vector<float> jet_pfcand_dxy;
  std::vector<float> jet_pfcand_dxysig;
  std::vector<float> jet_pfcand_normchi2;
  std::vector<float> jet_pfcand_quality;
  std::vector<float> jet_pfcand_btagEtaRel;
  std::vector<float> jet_pfcand_btagPtRatio;
  std::vector<float> jet_pfcand_btagPParRatio;
  std::vector<float> jet_pfcand_btagSip2dVal;
  std::vector<float> jet_pfcand_btagSip2dSig;
  std::vector<float> jet_pfcand_btagSip3dVal;
  std::vector<float> jet_pfcand_btagSip3dSig;
  std::vector<float> jet_pfcand_btagJetDistVal;

};

#endif

