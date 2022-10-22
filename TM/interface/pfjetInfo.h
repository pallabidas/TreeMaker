#ifndef __PFJET_INFO_H_
#define __PFJET_INFO_H_

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
//#include "CMGTools/External/interface/PileupJetIdentifier.h"
//#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "FWCore/Utilities/interface/isFinite.h"

using namespace btagbtvdeep;

class pfjetInfo : public baseTree{

 public:
  pfjetInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~pfjetInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();

 private:
  pfjetInfo(){};
  edm::InputTag pfjetLabel_;

  //variables which would become branches
  int PFJet_n;
  TClonesArray *PFJet_4Momentum, *PFJet_Vposition, *UcPFJet_4Momentum;

  std::vector<float> PFJet_CEMF;    		       
  std::vector<float> PFJet_NEMF;  		       
  std::vector<float> PFJet_CHF;  		       
  std::vector<float> PFJet_NHF;  		       
  std::vector<float> PFJet_MUF;		       
  std::vector<int> PFJet_NumNeutralParticles;		       
  std::vector<int>   PFJet_CHM;  		       
  std::vector<int>   PFJet_NumConst;            
  std::vector<float> PFJet_jecUncer		    ;     
  std::vector<float> PFJet_jecCorr                  ;
};

#endif

