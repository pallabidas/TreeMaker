#ifndef  TM_H
#define  TM_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "TreeMaker/TM/interface/genparticleInfo.h"
#include "TreeMaker/TM/interface/eventInfo.h"
#include "TreeMaker/TM/interface/pileUpInfo.h"
#include "TreeMaker/TM/interface/triggerInfo.h"
#include "TreeMaker/TM/interface/vertexInfo.h"
#include "TreeMaker/TM/interface/muonInfo.h"
#include "TreeMaker/TM/interface/photonInfo.h"
#include "TreeMaker/TM/interface/electronInfo.h"
#include "TreeMaker/TM/interface/rhoInfo.h"
#include "TreeMaker/TM/interface/pfmetInfo.h"
#include "TreeMaker/TM/interface/pfjetInfo.h"
#include "TreeMaker/TM/interface/tauInfo.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<std::string> >+;
#pragma link C++ class std::vector<std::vector<TString> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<bool> >+;
#pragma extra_include "std::vector";
#endif

class TM : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit TM(const edm::ParameterSet&);
      ~TM();


   private:
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run& , const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      TFile* file;
      TTree* tree_;
      const size_t MaxN;

      bool debug_;
      bool fillPhotInfo_;
      bool fillgenparticleInfo_;
      bool filleventInfo_;
      bool fillpileUpInfo_;
      bool filltriggerInfo_;
      bool fillvertexInfo_;
      bool fillmuonInfo_;
      bool fillelectronInfo_;
      bool fillrhoInfo_;
      bool fillpfmetInfo_;
      bool fillpfjetInfo_;
      bool filltauInfo_;

      genparticleInfo *genparticleInfo_;
      photonInfo      *photonInfo_;
      eventInfo       *eventInfo_;
      pileUpInfo      *pileUpInfo_;
      triggerInfo     *triggerInfo_;
      vertexInfo      *vertexInfo_;
      muonInfo        *muonInfo_;
      electronInfo    *electronInfo_;
      rhoInfo         *rhoInfo_;
      pfmetInfo       *pfmetInfo_;
      pfjetInfo       *pfjetInfo_;
      tauInfo         *tauInfo_;

      float rhoLepton; 
      HLTConfigProvider hltConfig_;
      HLTPrescaleProvider hltPrescale_;
      std::string hltlabel_;
      std::vector<std::string> all_triggers;

      edm::EDGetTokenT<reco::VertexCollection> vtxToken;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>>  puToken;
      edm::EDGetTokenT<edm::View<pat::Muon>> muonToken; 
      edm::EDGetTokenT<pat::ElectronCollection> eleToken;
      edm::EDGetTokenT<edm::View<pat::MET>> metToken;
      edm::EDGetTokenT<edm::View<pat::MET>> puppimetToken;
      edm::EDGetTokenT<edm::View<pat::Jet>> jetToken;
      //edm::EDGetTokenT<JetCorrectorParametersCollection> jecToken;
      edm::EDGetTokenT<edm::View<pat::Photon>> phoToken;
      edm::EDGetTokenT<edm::TriggerResults> TRToken;
      edm::EDGetTokenT<edm::View<pat::Tau>> tauToken;
      edm::EDGetTokenT<reco::GenParticleCollection> genToken;
      edm::EDGetTokenT<GenEventInfoProduct> genEvtToken;
      edm::EDGetTokenT<LHEEventProduct> generatorlheToken;
      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<bool> ecalBadCalib_token;
      edm::EDGetTokenT<edm::TriggerResults> metfilterspatToken_;
      edm::EDGetTokenT<edm::TriggerResults> metfiltersrecoToken_;

      //edm::ESGetToken<JetCorrector, JetCorrectionsRecord> jetCorrector_;

};

#endif
