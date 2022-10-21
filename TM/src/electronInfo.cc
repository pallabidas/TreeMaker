#include "TreeMaker/TM/interface/electronInfo.h"

electronInfo::electronInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in electron constructor"<<std::endl;
  electronLabel_               = iConfig.getUntrackedParameter<edm::InputTag> ("electronLabel_");
  Electron_4Momentum           = new TClonesArray("TLorentzVector");
  if(debug) std::cout<<"in electron constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

electronInfo::~electronInfo(){
  delete tree_;
  delete Electron_4Momentum  ; 
}

void electronInfo::Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx, float& rhoLepton){
  Clear();
  if(debug_)    std::cout<<"getting electron info"<<std::endl;

  edm::Handle<std::vector<pat::Electron> > electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);
  if(not iEvent.getByLabel(electronLabel_,electronHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<electronLabel_<<std::endl; 
    exit(0);
  }  
  
  pat::ElectronCollection::const_iterator electron;
  for(electron=electronHandle->begin(); electron!=electronHandle->end(); electron++){
    idPreselection = 0.;
    TLorentzVector p4(electron->px(),electron->py(),electron->pz(),electron->energy());
    new( (*Electron_4Momentum)[Electron_n]) TLorentzVector(p4);

    if(electron->superCluster().isNonnull()) scEta = abs(electron->superCluster()->eta());
    dEtaInSeed = electron->superCluster().isNonnull() && electron->superCluster()->seed().isNonnull() ? electron->deltaEtaSuperClusterTrackAtVtx() - electron->superCluster()->eta() + electron->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
    if(electron->gsfTrack().isNonnull()){
      dxy = electron->gsfTrack()->dxy(pv);
      dz = electron->gsfTrack()->dz(pv);
      numMissingHits = electron->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    // Find isolation (Fall17)
    pfIsoCharged = electron->pfIsolationVariables().sumChargedHadronPt;
    pfIsoNeutral = electron->pfIsolationVariables().sumNeutralHadronEt + electron->pfIsolationVariables().sumPhotonEt;
    if (abs(electron->eta()) >= 0. && abs(electron->eta()) < 1.0) EffArea = 0.1440;
    else if (abs(electron->eta()) >= 1.0 && abs(electron->eta()) < 1.479) EffArea = 0.1562;
    else if (abs(electron->eta()) >= 1.479 && abs(electron->eta()) < 2.0) EffArea = 0.1032;
    else if (abs(electron->eta()) >= 2.0 && abs(electron->eta()) < 2.2) EffArea = 0.0859;
    else if (abs(electron->eta()) >= 2.2 && abs(electron->eta()) < 2.3) EffArea = 0.1116;
    else if (abs(electron->eta()) >= 2.3 && abs(electron->eta()) < 2.4) EffArea = 0.1321;
    else if (abs(electron->eta()) >= 2.4 && abs(electron->eta()) < 2.5) EffArea = 0.1654;
    correction = rhoLepton*EffArea;
    pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
    relIso = (pfIsoCharged + pfIsoPUSubtracted)/electron->pt();

    // Selection of veto electrions (Fall17)
    if (electron->pt() > 10 && scEta < 1.479) {
      idPreselection = electron->full5x5_sigmaIetaIeta() < 0.0126 &&
         abs(dEtaInSeed) < 0.00463 &&
         abs(electron->deltaPhiSuperClusterTrackAtVtx()) < 0.148 &&
         electron->hcalOverEcal() < 0.05 + (1.16/electron->ecalEnergy()) + (0.0324*rhoLepton/electron->ecalEnergy()) &&
         relIso < 0.198+(0.506/electron->pt()) &&
         abs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()) < 0.209 &&
         numMissingHits <= 2 &&
         electron->passConversionVeto() &&
         abs(dxy) < 0.05 &&
         abs(dz) < 0.10;
    }
    else if (electron->pt() > 10 && scEta < 2.5) {
      idPreselection = electron->full5x5_sigmaIetaIeta() < 0.0457 &&
         abs(dEtaInSeed) <  0.00814 &&
         abs(electron->deltaPhiSuperClusterTrackAtVtx()) < 0.19 &&
         electron->hcalOverEcal() < 0.05 + (2.54/electron->ecalEnergy()) + (0.183*rhoLepton/electron->ecalEnergy()) &&
         relIso < 0.203+(0.963/electron->pt()) &&
         abs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()) < 0.132 &&
         numMissingHits <= 3 &&
         electron->passConversionVeto() &&
         abs(dxy) < 0.10 &&
         abs(dz) < 0.20;
    }
    Electron_idPreselection.push_back(idPreselection);
    Electron_n++;
  }
  if(debug_)    std::cout<<"got electron info"<<std::endl;
}

void electronInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Electron_n  ,"Electron_n");
  AddBranch(&Electron_4Momentum   ,"Electron_4Momentum");
  AddBranch(&Electron_idPreselection      ,"Electron_idPreselection");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void electronInfo::Clear(){
  if(debug_)std::cout<<"clearing Electron info"<<std::endl;
  Electron_n = 0;
  idPreselection = -99.;
  Electron_4Momentum->Clear();
  Electron_idPreselection.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
