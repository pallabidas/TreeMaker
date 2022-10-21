#include "TreeMaker/TM/interface/photonInfo.h"

photonInfo::photonInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in photon constructor"<<std::endl;
  photonLabel_               = iConfig.getUntrackedParameter<edm::InputTag>("photonLabel_");
  Photon_4Momentum           = new TClonesArray("TLorentzVector");
  Photon_4Momentum_corr      = new TClonesArray("TLorentzVector");
  Photon_Vposition           = new TClonesArray("TVector3");

  if(debug) std::cout<<"in photon constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

photonInfo::~photonInfo(){
  delete tree_;
  delete Photon_4Momentum         ;
  delete Photon_4Momentum_corr    ;
  delete Photon_Vposition         ;
}

void photonInfo::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  if(debug_)    std::cout<<"getting photon info"<<std::endl;
     
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByLabel(photonLabel_,photonHandle); 
  const edm::View<pat::Photon> & photons = *photonHandle;

  if(not iEvent.getByLabel(photonLabel_,photonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<photonLabel_<<std::endl;
    exit(0);
  }

  edm::View<pat::Photon>::const_iterator photon;

  for(photon=photons.begin(); photon!=photons.end(); photon++){
   
    double corrfactor = photon->userFloat("ecalEnergyPostCorr")/photon->energy();

    TLorentzVector p4(photon->px(),photon->py(),photon->pz(),photon->energy());
    TLorentzVector corrp4 = p4*corrfactor;

    new( (*Photon_4Momentum)[Photon_n]) TLorentzVector(p4);
    new( (*Photon_4Momentum_corr)[Photon_n]) TLorentzVector(corrp4);
    TVector3 v3(photon->vx(),photon->vy(),photon->vz());
    new( (*Photon_Vposition)[Photon_n]) TVector3(v3);
    Photon_HoverE.push_back((float) photon->hadronicOverEm());
    Photon_sigmaIetaIeta.push_back((float) photon->sigmaIetaIeta());
    Photon_chargedHadronIso.push_back((float) photon->chargedHadronIso());
    Photon_neutralHadronIso.push_back((float) photon->neutralHadronIso());
    Photon_photonIso.push_back((float) photon->photonIso());
    Photon_r9.push_back((float) photon->r9());
    Photon_hasPixelSeed.push_back((bool) photon->hasPixelSeed());
    Photon_passElectronVeto.push_back((bool) photon->passElectronVeto());
    Photon_hOVERe.push_back((float) photon->hadTowOverEm());
    Photon_full5x5_r9.push_back((float) photon->full5x5_r9());
    Photon_full5x5_sigmaIetaIeta.push_back((float) photon->full5x5_sigmaIetaIeta());
    Photon_full5x5_e5x5.push_back((float) photon->full5x5_e5x5());
    Photon_scEta.push_back((float) photon->superCluster()->eta());
    Photon_scEnergy.push_back((float) photon->superCluster()->energy());
 
    Photon_n++;	
  }//end of photonloop
  
  
  if(debug_)    std::cout<<"got photon info"<<std::endl;
}

void photonInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;

  AddBranch(&Photon_n, "Photon_n"); 
  AddBranch(&Photon_HoverE, "Photon_HoverE");
  AddBranch(&Photon_sigmaIetaIeta, "Photon_sigmaIetaIeta");
  AddBranch(&Photon_chargedHadronIso, "Photon_chargedHadronIso");
  AddBranch(&Photon_neutralHadronIso, "Photon_neutralHadronIso");
  AddBranch(&Photon_photonIso, "Photon_photonIso");
  AddBranch(&Photon_4Momentum , "Photon_4Momentum");
  AddBranch(&Photon_4Momentum_corr , "Photon_4Momentum_corr");
  AddBranch(&Photon_Vposition , "Photon_Vposition");
  AddBranch(&Photon_r9, "Photon_r9");
  AddBranch(&Photon_hasPixelSeed, "Photon_hasPixelSeed");
  AddBranch(&Photon_passElectronVeto, "Photon_passElectronVeto");
  AddBranch(&Photon_hOVERe, "Photon_hOVERe");
  AddBranch(&Photon_full5x5_r9, "Photon_full5x5_r9");
  AddBranch(&Photon_full5x5_sigmaIetaIeta, "Photon_full5x5_sigmaIetaIeta");
  AddBranch(&Photon_full5x5_e5x5, "Photon_full5x5_e5x5");
  AddBranch(&Photon_scEta, "Photon_scEta");
  AddBranch(&Photon_scEnergy, "Photon_scEnergy");
     
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void photonInfo::Clear(){
  if(debug_)std::cout<<"clearing Photon info"<<std::endl;
  Photon_n =0;
  Photon_4Momentum->Clear();
  Photon_4Momentum_corr->Clear();
  Photon_Vposition->Clear();
  Photon_HoverE.clear();
  Photon_sigmaIetaIeta.clear();
  Photon_chargedHadronIso.clear();
  Photon_neutralHadronIso.clear();
  Photon_photonIso.clear();
  Photon_r9.clear();
  Photon_hasPixelSeed.clear();
  Photon_passElectronVeto.clear();
  Photon_hOVERe.clear();
  Photon_full5x5_r9.clear();
  Photon_full5x5_sigmaIetaIeta.clear();
  Photon_full5x5_e5x5.clear();
  Photon_scEta.clear();
  Photon_scEnergy.clear();

  if(debug_) std::cout<<"cleared"<<std::endl;
}

