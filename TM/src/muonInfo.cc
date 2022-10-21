#include "TreeMaker/TM/interface/muonInfo.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TreeMaker/TM/interface/rhoInfo.h"

muonInfo::muonInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in muon constructor"<<std::endl;
  muonLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("muonLabel_");
  Muon_4Momentum                      = new TClonesArray("TLorentzVector");
  if(debug) std::cout<<"in muon constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

muonInfo::~muonInfo(){
  delete tree_;
  delete Muon_4Momentum;
}

void muonInfo::Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx, float& rhoLepton){
  Clear();
  if(debug_)    std::cout<<"getting muon info"<<std::endl;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;

  if(not iEvent.getByLabel(muonLabel_,muonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<muonLabel_<<std::endl; 
    exit(0);
  }  

  edm::View<pat::Muon>::const_iterator muon;
  for(muon = muons.begin(); muon!=muons.end(); muon++){
    idPreselection = 0.;
    TLorentzVector p4(muon->px(),muon->py(),muon->pz(),muon->energy());
    new( (*Muon_4Momentum)[Muon_n]) TLorentzVector(p4);
    Muon_isLooseMuon.push_back((bool)muon->isLooseMuon());

    if(muon->innerTrack().isNonnull()){
      dxy = muon->innerTrack()->dxy(pv);
      dz = muon->innerTrack()->dz(pv);
    }

    // Find isolation (Fall17)
    pfIsoCharged = muon->pfIsolationR03().sumChargedHadronPt;
    pfIsoNeutral = muon->pfIsolationR03().sumNeutralHadronEt + muon->pfIsolationR03().sumPhotonEt;
    if (abs(muon->eta()) < 0.8) EffArea = 0.0566;
    else if (abs(muon->eta()) < 1.3) EffArea = 0.0562;
    else if (abs(muon->eta()) < 2.0) EffArea = 0.0363;
    else if (abs(muon->eta()) < 2.2) EffArea = 0.0119;
    else EffArea = 0.0064;
    correction = rhoLepton*EffArea;
    pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
    relIso = (pfIsoCharged + pfIsoPUSubtracted)/muon->pt();

    if(muon->pt() > 10 && abs(muon->eta()) < 2.4 && abs(dxy) < 0.05 && abs(dz) < 0.1 && relIso < 0.4 && muon->isLooseMuon()) idPreselection = 1.;
    Muon_idPreselection.push_back(idPreselection);
    Muon_n++;
  }  
  if(debug_)    std::cout<<"got muon info"<<std::endl;
}

void muonInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Muon_n  ,"Muon_n");
  AddBranch(&Muon_4Momentum ,"Muon_4Momentum");
  AddBranch(&Muon_isLooseMuon         ,"Muon_isLooseMuon");
  AddBranch(&Muon_idPreselection      ,"Muon_idPreselection");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void muonInfo::Clear(){
  if(debug_)std::cout<<"clearing Muon info"<<std::endl;
  Muon_n = 0; 
  idPreselection = -99.;
  Muon_4Momentum->Clear();
  Muon_isLooseMuon.clear();
  Muon_idPreselection.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
