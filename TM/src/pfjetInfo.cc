#include "TreeMaker/TM/interface/pfjetInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

pfjetInfo::pfjetInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in pfjet constructor"<<std::endl;
  pfjetLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("pfjetLabel_");
  genpartLabel_   = iConfig.getUntrackedParameter<edm::InputTag> ("genParticleLabel_");
  genjetLabel_    = iConfig.getUntrackedParameter<edm::InputTag> ("genJetLabel_");
  //PFJet_4Momentum   = new TClonesArray("TLorentzVector");
  //UcPFJet_4Momentum = new TClonesArray("TLorentzVector");
  //PFJet_Vposition   = new TClonesArray("TVector3");
  if(debug) std::cout<<"in pfjet constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

pfjetInfo::~pfjetInfo(){
  delete tree_;
  //delete PFJet_4Momentum   ;
  //delete UcPFJet_4Momentum ;
  //delete PFJet_Vposition   ;
}

void pfjetInfo::Fill(const edm::Event& iEvent,const  edm::EventSetup& iSetup){
  Clear();
  if(debug_)    std::cout<<"getting pfjet info"<<std::endl;
  //for jec uncertainty
  //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *pfjecUnc = new JetCorrectionUncertainty(JetCorPar);

  edm::Handle<edm::View<pat::Jet> > pfjetHandle;
  iEvent.getByLabel(pfjetLabel_,pfjetHandle);
  const edm::View<pat::Jet> & pfjets = *pfjetHandle;
  if(not iEvent.getByLabel(pfjetLabel_,pfjetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<pfjetLabel_<<std::endl; 
    exit(0);
  }

  edm::Handle<edm::View<reco::GenJet> > genjetHandle;
  iEvent.getByLabel(genjetLabel_,genjetHandle);
  const edm::View<reco::GenJet> & genjets = *genjetHandle;
  if(not iEvent.getByLabel(genjetLabel_,genjetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<genjetLabel_<<std::endl;
    exit(0);
  }

  edm::View<pat::Jet>::const_iterator pfjet;
  for(pfjet = pfjets.begin(); pfjet != pfjets.end(); ++pfjet){
    if(pfjet->pt() < 20.) continue; 
    float deltar = 0.4;
    for(auto genjet = genjets.begin(); genjet != genjets.end(); ++genjet){
      if(reco::deltaR(*pfjet, *genjet) < deltar){
        deltar = reco::deltaR(*pfjet, *genjet);
        jet_genmatch_pt = genjet->pt();
        jet_genmatch_eta = genjet->eta();
        jet_genmatch_phi = genjet->phi();
        jet_genmatch_mass = genjet->mass();
      }
    }
    jet_pt = pfjet->pt();
    jet_eta = pfjet->eta();
    jet_phi = pfjet->phi();
    jet_mass = pfjet->mass();
    jet_ncand = pfjet->neutralMultiplicity() + pfjet->chargedMultiplicity();
    jet_nel = pfjet->electronMultiplicity();
    jet_nmu = pfjet->muonMultiplicity();
    jet_nbhad = pfjet->jetFlavourInfo().getbHadrons().size();
    jet_nchad = pfjet->jetFlavourInfo().getcHadrons().size();
    jet_hflav = pfjet->jetFlavourInfo().getHadronFlavour(); 
    jet_pflav = pfjet->jetFlavourInfo().getPartonFlavour();
    jet_lepflav = pfjet->jetFlavourInfo().getLeptons().size();

    //jet_elflav = 0; jet_muflav = 0; jet_tauflav = 0;
    //const reco::GenParticleRefVector& leptons = pfjet->jetFlavourInfo().getLeptons();
    //std::cout << "                      # of clustered leptons: " << leptons.size() << std::endl;
    //for (reco::GenParticleRefVector::const_iterator it = leptons.begin(); it != leptons.end(); ++it) {
    //    std::cout<<"inside loop"<<std::endl;
    //    std::cout<<(*it)->pdgId()<<std::endl;
    //  if(abs((*it)->pdgId() == 11)) jet_elflav++;
    //  if(abs((*it)->pdgId() == 13)) jet_muflav++;
    //  if(abs((*it)->pdgId() == 15)) jet_tauflav++;
    //}
    //std::cout << "jet_elflav: "<<jet_elflav<<"\t"<<"jet_muflav: "<<jet_muflav<<"\t"<<"jet_tauflav: "<<jet_tauflav<<std::endl;

    int nConstituents = pfjet->numberOfSourceCandidatePtrs();
    for(int k = 0; k < nConstituents; k++){
      reco::CandidatePtr pfcand = pfjet->sourceCandidatePtr(k);
      jet_pfcand_pt.push_back(pfcand->pt());
      jet_pfcand_eta.push_back(pfcand->eta());
      jet_pfcand_phi.push_back(pfcand->phi());
    }

/*
    TLorentzVector p4(pfjet->px(),pfjet->py(),pfjet->pz(),pfjet->energy());
    new( (*PFJet_4Momentum)[PFJet_n]) TLorentzVector(p4);
    TVector3 v3(pfjet->vx(),pfjet->vy(),pfjet->vz());
    new( (*PFJet_Vposition)[PFJet_n]) TVector3(v3);
    
    if(pfjet->jecFactor("Uncorrected")!= 0)
      PFJet_jecCorr.push_back(1.0/pfjet->jecFactor("Uncorrected")); 
    else PFJet_jecCorr.push_back(0.);
    
    PFJet_CEMF.push_back((float)pfjet->chargedEmEnergyFraction());
    PFJet_NEMF.push_back((float)pfjet->neutralEmEnergyFraction());
    PFJet_CHF.push_back((float)pfjet->chargedHadronEnergyFraction());
    PFJet_NHF.push_back((float)pfjet->neutralHadronEnergyFraction());
    PFJet_MUF.push_back((float)pfjet->muonEnergyFraction());
    PFJet_NumNeutralParticles.push_back((int)pfjet->neutralMultiplicity());
    PFJet_CHM.push_back((int)pfjet->chargedMultiplicity());
    PFJet_NumConst.push_back((int)pfjet->neutralMultiplicity()+(int)pfjet->chargedMultiplicity());
*/
     
    //jet energy uncertiany
    //pfjecUnc->setJetEta(pfjet->eta());
    //pfjecUnc->setJetPt(pfjet->pt());
    //double unc = pfjecUnc->getUncertainty(true);
    //PFJet_jecUncer.push_back(unc);
    
    //get the uncorrected jet and fill them
    //pat::Jet uncpfjet = pfjet->correctedJet("Uncorrected");
    //TLorentzVector ucp4(uncpfjet.px(),uncpfjet.py(),uncpfjet.pz(),uncpfjet.energy());
    //new( (*UcPFJet_4Momentum)[PFJet_n]) TLorentzVector(ucp4);
    
    //PFJet_n++;
    if(debug_)    std::cout<<"filling pfjet info"<<std::endl;
    tree_->Fill();
    Clear();
  }//end of for loop
  if(debug_)    std::cout<<"got pfjet info"<<std::endl;
}

void pfjetInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&jet_pt, "jet_pt");
  AddBranch(&jet_eta, "jet_eta");
  AddBranch(&jet_phi, "jet_phi");
  AddBranch(&jet_mass, "jet_mass");
  AddBranch(&jet_ncand, "jet_ncand");
  AddBranch(&jet_nel, "jet_nel");
  AddBranch(&jet_nmu, "jet_nmu");
  AddBranch(&jet_nbhad, "jet_nbhad");
  AddBranch(&jet_nchad, "jet_nchad");
  AddBranch(&jet_hflav, "jet_hflav");
  AddBranch(&jet_pflav, "jet_pflav");
  AddBranch(&jet_lepflav, "jet_lepflav");
  AddBranch(&jet_genmatch_pt, "jet_genmatch_pt");
  AddBranch(&jet_genmatch_eta, "jet_genmatch_eta");
  AddBranch(&jet_genmatch_phi, "jet_genmatch_phi");
  AddBranch(&jet_genmatch_mass, "jet_genmatch_mass");

  AddBranch(&jet_pfcand_pt , "jet_pfcand_pt");
  AddBranch(&jet_pfcand_eta , "jet_pfcand_eta");
  AddBranch(&jet_pfcand_phi , "jet_pfcand_phi");

/*
  AddBranch(&PFJet_n                     ,"PFJet_n");
  AddBranch(&PFJet_4Momentum             ,"PFJet_4Momentum");
  AddBranch(&UcPFJet_4Momentum           ,"UcPFJet_4Momentum");
  AddBranch(&PFJet_Vposition             ,"PFJet_Vposition");
  AddBranch(&PFJet_CEMF                  ,"PFJet_CEMF");
  AddBranch(&PFJet_NEMF                  ,"PFJet_NEMF");
  AddBranch(&PFJet_CHF                   ,"PFJet_CHF");
  AddBranch(&PFJet_NHF                   ,"PFJet_NHF");
  AddBranch(&PFJet_MUF                   ,"PFJet_MUF");
  AddBranch(&PFJet_NumNeutralParticles   ,"PFJet_NumNeutralParticles");
  AddBranch(&PFJet_CHM                   ,"PFJet_CHM");
  AddBranch(&PFJet_NumConst              ,"PFJet_NumConst");
  AddBranch(&PFJet_jecUncer              ,"PFJet_jecUncer");
  AddBranch(&PFJet_jecCorr               ,"PFJet_jecCorr");
*/

  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void pfjetInfo::Clear(){
  if(debug_)std::cout<<"clearing PFjet info"<<std::endl;
  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_mass = 0;
  jet_ncand = 0;
  jet_nel = 0;
  jet_nmu = 0;
  jet_nbhad = 0;
  jet_nchad = 0;
  jet_hflav = 0;
  jet_pflav = 0;
  jet_lepflav = 0;
  jet_genmatch_pt = 0;
  jet_genmatch_eta = 0;
  jet_genmatch_phi = 0;
  jet_genmatch_mass = 0;

  jet_pfcand_pt.clear();
  jet_pfcand_eta.clear();
  jet_pfcand_phi.clear();
/*
  PFJet_n = 0;  
  PFJet_4Momentum->Clear();
  UcPFJet_4Momentum->Clear();
  PFJet_Vposition->Clear();
  PFJet_CEMF.clear();    
  PFJet_NEMF.clear();
  PFJet_CHF.clear();
  PFJet_NHF.clear();
  PFJet_MUF.clear();
  PFJet_NumNeutralParticles.clear();
  PFJet_CHM.clear();
  PFJet_NumConst.clear();         
  PFJet_jecUncer.clear();
  PFJet_jecCorr.clear();           
*/
  
  if(debug_) std::cout<<"cleared"<<std::endl;
}
