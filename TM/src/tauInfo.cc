#include "TreeMaker/TM/interface/tauInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


tauInfo::tauInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in tau constructor"<<std::endl;
  tauLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("tauLabel_");
  Tau_4Momentum                      = new TClonesArray("TLorentzVector");
  if(debug) std::cout<<"in tau constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

tauInfo::~tauInfo(){
  delete tree_;
  delete Tau_4Momentum                      ;
}

void tauInfo::Fill(const edm::Event& iEvent, math::XYZPoint& pv, reco::Vertex& vtx){
  Clear();
  if(debug_)    std::cout<<"getting tau info"<<std::endl;
  float dxy_approx = -99.;
  float dxy = -99.;
  float dz = -99.;
  float z_impact = -99.;
  float dxy_new = -99.;
  float dz_new = -99.;

  edm::Handle<edm::View<pat::Tau>> tauHandle;
  iEvent.getByLabel(tauLabel_, tauHandle);
  const edm::View<pat::Tau> & taus = *tauHandle;

  if(not iEvent.getByLabel(tauLabel_, tauHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<tauLabel_<<std::endl; 
    exit(0);
  }  
  edm::View<pat::Tau>::const_iterator tau;
    for(tau = taus.begin(); tau !=taus.end(); tau++){
	if(!tau->tauID("decayModeFinding")) continue;
	TLorentzVector p4(tau->px(),tau->py(),tau->pz(),tau->energy());
	new( (*Tau_4Momentum)[Tau_n]) TLorentzVector(p4);
	math::XYZPoint tauVtx(tau->leadChargedHadrCand()->vertex().x(), tau->leadChargedHadrCand()->vertex().y(), tau->leadChargedHadrCand()->vertex().z());
	dxy_approx = (-(tauVtx.x()-pv.x())*(tau->py())+(tauVtx.y()-pv.y())*(tau->px()))/p4.Pt();
	if(abs(pv.z()-tau->vertex().z()) < 0.0001){
		dxy = tau->dxy();
	}
	else{
		dxy = dxy_approx;
	}
	dz = (tauVtx.z()-pv.z()) - ((tauVtx.x()-pv.x())*p4.X()+(tauVtx.y()-pv.y())*p4.Y())/ p4.Pt() * p4.Z()/ p4.Pt();
	z_impact = pv.z()+130./(tan(tau->theta()));

	auto packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());
	dz_new = packedLeadTauCand->dz();
	dxy_new = packedLeadTauCand->dxy();
	
	Tau_charge.push_back((int) tau->charge());
	Tau_decayModeNewDM.push_back((bool) tau->tauID("decayModeFindingNewDMs"));
	Tau_dxy.push_back((float) dxy);
	Tau_dz.push_back((float) dz);
	Tau_dxy_new.push_back((float) dxy_new);
	Tau_dz_new.push_back((float) dz_new);
	Tau_pvDiff.push_back((float) abs(pv.z()-tau->vertex().z()));
	Tau_z_impact.push_back((float) z_impact);

	Tau_chargedIsoPtSum.push_back((float) tau->tauID("chargedIsoPtSum"));
	Tau_chargedIsoPtSumdR03.push_back((float) tau->tauID("chargedIsoPtSumdR03"));
	Tau_footprintCorrection.push_back((float) tau->tauID("footprintCorrection"));
	Tau_footprintCorrectiondR03.push_back((float) tau->tauID("footprintCorrectiondR03"));
	Tau_neutralIsoPtSum.push_back((float) tau->tauID("neutralIsoPtSum"));
	Tau_neutralIsoPtSumWeight.push_back((float) tau->tauID("neutralIsoPtSumWeight"));
	Tau_neutralIsoPtSumWeightdR03.push_back((float) tau->tauID("neutralIsoPtSumWeightdR03"));
	Tau_neutralIsoPtSumdR03.push_back((float) tau->tauID("neutralIsoPtSumdR03"));
	Tau_photonPtSumOutsideSignalCone.push_back((float) tau->tauID("photonPtSumOutsideSignalCone"));
	Tau_photonPtSumOutsideSignalConedR03.push_back((float) tau->tauID("photonPtSumOutsideSignalConedR03"));
	Tau_puCorrPtSum.push_back((float) tau->tauID("puCorrPtSum"));	
	Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw.push_back((float) tau->tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw"));
	Tau_byIsolationMVArun2v1DBnewDMwLTraw.push_back((float) tau->tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
	Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back((float) tau->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
	Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back((float) tau->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
	Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back((float) tau->tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
	Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back((float) tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
	Tau_againstMuonLoose3.push_back((float) tau->tauID("againstMuonLoose3"));
	Tau_againstMuonTight3.push_back((float) tau->tauID("againstMuonTight3"));
	Tau_againstElectronLooseMVA6.push_back((float) tau->tauID("againstElectronLooseMVA6"));
	Tau_againstElectronMediumMVA6.push_back((float) tau->tauID("againstElectronMediumMVA6"));
	Tau_againstElectronTightMVA6.push_back((float) tau->tauID("againstElectronTightMVA6"));	

	Tau_n++;
  }  
  if(debug_)    std::cout<<"got tau info"<<std::endl;
}

void tauInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Tau_n  ,"Tau_n");
  AddBranch(&Tau_4Momentum ,"Tau_4Momentum");
  AddBranch(&Tau_charge, "Tau_charge"); 
  AddBranch(&Tau_decayModeNewDM, "Tau_decayModeNewDM");
  AddBranch(&Tau_dxy, "Tau_dxy");
  AddBranch(&Tau_dz, "Tau_dz");
  AddBranch(&Tau_dxy_new, "Tau_dxy_new");
  AddBranch(&Tau_dz_new, "Tau_dz_new");
  AddBranch(&Tau_pvDiff, "Tau_pvDiff");
  AddBranch(&Tau_z_impact, "Tau_z_impact");

  AddBranch(&Tau_chargedIsoPtSum, "Tau_chargedIsoPtSum");
  AddBranch(&Tau_chargedIsoPtSumdR03, "Tau_chargedIsoPtSumdR03");
  AddBranch(&Tau_footprintCorrection, "Tau_footprintCorrection");
  AddBranch(&Tau_footprintCorrectiondR03, "Tau_footprintCorrectiondR03");
  AddBranch(&Tau_neutralIsoPtSum, "Tau_neutralIsoPtSum");
  AddBranch(&Tau_neutralIsoPtSumWeight, "Tau_neutralIsoPtSumWeight");
  AddBranch(&Tau_neutralIsoPtSumWeightdR03, "Tau_neutralIsoPtSumWeightdR03");
  AddBranch(&Tau_neutralIsoPtSumdR03, "Tau_neutralIsoPtSumdR03");
  AddBranch(&Tau_photonPtSumOutsideSignalCone, "Tau_photonPtSumOutsideSignalCone");
  AddBranch(&Tau_photonPtSumOutsideSignalConedR03, "Tau_photonPtSumOutsideSignalConedR03");
  AddBranch(&Tau_puCorrPtSum, "Tau_puCorrPtSum");
  AddBranch(&Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw,"Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw");
  AddBranch(&Tau_byIsolationMVArun2v1DBnewDMwLTraw, "Tau_byIsolationMVArun2v1DBnewDMwLTraw");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, "Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, "Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, "Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, "Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_againstMuonLoose3, "Tau_againstMuonLoose3");
  AddBranch(&Tau_againstMuonTight3, "Tau_againstMuonTight3");
  AddBranch(&Tau_againstElectronLooseMVA6, "Tau_againstElectronLooseMVA6");
  AddBranch(&Tau_againstElectronMediumMVA6, "Tau_againstElectronMediumMVA6");
  AddBranch(&Tau_againstElectronTightMVA6, "Tau_againstElectronTightMVA6");

  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void tauInfo::Clear(){
  if(debug_)std::cout<<"clearing Tau info"<<std::endl;
  Tau_n =0;

  Tau_4Momentum->Clear();
  Tau_charge.clear();
  Tau_decayModeNewDM.clear();
  Tau_dxy.clear();
  Tau_dz.clear();
  Tau_dxy_new.clear();
  Tau_dz_new.clear();
  Tau_pvDiff.clear();
  Tau_z_impact.clear();

  Tau_chargedIsoPtSum.clear();
  Tau_chargedIsoPtSumdR03.clear();
  Tau_footprintCorrection.clear();
  Tau_footprintCorrectiondR03.clear();
  Tau_neutralIsoPtSum.clear();
  Tau_neutralIsoPtSumWeight.clear();
  Tau_neutralIsoPtSumWeightdR03.clear();
  Tau_neutralIsoPtSumdR03.clear();
  Tau_photonPtSumOutsideSignalCone.clear();
  Tau_photonPtSumOutsideSignalConedR03.clear();
  Tau_puCorrPtSum.clear();
  Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw.clear();
  Tau_byIsolationMVArun2v1DBnewDMwLTraw.clear();
  Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_againstMuonLoose3.clear();
  Tau_againstMuonTight3.clear();
  Tau_againstElectronLooseMVA6.clear();
  Tau_againstElectronMediumMVA6.clear();
  Tau_againstElectronTightMVA6.clear();

 
 if(debug_) std::cout<<"cleared"<<std::endl;
}


