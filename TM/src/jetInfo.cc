#include "TreeMaker/TM/interface/jetInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

jetInfo::jetInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in jet constructor"<<std::endl;
  jetLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("jetLabel_");
  genpartLabel_   = iConfig.getUntrackedParameter<edm::InputTag> ("genParticleLabel_");
  genjetLabel_    = iConfig.getUntrackedParameter<edm::InputTag> ("genJetLabel_");
  vtxLabel_       = iConfig.getUntrackedParameter<edm::InputTag> ("vertexLabel_");
  svLabel_        = iConfig.getUntrackedParameter<edm::InputTag> ("svLabel_");
  if(debug) std::cout<<"in jet constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

jetInfo::~jetInfo(){
  delete tree_;
}

void jetInfo::Fill(const edm::Event& iEvent,const  edm::EventSetup& iSetup){
  Clear();
  if(debug_)    std::cout<<"getting jet info"<<std::endl;

  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_, jetHandle);
  const edm::View<pat::Jet> & jets = *jetHandle;
  if(not iEvent.getByLabel(jetLabel_, jetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<jetLabel_<<std::endl; 
    exit(0);
  }

  edm::Handle<edm::View<reco::GenJet> > genjetHandle;
  iEvent.getByLabel(genjetLabel_, genjetHandle);
  const edm::View<reco::GenJet> & genjets = *genjetHandle;
  if(not iEvent.getByLabel(genjetLabel_, genjetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<genjetLabel_<<std::endl;
    exit(0);
  }

  const reco::Vertex *pv = nullptr;
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByLabel(vtxLabel_, vtxHandle);
  if(not iEvent.getByLabel(vtxLabel_, vtxHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<vtxLabel_<<std::endl;
    exit(0);
  }
  pv = &vtxHandle->at(0);

  edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > svHandle;
  iEvent.getByLabel(svLabel_, svHandle);
  if(not iEvent.getByLabel(svLabel_, svHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<svLabel_<<std::endl;
    exit(0);
  }

  /////// signal jet ///////

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genpartLabel_, genParticles);
  if(not iEvent.getByLabel(genpartLabel_, genParticles )){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<genpartLabel_<<std::endl;
    exit(0);
  }

  TLorentzVector v1;
  TLorentzVector v2;
  TLorentzVector higgs;
  reco::GenParticleCollection genColl(*(genParticles.product()));
  std::sort(genColl.begin(),genColl.end(),PtGreater());
  for(reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
    if(genparticle->pdgId() == 25 && genparticle->status() == 62) higgs.SetPtEtaPhiM(genparticle->pt(), genparticle->eta(), genparticle->phi(), 125.09);
    if(genparticle->pdgId() == 5 && genparticle->status() == 23) v1.SetPtEtaPhiM(genparticle->pt(), genparticle->eta(), genparticle->phi(), 4.18);
    if(genparticle->pdgId() == -5 && genparticle->status() == 23) v2.SetPtEtaPhiM(genparticle->pt(), genparticle->eta(), genparticle->phi(), 4.18);
  }
  float dR_genBs = v1.DeltaR(v2);
  //std::cout<<"dR_genBs: "<<dR_genBs<<std::endl;
  TLorentzVector diB = v1 + v2;

  edm::View<pat::Jet>::const_iterator jet;
  pat::Jet matchedJet;
  float deltar = 0.4;
  for(jet = jets.begin(); jet != jets.end(); ++jet){
    if(reco::deltaR(jet->eta(), jet->phi(), diB.Eta(), diB.Phi()) < deltar){
      deltar = reco::deltaR(jet->eta(), jet->phi(), diB.Eta(), diB.Phi());
      matchedJet = *jet;
    }
    //std::cout<<"deltar: "<<deltar<<std::endl;
  }
  //std::cout<<"matchedJet: "<<matchedJet.pt()<<"\t"<<matchedJet.eta()<<"\t"<<matchedJet.phi()<<std::endl;

  // select signal events where bi-b system is boosted and match to a reco jet
  if(dR_genBs < 0.8 && matchedJet.pt() > 20.){
    deltar = 0.4;
    for(auto genjet = genjets.begin(); genjet != genjets.end(); ++genjet){
      if(reco::deltaR(matchedJet, *genjet) < deltar){
        deltar = reco::deltaR(matchedJet, *genjet);
        jet_genmatch_pt = genjet->pt();
        jet_genmatch_eta = genjet->eta();
        jet_genmatch_phi = genjet->phi();
        jet_genmatch_mass = genjet->mass();
      }
    }
    // refer: https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/plugins/DeepBoostedJetTagInfoProducer.cc
    // find SVs associated with the jet
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    for(const auto &sv : *svHandle){
      if(reco::deltaR(matchedJet, sv) < 0.4){
        jetSVs.push_back(&sv);
      }
    }
    // sort by dxy significance
    std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) { return btagbtvdeep::sv_vertex_comparator(*sva, *svb, *pv); });
    nsv = jetSVs.size();
    float etasign = matchedJet.eta() > 0 ? 1 : -1;

    // fill associated secondar vertices information
    for (const auto *jetsv : jetSVs) {
      jet_sv_pt.push_back(jetsv->pt());
      jet_sv_pt_log.push_back(std::log(jetsv->pt()));
      jet_sv_ptrel.push_back(jetsv->pt() / matchedJet.pt());
      jet_sv_ptrel_log.push_back(std::log(jetsv->pt() / matchedJet.pt()));
      jet_sv_eta.push_back(jetsv->eta());
      jet_sv_phi.push_back(jetsv->phi());
      jet_sv_mass.push_back(jetsv->mass());
      jet_sv_energy.push_back(jetsv->energy());
      jet_sv_energy_log.push_back(std::log(jetsv->energy()));
      jet_sv_erel.push_back(jetsv->energy() / matchedJet.energy());
      jet_sv_erel_log.push_back(std::log(jetsv->energy() / matchedJet.energy()));
      jet_sv_deta.push_back(etasign * (jetsv->eta() - matchedJet.eta()));
      jet_sv_dphi.push_back(reco::deltaPhi(*jetsv, matchedJet));
      jet_sv_chi2.push_back(jetsv->vertexChi2());
      const auto &dxy = btagbtvdeep::vertexDxy(*jetsv, *pv);
      jet_sv_dxy.push_back(dxy.value());
      jet_sv_dxysig.push_back(dxy.significance());
      const auto &d3d = btagbtvdeep::vertexD3d(*jetsv, *pv);
      jet_sv_d3d.push_back(d3d.value());
      jet_sv_d3dsig.push_back(d3d.significance());
      jet_sv_ntrack.push_back(jetsv->numberOfDaughters());
    }

    jet_pt = matchedJet.pt();
    jet_eta = matchedJet.eta();
    jet_phi = matchedJet.phi();
    jet_mass = matchedJet.mass();
    jet_ncand = matchedJet.neutralMultiplicity() + matchedJet.chargedMultiplicity();
    jet_nel = matchedJet.electronMultiplicity();
    jet_nmu = matchedJet.muonMultiplicity();
    jet_nbhad = matchedJet.jetFlavourInfo().getbHadrons().size();
    jet_nchad = matchedJet.jetFlavourInfo().getcHadrons().size();
    jet_hflav = matchedJet.jetFlavourInfo().getHadronFlavour();
    jet_pflav = matchedJet.jetFlavourInfo().getPartonFlavour();
    jet_lepflav = matchedJet.jetFlavourInfo().getLeptons().size();
    jet_deepcsv_probb = matchedJet.bDiscriminator("pfDeepCSVJetTags:probb");
    jet_deepcsv_probc = matchedJet.bDiscriminator("pfDeepCSVJetTags:probc");
    jet_deepcsv_probudsg = matchedJet.bDiscriminator("pfDeepCSVJetTags:probudsg");
    jet_deepjet_probb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probb");
    jet_deepjet_probc = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probc");
    jet_deepjet_probuds = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probuds");
    jet_deepjet_probg = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probg");
    jet_deepjet_probbb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probbb");
    jet_deepjet_problepb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:problepb");

    // fill clustered particles information

    //for (const auto &dau : matchedJet.daughterPtrVector()) {
    //
    //}
    //std::cout<<"nC: "<< matchedJet.daughterPtrVector().size()<<"\t"<<matchedJet.numberOfSourceCandidatePtrs()<<std::endl;

    int nConstituents = matchedJet.numberOfSourceCandidatePtrs();
    for(int k = 0; k < nConstituents; k++){
      reco::CandidatePtr pfcand = matchedJet.sourceCandidatePtr(k);
      //jet_pfcand_pt.push_back(pfcand->pt());
      //jet_pfcand_eta.push_back(pfcand->eta());
      //jet_pfcand_phi.push_back(pfcand->phi());

      const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*pfcand)); // from https://github.com/hqucms/DNNTuples/blob/prod/UL/10_6_X_MiniAODv2/Ntupler/src/PFCompleteFiller.cc
      jet_pfcand_px.push_back(packed_cand->px());
      jet_pfcand_py.push_back(packed_cand->py());
      jet_pfcand_pz.push_back(packed_cand->pz());
      jet_pfcand_energy.push_back(packed_cand->energy());
      jet_pfcand_pt.push_back(packed_cand->pt());
      jet_pfcand_pt_log.push_back(btagbtvdeep::catch_infs(std::log(packed_cand->pt()), -99));
      jet_pfcand_e_log.push_back(btagbtvdeep::catch_infs(std::log(packed_cand->energy()), -99));
      jet_pfcand_phirel.push_back(reco::deltaPhi(*packed_cand, matchedJet));
      jet_pfcand_etarel.push_back(etasign * (packed_cand->eta() - matchedJet.eta()));
      jet_pfcand_abseta.push_back(std::abs(packed_cand->eta()));
      jet_pfcand_puppiw.push_back(packed_cand->puppiWeight()); //FIXME
      jet_pfcand_charge.push_back(packed_cand->charge());
      jet_pfcand_isEl.push_back(std::abs(packed_cand->pdgId()) == 11);
      jet_pfcand_isMu.push_back(std::abs(packed_cand->pdgId()) == 13);
      jet_pfcand_isChargedHad.push_back(std::abs(packed_cand->pdgId()) == 211);
      jet_pfcand_isGamma.push_back(std::abs(packed_cand->pdgId()) == 22);
      jet_pfcand_isNeutralHad.push_back(std::abs(packed_cand->pdgId()) == 130);

      // for neutral
      float hcal_fraction = 0.;
      if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
        hcal_fraction = packed_cand->hcalFraction();
      }
      else if (packed_cand->isIsolatedChargedHadron()) {
        hcal_fraction = packed_cand->rawHcalFraction();
      }
      jet_pfcand_hcalFrac.push_back(hcal_fraction);
      jet_pfcand_hcalFracCalib.push_back(packed_cand->hcalFraction());

      // for charged
      jet_pfcand_VTX_ass.push_back(packed_cand->pvAssociationQuality());
      jet_pfcand_fromPV.push_back(packed_cand->fromPV());
      jet_pfcand_lostInnerHits.push_back(packed_cand->lostInnerHits());
      jet_pfcand_trackHighPurity.push_back(packed_cand->trackHighPurity());

      // impact parameters
      jet_pfcand_dz.push_back(btagbtvdeep::catch_infs(packed_cand->dz()));
      jet_pfcand_dzsig.push_back(packed_cand->bestTrack() ? btagbtvdeep::catch_infs(packed_cand->dz() / packed_cand->dzError()) : 0);
      jet_pfcand_dxy.push_back(btagbtvdeep::catch_infs(packed_cand->dxy()));
      jet_pfcand_dxysig.push_back(packed_cand->bestTrack() ? btagbtvdeep::catch_infs(packed_cand->dxy() / packed_cand->dxyError()) : 0);

      if (packed_cand->bestTrack()) {
        const auto* trk = packed_cand->bestTrack();
        jet_pfcand_normchi2.push_back(btagbtvdeep::catch_infs(trk->normalizedChi2()));
        jet_pfcand_quality.push_back(trk->qualityMask());
      }
      else {
        jet_pfcand_normchi2.push_back(999);
        jet_pfcand_quality.push_back(0);
      }


      //const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*pfcand));
      //if (not packed_cand and not reco_cand) throw edm::Exception(edm::errors::InvalidReference) << "Cannot convert to either reco::PFCandidate or pat::PackedCandidate";
      //if(packed_cand){
      //}
      //else if(reco_cand){
      //}
    }
    if(debug_)    std::cout<<"filling jet info"<<std::endl;
    tree_->Fill();
    Clear();
  }

  /////// background jets ///////
 /*
  edm::View<pat::Jet>::const_iterator jet;
  for(jet = jets.begin(); jet != jets.end(); ++jet){
    if(!(dR_genBs < 0.8 && jet->pt() > 20.)) continue;
    float deltar = 0.4;
    for(auto genjet = genjets.begin(); genjet != genjets.end(); ++genjet){
      if(reco::deltaR(*jet, *genjet) < deltar){
        deltar = reco::deltaR(*jet, *genjet);
        jet_genmatch_pt = genjet->pt();
        jet_genmatch_eta = genjet->eta();
        jet_genmatch_phi = genjet->phi();
        jet_genmatch_mass = genjet->mass();
      }
    }
    jet_pt = jet->pt();
    jet_eta = jet->eta();
    jet_phi = jet->phi();
    jet_mass = jet->mass();
    jet_ncand = jet->neutralMultiplicity() + jet->chargedMultiplicity();
    jet_nel = jet->electronMultiplicity();
    jet_nmu = jet->muonMultiplicity();
    jet_nbhad = jet->jetFlavourInfo().getbHadrons().size();
    jet_nchad = jet->jetFlavourInfo().getcHadrons().size();
    jet_hflav = jet->jetFlavourInfo().getHadronFlavour(); 
    jet_pflav = jet->jetFlavourInfo().getPartonFlavour();
    jet_lepflav = jet->jetFlavourInfo().getLeptons().size();
    jet_deepcsv_probb = matchedJet.bDiscriminator("pfDeepCSVJetTags:probb");
    jet_deepcsv_probc = matchedJet.bDiscriminator("pfDeepCSVJetTags:probc");
    jet_deepcsv_probudsg = matchedJet.bDiscriminator("pfDeepCSVJetTags:probudsg");
    jet_deepjet_probb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probb");
    jet_deepjet_probc = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probc");
    jet_deepjet_probuds = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probuds");
    jet_deepjet_probg = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probg");
    jet_deepjet_probbb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:probbb");
    jet_deepjet_problepb = matchedJet.bDiscriminator("pfDeepFlavourJetTags:problepb");

    //jet_elflav = 0; jet_muflav = 0; jet_tauflav = 0;
    //const reco::GenParticleRefVector& leptons = jet->jetFlavourInfo().getLeptons();
    //std::cout << "                      # of clustered leptons: " << leptons.size() << std::endl;
    //for (reco::GenParticleRefVector::const_iterator it = leptons.begin(); it != leptons.end(); ++it) {
    //    std::cout<<"inside loop"<<std::endl;
    //    std::cout<<(*it)->pdgId()<<std::endl;
    //  if(abs((*it)->pdgId() == 11)) jet_elflav++;
    //  if(abs((*it)->pdgId() == 13)) jet_muflav++;
    //  if(abs((*it)->pdgId() == 15)) jet_tauflav++;
    //}
    //std::cout << "jet_elflav: "<<jet_elflav<<"\t"<<"jet_muflav: "<<jet_muflav<<"\t"<<"jet_tauflav: "<<jet_tauflav<<std::endl;

    int nConstituents = jet->numberOfSourceCandidatePtrs();
    for(int k = 0; k < nConstituents; k++){
      reco::CandidatePtr pfcand = jet->sourceCandidatePtr(k);
      jet_pfcand_pt.push_back(pfcand->pt());
      jet_pfcand_eta.push_back(pfcand->eta());
      jet_pfcand_phi.push_back(pfcand->phi());
    }
    
    if(debug_)    std::cout<<"filling jet info"<<std::endl;
    tree_->Fill();
    Clear();
  }//end of for loop
  */

  if(debug_)    std::cout<<"got jet info"<<std::endl;
}

void jetInfo::SetBranches(){
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
  AddBranch(&jet_deepcsv_probb, "jet_deepcsv_probb");
  AddBranch(&jet_deepcsv_probc, "jet_deepcsv_probc");
  AddBranch(&jet_deepcsv_probudsg, "jet_deepcsv_probudsg");
  AddBranch(&jet_deepjet_probb, "jet_deepjet_probb");
  AddBranch(&jet_deepjet_probc, "jet_deepjet_probc");
  AddBranch(&jet_deepjet_probuds, "jet_deepjet_probuds");
  AddBranch(&jet_deepjet_probg, "jet_deepjet_probg");
  AddBranch(&jet_deepjet_probbb, "jet_deepjet_probbb");
  AddBranch(&jet_deepjet_problepb, "jet_deepjet_problepb");
  AddBranch(&jet_genmatch_pt, "jet_genmatch_pt");
  AddBranch(&jet_genmatch_eta, "jet_genmatch_eta");
  AddBranch(&jet_genmatch_phi, "jet_genmatch_phi");
  AddBranch(&jet_genmatch_mass, "jet_genmatch_mass");
  AddBranch(&nsv, "nsv");
  AddBranch(&jet_sv_pt, "jet_sv_pt");
  AddBranch(&jet_sv_pt_log, "jet_sv_pt_log");
  AddBranch(&jet_sv_ptrel, "jet_sv_ptrel");
  AddBranch(&jet_sv_ptrel_log, "jet_sv_ptrel_log");
  AddBranch(&jet_sv_eta, "jet_sv_eta");
  AddBranch(&jet_sv_phi, "jet_sv_phi");
  AddBranch(&jet_sv_mass, "jet_sv_mass");
  AddBranch(&jet_sv_energy, "jet_sv_energy");
  AddBranch(&jet_sv_energy_log, "jet_sv_energy_log");
  AddBranch(&jet_sv_erel, "jet_sv_erel");
  AddBranch(&jet_sv_erel_log, "jet_sv_erel_log");
  AddBranch(&jet_sv_deta, "jet_sv_deta");
  AddBranch(&jet_sv_dphi, "jet_sv_dphi");
  AddBranch(&jet_sv_chi2, "jet_sv_chi2");
  AddBranch(&jet_sv_dxy, "jet_sv_dxy");
  AddBranch(&jet_sv_dxysig, "jet_sv_dxysig");
  AddBranch(&jet_sv_d3d, "jet_sv_d3d");
  AddBranch(&jet_sv_d3d, "jet_sv_d3d");
  AddBranch(&jet_sv_d3d, "jet_sv_d3d");
  AddBranch(&jet_sv_ntrack, "jet_sv_ntrack");
  AddBranch(&jet_sv_chi2, "jet_sv_chi2");
  AddBranch(&jet_sv_dxy, "jet_sv_dxy");

  AddBranch(&jet_pfcand_px , "jet_pfcand_px");
  AddBranch(&jet_pfcand_py , "jet_pfcand_py");
  AddBranch(&jet_pfcand_pz , "jet_pfcand_pz");
  AddBranch(&jet_pfcand_energy, "jet_pfcand_energy");
  AddBranch(&jet_pfcand_pt, "jet_pfcand_pt");
  AddBranch(&jet_pfcand_pt_log, "jet_pfcand_pt_log");
  AddBranch(&jet_pfcand_e_log, "jet_pfcand_e_log");
  AddBranch(&jet_pfcand_phirel, "jet_pfcand_phirel");
  AddBranch(&jet_pfcand_etarel, "jet_pfcand_etarel");
  AddBranch(&jet_pfcand_abseta, "jet_pfcand_abseta");
  AddBranch(&jet_pfcand_puppiw, "jet_pfcand_puppiw");
  AddBranch(&jet_pfcand_charge, "jet_pfcand_charge");
  AddBranch(&jet_pfcand_isEl, "jet_pfcand_isEl");
  AddBranch(&jet_pfcand_isMu, "jet_pfcand_isMu");
  AddBranch(&jet_pfcand_isChargedHad, "jet_pfcand_isChargedHad");
  AddBranch(&jet_pfcand_isGamma, "jet_pfcand_isGamma");
  AddBranch(&jet_pfcand_isNeutralHad, "jet_pfcand_isNeutralHad");
  AddBranch(&jet_pfcand_hcalFrac, "jet_pfcand_hcalFrac");
  AddBranch(&jet_pfcand_hcalFracCalib, "jet_pfcand_hcalFracCalib");
  AddBranch(&jet_pfcand_VTX_ass, "jet_pfcand_VTX_ass");
  AddBranch(&jet_pfcand_fromPV, "jet_pfcand_fromPV");
  AddBranch(&jet_pfcand_lostInnerHits, "jet_pfcand_lostInnerHits");
  AddBranch(&jet_pfcand_trackHighPurity, "jet_pfcand_trackHighPurity");
  AddBranch(&jet_pfcand_dz, "jet_pfcand_dz");
  AddBranch(&jet_pfcand_dzsig, "jet_pfcand_dzsig");
  AddBranch(&jet_pfcand_dxy, "jet_pfcand_dxy");
  AddBranch(&jet_pfcand_dxysig, "jet_pfcand_dxysig");
  AddBranch(&jet_pfcand_normchi2, "jet_pfcand_normchi2");
  AddBranch(&jet_pfcand_quality, "jet_pfcand_quality");

  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void jetInfo::Clear(){
  if(debug_)std::cout<<"clearing jet info"<<std::endl;
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
  jet_deepcsv_probb = 0;
  jet_deepcsv_probc = 0;
  jet_deepcsv_probudsg = 0;
  jet_deepjet_probb = 0;
  jet_deepjet_probc = 0;
  jet_deepjet_probuds = 0;
  jet_deepjet_probg = 0;
  jet_deepjet_probbb = 0;
  jet_deepjet_problepb = 0;
  jet_genmatch_pt = 0;
  jet_genmatch_eta = 0;
  jet_genmatch_phi = 0;
  jet_genmatch_mass = 0;

  nsv = 0;
  jet_sv_pt.clear();
  jet_sv_pt_log.clear();
  jet_sv_ptrel.clear();
  jet_sv_ptrel_log.clear();
  jet_sv_eta.clear();
  jet_sv_phi.clear();
  jet_sv_mass.clear();
  jet_sv_energy.clear();
  jet_sv_energy_log.clear();
  jet_sv_erel.clear();
  jet_sv_erel_log.clear();
  jet_sv_deta.clear();
  jet_sv_dphi.clear();
  jet_sv_chi2.clear();
  jet_sv_dxy.clear();
  jet_sv_dxysig.clear();
  jet_sv_d3d.clear();
  jet_sv_d3dsig.clear();
  jet_sv_ntrack.clear();

  jet_pfcand_px.clear();
  jet_pfcand_py.clear();
  jet_pfcand_pz.clear();
  jet_pfcand_energy.clear();
  jet_pfcand_pt.clear();
  jet_pfcand_pt_log.clear();
  jet_pfcand_e_log.clear();
  jet_pfcand_phirel.clear();
  jet_pfcand_etarel.clear();
  jet_pfcand_abseta.clear();
  jet_pfcand_puppiw.clear();
  jet_pfcand_charge.clear();
  jet_pfcand_isEl.clear();
  jet_pfcand_isMu.clear();
  jet_pfcand_isChargedHad.clear();
  jet_pfcand_isGamma.clear();
  jet_pfcand_isNeutralHad.clear();
  jet_pfcand_hcalFrac.clear();
  jet_pfcand_hcalFracCalib.clear();
  jet_pfcand_VTX_ass.clear();
  jet_pfcand_fromPV.clear();
  jet_pfcand_lostInnerHits.clear();
  jet_pfcand_trackHighPurity.clear();
  jet_pfcand_dz.clear();
  jet_pfcand_dzsig.clear();
  jet_pfcand_dxy.clear();
  jet_pfcand_dxysig.clear();
  jet_pfcand_normchi2.clear();
  jet_pfcand_quality.clear();

  if(debug_) std::cout<<"cleared"<<std::endl;
}
