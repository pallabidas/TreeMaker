// system include files
#include <memory>
#include <string>

#include "TreeMaker/TM/interface/TM.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"

TM::TM(const edm::ParameterSet& iConfig):
  MaxN(200),
  hltPrescale_(iConfig, consumesCollector(), *this)
{
  debug_                       = iConfig.getUntrackedParameter<bool>("debug_");
  fillgenparticleInfo_         = iConfig.getUntrackedParameter<bool>("fillgenparticleInfo_");
  filleventInfo_               = iConfig.getUntrackedParameter<bool>("filleventInfo_");
  fillpileUpInfo_              = iConfig.getUntrackedParameter<bool>("fillpileUpInfo_");
  filltriggerInfo_             = iConfig.getUntrackedParameter<bool>("filltriggerInfo_");
  hltlabel_                    = iConfig.getUntrackedParameter<string>("hltlabel_");
  fillvertexInfo_              = iConfig.getUntrackedParameter<bool>("fillvertexInfo_");
  fillmuonInfo_                = iConfig.getUntrackedParameter<bool>("fillmuonInfo_");
  fillPhotInfo_                = iConfig.getUntrackedParameter<bool>("fillPhotInfo_");
  fillelectronInfo_            = iConfig.getUntrackedParameter<bool>("fillelectronInfo_");
  fillrhoInfo_                 = iConfig.getUntrackedParameter<bool>("fillrhoInfo_");
  fillpfmetInfo_               = iConfig.getUntrackedParameter<bool>("fillpfmetInfo_");
  fillpfjetInfo_               = iConfig.getUntrackedParameter<bool>("fillpfjetInfo_");
  filljetInfo_                 = iConfig.getUntrackedParameter<bool>("filljetInfo_");
  filltauInfo_                 = iConfig.getUntrackedParameter<bool>("filltauInfo_");
  

  vtxToken = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexLabel_")); 
  puToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getUntrackedParameter<edm::InputTag>("pileUpLabel_"));
  muonToken = consumes<edm::View<pat::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("muonLabel_"));
  eleToken = consumes<pat::ElectronCollection>(iConfig.getUntrackedParameter<edm::InputTag>("electronLabel_"));
  metToken = consumes<edm::View<pat::MET>>(iConfig.getUntrackedParameter<edm::InputTag>("pfmetLabel_"));
  puppimetToken = consumes<edm::View<pat::MET>>(iConfig.getUntrackedParameter<edm::InputTag>("puppimetLabel_"));
  pfjetToken = consumes<edm::View<pat::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("pfjetLabel_"));
  jetToken = consumes<edm::View<pat::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("jetLabel_"));
  //jecToken = consumes<JetCorrectorParametersCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jetCorrectorLabel_")); //getParameter<std::string>
  phoToken = consumes<edm::View<pat::Photon>>(iConfig.getUntrackedParameter<edm::InputTag>("photonLabel_"));
  TRToken = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults_"));
  tauToken = consumes<edm::View<pat::Tau>>(iConfig.getUntrackedParameter<edm::InputTag>("tauLabel_"));
  genToken = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticleLabel_"));
  genEvtToken = consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("genEventLabel_"));
  generatorlheToken = consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorlheLabel_"));
  genJetToken = consumes<edm::View<reco::GenJet>>(iConfig.getUntrackedParameter<edm::InputTag>("genJetLabel_"));
  svToken = consumes<edm::View<reco::VertexCompositePtrCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("svLabel_"));
  rhoToken = consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoLabel_"));
  ecalBadCalib_token = consumes<bool>(iConfig.getUntrackedParameter<edm::InputTag>("ecalBadCalibLabel_")); 
  metfilterspatToken_ = consumes<TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("metfilterspatLabel_"));
  metfiltersrecoToken_ = consumes<TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("metfiltersrecoLabel_"));

  //jetCorrector_ = iC.esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("jetCorrectorLabel_")));
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");

  if( filltriggerInfo_)       triggerInfo_     = new triggerInfo    ("reco", tree_, debug_, iConfig);
  if( fillgenparticleInfo_)   genparticleInfo_ = new genparticleInfo("reco", tree_, debug_, iConfig);
  if( filleventInfo_)         eventInfo_       = new eventInfo      ("reco", tree_, debug_, iConfig);
  if( fillpileUpInfo_)        pileUpInfo_      = new pileUpInfo     ("reco", tree_, debug_, iConfig);
  if( fillvertexInfo_)        vertexInfo_      = new vertexInfo     ("reco", tree_, debug_, iConfig);
  if( fillmuonInfo_)          muonInfo_        = new muonInfo       ("reco", tree_, debug_, iConfig);
  if( fillPhotInfo_)          photonInfo_      = new photonInfo     ("reco", tree_, debug_, iConfig); 
  if( fillelectronInfo_)      electronInfo_    = new electronInfo   ("reco", tree_, debug_, iConfig); 
  if( fillrhoInfo_)           rhoInfo_         = new rhoInfo        ("reco", tree_, debug_, iConfig);  
  if( fillpfmetInfo_)         pfmetInfo_       = new pfmetInfo      ("pat", tree_, debug_, iConfig);
  if( fillpfjetInfo_)         pfjetInfo_       = new pfjetInfo      ("pat", tree_, debug_, iConfig);
  if( filljetInfo_)           jetInfo_         = new jetInfo        ("pat", tree_, debug_, iConfig);
  if( filltauInfo_)           tauInfo_         = new tauInfo        ("pat", tree_, debug_, iConfig);
  if(debug_) std::cout<<"got all the objects and branches in the tree"<<std::endl;
}


TM::~TM()
{
}

void
TM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vtxToken, vtxHandle);

  reco::Vertex vtx;
  //best known primary vertex
  math::XYZPoint pv(0, 0, 0);
  bool isGoodVertex = false;
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
      isGoodVertex = !(v->isFake()) && (v->ndof()>=4) && (abs(v->z())<=24.) && (abs((v->position().Rho())<=2.));
      if(isGoodVertex) {
          pv.SetXYZ(v->x(), v->y(), v->z());
          vtx = *v;
          break;
      }
  }    

  edm::Handle<double> rhoLeptonHandle;
  iEvent.getByToken(rhoToken, rhoLeptonHandle);
  rhoLepton = 0.;
  if(rhoLeptonHandle.isValid()) {
    rhoLepton = *(rhoLeptonHandle.product());
  }

  if( filleventInfo_)       eventInfo_      ->Fill(iEvent);
  if( fillgenparticleInfo_) genparticleInfo_->Fill(iEvent);
  if( fillpileUpInfo_)      pileUpInfo_     ->Fill(iEvent);
  if( filltriggerInfo_)     triggerInfo_    ->Fill(iEvent, iSetup, all_triggers, hltConfig_, hltPrescale_, hltlabel_, MaxN);
  if( fillvertexInfo_)      vertexInfo_     ->Fill(iEvent);
  if( fillmuonInfo_)        muonInfo_       ->Fill(iEvent, pv, vtx, rhoLepton);
  if( fillPhotInfo_)        photonInfo_     ->Fill(iEvent, iSetup);
  if( fillelectronInfo_)    electronInfo_   ->Fill(iEvent, pv, vtx, rhoLepton);
  if( fillrhoInfo_)         rhoInfo_        ->Fill(iEvent);
  if( fillpfmetInfo_)       pfmetInfo_      ->Fill(iEvent);
  if( fillpfjetInfo_)       pfjetInfo_      ->Fill(iEvent, iSetup);
  if( filljetInfo_)         jetInfo_        ->Fill(iEvent, iSetup);
  if( filltauInfo_)         tauInfo_        ->Fill(iEvent, pv, vtx);

  //if(debug_) std::cout<<"Filling tree"<<std::endl;
  //tree_->Fill();
  //if(debug_) std::cout<<"Filled tree"<<std::endl;
}

void
TM::beginJob(){

}

void
TM::endJob() {
  //tree_->Print();
}

void TM::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){
  if(debug_){
    cout<<"\n----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"BEGIN NEW RUN: "<<iRun.run()<<endl;
  }
  
  cout<<"running"<<endl;

  bool changed(true);
  if(hltPrescale_.init(iRun,iSetup,hltlabel_,changed)){
    HLTConfigProvider const&  hltConfig_ = hltPrescale_.hltConfigProvider();
    if(debug_) std::cout << "Initalizing HLTConfigProvider"  << std::endl;
    if(changed){
      if(debug_) cout<< "HLT config has changed wrt the previous Run"  << std::endl;
      std::vector<std::string> photon_triggers_in_run;
      photon_triggers_in_run.clear();
      if(debug_){
        cout<<" Trigger Table : "<<hltConfig_.tableName()<<endl;
      }
      unsigned int ntriggers = hltConfig_.size();
      cout<<"  "<<ntriggers <<" Total HLT path in this run\t"<<endl; 

      // Loop over all available triggers
      for(unsigned int t=0;t<ntriggers;++t){
        std::string hltname(hltConfig_.triggerName(t));
        //string string_search1 ("HLT_Photon30_R9Id90_HE10_IsoM_v"); //not in the trigger menu since 2017
        string string_search1 ("HLT_Photon30_HoverELoose_v");
        string string_search2 ("HLT_Photon50_R9Id90_HE10_IsoM_v");
        string string_search3 ("HLT_Photon75_R9Id90_HE10_IsoM_v");
        string string_search4 ("HLT_Photon90_R9Id90_HE10_IsoM_v");
        string string_search5 ("HLT_Photon120_R9Id90_HE10_IsoM_v");
	string string_search6 ("HLT_Photon165_R9Id90_HE10_IsoM_v");
        
        //search the trigger name for string_search. 
        size_t found1 = hltname.find(string_search1);
        size_t found2 = hltname.find(string_search2);
        size_t found3 = hltname.find(string_search3);
	size_t found4 = hltname.find(string_search4);
	size_t found5 = hltname.find(string_search5);
	size_t found6 = hltname.find(string_search6); 

	if(found1!=string::npos || found2!=string::npos || found3!=string::npos || found4!=string::npos || found5!=string::npos || found6!=string::npos){
	  photon_triggers_in_run.push_back(hltname);
	}
        
      }//loop over ntriggers

      //This has to be clean for every run as Menu get changed frequently
      all_triggers.clear(); 
  
      for(int x = 0; x< (int)photon_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==photon_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found)
          all_triggers.push_back(photon_triggers_in_run[x]);
      }//loop photon triggers
    }//loop over all available triggers

  }
  else{
    std::cout << " HLT config extraction failure with name " << "HLTConfigProvider" << std::endl;
  }
    
 if(debug_){
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"the triggers in HLT list till now are:"<<endl;
    for(int i = 0; i< (int)all_triggers.size();i++) cout<<"\t"<<all_triggers[i]<<endl;
  }    
}

DEFINE_FWK_MODULE(TM);
