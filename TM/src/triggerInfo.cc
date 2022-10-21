#include "TreeMaker/TM/interface/triggerInfo.h"

triggerInfo::triggerInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in triggerInfo constructor"<<std::endl;
  HLTriggerResults_ = iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults_");
  ntriggers = 0;
  if(debug) std::cout<<"in trigger constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

triggerInfo::~triggerInfo(){
  delete tree_;
}


void triggerInfo::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<std::string>& all_triggers, HLTConfigProvider& hltConfig_, HLTPrescaleProvider& hltPrescale_, std::string& hltlabel_, const size_t& MaxN =200 ){

  if(debug_)    std::cout<<"getting HL Trigger info"<<std::endl;

  Clear();
  edm::Handle<TriggerResults> HLTR;
  iEvent.getByLabel(HLTriggerResults_,HLTR);
  if(debug_) std::cout<<"prescale set: "<<hltPrescale_.prescaleSet( iEvent, iSetup)<<std::endl;

   if (HLTR.isValid())
     {
       const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*HLTR);
       
       vector<int> idx;
      Int_t hsize = Int_t(HLTR->size());


      ntriggers = all_triggers.size();
      for(int i = 0; i< ntriggers;i++){

  	all_triggerprescales.push_back(0);
  	all_ifTriggerpassed.push_back(0);
  	idx.push_back(triggerNames_.triggerIndex(all_triggers[i]));
 	
  	if(idx[i] < hsize){
          all_triggernames.push_back(all_triggers[i]);
  	  all_ifTriggerpassed[i]=HLTR->accept(idx[i]);
          all_triggerprescales[i] = hltPrescale_.prescaleValue( iEvent, iSetup, all_triggers[i]);
 	}	
 	
  	if(debug_){
  	  std::cout<<"prescale for "<<all_triggers[i]<<" is: "<< all_triggerprescales[i]<<std::endl;
  	  std::cout<<"if triggger passed for "<<all_triggers[i]<<" : "<<all_ifTriggerpassed[i]<<std::endl;
  	}
	
      }//loop over trigger
      
    }//HLT is valid collection
  
  if(debug_)    std::cout<<"got trigger info"<<std::endl;
}

void triggerInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&all_triggernames,"all_triggernames");
  AddBranch(&all_triggerprescales,"all_triggerprescales");
  AddBranch(&all_ifTriggerpassed ,"all_ifTriggerpassed" );
    if(debug_)    std::cout<<"set branches"<<std::endl;
}

void triggerInfo::Clear(){
  if(debug_)std::cout<<"clearing trigger info"<<std::endl;
  all_triggernames.clear();
  all_triggerprescales.clear();
  all_ifTriggerpassed.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
