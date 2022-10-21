#include "TreeMaker/TM/interface/eventInfo.h"

eventInfo::eventInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in eventInfo constructor"<<std::endl;
  //Event_n                = 0;
  if(debug) std::cout<<"in eventInfo constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

eventInfo::~eventInfo(){
  delete tree_;
}

void eventInfo::Fill(const edm::Event& iEvent){

  if(debug_)    std::cout<<"getting event info"<<std::endl;
  RunNumber   = iEvent.id().run();
  EventNumber = iEvent.id().event();
  LumiNumber  = iEvent.id().luminosityBlock();
  BXNumber = iEvent.bunchCrossing();
  //Event_n++;
  if(debug_)    std::cout<<"got event info"<<std::endl;
}

void eventInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&RunNumber,"RunNumber");
  AddBranch(&EventNumber,"EventNumber");
  AddBranch(&LumiNumber,"LumiNumber");
  AddBranch(&BXNumber,"BXNumber");
  //AddBranch(&Event_n,"Event_n");
  
  if(debug_)    std::cout<<"set branches"<<std::endl;
}
