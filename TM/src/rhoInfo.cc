#include "TreeMaker/TM/interface/rhoInfo.h"

rhoInfo::rhoInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in rho constructor"<<std::endl;
  rhoLabel_       = iConfig.getUntrackedParameter<edm::InputTag> ("rhoLabel_");
  if(debug) std::cout<<"in rho constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

rhoInfo::~rhoInfo(){
  delete tree_;
}

void rhoInfo::Fill(const edm::Event& iEvent){
  if(debug_)    std::cout<<"getting rho,sigma info"<<std::endl;
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoLabel_, rhoHandle);
  rho=0.;
  if(rhoHandle.isValid()) {
    rho= *(rhoHandle.product());
  }
  if(debug_)    std::cout<<"got rho info"<<std::endl;
}

void rhoInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&rho  ,"rho");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

