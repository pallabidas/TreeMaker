#include "TreeMaker/TM/interface/vertexInfo.h"

vertexInfo::vertexInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in vertex constructor"<<std::endl;
  vertexLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("vertexLabel_");
  Vertex_position  = new TClonesArray("TVector3");
  if(debug) std::cout<<"in vertex constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

vertexInfo::~vertexInfo(){
  delete tree_;
  delete Vertex_position;
}

void vertexInfo::Fill(const edm::Event& iEvent){
  Clear();
  if(debug_)    std::cout<<"getting vertex info"<<std::endl;

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(vertexLabel_, recVtxs);
  if(not iEvent.getByLabel(vertexLabel_,recVtxs)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<vertexLabel_<<std::endl; 
    exit(0);
  }  

  reco::VertexCollection::const_iterator vertex;
  for(vertex=recVtxs->begin();vertex!=recVtxs->end(); vertex++){
    TVector3 v3(vertex->x(),vertex->y(),vertex->z());
    new( (*Vertex_position)[Vertex_n]) TVector3(v3);
    if(debug_)   {
      std::cout<<"vertex:"<<Vertex_n<<std::endl;
      std::cout<<"vertex position:"<<vertex->x()<<"\t"<<vertex->y()<<"\t"<<vertex->z()<<std::endl;
    }
    Vertex_chi2.push_back((float)vertex->chi2());
    Vertex_d0.push_back((float)vertex->position().rho());
    Vertex_tracksize.push_back((int)vertex->tracksSize());
    Vertex_ndof.push_back((int)vertex->ndof());
    Vertex_isFake.push_back((bool)vertex->isFake());
    Vertex_n++;
  }    
  if(debug_)    std::cout<<"got vertices: "<<Vertex_n<<std::endl;
  if(debug_)    std::cout<<"got vertex info"<<std::endl;
}

void vertexInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Vertex_n,         "Vertex_n");
  AddBranch(&Vertex_chi2,      "Vertex_chi2");
  AddBranch(&Vertex_d0,        "Vertex_d0");
  AddBranch(&Vertex_tracksize, "Vertex_tracksize");
  AddBranch(&Vertex_ndof,      "Vertex_ndof");
  AddBranch(&Vertex_isFake,    "Vertex_isFake");
  AddBranch(&Vertex_position,  "Vertex_position");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void vertexInfo::Clear(){
  if(debug_)std::cout<<"clearing vertex info"<<std::endl;
  Vertex_n =0;
  Vertex_position->Clear();
  Vertex_chi2.clear();
  Vertex_d0.clear();
  Vertex_tracksize.clear();
  Vertex_ndof.clear();
  Vertex_isFake.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}
