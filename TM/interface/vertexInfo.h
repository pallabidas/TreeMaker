#ifndef __VERTEX_INFO_H_
#define __VERTEX_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"


class vertexInfo : public baseTree{

 public:
  vertexInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~vertexInfo();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
  vertexInfo(){};
  edm::InputTag vertexLabel_;

  //variables which would become branches
  int Vertex_n;
  TClonesArray *Vertex_position; 
  std::vector<float> Vertex_chi2, Vertex_d0;
  std::vector<int>   Vertex_tracksize, Vertex_ndof;
  std::vector<bool>  Vertex_isFake;
  
};

#endif

