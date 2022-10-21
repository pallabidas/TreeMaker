#ifndef __PILEUP_INFO_H_
#define __PILEUP_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
using namespace std;
using namespace edm;


class pileUpInfo : public baseTree{

 public:
  pileUpInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~pileUpInfo();
  void Fill(const edm::Event& iEvent);
  void SetBranches();


 private:
  pileUpInfo(){};
  edm::InputTag pileUpLabel_;
  //variables which would become branches
  int ootnpuVertices,npuVertices, npuVerticesm1, npuVerticesp1;
  float trueInteractions;
  
};

#endif

