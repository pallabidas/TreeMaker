#ifndef __EVENT_INFO_H_
#define __EVENT_INFO_H_


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
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

using namespace std;
using namespace edm;


class eventInfo : public baseTree{

 public:
  eventInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~eventInfo();
  void Fill(const edm::Event& iEvent);
  void SetBranches();


 private:
  eventInfo(){};
  //variables which would become branches
  unsigned int RunNumber, LumiNumber, BXNumber;
  ULong64_t EventNumber;
  //int Event_n;
  
};

#endif

