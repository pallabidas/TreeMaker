#ifndef __TRIGGER_INFO_H_
#define __TRIGGER_INFO_H_


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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
using namespace std;
using namespace edm;


class triggerInfo : public baseTree{

 public:
  triggerInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~triggerInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
	    std::vector<std::string>& all_triggers, HLTConfigProvider& hltConfig_, HLTPrescaleProvider& hltPrescale_,
	    std::string& hltlabel_, const size_t& MaxN );
  void SetBranches();
  void Clear();

 private:
  triggerInfo(){};
  edm::InputTag HLTriggerResults_;
  
  int ntriggers;
  std::vector<int>  all_triggerprescales;
  std::vector<bool> all_ifTriggerpassed;
  std::vector<std::string> all_triggernames; 
};

#endif

