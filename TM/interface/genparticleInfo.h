#ifndef __GENPARTICLE_INFO_H_
#define __GENPARTICLE_INFO_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TreeMaker/TM/interface/utils.h"
#include "TreeMaker/TM/interface/baseTree.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

using namespace std;

class genparticleInfo : public baseTree{

 public:
  genparticleInfo(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~genparticleInfo();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
 private:
  genparticleInfo(){};
  edm::InputTag genParticleLabel_;
  edm::InputTag genEventLabel_;
  edm::InputTag generatorlheLabel_;

  int Gen_Photon_n;
  TClonesArray *Gen_Photon_4Momentum;
  std::vector<bool> Gen_Photon_isPrompt; 
};

#endif

