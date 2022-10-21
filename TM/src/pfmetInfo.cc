#include "TreeMaker/TM/interface/pfmetInfo.h"

pfmetInfo::pfmetInfo(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in pfmet constructor"<<std::endl;
  pfmetLabel_     = iConfig.getUntrackedParameter<edm::InputTag> ("pfmetLabel_");
  puppimetLabel_  = iConfig.getUntrackedParameter<edm::InputTag> ("puppimetLabel_");
  ecalBadCalibLabel_ = iConfig.getUntrackedParameter<edm::InputTag> ("ecalBadCalibLabel_"); 
  metfilterspatLabel_  = iConfig.getUntrackedParameter<edm::InputTag> ("metfilterspatLabel_");
  metfiltersrecoLabel_ = iConfig.getUntrackedParameter<edm::InputTag> ("metfiltersrecoLabel_");
  isMC_         = iConfig.getUntrackedParameter<bool> ("isMC_");
  if(debug) std::cout<<"in pfmet constructor: calling SetBrances()"<<std::endl;
  SetBranches();
}

pfmetInfo::~pfmetInfo(){
  delete tree_;
}

void pfmetInfo::Fill(const edm::Event& iEvent){
  Clear();
  if(debug_)    std::cout<<"getting pfmet info"<<std::endl;

  //To directly access the filter decisions stored in miniaod, not running the filter modules:
  edm::Handle<edm::TriggerResults> METFilterResults;
  iEvent.getByLabel(metfilterspatLabel_, METFilterResults);
  if(!(METFilterResults.isValid())) iEvent.getByLabel(metfiltersrecoLabel_, METFilterResults);

  const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);

  unsigned int goodVerticesIndex = metfilterName.triggerIndex("Flag_goodVertices");
  filter_goodVertices.push_back(METFilterResults.product()->accept(goodVerticesIndex));

  unsigned int globalTightHalo2016FilterIndex = metfilterName.triggerIndex("Flag_globalTightHalo2016Filter");
  filter_globaltighthalo2016.push_back(METFilterResults.product()->accept(globalTightHalo2016FilterIndex));

  unsigned int globalSuperTightHalo2016FilterIndex = metfilterName.triggerIndex("Flag_globalSuperTightHalo2016Filter");
  filter_globalsupertighthalo2016.push_back(METFilterResults.product()->accept(globalSuperTightHalo2016FilterIndex));

  unsigned int HBHENoiseFilterIndex = metfilterName.triggerIndex("Flag_HBHENoiseFilter");
  filter_hbher2t.push_back(METFilterResults.product()->accept(HBHENoiseFilterIndex));

  unsigned int HBHENoiseIsoFilterIndex = metfilterName.triggerIndex("Flag_HBHENoiseIsoFilter");
  filter_hbheiso.push_back(METFilterResults.product()->accept(HBHENoiseIsoFilterIndex));

  unsigned int EcalDeadCellTriggerPrimitiveFilterIndex = metfilterName.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filter_ecaltp.push_back(METFilterResults.product()->accept(EcalDeadCellTriggerPrimitiveFilterIndex));

  unsigned int BadPFMuonFilterIndex = metfilterName.triggerIndex("Flag_BadPFMuonFilter");
  filter_badPFMuon.push_back(METFilterResults.product()->accept(BadPFMuonFilterIndex));

  unsigned int eeBadScFilterIndex = metfilterName.triggerIndex("Flag_eeBadScFilter");
  filter_ecalsc.push_back(METFilterResults.product()->accept(eeBadScFilterIndex));

  edm::Handle<bool> ifilterecalBadCalib;
  iEvent.getByLabel(ecalBadCalibLabel_, ifilterecalBadCalib);
  filter_ecalBadCalib.push_back(*ifilterecalBadCalib);

  edm::Handle<edm::View<pat::MET> > metPFHandle;
  iEvent.getByLabel(pfmetLabel_,metPFHandle);
  const edm::View<pat::MET> & metsPF = *metPFHandle;

  edm::Handle<edm::View<pat::MET> > metPuppiHandle;
  iEvent.getByLabel(puppimetLabel_,metPuppiHandle);
  const edm::View<pat::MET> & metsPuppi = *metPuppiHandle;

  if(not iEvent.getByLabel(pfmetLabel_,metPFHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<pfmetLabel_<<std::endl; 
    exit(0);
  }  
  for(int i = 0; i<1; i++)
    {
      PFMetPt.push_back(-99);
      PFMetPx.push_back(-99);
      PFMetPy.push_back(-99);
      PFMetPhi.push_back(-99);
      if(isMC_){
        PFMetPt_JetResUp.push_back(-99);
        PFMetPt_JetResDown.push_back(-99);
        PFMetPt_JetEnUp.push_back(-99);
        PFMetPt_JetEnDown.push_back(-99);
        PFMetPt_UnclusteredEnUp.push_back(-99);
        PFMetPt_UnclusteredEnDown.push_back(-99);
        PFMetPx_JetResUp.push_back(-99);
        PFMetPx_JetResDown.push_back(-99);
        PFMetPx_JetEnUp.push_back(-99);
        PFMetPx_JetEnDown.push_back(-99);
        PFMetPx_UnclusteredEnUp.push_back(-99);
        PFMetPx_UnclusteredEnDown.push_back(-99);
        PFMetPy_JetResUp.push_back(-99);
        PFMetPy_JetResDown.push_back(-99);
        PFMetPy_JetEnUp.push_back(-99);
        PFMetPy_JetEnDown.push_back(-99);
        PFMetPy_UnclusteredEnUp.push_back(-99);
        PFMetPy_UnclusteredEnDown.push_back(-99);
        PFMetPhi_JetResUp.push_back(-99);
        PFMetPhi_JetResDown.push_back(-99);
        PFMetPhi_JetEnUp.push_back(-99);
        PFMetPhi_JetEnDown.push_back(-99);
        PFMetPhi_UnclusteredEnUp.push_back(-99);
        PFMetPhi_UnclusteredEnDown.push_back(-99);
      }
      RawMetPt.push_back(-99);
      RawMetPx.push_back(-99);
      RawMetPy.push_back(-99);
      RawMetPhi.push_back(-99);
      RawChsMetPt.push_back(-99);
      RawChsMetPx.push_back(-99);
      RawChsMetPy.push_back(-99);
      RawChsMetPhi.push_back(-99);
      RawTrkMetPt.push_back(-99);
      RawTrkMetPx.push_back(-99);
      RawTrkMetPy.push_back(-99);
      RawTrkMetPhi.push_back(-99);
    }
  
  if ( metPFHandle.isValid() ){

     PFMetPt[0]                          = (float) metsPF[0].pt();
     PFMetPx[0]                          = (float) metsPF[0].px();
     PFMetPy[0]                          = (float) metsPF[0].py();
     //PFMetPhi[0]                         = (float) correct_phi(metsPF[0].phi());
     PFMetPhi[0]                         = (float) metsPF[0].phi();
     if(isMC_){
       PFMetPt_JetResUp[0]                 = (float) metsPF[0].shiftedPt(pat::MET::JetResUp); 
       PFMetPt_JetResDown[0]               = (float) metsPF[0].shiftedPt(pat::MET::JetResDown);
       PFMetPt_JetEnUp[0]                  = (float) metsPF[0].shiftedPt(pat::MET::JetEnUp);
       PFMetPt_JetEnDown[0]                = (float) metsPF[0].shiftedPt(pat::MET::JetEnDown);
       PFMetPt_UnclusteredEnUp[0]          = (float) metsPF[0].shiftedPt(pat::MET::UnclusteredEnUp);
       PFMetPt_UnclusteredEnDown[0]        = (float) metsPF[0].shiftedPt(pat::MET::UnclusteredEnDown);
       PFMetPx_JetResUp[0]                 = (float) metsPF[0].shiftedPx(pat::MET::JetResUp);
       PFMetPx_JetResDown[0]               = (float) metsPF[0].shiftedPx(pat::MET::JetResDown);
       PFMetPx_JetEnUp[0]                  = (float) metsPF[0].shiftedPx(pat::MET::JetEnUp);
       PFMetPx_JetEnDown[0]                = (float) metsPF[0].shiftedPx(pat::MET::JetEnDown);
       PFMetPx_UnclusteredEnUp[0]          = (float) metsPF[0].shiftedPx(pat::MET::UnclusteredEnUp);
       PFMetPx_UnclusteredEnDown[0]        = (float) metsPF[0].shiftedPx(pat::MET::UnclusteredEnDown);
       PFMetPy_JetResUp[0]                 = (float) metsPF[0].shiftedPy(pat::MET::JetResUp);
       PFMetPy_JetResDown[0]               = (float) metsPF[0].shiftedPy(pat::MET::JetResDown);
       PFMetPy_JetEnUp[0]                  = (float) metsPF[0].shiftedPy(pat::MET::JetEnUp);
       PFMetPy_JetEnDown[0]                = (float) metsPF[0].shiftedPy(pat::MET::JetEnDown);
       PFMetPy_UnclusteredEnUp[0]          = (float) metsPF[0].shiftedPy(pat::MET::UnclusteredEnUp);
       PFMetPy_UnclusteredEnDown[0]        = (float) metsPF[0].shiftedPy(pat::MET::UnclusteredEnDown);
       PFMetPhi_JetResUp[0]                 = (float) metsPF[0].shiftedPhi(pat::MET::JetResUp);
       PFMetPhi_JetResDown[0]               = (float) metsPF[0].shiftedPhi(pat::MET::JetResDown);
       PFMetPhi_JetEnUp[0]                  = (float) metsPF[0].shiftedPhi(pat::MET::JetEnUp);
       PFMetPhi_JetEnDown[0]                = (float) metsPF[0].shiftedPhi(pat::MET::JetEnDown);
       PFMetPhi_UnclusteredEnUp[0]          = (float) metsPF[0].shiftedPhi(pat::MET::UnclusteredEnUp);
       PFMetPhi_UnclusteredEnDown[0]        = (float) metsPF[0].shiftedPhi(pat::MET::UnclusteredEnDown);
     }
     RawMetPt[0]                         = (float) metsPF[0].corPt(pat::MET::Raw);
     RawMetPx[0]                         = (float) metsPF[0].corPx(pat::MET::Raw);
     RawMetPy[0]                         = (float) metsPF[0].corPy(pat::MET::Raw);
     RawMetPhi[0]                        = (float) metsPF[0].corPhi(pat::MET::Raw);
     RawChsMetPt[0]                      = (float) metsPF[0].corPt(pat::MET::RawChs);
     RawChsMetPx[0]                      = (float) metsPF[0].corPx(pat::MET::RawChs);
     RawChsMetPy[0]                      = (float) metsPF[0].corPy(pat::MET::RawChs);
     RawChsMetPhi[0]                      = (float) metsPF[0].corPhi(pat::MET::RawChs);
     RawTrkMetPt[0]                      = (float) metsPF[0].corPt(pat::MET::RawTrk);
     RawTrkMetPx[0]                      = (float) metsPF[0].corPx(pat::MET::RawTrk);
     RawTrkMetPy[0]                      = (float) metsPF[0].corPy(pat::MET::RawTrk);
     RawTrkMetPhi[0]                      = (float) metsPF[0].corPhi(pat::MET::RawTrk);                      
   }
    if(debug_)    std::cout<<"got pfmet info"<<std::endl;

  if(not iEvent.getByLabel(puppimetLabel_,metPuppiHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<puppimetLabel_<<std::endl;
    exit(0);
  }
  for(int i = 0; i<1; i++)
    {
      PuppiMetPt.push_back(-99);
      PuppiMetPx.push_back(-99);
      PuppiMetPy.push_back(-99);
      PuppiMetPhi.push_back(-99);
      if(isMC_){
        PuppiMetPt_JetResUp.push_back(-99);
        PuppiMetPt_JetResDown.push_back(-99);
        PuppiMetPt_JetEnUp.push_back(-99);
        PuppiMetPt_JetEnDown.push_back(-99);
        PuppiMetPt_UnclusteredEnUp.push_back(-99);
        PuppiMetPt_UnclusteredEnDown.push_back(-99);
        PuppiMetPx_JetResUp.push_back(-99);
        PuppiMetPx_JetResDown.push_back(-99);
        PuppiMetPx_JetEnUp.push_back(-99);
        PuppiMetPx_JetEnDown.push_back(-99);
        PuppiMetPx_UnclusteredEnUp.push_back(-99);
        PuppiMetPx_UnclusteredEnDown.push_back(-99);
        PuppiMetPy_JetResUp.push_back(-99);
        PuppiMetPy_JetResDown.push_back(-99);
        PuppiMetPy_JetEnUp.push_back(-99);
        PuppiMetPy_JetEnDown.push_back(-99);
        PuppiMetPy_UnclusteredEnUp.push_back(-99);
        PuppiMetPy_UnclusteredEnDown.push_back(-99);
        PuppiMetPhi_JetResUp.push_back(-99);
        PuppiMetPhi_JetResDown.push_back(-99);
        PuppiMetPhi_JetEnUp.push_back(-99);
        PuppiMetPhi_JetEnDown.push_back(-99);
        PuppiMetPhi_UnclusteredEnUp.push_back(-99);
        PuppiMetPhi_UnclusteredEnDown.push_back(-99);
      }
      PuppiMetcorPt.push_back(-99);
      PuppiMetcorPx.push_back(-99);
      PuppiMetcorPy.push_back(-99);
      PuppiMetcorPhi.push_back(-99);
    }

  if ( metPuppiHandle.isValid() ){

     PuppiMetPt[0]                          = (float) metsPuppi[0].pt();
     PuppiMetPx[0]                          = (float) metsPuppi[0].px();
     PuppiMetPy[0]                          = (float) metsPuppi[0].py();
     PuppiMetPhi[0]                         = (float) metsPuppi[0].phi();
     if(isMC_){
       PuppiMetPt_JetResUp[0]                 = (float) metsPuppi[0].shiftedPt(pat::MET::JetResUp);
       PuppiMetPt_JetResDown[0]               = (float) metsPuppi[0].shiftedPt(pat::MET::JetResDown);
       PuppiMetPt_JetEnUp[0]                  = (float) metsPuppi[0].shiftedPt(pat::MET::JetEnUp);
       PuppiMetPt_JetEnDown[0]                = (float) metsPuppi[0].shiftedPt(pat::MET::JetEnDown);
       PuppiMetPt_UnclusteredEnUp[0]          = (float) metsPuppi[0].shiftedPt(pat::MET::UnclusteredEnUp);
       PuppiMetPt_UnclusteredEnDown[0]        = (float) metsPuppi[0].shiftedPt(pat::MET::UnclusteredEnDown);
       PuppiMetPx_JetResUp[0]                 = (float) metsPuppi[0].shiftedPx(pat::MET::JetResUp);
       PuppiMetPx_JetResDown[0]               = (float) metsPuppi[0].shiftedPx(pat::MET::JetResDown);
       PuppiMetPx_JetEnUp[0]                  = (float) metsPuppi[0].shiftedPx(pat::MET::JetEnUp);
       PuppiMetPx_JetEnDown[0]                = (float) metsPuppi[0].shiftedPx(pat::MET::JetEnDown);
       PuppiMetPx_UnclusteredEnUp[0]          = (float) metsPuppi[0].shiftedPx(pat::MET::UnclusteredEnUp);
       PuppiMetPx_UnclusteredEnDown[0]        = (float) metsPuppi[0].shiftedPx(pat::MET::UnclusteredEnDown);
       PuppiMetPy_JetResUp[0]                 = (float) metsPuppi[0].shiftedPy(pat::MET::JetResUp);
       PuppiMetPy_JetResDown[0]               = (float) metsPuppi[0].shiftedPy(pat::MET::JetResDown);
       PuppiMetPy_JetEnUp[0]                  = (float) metsPuppi[0].shiftedPy(pat::MET::JetEnUp);
       PuppiMetPy_JetEnDown[0]                = (float) metsPuppi[0].shiftedPy(pat::MET::JetEnDown);
       PuppiMetPy_UnclusteredEnUp[0]          = (float) metsPuppi[0].shiftedPy(pat::MET::UnclusteredEnUp);
       PuppiMetPy_UnclusteredEnDown[0]        = (float) metsPuppi[0].shiftedPy(pat::MET::UnclusteredEnDown);
       PuppiMetPhi_JetResUp[0]                 = (float) metsPuppi[0].shiftedPhi(pat::MET::JetResUp);
       PuppiMetPhi_JetResDown[0]               = (float) metsPuppi[0].shiftedPhi(pat::MET::JetResDown);
       PuppiMetPhi_JetEnUp[0]                  = (float) metsPuppi[0].shiftedPhi(pat::MET::JetEnUp);
       PuppiMetPhi_JetEnDown[0]                = (float) metsPuppi[0].shiftedPhi(pat::MET::JetEnDown);
       PuppiMetPhi_UnclusteredEnUp[0]          = (float) metsPuppi[0].shiftedPhi(pat::MET::UnclusteredEnUp);
       PuppiMetPhi_UnclusteredEnDown[0]        = (float) metsPuppi[0].shiftedPhi(pat::MET::UnclusteredEnDown);
     }
     PuppiMetcorPt[0]                       = (float) metsPuppi[0].corPt(pat::MET::Raw);
     PuppiMetcorPx[0]                       = (float) metsPuppi[0].corPx(pat::MET::Raw);
     PuppiMetcorPy[0]                       = (float) metsPuppi[0].corPy(pat::MET::Raw);
     PuppiMetcorPhi[0]                      = (float) metsPuppi[0].corPhi(pat::MET::Raw);
   }

    if(debug_)    std::cout<<"got puppimet info"<<std::endl;

}

void pfmetInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&PFMetPt,"PFMetPt");
  AddBranch(&PFMetPx,"PFMetPx");
  AddBranch(&PFMetPy,"PFMetPy");
  AddBranch(&PFMetPhi,"PFMetPhi");
  if(isMC_){
    AddBranch(&PFMetPt_JetResUp,"PFMetPt_JetResUp");
    AddBranch(&PFMetPt_JetResDown,"PFMetPt_JetResDown");
    AddBranch(&PFMetPt_JetEnUp,"PFMetPt_JetEnUp");
    AddBranch(&PFMetPt_JetEnDown,"PFMetPt_JetEnDown");
    AddBranch(&PFMetPt_UnclusteredEnUp,"PFMetPt_UnclusteredEnUp");
    AddBranch(&PFMetPt_UnclusteredEnDown,"PFMetPt_UnclusteredEnDown");
    AddBranch(&PFMetPx_JetResUp,"PFMetPx_JetResUp");
    AddBranch(&PFMetPx_JetResDown,"PFMetPx_JetResDown");
    AddBranch(&PFMetPx_JetEnUp,"PFMetPx_JetEnUp");
    AddBranch(&PFMetPx_JetEnDown,"PFMetPx_JetEnDown");
    AddBranch(&PFMetPx_UnclusteredEnUp,"PFMetPx_UnclusteredEnUp");
    AddBranch(&PFMetPx_UnclusteredEnDown,"PFMetPx_UnclusteredEnDown");
    AddBranch(&PFMetPy_JetResUp,"PFMetPy_JetResUp");
    AddBranch(&PFMetPy_JetResDown,"PFMetPy_JetResDown");
    AddBranch(&PFMetPy_JetEnUp,"PFMetPy_JetEnUp");
    AddBranch(&PFMetPy_JetEnDown,"PFMetPy_JetEnDown");
    AddBranch(&PFMetPy_UnclusteredEnUp,"PFMetPy_UnclusteredEnUp");
    AddBranch(&PFMetPy_UnclusteredEnDown,"PFMetPy_UnclusteredEnDown");
    AddBranch(&PFMetPhi_JetResUp,"PFMetPhi_JetResUp");
    AddBranch(&PFMetPhi_JetResDown,"PFMetPhi_JetResDown");
    AddBranch(&PFMetPhi_JetEnUp,"PFMetPhi_JetEnUp");
    AddBranch(&PFMetPhi_JetEnDown,"PFMetPhi_JetEnDown");
    AddBranch(&PFMetPhi_UnclusteredEnUp,"PFMetPhi_UnclusteredEnUp");
    AddBranch(&PFMetPhi_UnclusteredEnDown,"PFMetPhi_UnclusteredEnDown");
  }
  AddBranch(&RawMetPt,"RawMetPt");
  AddBranch(&RawMetPx,"RawMetPx");
  AddBranch(&RawMetPy,"RawMetPy");
  AddBranch(&RawMetPhi,"RawMetPhi");
  AddBranch(&RawChsMetPt,"RawChsMetPt");
  AddBranch(&RawChsMetPx,"RawChsMetPx");
  AddBranch(&RawChsMetPy,"RawChsMetPy");
  AddBranch(&RawChsMetPhi,"RawChsMetPhi");
  AddBranch(&RawTrkMetPt,"RawTrkMetPt");
  AddBranch(&RawTrkMetPx,"RawTrkMetPx");
  AddBranch(&RawTrkMetPy,"RawTrkMetPy");
  AddBranch(&RawTrkMetPhi,"RawTrkMetPhi");
  AddBranch(&PuppiMetPt,"PuppiMetPt");
  AddBranch(&PuppiMetPx,"PuppiMetPx");
  AddBranch(&PuppiMetPy,"PuppiMetPy");
  AddBranch(&PuppiMetPhi,"PuppiMetPhi");
  if(isMC_){
    AddBranch(&PuppiMetPt_JetResUp,"PuppiMetPt_JetResUp");
    AddBranch(&PuppiMetPt_JetResDown,"PuppiMetPt_JetResDown");
    AddBranch(&PuppiMetPt_JetEnUp,"PuppiMetPt_JetEnUp");
    AddBranch(&PuppiMetPt_JetEnDown,"PuppiMetPt_JetEnDown");
    AddBranch(&PuppiMetPt_UnclusteredEnUp,"PuppiMetPt_UnclusteredEnUp");
    AddBranch(&PuppiMetPt_UnclusteredEnDown,"PuppiMetPt_UnclusteredEnDown");
    AddBranch(&PuppiMetPx_JetResUp,"PuppiMetPx_JetResUp");
    AddBranch(&PuppiMetPx_JetResDown,"PuppiMetPx_JetResDown");
    AddBranch(&PuppiMetPx_JetEnUp,"PuppiMetPx_JetEnUp");
    AddBranch(&PuppiMetPx_JetEnDown,"PuppiMetPx_JetEnDown");
    AddBranch(&PuppiMetPx_UnclusteredEnUp,"PuppiMetPx_UnclusteredEnUp");
    AddBranch(&PuppiMetPx_UnclusteredEnDown,"PuppiMetPx_UnclusteredEnDown");
    AddBranch(&PuppiMetPy_JetResUp,"PuppiMetPy_JetResUp");
    AddBranch(&PuppiMetPy_JetResDown,"PuppiMetPy_JetResDown");
    AddBranch(&PuppiMetPy_JetEnUp,"PuppiMetPy_JetEnUp");
    AddBranch(&PuppiMetPy_JetEnDown,"PuppiMetPy_JetEnDown");
    AddBranch(&PuppiMetPy_UnclusteredEnUp,"PuppiMetPy_UnclusteredEnUp");
    AddBranch(&PuppiMetPy_UnclusteredEnDown,"PuppiMetPy_UnclusteredEnDown");
    AddBranch(&PuppiMetPhi_JetResUp,"PuppiMetPhi_JetResUp");
    AddBranch(&PuppiMetPhi_JetResDown,"PuppiMetPhi_JetResDown");
    AddBranch(&PuppiMetPhi_JetEnUp,"PuppiMetPhi_JetEnUp");
    AddBranch(&PuppiMetPhi_JetEnDown,"PuppiMetPhi_JetEnDown");
    AddBranch(&PuppiMetPhi_UnclusteredEnUp,"PuppiMetPhi_UnclusteredEnUp");
    AddBranch(&PuppiMetPhi_UnclusteredEnDown,"PuppiMetPhi_UnclusteredEnDown");
  }
  AddBranch(&PuppiMetcorPt, "PuppiMetcorPt");
  AddBranch(&PuppiMetcorPx, "PuppiMetcorPx");
  AddBranch(&PuppiMetcorPy, "PuppiMetcorPy");
  AddBranch(&PuppiMetcorPhi, "PuppiMetcorPhi");
  AddBranch(&filter_goodVertices, "filter_goodVertices");
  AddBranch(&filter_globaltighthalo2016, "filter_globaltighthalo2016");
  AddBranch(&filter_globalsupertighthalo2016, "filter_globalsupertighthalo2016");
  AddBranch(&filter_hbher2t, "filter_hbher2t");
  AddBranch(&filter_hbheiso, "filter_hbheiso");
  AddBranch(&filter_ecaltp, "filter_ecaltp");
  AddBranch(&filter_ecalsc, "filter_ecalsc");
  AddBranch(&filter_badPFMuon, "filter_badPFMuon");
  AddBranch(&filter_ecalBadCalib, "filter_ecalBadCalib");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void pfmetInfo::Clear(){
  if(debug_)std::cout<<"clearing pfMet info"<<std::endl;
  PFMetPt.clear();
  PFMetPx.clear();
  PFMetPy.clear();
  PFMetPhi.clear();
  if(isMC_){
    PFMetPt_JetResUp.clear();
    PFMetPt_JetResDown.clear();
    PFMetPt_JetEnUp.clear();
    PFMetPt_JetEnDown.clear();
    PFMetPt_UnclusteredEnUp.clear();
    PFMetPt_UnclusteredEnDown.clear();
    PFMetPx_JetResUp.clear();
    PFMetPx_JetResDown.clear();
    PFMetPx_JetEnUp.clear();
    PFMetPx_JetEnDown.clear();
    PFMetPx_UnclusteredEnUp.clear();
    PFMetPx_UnclusteredEnDown.clear();
    PFMetPy_JetResUp.clear();
    PFMetPy_JetResDown.clear();
    PFMetPy_JetEnUp.clear();
    PFMetPy_JetEnDown.clear();
    PFMetPy_UnclusteredEnUp.clear();
    PFMetPy_UnclusteredEnDown.clear();
    PFMetPhi_JetResUp.clear();
    PFMetPhi_JetResDown.clear();
    PFMetPhi_JetEnUp.clear();
    PFMetPhi_JetEnDown.clear();
    PFMetPhi_UnclusteredEnUp.clear();
    PFMetPhi_UnclusteredEnDown.clear();
  }
  RawMetPt.clear();
  RawMetPx.clear();
  RawMetPy.clear();
  RawMetPhi.clear();
  RawChsMetPt.clear();
  RawChsMetPx.clear();
  RawChsMetPy.clear();
  RawChsMetPhi.clear();
  RawTrkMetPt.clear();
  RawTrkMetPx.clear();
  RawTrkMetPy.clear();
  RawTrkMetPhi.clear();
  PuppiMetPt.clear();
  PuppiMetPx.clear();
  PuppiMetPy.clear();
  PuppiMetPhi.clear();
  if(isMC_){
    PuppiMetPt_JetResUp.clear();
    PuppiMetPt_JetResDown.clear();
    PuppiMetPt_JetEnUp.clear();
    PuppiMetPt_JetEnDown.clear();
    PuppiMetPt_UnclusteredEnUp.clear();
    PuppiMetPt_UnclusteredEnDown.clear();
    PuppiMetPx_JetResUp.clear();
    PuppiMetPx_JetResDown.clear();
    PuppiMetPx_JetEnUp.clear();
    PuppiMetPx_JetEnDown.clear();
    PuppiMetPx_UnclusteredEnUp.clear();
    PuppiMetPx_UnclusteredEnDown.clear();
    PuppiMetPy_JetResUp.clear();
    PuppiMetPy_JetResDown.clear();
    PuppiMetPy_JetEnUp.clear();
    PuppiMetPy_JetEnDown.clear();
    PuppiMetPy_UnclusteredEnUp.clear();
    PuppiMetPy_UnclusteredEnDown.clear();
    PuppiMetPhi_JetResUp.clear();
    PuppiMetPhi_JetResDown.clear();
    PuppiMetPhi_JetEnUp.clear();
    PuppiMetPhi_JetEnDown.clear();
    PuppiMetPhi_UnclusteredEnUp.clear();
    PuppiMetPhi_UnclusteredEnDown.clear();
  }
  PuppiMetcorPt.clear();
  PuppiMetcorPx.clear();
  PuppiMetcorPy.clear();
  PuppiMetcorPhi.clear();
  filter_goodVertices.clear();
  filter_globaltighthalo2016.clear();
  filter_globalsupertighthalo2016.clear();
  filter_hbher2t.clear();
  filter_hbheiso.clear();
  filter_ecaltp.clear();
  filter_ecalsc.clear();
  filter_badPFMuon.clear();
  filter_ecalBadCalib.clear();
  if(debug_) std::cout<<"cleared"<<std::endl;
}

float pfmetInfo::correct_phi(float phi){
	return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}

