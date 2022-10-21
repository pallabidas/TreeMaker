import FWCore.ParameterSet.Config as cms

process = cms.Process("NTuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)
process.updatedPatJetsUpdatedJEC.userData.userFloats.src = []

#Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles.extend( [
	#"root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/FB4E61D1-33FC-C248-B47C-EEB5334F6331.root"
        #"root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/362B457C-213F-EE49-92B3-28063639CBBA.root"
#        "file:/eos/user/p/pdas/lowmass_bbtautau/signal_ma12/miniaod/277A52E2-73EE-1247-AAA2-777FD6CB7509.root"
        "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/BEA0934D-C518-9242-8390-9FBF304CF978.root"
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

#closes files after code is done running on that file
process.options = cms.untracked.PSet(
    fileMode = cms.untracked.string('NOMERGE')
)

process.load("RecoMET.METFilters.ecalBadCalibFilter_cfi")
baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313,
     872438182,872438951,872439990,872439864,
     872439609,872437181,872437182,872437053,
     872436794,872436667,872436536,872421541,
     872421413,872421414,872421031,872423083,872421439])
process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal       = baddetEcallist,
    taggingMode      = cms.bool(True),
    debug            = cms.bool(False)
    )

#input to analyzer
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.demo = cms.EDAnalyzer('TM',
                              debug_               = cms.untracked.bool(True),
                              isMC_               = cms.untracked.bool(True),
                              fillgenparticleInfo_ = cms.untracked.bool(False),  #false for data
                              genParticleLabel_    = cms.untracked.InputTag("prunedGenParticles"),
                              genEventLabel_       = cms.untracked.InputTag("generator"),
                              generatorlheLabel_   = cms.untracked.InputTag("externalLHEProducer"),
                              genJetLabel_         = cms.untracked.InputTag("slimmedGenJets"),
                              filleventInfo_       = cms.untracked.bool(True), 
                              fillpileUpInfo_      = cms.untracked.bool(False),   #false for data
                              pileUpLabel_         = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              filltriggerInfo_     = cms.untracked.bool(False),  # FIXME
                              hltlabel_            = cms.untracked.string("HLT"),
                              HLTriggerResults_    = cms.untracked.InputTag("TriggerResults","","HLT"),
                              triggerEventTag_     = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
                              fillvertexInfo_      = cms.untracked.bool(False),
                              vertexLabel_         = cms.untracked.InputTag("offlineSlimmedPrimaryVertices","",""),
                              fillmuonInfo_        = cms.untracked.bool(False),
                              muonLabel_           = cms.untracked.InputTag("slimmedMuons"),
                              fillPhotInfo_        = cms.untracked.bool(False), # FIXME
                              photonLabel_         = cms.untracked.InputTag('slimmedPhotons'),
                              fillelectronInfo_    = cms.untracked.bool(False),
                              electronLabel_       = cms.untracked.InputTag("slimmedElectrons"),
                              fillrhoInfo_         = cms.untracked.bool(False),
                              rhoLabel_            = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              fillpfmetInfo_       = cms.untracked.bool(False),
                              pfmetLabel_          = cms.untracked.InputTag("slimmedMETs"),
                              puppimetLabel_       = cms.untracked.InputTag("slimmedMETsPuppi"),  
                              metfilterspatLabel_  = cms.untracked.InputTag("TriggerResults::PAT"),
                              metfiltersrecoLabel_ = cms.untracked.InputTag("TriggerResults::RECO"),
                              ecalBadCalibLabel_   = cms.untracked.InputTag("ecalBadCalibReducedMINIAODFilter"),
                              fillpfjetInfo_       = cms.untracked.bool(True),
                              #jetCorrectorLabel_   = cms.untracked.InputTag("AK4PFchs"),
                              pfjetLabel_          = cms.untracked.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                              filltauInfo_         = cms.untracked.bool(False),
                              tauLabel_            = cms.untracked.InputTag("slimmedTaus"),                              
                              stageL1Trigger       = cms.uint32(1),
                              )

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histo.root"),
      closeFileFast = cms.untracked.bool(True)
  )

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=True,
                           )

from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True )
runMetCorAndUncFromMiniAOD(process,
                           isData=True,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )
process.puppi.useExistingWeights = True

#All paths are here
process.p = cms.Path(
#    process.egammaPostRecoSeq*
    process.ecalBadCalibReducedMINIAODFilter*
    process.patJetCorrFactorsUpdatedJEC*
    process.updatedPatJetsUpdatedJEC*
    process.selectedUpdatedPatJetsUpdatedJEC*
    process.fullPatMetSequence*
    process.puppiMETSequence*
    process.fullPatMetSequencePuppi*
    process.demo
    )


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

# process number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )

process.schedule=cms.Schedule(process.p)

process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

#print process.dumpPython()
