import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("test")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'rfio:/castor/cern.ch/user/s/slehti/testData/Ztautau_GEN_SIM_RECO_MC_31X_V2_preproduction_311_v1.root'
    )
)

process.load("FWCore/MessageService/MessageLogger_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR09_R_34X_V2::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
process.JPTZSPCorrectorICone5.ResponseMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_31X_resptowers.txt")
process.JPTZSPCorrectorICone5.EfficiencyMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackNonEff_one.txt")
process.JPTZSPCorrectorICone5.LeakageMap = cms.string("JetMETCorrections/Configuration/data/CMSSW_167_TrackLeakage_one.txt")
process.ZSPJetCorJetIcone5.src = cms.InputTag("iterativeCone5CaloJets")

from JetMETCorrections.TauJet.TCTauAlgoParameters_cfi import *
from JetMETCorrections.TauJet.TauJetAlgoParameters_cfi import *
process.tcTauCorrectorTest = cms.EDAnalyzer("TCTauCorrectorTest",
        tauJetAlgoParameters,                                    
        tcTauAlgoParameters,
        MCTauCollection = cms.InputTag("TauMCProducer:HadronicTauOneAndThreeProng"),                                    
        ProngSelection  = cms.string("1prong"),
        TauJetEt        = cms.double(20.),
        TauJetEta       = cms.double(2.1)                               
)

#TauMCProducer
process.load("HLTriggerOffline.Tau.Validation.HLTTauReferences_cfi")

process.runEDAna = cms.Path(
        process.TauMCProducer *
        process.ZSPJetCorrections *
        process.ZSPrecoJetAssociations *
        process.JetPlusTrackCorrections *
	process.tcTauCorrectorTest
)
