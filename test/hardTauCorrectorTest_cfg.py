import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("test")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_2_1_9/RelValHiggsChargedTausM200/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_Tauola_v1/0002/00D4B0FC-608E-DD11-BCBA-000423D6A6F4.root'
    )
)

process.load("FWCore/MessageService/MessageLogger_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.hardTauCorrectorTest = cms.EDAnalyzer("HardTauCorrectorTest",
	## Optional input for the HardTauAlgorithm
        ## uncomment here any line which needs changing.
        # EtCaloOverTrackMin 	= cms.untracked.double(-0.9),
	# EtCaloOverTrackMax	= cms.double(0.0),
	# EtHcalOverTrackMin	= cms.double(-0.5),
	# EtHcalOverTrackMax	= cms.double(0.5),
	# SignalConeSize	= cms.double(0.2),
	# EcalConeSize		= cms.double(0.5),
	# MatchingConeSize	= cms.double(0.1),
	# Track_minPt		= cms.double(1.0),
	# tkmaxipt		= cms.double(0.03),
	# tkmaxChi2		= cms.double(100.),
	# tkminPixelHitsn	= cms.int32(2),
	# tkminTrackerHitsn	= cms.int32(8),
	# TrackCollection	= cms.InputTag("generalTracks"),
	# PVProducer		= cms.InputTag("offlinePrimaryVertices"),
	# EBRecHitCollection	= cms.InputTag("ecalRecHit:EcalRecHitsEB"),
	# EERecHitCollection	= cms.InputTag("ecalRecHit:EcalRecHitsEE"),
	# HBHERecHitCollection	= cms.InputTag("hbhereco"),
	# HORecHitCollection	= cms.InputTag("horeco"),
	# HFRecHitCollection	= cms.InputTag("hfreco"),

        ## for the test program: ProngSelection = "1prong","3prong","any" (any = 1 or 3 prong)
        ProngSelection = cms.string("1prong"),
        ## for the test program: TauJet jet energy correction parameters
        src            = cms.InputTag("iterativeCone5CaloJets"),
        tagName        = cms.string("IterativeCone0.4_EtScheme_TowerEt0.5_E0.8_Jets871_2x1033PU_tau"),
        TauTriggerType = cms.int32(1)
)

process.runEDAna = cms.Path(
	process.hardTauCorrectorTest
)
