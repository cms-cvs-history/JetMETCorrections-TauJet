import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("test")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:/tmp/slehti/A161_RECO.root'
    )
)

process.load("FWCore/MessageService/MessageLogger_cfi")

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'IDEAL_V1::All'

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.tcTauCorrectorTest = cms.EDAnalyzer("TCTauCorrectorTest",
	## Optional input for the TCTauAlgorithm
        ## uncomment here any line which needs changing.
	# EtCaloOverTrackMin 	= cms.untracked.double(-0.9),
	# EtCaloOverTrackMax	= cms.untracked.double(0.0),
	# EtHcalOverTrackMin	= cms.untracked.double(-0.3),
	# EtHcalOverTrackMax	= cms.untracked.double(1.0),
	# SignalConeSize	= cms.untracked.double(0.2),
	# EcalConeSize		= cms.untracked.double(0.5),
	# MatchingConeSize	= cms.untracked.double(0.1),
	# Track_minPt		= cms.untracked.double(1.0),
	# tkmaxipt		= cms.untracked.double(0.03),
	# tkmaxChi2		= cms.untracked.double(100.),
	# tkminPixelHitsn	= cms.untracked.int32(2),
	# tkminTrackerHitsn	= cms.untracked.int32(8),
	# TrackCollection	= cms.untracked.InputTag("generalTracks"),
	# PVProducer		= cms.untracked.InputTag("offlinePrimaryVertices"),
	# EBRecHitCollection	= cms.untracked.InputTag("ecalRecHit:EcalRecHitsEB"),
	# EERecHitCollection	= cms.untracked.InputTag("ecalRecHit:EcalRecHitsEE"),
	# HBHERecHitCollection	= cms.untracked.InputTag("hbhereco"),
	# HORecHitCollection	= cms.untracked.InputTag("horeco"),
	# HFRecHitCollection	= cms.untracked.InputTag("hfreco"),

        ## for the test program: ProngSelection = "1prong","3prong","any" (any = 1 or 3 prong)
        ProngSelection = cms.string("1prong"),
        ## for the test program: TauJet jet energy correction parameters
        src            = cms.InputTag("iterativeCone5CaloJets"),
        tagName        = cms.string("IterativeCone0.4_EtScheme_TowerEt0.5_E0.8_Jets871_2x1033PU_tau"),
        TauTriggerType = cms.int32(1)
)

process.runEDAna = cms.Path(
	process.tcTauCorrectorTest
)
