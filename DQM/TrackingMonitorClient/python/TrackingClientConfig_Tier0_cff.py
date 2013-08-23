import FWCore.ParameterSet.Config as cms

#  TrackingOfflineDQM (for Tier0 Harvesting Step) ####
trackingOfflineAnalyser = cms.EDAnalyzer("TrackingOfflineDQM",
    GlobalStatusFilling      = cms.untracked.int32(2),
    UsedWithEDMtoMEConverter = cms.untracked.bool(True),
    TopFolderName            = cms.untracked.string("Tracking"),                                     
    MEsDirectories           = cms.vstring(
         "TrackParameters/GeneralProperties/GoodTracks",
         "TrackParameters/HitProperties/GoodTracks"
    ),
    atLumiMEsDirectories     = cms.vstring(
         "TrackParameters/LSanalysis"
    ),
    TrackRatePSet            = cms.PSet(
         Name     = cms.string("NumberOfGoodTracks_"),
         LowerCut = cms.double(1.0),
         UpperCut = cms.double(1000.0),
    ),
    TrackChi2PSet            = cms.PSet(
         Name     = cms.string("GoodTrackChi2oNDF_"),
         LowerCut = cms.double(0.0),
         UpperCut = cms.double(25.0),
    ),
    TrackHitPSet            = cms.PSet(
         Name     = cms.string("GoodTrackNumberOfRecHitsPerTrack_"),
         LowerCut = cms.double(5.0),
         UpperCut = cms.double(20.0),
    ),
    GoodTrackFractionPSet   = cms.PSet(
         Name     = cms.string("FractionOfGoodTracks_"),
         LowerCut = cms.double(0.85),
         UpperCut = cms.double(1.1),
    )
)

trackingQTester = cms.EDAnalyzer("QualityTester",
    qtList = cms.untracked.FileInPath('DQM/TrackingMonitorClient/data/tracking_qualitytest_config_tier0.xml'),
    prescaleFactor = cms.untracked.int32(1),                               
    getQualityTestsFromFile = cms.untracked.bool(True)
)


# Sequence
TrackingOfflineDQMClient = cms.Sequence(trackingQTester*trackingOfflineAnalyser)


