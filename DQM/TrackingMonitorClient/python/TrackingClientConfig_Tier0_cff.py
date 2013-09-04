import FWCore.ParameterSet.Config as cms

#  TrackingOfflineDQM (for Tier0 Harvesting Step) ####
trackingOfflineAnalyser = cms.EDAnalyzer("TrackingOfflineDQM",
    GlobalStatusFilling        = cms.untracked.int32(2),
    UsedWithEDMtoMEConverter   = cms.untracked.bool(True),
    TopFolderName              = cms.untracked.string("Tracking"),                                     
    TrackingGlobalQualityPSets = cms.VPSet(
         cms.PSet(
             QT         = cms.string("Rate"),
             dir        = cms.string("TrackParameters/GeneralProperties/GoodTracks"),
             name       = cms.string("NumberOfGoodTracks_"),
             LSdir      = cms.string("TrackParameters/LSanalysis"),
             LSname     = cms.string("NumberOfGoodTracks_lumiFlag_"),
             LSlowerCut = cms.double(    1.0 ),
             LSupperCut = cms.double( 1000.0 )    
         ),
         cms.PSet(
             QT         = cms.string("Chi2"),
             dir        = cms.string("TrackParameters/GeneralProperties/GoodTracks"),
             name       = cms.string("GoodTrackChi2oNDF_"),
             LSdir      = cms.string("TrackParameters/LSanalysis"),
             LSname     = cms.string("GoodTrackChi2oNDF_lumiFlag_"),
             LSlowerCut = cms.double(  0.0 ),
             LSupperCut = cms.double( 25.0 )
         ),
         cms.PSet(
             QT         = cms.string("RecHits"),
             dir        = cms.string("TrackParameters/HitProperties/GoodTracks"),
             name       = cms.string("GoodTrackNumberOfRecHitsPerTrack_"),
             LSdir      = cms.string("TrackParameters/LSanalysis"),
             LSname     = cms.string("GoodTrackNumberOfRecHitsPerTrack_lumiFlag_"),
             LSlowerCut = cms.double(  5.0 ),
             LSupperCut = cms.double( 20.0 )
         ),
#         cms.PSet(
#             QT       = cms.string("goodFrac"),
#             dir      = cms.string("TrackParameters/GeneralProperties/GoodTracks"),
#             name     = cms.string("FractionOfGoodTracks_"),
#             lowerCut = cms.double( 0.85 ),
#             upperCut = cms.double( 1.1  )
#         )
    )
)

trackingQTester = cms.EDAnalyzer("QualityTester",
    qtList = cms.untracked.FileInPath('DQM/TrackingMonitorClient/data/tracking_qualitytest_config_tier0.xml'),
    prescaleFactor = cms.untracked.int32(1),                               
    getQualityTestsFromFile = cms.untracked.bool(True)
)


# Sequence
TrackingOfflineDQMClient = cms.Sequence(trackingQTester*trackingOfflineAnalyser)


