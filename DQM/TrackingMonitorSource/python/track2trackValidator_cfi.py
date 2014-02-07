import FWCore.ParameterSet.Config as cms

from DQM.TrackingMonitorSource.histoHelper4hltTracking_cfi import *

track2trackValidator = cms.EDAnalyzer("Track2TrackValidator",
    numTrack           = cms.InputTag("hltIter4Merged"),
    denTrack           = cms.InputTag("generalTracks"),
    numBeamSpot        = cms.InputTag("hltOnlineBeamSpot"),
    denBeamSpot        = cms.InputTag("offlineBeamSpot"),
    numPrimaryVertices = cms.InputTag("hltPixelVertices"),
    denPrimaryVertices = cms.InputTag("pixelVertices"),
    topDirName         = cms.string("HLT/Tracking/ValidationWRTreco"),

    dRmin         = cms.double(0.002),

    # HistoProducerAlgo. Defines the set of plots to be booked and filled
    histoPSet = histoHelper4hltTracking,

)
