import FWCore.ParameterSet.Config as cms

from DQM.TrackingMonitorSource.histoHelper4hltTracking_cfi import *

vertex2vertexValidator = cms.EDAnalyzer("Vertex2VertexValidator",
    numVertex  = cms.InputTag("hltPixelVertices"),
    denVertex  = cms.InputTag("pixelVertices"),
    topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco"),

    dzmin         = cms.double(1.),

    # HistoProducerAlgo. Defines the set of plots to be booked and filled
    histoPSet = histoHelper4hltTracking,

)
