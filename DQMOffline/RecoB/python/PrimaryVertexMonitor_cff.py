import FWCore.ParameterSet.Config as cms


pvMonitor = cms.EDAnalyzer("PrimaryVertexMonitor",
   TopFolderName  = cms.string("OfflinePV"),
   AlignmentLabel = cms.string("Alignment"),                           
   vertexLabel    = cms.InputTag("offlinePrimaryVertices"),
   beamSpotLabel  = cms.InputTag("offlineBeamSpot")
)
