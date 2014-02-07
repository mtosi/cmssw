import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQM_cfg import *
DQMStore.collateHistograms =cms.untracked.bool(True)
from DQM.TrackingMonitorSource.vertex2vertexValidator_cfi import *

from RecoPixelVertexing.PixelVertexFinding.PixelVertexes_cfi import *
#pixelVertices.WtAverage = cms.bool(True),
#pixelVertices.ZOffset = cms.double(5.0),
#pixelVertices.beamSpot = cms.InputTag("offlineBeamSpot"),
#pixelVertices.Verbosity = cms.int32(0),
#pixelVertices.UseError = cms.bool(True),
#pixelVertices.TrackCollection = cms.InputTag("pixelTracks"),
#pixelVertices.ZSeparation = cms.double(0.05),
#pixelVertices.NTrkMin = cms.int32(2),
#pixelVertices.Method2 = cms.bool(True),
#pixelVertices.Finder = cms.string('DivisiveVertexFinder'),
#pixelVertices.PtMin = cms.double(1.0)

onlinePixelVertices = pixelVertices.clone()
onlinePixelVertices.beamSpot        = cms.InputTag("hltOnlineBeamSpot")
onlinePixelVertices.TrackCollection = cms.InputTag("hltPixelTracks")

hltPixelVerticesPt0p5 = onlinePixelVertices.clone()
hltPixelVerticesPt0p5.PtMin = cms.double(0.5)

hltPixelVerticesPt1p5 = onlinePixelVertices.clone()
hltPixelVerticesPt1p5.PtMin = cms.double(1.5)

hltPixelVerticesPt2p0 = onlinePixelVertices.clone()
hltPixelVerticesPt2p0.PtMin = cms.double(2.0)

hltPixelVerticesPt2p5 = onlinePixelVertices.clone()
hltPixelVerticesPt2p5.PtMin = cms.double(2.5)

hltPixelVerticesPt3p0 = onlinePixelVertices.clone()
hltPixelVerticesPt3p0.PtMin = cms.double(3.0)

hltPixelVerticesPt3p5 = onlinePixelVertices.clone()
hltPixelVerticesPt3p5.PtMin = cms.double(3.5)

hltPixelVerticesPt4p0 = onlinePixelVertices.clone()
hltPixelVerticesPt4p0.PtMin = cms.double(4.0)

hltPixelVerticesPt4p5 = onlinePixelVertices.clone()
hltPixelVerticesPt4p5.PtMin = cms.double(4.5)

hltPixelVertices2offlinePixelVertices = vertex2vertexValidator.clone()
hltPixelVertices2offlinePixelVertices.numVertex  = cms.InputTag("hltPixelVertices")
hltPixelVertices2offlinePixelVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVertices2offlinePixelVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVertices/WRTofflinePixelVertices")

hltPixelVertices2offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVertices2offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVertices")
hltPixelVertices2offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVertices2offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVertices/WRTofflinePrimaryVertices") 

hltPixelVerticesPt0p52offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt0p52offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt0p5")
hltPixelVerticesPt0p52offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt0p52offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt0p5/WRTofflinePrimaryVertices") 


hltPixelVerticesPt1p52offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt1p52offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt1p5")
hltPixelVerticesPt1p52offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt1p52offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt1p5/WRTofflinePrimaryVertices") 


hltPixelVerticesPt2p02offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt2p02offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt2p0")
hltPixelVerticesPt2p02offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt2p02offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt2p0/WRTofflinePrimaryVertices") 

hltPixelVerticesPt2p52offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt2p52offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt2p5")
hltPixelVerticesPt2p52offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt2p52offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt2p5/WRTofflinePrimaryVertices") 

hltPixelVerticesPt3p02offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt3p02offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt3p0")
hltPixelVerticesPt3p02offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt3p02offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt3p0/WRTofflinePrimaryVertices") 

hltPixelVerticesPt3p52offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt3p52offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt3p5")
hltPixelVerticesPt3p52offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt3p52offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt3p5/WRTofflinePrimaryVertices") 

hltPixelVerticesPt4p02offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt4p02offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt4p0")
hltPixelVerticesPt4p02offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt4p02offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt4p0/WRTofflinePrimaryVertices") 

hltPixelVerticesPt4p52offlinePrimaryVertices = vertex2vertexValidator.clone()
hltPixelVerticesPt4p52offlinePrimaryVertices.numVertex  = cms.InputTag("hltPixelVerticesPt4p5")
hltPixelVerticesPt4p52offlinePrimaryVertices.denVertex  = cms.InputTag("pixelVertices")
hltPixelVerticesPt4p52offlinePrimaryVertices.topDirName = cms.string("HLT/Tracking/vertex/ValidationWRTreco/hltPixelVerticesPt4p5/WRTofflinePrimaryVertices") 


vertex2vertexValidatorSequence = cms.Sequence(
    hltPixelVertices2offlinePixelVertices
    + hltPixelVertices2offlinePrimaryVertices
    + hltPixelVerticesPt0p5 + hltPixelVerticesPt0p52offlinePrimaryVertices 
    + hltPixelVerticesPt1p5 + hltPixelVerticesPt1p52offlinePrimaryVertices 
    + hltPixelVerticesPt2p0 + hltPixelVerticesPt2p02offlinePrimaryVertices
    + hltPixelVerticesPt2p5 + hltPixelVerticesPt2p52offlinePrimaryVertices
    + hltPixelVerticesPt3p0 + hltPixelVerticesPt3p02offlinePrimaryVertices
    + hltPixelVerticesPt3p5 + hltPixelVerticesPt3p52offlinePrimaryVertices
    + hltPixelVerticesPt4p0 + hltPixelVerticesPt4p02offlinePrimaryVertices
    + hltPixelVerticesPt4p5 + hltPixelVerticesPt4p52offlinePrimaryVertices

    
)    
