import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQM_cfg import *
DQMStore.collateHistograms =cms.untracked.bool(True)
from DQM.TrackingMonitorSource.track2trackValidator_cfi import *

trackSelector = cms.EDFilter('TrackSelector',
    src = cms.InputTag('generalTracks'),
    cut = cms.string("")
)
highPurityTracks = trackSelector.clone()
highPurityTracks.cut = cms.string("quality('highPurity')")

from PhysicsTools.RecoAlgos.recoTrackSelector_cfi import recoTrackSelector

iter0highPurityTracks = recoTrackSelector.clone()
iter0highPurityTracks.algorithm=cms.vstring("iter0")
iter0highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter0highPurityTracks.src = cms.InputTag("highPurityTracks")

iter1highPurityTracks = recoTrackSelector.clone()
iter1highPurityTracks.algorithm=cms.vstring("iter1")
iter1highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter1highPurityTracks.src = cms.InputTag("highPurityTracks")

iter2highPurityTracks = recoTrackSelector.clone()
iter2highPurityTracks.algorithm=cms.vstring("iter2")
iter2highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter2highPurityTracks.src = cms.InputTag("highPurityTracks")

iter3highPurityTracks = recoTrackSelector.clone()
iter3highPurityTracks.algorithm=cms.vstring("iter3")
iter3highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter3highPurityTracks.src = cms.InputTag("highPurityTracks")

iter4highPurityTracks = recoTrackSelector.clone()
iter4highPurityTracks.algorithm=cms.vstring("iter4")
iter4highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter4highPurityTracks.src = cms.InputTag("highPurityTracks")

iter5highPurityTracks = recoTrackSelector.clone()
iter5highPurityTracks.algorithm=cms.vstring("iter5")
iter5highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter5highPurityTracks.src = cms.InputTag("highPurityTracks")

iter6highPurityTracks = recoTrackSelector.clone()
iter6highPurityTracks.algorithm=cms.vstring("iter6")
iter6highPurityTracks.beamSpot = cms.InputTag("offlineBeamSpot")
iter6highPurityTracks.src = cms.InputTag("highPurityTracks")

hltIter4Merged2highPurityTracks = track2trackValidator.clone()
hltIter4Merged2highPurityTracks.numTrack      = cms.InputTag("hltIter4Merged")
hltIter4Merged2highPurityTracks.denTrack      = cms.InputTag("highPurityTracks")
hltIter4Merged2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter4Merged2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter4Merged2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter4Merged")

hltIter3Merged2highPurityTracks = track2trackValidator.clone()
hltIter3Merged2highPurityTracks.numTrack      = cms.InputTag("hltIter3Merged")
hltIter3Merged2highPurityTracks.denTrack      = cms.InputTag("highPurityTracks")
hltIter3Merged2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter3Merged2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter3Merged2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter3Merged")

hltIter2Merged2highPurityTracks = track2trackValidator.clone()
hltIter2Merged2highPurityTracks.numTrack      = cms.InputTag("hltIter2Merged")
hltIter2Merged2highPurityTracks.denTrack      = cms.InputTag("highPurityTracks")
hltIter2Merged2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter2Merged2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter2Merged2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter2Merged")

hltIter1Merged2highPurityTracks = track2trackValidator.clone()
hltIter1Merged2highPurityTracks.numTrack      = cms.InputTag("hltIter1Merged")
hltIter1Merged2highPurityTracks.denTrack      = cms.InputTag("highPurityTracks")
hltIter1Merged2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter1Merged2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter1Merged2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter1Merged")
hltIter1Merged2highPurityTracks.minEta = cms.double(-2.5)
    
hltIter0Merged2highPurityTracks = track2trackValidator.clone()
hltIter0Merged2highPurityTracks.numTrack      = cms.InputTag("hltPFlowTrackSelectionHighPurity")
hltIter0Merged2highPurityTracks.denTrack      = cms.InputTag("highPurityTracks")
hltIter0Merged2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter0Merged2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter0Merged2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter0Merged")

hltIter0Merged2iter0highPurityTracks = track2trackValidator.clone()
hltIter0Merged2iter0highPurityTracks.numTrack      = cms.InputTag("hltPFlowTrackSelectionHighPurity")
hltIter0Merged2iter0highPurityTracks.denTrack      = cms.InputTag("iter0highPurityTracks")
hltIter0Merged2iter0highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter0Merged2iter0highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter0Merged2iter0highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter0Merged/WRTiter0highPurityTracks")

hltIter1Merged2iter1highPurityTracks = track2trackValidator.clone()
hltIter1Merged2iter1highPurityTracks.numTrack      = cms.InputTag("hltIter1PFlowTrackSelectionHighPurity")
hltIter1Merged2iter1highPurityTracks.denTrack      = cms.InputTag("iter1highPurityTracks")
hltIter1Merged2iter1highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter1Merged2iter1highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter1Merged2iter1highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter1Merged/WRTiter1highPurityTracks")

hltIter2Merged2iter2highPurityTracks = track2trackValidator.clone()
hltIter2Merged2iter2highPurityTracks.numTrack      = cms.InputTag("hltIter2PFlowTrackSelectionHighPurity")
hltIter2Merged2iter2highPurityTracks.denTrack      = cms.InputTag("iter2highPurityTracks")
hltIter2Merged2iter2highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter2Merged2iter2highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter2Merged2iter2highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter2Merged/WRTiter2highPurityTracks")

hltIter3Merged2iter4highPurityTracks = track2trackValidator.clone()
hltIter3Merged2iter4highPurityTracks.numTrack      = cms.InputTag("hltIter3PFlowTrackSelectionHighPurity")
hltIter3Merged2iter4highPurityTracks.denTrack      = cms.InputTag("iter4highPurityTracks")
hltIter3Merged2iter4highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter3Merged2iter4highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter3Merged2iter4highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter3Merged/WRTiter4highPurityTracks")

hltIter4Merged2iter5highPurityTracks = track2trackValidator.clone()
hltIter4Merged2iter5highPurityTracks.numTrack      = cms.InputTag("hltIter4PFlowTrackSelectionHighPurity")
hltIter4Merged2iter5highPurityTracks.denTrack      = cms.InputTag("iter5highPurityTracks")
hltIter4Merged2iter5highPurityTracks.numBeamSpot   = cms.InputTag("hltOnlineBeamSpot")
hltIter4Merged2iter5highPurityTracks.denBeamSpot   = cms.InputTag("offlineBeamSpot")
hltIter4Merged2iter5highPurityTracks.topDirName    = cms.string("HLT/Tracking/ValidationWRTreco/hltIter4Merged/WRTiter4highPurityTracks")



track2trackValidatorSequence = cms.Sequence(
    highPurityTracks+
    hltIter0Merged2highPurityTracks+
    hltIter1Merged2highPurityTracks+
    hltIter2Merged2highPurityTracks+
    hltIter3Merged2highPurityTracks+
    hltIter4Merged2highPurityTracks+
    iter0highPurityTracks + hltIter0Merged2iter0highPurityTracks+
    iter1highPurityTracks + hltIter1Merged2iter1highPurityTracks+
    iter2highPurityTracks + hltIter2Merged2iter2highPurityTracks+
    iter4highPurityTracks + hltIter3Merged2iter4highPurityTracks+
    iter5highPurityTracks + hltIter4Merged2iter5highPurityTracks
)
