import FWCore.ParameterSet.Config as cms

import DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi
pixelTrackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
pixelTrackMonHLT.FolderName    = 'HLT/Tracking/pixelTracks'
pixelTrackMonHLT.TrackProducer    = 'hltPixelTracks'

iter0trackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
iter0trackMonHLT.FolderName    = 'HLT/Tracking/iter0tracks'
iter0trackMonHLT.TrackProducer = 'hltPFlowTrackSelectionHighPurity'

iter1trackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
iter1trackMonHLT.FolderName    = 'HLT/Tracking/iter1tracks'
iter1trackMonHLT.TrackProducer = 'hltIter1Merged'

iter2trackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
iter2trackMonHLT.FolderName    = 'HLT/Tracking/iter2tracks'
iter2trackMonHLT.TrackProducer = 'hltIter2Merged'

iter3trackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
iter3trackMonHLT.FolderName    = 'HLT/Tracking/iter3tracks'
iter3trackMonHLT.TrackProducer = 'hltIter3Merged'

iter4trackMonHLT = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
iter4trackMonHLT.FolderName    = 'HLT/Tracking/iter4tracks'
iter4trackMonHLT.TrackProducer = 'hltIter4Merged'

trackingHLTofflineMonSequence = cms.Sequence(
    pixelTrackMonHLT *
    iter0trackMonHLT *
    iter1trackMonHLT *
    iter2trackMonHLT *
    iter3trackMonHLT *
    iter4trackMonHLT 
)

