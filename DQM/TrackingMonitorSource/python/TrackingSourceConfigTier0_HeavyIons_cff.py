import FWCore.ParameterSet.Config as cms

# TrackingMonitor ####
import DQM.TrackingMonitor.TrackerHeavyIonTrackingMonitor_cfi
TrackMon_hi = DQM.TrackingMonitor.TrackerHeavyIonTrackingMonitor_cfi.TrackerHeavyIonTrackMon.clone()
TrackMon_hi.FolderName          = 'Tracking/TrackParameters'
TrackMon_hi.BSFolderName        = 'Tracking/TrackParameters/BeamSpotParameters'

SiStripDQMTier0_hi = cms.Sequence(TrackMon_hi)
