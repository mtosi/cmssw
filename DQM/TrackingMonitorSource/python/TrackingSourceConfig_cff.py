import FWCore.ParameterSet.Config as cms

# Tracking Monitor ####
# Clone for Sim Data
import DQM.TrackingMonitor.TrackingMonitor_cfi
TrackMonSim = DQM.TrackingMonitor.TrackingMonitor_cfi.TrackMon.clone()
TrackMonSim.FolderName = 'Tracking/TrackParameters'

# Clone for Real Data (Cosmic)
import DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi
TrackMonReal = DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi.TrackerCosmicTrackMon.clone()
TrackMonReal.TrackProducer = 'ctfWithMaterialTracksP5'
TrackMonReal.FolderName = 'Tracking/TrackParameters'
TrackMonReal.AlgoName = 'CKFTk'
TrackMonReal.doSeedParameterHistos = True

# Clone for Real Data (Collison)
import DQM.TrackingMonitor.TrackingMonitor_cfi
TrackMonColl = DQM.TrackingMonitor.TrackingMonitor_cfi.TrackMon.clone()
TrackMonColl.TrackProducer = 'generalTracks'
TrackMonColl.FolderName = 'Tracking/TrackParameters'
TrackMonColl.AlgoName = 'CKFTk'
TrackMonColl.doSeedParameterHistos = True

# Sequences
TrackingSourcesSimData           = cms.Sequence(TrackMonSim )
TrackingSourcesRealData          = cms.Sequence(TrackMonReal)
TrackingSourcesRealDataCollision = cms.Sequence(TrackMonColl)




