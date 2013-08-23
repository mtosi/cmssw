import FWCore.ParameterSet.Config as cms

# Clone for TrackingMonitor for all PDs but MinBias and Jet ####
import DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi
TrackerCollisionTrackMonCommon = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
TrackerCollisionTrackMonCommon.FolderName    = 'Tracking/TrackParameters'
TrackerCollisionTrackMonCommon.andOr         = cms.bool( False )
TrackerCollisionTrackMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackerCollisionTrackMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackerCollisionTrackMonCommon.andOrDcs      = cms.bool( False )
TrackerCollisionTrackMonCommon.errorReplyDcs = cms.bool( True )

# Clone for TrackingMonitor for MinBias ###
import DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi
TrackerCollisionTrackMonMB = DQM.TrackingMonitor.TrackerCollisionTrackingMonitor_cfi.TrackerCollisionTrackMon.clone()
TrackerCollisionTrackMonMB.FolderName    = 'Tracking/TrackParameters'
TrackerCollisionTrackMonMB.andOr         = cms.bool( False )
TrackerCollisionTrackMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackerCollisionTrackMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackerCollisionTrackMonMB.andOrDcs      = cms.bool( False )
TrackerCollisionTrackMonMB.errorReplyDcs = cms.bool( True )
TrackerCollisionTrackMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
TrackerCollisionTrackMonMB.hltInputTag = cms.InputTag( "TriggerResults::HLT" )
TrackerCollisionTrackMonMB.hltPaths = cms.vstring("HLT_ZeroBias_*")
TrackerCollisionTrackMonMB.hltDBKey = cms.string("Tracker_MB")
TrackerCollisionTrackMonMB.errorReplyHlt  = cms.bool( False )
TrackerCollisionTrackMonMB.andOrHlt = cms.bool(True) 


from DQM.TrackingMonitor.TrackingMonitorSeedNumber_cff import *

# LogMessageMonitor ####
from DQM.TrackingMonitor.LogMessageMonitor_cff import *
### LocalReco
# Clone for all PDs but MinBias ####
LocalRecoLogMessageMonCommon = DQM.TrackingMonitor.LogMessageMonitor_cff.LocalRecoLogMessageMon.clone()
LocalRecoLogMessageMonCommon.andOr         = cms.bool( False )
LocalRecoLogMessageMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
LocalRecoLogMessageMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
LocalRecoLogMessageMonCommon.andOrDcs      = cms.bool( False )
LocalRecoLogMessageMonCommon.errorReplyDcs = cms.bool( True )

# Clone for MinBias ###
LocalRecoLogMessageMonMB = DQM.TrackingMonitor.LogMessageMonitor_cff.LocalRecoLogMessageMon.clone()
LocalRecoLogMessageMonMB.andOr         = cms.bool( False )
LocalRecoLogMessageMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
LocalRecoLogMessageMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
LocalRecoLogMessageMonMB.andOrDcs      = cms.bool( False )
LocalRecoLogMessageMonMB.errorReplyDcs = cms.bool( True )
LocalRecoLogMessageMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
LocalRecoLogMessageMonMB.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
LocalRecoLogMessageMonMB.hltPaths      = cms.vstring("HLT_ZeroBias_*")
LocalRecoLogMessageMonMB.hltDBKey      = cms.string("Tracker_MB")
LocalRecoLogMessageMonMB.errorReplyHlt = cms.bool( False )
LocalRecoLogMessageMonMB.andOrHlt      = cms.bool(True) 

### Clusterizer
# Clone for all PDs but MinBias ####
ClusterizerLogMessageMonCommon = DQM.TrackingMonitor.LogMessageMonitor_cff.ClusterizerLogMessageMon.clone()
ClusterizerLogMessageMonCommon.andOr         = cms.bool( False )
ClusterizerLogMessageMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
ClusterizerLogMessageMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
ClusterizerLogMessageMonCommon.andOrDcs      = cms.bool( False )
ClusterizerLogMessageMonCommon.errorReplyDcs = cms.bool( True )

# Clone for MinBias ###
ClusterizerLogMessageMonMB = DQM.TrackingMonitor.LogMessageMonitor_cff.ClusterizerLogMessageMon.clone()
ClusterizerLogMessageMonMB.andOr         = cms.bool( False )
ClusterizerLogMessageMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
ClusterizerLogMessageMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
ClusterizerLogMessageMonMB.andOrDcs      = cms.bool( False )
ClusterizerLogMessageMonMB.errorReplyDcs = cms.bool( True )
ClusterizerLogMessageMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
ClusterizerLogMessageMonMB.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
ClusterizerLogMessageMonMB.hltPaths      = cms.vstring("HLT_ZeroBias_*")
ClusterizerLogMessageMonMB.hltDBKey      = cms.string("Tracker_MB")
ClusterizerLogMessageMonMB.errorReplyHlt = cms.bool( False )
ClusterizerLogMessageMonMB.andOrHlt      = cms.bool(True) 

### Seeding
# Clone for all PDs but MinBias ####
SeedingLogMessageMonCommon = DQM.TrackingMonitor.LogMessageMonitor_cff.SeedingLogMessageMon.clone()
SeedingLogMessageMonCommon.andOr         = cms.bool( False )
SeedingLogMessageMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
SeedingLogMessageMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
SeedingLogMessageMonCommon.andOrDcs      = cms.bool( False )
SeedingLogMessageMonCommon.errorReplyDcs = cms.bool( True )

# Clone for MinBias ###
SeedingLogMessageMonMB = DQM.TrackingMonitor.LogMessageMonitor_cff.SeedingLogMessageMon.clone()
SeedingLogMessageMonMB.andOr         = cms.bool( False )
SeedingLogMessageMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
SeedingLogMessageMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
SeedingLogMessageMonMB.andOrDcs      = cms.bool( False )
SeedingLogMessageMonMB.errorReplyDcs = cms.bool( True )
SeedingLogMessageMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
SeedingLogMessageMonMB.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
SeedingLogMessageMonMB.hltPaths      = cms.vstring("HLT_ZeroBias_*")
SeedingLogMessageMonMB.hltDBKey      = cms.string("Tracker_MB")
SeedingLogMessageMonMB.errorReplyHlt = cms.bool( False )
SeedingLogMessageMonMB.andOrHlt      = cms.bool(True) 

### TrackCandidate
# Clone for all PDs but MinBias ####
TrackCandidateLogMessageMonCommon = DQM.TrackingMonitor.LogMessageMonitor_cff.TrackCandidateLogMessageMon.clone()
TrackCandidateLogMessageMonCommon.andOr         = cms.bool( False )
TrackCandidateLogMessageMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackCandidateLogMessageMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackCandidateLogMessageMonCommon.andOrDcs      = cms.bool( False )
TrackCandidateLogMessageMonCommon.errorReplyDcs = cms.bool( True )

# Clone for MinBias ###
TrackCandidateLogMessageMonMB = DQM.TrackingMonitor.LogMessageMonitor_cff.TrackCandidateLogMessageMon.clone()
TrackCandidateLogMessageMonMB.andOr         = cms.bool( False )
TrackCandidateLogMessageMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackCandidateLogMessageMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackCandidateLogMessageMonMB.andOrDcs      = cms.bool( False )
TrackCandidateLogMessageMonMB.errorReplyDcs = cms.bool( True )
TrackCandidateLogMessageMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
TrackCandidateLogMessageMonMB.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
TrackCandidateLogMessageMonMB.hltPaths      = cms.vstring("HLT_ZeroBias_*")
TrackCandidateLogMessageMonMB.hltDBKey      = cms.string("Tracker_MB")
TrackCandidateLogMessageMonMB.errorReplyHlt = cms.bool( False )
TrackCandidateLogMessageMonMB.andOrHlt      = cms.bool(True) 

### TrackFinder
# Clone for all PDs but MinBias ####
TrackFinderLogMessageMonCommon = DQM.TrackingMonitor.LogMessageMonitor_cff.TrackFinderLogMessageMon.clone()
TrackFinderLogMessageMonCommon.andOr         = cms.bool( False )
TrackFinderLogMessageMonCommon.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackFinderLogMessageMonCommon.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackFinderLogMessageMonCommon.andOrDcs      = cms.bool( False )
TrackFinderLogMessageMonCommon.errorReplyDcs = cms.bool( True )

# Clone for MinBias ###
TrackFinderLogMessageMonMB = DQM.TrackingMonitor.LogMessageMonitor_cff.TrackFinderLogMessageMon.clone()
TrackFinderLogMessageMonMB.andOr         = cms.bool( False )
TrackFinderLogMessageMonMB.dcsInputTag   = cms.InputTag( "scalersRawToDigi" )
TrackFinderLogMessageMonMB.dcsPartitions = cms.vint32 ( 24, 25, 26, 27, 28, 29)
TrackFinderLogMessageMonMB.andOrDcs      = cms.bool( False )
TrackFinderLogMessageMonMB.errorReplyDcs = cms.bool( True )
TrackFinderLogMessageMonMB.dbLabel       = cms.string("SiStripDQMTrigger")
TrackFinderLogMessageMonMB.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
TrackFinderLogMessageMonMB.hltPaths      = cms.vstring("HLT_ZeroBias_*")
TrackFinderLogMessageMonMB.hltDBKey      = cms.string("Tracker_MB")
TrackFinderLogMessageMonMB.errorReplyHlt = cms.bool( False )
TrackFinderLogMessageMonMB.andOrHlt      = cms.bool(True) 

# dEdx monitor ####
from DQM.TrackingMonitor.dEdxAnalyzer_cff import *
import DQM.TrackingMonitor.dEdxAnalyzer_cfi
# Clone for all PDs but MinBias ####
dEdxMonCommon = DQM.TrackingMonitor.dEdxAnalyzer_cfi.dEdxAnalyzer.clone()

# Clone for MinBias ####
dEdxMonMB = DQM.TrackingMonitor.dEdxAnalyzer_cfi.dEdxAnalyzer.clone()
dEdxMonMB.dEdxParameters.andOr         = cms.bool( False )
dEdxMonMB.dEdxParameters.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
dEdxMonMB.dEdxParameters.hltPaths      = cms.vstring("HLT_ZeroBias_*")
dEdxMonMB.dEdxParameters.hltDBKey      = cms.string("Tracker_MB")
dEdxMonMB.dEdxParameters.errorReplyHlt = cms.bool( False )
dEdxMonMB.dEdxParameters.andOrHlt      = cms.bool(True) 

# Clone for SingleMu ####
dEdxMonMU = DQM.TrackingMonitor.dEdxAnalyzer_cfi.dEdxAnalyzer.clone()
dEdxMonMU.dEdxParameters.andOr         = cms.bool( False )
dEdxMonMU.dEdxParameters.hltInputTag   = cms.InputTag( "TriggerResults::HLT" )
dEdxMonMU.dEdxParameters.hltPaths      = cms.vstring("HLT_SingleMu40_Eta2p1_*")
dEdxMonMU.dEdxParameters.errorReplyHlt = cms.bool( False )
dEdxMonMU.dEdxParameters.andOrHlt      = cms.bool(True) 



# DQM Services
dqmInfoTracking = cms.EDAnalyzer("DQMEventInfo",
    subSystemFolder = cms.untracked.string('Tracking')
)

## Services needed for TkHistoMap
#TkDetMap = cms.Service("TkDetMap")
#SiStripDetInfoFileReade = cms.Service("SiStripDetInfoFileReader")

# Event History Producer
from  DPGAnalysis.SiStripTools.eventwithhistoryproducerfroml1abc_cfi import *


# temporary patch in order to have BXlumi 
from RecoLuminosity.LumiProducer.lumiProducer_cff import *

# temporary test in order to temporary produce the "goodPrimaryVertexCollection"
# define with a new name if changes are necessary, otherwise simply include
# it from CommonTools/ParticleFlow/python/goodOfflinePrimaryVertices_cfi.py
# uncomment when necessary
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
trackingDQMgoodOfflinePrimaryVertices = goodOfflinePrimaryVertices.clone()
trackingDQMgoodOfflinePrimaryVertices.filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
trackingDQMgoodOfflinePrimaryVertices.src=cms.InputTag('offlinePrimaryVertices')
trackingDQMgoodOfflinePrimaryVertices.filter = cms.bool(False)

# Sequence
TrackingDQMSourceTier0 = cms.Sequence(
    # dEdx monitoring
    dedxDQMHarm2SP * dedxDQMHarm2SO * dedxDQMHarm2PO * dEdxMonCommon    
#    # temporary patch in order to have BXlumi
#    * lumiProducer
#    # temporary test in order to have the "goodPrimaryVertexCollection"
#    * trackingDQMgoodOfflinePrimaryVertices
   # tracks monitoring
    *TrackerCollisionTrackMonCommon
   # seeding monitoring [one per iteration] 
    *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
   # MessageLog
    * LocalRecoLogMessageMonCommon * ClusterizerLogMessageMonCommon * SeedingLogMessageMonCommon * TrackCandidateLogMessageMonCommon * TrackFinderLogMessageMonCommon
    *dqmInfoTracking)

TrackingDQMSourceTier0Common = cms.Sequence(
    # dEdx monitoring
    dedxDQMHarm2SP * dedxDQMHarm2SO * dedxDQMHarm2PO * dEdxMonCommon    
#    # temporary patch in order to have BXlumi
#    * lumiProducer
#    # temporary test in order to have the "goodPrimaryVertexCollection"
#    * trackingDQMgoodOfflinePrimaryVertices
    # tracks monitoring
     *TrackerCollisionTrackMonCommon
    # seeding monitoring [one per iteration]
     *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
    # MessageLog
     * LocalRecoLogMessageMonCommon * ClusterizerLogMessageMonCommon * SeedingLogMessageMonCommon * TrackCandidateLogMessageMonCommon * TrackFinderLogMessageMonCommon
     *dqmInfoTracking)

TrackingDQMSourceTier0MinBias = cms.Sequence(
    # dEdx monitoring
    dedxDQMHarm2SP * dedxDQMHarm2SO * dedxDQMHarm2PO * dEdxMonMB    
#    # temporary patch in order to have BXlumi
#    * lumiProducer
#    # temporary test in order to have the "goodPrimaryVertexCollection"
#    * trackingDQMgoodOfflinePrimaryVertices
    # tracks monitoring
     *TrackerCollisionTrackMonMB
    # seeding monitoring [one per iteration]     
     *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
    # MessageLog
     * LocalRecoLogMessageMonMB * ClusterizerLogMessageMonMB * SeedingLogMessageMonMB * TrackCandidateLogMessageMonMB * TrackFinderLogMessageMonMB
     *dqmInfoTracking)

