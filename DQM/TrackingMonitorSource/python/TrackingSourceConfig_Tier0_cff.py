import FWCore.ParameterSet.Config as cms

### load which are the tracks collection 2 be monitored
from DQM.TrackingMonitorSource.TrackCollections2monitor_cff import *

### load the different flavour of settings of the TrackingMonitor module
from DQM.TrackingMonitorSource.TrackerCollisionTrackingMonitor_cff import *


### define one EDAnalyzer per each track collection
### following suggestion 2. in
### https://hypernews.cern.ch/HyperNews/CMS/get/sw-develtools/1908/1.html
for tracks in selectedTracks :
    label = 'TrackerCollisionSelectedTrackMonCommon' + str(tracks)
    locals()[label] = TrackerCollisionTrackMonCommon.clone()
    locals()[label].TrackProducer = cms.InputTag(tracks)
    locals()[label].FolderName    = cms.string(mainfolderName[tracks])
    locals()[label].PVFolderName  = cms.string(vertexfolderName[tracks])
    locals()[label].setLabel(label)

    label = 'TrackerCollisionSelectedTrackMonMB' + str(tracks)                       
    locals()[label] = TrackerCollisionTrackMonMB.clone()
    locals()[label].TrackProducer = cms.InputTag(tracks)
    locals()[label].FolderName    = cms.string(mainfolderName[tracks])
    locals()[label].PVFolderName  = cms.string(vertexfolderName[tracks])
    locals()[label].setLabel(label)


### define the sequence w/ all the EDAnalyzers
TrackerCollisionSelectedTrackMonCommonSequence = cms.Sequence()
TrackerCollisionSelectedTrackMonMBSequence     = cms.Sequence()
for tracks in selectedTracks :
    label = 'TrackerCollisionSelectedTrackMonCommon' + str(tracks)
    TrackerCollisionSelectedTrackMonCommonSequence+=cms.Sequence(locals()[label])
    label = 'TrackerCollisionSelectedTrackMonMB' + str(tracks)
    TrackerCollisionSelectedTrackMonMBSequence+=cms.Sequence(locals()[label])
    
    
#-------------------------------------------------
# Tracking Monitor 
#-------------------------------------------------
import DQM.TrackingMonitor.TrackingMonitorSeed_cfi

from DQM.TrackingMonitorSource.IterTrackingModules4seedMonitoring_cfi import *
for step in selectedIterTrackingStep :
    label = 'TrackSeedMon'+str(step)
    locals()[label] = DQM.TrackingMonitor.TrackingMonitorSeed_cfi.TrackMonSeed.clone()
    locals()[label].TrackProducer = cms.InputTag("generalTracks")
    locals()[label].SeedProducer  = seedInputTag[step]
    locals()[label].TCProducer    = trackCandInputTag[step]
    locals()[label].AlgoName      = cms.string( str(step) )
    locals()[label].TkSeedSizeBin = trackSeedSizeBin[step]
    locals()[label].TkSeedSizeMin = trackSeedSizeMin[step]
    locals()[label].TkSeedSizeMax = trackSeedSizeMax[step]
    locals()[label].ClusterLabels = clusterLabel[step]
    if clusterLabel[step] == cms.vstring('Pix') :
        locals()[label].NClusPxBin = clusterBin[step]
        locals()[label].NClusPxMax = clusterMax[step]
    elif clusterLabel[step] == cms.vstring('Strip') or clusterLabel[step] == cms.vstring('Tot') :
        locals()[label].NClusStrBin = clusterBin[step]
        locals()[label].NClusStrMax = clusterMax[step]

seedingMonitorSequence = cms.Sequence()
for step in selectedIterTrackingStep :
    label = 'TrackSeedMon'+str(step)
    seedingMonitorSequence+=cms.Sequence(locals()[label])

# DQM Services
dqmInfoTracking = cms.EDAnalyzer("DQMEventInfo",
    subSystemFolder = cms.untracked.string('Tracking')
)

# LogMessageMonitor ####
### load which are the module to monitor
from DQM.TrackingMonitorSource.EDModules2monitor_cfi import *

### load the different flavour of settings of the LogMessageMonitor module
from DQM.TrackingMonitorSource.LogMessageMonitor_cff import *

for module in selectedModules :
    label = str(module)+'LogMessageMonCommon'
    locals()[label] = LogMessageMonCommon.clone()
    locals()[label].pluginsMonName = pluginsMonName[module]
    locals()[label].modules        = modulesLabel[module]
    locals()[label].categories     = categories[module]
    locals()[label].setLabel(label)

    label = str(module)+'LogMessageMonMB'
    locals()[label] = LogMessageMonMB.clone()
    locals()[label].pluginsMonName = pluginsMonName[module]
    locals()[label].modules        = modulesLabel[module]
    locals()[label].categories     = categories[module]
    locals()[label].setLabel(label)

logMessageMonCommonSequence = cms.Sequence()
logMessageMonMBSequence     = cms.Sequence()
for module in selectedModules :
    label = str(module)+'LogMessageMonCommon'
    logMessageMonCommonSequence+=cms.Sequence(locals()[label])
    label = str(module)+'LogMessageMonMB'
    logMessageMonMBSequence+=cms.Sequence(locals()[label])


# dEdx monitor ####
### load which dedx
from DQM.TrackingMonitorSource.dedxHarmonic2monitor_cfi import *

### load the different flavour of settings of the dEdxAnalyzer module
from DQM.TrackingMonitorSource.dEdxAnalyzer_cff import *


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
    dedxHarmonicSequence * dEdxMonCommon    
#    # temporary patch in order to have BXlumi
#    * lumiProducer
    # temporary test in order to have the "goodPrimaryVertexCollection"
    ## monitor track collections
    * selectedTracks2runSequence * TrackerCollisionSelectedTrackMonCommonSequence
    # seeding monitoring
    * seedingMonitorSequence
#    *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
     # MessageLog
    * logMessageMonCommonSequence
#    * LocalRecoLogMessageMonCommon * ClusterizerLogMessageMonCommon * SeedingLogMessageMonCommon * TrackCandidateLogMessageMonCommon * TrackFinderLogMessageMonCommon
    *dqmInfoTracking)


TrackingDQMSourceTier0Common = cms.Sequence(
    # dEdx monitoring
    dedxHarmonicSequence * dEdxMonCommon    
#    # temporary patch in order to have BXlumi
#    * lumiProducer
#    # temporary test in order to have the "goodPrimaryVertexCollection"
#    * trackingDQMgoodOfflinePrimaryVertices
      ## monitor track collections
    * selectedTracks2runSequence * TrackerCollisionSelectedTrackMonCommonSequence
      # seeding monitoring
    * seedingMonitorSequence
#      *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
    # MessageLog
      * logMessageMonCommonSequence
#    * LocalRecoLogMessageMonCommon * ClusterizerLogMessageMonCommon * SeedingLogMessageMonCommon * TrackCandidateLogMessageMonCommon * TrackFinderLogMessageMonCommon
    *dqmInfoTracking)

TrackingDQMSourceTier0MinBias = cms.Sequence(
    # dEdx monitoring
    dedxHarmonicSequence * dEdxMonCommon    
#    * lumiProducer
#    # temporary test in order to have the "goodPrimaryVertexCollection"
#    * trackingDQMgoodOfflinePrimaryVertices
    ## monitor track collections
    * selectedTracks2runSequence * TrackerCollisionSelectedTrackMonMBSequence
     # seeding monitoring
    * seedingMonitorSequence
#    *TrackMonStep0*TrackMonStep1*TrackMonStep2*TrackMonStep3*TrackMonStep4*TrackMonStep5*TrackMonStep6*TrackMonStep9*TrackMonStep10
    # MessageLog
    * logMessageMonMBSequence
#    * LocalRecoLogMessageMonMB * ClusterizerLogMessageMonMB * SeedingLogMessageMonMB * TrackCandidateLogMessageMonMB * TrackFinderLogMessageMonMB
    *dqmInfoTracking)

