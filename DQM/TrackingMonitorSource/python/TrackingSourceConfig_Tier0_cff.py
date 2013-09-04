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


############################################################################################################
############################################################################################################


# tracking monitor Sequence: VALIDATION
TrackingDQMSourceTier0 = cms.Sequence()
# dEdx monitoring
TrackingDQMSourceTier0+=(dedxHarmonicSequence * dEdxMonCommon)
## temporary patch in order to have BXlumi
#TrackingDQMSourceTier0+=(lumiProducer)
## temporary test in order to have the "goodPrimaryVertexCollection"
#TrackingDQMSourceTier0+=(trackingDQMgoodOfflinePrimaryVertices)
## monitor track collections
for tracks in selectedTracks :
    if tracks != 'generalTracks':
        TrackingDQMSourceTier0+=sequenceName[tracks]
    label = 'TrackerCollisionSelectedTrackMonCommon' + str(tracks)
    TrackingDQMSourceTier0+=(locals()[label])
# seeding monitoring
for step in selectedIterTrackingStep :
    label = 'TrackSeedMon'+str(step)
    TrackingDQMSourceTier0+=(locals()[label])
# MessageLog
for module in selectedModules :
    label = str(module)+'LogMessageMonCommon'
    TrackingDQMSourceTier0+=cms.Sequence(locals()[label])
TrackingDQMSourceTier0+=(dqmInfoTracking)
############################################################################################################


# tracking monitor Sequence: DQM all PD, but MiminumBias
TrackingDQMSourceTier0Common = cms.Sequence()
# dEdx monitoring
TrackingDQMSourceTier0Common+=(dedxHarmonicSequence * dEdxMonCommon)
## temporary patch in order to have BXlumi
#TrackingDQMSourceTier0Common+=(lumiProducer)
## temporary test in order to have the "goodPrimaryVertexCollection"
#TrackingDQMSourceTier0Common+=(trackingDQMgoodOfflinePrimaryVertices)

## monitor track collections
for tracks in selectedTracks :
    if tracks != 'generalTracks':
        TrackingDQMSourceTier0Common+=sequenceName[tracks]
    label = 'TrackerCollisionSelectedTrackMonCommon' + str(tracks)
    TrackingDQMSourceTier0Common+=(locals()[label])
# seeding monitoring
for step in selectedIterTrackingStep :
    label = 'TrackSeedMon'+str(step)
    TrackingDQMSourceTier0Common+=(locals()[label])
# MessageLog
for module in selectedModules :
    label = str(module)+'LogMessageMonCommon'
    TrackingDQMSourceTier0Common+=(locals()[label])
TrackingDQMSourceTier0Common+=(dqmInfoTracking)
############################################################################################################

# tracking monitor Sequence: DQM only MiminumBias PD
TrackingDQMSourceTier0MinBias = cms.Sequence()
# dEdx monitoring
TrackingDQMSourceTier0MinBias+=(dedxHarmonicSequence * dEdxMonCommon)    
## temporary patch in order to have BXlumi
#TrackingDQMSourceTier0MinBias+=(lumiProducer)
## temporary test in order to have the "goodPrimaryVertexCollection"
#TrackingDQMSourceTier0MinBias+=(trackingDQMgoodOfflinePrimaryVertices)
## monitor track collections
for tracks in selectedTracks :
    if tracks != 'generalTracks':
        TrackingDQMSourceTier0MinBias+=sequenceName[tracks]
    label = 'TrackerCollisionSelectedTrackMonMB' + str(tracks)
    TrackingDQMSourceTier0MinBias+=(locals()[label])
# seeding monitoring
for step in selectedIterTrackingStep :
    label = 'TrackSeedMon'+str(step)
    TrackingDQMSourceTier0MinBias+=(locals()[label])
# MessageLog
for module in selectedModules :
    label = str(module)+'LogMessageMonMB'
    TrackingDQMSourceTier0MinBias+=cms.Sequence(locals()[label])
TrackingDQMSourceTier0MinBias+=(dqmInfoTracking)
############################################################################################################

