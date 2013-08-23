import FWCore.ParameterSet.Config as cms

# Hardware Monitor ###
from DQM.SiStripMonitorHardware.siStripFEDMonitor_P5_cff import *

# Pedestal Monitor ###
from DQM.SiStripMonitorPedestals.SiStripMonitorPedestals_cfi import *
PedsMon.OutputMEsInRootFile = False
PedsMon.StripQualityLabel = ''
PedsMon.RunTypeFlag = 'CalculatedPlotsOnly'

# Digi Monitor #####
from DQM.SiStripMonitorDigi.SiStripMonitorDigi_cfi import *
SiStripMonitorDigi.SelectAllDetectors = True

# Cluster Monitor ####
from DQM.SiStripMonitorCluster.SiStripMonitorCluster_cfi import *
SiStripMonitorCluster.OutputMEsInRootFile = False
SiStripMonitorCluster.SelectAllDetectors = True
SiStripMonitorCluster.StripQualityLabel = ''

# On/Off Track Cluster Monitor ####
# Clone for Sim data
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrackSim = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrackSim.TrackProducer = 'TrackRefitter'
SiStripMonitorTrackSim.TrackLabel    = ''
SiStripMonitorTrackSim.Cluster_src   = 'siStripClusters'
SiStripMonitorTrackSim.Mod_On        = True

# Clone for Real Data
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrackReal = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrackReal.TrackProducer = 'ctfWithMaterialTracksP5'
SiStripMonitorTrackReal.TrackLabel    = ''
SiStripMonitorTrackReal.Cluster_src   = 'siStripClusters'
SiStripMonitorTrackReal.Mod_On        = True

# Clone for Real Data (Collision)
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrackColl = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrackColl.TrackProducer = 'generalTracks'
SiStripMonitorTrackColl.TrackLabel    = ''
SiStripMonitorTrackColl.Cluster_src   = 'siStripClusters'
SiStripMonitorTrackColl.Mod_On        = True


# Residual Monitor ####
# Clone for Sim Data
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResidualsSim = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
# Clone for Real Data
MonitorTrackResidualsReal = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResidualsReal.Tracks              = 'ctfWithMaterialTracksP5'
MonitorTrackResidualsReal.trajectoryInput     = 'ctfWithMaterialTracksP5'
MonitorTrackResidualsReal.OutputMEsInRootFile = False
# Clone for Real Data
MonitorTrackResidualsColl = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResidualsColl.Tracks              = 'generalTracks'
MonitorTrackResidualsColl.trajectoryInput     = 'generalTracks'
MonitorTrackResidualsColl.OutputMEsInRootFile = False


# Sequences
SiStripSourcesSimData = cms.Sequence(SiStripMonitorDigi*SiStripMonitorCluster*SiStripMonitorTrackSim*MonitorTrackResidualsSim)
SiStripSourcesRealData = cms.Sequence(SiStripMonitorDigi*SiStripMonitorCluster*SiStripMonitorTrackReal*MonitorTrackResidualsReal)
SiStripSourcesRealDataCollision = cms.Sequence(SiStripMonitorDigi*SiStripMonitorCluster*SiStripMonitorTrackColl*MonitorTrackResidualsColl)




