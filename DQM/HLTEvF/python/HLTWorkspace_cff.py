import FWCore.ParameterSet.Config as cms

#commented out in 74X
#from DQM.HLTEvF.FourVectorHLTOnline_cfi import *

#from DQM.HLTEvF.OccupancyPlotter_cfi import *

from DQM.HLTEvF.HLTWorkspace_cfi import *
# strip cluster monitor
from DQMOffline.Trigger.SiStrip_OfflineMonitoring_cff import *
# tracking monitor
from DQMOffline.Trigger.TrackingMonitoring_cff import *
iterHLTTracksMonitoringHLT.doProfilesVsLS   = cms.bool(True)

hlt4vector = cms.Path(
    hltWorkspace
#    * sistripMonitorHLTsequence # strip cluster monitoring
    * pixelTracksMonitoringHLT # hltPixel tracks monitoring
    * iterHLTTracksMonitoringHLT # hltIter2Merged tracks monitoring
)

#hlt4vector = cms.Path(onlineOccPlot * hltWorkspace)
#hlt4vector = cms.Path(onlineOccPlot)

