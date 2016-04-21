import FWCore.ParameterSet.Config as cms

from DQM.HLTEvF.triggerBxMonitor_cfi import *
from DQM.HLTEvF.triggerRatesMonitor_cfi import *

allTriggerPathMonitoringVsBX = triggerBxMonitor.clone(
  doPlots4l1 = cms.untracked.bool(False),
  doPlots4hlt = cms.untracked.bool(True),
  l1tResults = cms.untracked.InputTag('hltGtDigis', '', 'HLT'),
  hltResults = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  dqmPath = cms.untracked.string('HLT/TriggerBx'),
  hltTriggerPaths = cms.untracked.vstring()
)

zeroBiasMonitoringVsBX = triggerBxMonitor.clone(
  doPlots4l1      = cms.untracked.bool(False),
  doPlots4hlt     = cms.untracked.bool(True),
  l1tResults      = cms.untracked.InputTag('hltGtDigis', '', 'HLT'),
  hltResults      = cms.untracked.InputTag('TriggerResults', '', 'HLT'),
  dqmPath         = cms.untracked.string('HLT/TriggerBx'),
  hltTriggerPaths = cms.untracked.vstring( 'HLT_ZeroBias_v', 'HLT_ZeroBias_IsolatedBunches_v', 'HLT_ZeroBias_FirstCollisionAfterAbortGap' )
)

