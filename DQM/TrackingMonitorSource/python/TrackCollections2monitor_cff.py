import FWCore.ParameterSet.Config as cms

selectedTracks = [ 'generalTracks' ]

mainfolderName   = {}
vertexfolderName = {}
sequenceName = {}

mainfolderName  ['generalTracks'] = 'Tracking/TrackParameters'
vertexfolderName['generalTracks'] = 'Tracking/PrimaryVertices'

 
trackSelector = cms.EDFilter('TrackSelector',
    src = cms.InputTag('generalTracks'),
    cut = cms.string("")
)

highPurityPtRange0to1 = trackSelector.clone()
highPurityPtRange0to1.cut = cms.string("quality('highPurity') & pt >= 0 & pt < 1 ")

sequenceName    ['highPurityPtRange0to1'] = cms.Sequence(highPurityPtRange0to1)
mainfolderName  ['highPurityPtRange0to1'] = 'Tracking/highPurityTracks/pt_0to1/TrackParameters'
vertexfolderName['highPurityPtRange0to1'] = 'Tracking/highPurityTracks/pt_0to1/PrimaryVertices'


highPurityPtRange1to10 = trackSelector.clone()
highPurityPtRange1to10.cut = cms.string("quality('highPurity') & pt >= 1 & pt < 10 ")

sequenceName     ['highPurityPtRange1to10'] = cms.Sequence( highPurityPtRange1to10 )
mainfolderName  ['highPurityPtRange1to10'] = 'Tracking/highPurityTracks/pt_1to10/TrackParameters'
vertexfolderName['highPurityPtRange1to10'] = 'Tracking/highPurityTracks/pt_1to10/PrimaryVertices'



highPurityPt10 = trackSelector.clone()
highPurityPt10.cut = cms.string("quality('highPurity') & pt >= 10")

sequenceName     ['highPurityPt10'] = cms.Sequence( highPurityPt10 )
mainfolderName  ['highPurityPt10'] = 'Tracking/highPurityTracks/pt_10/TrackParameters'
vertexfolderName['highPurityPt10'] = 'Tracking/highPurityTracks/pt_10/PrimaryVertices'


selectedTracks.extend( ['highPurityPtRange0to1']  )
selectedTracks.extend( ['highPurityPtRange1to10'] )
selectedTracks.extend( ['highPurityPt10']         )

selectedTracks2runSequence=cms.Sequence()
for tracks in selectedTracks :
    if tracks != 'generalTracks':
        selectedTracks2runSequence+=sequenceName[tracks]

