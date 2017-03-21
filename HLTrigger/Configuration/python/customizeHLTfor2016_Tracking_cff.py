import FWCore.ParameterSet.Config as cms

import itertools

# customisation functions for the HLT configuration
from HLTrigger.Configuration.common import *

# import the relevant eras from Configuration.Eras.*
from Configuration.Eras.Modifier_HLT_2016_cff import HLT_2016

# modify the HLT configuration to run the 2016 tracking in the particle flow sequence
def customizeHLTfor2016_Tracking(process):

    ########    ########    ########    ########
    #                pixelTracks               #
    ########    ########    ########    ########
    if 'hltPixelLayerTriplets' in process.__dict__:
        process.hltPixelLayerTriplets.layerList = cms.vstring(
            'BPix1+BPix2+BPix3',
            'BPix1+BPix2+FPix1_pos',
            'BPix1+BPix2+FPix1_neg',
            'BPix1+FPix1_pos+FPix2_pos',
            'BPix1+FPix1_neg+FPix2_neg'
        )

    if 'hltPixelTracksTrackingRegions' in process.__dict__:    
        process.hltPixelTracksTrackingRegions.RegionPSet.nSigmaZ      = cms.double( 0.0 ) ## 4.0 in 2017
        process.hltPixelTracksTrackingRegions.RegionPSet.ptMin        = cms.double( 0.9 ) ## 0.8 in 2017
        process.hltPixelTracksTrackingRegions.RegionPSet.originRadius = cms.double( 0.2 ) ## 0.02 in 2017

    if 'hltPixelTracksHitDoublets' in process.__dict__:    
        process.hltPixelTracksHitDoublets.seedingLayers = "hltPixelLayerTriplets" ## 
        process.hltPixelTracksHitDoublets.layerPairs = [0] ## [0,1,2] in 2017 ( BUT I'M CONFUSED ! [0] or [0,1] ?!?!? )

    process.hltPixelTracksHitTriplets = cms.EDProducer( "PixelTripletHLTEDProducer",
        useBending = cms.bool( True ),
        useFixedPreFiltering = cms.bool( False ),
        produceIntermediateHitTriplets = cms.bool( False ),
        maxElement = cms.uint32( 100000 ),
        SeedComparitorPSet = cms.PSet( 
          clusterShapeHitFilter = cms.string( "ClusterShapeHitFilter" ),
          ComponentName = cms.string( "LowPtClusterShapeSeedComparitor" ),
          clusterShapeCacheSrc = cms.InputTag( "hltSiPixelClustersCache" )
        ),
        extraHitRPhitolerance = cms.double( 0.06 ),
        produceSeedingHitSets = cms.bool( True ),
        doublets = cms.InputTag( "hltPixelTracksHitDoublets" ),
        useMultScattering = cms.bool( True ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRZtolerance = cms.double( 0.06 )
    )

    if 'hltPixelTracks' in process.__dict__:
        process.hltPixelTracks.SeedingHitSets = "hltPixelTracksHitTriplets"

    
    ########    ########    ########
    #             iter0            #
    ########    ########    ########
    if 'HLTIter0PSetTrajectoryFilterIT' in process.__dict__:
        process.HLTIter0PSetTrajectoryFilterIT.minimumNumberOfHits = cms.int32( 3 ) ## 4 in 2017
        process.HLTIter0PSetTrajectoryFilterIT.minHitsMinPt        = cms.int32( 3 ) ## 4 in 2017
        process.hltIter0PFlowTrackCutClassifier.mva.minLayers    = cms.vint32( 3, 3, 3 ) ## 3, 3, 4 in 2017
        process.hltIter0PFlowTrackCutClassifier.mva.min3DLayers  = cms.vint32( 0, 0, 0 ) ## 0, 3, 4 in 2017
        process.hltIter0PFlowTrackCutClassifier.mva.minPixelHits = cms.vint32( 0, 0, 0 ) ## 0, 3, 4 in 2017

    process.HLTIter0PSetTrajectoryBuilderIT = cms.PSet( 
         propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
         trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
         maxCand = cms.int32( 2 ),
         ComponentType = cms.string( "CkfTrajectoryBuilder" ),
         propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
         estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
         TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
         updator = cms.string( "hltESPKFUpdator" ),
         alwaysUseInvalidHits = cms.bool( False ),
         intermediateCleaning = cms.bool( True ),
         lostHitPenalty = cms.double( 30.0 )
    )

    if 'hltIter0PFlowCkfTrackCandidates' in process.__dict__:
        process.hltIter0PFlowCkfTrackCandidates.TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('HLTIter0PSetTrajectoryBuilderIT'))


    ########    ########    ########
    #             iter1            #
    ########    ########    ########
    process.hltIter1PixelLayerTriplets = cms.EDProducer( "SeedingLayersEDProducer",
        layerList = cms.vstring( 
          'BPix1+BPix2+BPix3',
          'BPix1+BPix2+FPix1_pos',
          'BPix1+BPix2+FPix1_neg',
          'BPix1+FPix1_pos+FPix2_pos',
          'BPix1+FPix1_neg+FPix2_neg' ),
        MTOB = cms.PSet(  ),
        TEC = cms.PSet(  ),
        MTID = cms.PSet(  ),
        FPix = cms.PSet( 
          hitErrorRPhi = cms.double( 0.0051 ),
          TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
          skipClusters = cms.InputTag( "hltIter1ClustersRefRemoval" ),
          useErrorsFromParam = cms.bool( True ),
          hitErrorRZ = cms.double( 0.0036 ),
          HitProducer = cms.string( "hltSiPixelRecHits" )
        ),
        MTEC = cms.PSet(  ),
        MTIB = cms.PSet(  ),
        TID = cms.PSet(  ),
        TOB = cms.PSet(  ),
        BPix = cms.PSet( 
          hitErrorRPhi = cms.double( 0.0027 ),
          TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
          skipClusters = cms.InputTag( "hltIter1ClustersRefRemoval" ),
          useErrorsFromParam = cms.bool( True ),
          hitErrorRZ = cms.double( 0.006 ),
          HitProducer = cms.string( "hltSiPixelRecHits" )
        ),
        TIB = cms.PSet(  )
    )
    
    if 'hltIter1PFlowPixelHitDoublets' in process.__dict__:
        process.hltIter1PFlowPixelHitDoublets.layerPairs = [0] ## [0,1,2] in 2017 ( BUT I'M CONFUSED ! [0] or [0,1] ?!?!? )
        process.hltIter1PFlowPixelHitDoublets.seedingLayers = "hltIter1PixelLayerTriplets"

    if 'hltIter1PFlowPixelHitTriplets' in process.__dict__:
        process.hltIter1PFlowPixelHitTriplets = cms.EDProducer( "PixelTripletHLTEDProducer",
        useBending = cms.bool( True ),
        useFixedPreFiltering = cms.bool( False ),
        produceIntermediateHitTriplets = cms.bool( False ),
        maxElement = cms.uint32( 100000 ),
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        produceSeedingHitSets = cms.bool( True ),
        doublets = cms.InputTag( "hltIter1PFlowPixelHitDoublets" ),
        useMultScattering = cms.bool( True ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRZtolerance = cms.double( 0.037 )
    )

    if 'hltIter1PFlowPixelTrackingRegions' in process.__dict__:
        process.hltIter1PFlowPixelTrackingRegions.RegionPSet.nSigmaZVertex   = cms.double( 3.0 ) ## 4 in 2017
        process.hltIter1PFlowPixelTrackingRegions.RegionPSet.nSigmaZBeamSpot = cms.double( 3.0 ) ## 4 in 2017      
        process.hltIter1PFlowPixelTrackingRegions.RegionPSet.ptMin           = cms.double( 0.5 ) ## 0.3 in 2017

    from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsEDProducer_cff import seedCreatorFromRegionConsecutiveHitsEDProducer as _seedCreatorFromRegionConsecutiveHitsEDProducer
    if 'hltIter1PFlowPixelSeeds' in process.__dict__:
        replace_with(process.hltIter1PFlowPixelSeeds,_seedCreatorFromRegionConsecutiveHitsEDProducer.clone(
                seedingHitSets = "hltIter1PFlowPixelHitTriplets",
                TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        ))

    process.HLTIter1PSetTrajectoryBuilderIT = cms.PSet( 
	propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
	trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter1PSetTrajectoryFilterIT" ) ),
	maxCand = cms.int32( 2 ),
	ComponentType = cms.string( "CkfTrajectoryBuilder" ),
	propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
	MeasurementTrackerName = cms.string( "hltIter1ESPMeasurementTracker" ),
	estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
	updator = cms.string( "hltESPKFUpdator" ),
	alwaysUseInvalidHits = cms.bool( False ),
	intermediateCleaning = cms.bool( True ),
	lostHitPenalty = cms.double( 30.0 ),
	useSameTrajFilter = cms.bool(True) 
    )

    if 'hltIter1PFlowCkfTrackCandidates' in process.__dict__:
        process.hltIter1PFlowCkfTrackCandidates.TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('HLTIter1PSetTrajectoryBuilderIT'))

    if 'HLTIterativeTrackingIteration1' in process.__dict__:
        replace_with(process.HLTIterativeTrackingIteration1 , cms.Sequence( process.hltIter1ClustersRefRemoval + process.hltIter1MaskedMeasurementTrackerEvent + process.hltIter1PixelLayerTriplets + process.hltIter1PFlowPixelTrackingRegions + process.hltIter1PFlowPixelClusterCheck + process.hltIter1PFlowPixelHitDoublets + process.hltIter1PFlowPixelHitTriplets + process.hltIter1PFlowPixelSeeds + process.hltIter1PFlowCkfTrackCandidates + process.hltIter1PFlowCtfWithMaterialTracks + process.hltIter1PFlowTrackCutClassifierPrompt + process.hltIter1PFlowTrackCutClassifierDetached + process.hltIter1PFlowTrackCutClassifierMerged + process.hltIter1PFlowTrackSelectionHighPurity ) )

    ########    ########    ########
    #             iter2            #
    ########    ########    ########
    process.hltIter2PixelLayerPairs = cms.EDProducer( "SeedingLayersEDProducer",
        layerList = cms.vstring( 
          'BPix1+BPix2',
          'BPix1+BPix3',
          'BPix2+BPix3',
          'BPix1+FPix1_pos',
          'BPix1+FPix1_neg',
          'BPix1+FPix2_pos',
          'BPix1+FPix2_neg',
          'BPix2+FPix1_pos',
          'BPix2+FPix1_neg',
          'BPix2+FPix2_pos',
          'BPix2+FPix2_neg',
          'FPix1_pos+FPix2_pos',
          'FPix1_neg+FPix2_neg' ),
        MTOB = cms.PSet(  ),
        TEC = cms.PSet(  ),
        MTID = cms.PSet(  ),
        FPix = cms.PSet( 
          hitErrorRPhi = cms.double( 0.0051 ),
          TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
          skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
          useErrorsFromParam = cms.bool( True ),
          hitErrorRZ = cms.double( 0.0036 ),
          HitProducer = cms.string( "hltSiPixelRecHits" )
        ),
        MTEC = cms.PSet(  ),
        MTIB = cms.PSet(  ),
        TID = cms.PSet(  ),
        TOB = cms.PSet(  ),
        BPix = cms.PSet( 
          hitErrorRPhi = cms.double( 0.0027 ),
          TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
          skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
          useErrorsFromParam = cms.bool( True ),
          hitErrorRZ = cms.double( 0.006 ),
          HitProducer = cms.string( "hltSiPixelRecHits" )
        ),
        TIB = cms.PSet(  )
    )
 
    if 'hltIter2PFlowPixelTrackingRegions' in process.__dict__:
        process.hltIter2PFlowPixelTrackingRegions.RegionPSet.ptMin         = cms.double(1.2) ## 0.8 in 2017
        process.hltIter2PFlowPixelTrackingRegions.RegionPSet.nSigmaZVertex = cms.double(3.0) ## 4.0 in 2017

    if 'hltIter2PFlowPixelHitDoublets' in process.__dict__:
        process.hltIter2PFlowPixelHitDoublets.seedingLayers                  = "hltIter2PixelLayerPairs"
        process.hltIter2PFlowPixelHitDoublets.produceIntermediateHitDoublets = False ## True in 2017
        process.hltIter2PFlowPixelHitDoublets.produceSeedingHitSets          = True ## False in 2017
        process.hltIter2PFlowPixelHitDoublets.layerPairs                     = [0] ## [0,1] in 2017 ( BUT I'M CONFUSED ! [0] or [0,1] ?!?!? )

    def _copy(old, new, skip=[]):
        skipSet = set(skip)
        for key in old.parameterNames_():
            if key not in skipSet:
                setattr(new, key, getattr(old, key))
    from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsEDProducer_cfi import seedCreatorFromRegionConsecutiveHitsEDProducer as _seedCreatorFromRegionConsecutiveHitsEDProducer
    if 'hltIter2PFlowPixelSeeds' in process.__dict__:
        replace_with(process.hltIter2PFlowPixelSeeds, _seedCreatorFromRegionConsecutiveHitsEDProducer.clone(seedingHitSets="hltIter2PFlowPixelHitDoublets"))
    if 'HLTSeedFromConsecutiveHitsCreatorIT' in process.__dict__:
        _copy(process.HLTSeedFromConsecutiveHitsCreatorIT, process.hltIter2PFlowPixelSeeds, skip=["ComponentName"])

    process.HLTIter2PSetTrajectoryBuilderIT = cms.PSet( 
      propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
      trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
      maxCand = cms.int32( 2 ),
      ComponentType = cms.string( "CkfTrajectoryBuilder" ),
      propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
      MeasurementTrackerName = cms.string( "hltIter2ESPMeasurementTracker" ),
      estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
      TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      updator = cms.string( "hltESPKFUpdator" ),
      alwaysUseInvalidHits = cms.bool( False ),
      intermediateCleaning = cms.bool( True ),
      lostHitPenalty = cms.double( 30.0 )
    )

    if 'hltIter2PFlowCkfTrackCandidates' in process.__dict__:
        process.hltIter2PFlowCkfTrackCandidates.TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('HLTIter2PSetTrajectoryBuilderIT'))

    if 'HLTIterativeTrackingIteration2' in process.__dict__:
        replace_with(process.HLTIterativeTrackingIteration2 , cms.Sequence( process.hltIter2ClustersRefRemoval + process.hltIter2MaskedMeasurementTrackerEvent + process.hltIter2PixelLayerPairs + process.hltIter2PFlowPixelTrackingRegions + process.hltIter2PFlowPixelClusterCheck + process.hltIter2PFlowPixelHitDoublets + process.hltIter2PFlowPixelSeeds + process.hltIter2PFlowCkfTrackCandidates + process.hltIter2PFlowCtfWithMaterialTracks + process.hltIter2PFlowTrackCutClassifier + process.hltIter2PFlowTrackSelectionHighPurity ))




    ########    ########    ########
    # replace hltPixelLayerTriplets and hltPixelTracksHitTriplets with hltPixelLayerQuadruplets and hltPixelTracksHitQuadruplets
    # in any Sequence, Paths or EndPath that contains the former and not the latter
    from FWCore.ParameterSet.SequenceTypes import ModuleNodeVisitor
    for sequence in itertools.chain(
        process._Process__sequences.itervalues(),
        process._Process__paths.itervalues(),
        process._Process__endpaths.itervalues()
    ):
        modules = list()
        sequence.visit(ModuleNodeVisitor(modules))

        if process.hltPixelTracks in modules and not process.hltPixelLayerTriplets in modules:
            # note that this module does not necessarily exist in sequence 'sequence', if it doesn't, it does not get removed
            sequence.remove(process.hltPixelLayerQuadruplets)
            index = sequence.index(process.hltPixelTracksHitDoublets)
            sequence.insert(index,process.hltPixelLayerTriplets)
            index = sequence.index(process.hltPixelTracksHitQuadruplets)
            sequence.remove(process.hltPixelTracksHitQuadruplets)
            sequence.insert(index, process.hltPixelTracksHitTriplets)

	if process.hltIter1PFlowPixelHitQuadruplets in modules and not process.hltIter1PFlowPixelHitTriplets in modules:
            index = sequence.index(process.hltIter1PFlowPixelHitQuadruplets)
            sequence.insert(index, process.hltIter1PixelTracks)
            sequence.insert(index, process.hltIter1PFlowPixelHitTriplets)
            sequence.remove(process.hltIter1PFlowPixelHitQuadruplets)


    # Remove entirely to avoid warning from the early deleter
    if 'hltPixelLayerQuadruplets' in process.__dict__:
        del process.hltPixelLayerQuadruplets
    if 'hltPixelTracksHitQuadruplets' in process.__dict__:
        del process.hltPixelTracksHitQuadruplets
    if 'hltIter1PixelLayerQuadruplets' in process.__dict__:
        del process.hltIter1PixelLayerQuadruplets
    if 'hltIter1PFlowPixelHitQuadruplets' in process.__dict__:
        del process.hltIter1PFlowPixelHitQuadruplets
    if 'hltIter2PixelLayerTriplets' in process.__dict__:
        del process.hltIter2PixelLayerTriplets
    if 'hltIter2PFlowPixelHitTriplets' in process.__dict__:
        del process.hltIter2PFlowPixelHitTriplets

    if 'HLTIter0GroupedCkfTrajectoryBuilderIT' in process.__dict__:
        del process.HLTIter0GroupedCkfTrajectoryBuilderIT
    if 'HLTIter1GroupedCkfTrajectoryBuilderIT' in process.__dict__:
        del process.HLTIter1GroupedCkfTrajectoryBuilderIT
    if 'HLTIter2GroupedCkfTrajectoryBuilderIT' in process.__dict__:
        del process.HLTIter2GroupedCkfTrajectoryBuilderIT

    return process

# attach `customizeHLTForPFTrackingPhaseI2017` to the `phase1Pixel` era
def modifyHLTfor2016_Tracking(process):
    HLT_2016.toModify(process, customizeHLTfor2016_Tracking)
