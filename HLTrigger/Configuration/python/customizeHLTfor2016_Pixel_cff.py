import FWCore.ParameterSet.Config as cms

import itertools

# customisation functions for the HLT configuration
from HLTrigger.Configuration.common import *

# import the relevant eras from Configuration.Eras.*
from Configuration.Eras.Modifier_HLT_2016_cff import HLT_2016

# modify the HLT configuration for the Phase I pixel geometry
def customizeHLTfor2016_Pixel(process):

    for esproducer in esproducers_by_type(process,"ClusterShapeHitFilterESProducer"):
         esproducer.PixelShapeFile = 'RecoPixelVertexing/PixelLowPtUtilities/data/pixelShape.par'
    for producer in producers_by_type(process,"SiPixelRawToDigi"):
        if "hlt" in producer.label():
            producer.UsePhase1 = cms.bool( False )
    return process

# attach `modifyHLTPhaseIPixelGeom' to the `2016` era
def modifyHLTfor2016_Pixel(process):
    HLT_2016.toModify(process, customizeHLTfor2016_Pixel)
