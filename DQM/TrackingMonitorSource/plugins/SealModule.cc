#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DQM/TrackingMonitorSource/interface/Track2TrackValidator.h"
DEFINE_FWK_MODULE(Track2TrackValidator);

#include "DQM/TrackingMonitorSource/interface/Vertex2VertexValidator.h"
DEFINE_FWK_MODULE(Vertex2VertexValidator);
