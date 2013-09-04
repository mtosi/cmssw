#ifndef _TrackingQualityChecker_h_
#define _TrackingQualityChecker_h_

#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

class DQMStore;
class MonitorElement;
class TkDetMap;
class TrackingDetCabling;

class TrackingQualityChecker {

 public:


  TrackingQualityChecker(edm::ParameterSet const& ps);
  virtual ~TrackingQualityChecker();


  void bookStatus(DQMStore* dqm_store);     
  void resetStatus();
  void fillDummyStatus();
  void fillStatus(DQMStore* dqm_store);
  void fillStatusAtLumi(DQMStore* dqm_store);
  void fillFaultyModuleStatus(DQMStore* dqm_store, const edm::EventSetup& eSetup);
  
 private:

  struct TrackingMEs{
    MonitorElement* TrackingFlag;
    std::string     HistoDir;
    std::string     HistoName;
    std::string     HistoLSDir;
    std::string     HistoLSName;
    float           HistoLSLowerCut;
    float           HistoLSUpperCut; 
  };

  void fillTrackingStatus(DQMStore* dqm_store);
  void getModuleStatus(DQMStore* dqm_store, std::vector<MonitorElement*>& layer_mes, int& errdet);

  void fillStatusHistogram(MonitorElement*, int xbin, int ybin, float val);

  void fillTrackingStatusAtLumi(DQMStore* dqm_store);
  
  std::map<std::string, TrackingMEs> TrackingMEsMap;
  std::map<std::string, TrackingMEs> TrackingMEsLSMap;
  
  MonitorElement* TrackSummaryReportMap;
  MonitorElement* TrackSummaryReportGlobal;

  edm::ParameterSet pSet_;

  bool bookedTrackingStatus_;
  int globalStatusFilling_;
  bool useGoodTracks_;

  TkDetMap* tkDetMap_;
 
  float cutoffTrackRate_;
  float cutoffChi2overDoF_;
  float cutoffRecHits_;

  std::string TopFolderName_;

};
#endif
