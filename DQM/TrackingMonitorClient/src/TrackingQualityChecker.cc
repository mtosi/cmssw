#include "DQM/TrackingMonitorClient/interface/TrackingQualityChecker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/QReport.h"

#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"

#include "CalibTracker/SiStripCommon/interface/TkDetMap.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "DQM/SiStripMonitorClient/interface/SiStripUtility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include <iomanip>
//
// -- Constructor
// 
TrackingQualityChecker::TrackingQualityChecker(edm::ParameterSet const& ps) :
  pSet_(ps)
{
  edm::LogInfo("TrackingQualityChecker") << " Creating TrackingQualityChecker " << "\n" ;

  bookedTrackingStatus_ = false;

  if(!edm::Service<TkDetMap>().isAvailable()){
    edm::LogError("TkHistoMap") <<
      "\n------------------------------------------"
      "\nUnAvailable Service TkHistoMap: please insert in the configuration file an instance like"
      "\n\tprocess.TkDetMap = cms.Service(\"TkDetMap\")"
      "\n------------------------------------------";
  }
  tkDetMap_=edm::Service<TkDetMap>().operator->();

  TopFolderName_ = pSet_.getUntrackedParameter<std::string>("TopFolderName","Tracking");

  TrackingMEs tracking_mes;

  edm::ParameterSet TkPSet;
 
  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackChi2PSet"); 
  tracking_mes.TrackingFlag = 0;
  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
  tracking_mes.LowerCut     = TkPSet.getParameter<double>("LowerCut");
  tracking_mes.UpperCut     = TkPSet.getParameter<double>("UpperCut");
  TrackingMEsMap.insert(std::pair<std::string, TrackingMEs>("Chi2", tracking_mes));

  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackRatePSet"); 
  tracking_mes.TrackingFlag = 0;
  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
  tracking_mes.LowerCut     = TkPSet.getParameter<double>("LowerCut");
  tracking_mes.UpperCut     = TkPSet.getParameter<double>("UpperCut");
  TrackingMEsMap.insert(std::pair<std::string, TrackingMEs>("Rate", tracking_mes));

  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackHitPSet"); 
  tracking_mes.TrackingFlag = 0;
  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
  tracking_mes.LowerCut     = TkPSet.getParameter<double>("LowerCut");
  tracking_mes.UpperCut     = TkPSet.getParameter<double>("UpperCut");
  TrackingMEsMap.insert(std::pair<std::string, TrackingMEs>("RecHits", tracking_mes));

//  // LS analysis
//  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackRateLSPSet"); 
//  tracking_mes.TrackingFlag = 0;
//  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
//  tracking_mes.LowerCut = TkPSet.getParameter<double>("LowerCut");
//  tracking_mes.UpperCut = TkPSet.getParameter<double>("UpperCut");
//  TrackingMEsLSMap.insert(std::pair<std::string, TrackingMEs>("Rate", tracking_mes));
//
//  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackChi2LSPSet"); 
//  tracking_mes.TrackingFlag = 0;
//  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
//  tracking_mes.LowerCut = TkPSet.getParameter<double>("LowerCut");
//  tracking_mes.UpperCut = TkPSet.getParameter<double>("UpperCut");
//  TrackingMEsLSMap.insert(std::pair<std::string, TrackingMEs>("Chi2", tracking_mes));
//
//  TkPSet = pSet_.getParameter<edm::ParameterSet>("TrackHitLSPSet"); 
//  tracking_mes.TrackingFlag = 0;
//  tracking_mes.HistoName    = TkPSet.getParameter<std::string>("Name");
//  tracking_mes.LowerCut = TkPSet.getParameter<double>("LowerCut");
//  tracking_mes.UpperCut = TkPSet.getParameter<double>("UpperCut");
//  TrackingMEsLSMap.insert(std::pair<std::string, TrackingMEs>("RecHits", tracking_mes));

  useGoodTracks_  = pSet_.getUntrackedParameter<bool>("UseGoodTracks", false);
  if (useGoodTracks_) edm::LogInfo("TrackingQualityChecker") <<  " use GoodTrack histograms for certification " << "\n" ;
}
//
// --  Destructor
// 
TrackingQualityChecker::~TrackingQualityChecker() {
  edm::LogInfo("TrackingQualityChecker") << " Deleting TrackingQualityChecker " << "\n" ;
}
//
// -- create reportSummary MEs
//
void TrackingQualityChecker::bookStatus(DQMStore* dqm_store) {

  if (!bookedTrackingStatus_) {
    dqm_store->cd();     
    edm::LogInfo("TrackingQualityChecker") << " booking TrackingQualityStatus" << "\n";

    std::string tracking_dir = "";
    SiStripUtility::getTopFolderPath(dqm_store, TopFolderName_, tracking_dir);
    if (tracking_dir.size() ==  0) tracking_dir = "Tracking"; // is it necessary !?!?
    dqm_store->setCurrentFolder(tracking_dir+"/EventInfo"); 
    
    TrackSummaryReportGlobal = dqm_store->bookFloat("reportSummary");
    
    std::string hname, htitle;
    hname  = "reportSummaryMap";
    htitle = "Tracking Report Summary Map";
    
    size_t nQT = TrackingMEsMap.size();
    TrackSummaryReportMap    = dqm_store->book2D(hname, htitle, nQT,0.5,double(nQT)+0.5,1,0.5,1.5);
    TrackSummaryReportMap->setAxisTitle("Track Quality Type", 1);
    TrackSummaryReportMap->setAxisTitle("QTest Flag", 2);
    size_t ibin =0;
    for (auto me : TrackingMEsMap) {
      ibin++;
      TrackSummaryReportMap->setBinLabel(ibin+1,me.first);
    }

    dqm_store->setCurrentFolder(tracking_dir+"/EventInfo/reportSummaryContents");  
    ibin = 0;
    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      ibin++;
      std::string MEname = it->first;
      it->second.TrackingFlag = dqm_store->bookFloat("Track"+MEname);
      TrackSummaryReportMap->setBinLabel(ibin,MEname);
    }
    bookedTrackingStatus_ = true;
    dqm_store->cd();
  }
}
//
// -- Fill Dummy  Status
//
void TrackingQualityChecker::fillDummyStatus(){
 
  resetStatus();
  if (bookedTrackingStatus_) {  
    TrackSummaryReportGlobal->Fill(-1.0);
    for (int xbin = 1; xbin < TrackSummaryReportMap->getNbinsX()+1; xbin++) {
      for (int ybin = 1; ybin < TrackSummaryReportMap->getNbinsY()+1; ybin++) {
        TrackSummaryReportMap->Fill(xbin, ybin, -1.0);
      }
    }
    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      it->second.TrackingFlag->Fill(-1.0);
    }
    /*
    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsLSMap.begin();
         it != TrackingMEsLSMap.end(); it++) {
      it->second.TrackingFlag->Fill(-1.0);
    }
    */
  }
}
//
// -- Reset Status
//
void TrackingQualityChecker::resetStatus() {
  if (bookedTrackingStatus_) {  
    TrackSummaryReportGlobal->Reset();
    TrackSummaryReportMap->Reset();
    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      it->second.TrackingFlag->Reset();
    }
    /*
    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsLSMap.begin();
         it != TrackingMEsLSMap.end(); it++) {
      it->second.TrackingFlag->Reset();
    }
    */
  }
}
//
// -- Fill Status
//
void TrackingQualityChecker::fillStatus(DQMStore* dqm_store) {

  if (!bookedTrackingStatus_) bookStatus(dqm_store);

  fillDummyStatus();
  fillTrackingStatus(dqm_store);

}

//
// -- Fill Tracking Status
//
void TrackingQualityChecker::fillTrackingStatus(DQMStore* dqm_store) {

  dqm_store->cd();
  //  std::string dir = "Tracking"; 
  if (!SiStripUtility::goToDir(dqm_store, TopFolderName_)) return;
  //  dir = "TrackParameters"; 
  //  if (!SiStripUtility::goToDir(dqm_store, dir)) return;
  
  /*
  std::vector<MonitorElement*> meVec1;
  std::vector<MonitorElement*> meVec2;
  if (useGoodTracks_){
    meVec1 = dqm_store->getContents(dqm_store->pwd()+"/GeneralProperties/GoodTracks");
    meVec2 = dqm_store->getContents(dqm_store->pwd()+"/HitProperties/GoodTracks");
  }else{
    meVec1 = dqm_store->getContents(dqm_store->pwd()+"/GeneralProperties");
    meVec2 = dqm_store->getContents(dqm_store->pwd()+"/HitProperties");
  }
  std::vector<MonitorElement*> meVec(meVec1.size() + meVec2.size()); 
  std::merge(meVec1.begin(), meVec1.end(), meVec2.begin(), meVec2.end(), meVec.begin());
  */
  std::vector<std::string> MEsDirectories = pSet_.getParameter<std::vector<std::string> >("MEsDirectories");
  std::vector<MonitorElement*> MEs;
  for ( auto medir : MEsDirectories ) {
    std::vector<MonitorElement*> tmpMEs = dqm_store->getContents(dqm_store->pwd()+"/"+medir);
    MEs.insert(MEs.end(),tmpMEs.begin(),tmpMEs.end());
  }

  float gstatus = 1.0;
  for (std::vector<MonitorElement*>::const_iterator itME = MEs.begin(); itME != MEs.end(); itME++) {
    MonitorElement * me = (*itME);     
    if (!me) continue;     
    std::vector<QReport *> qt_reports = me->getQReports();          
    if (qt_reports.size() == 0) continue;
    std::string name = me->getName();

    float status = 1.0; 

    int ibin = 0;
    for (std::map<std::string, TrackingMEs>::const_iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      ibin++;
      std::string hname = it->second.HistoName;
      if (name.find(hname) != std::string::npos) {

	status = qt_reports[0]->getQTresult();
	
	it->second.TrackingFlag->Fill(status);
	fillStatusHistogram(TrackSummaryReportMap, ibin, 1, status);
        break;
      }
    }
    if ( status < 0. ) gstatus = -1.;
    else gstatus = gstatus * status; 
  }

  TrackSummaryReportGlobal->Fill(gstatus);
  dqm_store->cd();
}

//
// -- Fill Report Summary Map
//
 void TrackingQualityChecker::fillStatusHistogram(MonitorElement* me, int xbin, int ybin, float val){
   if (me &&  me->kind() == MonitorElement::DQM_KIND_TH2F) {
     TH2F*  th2d = me->getTH2F();
     th2d->SetBinContent(xbin, ybin, val);
   }
 }

//
// -- Fill Status information and the lumi block
//
void TrackingQualityChecker::fillStatusAtLumi(DQMStore* dqm_store){
  if (!bookedTrackingStatus_) bookStatus(dqm_store);
  fillDummyStatus();
  fillTrackingStatusAtLumi(dqm_store);
}
//
// Fill Tracking Status MEs at the Lumi block
// 
void TrackingQualityChecker::fillTrackingStatusAtLumi(DQMStore* dqm_store){
  dqm_store->cd();
  //  std::string dir = "Tracking"; 
  if (!SiStripUtility::goToDir(dqm_store, TopFolderName_)) return;
  //  dir = "TrackParameters"; 
  //  if (!SiStripUtility::goToDir(dqm_store, dir)) return;

  std::vector<std::string> MEsDirectories = pSet_.getParameter<std::vector<std::string> >("atLumiMEsDirectories");
  std::vector<MonitorElement*> MEs;
  for ( auto medir : MEsDirectories ) {
    std::vector<MonitorElement*> tmpMEs = dqm_store->getContents(dqm_store->pwd()+"/"+medir);
    MEs.insert(MEs.end(),tmpMEs.begin(),tmpMEs.end());
  }

  /*
  std::vector<MonitorElement*> meVec1;
  std::vector<MonitorElement*> meVec2;
  if (useGoodTracks_){
    meVec1 = dqm_store->getContents(dqm_store->pwd()+"/LSanalysis");
  }else{
    meVec1 = dqm_store->getContents(dqm_store->pwd()+"/GeneralProperties");
    meVec2 = dqm_store->getContents(dqm_store->pwd()+"/HitProperties");
  }
  std::vector<MonitorElement*> meVec(meVec1.size() + meVec2.size()); 
  std::merge(meVec1.begin(), meVec1.end(), meVec2.begin(), meVec2.end(), meVec.begin());
  */

  float gstatus = 1.0;
  for (std::vector<MonitorElement*>::const_iterator itME = MEs.begin(); itME != MEs.end(); itME++) {
    MonitorElement * me = (*itME);     
    if (!me) continue;     
    std::string name = me->getName();

    float status = -1.0; 
    int ibin = 0;
    for (std::map<std::string, TrackingMEs>::const_iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      ibin++;
      std::string hname = it->second.HistoName+"lumiFlag_";
      float lower_cut = it->second.LowerCut; 
      float upper_cut = it->second.UpperCut; 
      if (name.find(hname) != std::string::npos) {
        if (me->getMean() <= lower_cut || me->getMean() > upper_cut) status = 0.0;
        else status = 1.0; 

	it->second.TrackingFlag->Fill(status);
	fillStatusHistogram(TrackSummaryReportMap, ibin, 1, status);
        break;
      } else {
      }
    }
    if (status == -1.0) gstatus = -1.0;
    else gstatus = gstatus * status; 
  }
  TrackSummaryReportGlobal->Fill(gstatus);
  dqm_store->cd();
}
