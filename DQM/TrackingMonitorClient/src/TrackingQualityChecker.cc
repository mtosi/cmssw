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

  bookedTrackingStatus_       = false;

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
  std::vector<edm::ParameterSet> TrackingGlobalQualityMEs = pSet_.getParameter< std::vector<edm::ParameterSet> >("TrackingGlobalQualityPSets" );
  for ( auto meQTset : TrackingGlobalQualityMEs ) {

    std::string QTname           = meQTset.getParameter<std::string>("QT");
    tracking_mes.HistoDir        = meQTset.getParameter<std::string>("dir");
    tracking_mes.HistoName       = meQTset.getParameter<std::string>("name");
    tracking_mes.HistoLSDir      = meQTset.getParameter<std::string>("LSdir");
    tracking_mes.HistoLSName     = meQTset.getParameter<std::string>("LSname");
    tracking_mes.HistoLSLowerCut = meQTset.getParameter<double>("LSlowerCut");
    tracking_mes.HistoLSUpperCut = meQTset.getParameter<double>("LSupperCut");
    tracking_mes.TrackingFlag = 0;

    std::cout << "[TrackingQualityChecker::TrackingQualityChecker] inserting " << QTname << " in TrackingMEsMap" << std::endl;
    TrackingMEsMap.insert(std::pair<std::string, TrackingMEs>(QTname, tracking_mes));
  }

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
  
  std::cout << "[TrackingQualityChecker::bookStatus] already booked ? " << (bookedTrackingStatus_ ? "yes" : "nope") << std::endl;

  if (!bookedTrackingStatus_) {
    dqm_store->cd();     
    edm::LogInfo("TrackingQualityChecker") << " booking TrackingQualityStatus" << "\n";

    std::string tracking_dir = "";
    SiStripUtility::getTopFolderPath(dqm_store, TopFolderName_, tracking_dir);
    dqm_store->setCurrentFolder(TopFolderName_+"/EventInfo"); 
    
    TrackSummaryReportGlobal = dqm_store->bookFloat("reportSummary");
    
    std::string hname, htitle;
    hname  = "reportSummaryMap";
    htitle = "Tracking Report Summary Map";
    
    size_t nQT = TrackingMEsMap.size();
    std::cout << "[TrackingQualityChecker::bookStatus] nQT: " << nQT << std::endl;
    TrackSummaryReportMap    = dqm_store->book2D(hname, htitle, nQT,0.5,double(nQT)+0.5,1,0.5,1.5);
    TrackSummaryReportMap->setAxisTitle("Track Quality Type", 1);
    TrackSummaryReportMap->setAxisTitle("QTest Flag", 2);
    size_t ibin =0;
    for ( auto meQTset : TrackingMEsMap ) {
      TrackSummaryReportMap->setBinLabel(ibin+1,meQTset.first);
      ibin++;
    }

    dqm_store->setCurrentFolder(TopFolderName_+"/EventInfo/reportSummaryContents");  

    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      std::string meQTname = it->first;
      it->second.TrackingFlag = dqm_store->bookFloat("Track"+meQTname);
      std::cout << "[TrackingQualityChecker::bookStatus] " << it->first << " exists ? " << it->second.TrackingFlag << std::endl;      
      std::cout << "[TrackingQualityChecker::bookStatus] DONE w/ TrackingMEsMap" << std::endl;
    }

    bookedTrackingStatus_ = true;
    dqm_store->cd();
  }
}

//
// -- Fill Dummy  Status
//
void TrackingQualityChecker::fillDummyStatus(){
  std::cout << "[TrackingQualityChecker::fillDummyStatus] starting ..." << std::endl;

  resetStatus();
  std::cout << "[TrackingQualityChecker::fillDummyStatus] already booked ? " << (bookedTrackingStatus_ ? "yes" : "nope") << std::endl;
  if (bookedTrackingStatus_) {  

    TrackSummaryReportGlobal->Fill(-1.0);

    for (int ibin = 1; ibin < TrackSummaryReportMap->getNbinsX()+1; ibin++) {
      fillStatusHistogram(TrackSummaryReportMap, ibin, 1, -1.0);
    }

    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++)
      it->second.TrackingFlag->Fill(-1.0);
    std::cout << "[TrackingQualityChecker::fillDummyStatus] DONE w/ TrackingMEsMap" << std::endl;

  }
}

//
// -- Reset Status
//
void TrackingQualityChecker::resetStatus() {

  std::cout << "[TrackingQualityChecker::resetStatus] already booked ? " << (bookedTrackingStatus_ ? "yes" : "nope") << std::endl;
  if (bookedTrackingStatus_) {  

    TrackSummaryReportGlobal -> Reset();
    TrackSummaryReportMap    -> Reset();

    for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
         it != TrackingMEsMap.end(); it++) {
      MonitorElement* me = it->second.TrackingFlag;
      std::cout << "[TrackingQualityChecker::resetStatus] " << it->second.HistoName << " exist ? " << ( it->second.TrackingFlag == NULL ? "nope" : "yes" ) << " ---> " << me << std::endl;      
      it->second.TrackingFlag->Reset();
    }
    std::cout << "[TrackingQualityChecker::resetStatus] DONE w/ TrackingMEsMap" << std::endl;

  }
}

//
// -- Fill Status
//
void TrackingQualityChecker::fillStatus(DQMStore* dqm_store) {

  std::cout << "[TrackingQualityChecker::fillStatus] already booked ? " << (bookedTrackingStatus_ ? "yes" : "nope") << std::endl;
  if (!bookedTrackingStatus_) bookStatus(dqm_store);

  fillDummyStatus();
  fillTrackingStatus(dqm_store);
  std::cout << "[TrackingQualityChecker::fillStatus] DONE" << std::endl;
  dqm_store->cd();
}

//
// -- Fill Tracking Status
//
void TrackingQualityChecker::fillTrackingStatus(DQMStore* dqm_store) {

  float gstatus = 0.0;

  dqm_store->cd();
  if (!SiStripUtility::goToDir(dqm_store, TopFolderName_)) return;
  
    
  int ibin = 0;
  for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
       it != TrackingMEsMap.end(); it++) {

    std::cout << "[TrackingQualityChecker::fillTrackingStatus] ME: " << it->first << " [" << it->second.TrackingFlag->getFullname() << "] flag: " << it->second.TrackingFlag->getFloatValue() << std::endl;

    ibin++;
    
    std::string localMEdirpath = it->second.HistoDir;
    std::string MEname         = it->second.HistoName;

    std::vector<MonitorElement*> tmpMEvec = dqm_store->getContents(dqm_store->pwd()+"/"+localMEdirpath);
    MonitorElement* me = NULL;
    float status = 0.;
    for ( auto ime : tmpMEvec ) {
      std::string name = ime->getName();
      if ( name.find(MEname) != std::string::npos) {
	me = ime;
	break;
      }
    }
    if (!me) continue;
    std::cout << "[TrackingQualityChecker::fillTrackingStatus] status: " << status << std::endl;
    std::vector<QReport *> qt_reports = me->getQReports();          
    size_t nQTme = qt_reports.size();
    if (nQTme != 0) {
      std::cout << "[TrackingQualityChecker::fillTrackingStatus] qt_reports: " << qt_reports.size() << std::endl;
      for ( auto iQT : qt_reports ) {
	status += iQT->getQTresult();
	std::cout << "[TrackingQualityChecker::fillTrackingStatus] iQT: " << iQT->getQRName() << std::endl;
	std::cout << "[TrackingQualityChecker::fillTrackingStatus] MEname: " << MEname << " status: " << iQT->getQTresult() << " exists ? " << (it->second.TrackingFlag ? "yes " : "no ") << it->second.TrackingFlag << std::endl;
      }
      status = status/float(nQTme);
      std::cout << "[TrackingQualityChecker::fillTrackingStatus] MEname: " << MEname << " status: " << status << std::endl;
      it->second.TrackingFlag->Fill(status);
      std::cout << "[TrackingQualityChecker::fillTrackingStatus] TrackSummaryReportMap: " << TrackSummaryReportMap << std::endl;
      fillStatusHistogram(TrackSummaryReportMap, ibin, 1, status);
    }
   
    std::cout << "[TrackingQualityChecker::fillTrackingStatus] gstatus: " << gstatus << " x status: " << status << std::endl;
    if ( status < 0. ) gstatus = -1.;
    else gstatus += status; 
    std::cout << "[TrackingQualityChecker::fillTrackingStatus] ===> gstatus: " << gstatus << std::endl;
    std::cout << "[TrackingQualityChecker::fillTrackingStatus] ME: " << it->first << " [" << it->second.TrackingFlag->getFullname() << "] flag: " << it->second.TrackingFlag->getFloatValue() << std::endl;
  }
  
  size_t nQT = TrackingMEsMap.size();
  if (gstatus < 1.) gstatus = -1.;
  else gstatus = gstatus/float(nQT);

  std::cout << "[TrackingQualityChecker::fillTrackingStatus] ===> gstatus: " << gstatus << std::endl;
  TrackSummaryReportGlobal->Fill(gstatus);
  dqm_store->cd();

  std::cout << "[TrackingQualityChecker::fillTrackingStatus] DONE" << std::endl;

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
  std::cout << "[TrackingQualityChecker::fillStatusAtLumi] starting .." << std::endl;
  if (!bookedTrackingStatus_) bookStatus(dqm_store);
  fillDummyStatus();
  fillTrackingStatusAtLumi(dqm_store);
  std::cout << "[TrackingQualityChecker::fillStatusAtLumi] DONE" << std::endl; 

  for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
       it != TrackingMEsMap.end(); it++) {
    std::cout << "[TrackingQualityChecker::fillStatusAtLumi] ME: " << it->first << " [" << it->second.TrackingFlag->getFullname() << "] flag: " << it->second.TrackingFlag->getFloatValue() << std::endl;
  }
}
//
// Fill Tracking Status MEs at the Lumi block
// 
void TrackingQualityChecker::fillTrackingStatusAtLumi(DQMStore* dqm_store){

  std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] starting .. " << std::endl;
  float gstatus = 1.0;

  dqm_store->cd();
  if (!SiStripUtility::goToDir(dqm_store, TopFolderName_)) return;


  int ibin = 0;
  for (std::map<std::string, TrackingMEs>::iterator it = TrackingMEsMap.begin();
       it != TrackingMEsMap.end(); it++) {
    
    std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] ME: " << it->first << " [" << it->second.TrackingFlag->getFullname() << "] flag: " << it->second.TrackingFlag->getFloatValue() << std::endl;

    ibin++;
  
    std::string localMEdirpath = it->second.HistoLSDir;
    std::string MEname         = it->second.HistoLSName;
    float lower_cut            = it->second.HistoLSLowerCut; 
    float upper_cut            = it->second.HistoLSUpperCut; 

    float status = 1.0; 

    std::vector<MonitorElement*> tmpMEvec = dqm_store->getContents(dqm_store->pwd()+"/"+localMEdirpath);
    MonitorElement* me = NULL;
    for ( auto ime : tmpMEvec ) {
      std::string name = ime->getName();
      std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] name: " << name << std::endl;
      if (name.find(MEname) != std::string::npos) {
	me = ime;
	break;
      }
    }
    if (!me) continue;     

    if (me->kind() == MonitorElement::DQM_KIND_TH1F) {
      double x_mean = me->getMean();
      std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] MEname: " << MEname << " x_mean: " << x_mean << std::endl;
      if (x_mean <= lower_cut || x_mean > upper_cut) status = 0.0;
      else status = 1.0; 
      
      it->second.TrackingFlag->Fill(status);
      fillStatusHistogram(TrackSummaryReportMap, ibin, 1, status);
      std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] ===> status: " << status << " [" << gstatus << "]" << std::endl;
    }
    if (status == 0.0) gstatus = -1.0;
    else gstatus = gstatus * status; 
    std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] ===> gstatus: " << gstatus << std::endl;
    std::cout << "[TrackingQualityChecker::fillTrackingStatusAtLumi] ME: " << it->first << " [" << it->second.TrackingFlag->getFullname() << "] flag: " << it->second.TrackingFlag->getFloatValue() << std::endl;
  }
  TrackSummaryReportGlobal->Fill(gstatus);
  dqm_store->cd();
}
