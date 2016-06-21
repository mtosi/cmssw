// C++ headers
#include <string>
#include <cstring>

// boost headers
#include <boost/regex.hpp>

// Root headers
#include <TH1F.h>

// CMSSW headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/MonitorElement.h"

class FastTimerServiceClient : public DQMEDHarvester {
public:
  explicit FastTimerServiceClient(edm::ParameterSet const &);
  ~FastTimerServiceClient();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  static void fillLumiMePSetDescription(edm::ParameterSetDescription & pset);

private:
  std::string m_dqm_path;

  void dqmEndLuminosityBlock(DQMStore::IBooker & booker, DQMStore::IGetter & getter, edm::LuminosityBlock const &, edm::EventSetup const&) override;
  void dqmEndJob(DQMStore::IBooker & booker, DQMStore::IGetter & getter) override;

  void fillSummaryPlots(        DQMStore::IBooker & booker, DQMStore::IGetter & getter);
  void fillProcessSummaryPlots( DQMStore::IBooker & booker, DQMStore::IGetter & getter, std::string const & path);
  void fillPathSummaryPlots(    DQMStore::IBooker & booker, DQMStore::IGetter & getter, double events, std::string const & path);
  void fillPlotsVsLumi(DQMStore::IBooker & booker, DQMStore::IGetter & getter, std::string const & current_path, std::string const & suffix, edm::ParameterSet conf);

  bool doPlotsVsSCALlumi_;
  bool doPlotsVsPIXELlumi_;

  edm::ParameterSet scallumiMEPSet_; 
  edm::ParameterSet pixellumiMEPSet_;

};


FastTimerServiceClient::FastTimerServiceClient(edm::ParameterSet const & config) :
  m_dqm_path( config.getUntrackedParameter<std::string>( "dqmPath" ) )
  , doPlotsVsSCALlumi_ ( config.getParameter<bool>( "doPlotsVsSCALlumi" )  )
  , doPlotsVsPIXELlumi_( config.getParameter<bool>( "doPlotsVsPIXELlumi" ) )
  , scallumiMEPSet_ ( config.getParameter<edm::ParameterSet>("scallumiME")  )
  , pixellumiMEPSet_( config.getParameter<edm::ParameterSet>("pixellumiME") )
{
}

FastTimerServiceClient::~FastTimerServiceClient()
{
}

void
FastTimerServiceClient::dqmEndJob(DQMStore::IBooker & booker, DQMStore::IGetter & getter)
{
  fillSummaryPlots(booker, getter);
}

void
FastTimerServiceClient::dqmEndLuminosityBlock(DQMStore::IBooker & booker, DQMStore::IGetter & getter, edm::LuminosityBlock const & lumi, edm::EventSetup const & setup)
{
  fillSummaryPlots(booker, getter);
}

void
FastTimerServiceClient::fillSummaryPlots(DQMStore::IBooker & booker, DQMStore::IGetter & getter)
{
  if (getter.get(m_dqm_path + "/event")) {
    // the plots are directly in the configured folder
    fillProcessSummaryPlots(booker, getter, m_dqm_path);
  } else {
    static const boost::regex running_n_processes(".*/Running [0-9]+ processes");
    booker.setCurrentFolder(m_dqm_path);
    std::vector<std::string> subdirs = getter.getSubdirs();
    for (auto const & subdir: subdirs) {
      if (boost::regex_match(subdir, running_n_processes)) {
        // the plots are in a per-number-of-processes folder
        if (getter.get(subdir + "/event"))
          fillProcessSummaryPlots(booker, getter, subdir);
      }
    }
  }
}



void
FastTimerServiceClient::fillProcessSummaryPlots(DQMStore::IBooker & booker, DQMStore::IGetter & getter, std::string const & current_path) {

  MonitorElement * me = getter.get(current_path + "/event");
  if (me == 0)
    // no FastTimerService DQM information
    return;

  if ( doPlotsVsSCALlumi_ ) {
    fillPlotsVsLumi( booker,getter, current_path, "VsSCALlumi", scallumiMEPSet_ );
  }
  if ( doPlotsVsPIXELlumi_ ) {
    fillPlotsVsLumi( booker,getter, current_path, "VsPIXELlumi", pixellumiMEPSet_ );
  }

  getter.setCurrentFolder(current_path);

  double events = me->getTH1F()->GetEntries();

  // look for per-process directories
  static const boost::regex process_name(".*/process .*");
  booker.setCurrentFolder(current_path);
  std::vector<std::string> subdirs = getter.getSubdirs();
  for (auto const & subdir: subdirs) {
    if (boost::regex_match(subdir, process_name)) {
      // look for per-path plots inside each per-process directory
      if (getter.dirExists(subdir + "/Paths"))
        fillPathSummaryPlots(booker, getter, events, subdir);
    }
  }

  // look for per-path plots inside the current directory
  if (getter.dirExists(current_path + "/Paths"))
    fillPathSummaryPlots(booker, getter, events, current_path);
}

void
FastTimerServiceClient::fillPathSummaryPlots(DQMStore::IBooker & booker, DQMStore::IGetter & getter, double events, std::string const & current_path) {
  // note: the following checks need to be kept separate, as any of these histograms might be missing
  // if any of them is filled, size will have the total number of paths, and "paths" can be used to extract the list of labels
  MonitorElement * me;
  TProfile const * paths = nullptr;
  uint32_t         size  = 0;

  // extract the list of Paths and EndPaths from the summary plots
  if (( me = getter.get(current_path + "/paths_active_time") )) {
    paths = me->getTProfile();
    size  = paths->GetXaxis()->GetNbins();
  } else
  if (( me = getter.get(current_path + "/paths_total_time") )) {
    paths = me->getTProfile();
    size  = paths->GetXaxis()->GetNbins();
  } else
  if (( me = getter.get(current_path + "/paths_exclusive_time") )) {
    paths = me->getTProfile();
    size  = paths->GetXaxis()->GetNbins();
  }
  if (paths == nullptr)
    return;

  // for each path, fill histograms with
  //  - the average time spent in each module (total time spent in that module, averaged over all events)
  //  - the running time spent in each module (total time spent in that module, averaged over the events where that module actually ran)
  //  - the "efficiency" of each module (number of time a module succeded divided by the number of times the has run)
  booker.setCurrentFolder(current_path + "/Paths");
  for (uint32_t p = 1; p <= size; ++p) {
    // extract the list of Paths and EndPaths from the bin labels of one of the summary plots
    std::string label = paths->GetXaxis()->GetBinLabel(p);
    MonitorElement * me_counter = getter.get( current_path + "/Paths/" + label + "_module_counter" );
    MonitorElement * me_total   = getter.get( current_path + "/Paths/" + label + "_module_total" );
    if (me_counter == 0 or me_total == 0)
      continue;
    TH1F * counter = me_counter->getTH1F();
    TH1F * total   = me_total  ->getTH1F();
    uint32_t bins = counter->GetXaxis()->GetNbins() - 1;
    double   min  = counter->GetXaxis()->GetXmin();
    double   max  = counter->GetXaxis()->GetXmax() - 1;
    booker.setCurrentFolder(current_path + "/Paths");

    TH1F * average;
    TH1F * running;
    TH1F * efficiency;
    MonitorElement * me;

    me = getter.get( current_path + "/Paths/" + label + "_module_average" );
    if (me) {
      average = me->getTH1F();
      //assert( me->getTH1F()->GetXaxis()->GetNbins() == (int) bins );
      assert( me->getTH1F()->GetXaxis()->GetXmin()  == min );
      assert( me->getTH1F()->GetXaxis()->GetXmax()  == max );
      average->Reset();
    } else {
      average = booker.book1D(label + "_module_average", label + " module average", bins, min, max)->getTH1F();
      average->SetYTitle("processing time [ms]");
      for (uint32_t i = 1; i <= bins; ++i) {
        const char * module = counter->GetXaxis()->GetBinLabel(i);
        average->GetXaxis()->SetBinLabel(i, module);
      }
    }

    me = getter.get( current_path + "/Paths/" + label + "_module_running" );
    if (me) {
      running = me->getTH1F();
      //assert( me->getTH1F()->GetXaxis()->GetNbins() == (int) bins );
      assert( me->getTH1F()->GetXaxis()->GetXmin()  == min );
      assert( me->getTH1F()->GetXaxis()->GetXmax()  == max );
      running->Reset();
    } else {
      running = booker.book1D(label + "_module_running", label + " module running", bins, min, max)->getTH1F();
      running->SetYTitle("processing time [ms]");
      for (uint32_t i = 1; i <= bins; ++i) {
        const char * module = counter->GetXaxis()->GetBinLabel(i);
        running->GetXaxis()->SetBinLabel(i, module);
      }
    }

    me = getter.get( current_path + "/Paths/" + label + "_module_efficiency" );
    if (me) {
      efficiency = me->getTH1F();
      //assert( me->getTH1F()->GetXaxis()->GetNbins() == (int) bins );
      assert( me->getTH1F()->GetXaxis()->GetXmin()  == min );
      assert( me->getTH1F()->GetXaxis()->GetXmax()  == max );
      efficiency->Reset();
    } else {
      efficiency = booker.book1D(label + "_module_efficiency", label + " module efficiency", bins, min, max)->getTH1F();
      efficiency->SetYTitle("filter efficiency");
      efficiency->SetMaximum(1.05);
      for (uint32_t i = 1; i <= bins; ++i) {
        const char * module = counter->GetXaxis()->GetBinLabel(i);
        efficiency->GetXaxis()->SetBinLabel(i, module);
      }
    }

    for (uint32_t i = 1; i <= bins; ++i) {
      double t = total  ->GetBinContent(i);
      double n = counter->GetBinContent(i);
      double p = counter->GetBinContent(i+1);
      average ->SetBinContent(i, t / events);
      if (n) {
        running   ->SetBinContent(i, t / n);
        efficiency->SetBinContent(i, p / n);
      }
    }
  }

}

#include "FWCore/MessageLogger/interface/MessageLogger.h"
void
FastTimerServiceClient::fillPlotsVsLumi(DQMStore::IBooker & booker, DQMStore::IGetter & getter, std::string const & current_path, std::string const & suffix, edm::ParameterSet conf ) {

  std::vector<std::string> menames;

  static const boost::regex byls(".*byls");
  static const boost::regex test(suffix);
  // get all MEs in the current_path
  getter.setCurrentFolder(current_path);
  std::vector<std::string> allmenames = getter.getMEs();
  for ( auto m : allmenames )
    // get only MEs vs LS
    if (boost::regex_match(m, byls))
      menames.push_back(m);

  // if no MEs available, return
  if ( menames.size() == 0 )
    return;
  

  // get info for getting the lumi VS LS histogram
  std::string folder = conf.getParameter<std::string>("folder");
  std::string name   = conf.getParameter<std::string>("name");
  int         nbins  = conf.getParameter<int>("nbins");
  double      xmin   = conf.getParameter<int>("xmin");
  double      xmax   = conf.getParameter<int>("xmax");

  // get lumi VS LS ME
  getter.setCurrentFolder(folder);
  MonitorElement* lumiVsLS = getter.get(folder+"/"+name);
  // if no ME available, return
  if ( !lumiVsLS ) {
    edm::LogWarning("FastTimerServiceClient") << "no " << name << " ME is available in " << folder << std::endl;
    return;
  }

  // get range and binning for new MEs x-axis
  size_t size        = lumiVsLS->getTProfile()->GetXaxis()->GetNbins();
  //  double xmin        = lumiVsLS->getTProfile()->GetMinimum();
  //  double xmax        = lumiVsLS->getTProfile()->GetMaximum()*1.1;
  std::string xtitle = lumiVsLS->getTProfile()->GetYaxis()->GetTitle();

  std::vector<double> lumi;
  std::vector<int> LS;
  for ( size_t ibin=1; ibin <= size; ++ibin ) {
    // avoid to store points w/ no info
    if ( lumiVsLS->getTProfile()->GetBinContent(ibin) == 0. ) continue;

    lumi.push_back( lumiVsLS->getTProfile()->GetBinContent(ibin) );
    LS.push_back  ( lumiVsLS->getTProfile()->GetXaxis()->GetBinCenter(ibin) );
  }

  booker.setCurrentFolder(current_path);
  getter.setCurrentFolder(current_path);
  for ( auto m : menames ) {
    std::string label = m;
    label.erase(label.find("_byls"));

    MonitorElement* me = getter.get(current_path + "/" + m);
    double ymin        = me->getTProfile()->GetMinimum();
    double ymax        = me->getTProfile()->GetMaximum();
    std::string ytitle = me->getTProfile()->GetYaxis()->GetTitle();

    MonitorElement* meVsLumi = getter.get( current_path + "/" + label + "_" + suffix );
    if (meVsLumi) {
      assert( meVsLumi->getTProfile()->GetXaxis()->GetXmin()  == xmin );
      assert( meVsLumi->getTProfile()->GetXaxis()->GetXmax()  == xmax );
      meVsLumi->Reset(); // do I have to do it ?!?!?
    } else {
      meVsLumi = booker.bookProfile(label + "_" + suffix, label + "_" + suffix, nbins, xmin, xmax, ymin, ymax);
      //    TProfile* meVsLumi_p = meVsLumi->getTProfile();
      meVsLumi->getTProfile()->GetXaxis()->SetTitle(xtitle.c_str());
      meVsLumi->getTProfile()->GetYaxis()->SetTitle(ytitle.c_str());
    }
    for ( size_t ils=0; ils < LS.size(); ++ils ) {
      int ibin = me->getTProfile()->GetXaxis()->FindBin(LS[ils]);
      double y = me->getTProfile()->GetBinContent(ibin);

      meVsLumi->Fill(lumi[ils],y);
    }
  }
  
  
}

void 
FastTimerServiceClient::fillLumiMePSetDescription(edm::ParameterSetDescription & pset) {
  pset.add<std::string>("folder", "HLT/LumiMonitoring");
  pset.add<std::string>("name"  , "lumiVsLS");
  pset.add<int>   ("nbins",  6500 );
  pset.add<double>("xmin",      0.);
  pset.add<double>("xmax",  13000.);
}


void
FastTimerServiceClient::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>( "dqmPath", "HLT/TimerService");
  desc.add<bool>( "doPlotsVsSCALlumi",  true  );
  desc.add<bool>( "doPlotsVsPIXELlumi", false );

  edm::ParameterSetDescription scallumiMEPSet;
  fillLumiMePSetDescription(scallumiMEPSet);
  desc.add<edm::ParameterSetDescription>("scallumiME", scallumiMEPSet);
  
  edm::ParameterSetDescription pixellumiMEPSet;
  fillLumiMePSetDescription(pixellumiMEPSet);
  desc.add<edm::ParameterSetDescription>("pixellumiME", pixellumiMEPSet);

  descriptions.add("fastTimerServiceClient", desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FastTimerServiceClient);
