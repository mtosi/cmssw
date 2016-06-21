#include "DQM/HLTEvF/plugins/LumiMonitor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQM/TrackingMonitor/interface/GetLumi.h"


// -----------------------------
//  constructors and destructor
// -----------------------------

LumiMonitor::LumiMonitor( const edm::ParameterSet& iConfig ) : 
  folderName_                     ( iConfig.getParameter<std::string>("FolderName")                  )
  , pixelClustersToken_           ( mayConsume<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("pixelClusters") ) )
  , lumiScalersToken_             ( consumes<LumiScalersCollection>                  (iConfig.getParameter<edm::InputTag>("scalers")       ) ) 
  , doPixelLumi_                  ( iConfig.getParameter<bool>       ("doPixelLumi")                 )
  , useBPixLayer1_                ( iConfig.getParameter<bool>       ("useBPixLayer1")               )
  , minNumberOfPixelsPerCluster_  ( iConfig.getParameter<int>        ("minNumberOfPixelsPerCluster") )
  , minPixelClusterCharge_        ( iConfig.getParameter<double>     ("minPixelClusterCharge")       )
{

  numberOfPixelClustersVsLS_   = nullptr;
  numberOfPixelClustersVsLumi_ = nullptr;
  lumiVsLS_                    = nullptr;
  pixelLumiVsLS_               = nullptr;
  pixelLumiVsLumi_             = nullptr;

  if(useBPixLayer1_) 
    lumi_factor_per_bx_ = GetLumi::FREQ_ORBIT * GetLumi::SECONDS_PER_LS / GetLumi::XSEC_PIXEL_CLUSTER  ;
  else
    lumi_factor_per_bx_ = GetLumi::FREQ_ORBIT * GetLumi::SECONDS_PER_LS / GetLumi::rXSEC_PIXEL_CLUSTER  ;

  edm::ParameterSet histoPSet        = iConfig.getParameter<edm::ParameterSet>("histoPSet");
  edm::ParameterSet pixelClusterPSet = histoPSet.getParameter<edm::ParameterSet>("pixelClusterPSet");
  edm::ParameterSet lumiPSet         = histoPSet.getParameter<edm::ParameterSet>("lumiPSet");
  edm::ParameterSet pixellumiPSet    = histoPSet.getParameter<edm::ParameterSet>("pixellumiPSet");
  edm::ParameterSet lsPSet           = histoPSet.getParameter<edm::ParameterSet>("lsPSet");

  /*
  pixelCluster_binning_(pixelClusterPSet.getParameter<int32_t>("nbins"), pixelClusterPSet.getParameter<double>("xmin"), pixelClusterPSet.getParameter<double>("xmax"));
  lumi_binning_        (lumiPSet.getParameter<int32_t>        ("nbins"), lumiPSet.getParameter<double>        ("xmin"), lumiPSet.getParameter<double>        ("xmax"));
  pixellumi_binning_   (pixellumiPSet.getParameter<int32_t>   ("nbins"), pixellumiPSet.getParameter<double>   ("xmin"), pixellumiPSet.getParameter<double>   ("xmax"));
  ls_binning_          (lsPSet.getParameter<int32_t>          ("nbins"), lsPSet.getParameter<double>          ("xmin"), lsPSet.getParameter<double>          ("xmax"));
  */

  getHistoPSet  (pixellumiPSet,   pixellumi_binning_);
  getHistoPSet  (pixelClusterPSet,pixelCluster_binning_);
  getHistoPSet  (lumiPSet,        lumi_binning_);
  getHistoLSPSet(lsPSet,          ls_binning_);

  /*
  pixellumi_binning_.nbins = pixellumiPSet.getParameter<int32_t>("nbins");
  pixellumi_binning_.xmin  = pixellumiPSet.getParameter<double>("xmin");
  pixellumi_binning_.xmax  = pixellumiPSet.getParameter<double>("xmax");

  pixelCluster_binning_.nbins = pixelClusterPSet.getParameter<int32_t>("nbins");
  pixelCluster_binning_.xmin  = pixelClusterPSet.getParameter<double>("xmin");
  pixelCluster_binning_.xmax  = pixelClusterPSet.getParameter<double>("xmax");

  lumi_binning_.nbins = lumiPSet.getParameter<int32_t>("nbins");
  lumi_binning_.xmin  = lumiPSet.getParameter<double>("xmin");
  lumi_binning_.xmax  = lumiPSet.getParameter<double>("xmax");

  ls_binning_.nbins = lsPSet.getParameter<int32_t>("nbins");
  ls_binning_.xmin  = lsPSet.getParameter<double>("xmin");
  ls_binning_.xmax  = lsPSet.getParameter<double>("xmax");
  */
}

void LumiMonitor::getHistoPSet(edm::ParameterSet& pset, MEbinning& mebinning)
{
  mebinning.nbins = pset.getParameter<int32_t>("nbins");
  mebinning.xmin  = pset.getParameter<double>("xmin");
  mebinning.xmax  = pset.getParameter<double>("xmax");
}

void LumiMonitor::getHistoLSPSet(edm::ParameterSet& pset, MEbinning& mebinning)
{
  mebinning.nbins = pset.getParameter<int32_t>("nbins");
  mebinning.xmin  = 0.;
  mebinning.xmax  = double(pset.getParameter<double>("nbins"));
}

void LumiMonitor::bookHistograms(DQMStore::IBooker     & ibooker,
				 edm::Run const        & iRun,
				 edm::EventSetup const & iSetup) 
{  
  
  std::string histname, histtitle;

  std::string currentFolder = folderName_ ;
  ibooker.setCurrentFolder(currentFolder.c_str());

  if ( doPixelLumi_ ) {
    histname = "numberOfPixelClustersVsLS"; histtitle = "number of pixel clusters vs LS";
    numberOfPixelClustersVsLS_ = ibooker.book1D(histname, histtitle, 
						ls_binning_.nbins, ls_binning_.xmin, ls_binning_.xmax);
    //    numberOfPixelClustersVsLS_->getTH1()->SetCanExtend(TH1::kAllAxes);
    numberOfPixelClustersVsLS_->setAxisTitle("LS",1);
    numberOfPixelClustersVsLS_->setAxisTitle("number of pixel clusters",2);
    
    histname = "numberOfPixelClustersVsLumi"; histtitle = "number of pixel clusters vs scal lumi";
    numberOfPixelClustersVsLumi_ = ibooker.bookProfile(histname, histtitle, 
						       lumi_binning_.nbins, lumi_binning_.xmin, lumi_binning_.xmax,
						       pixelCluster_binning_.xmin,pixelCluster_binning_.xmax);
    numberOfPixelClustersVsLumi_->setAxisTitle("scal inst lumi E30 [Hz cm^{-2}]",1);
    numberOfPixelClustersVsLumi_->setAxisTitle("number of pixel clusters",2);

    histname = "pixelLumiVsLS"; histtitle = "pixel-lumi vs LS";
    pixelLumiVsLS_ = ibooker.bookProfile(histname, histtitle, 
					 ls_binning_.nbins, ls_binning_.xmin, ls_binning_.xmax,
					 pixellumi_binning_.xmin,pixellumi_binning_.xmax);
    //    pixelLumiVsLS_->getTH1()->SetCanExtend(TH1::kAllAxes);
    pixelLumiVsLS_->setAxisTitle("LS",1);
    pixelLumiVsLS_->setAxisTitle("pixel-based inst lumi E30 [Hz cm^{-2}]",2);
    
    histname = "pixelLumiVsLumi"; histtitle = "pixel-lumi vs scal lumi";
    pixelLumiVsLumi_ = ibooker.bookProfile(histname, histtitle, 
					   lumi_binning_.nbins,lumi_binning_.xmin,lumi_binning_.xmax,
					   pixellumi_binning_.xmin,lumi_binning_.xmax);
    pixelLumiVsLumi_->setAxisTitle("scal inst lumi E30 [Hz cm^{-2}]",1);
    pixelLumiVsLumi_->setAxisTitle("pixel-based inst lumi E30 [Hz cm^{-2}]",2);
  }

  histname = "lumiVsLS"; histtitle = "scal lumi vs LS";
  lumiVsLS_ = ibooker.bookProfile(histname, histtitle, 
				  ls_binning_.nbins, ls_binning_.xmin, ls_binning_.xmax,
				  lumi_binning_.xmin, lumi_binning_.xmax);
  //  lumiVsLS_->getTH1()->SetCanExtend(TH1::kAllAxes);
  lumiVsLS_->setAxisTitle("LS",1);
  lumiVsLS_->setAxisTitle("scal inst lumi E30 [Hz cm^{-2}]",2);

}

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
void LumiMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)  {

  //  int bx = iEvent.bunchCrossing();
  int ls = iEvent.id().luminosityBlock();

  float scal_lumi = -1.;
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(lumiScalersToken_, lumiScalers);
  if ( lumiScalers.isValid() && lumiScalers->size() ) {
    LumiScalersCollection::const_iterator scalit = lumiScalers->begin();
    scal_lumi = scalit->instantLumi();
  } else {
    scal_lumi = -1.;
  }
  lumiVsLS_ -> Fill(ls, scal_lumi);

  if ( doPixelLumi_ ) {
    size_t pixel_clusters = 0;
    float  pixel_lumi = -1.;
    edm::Handle< edmNew::DetSetVector<SiPixelCluster> > pixelClusters;
    iEvent.getByToken(pixelClustersToken_, pixelClusters);
    if ( pixelClusters.isValid() ) {
      
      edm::ESHandle<TrackerTopology> tTopoHandle;
      iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
      const TrackerTopology* const tTopo = tTopoHandle.product();
      
      // Count the number of clusters with at least a minimum
      // number of pixels per cluster and at least a minimum charge.
      size_t tot = 0;
      edmNew::DetSetVector<SiPixelCluster>::const_iterator  pixCluDet = pixelClusters->begin();
      for ( ; pixCluDet!=pixelClusters->end(); ++pixCluDet) {
	
	DetId detid = pixCluDet->detId();
	size_t subdetid = detid.subdetId();
	if ( subdetid == (int) PixelSubdetector::PixelBarrel ) 
	  if ( tTopo->layer(detid)==1 ) 
	    continue;
	
	edmNew::DetSet<SiPixelCluster>::const_iterator  pixClu = pixCluDet->begin();    
	for ( ; pixClu != pixCluDet->end(); ++pixClu ) {
	  ++tot;
	  if ( (pixClu->size()   >= minNumberOfPixelsPerCluster_) &&
	       (pixClu->charge() >= minPixelClusterCharge_      ) ) {
	    ++pixel_clusters;
	  }
	}
      }
      pixel_lumi = lumi_factor_per_bx_ * pixel_clusters / GetLumi::CM2_TO_NANOBARN ; // ?!?!
    } else
      pixel_lumi = -1.;

    numberOfPixelClustersVsLS_   -> Fill(ls,       pixel_clusters);
    numberOfPixelClustersVsLumi_ -> Fill(scal_lumi,pixel_clusters);
    pixelLumiVsLS_               -> Fill(ls,       pixel_lumi);
    pixelLumiVsLumi_             -> Fill(scal_lumi,pixel_lumi);
  }

}

void LumiMonitor::fillHistoPSetDescription(edm::ParameterSetDescription & pset)
{
  pset.add<int>   ( "nbins");
  pset.add<double>( "xmin" );
  pset.add<double>( "xmax" );
}

void LumiMonitor::fillHistoLSPSetDescription(edm::ParameterSetDescription & pset)
{
  pset.add<int>   ( "nbins", 2500);
}

void LumiMonitor::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>( "pixelClusters",               edm::InputTag("hltSiPixelClusters") );
  desc.add<edm::InputTag>( "scalers",                     edm::InputTag("hltScalersRawToDigi"));
  desc.add<std::string>  ( "FolderName",                  "HLT/LumiMonitoring"                );
  desc.add<bool>         ( "doPixelLumi",                 false                               );
  desc.add<bool>         ( "useBPixLayer1",               false                               );
  desc.add<int>          ( "minNumberOfPixelsPerCluster", 2                                   ); // from DQM/PixelLumi/python/PixelLumiDQM_cfi.py
  desc.add<double>       ( "minPixelClusterCharge",       15000.                              );

  edm::ParameterSetDescription histoPSet;
  edm::ParameterSetDescription pixelClusterPSet;
  LumiMonitor::fillHistoPSetDescription(pixelClusterPSet);
  histoPSet.add("pixelClusterPSet", pixelClusterPSet);

  edm::ParameterSetDescription lumiPSet;
  fillHistoPSetDescription(lumiPSet);
  histoPSet.add<edm::ParameterSetDescription>("lumiPSet", lumiPSet);

  edm::ParameterSetDescription pixellumiPSet;
  fillHistoPSetDescription(pixellumiPSet);
  histoPSet.add<edm::ParameterSetDescription>("pixellumiPSet", pixellumiPSet);

  edm::ParameterSetDescription lsPSet;
  fillHistoLSPSetDescription(lsPSet);
  histoPSet.add<edm::ParameterSetDescription>("lsPSet", lsPSet);

  desc.add<edm::ParameterSetDescription>("histoPSet",histoPSet);

  descriptions.add("lumiMonitor", desc);
}

// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LumiMonitor);
