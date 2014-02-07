#include "DQM/TrackingMonitorSource/interface/HistoHelper.h"

#include "DQMServices/ClientConfig/interface/FitSlicesYTool.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

HistoHelper::HistoHelper(const edm::ParameterSet& iConfig) 
  : pset_ (iConfig)
{
  //parameters for _vs_eta plots
  minEta     = iConfig.getParameter<double>("minEta");
  maxEta     = iConfig.getParameter<double>("maxEta");
  nintEta    = iConfig.getParameter<int>("nintEta");
  useFabsEta = iConfig.getParameter<bool>("useFabsEta");

  //parameters for _vs_pt plots
  minPt    = iConfig.getParameter<double>("minPt");
  maxPt    = iConfig.getParameter<double>("maxPt");
  nintPt   = iConfig.getParameter<int>("nintPt");
  useInvPt = iConfig.getParameter<bool>("useInvPt");
  useLogPt = iConfig.getUntrackedParameter<bool>("useLogPt",false);

  //parameters for _vs_Hit plots
  minHit  = iConfig.getParameter<double>("minHit");
  maxHit  = iConfig.getParameter<double>("maxHit");
  nintHit = iConfig.getParameter<int>("nintHit");

  //parameters for _vs_Layer plots
  minLayers  = iConfig.getParameter<double>("minLayers");
  maxLayers  = iConfig.getParameter<double>("maxLayers");
  nintLayers = iConfig.getParameter<int>("nintLayers");

  //parameters for _vs_phi plots
  minPhi  = iConfig.getParameter<double>("minPhi");
  maxPhi  = iConfig.getParameter<double>("maxPhi");
  nintPhi = iConfig.getParameter<int>("nintPhi");

  //parameters for _vs_Dxy plots
  minDxy  = iConfig.getParameter<double>("minDxy");
  maxDxy  = iConfig.getParameter<double>("maxDxy");
  nintDxy = iConfig.getParameter<int>("nintDxy");

  //parameters for _vs_Dz plots
  minDz  = iConfig.getParameter<double>("minDz");
  maxDz  = iConfig.getParameter<double>("maxDz");
  nintDz = iConfig.getParameter<int>("nintDz");

  //parameters for _vs_ProductionVertexTransvPosition plots
  minVertpos  = iConfig.getParameter<double>("minVertpos");
  maxVertpos  = iConfig.getParameter<double>("maxVertpos");
  nintVertpos = iConfig.getParameter<int>("nintVertpos");

  //parameters for _vs_ProductionVertexZPosition plots
  minZpos  = iConfig.getParameter<double>("minZpos");
  maxZpos  = iConfig.getParameter<double>("maxZpos");
  nintZpos = iConfig.getParameter<int>("nintZpos");

  //parameters for resolution plots
  ptRes_rangeMin = iConfig.getParameter<double>("ptRes_rangeMin");
  ptRes_rangeMax = iConfig.getParameter<double>("ptRes_rangeMax");
  ptRes_nbin = iConfig.getParameter<int>("ptRes_nbin");

  phiRes_rangeMin = iConfig.getParameter<double>("phiRes_rangeMin");
  phiRes_rangeMax = iConfig.getParameter<double>("phiRes_rangeMax");
  phiRes_nbin = iConfig.getParameter<int>("phiRes_nbin");

  cotThetaRes_rangeMin = iConfig.getParameter<double>("cotThetaRes_rangeMin");
  cotThetaRes_rangeMax = iConfig.getParameter<double>("cotThetaRes_rangeMax");
  cotThetaRes_nbin = iConfig.getParameter<int>("cotThetaRes_nbin");

  dxyRes_rangeMin = iConfig.getParameter<double>("dxyRes_rangeMin");
  dxyRes_rangeMax = iConfig.getParameter<double>("dxyRes_rangeMax");
  dxyRes_nbin = iConfig.getParameter<int>("dxyRes_nbin");

  dzRes_rangeMin = iConfig.getParameter<double>("dzRes_rangeMin");
  dzRes_rangeMax = iConfig.getParameter<double>("dzRes_rangeMax");
  dzRes_nbin = iConfig.getParameter<int>("dzRes_nbin");

  // fix for the LogScale
  if(useLogPt){
    maxPt=log10(maxPt);
    if(minPt > 0){
      minPt=log10(minPt);
    }
    else{
    }
  }

}

void HistoHelper::setUpVectors(){

  double step=(maxEta-minEta)/nintEta;
  //std::ostringstream title,name; ///BM, what is this?
  ETAintervals.push_back(minEta);
  for (int k=1;k<nintEta+1;k++) {
    double d=minEta+k*step;
    ETAintervals.push_back(d);
  }

  double stepPt = (maxPt-minPt)/nintPt;
  PTintervals.push_back(minPt);
  for (int k=1;k<nintPt+1;k++) {
    double d=0;
    if(useLogPt)d=pow(10,minPt+k*stepPt);
    else d=minPt+k*stepPt;
    PTintervals.push_back(d);
  }

  //  for (int k=1;k<nintHit+1;k++) {
  //    totSIMv_hit.push_back(0);
  //  }

  double stepPhi = (maxPhi-minPhi)/nintPhi;
  PHIintervals.push_back(minPhi);
  for (int k=1;k<nintPhi+1;k++) {
    double d=minPhi+k*stepPhi;
    PHIintervals.push_back(d);
  }

  double stepDxy = (maxDxy-minDxy)/nintDxy;
  DXYintervals.push_back(minDxy);
  for (int k=1;k<nintDxy+1;k++) {
    double d=minDxy+k*stepDxy;
    DXYintervals.push_back(d);
  }

  double stepDz = (maxDz-minDz)/nintDz;
  DZintervals.push_back(minDz);
  for (int k=1;k<nintDz+1;k++) {
    double d=minDz+k*stepDz;
    DZintervals.push_back(d);
  }

  double stepVertpos = (maxVertpos-minVertpos)/nintVertpos;
  VERTEXPOSintervals.push_back(minVertpos);
  for (int k=1;k<nintVertpos+1;k++) {
    double d=minVertpos+k*stepVertpos;
    VERTEXPOSintervals.push_back(d);
  }

  double stepZpos = (maxZpos-minZpos)/nintZpos;
  ZPOSintervals.push_back(minZpos);
  for (int k=1;k<nintZpos+1;k++) {
    double d=minZpos+k*stepZpos;
    ZPOSintervals.push_back(d);
  }

}

void HistoHelper::doProfileX(TH2 * th2, MonitorElement* me){
  if (th2->GetNbinsX()==me->getNbinsX()){
    TProfile * p1 = (TProfile*) th2->ProfileX();
    p1->Copy(*me->getTProfile());
    delete p1;
  } else {
    throw cms::Exception("MultiTrackValidator") << "Different number of bins!";
  }
}


void HistoHelper::fillPlotFromVector(MonitorElement* h, std::vector<int>& vec) {
  int i = 0;
  for (auto j : vec ) {
    i++;
    h->setBinContent(i, j);
  }
}

void HistoHelper::fillPlotFromVectors(MonitorElement* h, 
				      std::vector<int>& numerator, 
				      std::vector<int>& denominator,
				      std::string type){
  double value = 0.;
  double err = 0.;
  for ( size_t j=0, n=numerator.size(); j<n-1; j++ ) {
    if (denominator[j]!=0){

      if (type=="effic"){
	value = ((double) numerator[j]) / ((double) denominator[j]);
        err = sqrt( value*(1.-value)/(double) denominator[j] );

      } else if (type=="fakerate"){
	value = 1.-((double) numerator[j]) / ((double) denominator[j]);
        err = sqrt( value*(1.-value)/(double) denominator[j] );

      } else if (type=="pileup"){
	value = ((double) numerator[j]) / ((double) denominator[j]);
        err = sqrt( value*(1.+value)/(double) denominator[j] );

      } else return;

      h->setBinContent(j+1, value);
      h->setBinError(j+1, err);
    }
    else {
      h->setBinContent(j+1, 0.);
      h->setBinError(j+1, 0.);
    }
  }

}



#include <TMath.h>

void HistoHelper::BinLogX(TH1*h){  
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();
  
  float from = axis->GetXmin();
  float to = axis->GetXmax();
  float width = (to - from) / bins;
  float *new_bins = new float[bins + 1];
  
  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
    
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}


