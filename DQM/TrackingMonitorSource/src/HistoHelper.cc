#include "DQM/TrackingMonitorSource/interface/HistoHelper.h"

#include "DQMServices/ClientConfig/interface/FitSlicesYTool.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

HistoHelper::HistoHelper(const edm::ParameterSet& iConfig) 
  : pset_ (iConfig)
{
  //parameters for _vs_eta plots
  minEta     = pset_.getParameter<double>("minEta");
  maxEta     = pset_.getParameter<double>("maxEta");
  nintEta    = pset_.getParameter<int>("nintEta");
  useFabsEta = pset_.getParameter<bool>("useFabsEta");

  //parameters for _vs_pt plots
  minPt    = pset_.getParameter<double>("minPt");
  maxPt    = pset_.getParameter<double>("maxPt");
  nintPt   = pset_.getParameter<int>("nintPt");
  useInvPt = pset_.getParameter<bool>("useInvPt");
  useLogPt = pset_.getUntrackedParameter<bool>("useLogPt",false);

  //parameters for _vs_Hit plots
  minHit  = pset_.getParameter<double>("minHit");
  maxHit  = pset_.getParameter<double>("maxHit");
  nintHit = pset_.getParameter<int>("nintHit");

  //parameters for _vs_Layer plots
  minLayers  = pset_.getParameter<double>("minLayers");
  maxLayers  = pset_.getParameter<double>("maxLayers");
  nintLayers = pset_.getParameter<int>("nintLayers");

  //parameters for _vs_phi plots
  minPhi  = pset_.getParameter<double>("minPhi");
  maxPhi  = pset_.getParameter<double>("maxPhi");
  nintPhi = pset_.getParameter<int>("nintPhi");

  //parameters for _vs_Dxy plots
  minDxy  = pset_.getParameter<double>("minDxy");
  maxDxy  = pset_.getParameter<double>("maxDxy");
  nintDxy = pset_.getParameter<int>("nintDxy");

  //parameters for _vs_Dz plots
  minDz  = pset_.getParameter<double>("minDz");
  maxDz  = pset_.getParameter<double>("maxDz");
  nintDz = pset_.getParameter<int>("nintDz");

  //parameters for _vs_ProductionVertexTransvPosition plots
  minVertpos  = pset_.getParameter<double>("minVertpos");
  maxVertpos  = pset_.getParameter<double>("maxVertpos");
  nintVertpos = pset_.getParameter<int>("nintVertpos");

  //parameters for _vs_ProductionVertexZPosition plots
  minZpos  = pset_.getParameter<double>("minZpos");
  maxZpos  = pset_.getParameter<double>("maxZpos");
  nintZpos = pset_.getParameter<int>("nintZpos");

  //parameters for resolution plots
  ptRes_rangeMin = pset_.getParameter<double>("ptRes_rangeMin");
  ptRes_rangeMax = pset_.getParameter<double>("ptRes_rangeMax");
  ptRes_nbin = pset_.getParameter<int>("ptRes_nbin");

  phiRes_rangeMin = pset_.getParameter<double>("phiRes_rangeMin");
  phiRes_rangeMax = pset_.getParameter<double>("phiRes_rangeMax");
  phiRes_nbin = pset_.getParameter<int>("phiRes_nbin");

  cotThetaRes_rangeMin = pset_.getParameter<double>("cotThetaRes_rangeMin");
  cotThetaRes_rangeMax = pset_.getParameter<double>("cotThetaRes_rangeMax");
  cotThetaRes_nbin = pset_.getParameter<int>("cotThetaRes_nbin");

  dxyRes_rangeMin = pset_.getParameter<double>("dxyRes_rangeMin");
  dxyRes_rangeMax = pset_.getParameter<double>("dxyRes_rangeMax");
  dxyRes_nbin = pset_.getParameter<int>("dxyRes_nbin");

  dzRes_rangeMin = pset_.getParameter<double>("dzRes_rangeMin");
  dzRes_rangeMax = pset_.getParameter<double>("dzRes_rangeMax");
  dzRes_nbin = pset_.getParameter<int>("dzRes_nbin");

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
void 
HistoHelper::bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

  book_generic_tracks_histos(ibooker,mes,label,dir);

}

void 
HistoHelper::book_generic_tracks_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

  ibooker.cd();
  ibooker.setCurrentFolder(dir);

  (mes.h_pt )      = ibooker.book1D(label+"_pt",       "track p_{T}",                               nintPt,  minPt,  maxPt   );
  (mes.h_eta)      = ibooker.book1D(label+"_eta",      "track pseudorapidity",                     nintEta, minEta, maxEta   );
  (mes.h_phi)      = ibooker.book1D(label+"_phi",      "track #phi",                               nintPhi, minPhi, maxPhi   );
  (mes.h_dxy)      = ibooker.book1D(label+"_dxy",      "track transverse dca to beam spot",        nintDxy, minDxy, maxDxy   );
  (mes.h_dz )      = ibooker.book1D(label+"_dz",       "track longitudinal dca to beam spot",       nintDz,  minDz,  maxDz   );
  (mes.h_dxyWRTpv) = ibooker.book1D(label+"_dxyWRTpv", "track transverse dca to primary vertex",   nintDxy, minDxy, maxDxy   );
  (mes.h_dzWRTpv)  = ibooker.book1D(label+"_dzWRTpv",  "track longitudinal dca to primary vertex",  nintDz,  minDz,  maxDz   );
  (mes.h_charge)   = ibooker.book1D(label+"_charge",   "track charge",                                   5,   -2,        2   );
  (mes.h_hits  )   = ibooker.book1D(label+"_hits",     "track number of hits",                          35,   -0.5,     34.5 );
  (mes.h_dRmin)    = ibooker.book1D(label+"_dRmin",    "track min dR",                                 100,    0.,       0.01); 

  (mes.h_pt_vs_eta)  = ibooker.book2D(label+"_ptVSeta","track p_{T} vs #eta", nintEta, minEta, maxEta, nintPt, minPt, maxPt);

}


void
HistoHelper::fill_generic_tracks_histos(generalME& mes, reco::Track* trk, reco::BeamSpot* bs, reco::Vertex* pv) {

  float pt       = trk->pt();
  float eta      = trk->eta();
  float phi      = trk->phi();
  float dxy      = trk->dxy(bs->position());
  float dz       = trk->dz(bs->position());
  float dxyWRTpv = trk->dxy(pv->position());
  float dzWRTpv  = trk->dz(pv->position());
  float charge   = trk->charge();
  float nhits    = trk->hitPattern().numberOfValidHits();

  (mes.h_pt      ) -> Fill(pt);
  (mes.h_eta     ) -> Fill(eta);
  (mes.h_phi     ) -> Fill(phi);
  (mes.h_dxy     ) -> Fill(dxy);
  (mes.h_dz      ) -> Fill(dz);
  (mes.h_dxyWRTpv) -> Fill(dxyWRTpv);
  (mes.h_dzWRTpv ) -> Fill(dzWRTpv);
  (mes.h_charge  ) -> Fill(charge);
  (mes.h_hits    ) -> Fill(nhits);

  (mes.h_pt_vs_eta) -> Fill(eta,pt);
  
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


