#ifndef HistoHelper_h
#define HistoHelper_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"


class HistoHelper {
 public:

  HistoHelper(const edm::ParameterSet& iConfig);
  virtual ~HistoHelper() {}

  // to be implemented in the concrete classes
  void initialize();
  void setUpVectors();

  struct generalME {
    std::string label;
    MonitorElement *h_tracks, *h_pt, *h_eta, *h_phi, *h_dxy, *h_dz, *h_charge, *h_hits;
    MonitorElement *h_dRmin;
    MonitorElement *h_pt_vs_eta;
  };
  
  void bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);
  void book_generic_tracks_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);

  void fill_generic_tracks_histos(generalME& mes, reco::Track* trk, reco::BeamSpot* bs);
  void fillNumeratorHistos(const reco::Track::Vector&,const reco::Track::Point& vertex, int bx);

  void fillDenominatorHistos(
			     const reco::Track& tp,
			     const reco::Track::Vector& momentumTP, const reco::Track::Point& vertexTP,
			     double dxy, double dz, int nSimHits,
			     const reco::Track* track,
			     int numVertices, double vertz);

  void fillAssociatedNumeratorHistos(
					     const reco::Track& tp,
					     const reco::Track::Vector & momentumTP, const reco::Track::Point & vertexTP,
					     double dxy, double dz, int nSimHits,
					     const reco::Track* track,
					     int numVertices, double vertz);

  void fill_generic_recoTrack_histos(
				     	     const reco::Track& track,
				     	     const math::XYZPoint& bsPosition,
				     	     bool isMatched,
				     	     bool isSigMatched,
				     	     bool isChargeMatched,
					     int numAssocRecoTracks,
                         	             int numVertices,
                         		     int tpbunchcrossing,
				             int nSimHits,
   					     double sharedFraction);

  void fillAssociatedDenominatorHistos(
					       const reco::Track& track);

  void fillGeneralHistos(
				 int assTracks,
				 int numRecoTracks,
				 int numSimTracks);

  void fillResolutionHistos(
				    const reco::Track::Vector& momentumTP,
				    const reco::Track::Point& vertexTP,
				    int chargeTP,
				    const reco::Track& track,
				    const math::XYZPoint& bsPosition);

  void fillPullHistos(
			      const reco::Track::Vector& momentumTP,
			      const reco::Track::Point& vertexTP,
			      int chargeTP,
			      const reco::Track& track,
			      const math::XYZPoint& bsPosition);

  void finalHistoFits();


  void fillHistosFromVectors();
  void fillProfileHistosFromVectors();


 protected:
  //protected functions
  
  
  double getEta(double eta);

  double getPt(double pt);

  void doProfileX(TH2 * th2, MonitorElement* me);

  void doProfileX(MonitorElement * th2m, MonitorElement* me) {
    doProfileX(th2m->getTH2F(), me);
  }

  void fillPlotFromVector(MonitorElement* h, std::vector<int>& vec);

  void fillPlotFromVectors(MonitorElement* h,
			   std::vector<int>& numerator,
			   std::vector<int>& denominator,
			   std::string type);

  void BinLogX(TH1*h);

 private:
  //private data members
  const edm::ParameterSet& pset_;


  double minEta, maxEta;  int nintEta;  bool useFabsEta;
  double minPt, maxPt;  int nintPt;   bool useInvPt;   bool useLogPt;
  double minHit, maxHit;  int nintHit;
  double minLayers, maxLayers;  int nintLayers;
  double minPhi, maxPhi;  int nintPhi;
  double minDxy, maxDxy;  int nintDxy;
  double minDz, maxDz;  int nintDz;
  double minVertpos, maxVertpos;  int nintVertpos;
  double minZpos, maxZpos;  int nintZpos;
  double minDeDx, maxDeDx;  int nintDeDx;
  double minVertcount, maxVertcount;  int nintVertcount;

  //
  double ptRes_rangeMin,ptRes_rangeMax; int ptRes_nbin;
  double phiRes_rangeMin,phiRes_rangeMax; int phiRes_nbin;
  double cotThetaRes_rangeMin,cotThetaRes_rangeMax; int cotThetaRes_nbin;
  double dxyRes_rangeMin,dxyRes_rangeMax; int dxyRes_nbin;
  double dzRes_rangeMin,dzRes_rangeMax; int dzRes_nbin;


  // denominator

  //1D
  MonitorElement* h_eff_vs_eta, h_fake_vs_eta;
  MonitorElement* h_eff_vs_pt,  h_fake_vs_pt;
  MonitorElement* h_eff_vs_hit, h_fake_vs_hit;
  MonitorElement* h_eff_vs_phi, h_fake_vs_phi;
  MonitorElement* h_eff_vs_dxy, h_fake_vs_dxy;
  MonitorElement* h_eff_vs_dz,  h_fake_vs_dz;

  
  MonitorElement* nrec_vs_nsim;
  MonitorElement* nrecHit_vs_nsimHit_sim2rec;
  MonitorElement* nrecHit_vs_nsimHit_rec2sim;

  //assoc hits
  MonitorElement* h_assocFraction, h_assocSharedHit;

  //pulls of track params vs eta: to be used with fitslicesytool
  MonitorElement* pull_dxy_vs_eta, ptpull_vs_eta, dzpull_vs_eta, phipull_vs_eta, thetapull_vs_eta;

  std::vector<double> ETAintervals;
  std::vector<double> PTintervals;
  std::vector<double> PHIintervals;
  std::vector<double> DXYintervals;
  std::vector<double> DZintervals;
  std::vector<double> VERTEXPOSintervals;
  std::vector<double> ZPOSintervals;
  std::vector<int> totSIMeta,totRECeta,totASSeta,totASS2eta,totloopeta,totmisideta,totASS2etaSig;
  std::vector<int> totSIMpT,totRECpT,totASSpT,totASS2pT,totlooppT,totmisidpT;
  std::vector<int> totSIM_hit,totREC_hit,totASS_hit,totASS2_hit,totloop_hit,totmisid_hit;
  std::vector<int> totSIM_phi,totREC_phi,totASS_phi,totASS2_phi,totloop_phi,totmisid_phi;
  std::vector<int> totSIM_dxy,totREC_dxy,totASS_dxy,totASS2_dxy,totloop_dxy,totmisid_dxy;
  std::vector<int> totSIM_dz,totREC_dz,totASS_dz,totASS2_dz,totloop_dz,totmisid_dz;

  std::vector<int> totSIM_vertpos,totASS_vertpos,totSIM_zpos,totASS_zpos;
  std::vector<int> totSIM_vertcount_entire,totASS_vertcount_entire,totREC_vertcount_entire,totASS2_vertcount_entire,totASS2_vertcount_entire_signal;
  std::vector<int> totSIM_vertcount_barrel,totASS_vertcount_barrel,totREC_vertcount_barrel,totASS2_vertcount_barrel;
  std::vector<int> totSIM_vertcount_fwdpos,totASS_vertcount_fwdpos,totREC_vertcount_fwdpos,totASS2_vertcount_fwdpos;
  std::vector<int> totSIM_vertcount_fwdneg,totASS_vertcount_fwdneg,totREC_vertcount_fwdneg,totASS2_vertcount_fwdneg;
  std::vector<int> totSIM_vertz_entire,totASS_vertz_entire;
  std::vector<int> totSIM_vertz_barrel,totASS_vertz_barrel;
  std::vector<int> totSIM_vertz_fwdpos,totASS_vertz_fwdpos;
  std::vector<int> totSIM_vertz_fwdneg,totASS_vertz_fwdneg;
  std::vector<int> totREC_algo;
  std::vector<int> totREC_ootpu_entire, totASS2_ootpu_entire;
  std::vector<int> totREC_ootpu_barrel, totASS2_ootpu_barrel;
  std::vector<int> totREC_ootpu_fwdpos, totASS2_ootpu_fwdpos;
  std::vector<int> totREC_ootpu_fwdneg, totASS2_ootpu_fwdneg;
  std::vector<int> totREC_ootpu_eta_entire, totASS2_ootpu_eta_entire;
  std::vector<int> totASS2_itpu_eta_entire, totASS2_itpu_eta_entire_signal, totASS2_itpu_vertcount_entire, totASS2_itpu_vertcount_entire_signal;
  std::vector<int> totFOMT_eta, totFOMT_vertcount;
  std::vector<int> totCONeta, totCONvertcount, totCONzpos;

};

#endif
