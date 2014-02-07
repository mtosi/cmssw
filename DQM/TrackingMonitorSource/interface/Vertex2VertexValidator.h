//
// Original Author:  Mia Tosi
//         Created:  Sat, 01 Feb 2014 15:47:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

#include <DQMServices/Core/interface/DQMEDAnalyzer.h>

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include <iostream>
#include <sstream>
#include <string>

//
// class declaration
//
class DQMStore;
namespace reco {
  class Track;
  class BeamSpot;
  class Vertex;
}
class DQMStore;

class Vertex2VertexValidator : public DQMEDAnalyzer {

  /*
  vector<pair<int, map<double, int> > > res
  for j in loop_off_trk
     map<double, int> tmp
     for i in loop_online_trk
       tmp[dR(j,i)] = i
     res.push_back(make_pair<j, tmp>)
  */     

   public:

  MonitorElement* h_nvertices_vs_nvertices;
  
  struct generalME {
    std::string label;
    MonitorElement *h_ntracks, *h_sumpt2;
    MonitorElement *h_dzmin;
  };
  
  struct assME {
    std::string label;
    MonitorElement *h_ntracks_vs_ntracks, *h_sumpt2_vs_sumpt2;
  };
  
  typedef std::vector<std::pair<int, std::map<double, int> > > idx2idxByDoubleColl;

      explicit Vertex2VertexValidator(const edm::ParameterSet&);
      ~Vertex2VertexValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      void beginJob(const edm::EventSetup& iSetup);
      void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void bookHistograms(DQMStore::IBooker & iBooker, edm::Run const & iRun, edm::EventSetup const & iSetup) override;
      void endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup);
      void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
      // ----------member data ---------------------------
 protected:
  
  void fillMap(reco::VertexCollection tracks1, reco::VertexCollection tracks2, idx2idxByDoubleColl& map);

  void initialize_parameter(const edm::ParameterSet& iConfig);
  void bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);
  void book_generic_vertex_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);
  void book_associate_vertex_histos(DQMStore::IBooker & ibooker, assME& mes, TString label, std::string & dir);

  void fill_generic_vertex_histos(generalME& mes, reco::Vertex* trk);
  void fill_associate_vertex_histos(assME& mes, reco::Vertex* mon, reco::Vertex* ref);


  DQMStore* dqmStore_;

  edm::InputTag numVertexInputTag_;
  edm::InputTag denVertexInputTag_;

  //these are used by MTVGenPs
  edm::EDGetTokenT<reco::VertexCollection> numpvToken_;
  edm::EDGetTokenT<reco::VertexCollection> denpvToken_;

  std::string out;

 private:
  //  edm::ParameterSet conf_;
  std::string topDirName_;
  double dzmin_;

  generalME numVertexMEs_;  
  generalME numassVertexMEs_;  
  generalME denVertexMEs_;  
  generalME denassVertexMEs_;  

  assME assVertexMEs_;

  double minNtracks, maxNtracks; int ntracks;
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


  MonitorElement* nrec_vs_nsim;
  MonitorElement* nrecHit_vs_nsimHit_sim2rec;
  MonitorElement* nrecHit_vs_nsimHit_rec2sim;

  //assoc hits
  MonitorElement* h_assocFraction, h_assocSharedHit;

  //pulls of track params vs eta: to be used with fitslicesytool
  MonitorElement* pull_dxy_vs_eta, ptpull_vs_eta, dzpull_vs_eta, phipull_vs_eta, thetapull_vs_eta;

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

