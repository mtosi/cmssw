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

class Track2TrackValidator : public DQMEDAnalyzer {

  /*
  vector<pair<int, map<double, int> > > res
  for j in loop_off_trk
     map<double, int> tmp
     for i in loop_online_trk
       tmp[dR(j,i)] = i
     res.push_back(make_pair<j, tmp>)
  */     

   public:

  struct generalME {
    std::string label;
    MonitorElement *h_tracks, *h_pt, *h_eta, *h_phi, *h_dxy, *h_dz, *h_dxyWRTpv, *h_dzWRTpv, *h_charge, *h_hits;
    MonitorElement *h_dRmin;
    MonitorElement *h_pt_vs_eta;
  };
  
  typedef std::vector<std::pair<int, std::map<double, int> > > idx2idxByDoubleColl;

      explicit Track2TrackValidator(const edm::ParameterSet&);
      ~Track2TrackValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      void beginJob(const edm::EventSetup& iSetup);
      void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void bookHistograms(DQMStore::IBooker & iBooker, edm::Run const & iRun, edm::EventSetup const & iSetup) override;
      void endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup);
      void endRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
      // ----------member data ---------------------------
 protected:
  
  void fillMap(reco::TrackCollection tracks1, reco::TrackCollection tracks2, idx2idxByDoubleColl& map);

  void initialize_parameter(const edm::ParameterSet& iConfig);
  void bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);
  void book_generic_tracks_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir);

  void fill_generic_tracks_histos(generalME& mes, reco::Track* trk, reco::BeamSpot* bs, reco::Vertex* pv);


  DQMStore* dqmStore_;

  edm::InputTag numTrackInputTag_;
  edm::InputTag denTrackInputTag_;

  //these are used by MTVGenPs
  edm::EDGetTokenT<reco::TrackCollection> numTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> denTrackToken_;
  edm::EDGetTokenT<reco::BeamSpot> numbsToken_;
  edm::EDGetTokenT<reco::BeamSpot> denbsToken_;
  edm::EDGetTokenT<reco::VertexCollection> numpvToken_;
  edm::EDGetTokenT<reco::VertexCollection> denpvToken_;

  std::string out;

 private:
  //  edm::ParameterSet conf_;
  std::string topDirName_;
  double dRmin_;

  generalME numTracksMEs_;  
  generalME numassTracksMEs_;  
  generalME denTracksMEs_;  
  generalME denassTracksMEs_;  

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

#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Track2TrackValidator::Track2TrackValidator(const edm::ParameterSet& iConfig)
  : numTrackInputTag_ ( iConfig.getParameter<edm::InputTag>("numTrack")       )
  , denTrackInputTag_ ( iConfig.getParameter<edm::InputTag>("denTrack")       )
  , topDirName_       ( iConfig.getParameter<std::string>  ("topDirName")     )
  , dRmin_            ( iConfig.getParameter<double>("dRmin")                 )
{
   //now do what ever initialization is needed
  numTrackToken_      = consumes<reco::TrackCollection>(numTrackInputTag_);
  denTrackToken_      = consumes<reco::TrackCollection>(denTrackInputTag_);
  numbsToken_         = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("numBeamSpot"));
  denbsToken_         = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("denBeamSpot"));
  numpvToken_         = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("numPrimaryVertices"));
  denpvToken_         = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("denPrimaryVertices"));
  
  numTracksMEs_.label    = numTrackInputTag_.label();
  numassTracksMEs_.label = numTrackInputTag_.label()+"_ass";
  denTracksMEs_.label    = denTrackInputTag_.label();
  denassTracksMEs_.label = denTrackInputTag_.label()+"_ass";
}


Track2TrackValidator::~Track2TrackValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void
Track2TrackValidator::beginJob(const edm::EventSetup& iSetup) {
}

//
// member functions
//
// ------------ method called for each event  ------------
void
Track2TrackValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle<reco::TrackCollection> numTracksHandle;
  iEvent.getByToken(numTrackToken_,numTracksHandle);
  reco::TrackCollection numTracks = *numTracksHandle;

  edm::Handle<reco::TrackCollection> denTracksHandle;
  iEvent.getByToken(denTrackToken_,denTracksHandle);
  reco::TrackCollection denTracks = *denTracksHandle;

  edm::LogVerbatim("Track2TrackValidator") << "analyzing "
					   << numTrackInputTag_.process()  << ":"
					   << numTrackInputTag_.label()    << ":"
					   << numTrackInputTag_.instance() << " w.r.t. "
					   << denTrackInputTag_.process()  << ":"
					   << denTrackInputTag_.label()    << ":"
					   << denTrackInputTag_.instance() << " \n";
  
  std::cout << "analyzing "
	    << numTrackInputTag_.process()  << ":"
	    << numTrackInputTag_.label()    << ":"
	    << numTrackInputTag_.instance() << " w.r.t. "
	    << denTrackInputTag_.process()  << ":"
	    << denTrackInputTag_.label()    << ":"
	    << denTrackInputTag_.instance() << " \n";
  
  idx2idxByDoubleColl num2denColl;
  fillMap(numTracks,denTracks,num2denColl);

  idx2idxByDoubleColl den2numColl;
  fillMap(denTracks,numTracks,den2numColl);

  // ########################################################
  // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
  // ########################################################
  
  //compute number of tracks per eta interval
  //
  edm::LogVerbatim("Track2TrackValidator") << "\n# of tracks (denominator): " << denTracks.size() << "\n";
  int st(0);  	     //This counter counts the number of tracks
  int ast(0);  	     //This counter counts the number of tracks that are "associated" to recoTracks
  unsigned asts(0);  //This counter counts the number of tracks that are "associated" to recoTracks surviving the bunchcrossing cut

  edm::Handle<reco::BeamSpot> denbsHandle;
  iEvent.getByToken(denbsToken_,denbsHandle);
  reco::BeamSpot denbs = *denbsHandle;
  
  edm::Handle<reco::VertexCollection> denpvHandle;
  iEvent.getByToken(denpvToken_,denpvHandle);
  reco::Vertex denpv = denpvHandle->at(0);
  // loop over tracks (denominator)
  std::cout << "[Track2TrackValidator::analyze] starting loop over denominator" << std::endl;
  for (idx2idxByDoubleColl::const_iterator pItr = den2numColl.begin(), eItr = den2numColl.end();
       pItr != eItr; ++pItr) {
    
    int trackIdx = pItr->first;
    reco::Track track = denTracks.at(trackIdx);

    fill_generic_tracks_histos(*&denTracksMEs_,&track,&denbs,&denpv);
    st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
    
    std::map<double,int> trackDRmap = pItr->second;
    double dRmin = trackDRmap.begin()->first;
    (denTracksMEs_.h_dRmin)->Fill(dRmin);
    bool matched = false;
    if ( dRmin < dRmin_ ) matched = true;
    if ( matched ) {
      fill_generic_tracks_histos(*&denassTracksMEs_,&track,&denbs,&denpv);
      (denassTracksMEs_.h_dRmin)->Fill(dRmin);
    }
    for (std::map<double,int>::const_iterator mItr = trackDRmap.begin(), meItr = trackDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dR: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
  }

  edm::Handle<reco::BeamSpot> numbsHandle;
  iEvent.getByToken(numbsToken_,numbsHandle);
  reco::BeamSpot numbs = *numbsHandle;

  edm::Handle<reco::VertexCollection> numpvHandle;
  iEvent.getByToken(numpvToken_,numpvHandle);
  reco::Vertex numpv = numpvHandle->at(0);
  // loop over tracks (numerator)
  for (idx2idxByDoubleColl::const_iterator pItr = num2denColl.begin(), eItr = num2denColl.end();
       pItr != eItr; ++pItr) {

    int trackIdx = pItr->first;
    reco::Track track = numTracks.at(trackIdx);
    fill_generic_tracks_histos(*&numTracksMEs_,&track,&numbs,&numpv);
   
    st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
    
    int matchedTrackIdx = -1;
    std::map<double,int> trackDRmap = pItr->second;
    double dRmin = trackDRmap.begin()->first;
    (numTracksMEs_.h_dRmin)->Fill(dRmin);
    std::cout << "num2denColl[" <<  trackIdx << "]" << " map size: " << trackDRmap.size() << " dRmin: " << dRmin << std::endl;
    bool matched = false;
    if ( dRmin < dRmin_ ) matched = true;
    if ( !matched ) {
      (numassTracksMEs_.h_dRmin)->Fill(dRmin);
      fill_generic_tracks_histos(*&numassTracksMEs_,&track,&numbs,&numpv);
    }
    for (std::map<double,int>::const_iterator mItr = trackDRmap.begin(), meItr = trackDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dR: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
    //    fill_recoAssociated_simTrack_histos(w,*tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,matchedTrackPointer,puinfo.getPU_NumInteractions(), vtx_z_PU);
    if (matchedTrackIdx >= 0) asts++;

  } // end loop over tracks (denominator)

  // ##############################################
  // fill reco Tracks histograms (LOOP OVER TRACKS)
  // ##############################################
  edm::LogVerbatim("Track2TrackValidator") << "\n# of reco::Tracks with "
					   << numTrackInputTag_.process()<<":"
					   << numTrackInputTag_.label()<<":"
					   << numTrackInputTag_.instance()
					   << ": " << numTracks.size() << "\n";
  
  int at(0);  //This counter counts the number of recoTracks that are associated to Tracks
  int rT(0);  //This counter counts the number of recoTracks in general
  
  edm::LogVerbatim("Track2TrackValidator") << "Total tracks: "        << st << "\n"
					   << "Total Associated: "    << ast << "\n"
					   << "Total Reconstructed: " << rT << "\n"
					   << "Total Associated: "    << at << "\n"
					   << "Total Fakes: "         << rT-at << "\n";
  
}

void 
Track2TrackValidator::bookHistograms(DQMStore::IBooker & ibooker,
				     edm::Run const & iRun,
				     edm::EventSetup const & iSetup)
{
  
  std::cout << "[Track2TrackValidator::bookHistograms]" << std::endl;

  std::string dir = topDirName_;

  bookHistos(ibooker,numTracksMEs_,   "mon",    dir);
  bookHistos(ibooker,numassTracksMEs_,"mon_ass",dir);
  bookHistos(ibooker,denTracksMEs_,   "ref",    dir);
  bookHistos(ibooker,denassTracksMEs_,"ref_ass",dir);

}

void 
Track2TrackValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
Track2TrackValidator::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  //  fillHistosFromVectors();
  if ( out.size() != 0 && dqmStore_ ) dqmStore_->save(out);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Track2TrackValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
Track2TrackValidator::fillMap(reco::TrackCollection tracks1, reco::TrackCollection tracks2, idx2idxByDoubleColl& map)
{
  int i = 0;
  for ( auto track1 : tracks1 ) {
    std::map<double,int> tmp;
    int j = 0;
    for ( auto track2 : tracks2 ) {
      double dR = reco::deltaR(track1.eta(),track1.phi(),track2.eta(),track2.phi());
      tmp[dR] = j;
      j++;
    }
    map.push_back(std::make_pair(i,tmp));
    i++;
  }
  std::cout << "map: " << map.size() << "[tracks1: " << tracks1.size() << ", tracks2: " << tracks2.size() << "]" << std::endl;

}

void 
Track2TrackValidator::bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

  book_generic_tracks_histos(ibooker,mes,label,dir);

}

void 
Track2TrackValidator::book_generic_tracks_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

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
Track2TrackValidator::fill_generic_tracks_histos(generalME& mes, reco::Track* trk, reco::BeamSpot* bs, reco::Vertex* pv) {

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

void
Track2TrackValidator::initialize_parameter(const edm::ParameterSet& iConfig)
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

//define this as a plug-in
DEFINE_FWK_MODULE(Track2TrackValidator);
