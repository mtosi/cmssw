#include "DQM/TrackingMonitorSource/interface/Track2TrackValidator.h"

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
  initialize_parameter(iConfig);

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

  assTracksMEs_.label = "associate";
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
  
  /*
  std::cout << "analyzing "
	    << numTrackInputTag_.process()  << ":"
	    << numTrackInputTag_.label()    << ":"
	    << numTrackInputTag_.instance() << " w.r.t. "
	    << denTrackInputTag_.process()  << ":"
	    << denTrackInputTag_.label()    << ":"
	    << denTrackInputTag_.instance() << " \n";
  */
  
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
  //  std::cout << "[Track2TrackValidator::analyze] starting loop over denominator" << std::endl;
  for (idx2idxByDoubleColl::const_iterator pItr = den2numColl.begin(), eItr = den2numColl.end();
       pItr != eItr; ++pItr) {
    
    int trackIdx = pItr->first;
    reco::Track track = denTracks.at(trackIdx);

    fill_generic_tracks_histos(*&denTracksMEs_,&track,&denbs,&denpv);
    st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
    
    std::map<double,int> trackDRmap = pItr->second;
    if (trackDRmap.size() == 0) {
      (denassTracksMEs_.h_dRmin)->Fill(-1.);
      continue;
    }
    double dRmin = trackDRmap.begin()->first;
    //    std::cout << "[Track2TrackValidator::analyze] den2numColl loop: dRmin: " << dRmin << " [map size: " << trackDRmap.size() << "]" << std::endl;
    (denTracksMEs_.h_dRmin)->Fill(dRmin);
    bool matched = false;
    if ( dRmin < dRmin_ ) matched = true;
    if ( matched ) {
      fill_generic_tracks_histos(*&denassTracksMEs_,&track,&denbs,&denpv);
      (denassTracksMEs_.h_dRmin)->Fill(dRmin);
      int matchedTrackIndex = trackDRmap[dRmin];
      //      std::cout << " --> matchedTrackIndex: " << matchedTrackIndex << std::endl;
      reco::Track matchedTrack = numTracks.at(matchedTrackIndex);
      fill_associate_tracks_histos(*&assTracksMEs_,&track,&matchedTrack,&denbs,&denpv);
    }
    for (std::map<double,int>::const_iterator mItr = trackDRmap.begin(), meItr = trackDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      //      std::cout << " --> dR: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
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
    if (trackDRmap.size() == 0) {
      (numassTracksMEs_.h_dRmin)->Fill(-1.);
      continue;
    }
    double dRmin = trackDRmap.begin()->first;
    (numTracksMEs_.h_dRmin)->Fill(dRmin);
    //    //    std::cout << "num2denColl[" <<  trackIdx << "]" << " map size: " << trackDRmap.size() << " dRmin: " << dRmin << std::endl;
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
  
  //  std::cout << "[Track2TrackValidator::bookHistograms]" << std::endl;

  std::string dir = topDirName_;

  bookHistos(ibooker,numTracksMEs_,   "mon",    dir);
  bookHistos(ibooker,numassTracksMEs_,"mon_ass",dir);
  bookHistos(ibooker,denTracksMEs_,   "ref",    dir);
  bookHistos(ibooker,denassTracksMEs_,"ref_ass",dir);

  book_associate_tracks_histos(ibooker,assTracksMEs_,"associate",dir);

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
  //  std::cout << "map: " << map.size() << "[tracks1: " << tracks1.size() << ", tracks2: " << tracks2.size() << "]" << std::endl;

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
Track2TrackValidator::book_associate_tracks_histos(DQMStore::IBooker & ibooker, assME& mes, TString label, std::string & dir){

  ibooker.cd();
  ibooker.setCurrentFolder(dir);

  (mes.h_hits_vs_hits)  = ibooker.book2D(label+"_hits_vs_hits","monitored track # hits vs reference track # hits", 35, -0.5, 34.5, 35,-0.5, 34.5);
  (mes.h_pt_vs_pt)      = ibooker.book2D(label+"_pt_vs_pt",    "monitored track p_{T} vs reference track p_{T}",    nintPt,  minPt,  maxPt,  nintPt,  minPt,  maxPt);
  (mes.h_eta_vs_eta)    = ibooker.book2D(label+"_eta_vs_eta",  "monitored track #eta vs reference track #eta",     nintEta, minEta, maxEta, nintEta, minEta, maxEta);
  (mes.h_phi_vs_phi)    = ibooker.book2D(label+"_phi_vs_phi",  "monitored track #phi vs reference track #phi",     nintPhi, minPhi, maxPhi, nintPhi, minPhi, maxPhi);

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
Track2TrackValidator::fill_associate_tracks_histos(assME& mes, reco::Track* mon, reco::Track* ref, reco::BeamSpot* bs, reco::Vertex* pv) {

  float mon_pt       = mon->pt();
  float mon_eta      = mon->eta();
  float mon_phi      = mon->phi();
  //  float mon_dxy      = mon->dxy(bs->position());
  //  float mon_dz       = mon->dz(bs->position());
  //  float mon_dxyWRTpv = mon->dxy(pv->position());
  //  float mon_dzWRTpv  = mon->dz(pv->position());
  //  float mon_charge   = mon->charge();
  float mon_nhits    = mon->hitPattern().numberOfValidHits();

  float ref_pt       = ref->pt();
  float ref_eta      = ref->eta();
  float ref_phi      = ref->phi();
  //  float ref_dxy      = ref->dxy(bs->position());
  //  float ref_dz       = ref->dz(bs->position());
  //  float ref_dxyWRTpv = ref->dxy(pv->position());
  //  float ref_dzWRTpv  = ref->dz(pv->position());
  //  float ref_charge   = ref->charge();
  float ref_nhits    = ref->hitPattern().numberOfValidHits();

  (mes.h_hits_vs_hits) -> Fill(ref_nhits,mon_nhits);
  (mes.h_pt_vs_pt    ) -> Fill(ref_pt,   mon_pt);
  (mes.h_eta_vs_eta  ) -> Fill(ref_eta,  mon_eta);
  (mes.h_phi_vs_phi  ) -> Fill(ref_phi,  mon_phi);
  
}

void
Track2TrackValidator::initialize_parameter(const edm::ParameterSet& iConfig)
{

  const edm::ParameterSet& pset = iConfig.getParameter<edm::ParameterSet>("histoPSet");

  //parameters for _vs_eta plots
  minEta     = pset.getParameter<double>("minEta");
  maxEta     = pset.getParameter<double>("maxEta");
  nintEta    = pset.getParameter<int>("nintEta");
  useFabsEta = pset.getParameter<bool>("useFabsEta");

  //parameters for _vs_pt plots
  minPt    = pset.getParameter<double>("minPt");
  maxPt    = pset.getParameter<double>("maxPt");
  nintPt   = pset.getParameter<int>("nintPt");
  useInvPt = pset.getParameter<bool>("useInvPt");
  useLogPt = pset.getUntrackedParameter<bool>("useLogPt",false);

  //parameters for _vs_Hit plots
  minHit  = pset.getParameter<double>("minHit");
  maxHit  = pset.getParameter<double>("maxHit");
  nintHit = pset.getParameter<int>("nintHit");

  //parameters for _vs_Layer plots
  minLayers  = pset.getParameter<double>("minLayers");
  maxLayers  = pset.getParameter<double>("maxLayers");
  nintLayers = pset.getParameter<int>("nintLayers");

  //parameters for _vs_phi plots
  minPhi  = pset.getParameter<double>("minPhi");
  maxPhi  = pset.getParameter<double>("maxPhi");
  nintPhi = pset.getParameter<int>("nintPhi");

  //parameters for _vs_Dxy plots
  minDxy  = pset.getParameter<double>("minDxy");
  maxDxy  = pset.getParameter<double>("maxDxy");
  nintDxy = pset.getParameter<int>("nintDxy");

  //parameters for _vs_Dz plots
  minDz  = pset.getParameter<double>("minDz");
  maxDz  = pset.getParameter<double>("maxDz");
  nintDz = pset.getParameter<int>("nintDz");

  //parameters for _vs_ProductionVertexTransvPosition plots
  minVertpos  = pset.getParameter<double>("minVertpos");
  maxVertpos  = pset.getParameter<double>("maxVertpos");
  nintVertpos = pset.getParameter<int>("nintVertpos");

  //parameters for _vs_ProductionVertexZPosition plots
  minZpos  = pset.getParameter<double>("minZpos");
  maxZpos  = pset.getParameter<double>("maxZpos");
  nintZpos = pset.getParameter<int>("nintZpos");

  //parameters for resolution plots
  ptRes_rangeMin = pset.getParameter<double>("ptRes_rangeMin");
  ptRes_rangeMax = pset.getParameter<double>("ptRes_rangeMax");
  ptRes_nbin = pset.getParameter<int>("ptRes_nbin");

  phiRes_rangeMin = pset.getParameter<double>("phiRes_rangeMin");
  phiRes_rangeMax = pset.getParameter<double>("phiRes_rangeMax");
  phiRes_nbin = pset.getParameter<int>("phiRes_nbin");

  cotThetaRes_rangeMin = pset.getParameter<double>("cotThetaRes_rangeMin");
  cotThetaRes_rangeMax = pset.getParameter<double>("cotThetaRes_rangeMax");
  cotThetaRes_nbin = pset.getParameter<int>("cotThetaRes_nbin");

  dxyRes_rangeMin = pset.getParameter<double>("dxyRes_rangeMin");
  dxyRes_rangeMax = pset.getParameter<double>("dxyRes_rangeMax");
  dxyRes_nbin = pset.getParameter<int>("dxyRes_nbin");

  dzRes_rangeMin = pset.getParameter<double>("dzRes_rangeMin");
  dzRes_rangeMax = pset.getParameter<double>("dzRes_rangeMax");
  dzRes_nbin = pset.getParameter<int>("dzRes_nbin");

  useLogPt = pset.getUntrackedParameter<bool>("useLogPt",false);

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
//DEFINE_FWK_MODULE(Track2TrackValidator);
