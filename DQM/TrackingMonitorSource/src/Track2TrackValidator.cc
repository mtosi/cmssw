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

#include "DQM/TrackingMonitorSource/interface/HistoHelper.h"

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

  DQMStore* dqmStore_;

  edm::InputTag numTrackInputTag_;
  edm::InputTag denTrackInputTag_;

  //these are used by MTVGenPs
  edm::EDGetTokenT<reco::TrackCollection> numTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> denTrackToken_;
  edm::EDGetTokenT<reco::BeamSpot> numbsToken_;
  edm::EDGetTokenT<reco::BeamSpot> denbsToken_;

  std::string out;

 private:
  //  edm::ParameterSet conf_;
  std::string topDirName_;
  double dRmin_;

  HistoHelper* histoHelper_;
  HistoHelper::generalME numTracksMEs_;  
  HistoHelper::generalME numassTracksMEs_;  
  HistoHelper::generalME denTracksMEs_;  
  HistoHelper::generalME denassTracksMEs_;  
  
};

#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
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
  
  histoHelper_ = new HistoHelper(iConfig.getParameter<edm::ParameterSet>("histoPSet"));
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
  
  edm::Handle<reco::BeamSpot> numbsHandle;
  iEvent.getByToken(numbsToken_,numbsHandle);
  reco::BeamSpot numbs = *numbsHandle;

  edm::Handle<reco::BeamSpot> denbsHandle;
  iEvent.getByToken(denbsToken_,denbsHandle);
  reco::BeamSpot denbs = *denbsHandle;
  std::cout << "denbs: " << denbs << std::endl;

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

  // loop over tracks (denominator)

  size_t i = 0;

  std::cout << "[Track2TrackValidator::analyze] starting loop over denominator" << std::endl;
  for (idx2idxByDoubleColl::const_iterator pItr = den2numColl.begin(), eItr = den2numColl.end();
       pItr != eItr; ++pItr) {
    
    int trackIdx = pItr->first;
    reco::Track track = denTracks.at(trackIdx);

    histoHelper_->fill_generic_tracks_histos(*&denTracksMEs_,&track,&denbs);
    st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
    
    std::map<double,int> trackDRmap = pItr->second;
    double dRmin = trackDRmap.begin()->first;
    (denTracksMEs_.h_dRmin)->Fill(dRmin);
    bool matched = false;
    if ( dRmin < dRmin_ ) matched = true;
    if ( matched ) {
      histoHelper_->fill_generic_tracks_histos(*&denassTracksMEs_,&track,&numbs);
      (denassTracksMEs_.h_dRmin)->Fill(dRmin);
    }
    for (std::map<double,int>::const_iterator mItr = trackDRmap.begin(), meItr = trackDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dR: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
  }

  // loop over tracks (numerator)
  for (idx2idxByDoubleColl::const_iterator pItr = num2denColl.begin(), eItr = num2denColl.end();
       pItr != eItr; ++pItr) {

    int trackIdx = pItr->first;
    reco::Track track = numTracks.at(trackIdx);
    histoHelper_->fill_generic_tracks_histos(*&numTracksMEs_,&track,&numbs);
   
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
      histoHelper_->fill_generic_tracks_histos(*&numassTracksMEs_,&track,&denbs);
    }
    for (std::map<double,int>::const_iterator mItr = trackDRmap.begin(), meItr = trackDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dR: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
    //    histoHelper_->fill_recoAssociated_simTrack_histos(w,*tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,matchedTrackPointer,puinfo.getPU_NumInteractions(), vtx_z_PU);
    if (matchedTrackIdx >= 0) asts++;

    i++;
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

  histoHelper_->bookHistos(ibooker,numTracksMEs_,   "mon",    dir);
  histoHelper_->bookHistos(ibooker,numassTracksMEs_,"mon_ass",dir);
  histoHelper_->bookHistos(ibooker,denTracksMEs_,   "ref",    dir);
  histoHelper_->bookHistos(ibooker,denassTracksMEs_,"ref_ass",dir);

}

void 
Track2TrackValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
Track2TrackValidator::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  //  histoHelper_->fillHistosFromVectors();
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



//define this as a plug-in
DEFINE_FWK_MODULE(Track2TrackValidator);
