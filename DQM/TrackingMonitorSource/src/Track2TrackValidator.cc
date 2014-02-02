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

#include "DQM/TrackingMonitorSource/interface/Track2TrackAssociatorBase.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelector.h"
#include "CommonTools/RecoAlgos/interface/CosmicTrackingParticleSelector.h"

#include <DQMServices/Core/interface/DQMEDAnalyzer.h>

#include <iostream>
#include <sstream>
#include <string>

//
// class declaration
//
class DQMStore;
class Track;
class DQMStore;
class HistoHelper;

class Track2TrackValidator : public DQMEDAnalyzer {
   public:
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
  
  DQMStore* dqmStore_;

  edm::InputTag numTrackInputTag_;
  edm::InputTag denTrackInputTag_;

  //these are used by MTVGenPs
  edm::EDGetTokenT<reco::Track> numTrackToken_;
  edm::EDGetTokenT<reco::Track> denTrackToken_;
  edm::EDGetTokenT<reco::RecoToRecoCollection> associatorMapToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;

  bool useAssociator_;
  std::string associators;
  edm::InputTag associatorMap_;
  const Track2TrackAssociatorBase* associator;

  std::string out;


  edm::ESHandle<MagneticField> theMF;

 private:
  std::string topDirName_;

  // select tracking particles 
  //(i.e. "denominator" of the efficiency ratio)
  //  TrackSelector tpSelector;				      
  edm::InputTag hitMapInputTag_;

  HistoHelper* histoHelper_;

};

#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DQM/TrackingMonitorSource/interface/HistoHelper.h"
//#include "SimTracker/Records/interface/TrackAssociatorRecord.h"


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
  , useAssociator_    ( iConfig.getParameter<bool>         ("useAssociator")  )
  , associatorMap_    ( iConfig.getParameter<edm::InputTag>("associatorMap")  )
  , topDirName_       ( iConfig.getParameter<std::string>  ("topDirName")     )
  , hitMapInputTag_   ( iConfig.getParameter<edm::InputTag>("hitMapInputTag") )
{
   //now do what ever initialization is needed
  numTrackToken_      = consumes<reco::Track>(numTrackInputTag_);
  denTrackToken_      = consumes<reco::Track>(denTrackInputTag_);
  associatorMapToken_ = consumes<reco::RecoToRecoCollection>(associatorMap_);
  bsToken_            = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  
  associators = associatorMap_.label();
  
  histoHelper_ = new HistoHelper(iConfig.getParameter<edm::ParameterSet>("histoPSet"));

}


Track2TrackValidator::~Track2TrackValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void
Track2TrackValidator::beginJob(const edm::EventSetup& iSetup) {
  if (useAssociator_) {
    //    edm::ESHandle<Track2TrackAssociatorBase> theAssociator;
    //    iSetup.get<TrackAssociatorRecord>().get(associators,theAssociator);
    //    associator = theAssociator.product();
  }
}

//
// member functions
//
// ------------ method called for each event  ------------
void
Track2TrackValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_,bsHandle);
  reco::BeamSpot bs = *bsHandle;

  edm::Handle<edm::View<reco::Track> > numTracks;
  iEvent.getByToken(numTrackToken_,numTracks);

  edm::Handle<edm::View<reco::Track> > denTracks;
  iEvent.getByToken(denTrackToken_,denTracks);

  reco::RecoToRecoCollection reco1TOreco2Coll;
  reco::RecoToRecoCollection reco2TOreco1Coll;
  if (useAssociator_) {
    edm::LogVerbatim("Track2TrackValidator") << "analyzing "
					     << numTrackInputTag_.process()  << ":"
					     << numTrackInputTag_.label()    << ":"
					     << numTrackInputTag_.instance() << " w.r.t. "
					     << denTrackInputTag_.process()  << ":"
					     << denTrackInputTag_.label()    << ":"
					     << denTrackInputTag_.instance() << " \n";

    reco1TOreco2Coll = associator->associateRecoToReco(numTracks,denTracks,&iEvent,&iSetup);
    reco2TOreco1Coll = associator->associateRecoToReco(denTracks,numTracks,&iEvent,&iSetup);
  } else {
    
    edm::Handle<reco::RecoToRecoCollection> reco1TOreco2CollHandle;
    iEvent.getByToken(associatorMapToken_,reco1TOreco2CollHandle);
    reco1TOreco2Coll = *(reco1TOreco2CollHandle.product());
  }


  // ########################################################
  // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
  // ########################################################
  
  //compute number of tracks per eta interval
  //
  edm::LogVerbatim("Track2TrackValidator") << "\n# of tracks (denominator): " << denTracks->size() << "\n";
  int st(0);  	     //This counter counts the number of tracks
  int ast(0);  	     //This counter counts the number of tracks that are "associated" to recoTracks
  unsigned sts(0);   //This counter counts the number of tracks surviving the bunchcrossing cut
  unsigned asts(0);  //This counter counts the number of tracks that are "associated" to recoTracks surviving the bunchcrossing cut

  // loop over tracks (denominator)
  edm::View<reco::Track>::size_type denTracksN = denTracks->size();
  for ( edm::View<reco::Track>::size_type i=0; i < denTracksN; i++ ) {

    edm::RefToBase<reco::Track> track(denTracks, i);
    double dxy(0);
    double dz(0);
    
    reco::Track::Vector momentumTP = track->momentum();
    reco::Track::Point  vertexTP   = track->vertex();
    //Calcualte the impact parameters w.r.t. PCA
    dxy = track->dxy(bs.position());
    dz  = track->dz(bs.position());
    std::cout << "dxy: " << dxy << " dz: " << dz << " vertexTP: " << vertexTP << std::endl;
   
    st++;   //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
    
    // in the coming lines, histos are filled using as input
    // - momentumTP
    // - vertexTP
    // - dxySim
    // - dzSim
    
    //    histoHelper_->fill_generic_simTrack_histos(momentumTP,vertexTP, tp->eventId().bunchCrossing());
    
    
    // ######################################
    // fill RecoAssociated Tracks' histograms
    // ######################################
    const reco::Track* matchedTrackPointer=0;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
    if(reco2TOreco1Coll.find(track) != reco2TOreco1Coll.end()){
      rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) reco2TOreco1Coll[track];

      if (rt.size()!=0) {
	ast++; //This counter counts the number of tracks that have a recoTrack associated
	matchedTrackPointer = rt.begin()->first.get();
	edm::LogVerbatim("Track2TrackValidator") << "Track #" << st
						 << " with pt=" << sqrt(momentumTP.perp2())
						 << " associated with quality:" << rt.begin()->second <<"\n";
      }
    } else{
      edm::LogVerbatim("TrackValidator")
	<< "Track #" << st
	<< " with pt,eta,phi: "
	<< sqrt(momentumTP.perp2()) << " , "
	<< momentumTP.eta() << " , "
	<< momentumTP.phi() << " , "
	<< " NOT associated to any reco::Track" << "\n";
    }
    
    int nRecoHits = track->hitPattern().numberOfHits();
    std::cout << "nRecoHits: " << nRecoHits << std::endl;

    //    histoHelper_->fill_recoAssociated_simTrack_histos(w,*tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,matchedTrackPointer,puinfo.getPU_NumInteractions(), vtx_z_PU);
    sts++;
    if (matchedTrackPointer) asts++;
    
  } // end loop over tracks (denominator)
  
  // ##############################################
  // fill reco Tracks histograms (LOOP OVER TRACKS)
  // ##############################################
  edm::LogVerbatim("Track2TrackValidator") << "\n# of reco::Tracks with "
					   << numTrackInputTag_.process()<<":"
					   << numTrackInputTag_.label()<<":"
					   << numTrackInputTag_.instance()
					   << ": " << numTracks->size() << "\n";
  
  //  int sat(0); //This counter counts the number of recoTracks that are associated to Tracks from Signal only
  int at(0);  //This counter counts the number of recoTracks that are associated to Tracks
  int rT(0);  //This counter counts the number of recoTracks in general
  
  
  edm::View<reco::Track>::size_type numTracksN = numTracks->size();
  for ( edm::View<reco::Track>::size_type i=0; i < numTracksN; i++ ) {

    edm::RefToBase<reco::Track> track(numTracks, i);

    rT++;
    
    //    bool isSigMatched(false);
    //    bool isMatched(false);
    //    bool isChargeMatched(true);

    int numAssocRecoTracks = 0;
    //    int tpbx = 0;
    int nHits = 0;
    double sharedFraction = 0.;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > tp;
    if(reco2TOreco1Coll.find(track) != reco2TOreco1Coll.end()){
      tp = reco1TOreco2Coll[track];
      if (tp.size()!=0) {
	nHits = tp[0].first->hitPattern().numberOfHits();
	sharedFraction = tp[0].second;
	std::cout << "nHits: " << nHits << " sharedFraction: " << sharedFraction << std::endl;
	//	isMatched = true;
	//        if (tp[0].first->charge() != track->charge()) 
	//	  isChargeMatched = false;
        if(reco1TOreco2Coll.find(tp[0].first) != reco1TOreco2Coll.end()) {
	  numAssocRecoTracks = reco1TOreco2Coll[tp[0].first].size();
	}
        std::cout << numAssocRecoTracks << std::endl;
	at++;
	for (size_t tp_ite=0;tp_ite<tp.size();++tp_ite){
	  reco::Track trackpart = *(tp[tp_ite].first);
	}
	edm::LogVerbatim("Track2TrackValidator") << "reco::Track #" << rT << " with pt=" << track->pt()
						 << " associated with quality:" << tp.begin()->second <<"\n";
      }
    } else {
      edm::LogVerbatim("Track2TrackValidator") << "reco::Track #" << rT << " with pt=" << track->pt()
					       << " NOT associated to any TrackingParticle" << "\n";
    }
    
    
    //    histoHelper_->fill_generic_recoTrack_histos(w,*track,bs.position(),isSimMatched,isSigSimMatched, isChargeMatched, numAssocRecoTracks, puinfo.getPU_NumInteractions(), tpbx, nSimHits, sharedFraction);
    

    //Fill other histos
    if (tp.size()==0) continue;
    
    //    histoHelper_->fill_simAssociated_recoTrack_histos(w,*track);
    
    edm::RefToBase<reco::Track> tpr = tp.begin()->first;
    
    //Get tracking particle parameters at point of closest approach to the beamline
    reco::Track::Vector momentumTP = tpr->momentum();
    reco::Track::Point vertexTP    = tpr->vertex();
    int chargeTP = tpr->charge();
    std::cout << "vertexTP: " << vertexTP << " chargeTP: " << chargeTP << " momentumTP: " << momentumTP << std::endl;

    //    histoHelper_->fill_ResoAndPull_recoTrack_histos(w,momentumTP,vertexTP,chargeTP,
    //							  *track,bs.position());
    
  } // end loop over tracks (numerator)
  
  //  histoHelper_->fill_trackBased_histos(w,at,rT,st);
  
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
  histoHelper_->bookHistos(ibooker);
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

//define this as a plug-in
DEFINE_FWK_MODULE(Track2TrackValidator);
