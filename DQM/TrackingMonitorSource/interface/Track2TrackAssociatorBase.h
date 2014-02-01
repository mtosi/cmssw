#ifndef Track2TrackAssociatorBase_h
#define Track2TrackAssociatorBase_h

/** \class Track2TrackAssociatorBase
 *  Base class for TrackAssociators. Methods take as input the handle of Track and Track collections and return an AssociationMap (oneToManyWithQuality)
 *
 *  \author tosi
 */

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

namespace reco {


  typedef edm::AssociationMap<edm::OneToManyWithQualityGeneric
    <TrackCandidateCollection, edm::View<TrajectorySeed>, double> >
    RecoToRecoCollectionSeed;  

  typedef edm::AssociationMap<edm::OneToManyWithQualityGeneric
    <TrackCandidateCollection, TrackCandidateCollection, double> >
    RecoToRecoCollectionTCandidate;  
}

class Track2TrackAssociatorBase {
 public:
  /// Constructor
  Track2TrackAssociatorBase() {;} 
  /// Destructor
  virtual ~Track2TrackAssociatorBase() {;}


  /// compare reco to reco the handle of reco::Track and TrackingParticle collections
  virtual reco::RecoToRecoCollection associateRecoToReco(edm::Handle<edm::View<reco::Track> >& t1CH, 
							 edm::Handle<edm::View<reco::Track> >& t2CH, 
							 const edm::Event * event ,
							 const edm::EventSetup * setup ) const {
    edm::RefToBaseVector<reco::Track> t1c(t1CH);
    for (unsigned int j=0; j<t1CH->size();j++)
      t1c.push_back(edm::RefToBase<reco::Track>(t1CH,j));
    
    edm::RefToBaseVector<reco::Track> t2c(t2CH);
    for (unsigned int j=0; j<t2CH->size();j++)
      t2c.push_back(edm::RefToBase<reco::Track>(t2CH,j));
    
    return associateRecoToReco(t1c,t2c,event,setup);
  }
  
  /// Association Reco To reco with Collections
  virtual  reco::RecoToRecoCollection associateRecoToReco(const edm::RefToBaseVector<reco::Track> & t1c,
							  const edm::RefToBaseVector<reco::Track> & t2c,
							  const edm::Event * event ,
							  const edm::EventSetup * setup  ) const = 0 ;
  
  //TrajectorySeed
  virtual reco::RecoToRecoCollectionSeed associateRecoToReco(edm::Handle<edm::View<TrajectorySeed> >&, 
							     edm::Handle<edm::View<reco::Track> >&,
							     const edm::Event * event ,
							     const edm::EventSetup * setup ) const {
    reco::RecoToRecoCollectionSeed empty;
    return empty;
  }
  
  //TrackCandidate
  virtual reco::RecoToRecoCollectionTCandidate associateRecoToReco(edm::Handle<TrackCandidateCollection>&, 
								   edm::Handle<TrackCandidateCollection>&, 
								   const edm::Event * event ,
								   const edm::EventSetup * setup ) const {
    reco::RecoToRecoCollectionTCandidate empty;
    return empty;
  }
  
};

#endif
