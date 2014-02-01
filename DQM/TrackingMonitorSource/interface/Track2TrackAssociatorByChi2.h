#ifndef Track2TrackAssociatorByChi2_h
#define Track2TrackAssociatorByChi2_h

/** \class Track2TrackAssociatorByChi2
 *  Class that performs the association of reco::Tracks and reco::Tracks evaluating the chi2 of reco tracks parameters and sim tracks parameters.
 *  The cut can be tuned from the config file: see data/Track2TrackAssociatorByChi2.cfi.
 *  Note that the Association Map is filled with -ch2 and not chi2 because it is ordered using std::greater: the track with the lowest association chi2 will be the first in the output map.
 *  It is possible to use only diagonal terms (associator by pulls) seeting onlyDiagonal = true in the PSet 
 *
 *  \author tosi
 */

#include "DQM/TrackingMonitorSource/interface/Track2TrackAssociatorBase.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include<map>

//Note that the Association Map is filled with -ch2 and not chi2 because it is ordered using std::greater:
//the track with the lowest association chi2 will be the first in the output map.

class Track2TrackAssociatorByChi2 : public Track2TrackAssociatorBase {

 public:
  typedef std::map<double,  reco::Track> Chi2RecoMap;
  typedef std::pair< reco::Track, Chi2RecoMap> RecoToRecoPair;
  typedef std::vector< RecoToRecoPair > RecoToRecoPairAssociation;

  /// Constructor with PSet
  Track2TrackAssociatorByChi2(const edm::ESHandle<MagneticField> mF, const edm::ParameterSet& conf)
    : theMF       (mF) 
    , chi2cut     (conf.getParameter<double>       ("chi2cut")     )
    , onlyDiagonal(conf.getParameter<bool>         ("onlyDiagonal"))
    , bsSrc       (conf.getParameter<edm::InputTag>("beamSpot")    )
    {
      if (onlyDiagonal)
	edm::LogInfo("Track2TrackAssociatorByChi2") << " ---- Using Off Diagonal Covariance Terms = 0 ---- " <<  "\n";
      else 
	edm::LogInfo("Track2TrackAssociatorByChi2") << " ---- Using Off Diagonal Covariance Terms != 0 ---- " <<  "\n";
    }

  /// Constructor with magnetic field, double, bool and InputTag
  Track2TrackAssociatorByChi2(const edm::ESHandle<MagneticField> mF, double chi2Cut, bool onlyDiag, const edm::InputTag& beamspotSrc)
    : theMF       (mF)
    , chi2cut     (chi2Cut)
    , onlyDiagonal(onlyDiag)
    , bsSrc       (beamspotSrc)
    {}

  /// Destructor
  ~Track2TrackAssociatorByChi2(){}

  /// compare reco::TrackCollection and reco::TrackCollection iterators: returns the chi2
  double compareTracksParam(reco::TrackCollection::const_iterator, 
			    reco::TrackCollection::const_iterator, 
			    const math::XYZTLorentzVectorD&, 
			    const GlobalVector&,
			    const reco::TrackBase::CovarianceMatrix&,
			    const reco::BeamSpot&) const;

  /// compare collections reco to reco
  RecoToRecoPairAssociation compareTracksParam(const reco::TrackCollection&, 
					       const reco::TrackCollection&, 
					       const reco::VertexCollection&,
					       const reco::BeamSpot&) const;

  /// basic method where chi2 is computed
  double getChi2(reco::TrackBase::ParameterVector& rParameters,
		 reco::TrackBase::CovarianceMatrix& recoTrackCovMatrix,
		 Basic3DVector<double>& momAtVtx,
		 Basic3DVector<double>& vert,
		 int& charge,
		 const reco::BeamSpot&) const;

  /// compare reco::TrackCollection and TrackingParticleCollection iterators: returns the chi2
  double associateRecoToReco(reco::TrackCollection::const_iterator,
			     reco::TrackCollection::const_iterator,
			     const reco::BeamSpot&) const;

  /// propagate the track parameters of TrackinParticle from production vertex to the point of closest approach to the beam line. 
  std::pair<bool,reco::TrackBase::ParameterVector> parametersAtClosestApproach(const Basic3DVector<double>&,// vertex
									       const Basic3DVector<double>&,// momAtVtx
									       float,// charge
									       const reco::BeamSpot&) const;//beam spot
  /// Association Reco To Reco with Collections
  virtual
    reco::RecoToRecoCollection associateRecoToReco(const edm::RefToBaseVector<reco::Track>&,
						   const edm::RefToBaseVector<reco::Track>&,
						   const edm::Event * event = 0,
						   const edm::EventSetup * setup = 0 ) const override;

  /// compare reco to sim the handle of reco::Track and reco::Track collections
  virtual
    reco::RecoToRecoCollection associateRecoToReco(edm::Handle<edm::View<reco::Track> >& t1CH, 
						   edm::Handle<edm::View<reco::Track> >& t2CH, 
						   const edm::Event * event = 0,
						   const edm::EventSetup * setup = 0) const override {
    return Track2TrackAssociatorBase::associateRecoToReco(t1CH,t2CH,event,setup);
  }
  
 private:
  edm::ESHandle<MagneticField> theMF;
  double chi2cut;
  bool onlyDiagonal;
  edm::InputTag bsSrc;
};

#endif
