#include "DQM/TrackingMonitorSource/interface/Track2TrackAssociatorByChi2.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

using namespace edm;
using namespace reco;
using namespace std;

double
Track2TrackAssociatorByChi2::compareTracksParam ( TrackCollection::const_iterator rt, 
						  TrackCollection::const_iterator st, 
						  const math::XYZTLorentzVectorD& vertexPosition, 
						  const GlobalVector& magField,
						  const TrackBase::CovarianceMatrix& invertedCovariance,
						  const reco::BeamSpot& bs) const{
  
  Basic3DVector<double> momAtVtx(st->momentum().x(),st->momentum().y(),st->momentum().z());
  Basic3DVector<double> vert = (Basic3DVector<double>) vertexPosition;

  std::pair<bool,reco::TrackBase::ParameterVector> params = parametersAtClosestApproach(vert, momAtVtx, st->charge(), bs);
  if (params.first){
    TrackBase::ParameterVector sParameters = params.second;
    TrackBase::ParameterVector rParameters = rt->parameters();

    TrackBase::ParameterVector diffParameters = rParameters - sParameters;
    diffParameters[2] = reco::deltaPhi(diffParameters[2],0.f);
    double chi2 = ROOT::Math::Dot(diffParameters * invertedCovariance, diffParameters);
    
    return chi2;
  } else {
    return 10000000000.;
  }
}


Track2TrackAssociatorByChi2::RecoToRecoPairAssociation 
Track2TrackAssociatorByChi2::compareTracksParam(const TrackCollection& rtColl,
						const TrackCollection& stColl,
						const VertexCollection& svColl,
						const reco::BeamSpot& bs) const{
  
  RecoToRecoPairAssociation outputVec;

  for (TrackCollection::const_iterator track=rtColl.begin(); track!=rtColl.end(); track++){
    Chi2RecoMap outMap;

    TrackBase::ParameterVector rParameters = track->parameters();

    TrackBase::CovarianceMatrix recoTrackCovMatrix = track->covariance();
    if (onlyDiagonal){
      for (unsigned int i=0;i<5;i++){
	for (unsigned int j=0;j<5;j++){
	  if (i!=j) recoTrackCovMatrix(i,j)=0;
	}
      }
    }
    recoTrackCovMatrix.Invert();

    for (TrackCollection::const_iterator st=stColl.begin(); st!=stColl.end(); st++){
      
      Basic3DVector<double> momAtVtx(st->momentum().x(),st->momentum().y(),st->momentum().z());
      //      Basic3DVector<double> vert = (Basic3DVector<double>)  svColl[st->vertIndex()].position();
      //      Basic3DVector<double> vert = (Basic3DVector<double>)  svColl[0].position();
      Basic3DVector<double> vert = (Basic3DVector<double>) st->vertex();

      std::pair<bool,reco::TrackBase::ParameterVector> params = parametersAtClosestApproach(vert, momAtVtx, st->charge(), bs);
      if (params.first){
	TrackBase::ParameterVector sParameters = params.second;
      
	TrackBase::ParameterVector diffParameters = rParameters - sParameters;
        diffParameters[2] = reco::deltaPhi(diffParameters[2],0.f);
	double chi2 = ROOT::Math::Dot(diffParameters * recoTrackCovMatrix, diffParameters);
	chi2/=5;
	if (chi2<chi2cut) outMap[chi2]=*st;
      }
    }
    outputVec.push_back(RecoToRecoPair(*track,outMap));
  }
  return outputVec;
}

double
Track2TrackAssociatorByChi2::getChi2(TrackBase::ParameterVector& rParameters,
				     TrackBase::CovarianceMatrix& recoTrackCovMatrix,
				     Basic3DVector<double>& momAtVtx,
				     Basic3DVector<double>& vert,
				     int& charge,
				     const reco::BeamSpot& bs) const{
  
  double chi2;
  
  std::pair<bool,reco::TrackBase::ParameterVector> params = parametersAtClosestApproach(vert, momAtVtx, charge, bs);
  if (params.first){
    TrackBase::ParameterVector sParameters=params.second;
    
    TrackBase::ParameterVector diffParameters = rParameters - sParameters;
    diffParameters[2] = reco::deltaPhi(diffParameters[2],0.f);
    chi2 = ROOT::Math::Dot(diffParameters * recoTrackCovMatrix, diffParameters);
    chi2 /= 5;
    
    LogDebug("Track2TrackAssociatorByChi2") << "====NEW RECO TRACK WITH PT=" << sin(rParameters[1])*float(charge)/rParameters[0] << "====\n" 
					    << "qoverp reco1: " << sParameters[0] << "\n" 
					    << "lambda reco1: " << sParameters[1] << "\n" 
					    << "phi    reco1: " << sParameters[2] << "\n" 
					    << "dxy    reco1: " << sParameters[3] << "\n" 
					    << "dsz    reco1: " << sParameters[4] << "\n" 
					    << ": " /*<< */ << "\n" 
					    << "qoverp reco2: " << rParameters[0] << "\n" 
					    << "lambda reco2: " << rParameters[1] << "\n" 
					    << "phi    reco2: " << rParameters[2] << "\n" 
					    << "dxy    reco2: " << rParameters[3] << "\n" 
					    << "dsz    reco2: " << rParameters[4] << "\n" 
					    << ": " /*<< */ << "\n" 
					    << "chi2: " << chi2 << "\n";
    
    return chi2;  
  } else {
    return 10000000000.;
  }
}


double
Track2TrackAssociatorByChi2::associateRecoToReco( TrackCollection::const_iterator rt, 
						  TrackCollection::const_iterator tp, 
						  const reco::BeamSpot& bs) const{  

  TrackBase::ParameterVector rParameters = rt->parameters();
  TrackBase::CovarianceMatrix recoTrackCovMatrix = rt->covariance();
  if (onlyDiagonal){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	if (i!=j) recoTrackCovMatrix(i,j)=0;
      }
    }
  } 
  
  recoTrackCovMatrix.Invert();
  Basic3DVector<double> momAtVtx(tp->momentum().x(),tp->momentum().y(),tp->momentum().z());
  Basic3DVector<double> vert(tp->vertex().x(),tp->vertex().y(),tp->vertex().z());
  int charge = tp->charge();
  return getChi2(rParameters,recoTrackCovMatrix,momAtVtx,vert,charge,bs);
}

pair<bool,TrackBase::ParameterVector> 
Track2TrackAssociatorByChi2::parametersAtClosestApproach(const Basic3DVector<double>& vertex,
							 const Basic3DVector<double>& momAtVtx,
							 float charge,
							 const BeamSpot& bs) const{
  
  TrackBase::ParameterVector sParameters;
  try {
    FreeTrajectoryState ftsAtProduction(GlobalPoint(vertex.x(),vertex.y(),vertex.z()),
					GlobalVector(momAtVtx.x(),momAtVtx.y(),momAtVtx.z()),
					TrackCharge(charge),
					theMF.product());
    TSCBLBuilderNoMaterial tscblBuilder;
    TrajectoryStateClosestToBeamLine tsAtClosestApproach = tscblBuilder(ftsAtProduction,bs);//as in TrackProducerAlgorithm
    
    GlobalPoint v = tsAtClosestApproach.trackStateAtPCA().position();
    GlobalVector p = tsAtClosestApproach.trackStateAtPCA().momentum();
    sParameters[0] = tsAtClosestApproach.trackStateAtPCA().charge()/p.mag();
    sParameters[1] = Geom::halfPi() - p.theta();
    sParameters[2] = p.phi();
    sParameters[3] = (-v.x()*sin(p.phi())+v.y()*cos(p.phi()));
    sParameters[4] = v.z()*p.perp()/p.mag() - (v.x()*p.x()+v.y()*p.y())/p.perp() * p.z()/p.mag();
    
    return pair<bool,TrackBase::ParameterVector>(true,sParameters);
  } catch ( ... ) {
    return pair<bool,TrackBase::ParameterVector>(false,sParameters);
  }
}

RecoToRecoCollection 
Track2TrackAssociatorByChi2::associateRecoToReco(const edm::RefToBaseVector<reco::Track>& t1C, 
						 const edm::RefToBaseVector<reco::Track>& t2C,
						 const edm::Event * e,
						 const edm::EventSetup *setup ) const{

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  e->getByLabel(bsSrc,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;      

  RecoToRecoCollection  outputCollection;

  int t1index=0;
  for (RefToBaseVector<reco::Track>::const_iterator r1t=t1C.begin(); r1t!=t1C.end(); r1t++, t1index++){
    
    LogDebug("TrackAssociator") << "=========LOOKING FOR ASSOCIATION===========" << "\n"
				<< "rec::Track #"<<t1index<<" with pt=" << (*r1t)->pt() <<  "\n"
				<< "===========================================" << "\n";
    
    TrackBase::ParameterVector rParameters = (*r1t)->parameters();

    TrackBase::CovarianceMatrix recoTrackCovMatrix = (*r1t)->covariance();
    if (onlyDiagonal){
      for (unsigned int i=0;i<5;i++){
	for (unsigned int j=0;j<5;j++){
	  if (i!=j) recoTrackCovMatrix(i,j)=0;
	}
      }
    } 

    recoTrackCovMatrix.Invert();

    int t2index =0;
  for (RefToBaseVector<reco::Track>::const_iterator r2t=t2C.begin(); r2t!=t2C.end(); r2t++, t2index++){
	
      //skip tps with a very small pt
      //if (sqrt(tp->momentum().perp2())<0.5) continue;
      int charge = (*r2t)->charge();
      if (charge==0) continue;
      Basic3DVector<double> momAtVtx((*r2t)->momentum().x(),(*r2t)->momentum().y(),(*r2t)->momentum().z());
      Basic3DVector<double> vert=(Basic3DVector<double>) (*r2t)->vertex();

      double chi2 = getChi2(rParameters,recoTrackCovMatrix,momAtVtx,vert,charge,bs);
      
      if (chi2<chi2cut) {
	outputCollection.insert(t1C[t1index], 
				std::make_pair(t2C[t2index],
					       -chi2));//-chi2 because the Association Map is ordered using std::greater
      }
    }
  }
  outputCollection.post_insert();
  return outputCollection;
}

