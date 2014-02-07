#ifndef HistoHelper_h
#define HistoHelper_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <TH1.h>
#include <TH2.h>

namespace reco{
  class Vertex;
}

class HistoHelper {
 public:

  HistoHelper(const edm::ParameterSet& iConfig);
  virtual ~HistoHelper() {}

  // to be implemented in the concrete classes
  void initialize();
  void setUpVectors();

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
  const edm::ParameterSet& pset_;

  //private data members
  std::vector<double> ETAintervals;
  std::vector<double> PTintervals;
  std::vector<double> PHIintervals;
  std::vector<double> DXYintervals;
  std::vector<double> DZintervals;
  std::vector<double> VERTEXPOSintervals;
  std::vector<double> ZPOSintervals;

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


};

#endif
