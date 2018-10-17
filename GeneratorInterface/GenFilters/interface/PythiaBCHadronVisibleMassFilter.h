#ifndef PYTHIABCHADRONVISIBLEMASSFILTER_h
#define PYTHIABCHADRONVISIBLEMASSFILTER_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "Pythia8/Pythia.h"

//
// class decleration
//

class PythiaBCHadronVisibleMassFilter : public edm::global::EDFilter<> {
 public:
  explicit PythiaBCHadronVisibleMassFilter(const edm::ParameterSet&);
  ~PythiaBCHadronVisibleMassFilter() override;
  
  
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  bool isBCHadron(const HepMC::GenParticle& gp) const;
  bool isBCMeson (const HepMC::GenParticle& gp) const;
  bool isBCBaryon(const HepMC::GenParticle& gp) const;
  bool isKs      (const HepMC::GenParticle& gp) const;
  bool isKL      (const HepMC::GenParticle& gp) const;
  bool isLambda0 (const HepMC::GenParticle& gp) const;
  bool isCharged (const HepMC::GenParticle& gp) const;
  bool isNeutrino(const HepMC::GenParticle& gp) const;

 private:
  const edm::EDGetTokenT<edm::HepMCProduct> token_;
  std::unique_ptr<Pythia8::Pythia> fLookupGen; // this instance is for accessing particleData information

  double minPt_;
  double minEta_;
  double maxEta_;
  double minVisMass_;
  double maxVisMass_;

};
#endif // PYTHIABCHADRONVISIBLEMASSFILTER_h
