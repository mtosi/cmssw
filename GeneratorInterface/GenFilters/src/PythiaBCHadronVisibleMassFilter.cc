#include "GeneratorInterface/GenFilters/interface/PythiaBCHadronVisibleMassFilter.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>

using namespace edm;
using namespace std;
using namespace Pythia8;


PythiaBCHadronVisibleMassFilter::PythiaBCHadronVisibleMassFilter(const edm::ParameterSet& iConfig) :
   token_   ( consumes<edm::HepMCProduct>( edm::InputTag( iConfig.getUntrackedParameter( "moduleLabel",std::string("generator")),"unsmeared")))
  , minPt_      (iConfig.getParameter<double>("minPt")   )
  , minEta_     (iConfig.getParameter<double>("minEta")  )
  , maxEta_     (iConfig.getParameter<double>("maxEta")  )
  , minVisMass_    (iConfig.getParameter<double>("minVisMass") )
  , maxVisMass_    (iConfig.getParameter<double>("maxVisMass") )
{
  //now do what ever initialization is needed

 // create pythia8 instance to access particle data
  edm::LogInfo("PythiaBCHadronVisibleMassFilter") << "Creating pythia8 instance for particle properties" << endl;
  std::cout  << "Creating pythia8 instance for particle properties" << endl;
  if(!fLookupGen.get()) fLookupGen.reset(new Pythia());
}


PythiaBCHadronVisibleMassFilter::~PythiaBCHadronVisibleMassFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool PythiaBCHadronVisibleMassFilter::filter(edm::StreamID,edm::Event& iEvent, const edm::EventSetup& iSetup) const {

  using namespace edm;

  bool accepted = false;

  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(token_, evt);

  HepMC::GenEvent *genEvent = new HepMC::GenEvent(*(evt->GetEvent()));
  
  for (HepMC::GenEvent::particle_iterator p = genEvent->particles_begin();
       p != genEvent->particles_end(); ++p) {
    
    if ( (*p)->pdg_id() == 2212 && (*p)->status() == 4 ) continue; // proton from beam
    
    //    std::cout << "particle ID: " << (*p)->pdg_id() << " status: " << (*p)->status() << std::endl;
    if ( !isBCHadron(*(*p)) ) continue;

    std::cout  << "found ID: " << (*p)->pdg_id() << " [mass: " << (*p)->momentum().m() << "] pT: " << (*p)->momentum().perp() << " eta: " << (*p)->momentum().eta() << endl;

    // -- check for daugthers
    TLorentzVector vis(0.,0.,0.,0.);

    int ndau = 0;     
    int nOK = 0;     
    if ((*p)->end_vertex()) {	
      for (HepMC::GenVertex::particle_iterator des=(*p)->end_vertex()->particles_begin(HepMC::children);
	   des != (*p)->end_vertex()->particles_end(HepMC::children);
	   ++des) {
	++ndau;       

	edm::LogInfo("PythiaBCHadronVisibleMassFilter") << " ID: " << (*des)->pdg_id()
						<< " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta()  << " mass = " << (*des)->momentum().m() << endl;
	std::cout  << " ID: " << (*des)->pdg_id() << " charge "
		   << " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta() << " mass = " << (*des)->momentum().m() << endl;
      	  
	if ( isNeutrino(*(*des)) ) continue;
	if ( isKL(*(*des)) ) continue;

	if ((*des)->momentum().perp() <  minPt_  ) continue;
	if ((*des)->momentum().eta()  <  minEta_ ) continue;
	if ((*des)->momentum().eta()  >  maxEta_ ) continue;

	//	if ( !isCharged(*(*des)) )
	//	  if( !isLambda0(*(*des)) && !isKs(*(*des)) ) continue;

	edm::LogInfo("PythiaBCHadronVisibleMassFilter") << "  accepted" << endl;
	std::cout  << "  accepted" << std::endl;
	  
	nOK++;
	TLorentzVector tmp((*des)->momentum().x(),(*des)->momentum().y(),(*des)->momentum().z(),(*des)->momentum().t());
	vis += tmp;
      }
      std::cout << "visible mass: " << vis.M() << std::endl;
      std::cout << "ndau: " << ndau << " nOK: " << nOK << std::endl;
      double deltaM = (*p)->momentum().m()-vis.M();
      if (deltaM<0) std::cout << "visible mass > original mass !!!" << std::endl;
      else if (deltaM < 1. && deltaM > 0.) {
	double relDeltaM = deltaM/(*p)->momentum().m();
	std::cout << "relative difference: " << relDeltaM << std::endl;
      }
      if (!accepted) 
	if (vis.M() >= minVisMass_ && vis.M() <= maxVisMass_) accepted = true;
    }
    if (accepted) break;
  }
  std::cout << "accepted: " << accepted << std::endl;

  return accepted;

}

bool PythiaBCHadronVisibleMassFilter::isBCHadron(const HepMC::GenParticle& gp) const {
  return isBCMeson(gp) || isBCBaryon(gp);
}

bool PythiaBCHadronVisibleMassFilter::isBCMeson(const HepMC::GenParticle& gp) const {
  
  uint32_t pdgid=abs(gp.pdg_id());
  uint32_t hundreds=pdgid%1000;
  if (hundreds>=400 && hundreds<600) {
    return true;
  } else {
    return false;
  }

}

bool PythiaBCHadronVisibleMassFilter::isBCBaryon(const HepMC::GenParticle& gp) const {
  
  uint32_t pdgid=abs(gp.pdg_id());
  uint32_t thousands=pdgid%10000;
  if (thousands>=4000 && thousands <6000) {
    return true;
  } else {
    return false;
  }

}

bool PythiaBCHadronVisibleMassFilter::isLambda0(const HepMC::GenParticle& gp) const {
  uint32_t pdgid=abs(gp.pdg_id());
  uint32_t hundreds=pdgid%1000;
  if (hundreds==122 || hundreds == 124) return true;
  else return false;
}

bool PythiaBCHadronVisibleMassFilter::isKs(const HepMC::GenParticle& gp) const {
  uint32_t pdgid=abs(gp.pdg_id());
  if (pdgid==310) return true;
  else return false;
}

bool PythiaBCHadronVisibleMassFilter::isKL(const HepMC::GenParticle& gp) const {
  uint32_t pdgid=abs(gp.pdg_id());
  if (pdgid==130) return true;
  else return false;
}

bool PythiaBCHadronVisibleMassFilter::isCharged(const HepMC::GenParticle& gp) const {
  uint32_t pdgid=abs(gp.pdg_id());  
  return true;
}

bool PythiaBCHadronVisibleMassFilter::isNeutrino(const HepMC::GenParticle& gp) const {
  uint32_t pdgid=abs(gp.pdg_id());
  if (pdgid == 12) return true;
  else if (pdgid == 14) return true;
  else if (pdgid == 16) return true;
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PythiaBCHadronVisibleMassFilter);
