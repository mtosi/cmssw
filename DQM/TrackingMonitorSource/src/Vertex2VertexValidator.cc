#include "DQM/TrackingMonitorSource/interface/Vertex2VertexValidator.h"
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
Vertex2VertexValidator::Vertex2VertexValidator(const edm::ParameterSet& iConfig)
  : numVertexInputTag_ ( iConfig.getParameter<edm::InputTag>("numVertex")       )
  , denVertexInputTag_ ( iConfig.getParameter<edm::InputTag>("denVertex")       )
  , topDirName_       ( iConfig.getParameter<std::string>  ("topDirName")     )
  , dzmin_            ( iConfig.getParameter<double>("dzmin")                 )
{
   //now do what ever initialization is needed
  numpvToken_         = consumes<reco::VertexCollection>(numVertexInputTag_);
  denpvToken_         = consumes<reco::VertexCollection>(denVertexInputTag_);
  
  numVertexMEs_.label    = numVertexInputTag_.label();
  numassVertexMEs_.label = numVertexInputTag_.label()+"_ass";
  denVertexMEs_.label    = denVertexInputTag_.label();
  denassVertexMEs_.label = denVertexInputTag_.label()+"_ass";

  assVertexMEs_.label = "associate";
}


Vertex2VertexValidator::~Vertex2VertexValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void
Vertex2VertexValidator::beginJob(const edm::EventSetup& iSetup) {
}

//
// member functions
//
// ------------ method called for each event  ------------
void
Vertex2VertexValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle<reco::VertexCollection> numpvHandle;
  iEvent.getByToken(numpvToken_,numpvHandle);
  reco::VertexCollection numVertex = *numpvHandle;

  edm::Handle<reco::VertexCollection> denpvHandle;
  iEvent.getByToken(denpvToken_,denpvHandle);
  reco::VertexCollection denVertex = *denpvHandle;

  float numVertexN = numVertex.size();
  float denVertexN = denVertex.size();
  h_nvertices_vs_nvertices -> Fill(denVertexN,numVertexN);

  edm::LogVerbatim("Vertex2VertexValidator") << "analyzing "
					   << numVertexInputTag_.process()  << ":"
					   << numVertexInputTag_.label()    << ":"
					   << numVertexInputTag_.instance() << " w.r.t. "
					   << denVertexInputTag_.process()  << ":"
					   << denVertexInputTag_.label()    << ":"
					   << denVertexInputTag_.instance() << " \n";
  
  /*
  std::cout << "analyzing "
	    << numVertexInputTag_.process()  << ":"
	    << numVertexInputTag_.label()    << ":"
	    << numVertexInputTag_.instance() << " w.r.t. "
	    << denVertexInputTag_.process()  << ":"
	    << denVertexInputTag_.label()    << ":"
	    << denVertexInputTag_.instance() << " \n";
  */
  
  idx2idxByDoubleColl num2denColl;
  fillMap(numVertex,denVertex,num2denColl);

  idx2idxByDoubleColl den2numColl;
  fillMap(denVertex,numVertex,den2numColl);

  // ########################################################

  //compute number of tracks per eta interval
  //
  edm::LogVerbatim("Vertex2VertexValidator") << "\n# of tracks (denominator): " << denVertex.size() << "\n";
  int st(0);  	     //This counter counts the number of tracks
  int ast(0);  	     //This counter counts the number of tracks that are "associated" to recoVertex
  unsigned asts(0);  //This counter counts the number of tracks that are "associated" to recoVertex surviving the bunchcrossing cut

  // loop over vertex (denominator)
  //  std::cout << "[Vertex2VertexValidator::analyze] starting loop over denominator" << std::endl;
  for (idx2idxByDoubleColl::const_iterator pItr = den2numColl.begin(), eItr = den2numColl.end();
       pItr != eItr; ++pItr) {
    
    int vertexIdx = pItr->first;
    reco::Vertex vertex = denVertex.at(vertexIdx);

    fill_generic_vertex_histos(*&denVertexMEs_,&vertex);
    st++;   //This counter counts the number of simulated vertex passing the MTV selection (i.e. tpSelector(tp) )
    
    std::map<double,int> vertexDRmap = pItr->second;
    double dzmin = vertexDRmap.begin()->first;
    (denVertexMEs_.h_dzmin)->Fill(dzmin);
    bool matched = false;
    if ( dzmin < dzmin_ ) matched = true;
    if ( matched ) {
      fill_generic_vertex_histos(*&denassVertexMEs_,&vertex);
      (denassVertexMEs_.h_dzmin)->Fill(dzmin);
      int matchedVertexIndex = vertexDRmap[dzmin];
      reco::Vertex matchedVertex = numVertex.at(matchedVertexIndex);
      fill_associate_vertex_histos(*&assVertexMEs_,&vertex,&matchedVertex);
    }
    for (std::map<double,int>::const_iterator mItr = vertexDRmap.begin(), meItr = vertexDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dz: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
  }

  // loop over vertex (numerator)
  for (idx2idxByDoubleColl::const_iterator pItr = num2denColl.begin(), eItr = num2denColl.end();
       pItr != eItr; ++pItr) {

    int vertexIdx = pItr->first;
    reco::Vertex vertex = numVertex.at(vertexIdx);
    fill_generic_vertex_histos(*&numVertexMEs_,&vertex);
   
    st++;   //This counter counts the number of simulated vertex passing the MTV selection (i.e. tpSelector(tp) )
    
    int matchedVertexIdx = -1;
    std::map<double,int> vertexDRmap = pItr->second;
    double dzmin = vertexDRmap.begin()->first;
    (numVertexMEs_.h_dzmin)->Fill(dzmin);
    //    std::cout << "num2denColl[" <<  vertexIdx << "]" << " map size: " << vertexDRmap.size() << " dzmin: " << dzmin << std::endl;
    bool matched = false;
    if ( dzmin < dzmin_ ) matched = true;
    if ( !matched ) {
      (numassVertexMEs_.h_dzmin)->Fill(dzmin);
      fill_generic_vertex_histos(*&numassVertexMEs_,&vertex);
    }
    for (std::map<double,int>::const_iterator mItr = vertexDRmap.begin(), meItr = vertexDRmap.end();
	 mItr != meItr; ++mItr ) {
      //      std::cout << " --> dz: " <<  mItr->first << " trkIdx: " << mItr->second << std::endl;
    }
    //    fill_recoAssociated_simVertex_histos(w,*tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,matchedVertexPointer,puinfo.getPU_NumInteractions(), vtx_z_PU);
    if (matchedVertexIdx >= 0) asts++;

  } // end loop over vertex (denominator)

  // ##############################################
  // fill reco Vertex histograms (LOOP OVER VERTEX)
  // ##############################################
  edm::LogVerbatim("Vertex2VertexValidator") << "\n# of reco::Vertex with "
					   << numVertexInputTag_.process()<<":"
					   << numVertexInputTag_.label()<<":"
					   << numVertexInputTag_.instance()
					   << ": " << numVertex.size() << "\n";
  
  int at(0);  //This counter counts the number of recoVertex that are associated to Vertex
  int rT(0);  //This counter counts the number of recoVertex in general
  
  edm::LogVerbatim("Vertex2VertexValidator") << "Total vertex: "        << st << "\n"
					   << "Total Associated: "    << ast << "\n"
					   << "Total Reconstructed: " << rT << "\n"
					   << "Total Associated: "    << at << "\n"
					   << "Total Fakes: "         << rT-at << "\n";
  
}

void 
Vertex2VertexValidator::bookHistograms(DQMStore::IBooker & ibooker,
				     edm::Run const & iRun,
				     edm::EventSetup const & iSetup)
{
  
  //  std::cout << "[Vertex2VertexValidator::bookHistograms]" << std::endl;

  std::string dir = topDirName_;

  bookHistos(ibooker,numVertexMEs_,   "mon",    dir);
  bookHistos(ibooker,numassVertexMEs_,"mon_ass",dir);
  bookHistos(ibooker,denVertexMEs_,   "ref",    dir);
  bookHistos(ibooker,denassVertexMEs_,"ref_ass",dir);

  book_associate_vertex_histos(ibooker,assVertexMEs_,"associate",dir);

}

void 
Vertex2VertexValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
Vertex2VertexValidator::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  //  fillHistosFromVectors();
  if ( out.size() != 0 && dqmStore_ ) dqmStore_->save(out);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Vertex2VertexValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
Vertex2VertexValidator::fillMap(reco::VertexCollection vertex1, reco::VertexCollection vertex2, idx2idxByDoubleColl& map)
{
  int i = 0;
  for ( auto vtx1 : vertex1 ) {
    std::map<double,int> tmp;
    int j = 0;
    for ( auto vtx2 : vertex2 ) {
      double dz = fabs(vtx1.z()-vtx2.z());
      tmp[dz] = j;
      j++;
    }
    map.push_back(std::make_pair(i,tmp));
    i++;
  }
  //  std::cout << "map: " << map.size() << "[vertex1: " << vertex1.size() << ", vertex2: " << vertex2.size() << "]" << std::endl;

}

void 
Vertex2VertexValidator::bookHistos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

  h_nvertices_vs_nvertices = ibooker.book2D("nvertices_vs_nvertices","monitored # vertex vs reference # vertex", 100, -0.5,   99.5,   100,-0.5,    99.5);
  book_generic_vertex_histos(ibooker,mes,label,dir);

}

void 
Vertex2VertexValidator::book_generic_vertex_histos(DQMStore::IBooker & ibooker, generalME& mes, TString label, std::string & dir){

  ibooker.cd();
  ibooker.setCurrentFolder(dir);

  (mes.h_ntracks ) = ibooker.book1D(label+"_ntracks", "number of track",       100, -0.5,    99.5);
  (mes.h_sumpt2)   = ibooker.book1D(label+"_sumpt2",  "tracks sum p_{T}^{2}", 1000,  0.1, 10000.);
  (mes.h_dzmin)    = ibooker.book1D(label+"_dzmin",   "delta z min",           300,  0.,     30.);

}

void 
Vertex2VertexValidator::book_associate_vertex_histos(DQMStore::IBooker & ibooker, assME& mes, TString label, std::string & dir){

  ibooker.cd();
  ibooker.setCurrentFolder(dir);

  (mes.h_ntracks_vs_ntracks) = ibooker.book2D(label+"_ntracks_vs_ntracks","monitored vertex # track vs reference vertex # track",                   100, -0.5,   99.5,   100,-0.5,    99.5);
  (mes.h_sumpt2_vs_sumpt2)   = ibooker.book2D(label+"_sumpt2_vs_sumpt2",  "monitored vertex #sum p_{T}^{2} vs reference vertex #sum p_{T}^{2}",    1000,  0.1, 10000.,  1000, 0.1, 10000.);

}


void
Vertex2VertexValidator::fill_generic_vertex_histos(generalME& mes, reco::Vertex* vtx)
{

  float ntracks  = vtx->nTracks();
  float sumpt2   = 0.;
  for(reco::Vertex::trackRef_iterator it = vtx->tracks_begin(),  et = vtx->tracks_end();
      it != et; ++it) {

    bool isHighPurity = (**it).quality(reco::TrackBase::highPurity);
    if ( !isHighPurity ) continue;
    float pt = (**it).pt();  
    sumpt2 += pt*pt;
  }

  (mes.h_ntracks) -> Fill(ntracks);
  (mes.h_sumpt2)  -> Fill(sumpt2);
  
}

void
Vertex2VertexValidator::fill_associate_vertex_histos(assME& mes, reco::Vertex* mon, reco::Vertex* ref)
{

  float mon_ntracks = mon->nTracks();
  float ref_ntracks = ref->nTracks();

  float mon_sumpt2   = 0.;
  for(reco::Vertex::trackRef_iterator it = mon->tracks_begin(),  et = mon->tracks_end();
      it != et; ++it) {

    bool isHighPurity = (**it).quality(reco::TrackBase::highPurity);
    if ( !isHighPurity ) continue;
    float pt = (**it).pt();  
    mon_sumpt2 += pt*pt;
  }
  float ref_sumpt2   = 0.;
  for(reco::Vertex::trackRef_iterator it = ref->tracks_begin(),  et = ref->tracks_end();
      it != et; ++it) {

    bool isHighPurity = (**it).quality(reco::TrackBase::highPurity);
    if ( !isHighPurity ) continue;
    float pt = (**it).pt();  
    ref_sumpt2 += pt*pt;
  }

  (mes.h_ntracks_vs_ntracks) -> Fill(ref_ntracks,mon_ntracks);
  (mes.h_sumpt2_vs_sumpt2)   -> Fill(ref_sumpt2, mon_sumpt2);
}

void
Vertex2VertexValidator::initialize_parameter(const edm::ParameterSet& iConfig)
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
