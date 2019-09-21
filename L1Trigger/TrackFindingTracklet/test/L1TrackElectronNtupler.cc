//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Phase2L1CaloTrig/interface/L1EGCrystalCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackElectronNtupler : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackElectronNtupler(const edm::ParameterSet& iConfig);
  virtual ~L1TrackElectronNtupler();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

protected:

private:

  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  int L1Tk_nPar;        // use 4 or 5 parameter track fit?

  edm::InputTag L1TrackInputTag;         // L1 track collection
  edm::InputTag MCTruthTrackInputTag;    // MC truth collection
  edm::InputTag GenParticleInputTag;     // GenParticles
  edm::InputTag L1EGammaCrystalInputTag; // egamma, barrel, crystals
  edm::InputTag L1EGammaHGCInputTag;     // egamma, HGC

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;  
  edm::EDGetTokenT< reco::GenParticleCollection > genToken_;
  edm::EDGetTokenT< l1t::EGammaBxCollection > egCrystalToken_;
  edm::EDGetTokenT< l1t::EGammaBxCollection > egHGCToken_;


  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;

  // Gen Particles
  std::vector<float>* m_gen_pt;
  std::vector<float>* m_gen_eta;
  std::vector<float>* m_gen_phi;
  std::vector<float>* m_gen_mass;
  std::vector<float>* m_gen_vx;
  std::vector<float>* m_gen_vy;
  std::vector<float>* m_gen_vz;

  // L1 tracks
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_d0;    // (filled if L1Tk_nPar==5, else 999)
  std::vector<float>* m_trk_z0;
  std::vector<int>*   m_trk_charge;
  std::vector<float>* m_trk_chi2;
  std::vector<float>* m_trk_bendchi2;
  std::vector<int>*   m_trk_nstub;
  std::vector<int>*   m_trk_lhits;
  std::vector<int>*   m_trk_dhits;
  std::vector<int>*   m_trk_seed;
  std::vector<int>*   m_trk_genuine;
  std::vector<int>*   m_trk_fake; //0 fake, 1 track from primary interaction, 2 secondary track
  std::vector<int>*   m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_dxy;

  // electron info
  std::vector<float>* m_elec_et;
  std::vector<float>* m_elec_eta;
  std::vector<float>* m_elec_phi;
  std::vector<float>* m_elec_etiso;
  std::vector<int>*   m_elec_hwQual;
  std::vector<int>*   m_elec_looseTkID;
  std::vector<int>*   m_elec_photonID;
  std::vector<int>*   m_elec_hgc;     //0: barrel, 1: HGC

};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackElectronNtupler::L1TrackElectronNtupler(edm::ParameterSet const& iConfig) :
  config(iConfig)
{
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");

  L1TrackInputTag      = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  GenParticleInputTag  = iConfig.getParameter<edm::InputTag>("GenParticleInputTag");
  L1EGammaCrystalInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaCrystalInputTag");
  L1EGammaHGCInputTag  = iConfig.getParameter<edm::InputTag>("L1EGammaHGCInputTag");

  ttTrackToken_        = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);
  genToken_            = consumes< reco::GenParticleCollection >(GenParticleInputTag);
  egCrystalToken_  = consumes< l1t::EGammaBxCollection >(L1EGammaCrystalInputTag);
  egHGCToken_      = consumes< l1t::EGammaBxCollection >(L1EGammaHGCInputTag);

}

/////////////
// DESTRUCTOR
L1TrackElectronNtupler::~L1TrackElectronNtupler()
{
}

//////////
// END JOB
void L1TrackElectronNtupler::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackElectronNtupler::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackElectronNtupler::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackElectronNtupler::beginJob" << endl;

  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;

  // initilize
  m_gen_pt   = new std::vector<float>;
  m_gen_eta  = new std::vector<float>;
  m_gen_phi  = new std::vector<float>;
  m_gen_mass = new std::vector<float>;
  m_gen_vx   = new std::vector<float>;
  m_gen_vy   = new std::vector<float>;
  m_gen_vz   = new std::vector<float>;

  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_charge    = new std::vector<int>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_bendchi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_lhits = new std::vector<int>;
  m_trk_dhits = new std::vector<int>;
  m_trk_seed  = new std::vector<int>;
  m_trk_genuine       = new std::vector<int>;
  m_trk_fake          = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt    = new std::vector<float>;
  m_trk_matchtp_eta   = new std::vector<float>;
  m_trk_matchtp_phi   = new std::vector<float>;
  m_trk_matchtp_z0    = new std::vector<float>;
  m_trk_matchtp_dxy   = new std::vector<float>;

  m_elec_et  = new std::vector<float>;
  m_elec_eta = new std::vector<float>;
  m_elec_phi = new std::vector<float>;
  m_elec_etiso  = new std::vector<float>;
  m_elec_hwQual = new std::vector<int>;
  m_elec_looseTkID = new std::vector<int>;
  m_elec_photonID  = new std::vector<int>;
  m_elec_hgc       = new std::vector<int>;


  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

  eventTree->Branch("gen_pt",    &m_gen_pt);
  eventTree->Branch("gen_eta",   &m_gen_eta);
  eventTree->Branch("gen_phi",   &m_gen_phi);
  eventTree->Branch("gen_mass",  &m_gen_mass);
  eventTree->Branch("gen_vx",    &m_gen_vx);
  eventTree->Branch("gen_vy",    &m_gen_vy);
  eventTree->Branch("gen_vz",    &m_gen_vz);

  eventTree->Branch("trk_pt",    &m_trk_pt);
  eventTree->Branch("trk_eta",   &m_trk_eta);
  eventTree->Branch("trk_phi",   &m_trk_phi);
  eventTree->Branch("trk_d0",    &m_trk_d0);
  eventTree->Branch("trk_z0",    &m_trk_z0);
  eventTree->Branch("trk_charge",    &m_trk_charge);
  eventTree->Branch("trk_chi2",  &m_trk_chi2);
  eventTree->Branch("trk_bendchi2",  &m_trk_bendchi2);
  eventTree->Branch("trk_nstub", &m_trk_nstub);
  eventTree->Branch("trk_lhits", &m_trk_lhits);
  eventTree->Branch("trk_dhits", &m_trk_dhits);
  eventTree->Branch("trk_seed",  &m_trk_seed);
  eventTree->Branch("trk_genuine",      &m_trk_genuine);
  eventTree->Branch("trk_fake",         &m_trk_fake);
  eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
  eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
  eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
  eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
  eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
  eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);

  eventTree->Branch("elec_et",        &m_elec_et);
  eventTree->Branch("elec_eta",       &m_elec_eta);
  eventTree->Branch("elec_phi",       &m_elec_phi);
  eventTree->Branch("elec_etiso",     &m_elec_etiso);
  eventTree->Branch("elec_hwQual",    &m_elec_hwQual);
  eventTree->Branch("elec_looseTkID", &m_elec_looseTkID);
  eventTree->Branch("elec_photonID",  &m_elec_photonID);
  eventTree->Branch("elec_hgc",       &m_elec_hgc);

}


//////////
// ANALYZE
void L1TrackElectronNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    return;
  }

  // clear variables
  m_gen_pt->clear();
  m_gen_eta->clear();
  m_gen_phi->clear();
  m_gen_mass->clear();
  m_gen_vx->clear();
  m_gen_vy->clear();
  m_gen_vz->clear();

  m_trk_pt->clear();
  m_trk_eta->clear();
  m_trk_phi->clear();
  m_trk_d0->clear();
  m_trk_z0->clear();
  m_trk_charge->clear();
  m_trk_chi2->clear();
  m_trk_bendchi2->clear();
  m_trk_nstub->clear();
  m_trk_lhits->clear();
  m_trk_dhits->clear();
  m_trk_seed->clear();
  m_trk_genuine->clear();
  m_trk_fake->clear();
  m_trk_matchtp_pdgid->clear();
  m_trk_matchtp_pt->clear();
  m_trk_matchtp_eta->clear();
  m_trk_matchtp_phi->clear();
  m_trk_matchtp_z0->clear();
  m_trk_matchtp_dxy->clear();

  m_elec_et->clear();
  m_elec_eta->clear();
  m_elec_phi->clear();
  m_elec_etiso->clear();
  m_elec_hwQual->clear();
  m_elec_looseTkID->clear();
  m_elec_photonID->clear();
  m_elec_hgc->clear();



  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  // MC truth association map
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // genparticles
  edm::Handle< std::vector<reco::GenParticle > > GenParticleHandle;
  iEvent.getByToken(genToken_, GenParticleHandle);

  // tracker geometry information
  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  // egamma objects
  edm::Handle<l1t::EGammaBxCollection> egCrystalHandle;
  iEvent.getByToken(egCrystalToken_, egCrystalHandle);

  edm::Handle<l1t::EGammaBxCollection> egHGCHandle;
  iEvent.getByToken(egHGCToken_, egHGCHandle);


  // ----------------------------------------------------------------------------------------------
  // loop over gen particles
  // ----------------------------------------------------------------------------------------------

  if (GenParticleHandle.isValid() ) {

    std::vector<reco::GenParticle>::const_iterator genpartIter ;
    for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {

      int genStatus = genpartIter->status() ;
      int genPdgId = genpartIter->pdgId() ;

      if (genStatus == 1 && abs(genPdgId)==11) {
	float tmp_gen_pt = genpartIter->pt();
	float tmp_gen_eta = genpartIter->eta();
	float tmp_gen_phi = genpartIter->phi();
	float tmp_gen_mass = genpartIter->mass();
	float tmp_gen_vx = genpartIter->vx();
	float tmp_gen_vy = genpartIter->vy();
	float tmp_gen_vz = genpartIter->vz();

	m_gen_pt ->push_back(tmp_gen_pt);
	m_gen_eta->push_back(tmp_gen_eta);
	m_gen_phi->push_back(tmp_gen_phi);
	m_gen_mass->push_back(tmp_gen_mass);
	m_gen_vx->push_back(tmp_gen_vx);
	m_gen_vy->push_back(tmp_gen_vy);
	m_gen_vz->push_back(tmp_gen_vz);

      }
    }
  }


  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  int this_l1track = 0;
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
  for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
    
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
    this_l1track++;
    
    float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
    float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
    float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
    float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm
    float tmp_trk_charge = 0;
    if(iterL1Track->getRInv(L1Tk_nPar)>0){ tmp_trk_charge = 1; }
    else{ tmp_trk_charge = -1; }
    
    float tmp_trk_d0 = -999;
    if (L1Tk_nPar == 5) {
      float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
      float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();
      tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
    }
    
    float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
    float tmp_trk_bendchi2 = iterL1Track->getStubPtConsistency(L1Tk_nPar);
    
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub  = (int) stubRefs.size();
    
    int tmp_trk_seed = (int) iterL1Track->getWedge();


    // ----------------------------------------------------------------------------------------------
    // loop over stubs on tracks to get layer / disk stubs

    int tmp_trk_dhits=0;
    int tmp_trk_lhits=0;

    for (int is=0; is<tmp_trk_nstub; is++) {

      //detID of stub
      DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();
      
      MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
      const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
      Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );
      
      double x=posStub.x();
      double y=posStub.y();
      double z=posStub.z();
      
      int layer=-999999;
      if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
	layer  = static_cast<int>(tTopo->layer(detIdStub));
	tmp_trk_lhits+=pow(10,layer-1);
      }
      else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	layer  = static_cast<int>(tTopo->layer(detIdStub));
	tmp_trk_dhits+=pow(10,layer-1);
      }

    }
 

    // ----------------------------------------------------------------------------------------------
    // MC truth information, is track genuinely matched to Tracking Particle?

    int tmp_trk_genuine = 0;
    if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
    

    // ----------------------------------------------------------------------------------------------
    // fill branches

    m_trk_pt ->push_back(tmp_trk_pt);
    m_trk_eta->push_back(tmp_trk_eta);
    m_trk_phi->push_back(tmp_trk_phi);
    m_trk_z0 ->push_back(tmp_trk_z0);
    m_trk_charge ->push_back(tmp_trk_charge);
    if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
    else m_trk_d0->push_back(999.);
    m_trk_chi2 ->push_back(tmp_trk_chi2);
    m_trk_bendchi2 ->push_back(tmp_trk_bendchi2);
    m_trk_nstub->push_back(tmp_trk_nstub);
    m_trk_dhits->push_back(tmp_trk_dhits);
    m_trk_lhits->push_back(tmp_trk_lhits);
    m_trk_seed->push_back(tmp_trk_seed);
    m_trk_genuine->push_back(tmp_trk_genuine);


    // ----------------------------------------------------------------------------------------------
    // get truth information about the tracking particled possibly matched to the track
    
    edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
    
    int myFake = 0;
    
    int myTP_pdgid = -999;
    float myTP_pt = -999;
    float myTP_eta = -999;
    float myTP_phi = -999;
    float myTP_z0 = -999;
    float myTP_dxy = -999;
    
    if (my_tp.isNull()) myFake = 0;
    else {
      int tmp_eventid = my_tp->eventId().event();
      
      if (tmp_eventid > 0) myFake = 2;
      else myFake = 1;
      
      myTP_pdgid = my_tp->pdgId();
      myTP_pt = my_tp->p4().pt();
      myTP_eta = my_tp->p4().eta();
      myTP_phi = my_tp->p4().phi();
      myTP_z0 = my_tp->vertex().z();
      
      float myTP_x0 = my_tp->vertex().x();
      float myTP_y0 = my_tp->vertex().y();
      myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

    }

    
    // ----------------------------------------------------------------------------------------------
    // fill more branches

    m_trk_fake->push_back(myFake);
    
    m_trk_matchtp_pdgid->push_back(myTP_pdgid);
    m_trk_matchtp_pt->push_back(myTP_pt);
    m_trk_matchtp_eta->push_back(myTP_eta);
    m_trk_matchtp_phi->push_back(myTP_phi);
    m_trk_matchtp_z0->push_back(myTP_z0);
    m_trk_matchtp_dxy->push_back(myTP_dxy);
 

  }//end track loop


  // ----------------------------------------------------------------------------------------------
  // loop over barrel crystal egamma objects
  // ----------------------------------------------------------------------------------------------

  l1t::EGammaBxCollection::const_iterator iterEgCrystal;
  for ( iterEgCrystal = egCrystalHandle->begin(0); iterEgCrystal != egCrystalHandle->end(0); iterEgCrystal++ ) {
    
    float tmp_elec_et = iterEgCrystal->et();
    float tmp_elec_eta = iterEgCrystal->eta();
    float tmp_elec_phi = iterEgCrystal->phi();
    float tmp_elec_hwQual = iterEgCrystal->hwQual();
    float tmp_elec_etiso = iterEgCrystal->isoEt();

    m_elec_et->push_back(tmp_elec_et);
    m_elec_eta->push_back(tmp_elec_eta);
    m_elec_phi->push_back(tmp_elec_phi);
    m_elec_hwQual->push_back(tmp_elec_hwQual);
    m_elec_etiso->push_back(tmp_elec_etiso);

    bool quality = ( ( iterEgCrystal->hwQual() >> 1 ) & 1 ) > 0;
    if (quality) m_elec_looseTkID->push_back(1);
    else m_elec_looseTkID->push_back(0);
    quality = ( ( iterEgCrystal->hwQual() >> 2 ) & 1 ) > 0;
    if (quality) m_elec_photonID->push_back(1);
    else m_elec_photonID->push_back(0);

    m_elec_hgc->push_back(0);
  }

  l1t::EGammaBxCollection::const_iterator iterEgHGC;
  for ( iterEgHGC = egHGCHandle->begin(0); iterEgHGC != egHGCHandle->end(0); iterEgHGC++ ) {
    
    float tmp_elec_et = iterEgHGC->et();
    float tmp_elec_eta = iterEgHGC->eta();
    float tmp_elec_phi = iterEgHGC->phi();
    float tmp_elec_hwQual = iterEgHGC->hwQual();
    float tmp_elec_etiso = iterEgHGC->isoEt();

    m_elec_et->push_back(tmp_elec_et);
    m_elec_eta->push_back(tmp_elec_eta);
    m_elec_phi->push_back(tmp_elec_phi);
    m_elec_hwQual->push_back(tmp_elec_hwQual);
    m_elec_etiso->push_back(tmp_elec_etiso);

    bool quality = (iterEgHGC->hwQual() == 3);
    if (quality){
      m_elec_looseTkID->push_back(1);
      m_elec_photonID->push_back(1);
    }
    else {
      m_elec_looseTkID->push_back(0);
      m_elec_photonID->push_back(0);
    }

    m_elec_hgc->push_back(1);
  }


  eventTree->Fill();


} // end of analyze()


///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackElectronNtupler);
