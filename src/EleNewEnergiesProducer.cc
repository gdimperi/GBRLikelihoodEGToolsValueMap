// -*- C++ -*-
//
// Package:    EleNewEnergiesProducer
// Class:      EleNewEnergiesProducer
// 
/**\class EleNewEnergiesProducer EleNewEnergiesProducer.cc Calibration/EleNewEnergiesProducer/src/EleNewEnergiesProducer.cc

 Description: [one line class summary]

 Implementation:
 TODO:
     - add the error on the  regression energy
     - define the name of the collection
     - implement a switch to select which energies to produce
     - remove unuseful variables


     [Notes on implementation]
*/
//
// Original Author:  Shervin Nourbakhsh,40 1-B24,+41227671643,
//         Created:  Thu Jul  5 20:17:56 CEST 2012
// $Id$
//
//
//#define REGRESSION 3

// system include files
#include <memory>
#include <TString.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
// for the output
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "HiggsAnalysis/GBRLikelihoodEGTools/interface/EGEnergyCorrectorSemiParm.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"

//
// class declaration
//

class EleNewEnergiesProducer : public edm::EDProducer {
  
  typedef float v_t;
  typedef edm::ValueMap<v_t> NewEnergyMap;
  typedef reco::Photon object_t;
  //typedef std::vector<pat::Photon>  object_t;                 
   public:
      explicit EleNewEnergiesProducer(const edm::ParameterSet&);
      ~EleNewEnergiesProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  // ----------member data ---------------------------

private:

  bool isMC;
  // input tag for primary vertex
  edm::InputTag vtxCollectionTAG;
  edm::InputTag BeamSpotTAG;
  
  /// input tag for electrons
  edm::InputTag electronsTAG;
  edm::InputTag recHitCollectionEBTAG;
  edm::InputTag recHitCollectionEETAG;

  /// input rho
  edm::InputTag rhoTAG;

  edm::InputTag conversionsProducerTAG;

  std::string foutName;
  TString r9weightsFilename;
  std::string puWeightFile;

private:
  // Handle to the electron collection
  edm::Handle<std::vector< object_t> > photonsHandle; //test giulia to be used for 1rst step reco::Photons
  
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitCollectionEBHandle;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitCollectionEEHandle;
  edm::Handle<reco::BeamSpot> bsHandle;
  edm::Handle<reco::VertexCollection> primaryVertexHandle;
  edm::ESHandle<CaloTopology> topologyHandle;
  edm::Handle<double> rhoHandle;
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  edm::Handle<reco::ConversionCollection> conversionsHandle;

  //------------------------------ Josh's regression (Hgg)


  EGEnergyCorrectorSemiParm corV5_ele;
  EGEnergyCorrectorSemiParm corV5_pho;
  std::string regrEleJoshV5_SemiParamFile;
  std::string regrPhoJoshV5_SemiParamFile;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
EleNewEnergiesProducer::EleNewEnergiesProducer(const edm::ParameterSet& iConfig):
  vtxCollectionTAG(iConfig.getParameter<edm::InputTag>("vertexCollection")),
  BeamSpotTAG(iConfig.getParameter<edm::InputTag>("BeamSpotCollection")),
  electronsTAG(iConfig.getParameter<edm::InputTag>("electronCollection")), 
  recHitCollectionEBTAG(iConfig.getParameter<edm::InputTag>("recHitCollectionEB")),
  recHitCollectionEETAG(iConfig.getParameter<edm::InputTag>("recHitCollectionEE")),
  rhoTAG(iConfig.getParameter<edm::InputTag>("rhoFastJet")),
  conversionsProducerTAG(iConfig.getParameter<edm::InputTag>("conversionCollection")),
  //  foutName(iConfig.getParameter<std::string>("foutName")),
  regrEleJoshV5_SemiParamFile(iConfig.getParameter<std::string>("regrEleJoshV5_SemiParamFile")),
  regrPhoJoshV5_SemiParamFile(iConfig.getParameter<std::string>("regrPhoJoshV5_SemiParamFile"))
  
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
  
//   produces< std::vector<float> >("energySCEleJoshEle");
//   produces< std::vector<float> >("energySCEleJoshPho");

  // this name are hard coded, should be put in the cfi
  //------------------------------
  produces< NewEnergyMap >("energySCEleJoshEleSemiParamV5ecorr");
					   
  produces< NewEnergyMap >("energySCEleJoshPhoSemiParamV5ecorr");
  //------------------------------


 //now do what ever other initialization is needed
  
}


EleNewEnergiesProducer::~EleNewEnergiesProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EleNewEnergiesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/


/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

   std::vector<v_t>  energySCEleJoshEleSemiParamV5_ecorr;
   std::vector<v_t>  energySCEleJoshPhoSemiParamV5_ecorr;


   std::auto_ptr<NewEnergyMap>  energySCEleJoshEleSemiParamV5_ecorr_Map(new NewEnergyMap());

   std::auto_ptr<NewEnergyMap>  energySCEleJoshPhoSemiParamV5_ecorr_Map(new NewEnergyMap());

   //------------------------------ ELECTRON
   iEvent.getByLabel(electronsTAG, photonsHandle);
   //iEvent.getByLabel(electronsTAG, electronsHandle);

   // altrimenti non tira eccezione
   //if(!electronsHandle.isValid()){
   //  std::cerr << "[ERROR] electron collection not found" << std::endl;
   //  return;
   //}
   
   //------------------------------  VERTEX
   iEvent.getByLabel(vtxCollectionTAG, primaryVertexHandle); 
   // for conversions with full vertex fit
   iEvent.getByLabel(BeamSpotTAG, bsHandle);
   //  bs = bsHandle.product();
   
   //------------------------------ RHO
   iEvent.getByLabel(rhoTAG,rhoHandle);
   
   //------------------------------ RECHIT
   //   iEvent.getByLabel(recHitCollectionEBTAG, recHitCollectionEBHandle);
   //   iEvent.getByLabel(recHitCollectionEETAG, recHitCollectionEEHandle);
   
   //------------------------------ CONVERSIONS
   iEvent.getByLabel(conversionsProducerTAG, conversionsHandle);
   
   
   iSetup.get<CaloTopologyRecord>().get(topologyHandle);
   //  if (! pTopology.isValid() ) return;
   //  topology = pTopology.product();
   


  EcalClusterLazyTools lazyTools(iEvent, iSetup, 
				 recHitCollectionEBTAG, 
                                 recHitCollectionEETAG);  



  for(std::vector<object_t>::const_iterator ele_itr = (photonsHandle)->begin(); 
      ele_itr != (photonsHandle)->end(); ele_itr++){
    // the new tree has one event per each electron
    // reject trackerDriven only electrons because I haven't saved the PFClusters                                                                                                 
    //    if(ele_itr->trackerDrivenSeed() && !ele_itr->ecalDrivenSeed()) continue;
    
    //energyScaleSmearer scaler;
    //float scaledEnergy = scaler.scale(*ele_itr);
    //float smearedEnergy = scaler.smear(*ele_itr);
    
    double ecor, sigmaEoverE, cbmean, sigma, alpha1, n1, alpha2, n2, pdfval;
    
    
    corV5_ele.CorrectedEnergyWithErrorV5(*ele_itr, *primaryVertexHandle, *rhoHandle, lazyTools, iSetup,ecor, sigma, alpha1, n1, alpha2, n2, pdfval);
    energySCEleJoshEleSemiParamV5_ecorr.push_back(ecor);
    
    
    corV5_pho.CorrectedEnergyWithErrorV5(*ele_itr, *primaryVertexHandle, *rhoHandle, lazyTools, iSetup,ecor, sigma, alpha1, n1, alpha2, n2, pdfval);
    energySCEleJoshPhoSemiParamV5_ecorr.push_back(ecor);
    
    //  std::cout<<corEle_fra.second<<std::endl;
  }
  
  //prepare product 
  // declare the filler of the ValueMap
  
  NewEnergyMap::Filler energySCEleJoshEleSemiParamV5_ecorr_filler(*energySCEleJoshEleSemiParamV5_ecorr_Map);
  				      						   
  NewEnergyMap::Filler energySCEleJoshPhoSemiParamV5_ecorr_filler (*energySCEleJoshPhoSemiParamV5_ecorr_Map);

  energySCEleJoshEleSemiParamV5_ecorr_filler.insert(photonsHandle, energySCEleJoshEleSemiParamV5_ecorr.begin(), energySCEleJoshEleSemiParamV5_ecorr.end());

  energySCEleJoshPhoSemiParamV5_ecorr_filler.insert(photonsHandle,  energySCEleJoshPhoSemiParamV5_ecorr.begin(), energySCEleJoshPhoSemiParamV5_ecorr.end());

  energySCEleJoshEleSemiParamV5_ecorr_filler.fill();
		 
  energySCEleJoshPhoSemiParamV5_ecorr_filler.fill();

  //------------------------------
  // put the ValueMap in the event
  
  iEvent.put(energySCEleJoshEleSemiParamV5_ecorr_Map,  "energySCEleJoshEleSemiParamV5ecorr");

  iEvent.put(energySCEleJoshPhoSemiParamV5_ecorr_Map,  "energySCEleJoshPhoSemiParamV5ecorr");
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
EleNewEnergiesProducer::beginJob()
{
  
  if (!corV5_ele.IsInitialized()) {
    std::cout << "[STATUS] Initializing V5 regrEle: " << regrEleJoshV5_SemiParamFile <<std::endl;
    corV5_ele.Initialize(regrEleJoshV5_SemiParamFile, 5); //"/afs/cern.ch/user/b/bendavid/CMSSWhgg/CMSSW_5_3_11_patch5/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v5_forest_ph.root");
  }
  
  if (!corV5_pho.IsInitialized()) {
    std::cout << "[STATUS] Initializing V5 regrPho: " << regrPhoJoshV5_SemiParamFile <<std::endl;
    corV5_pho.Initialize(regrPhoJoshV5_SemiParamFile, 5); //"/afs/cern.ch/user/b/bendavid/CMSSWhgg/CMSSW_5_3_11_patch5/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v5_forest_ph.root");
  }


}

// ------------ method called once each job just after ending the event loop  ------------
void 
EleNewEnergiesProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
EleNewEnergiesProducer::beginRun(edm::Run&, edm::EventSetup const& iSetup)
{

}

// ------------ method called when ending the processing of a run  ------------
void 
EleNewEnergiesProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EleNewEnergiesProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EleNewEnergiesProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EleNewEnergiesProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EleNewEnergiesProducer);
