// -*- C++ -*-
//
// Package:    QstartToGJ
// Class:      GENSIMAnalysis
// 
/**\class QstartToGJ GENSIMAnalysis.cc QstartToGJ/QstartToGJ/src/GENSIMAnalysis.cc

	Check singal MC of Qstar -> gamma + jet in GEN-SIM samples

 */
//
// Original Author:  Jui-FA Tsai,,067231805,0917970159
//

// system include files
#include <memory>
#include <cmath>                                                   

// user include files
#include "QstartToGJ/QstartToGJ/interface/GENSIMAnalysis.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"              
#include "CommonTools/UtilAlgos/interface/TFileService.h"          

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"     

//ROOT
#include <TFile.h>                                                
#include <TTree.h>                                                
#include <TH1D.h>                                                
#include <TH2D.h>                                               

// 
// class declaration
// 
GENSIMAnalysis::GENSIMAnalysis(const edm::ParameterSet& iConfig)
{
	h_pdgId              = tFileService->make<TH1D>("pdgId", "",   100, -50, 50);
	h_ndstar             = tFileService->make<TH1D>("NumDstar", "", 20, -10, 10);
	h_nustar             = tFileService->make<TH1D>("NumUstar", "", 20, -10, 10);
	h_mdstar             = tFileService->make<TH1D>("MassDstar", "", 2000, 0, 2000);
	h_mustar             = tFileService->make<TH1D>("MassUstar", "", 2000, 0, 2000);
	h_pTdstar            = tFileService->make<TH1D>("pTDstar", "", 1000, 0, 1000);
	h_pTustar            = tFileService->make<TH1D>("pTUstar", "", 1000, 0, 1000);
	h_etadstar           = tFileService->make<TH1D>("EtaDstar", "", 100, -5, 5);
	h_etaustar           = tFileService->make<TH1D>("EtaUstar", "", 100, -5, 5);
	h_phidstar           = tFileService->make<TH1D>("PhiDstar", "", 400, -4, 4);
	h_phiustar           = tFileService->make<TH1D>("PhiUstar", "", 400, -4, 4);
	h_dstarDecay         = tFileService->make<TH1D>("DstarDecay", "", 100, -50, 50);
	h_ustarDecay         = tFileService->make<TH1D>("UstarDecay", "", 100, -50, 50);
}

GENSIMAnalysis::~GENSIMAnalysis()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void GENSIMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByLabel("genParticles", genParticles);

	//Gen
	int ndstar_p, ndstar_m, nustar_p, nustar_m;
	ndstar_p=ndstar_m=nustar_p=nustar_m=0;

	for (reco::GenParticleCollection::const_iterator genit = genParticles->begin(); genit != genParticles->end();  ++genit) {
		// pythia only status 1
		h_pdgId->Fill(genit->pdgId());
		if( genit->pdgId() ==  4000001 ){ ndstar_p++; }
		if( genit->pdgId() == -4000001 ){ ndstar_m--; }
		if( genit->pdgId() ==  4000002 ){ nustar_p++; }
		if( genit->pdgId() == -4000002 ){ nustar_m--; }

		if(fabs(genit->pdgId())== 4000001 ){
			h_mdstar->Fill(genit->mass());
			h_pTdstar->Fill(genit->pt());
			h_etadstar->Fill(genit->eta());
			h_phidstar->Fill(genit->phi());
			for(unsigned int dit=0; dit<genit->numberOfDaughters();dit++){
				h_dstarDecay->Fill(genit->daughter(dit)->pdgId());				
			}
		}
		if(fabs(genit->pdgId())== 4000002 ){
			h_mustar->Fill(genit->mass());
			h_pTustar->Fill(genit->pt());
			h_etaustar->Fill(genit->eta());
			h_phiustar->Fill(genit->phi());
			for(unsigned int dit=0; dit<genit->numberOfDaughters();dit++){
				h_ustarDecay->Fill(genit->daughter(dit)->pdgId());				
			}
		}
	}
	if( ndstar_p != 0 ) h_ndstar->Fill(ndstar_p);
	if( ndstar_m != 0 ) h_ndstar->Fill(ndstar_m);
	if( nustar_p != 0 ) h_nustar->Fill(nustar_p);
	if( nustar_m != 0 ) h_nustar->Fill(nustar_m);
}


// ------------ method called once each job just before starting event loop  ------------
void GENSIMAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GENSIMAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void GENSIMAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void GENSIMAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void GENSIMAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GENSIMAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GENSIMAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GENSIMAnalysis);
