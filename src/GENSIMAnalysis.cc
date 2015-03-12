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
	edm::Service<TFileService> tFileService;
	m_Acceptance          = tFileService->make<TH1F>("Acceptance", "Acceptance;;number", 10, 0, 10);
	m_h0_pt               = tFileService->make<TH1F>("h0_pt", "GenHiggs_pt;pt [GeV/c];Events", 200, 0., 800.);
	m_h0_mass             = tFileService->make<TH1F>("h0_mass", "GenHiggs_mass; mass [GeV/c^{2}];Events", 20, 0., 200.);  
	m_h0gg_mass           = tFileService->make<TH1F>("h0gg_mass", "gg->h;mass [GeV/c^{2}];Events", 200, 0., 500.);
	m_Z_lep_pt            = tFileService->make<TH1F>("m_Z_lep_pt", "GenZ_lep_pt;pt [GeV/c];Events", 1000, 0., 1000.);
	m_Z_pt                = tFileService->make<TH1F>("Z_pt", "GenZ_pt;pt [GeV/c];Events", 1000, 0., 1000.);
	m_ele_pt              = tFileService->make<TH1F>("ele_pt", "Genele_pt;pt [GeV/c];Events", 1000, 0., 1000.);
	m_mu_pt               = tFileService->make<TH1F>("mu_pt", "Genmu_pt;pt [GeV/c];Events", 1000, 0., 1000.);
	m_pdgId               = tFileService->make<TH1F>("pdgId", "Gen_pdgId;pdgId;number", 100, -50, 50);
	m_ZmassEle            = tFileService->make<TH1F>("ZMassEle", "ZMassEle;mass [GeV/c^{2}];Events", 100, 80., 100.);
	m_ZmassMu             = tFileService->make<TH1F>("ZMassMu", "ZMassMu;mass [GeV/c^{2}];Events", 100, 80., 100.);
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

	// load collection from the Event
	//Handle<reco::GenParticleCollection> genParticles;
	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByLabel("genParticles", genParticles);

	//Gen
	for (reco::GenParticleCollection::const_iterator genit = genParticles->begin();
			genit != genParticles->end();  ++genit) {

		if(genit->status()==3)                
			m_pdgId->Fill(genit->pdgId());

		//Loop over Z's daughter
		if(fabs(genit->pdgId())==23 && genit->status()==3 ){

			for(unsigned int dit=0;dit<genit->numberOfDaughters();dit++){
				if(fabs(genit->daughter(dit)->pdgId())==11 || fabs(genit->daughter(dit)->pdgId())==13){
					m_Z_lep_pt->Fill(genit->daughter(dit)->pt());
				}
			}
		}
		//Gen Electrons
		if(fabs(genit->pdgId())==11 && genit->status()==3){
			m_ele_pt->Fill(genit->pt());
			m_Acceptance->Fill(0);
			if(genit->pt()>20. && fabs(genit->eta())<2.5){
				m_Acceptance->Fill(1);
			}
		}

		//Gen Muons
		if(fabs(genit->pdgId())==13 && genit->status()==3){
			m_mu_pt->Fill(genit->pt());
			m_Acceptance->Fill(4);
			if(genit->pt()>20. && fabs(genit->eta())<2.4){
				m_Acceptance->Fill(5);
			}
		}

		//Gen Higgs
		if(fabs(genit->pdgId())==25 && genit->status()==3){
			m_h0_pt->Fill(genit->pt());
			m_h0_mass->Fill(genit->mass());
		}
	}

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
