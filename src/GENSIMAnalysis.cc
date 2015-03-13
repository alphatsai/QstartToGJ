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
	numEvts_=0;
	printEvts_=2;

	genLists.open ("genLists.txt", std::ofstream::out | std::ofstream::app);

	h_numEvt           = tFileService->make<TH1D>("Num_Evt", "",  1, 1, 2);
	h_pdgId            = tFileService->make<TH1D>("PdgId", "",   100, -50, 50);
	h_nQstar           = tFileService->make<TH1D>("Num_Qstar", "", 50, -25, 25);
	h_ndstar           = tFileService->make<TH1D>("Num_Dstar", "", 50, -25, 25);
	h_nustar           = tFileService->make<TH1D>("Num_Ustar", "", 50, -25, 25);
	h_mdstar           = tFileService->make<TH1D>("Dstar_Mass", "", 2000, 0, 2000);
	h_mustar           = tFileService->make<TH1D>("Ustar_Mass", "", 2000, 0, 2000);
	h_pTdstar          = tFileService->make<TH1D>("Dstar_pT", "", 1000, 0, 1000);
	h_pTustar          = tFileService->make<TH1D>("Ustar_pT", "", 1000, 0, 1000);
	h_etadstar         = tFileService->make<TH1D>("Dstar_Eta", "", 100, -5, 5);
	h_etaustar         = tFileService->make<TH1D>("Ustar_Eta", "", 100, -5, 5);
	h_phidstar         = tFileService->make<TH1D>("Dstar_Phi", "", 400, -4, 4);
	h_phiustar         = tFileService->make<TH1D>("Ustar_Phi", "", 400, -4, 4);
	h_dstarDecay       = tFileService->make<TH1D>("Dstar_Decay", "", 100, -50, 50);
	h_ustarDecay       = tFileService->make<TH1D>("Ustar_Decay", "", 100, -50, 50);
	h_dstarStatus      = tFileService->make<TH1D>("Dstar_Status", "", 600, -300, 300);
	h_ustarStatus      = tFileService->make<TH1D>("Ustar_Status", "", 600, -300, 300);
	h_dstarStatus_phi0 = tFileService->make<TH1D>("Dstar_StatusInPhi0", "", 600, -300, 300);
	h_ustarStatus_phi0 = tFileService->make<TH1D>("Ustar_StatusInPhi0", "", 600, -300, 300);
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
	using namespace std;

	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByLabel("genParticles", genParticles);

	//Gen
	int nQstar, ndstar_p, ndstar_a, nustar_p, nustar_a;
	nQstar=ndstar_p=ndstar_a=nustar_p=nustar_a=0;

	if( numEvts_ <= printEvts_ ){ 
		genLists<<"\n=========================================="<<endl;	
		genLists<<"Event: "<<numEvts_<<endl;	
		genLists<<"PdgID\t"
			<<"Status\t"
			<<"Mother 0\t"
			<<"Mother 1\t"
			<<"# of Daughters"
			<<endl;	
	}
	for (reco::GenParticleCollection::const_iterator itGen=genParticles->begin(); itGen!=genParticles->end(); ++itGen){
		// Print out particel decay lists
		if( numEvts_ <= printEvts_ ){
			genLists<<"\n"<<itGen->pdgId()
				<<"\t"<<itGen->status()
				<<"\t"<<itGen->mother(0)->pdgId()
				<<"\t"<<itGen->mother(1)->pdgId()
				<<"\t"<<itGen->numberOfDaughters()
				<<endl;	
			for( unsigned int iDa=0; iDa<itGen->numberOfDaughters(); iDa++ ){
				genLists<<"|"<<endl;
				if( iDa < (itGen->numberOfDaughters()-1)){ 
					genLists<<"|-> "<<itGen->daughter(iDa)->pdgId()
						<<"\t"	<<itGen->daughter(iDa)->status()
						<<endl;
				}else{ 
					genLists<<"`-> "<<itGen->daughter(iDa)->pdgId()
						<<"\t"	<<itGen->daughter(iDa)->status()
						<<endl;
				}
			}
		}

		h_pdgId->Fill(itGen->pdgId());
		if( itGen->pdgId() ==  4000001 ){ ndstar_p++; nQstar++; }
		if( itGen->pdgId() == -4000001 ){ ndstar_a--; nQstar++; }
		if( itGen->pdgId() ==  4000002 ){ nustar_p++; nQstar++; }
		if( itGen->pdgId() == -4000002 ){ nustar_a--; nQstar++; }

		if(fabs(itGen->pdgId())== 4000001 ){
			h_mdstar->Fill(itGen->mass());
			h_pTdstar->Fill(itGen->pt());
			h_etadstar->Fill(itGen->eta());
			h_phidstar->Fill(itGen->phi());
			h_dstarStatus->Fill(itGen->status());
			if( itGen->phi() == 0 ) h_dstarStatus_phi0->Fill(itGen->status());
			for( unsigned int iDa=0; iDa<itGen->numberOfDaughters(); iDa++ ){
				h_dstarDecay->Fill(itGen->daughter(iDa)->pdgId());				
			}
		}
		if(fabs(itGen->pdgId())== 4000002 ){
			h_mustar->Fill(itGen->mass());
			h_pTustar->Fill(itGen->pt());
			h_etaustar->Fill(itGen->eta());
			h_phiustar->Fill(itGen->phi());
			h_ustarStatus->Fill(itGen->status());
			if( itGen->phi() == 0 ) h_ustarStatus_phi0->Fill(itGen->status());
			for( unsigned int iDa=0; iDa<itGen->numberOfDaughters(); iDa++ ){
				h_ustarDecay->Fill(itGen->daughter(iDa)->pdgId());				
			}
		}
	}
	if( ndstar_p != 0 ) h_ndstar->Fill(ndstar_p);
	if( ndstar_a != 0 ) h_ndstar->Fill(ndstar_a);
	if( nustar_p != 0 ) h_nustar->Fill(nustar_p);
	if( nustar_a != 0 ) h_nustar->Fill(nustar_a);
	h_nQstar->Fill(nQstar);	
	h_numEvt->Fill(1);
	numEvts_++;
}


// ------------ method called once each job just before starting event loop  ------------
void GENSIMAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GENSIMAnalysis::endJob() 
{
	genLists.close();	
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
