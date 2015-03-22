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
	std::cout<<"GENSIMAnalysis constructor..."<<std::endl;
	genInfoLabel_ = iConfig.getParameter<edm::InputTag>("genInfoLabel");	
	selectQstarStatus_ = iConfig.getParameter<int>("selectQstarStatus");	
	numEventListsPrint_ = iConfig.getParameter<int>("numEventListsPrint");
}

void GENSIMAnalysis::beginJob()
{
	std::cout<<"GENSIMAnalysis beginJob()..."<<std::endl;
	numEvts_=0;
	genLists.open ("genLists.txt", std::ofstream::out );

	h_numEvt           = tFileService->make<TH1D>("Num_Evt", "",  1, 1, 2);
	h_pdgId            = tFileService->make<TH1D>("PdgId", "",   100, -50, 50);
	h_nQstar           = tFileService->make<TH1D>("Num_Qstar", "", 50,   0, 50);
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
}

// ------------ My functions  ------------
void GENSIMAnalysis::fillUstarInfo( const reco::GenParticle* gen, int& nustar_p, int& nustar_a, int& nQstar )
{
	if(fabs(gen->pdgId())== 4000002 ){
		if( gen->pdgId() == 4000002 ) nustar_p++;
		else nustar_a--;
		h_mustar->Fill(gen->mass());
		h_pTustar->Fill(gen->pt());
		h_etaustar->Fill(gen->eta());
		h_phiustar->Fill(gen->phi());
		h_ustarStatus->Fill(gen->status());
		for( unsigned int iDa=0; iDa<gen->numberOfDaughters(); iDa++ ){
			h_ustarDecay->Fill(gen->daughter(iDa)->pdgId());				
		}
	}
	nQstar++;
}
void GENSIMAnalysis::fillDstarInfo( const reco::GenParticle* gen, int& ndstar_p, int& ndstar_a, int& nQstar )
{
	if(fabs(gen->pdgId())== 4000001 ){
		if( gen->pdgId() == 4000001 ) ndstar_p++;
		else ndstar_a--;
		h_mdstar->Fill(gen->mass());
		h_pTdstar->Fill(gen->pt());
		h_etadstar->Fill(gen->eta());
		h_phidstar->Fill(gen->phi());
		h_dstarStatus->Fill(gen->status());
		for( unsigned int iDa=0; iDa<gen->numberOfDaughters(); iDa++ ){
			h_dstarDecay->Fill(gen->daughter(iDa)->pdgId());				
		}
	}
	nQstar++;
}
// ------------ method called for each event  ------------
void GENSIMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	std::cout<<"GENSIMAnalysis analyze() in Event: "<<numEvts_<<std::endl;
	using namespace edm;
	using namespace std;

	//Handle objects
	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByLabel(genInfoLabel_, genParticles);

	int nQstar, ndstar_p, ndstar_a, nustar_p, nustar_a;
	nQstar=ndstar_p=ndstar_a=nustar_p=nustar_a=0;

	if( numEvts_ <= numEventListsPrint_ ){ 
		genLists<<"\n=========================================="<<endl;	
		genLists<<"Event: "<<numEvts_<<endl;	
	}
	for (reco::GenParticleCollection::const_iterator itGen=genParticles->begin(); itGen!=genParticles->end(); ++itGen){
		//** Print out particel decay lists
		if( numEvts_ <= numEventListsPrint_ ){
			genLists<<"\n\t\t"
				<<"PdgID\t"
				<<"Status\t"
				<<"Mother0\t"
				<<"Mother1\t"
				<<"#Daughters"
				<<endl;
			//PdgID	
			genLists<<"Mother\t\t"<<itGen->pdgId();
			//Status	
			genLists<<"\t"<<itGen->status();
			//Mother0 and Mother1
			if( itGen->numberOfMothers()==1 ){ 
				genLists<<"\t"<<itGen->mother(0)->pdgId();
				genLists<<"\tNULL";
			}else if( itGen->numberOfMothers()>1 ){
				genLists<<"\t"<<itGen->mother(0)->pdgId();
				genLists<<"\t"<<itGen->mother(1)->pdgId();
			}else{
				genLists<<"\tNULL";
				genLists<<"\tNULL";
			}
			//#Daughters
			genLists<<"\t"<<itGen->numberOfDaughters();
			genLists<<endl;
			// Daughter paticles	
			for( unsigned int iDa=0; iDa<itGen->numberOfDaughters(); iDa++ ){
				if( iDa < (itGen->numberOfDaughters()-1)){ 
					genLists<<"Daughter"<<iDa<<"\t"
						<<"|-> "<<itGen->daughter(iDa)->pdgId()
						<<"\t"	<<itGen->daughter(iDa)->status()
						<<endl;
				}else{ 
					genLists<<"Daughter"<<iDa<<"\t"
						<<"`-> "<<itGen->daughter(iDa)->pdgId()
						<<"\t"	<<itGen->daughter(iDa)->status()
						<<endl;
				}
			}
		}
		//** Fill histogram
		h_pdgId->Fill(itGen->pdgId());

		if( selectQstarStatus_ == -1 ){
			fillUstarInfo( &(*itGen), nustar_p, nustar_a, nQstar );
			fillDstarInfo( &(*itGen), ndstar_p, ndstar_a, nQstar );
		}else if( itGen->status() == selectQstarStatus_ ){
			fillUstarInfo( &(*itGen), nustar_p, nustar_a, nQstar );
			fillDstarInfo( &(*itGen), ndstar_p, ndstar_a, nQstar );
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

GENSIMAnalysis::~GENSIMAnalysis()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
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
