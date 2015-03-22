#ifndef QstartToGJ_GENSIMAnalysis_h
#define QstartToGJ_GENSIMAnalysis_h

// -*- C++ -*-
// //// -*- C++ -*-
// //
// // Package:    QstartToGJ 
// // Class:      GENSIMAnalysis
// //

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
class GENSIMAnalysis : public edm::EDAnalyzer {
	public:
		explicit GENSIMAnalysis(const edm::ParameterSet&);
		~GENSIMAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	
		void fillUstarInfo( const reco::GenParticle* gen, int& nustar_p, int& nustar_a, int& nQstar);	
		void fillDstarInfo( const reco::GenParticle* gen, int& ndstar_p, int& ndstar_a, int& nQstar);

		// ----------member data ---------------------------
		edm::InputTag genInfoLabel_;
		int selectQstarStatus_;
		int numEventListsPrint_;

		int numEvts_;

		std::ofstream genLists;
		edm::Service<TFileService> tFileService;

		TH1D* h_numEvt;
		TH1D* h_nQstar;
		TH1D* h_pdgId;
		TH1D* h_ndstar;     
		TH1D* h_nustar;     
		TH1D* h_mdstar;     
		TH1D* h_mustar;     
		TH1D* h_pTdstar;    
		TH1D* h_pTustar;    
		TH1D* h_etadstar;   
		TH1D* h_etaustar;   
		TH1D* h_phidstar;   
		TH1D* h_phiustar;   
		TH1D* h_dstarDecay; 
		TH1D* h_ustarDecay;
		TH1D* h_dstarStatus; 
		TH1D* h_ustarStatus; 

};





#endif
