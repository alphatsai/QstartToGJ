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
#include <cmath>                                                   

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"              
#include "CommonTools/UtilAlgos/interface/TFileService.h"          

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

		// ----------member data ---------------------------
		TH1D *m_Acceptance;
		TH1D *m_Z_pt;
		TH1D *m_h0_pt;
		TH1D *m_h0_mass;
		TH1D *m_h0gg_mass;
		TH1D *m_ele_pt;
		TH1D *m_mu_pt;
		TH1D *m_Z_lep_pt;
		TH1D *m_pdgId;
		TH1D *m_ZmassMu;
		TH1D *m_ZmassEle;
};





#endif
