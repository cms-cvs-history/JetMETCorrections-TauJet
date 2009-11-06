// system include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <iostream>
using namespace std;

#include "JetMETCorrections/TauJet/interface/TCTauCorrector.h"
#include "JetMETCorrections/TauJet/interface/TauJetCorrector.h"

#include "DataFormats/TauReco/interface/CaloTau.h"
#include "RecoTauTag/TauTagTools/interface/CaloTauElementsOperators.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "RecoTauTag/TauTagTools/interface/PFTauElementsOperators.h"

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"

#include "JetMETCorrections/TauJet/test/visibleTaus.h"

class TCTauCorrectorTest : public edm::EDAnalyzer {
  public:
  	explicit TCTauCorrectorTest(const edm::ParameterSet&);
  	~TCTauCorrectorTest();

  	virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  	virtual void beginJob(const edm::EventSetup& );
  	virtual void endJob();
  private:
	bool prongSelection(short int);

	TCTauCorrector* tcTauCorrector;
	TauJetCorrector* tauJetCorrector;

	TH1F* h_CaloTau_dEt;
	TH1F* h_CaloTau_caloTauCorrected_dEt;
	TH1F* h_CaloTau_TCTauCorrected_dEt;
	TH1F* h_CaloTau_doubleCorrected_dEt;
	TH1F* h_CaloTau_jptCorrected_dEt;
	TH1F* h_CaloTau_jptTCTauCorrected_dEt;
	TH1F* h_PFTau_dEt;
	TH1F* h_counters;

	int all,
	    caloTauIn01Counter,
	    caloTauTauJetCorrectedIn01Counter,
	    tcTauIn01Counter,
	    doubleCorrectedIn01Counter,
	    jptTauIn01Counter,
	    jptTCTauCorrectedIn01Counter,
	    pfAll,
	    pfTauIn01Counter;
	int prongs;
	double tauEtCut,tauEtaCut;
};

TCTauCorrectorTest::TCTauCorrectorTest(const edm::ParameterSet& iConfig){
	// this algo (TCTauAlgorithm)
	tcTauCorrector = new TCTauCorrector(iConfig);

	// TauJet correction from JetMETCorrections/TauJet
	tauJetCorrector = new TauJetCorrector(iConfig);

	all                               = 0;
	caloTauIn01Counter	          = 0;
	caloTauTauJetCorrectedIn01Counter = 0;
	tcTauIn01Counter                  = 0;
	doubleCorrectedIn01Counter        = 0;
	jptTauIn01Counter		  = 0;
	jptTCTauCorrectedIn01Counter	  = 0;
	pfAll				  = 0;
	pfTauIn01Counter	          = 0;

	string prongSele = iConfig.getParameter<string>("ProngSelection");
        prongs = -1;
        if(prongSele == "1prong") prongs = 1;
        if(prongSele == "3prong") prongs = 3;

	tauEtCut  = iConfig.getParameter<double>("TauJetEt");
	tauEtaCut = iConfig.getParameter<double>("TauJetEta");
}

TCTauCorrectorTest::~TCTauCorrectorTest(){

	cout << endl << endl;
	cout << "All tau candidates               " << tcTauCorrector->allTauCandidates() << endl;
	cout << "Number of taus passing selection " << tcTauCorrector->statistics() << endl;
	cout << "Algorithm efficiency             " << tcTauCorrector->efficiency() << endl;

	cout << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau                  " << double(caloTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TauJetCorrection " << double(caloTauTauJetCorrectedIn01Counter)/all << endl;
	cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TCTauCorrection  " << double(tcTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TauJet+TCTau     " << double(doubleCorrectedIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+JPTCorrection    " << double(jptTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+JPT+TCTau        " << double(jptTCTauCorrectedIn01Counter)/all << endl;
	cout << "Fraction of jets in abs(dEt) < 0.1, reco::PFTau                    " << double(pfTauIn01Counter)/pfAll << endl;
        cout << endl;

        h_counters->SetBinContent(1,tcTauCorrector->allTauCandidates());
        h_counters->SetBinContent(2,tcTauCorrector->statistics());
        h_counters->SetBinContent(4,double(all));
        h_counters->SetBinContent(5,caloTauIn01Counter);
        h_counters->SetBinContent(6,caloTauTauJetCorrectedIn01Counter);
        h_counters->SetBinContent(7,tcTauIn01Counter);
        h_counters->SetBinContent(8,doubleCorrectedIn01Counter);
        h_counters->SetBinContent(9,jptTauIn01Counter);
        h_counters->SetBinContent(10,jptTCTauCorrectedIn01Counter);
        h_counters->SetBinContent(11,pfAll);
        h_counters->SetBinContent(12,pfTauIn01Counter);


        TFile* oFILE = new TFile("histograms.root","RECREATE");

        oFILE->cd();

        h_CaloTau_jptTCTauCorrected_dEt->Write();
        h_CaloTau_jptCorrected_dEt->Write();
        h_CaloTau_doubleCorrected_dEt->Write();
        h_CaloTau_TCTauCorrected_dEt->Write();
        h_CaloTau_caloTauCorrected_dEt->Write();
        h_CaloTau_dEt->Write();
        h_PFTau_dEt->Write();

        h_counters->Write();

        oFILE->Close();

	delete tcTauCorrector;
}

void TCTauCorrectorTest::beginJob(const edm::EventSetup& iSetup){
	h_CaloTau_dEt                   = new TH1F("h_CaloTau_dEt","",100,-1.0,1.0);
	h_CaloTau_caloTauCorrected_dEt  = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_caloTauCorrected_dEt");
	h_CaloTau_TCTauCorrected_dEt    = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_TCTauCorrected_dEt");
	h_CaloTau_doubleCorrected_dEt   = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_doubleCorrected_dEt");
	h_CaloTau_jptCorrected_dEt      = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_jptCorrected_dEt");
	h_CaloTau_jptTCTauCorrected_dEt = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_jptTCTauCorrected_dEt");
	h_PFTau_dEt                     = (TH1F*)h_CaloTau_dEt->Clone("h_PFTau_dEt");

	h_counters = new TH1F("h_counters","",12,0,12);
}

void TCTauCorrectorTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

        double matchingConeSize         = 0.1,
               signalConeSize           = 0.07,
               isolationConeSize        = 0.4,
               ptLeadingTrackMin        = 6,
               ptOtherTracksMin         = 1;
        string metric = "DR"; // can be DR,angle,area
        unsigned int isolationAnnulus_Tracksmaxn = 0;

	vector<TLorentzVector> mcTaus = ::visibleTaus(iEvent,0);//37);

// CaloTaus

        Handle<CaloTauCollection> theCaloTauHandle;
        try{
          iEvent.getByLabel("caloRecoTauProducer",theCaloTauHandle);
        }catch(...) {;}

	tcTauCorrector->eventSetup(iEvent,iSetup);

        if(theCaloTauHandle.isValid()){

	  const CaloTauCollection & caloTaus = *(theCaloTauHandle.product());
	  //int nCaloTaus = caloTaus.size();
	  //cout << "calotau collection size " << nCaloTaus << endl;

	  // for JPT
	  edm::Handle<CaloJetCollection> zspjets;
	  iEvent.getByLabel(InputTag("ZSPJetCorJetIcone5"), zspjets);

	  string jpt = "JetPlusTrackZSPCorrectorIcone5";
	  const JetCorrector* correctorJPT = JetCorrector::getJetCorrector (jpt, iSetup);
	  //

          CaloTauCollection::const_iterator iTau;
          for(iTau = caloTaus.begin(); iTau != caloTaus.end(); iTau++){
                if(!iTau->leadTrack()) continue;

		if(iTau->pt() < tauEtCut || fabs(iTau->eta()) > tauEtaCut) continue;

		CaloTau theCaloTau = *iTau;
		CaloTauElementsOperators op(theCaloTau);
		double d_trackIsolation = op.discriminatorByIsolTracksN(
                                metric,
                                matchingConeSize,
                                ptLeadingTrackMin,
                                ptOtherTracksMin,
                                metric,
                                signalConeSize,
                                metric,
                                isolationConeSize,
                                isolationAnnulus_Tracksmaxn);
		if(d_trackIsolation == 0) continue;

		const TrackRef leadingTrack =op.leadTk(metric,matchingConeSize,ptLeadingTrackMin);
		if(leadingTrack.isNull()) continue;

		const TrackRefVector signalTracks = op.tracksInCone(leadingTrack->momentum(),metric,signalConeSize,ptOtherTracksMin);
		if(!prongSelection(signalTracks.size())) continue;
 


		double tcTauCorrection = tcTauCorrector->correction(theCaloTau);
		cout << "CaloTau Et = " << iTau->pt() << endl;
		cout << "TCTau corrected Et = " << iTau->pt()*tcTauCorrection << endl;

		double MC_Et = 0;
		vector<TLorentzVector>::const_iterator i;
        	for(i = mcTaus.begin(); i!= mcTaus.end(); ++i){
			XYZVector direction(i->Px(),i->Py(),i->Pz());
                	double DR = ROOT::Math::VectorUtil::DeltaR(direction,iTau->momentum());
                	if(DR < 0.5) MC_Et = i->Pt();
        	}

	        // Using TauJet Jet energy correction AND TCTau correction
	        double tauJetCorrection = tauJetCorrector->correction(iTau->p4());
		CaloTau tauJetCorrected = theCaloTau;
		tauJetCorrected.setP4(iTau->p4()*tauJetCorrection);

		double doubleCorrection = tcTauCorrector->correction(tauJetCorrected);
                cout << "CaloTau+TauJet Et = " << tauJetCorrected.pt() << endl;
		cout << "CaloTau+TauJet+TCTau Et = " << tauJetCorrected.pt()*doubleCorrection << endl;

		// JPT
		CLHEP::HepLorentzVector cjetc(iTau->px(), iTau->py(), iTau->pz(), iTau->energy());
	        CaloJetCollection::const_iterator zspjet;
         	for( zspjet = zspjets->begin();
                     zspjet != zspjets->end(); ++zspjet ){
           		CLHEP::HepLorentzVector zspjetc(zspjet->px(), zspjet->py(), zspjet->pz(), zspjet->energy());
           		double dr = zspjetc.deltaR(cjetc);
           		if(dr < 0.001) break;
         	}

		double jptCorrection = correctorJPT->correction (*zspjet,iEvent,iSetup);
		CaloTau jptCorrected = theCaloTau;
		jptCorrected.setP4(iTau->p4()*jptCorrection);

                double jptTCTauCorrection = tcTauCorrector->correction(jptCorrected);
                cout << "CaloTau+JPT Et = " << jptCorrected.pt() << endl;
		cout << "CaloTau+JPT+TCTau Et = " << jptCorrected.pt()*jptTCTauCorrection << endl;

		if(MC_Et > 0){
			double caloTau_dEt = (iTau->pt() - MC_Et)/MC_Et;
			h_CaloTau_dEt->Fill(caloTau_dEt);

			double caloTau_TauJetCorrected_dEt = (tauJetCorrected.pt() - MC_Et)/MC_Et;
			h_CaloTau_caloTauCorrected_dEt->Fill(caloTau_TauJetCorrected_dEt);

			double caloTau_TCTauCorrected_dEt = (iTau->pt()*tcTauCorrection - MC_Et)/MC_Et;
			h_CaloTau_TCTauCorrected_dEt->Fill(caloTau_TCTauCorrected_dEt);

			double caloTau_DoubleCorrected_dEt = (tauJetCorrected.pt()*doubleCorrection - MC_Et)/MC_Et;
			h_CaloTau_doubleCorrected_dEt->Fill(caloTau_DoubleCorrected_dEt);

			double caloTau_jptCorrected_dEt = (jptCorrected.pt() - MC_Et)/MC_Et;
			h_CaloTau_jptCorrected_dEt->Fill(caloTau_jptCorrected_dEt);

			double caloTau_jptTCTauCorrected_dEt = (jptCorrected.pt()*jptTCTauCorrection - MC_Et)/MC_Et;
			h_CaloTau_jptTCTauCorrected_dEt->Fill(caloTau_jptTCTauCorrected_dEt);

			all++;
			if(fabs(caloTau_dEt) < 0.1) caloTauIn01Counter++;
			if(fabs(caloTau_TauJetCorrected_dEt) < 0.1) caloTauTauJetCorrectedIn01Counter++;
			if(fabs(caloTau_TCTauCorrected_dEt) < 0.1) tcTauIn01Counter++;
			if(fabs(caloTau_DoubleCorrected_dEt) < 0.1) doubleCorrectedIn01Counter++;
			if(fabs(caloTau_jptCorrected_dEt) < 0.1) jptTauIn01Counter++;
			if(fabs(caloTau_jptTCTauCorrected_dEt) < 0.1) jptTCTauCorrectedIn01Counter++;
		}
          }
        }

// PFTau

        Handle<PFTauCollection> thePFTauHandle;
        iEvent.getByLabel("fixedConePFTauProducer",thePFTauHandle);

        if(thePFTauHandle.isValid()){

          const PFTauCollection & pfTaus = *(thePFTauHandle.product());

          //int nPFTaus = pfTaus.size();
          //cout << "pftau collection size " << nPFTaus << endl;

          PFTauCollection::const_iterator iTau;
          for(iTau = pfTaus.begin(); iTau != pfTaus.end(); iTau++){
                if(!iTau->leadTrack()) continue;

                PFTau theCaloTau = *iTau;
                PFTauElementsOperators op(theCaloTau);
                double d_trackIsolation = op.discriminatorByIsolTracksN(
                                metric,
                                matchingConeSize,
                                ptLeadingTrackMin,
                                ptOtherTracksMin,
                                metric,
                                signalConeSize,
                                metric,
                                isolationConeSize,
                                isolationAnnulus_Tracksmaxn);
                if(d_trackIsolation == 0) continue;

                const TrackRef leadingTrack =op.leadTk(metric,matchingConeSize,ptLeadingTrackMin);
                if(leadingTrack.isNull()) continue;

                const TrackRefVector signalTracks = op.tracksInCone(leadingTrack->momentum(),metric,signalConeSize,ptOtherTracksMin);
                if(!prongSelection(signalTracks.size())) continue;

                double MC_Et = 0;
                vector<TLorentzVector>::const_iterator i;
                for(i = mcTaus.begin(); i!= mcTaus.end(); ++i){
                        XYZVector direction(i->Px(),i->Py(),i->Pz());
                        double DR = ROOT::Math::VectorUtil::DeltaR(direction,iTau->momentum());
                        if(DR < 0.5) MC_Et = i->Pt();
                }

                if(MC_Et > 0){
                        double pfTau_dEt = (iTau->pt() - MC_Et)/MC_Et;
                        h_PFTau_dEt->Fill(pfTau_dEt);

                        pfAll++;
                        if(fabs(pfTau_dEt) < 0.1) pfTauIn01Counter++;
                }
          }
        }
}

void TCTauCorrectorTest::endJob(){
}

bool TCTauCorrectorTest::prongSelection(short int nTracks){
        if( nTracks != 1 && nTracks != 3) return false;
        if( (prongs == 1 && nTracks!= 1) || (prongs == 3 && nTracks!= 3)) return false;
	return true;
}

DEFINE_FWK_MODULE(TCTauCorrectorTest);
