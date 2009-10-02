// system include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <iostream>
using namespace std;

#include "JetMETCorrections/TauJet/interface/TCTauCorrector.h"
#include "JetMETCorrections/TauJet/interface/TauJetCorrector.h"

#include "DataFormats/TauReco/interface/CaloTau.h"
#include "RecoTauTag/TauTagTools/interface/CaloTauElementsOperators.h"

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

	int all,
	    caloTauIn01Counter,
	    caloTauTauJetCorrectedIn01Counter,
	    tcTauIn01Counter,
	    doubleCorrectedIn01Counter;
	int prongs;
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

	string prongSele = iConfig.getParameter<string>("ProngSelection");
        prongs = -1;
        if(prongSele == "1prong") prongs = 1;
        if(prongSele == "3prong") prongs = 3;

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
        cout << endl;

	delete tcTauCorrector;
}

void TCTauCorrectorTest::beginJob(const edm::EventSetup& iSetup){
	h_CaloTau_dEt                  = new TH1F("h_CaloTau_dEt","",100,-1.0,1.0);
	h_CaloTau_caloTauCorrected_dEt = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_caloTauCorrected_dEt");
	h_CaloTau_TCTauCorrected_dEt = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_TCTauCorrected_dEt");
	h_CaloTau_doubleCorrected_dEt  = (TH1F*)h_CaloTau_dEt->Clone("h_CaloTau_doubleCorrected_dEt");
}

void TCTauCorrectorTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

        Handle<CaloTauCollection> theCaloTauHandle;
        try{
          iEvent.getByLabel("caloRecoTauProducer",theCaloTauHandle);
        }catch(...) {;}

	tcTauCorrector->eventSetup(iEvent,iSetup);

        if(theCaloTauHandle.isValid()){

          const CaloTauCollection & caloTaus = *(theCaloTauHandle.product());

          //int nCaloTaus = caloTaus.size();
          //cout << "calotau collection size " << nCaloTaus << endl;

	  vector<TLorentzVector> mcTaus = ::visibleTaus(iEvent,0);//37);
////	  if(mcTaus.size() == 0) mcTaus = ::visibleTaus(iEvent,23);
//cout << "mcTaus size() " << mcTaus.size() << endl;
          double matchingConeSize         = 0.1,
                 signalConeSize           = 0.07,
                 isolationConeSize        = 0.4,
                 ptLeadingTrackMin        = 20,
                 ptOtherTracksMin         = 1;
          string metric = "DR"; // can be DR,angle,area
          unsigned int isolationAnnulus_Tracksmaxn = 0;

          CaloTauCollection::const_iterator iTau;
          for(iTau = caloTaus.begin(); iTau != caloTaus.end(); iTau++){
                if(!iTau->leadTrack()) continue;

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
		cout << "CaloTau Et = " << iTau->pt() << ", corrected Et = " << iTau->pt()*tcTauCorrection << endl;

		double MC_Et = 0;
		vector<TLorentzVector>::const_iterator i;
        	for(i = mcTaus.begin(); i!= mcTaus.end(); ++i){
			XYZVector direction(i->Px(),i->Py(),i->Pz());
                	double DR = ROOT::Math::VectorUtil::DeltaR(direction,iTau->momentum());
                	if(DR < 0.5) MC_Et = i->Pt();
        	}

	        // Using TauJet Jet energy correction AND TCTau correction
	        double tauJetCorrection = tauJetCorrector->correction(iTau->p4());
		theCaloTau.setP4(iTau->p4()*tauJetCorrection);

		double doubleCorrection = tcTauCorrector->correction(theCaloTau);
                cout << "CaloTau+TauJet Et = " << theCaloTau.pt() << ", doublecorrected Et = " << theCaloTau.pt()*doubleCorrection << endl;

		if(MC_Et > 0){
			double caloTau_dEt = (iTau->pt() - MC_Et)/MC_Et;
			h_CaloTau_dEt->Fill(caloTau_dEt);

			double caloTau_TauJetCorrected_dEt = (iTau->pt()*tauJetCorrection - MC_Et)/MC_Et;
			h_CaloTau_caloTauCorrected_dEt->Fill(caloTau_TauJetCorrected_dEt);

			double caloTau_TCTauCorrected_dEt = (iTau->pt()*tcTauCorrection - MC_Et)/MC_Et;
			h_CaloTau_TCTauCorrected_dEt->Fill(caloTau_TCTauCorrected_dEt);

			double caloTau_DoubleCorrected_dEt = (theCaloTau.pt()*doubleCorrection - MC_Et)/MC_Et;
			h_CaloTau_doubleCorrected_dEt->Fill(caloTau_DoubleCorrected_dEt);

			all++;
			if(fabs(caloTau_dEt) < 0.1) caloTauIn01Counter++;
			if(fabs(caloTau_TauJetCorrected_dEt) < 0.1) caloTauTauJetCorrectedIn01Counter++;
			if(fabs(caloTau_TCTauCorrected_dEt) < 0.1) tcTauIn01Counter++;
			if(fabs(caloTau_DoubleCorrected_dEt) < 0.1) doubleCorrectedIn01Counter++;
		}
          }
        }
}

void TCTauCorrectorTest::endJob(){
	TFile* oFILE = new TFile("histograms.root","RECREATE");
	
	oFILE->cd();

	h_CaloTau_doubleCorrected_dEt->Write();
	h_CaloTau_TCTauCorrected_dEt->Write();
	h_CaloTau_caloTauCorrected_dEt->Write();
	h_CaloTau_dEt->Write();

	oFILE->Close();
/*
	TCanvas* resolution = new TCanvas("resolution","",500,500);
	resolution->SetFillColor(0);

        resolution->cd();

        h_CaloTau_doubleCorrected_dEt->SetFillColor(2);
        h_CaloTau_doubleCorrected_dEt->SetLineWidth(3);
	h_CaloTau_doubleCorrected_dEt->SetStats(0);
	h_CaloTau_doubleCorrected_dEt->GetXaxis()->SetTitle("(E_{T} - E_{T}(MC))/E_{T}(MC)");
        h_CaloTau_doubleCorrected_dEt->GetYaxis()->SetTitle("Arbitrary units");
        h_CaloTau_doubleCorrected_dEt->DrawClone();

        h_CaloTau_TCTauCorrected_dEt->SetFillColor(5);
        h_CaloTau_TCTauCorrected_dEt->SetLineWidth(3);
        h_CaloTau_TCTauCorrected_dEt->DrawClone("same");

        h_CaloTau_doubleCorrected_dEt->SetFillColor(0);
        h_CaloTau_doubleCorrected_dEt->SetLineWidth(1);
        h_CaloTau_doubleCorrected_dEt->DrawClone("same");


        h_CaloTau_caloTauCorrected_dEt->SetLineStyle(3);
        h_CaloTau_caloTauCorrected_dEt->SetLineWidth(3);
        h_CaloTau_caloTauCorrected_dEt->DrawClone("same");


        h_CaloTau_dEt->SetLineStyle(2);
        h_CaloTau_dEt->SetLineWidth(3);
        h_CaloTau_dEt->DrawClone("same");

        double text_y1 = h_CaloTau_doubleCorrected_dEt->GetMaximum();
        TLatex* text1_1 = new TLatex(-0.8,0.8*text_y1,"reco::CaloTau");
	text1_1->SetTextSize(0.03);
        text1_1->Draw();
        TLatex* text1_2 = new TLatex(-0.8,0.7*text_y1,"reco::CaloTau+TauJet");
        text1_2->SetTextSize(0.03);
        text1_2->Draw();
        TLatex* text1_3 = new TLatex(-0.8,0.6*text_y1,"TCTau");
        text1_3->SetTextSize(0.03);
        text1_3->Draw();
        TLatex* text1_4 = new TLatex(-0.8,0.5*text_y1,"TCTau+TauJet");
        text1_4->SetTextSize(0.03);
        text1_4->Draw();

        TLine* line1_1 = new TLine(-0.92,0.82*text_y1,-0.82,0.82*text_y1);
        line1_1->SetLineWidth(3);
        line1_1->SetLineStyle(2);
        line1_1->Draw();

        TLine* line1_2 = new TLine(-0.92,0.72*text_y1,-0.82,0.72*text_y1);
        line1_2->SetLineWidth(3);
        line1_2->SetLineStyle(3);
        line1_2->Draw();


	TLine* line1_3 = new TLine(-0.92,0.62*text_y1,-0.82,0.62*text_y1);
	line1_3->SetLineWidth(3);
   	line1_3->DrawClone();
	line1_3->SetLineColor(5);
        line1_3->SetLineWidth(1);
        line1_3->DrawClone();

        TLine* line1_4 = new TLine(-0.92,0.52*text_y1,-0.82,0.52*text_y1);
        line1_4->SetLineWidth(3);
        line1_4->DrawClone();
        line1_4->SetLineColor(2);
        line1_4->SetLineWidth(1);
        line1_4->DrawClone();


	resolution->Print("resolution.C");
*/
}

bool TCTauCorrectorTest::prongSelection(short int nTracks){
        if( nTracks != 1 && nTracks != 3) return false;
        if( (prongs == 1 && nTracks!= 1) || (prongs == 3 && nTracks!= 3)) return false;
	return true;
}

DEFINE_FWK_MODULE(TCTauCorrectorTest);
