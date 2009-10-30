void tcTauCorrectorPlot(){

	TFile* inFILE = TFile::Open("histograms.root");

	TCanvas* tctau = new TCanvas("tctau","",700,500);
	tctau->SetFrameFillColor(0);
	tctau->SetFillColor(0);
	tctau->Divide(4,2);

	tctau->cd(1);
	gPad->SetFrameFillColor(0);
	h_CaloTau_dEt->Scale(1/h_CaloTau_dEt->GetMaximum());
	h_CaloTau_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_dEt->Draw();
	tex = new TLatex(-0.28,1.086937,"CaloTau");
   	tex->SetLineWidth(2);
   	tex->Draw();

        tctau->cd(2);
        gPad->SetFrameFillColor(0);
	h_CaloTau_caloTauCorrected_dEt->Scale(1/h_CaloTau_caloTauCorrected_dEt->GetMaximum());
	h_CaloTau_caloTauCorrected_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_caloTauCorrected_dEt->Draw();
        tex = new TLatex(-0.55,1.086937,"CaloTau+TauJet corr.");
        tex->SetLineWidth(2);
        tex->Draw();

        tctau->cd(3);
        gPad->SetFrameFillColor(0);
	h_CaloTau_TCTauCorrected_dEt->Scale(1/h_CaloTau_TCTauCorrected_dEt->GetMaximum());
	h_CaloTau_TCTauCorrected_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_TCTauCorrected_dEt->Draw();
        tex = new TLatex(-0.55,1.086937,"CaloTau+TCTau corr.");
        tex->SetLineWidth(2);
        tex->Draw();

        tctau->cd(4);
        gPad->SetFrameFillColor(0);
	h_CaloTau_doubleCorrected_dEt->Scale(1/h_CaloTau_doubleCorrected_dEt->GetMaximum());
	h_CaloTau_doubleCorrected_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_doubleCorrected_dEt->Draw();
        tex = new TLatex(-0.90,1.086937,"CaloTau+TauJet+TCTau corr");
        tex->SetLineWidth(2);
        tex->Draw();

        tctau->cd(5);
        gPad->SetFrameFillColor(0);
	h_CaloTau_jptCorrected_dEt->Scale(1/h_CaloTau_jptCorrected_dEt->GetMaximum());
	h_CaloTau_jptCorrected_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_jptCorrected_dEt->Draw();
        tex = new TLatex(-0.90,1.086937,"CaloTau+JPT corr");
        tex->SetLineWidth(2);
        tex->Draw();

        tctau->cd(6);
        gPad->SetFrameFillColor(0);
	h_CaloTau_jptTCTauCorrected_dEt->Scale(1/h_CaloTau_jptTCTauCorrected_dEt->GetMaximum());
	h_CaloTau_jptTCTauCorrected_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
	h_CaloTau_jptTCTauCorrected_dEt->Draw();
	tex = new TLatex(-0.90,1.086937,"CaloTau+JPT+TCTau corr");
        tex->SetLineWidth(2);
        tex->Draw();

        tctau->cd(8);
        gPad->SetFrameFillColor(0);
        h_PFTau_dEt->Scale(1/h_PFTau_dEt->GetMaximum());
        h_PFTau_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
        h_PFTau_dEt->Draw();
        tex = new TLatex(-0.90,1.086937,"PFTau");
        tex->SetLineWidth(2);
        tex->Draw();

	//

	int allTauCandidates 				= h_counters->GetBinContent(1),
		statistics				= h_counters->GetBinContent(2),
		all					= h_counters->GetBinContent(4),
		caloTauIn01Counter			= h_counters->GetBinContent(5),
		caloTauTauJetCorrectedIn01Counter	= h_counters->GetBinContent(6),
		tcTauIn01Counter			= h_counters->GetBinContent(7),
		doubleCorrectedIn01Counter		= h_counters->GetBinContent(8),
		jptTauIn01Counter			= h_counters->GetBinContent(9),
		jptTCTauCorrectedIn01Counter		= h_counters->GetBinContent(10),
		pfAll					= h_counters->GetBinContent(11),
		pfTauIn01Counter			= h_counters->GetBinContent(12);

        cout << endl << endl;
        cout << "All tau candidates               " << allTauCandidates << endl;
        cout << "Number of taus passing selection " << statistics << endl;
        double efficiency = 0;
	if(allTauCandidates > 0) efficiency = double(statistics)/allTauCandidates;
        cout << "Algorithm efficiency             " << efficiency << endl;

        cout << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau                  " << double(caloTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TauJetCorrection " << double(caloTauTauJetCorrectedIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TCTauCorrection  " << double(tcTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+TauJet+TCTau     " << double(doubleCorrectedIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+JPTCorrection    " << double(jptTauIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::CaloTau+JPT+TCTau        " << double(jptTCTauCorrectedIn01Counter)/all << endl;
        cout << "Fraction of jets in abs(dEt) < 0.1, reco::PFTau                    " << double(pfTauIn01Counter)/pfAll << endl;
        cout << endl;
}
