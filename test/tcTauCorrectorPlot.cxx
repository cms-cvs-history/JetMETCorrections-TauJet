void tcTauCorrectorPlot(){

	TFile* inFILE = TFile::Open("histograms.root");

	TCanvas* tctau = new TCanvas("tctau","",500,500);
	tctau->SetFrameFillColor(0);
	tctau->SetFillColor(0);
	tctau->Divide(3,2);

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

        tctau->cd(6);
        gPad->SetFrameFillColor(0);
        h_PFTau_dEt->Scale(1/h_PFTau_dEt->GetMaximum());
        h_PFTau_dEt->GetXaxis()->SetTitle("(Et(RECO)-Et(MC))/Et(MC)");
        h_PFTau_dEt->Draw();
        tex = new TLatex(-0.90,1.086937,"PFTau");
        tex->SetLineWidth(2);
        tex->Draw();
}
