#include "Xi0cAnaFunction.C"

void Xi0cAnaExecute(bool Show = false)
{
	const char* Config[][2] =
	{
		{"MB", "0to100"},
		{"MB", "0p1to30"},
		{"MB", "30to100"},
		{"HMV0", "0to0p1"},
		{"", ""}
	};
	const int nConfig = sizeof(Config)/sizeof(Config[0]) - 1;

	//Preparations
	//-------------------------------------------

	vector<float> pTvec = {0, 1, 2, 3, 4, 5, 6, 8, 12, 16, 20};

	const char* Cut1     = "eRec";
	const char* CutFlag1 = "stand";
	const char* Cut2     = "Bayes";
	const char* CutFlag2 = "stand3";
	const char* UFMethod = "Bayes"; //Bayes or Svd

	const int nIter = 3; //for unfolding
	bool IsWeighted = true; //for MC

	TFile* F_MC = TFile::Open(Form("out_MC_%s.root", IsWeighted?"wgt":"raw"));

	//Loop over configurations
	//-------------------------------------------

	TH1D* H1_xSec[nConfig];
	for (int a=0; a<nConfig; a++)
	{
		const char* TRIG = Config[a][0];
		const char* PERC = Config[a][1];
		const char* Suffix = Form("%s_%s", TRIG, PERC);
		cout <<Form("\nProcessing configuration %s...", Suffix) <<endl;

		TFile* F_data = TFile::Open(Form("out_data_%s_%s.root", TRIG, PERC));
		const double normF = GetNormFactorFromANC(TRIG, PERC); //Require original train output

		TH1D* H1_preFCorr = GetPreFilterCorrSpectrum(F_data, Cut1, CutFlag1, Suffix, Show);
		TH1D* H1_unfolded = GetUnfoldedSpectrum(F_MC,H1_preFCorr,pTvec, Cut2,CutFlag2,UFMethod,nIter, Suffix,Show);
		TH1D* H1_eff      = GetEfficiency(F_MC, pTvec, Cut2, CutFlag2, IsWeighted, Suffix, Show);

		H1_xSec[a] = GetXSection(H1_unfolded, H1_eff, pTvec, normF, Suffix, Show);

		/*
		//May 19, 2021: require overhaul from this point - binnings, etc modified quite a lot
		TH1D* H1_unfolded;
		if ( !strcmp(Config[a][0], "MB") && !strcmp(Config[a][1], "0.0to100.0") )
		{
			TFile* F_CMSLb = TFile::Open("HEPData-ins1113442-v1-Table_2.root");
			TH1D* H1_CMSLb = GetCMSLbSpectrum(F_CMSLb, Show);
			TH1D* H1_XibCorr = GetBottomCorrSpectrum(F_MC,H1_CMSLb,H1_preFCorr,pTvec, Cut1,CutFlag1, Suffix,Show);
			H1_unfolded = GetUnfoldedSpectrum(F_MC,H1_XibCorr,pTvec, Cut2,CutFlag2,UFMethod,nIter, Suffix,Show);
		}
		else H1_unfolded = GetUnfoldedSpectrum(F_MC,H1_preFCorr,pTvec, Cut2,CutFlag2,UFMethod,nIter, Suffix,Show);

		const double normF = GetNormFactor(F_data, TRIG, PERC, true);//, Show);
		TH1D* H1_eff = GetEfficiency(F_MC, pTvec, Cut2, CutFlag2, IsWeighted, Suffix, Show);
		H1_xSec[a] = GetXSection(H1_unfolded, H1_eff, pTvec, normF, Suffix, Show);
		*/
	}//a

	//Draw final xSection
	//-------------------------------------------

	#if 1
	gStyle->SetOptStat(0);
	const int COLOR[] = {1, 2, 210, 4, 6};
	const int MARKER[] = {20, 24, 24, 25, 26};
	TCanvas* c1 = new TCanvas("xSecFin", "", 800*2, 600*2);	c1->cd();
	TLegend* L1 = new TLegend(0.575, 0.85 - 0.0625*nConfig, 0.85, 0.85); L1->SetMargin(0.3);
	//TFile* Fout = new TFile("out_xSec.root", "recreate");
	for (int a=0; a<nConfig; a++)
	{
		if (a==0)
		{
			gPad->SetGridy();
			gPad->SetLogy();
			TH1F* H1 = new TH1F("H1_frame", "#Xi_{c}^{0} #rightarrow e#Xi#nu_{e}, p + p @ 13 TeV;", 20,0,20);
			H1->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
			H1->GetYaxis()->SetTitle("d^{2}#sigma/(dp_{T}dy) (#mub GeV^{-1}#it{c})");
			H1->GetXaxis()->SetTitleOffset(1.2);
			H1->GetYaxis()->SetTitleOffset(1.2);
			H1->GetXaxis()->SetRangeUser(0.0, 16.5);
			H1->GetYaxis()->SetRangeUser(2.e-3, 2.e3);
			H1->Draw();
		}

		H1_xSec[a]->SetLineColor(COLOR[a]);
		H1_xSec[a]->SetMarkerColor(COLOR[a]);
		H1_xSec[a]->SetMarkerStyle(MARKER[a]);
		H1_xSec[a]->SetMarkerSize(1.50);
		H1_xSec[a]->DrawCopy("pe same");

		TString outName = Form("Xic0_xSec_%s_%s", Config[a][0], Config[a][1]);
		H1_xSec[a]->SetName(outName);
		//H1_xSec[a]->Write();

		TString CONFIG = Form("%s, [%s]", Config[a][0], Config[a][1]);
		CONFIG.ReplaceAll("to", ", ");
		CONFIG.ReplaceAll("p1", ".1");
		L1->AddEntry(H1_xSec[a], (const char*)CONFIG, "lp");
		if (a==nConfig-1) L1->Draw("same");
	}
	c1->Print("c1z_xSecFin.png");
	//c1->Write();
	#endif

	return;
}//Main
