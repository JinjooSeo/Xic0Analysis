#include "Xi0cAnaFunction.C"

void Xi0cAnaExecute(bool Show = false)
{
	const char* Config[][2] =
	{
		{"MB", "0.0to100.0"},
		{"MB", "0.1to30.0"},
		{"MB", "30.0to100.0"},
		//{"HMV0", "0.0to0.1"},
		//{"HMSPD", "0.0to0.1"},
		{"", ""}
	};
	const int nConfig = sizeof(Config)/sizeof(Config[0]) - 1;

	//Preparations
	//-------------------------------------------
	
	vector<float> pTvec = {1, 2, 3, 4, 5, 6, 8, 12, 16, 20};

	TFile* F_MCraw = TFile::Open("out_MCraw_MB_0.0to100.0.root");
	TFile* F_MCwgt = TFile::Open("out_MCwgt_MB_0.0to100.0.root");
	TFile* F_CMSLb = TFile::Open("HEPData-ins1113442-v1-Table_2.root");
	TH1D* H1_CMSLb = GetCMSLbSpectrum(F_CMSLb, Show);

	const char* Cut1     = "eRec";
	const char* CutFlag1 = "stand";

	const char* Cut2     = "Bayes";
	const char* CutFlag2 = "stand3";
	const char* UFMethod = "Bayes"; //Bayes or Svd

	const int nIter = 3; //for unfolding
	bool IsWeighted = true;

	//Loop over configurations
	//-------------------------------------------

	TH1D* H1_xSec[nConfig];
	for (int a=0; a<nConfig; a++)
	{
		const char* TRIG = Config[a][0];
		const char* PERC = Config[a][1];

		TString SuffixOrig = Form("%s_%s", TRIG, PERC);
		SuffixOrig.ReplaceAll(".0", "");
		const char* Suffix = SuffixOrig;
		cout <<Form("\nProcessing configuration %s...", Suffix) <<endl;

		TFile* F_data = TFile::Open(Form("out_data_%s_%s.root", TRIG, PERC));
		TFile* F_MC   = IsWeighted?F_MCwgt:F_MCraw;

		TH1D* H1_preFCorr = GetPreFilterCorrSpectrum(F_data, Cut1, CutFlag1, Suffix, Show);
		TH1D* H1_XibCorr  = GetBottomCorrSpectrum(F_MC, H1_CMSLb, H1_preFCorr, pTvec, Cut1, CutFlag1, Suffix, Show);
		TH1D* H1_unfolded = GetUnfoldedSpectrum(F_MC,H1_XibCorr,pTvec, Cut2,CutFlag2,UFMethod,nIter, Suffix,Show);
		TH1D* H1_eff      = GetEfficiency(F_MC, pTvec, Cut2, CutFlag2, IsWeighted, Suffix, Show);

		const double normF = GetNormFactor(F_data, TRIG, PERC, Show);
		H1_xSec[a] = GetXSection(H1_unfolded, H1_eff, pTvec, normF, Suffix, Show);
		
		cout <<H1_preFCorr->GetName() <<endl;
		cout <<H1_XibCorr ->GetName() <<endl;
		cout <<H1_unfolded->GetName() <<endl;
		cout <<H1_eff     ->GetName() <<endl;
		cout <<H1_xSec[a] ->GetName() <<endl;
	}//a

	//Draw final xSection
	//-------------------------------------------

	#if 1
	gStyle->SetOptStat(0);
	TCanvas* c1 = new TCanvas("xSecFin", "", 800, 600);	c1->cd();
	TLegend* L1 = new TLegend(0.6, 0.6, 0.85, 0.85);
	L1->SetMargin(0.3);
	for (int a=0; a<nConfig; a++)
	{
		if (a==0)
		{
			gPad->SetLogy();
			gPad->SetGrid();
			TH1F* H1 = new TH1F("H1_frame", "#Xi_{c}^{0} cross section;p_{T} (GeV);", 20,0,20);
			H1->GetXaxis()->SetTitleOffset(1.2);
			H1->GetXaxis()->SetRangeUser(0.0, 16.5);
			H1->GetYaxis()->SetRangeUser(3.e-3, 5.e2);
			H1->Draw();
		}

		const int COLOR = a!=2?(a+1):(a+2);
		H1_xSec[a]->SetLineColor(COLOR);
		H1_xSec[a]->SetMarkerColor(COLOR);
		H1_xSec[a]->SetMarkerStyle(24);
		H1_xSec[a]->DrawCopy("pe same");

		TString CONFIG = Form("%s, [%s]", Config[a][0], Config[a][1]);
		CONFIG.ReplaceAll("to", ", ");
		CONFIG.ReplaceAll(".0", "");
		L1->AddEntry(H1_xSec[a], (const char*)CONFIG, "lp");
		if (a==nConfig-1) L1->Draw("same");
	}
	c1->Print("c1z_xSecFin.png");
	#endif

	return;
}//Main
