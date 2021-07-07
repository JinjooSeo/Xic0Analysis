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
	const int nSyst = 2; //# of systematic items: kinematic cut, unfolding

	bool IsWeighted = true; //for MC, use weighted MC
	bool INELLgt0 = false; //for data, enable INEL>0 flag

	//Loop over configurations
	//--------------------------------------------------------------------

	TH1D* H1_xSec[nConfig];
	TH1D* H1_Syst[nConfig];
	for (int a=0; a<nConfig; a++)
	{
		const char* TRIG = Config[a][0];
		const char* PERC = Config[a][1];
		const char* Suffix = Form("%s_%s", TRIG, PERC); //Suffix
		cout <<Form("\nProcessing configuration %s...", Suffix) <<endl;
		if (IsWeighted==false) cout <<"WARNING: MC being used w/o weighting factors!\n";
		if (INELLgt0==true) cout <<"WARNING: INEL>0 flag is on!\n";

		TFile* F_data = TFile::Open(Form("out_data_%s_%s.root", TRIG, PERC));
		TFile* F_MC   = TFile::Open(IsWeighted?Form("out_MC_%s_%s.root", TRIG, PERC):"out_MC_raw.root");

		const double NormF = GetNormFactor("./AnalysisResults_data.root", TRIG, PERC, INELLgt0);
		const double V0xSec = 57.8 * 1.E3; //V0AND xSec @ 13 TeV: 57.8 mb +- 2.2 (%)
		const double Xic0SemiLBR = 0.018; //Xic0 -> eXi: 1.8 +- 1.2 (%)

		for (int b=0; b<nSyst+1; b++)
		{
			//Default cuts
			const char* Cut1     = "eRec";
			const char* Cut1Flag = "stand";

			const char* Cut2     = "Bayes";
			const char* Cut2Flag = "stand3";
			const char* UFMethod = "Bayes"; //Bayes or Svd
			const int nIter = 3;

			//+++++++++++++++++++++++++++++++++++++++
		
			//Systematic variations: override default one
			if (b!=0) Suffix = Form("%s_Syst%i", Suffix, b);
			if (b==1) Cut1Flag = "vtight";
			if (b==2) UFMethod = "Svd";

			//+++++++++++++++++++++++++++++++++++++++
		
			TH1D* H1a = (TH1D*)F_data->Get(Form("hRawPt_%s_%s", Cut1, Cut1Flag))->Clone();
			TH1D* H1b = ApplyPreFilterEff(H1a, F_data, Cut1, Cut1Flag, Suffix, Show);
			TH1D* H1c = ApplyUnfolding(H1b, F_MC, Cut2, Cut2Flag, UFMethod, nIter, Suffix, Show);
			TH1D* H1d = ApplyXic0Eff(H1c, F_MC, IsWeighted, Cut2, Cut2Flag, Suffix, Show);
			TH1D* H1e = GetXS(H1d, NormF, V0xSec, Xic0SemiLBR, Suffix, 1, 16);//!strcmp(TRIG, "HMV0")?2:1, 16);

			if (b==0)
			{
				H1_xSec[a] = (TH1D*)H1e->Clone();
				H1_Syst[a] = (TH1D*)H1e->Clone(Form("%s_Syst", H1e->GetName()));
				H1_Syst[a]->Reset();
			}//Default XS
			else
			{
				for (int c=0; c<H1e->GetNbinsX(); c++)
				{
					if (H1_xSec[a]->GetBinContent(c+1) == 0) continue;
					const double SystErrOld = H1_Syst[a]->GetBinContent(c+1);
					const double SystErrNew = fabs(H1_xSec[a]->GetBinContent(c+1) - H1e->GetBinContent(c+1));
					H1_Syst[a]->SetBinContent(c+1, SystErrOld + SystErrNew);

					cout <<Form("Updating... %s_%s, SystErr%i, pT = %4.1f, %6.5f + %6.5f = %6.5f\n",
							TRIG, PERC, b, H1e->GetBinCenter(c+1), SystErrOld, SystErrNew, SystErrOld+SystErrNew);
				}//c, loop over pT bins
			}//Systematic

			H1a->Delete();
			H1b->Delete();
			H1c->Delete();
			H1d->Delete();
			H1e->Delete();
		}//b, check systematic errors

		F_data->Close();
		F_MC  ->Close();
	}//a

	//Draw final xSection
	//--------------------------------------------------------------------

	#if 1
	gStyle->SetOptStat(0);
	const int COLOR[] = {1, 2, 210, 4};
	const int MARKER[] = {20, 24, 24, 32};
	TCanvas* c1 = new TCanvas("xSecFin", "", 800*2, 600*2);	c1->cd();
	TLegend* L1 = new TLegend(0.575, 0.85 - 0.0625*nConfig, 0.85, 0.85); L1->SetMargin(0.3);
	for (int a=0; a<nConfig; a++)
	{
		if (a==0)
		{
			gPad->SetGridy();
			gPad->SetLogy();
			TH1F* H1 = new TH1F("H1_frame", "#Xi_{c}^{0} #rightarrow e#Xi#nu_{e}, pp 13 TeV;", 20,0,20);
			H1->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
			H1->GetYaxis()->SetTitle("d^{2}#sigma/(dp_{T}dy) (#mub GeV^{-1}#it{c})");
			H1->GetXaxis()->SetTitleOffset(1.2);
			H1->GetYaxis()->SetTitleOffset(1.2);
			H1->GetXaxis()->SetRangeUser(0.0, 16.5);
			H1->GetYaxis()->SetRangeUser(2.e-3, 2.e3);
			H1->DrawCopy("9");
		}

		H1_xSec[a]->SetLineColor(COLOR[a]);
		H1_xSec[a]->SetMarkerColor(COLOR[a]);
		H1_xSec[a]->SetMarkerStyle(MARKER[a]);
		H1_xSec[a]->SetMarkerSize(1.65);
		H1_xSec[a]->DrawCopy("pe 9 same");
		H1_xSec[a]->SetName(Form("Xic0_xSec_%s_%s", Config[a][0], Config[a][1]));

		for (int b=0; b<H1_xSec[a]->GetNbinsX(); b++)
		{
			if (H1_xSec[a]->GetBinContent(b+1) == 0) continue;
			if (H1_Syst[a]->GetBinContent(b+1) == 0) continue;
			const float x1 = H1_xSec[a]->GetBinCenter(b+1) - 0.25;
			const float x2 = H1_xSec[a]->GetBinCenter(b+1) + 0.25;
			const float y1 = H1_xSec[a]->GetBinContent(b+1) - H1_Syst[a]->GetBinContent(b+1);
			const float y2 = H1_xSec[a]->GetBinContent(b+1) + H1_Syst[a]->GetBinContent(b+1);
			TBox* SystBox = new TBox(x1, y1, x2, y2);
			SystBox->SetFillStyle(0);
			SystBox->SetLineColor(COLOR[a]);
			SystBox->Draw("9 same");
		}

		TString CONFIG = Form("%s, [%s]", Config[a][0], Config[a][1]);
		CONFIG.ReplaceAll("to", ", ");
		CONFIG.ReplaceAll("p1", ".1");
		L1->AddEntry(H1_xSec[a], (const char*)CONFIG, "lp");
		if (a==nConfig-1) L1->Draw("9 same");
	}
	c1->Print("c1z_xSecFin.png");

	//TFile* Fout = new TFile("out_xSec.root", "recreate");
	//for (int a=0; a<nConfig; a++) H1_xSec[a]->Write();
	//c1->Write();
	//Fout->Close();
	#endif

	return;
}//Main
