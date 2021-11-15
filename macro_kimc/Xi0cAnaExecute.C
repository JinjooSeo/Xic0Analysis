#include "Xi0cAnaFunction.C"

void Xi0cAnaExecute(bool Show = false)
{
	const char* Config[][2] =
	{
		{"MB", "0to100"},
		//{"MB", "0p1to30"},
		//{"MB", "30to100"},
		//{"HMV0", "0to0p1"},
		{"", ""}
	};
	const int nConfig = sizeof(Config)/sizeof(Config[0]) - 1;
	const int nSystErr = 0; //# of systematic items

	const double BRFrac = (3.9 * 1.E-3)/(5.8 * 1.E-4); //Xib fraction (?)
	const double BRXic0eXi = 0.018; //Xic0 -> eXi: 1.8 +- 1.2 (%)

	const bool Weighted = true; //for MC, use weighted MC
	const bool INELLgt0 = false; //for data, enable INEL>0 flag

	//Loop over configurations
	//--------------------------------------------------------------------

	TH1D* H1_xSec[nConfig];
	TH1D* H1_syst[nConfig][nSystErr+1]; 

	for (int a=0; a<nConfig; a++)
	{
		//Setup

		const char* INEL = INELLgt0?"INEL0":"";
		const char* TRIG = Config[a][0];
		const char* PERC = Config[a][1];
		const char* Suffix = Form("%s_%s_%s", INEL, TRIG, PERC);
		cout <<"----------------------------------------------\n";
		cout <<Form("\nProcessing configuration %s...\n", Suffix);
		if (Weighted==false) cout <<"WARNING: unweighted MC being used!\n";
		if (INELLgt0==true) cout <<"WARNING: INEL>0 flag is on!\n";

		TFile* F_data = TFile::Open(Form("out_data%s.root", Suffix));
		TFile* F_MC   = TFile::Open(Weighted?Form("out_MC%s.root", Suffix):"out_MC_raw.root");

		//Get normalization factor and year weighted V0xSec
		const double NormF = GetNormFac("AnalysisResults_data.root", Form("ANC%s", Suffix), true);
		const double Frac2016 = GetNormFac("AnalysisResults_data_2016.root", Form("ANC%s", Suffix))/NormF;
		const double Frac2017 = GetNormFac("AnalysisResults_data_2017.root", Form("ANC%s", Suffix))/NormF;
		const double Frac2018 = GetNormFac("AnalysisResults_data_2018.root", Form("ANC%s", Suffix))/NormF;

		const double xSec2016 = GetV0xSec(2016, "pp");
		const double xSec2017 = GetV0xSec(2017, "pp");
		const double xSec2018 = GetV0xSec(2018, "pp");
		const double xSecWgt = Frac2016*xSec2016 + Frac2017*xSec2017 + Frac2018*xSec2018;
		const double V0xSec = xSecWgt * 1.E3; // multiplying 1.E3 - to match the scale? (from legacy code)
		cout <<Form("Frac: %4.3f (2016), %4.3f (2017), %4.3f (2018) / Weighted V0 xSec: %5.3f\n\n",
				Frac2016, Frac2017, Frac2018, xSecWgt);

		//Analysis
		for (int b=0; b<nSystErr+1; b++)
		{
			//Prefilter
			const char* Cut1     = "eRec";//"Svd";
			const char* Cut1Flag = "stand";//"stand3";

			//Unfolding
			const char* Cut2     = "Bayes";//"Svd";
			const char* Cut2Flag = "stand3";
			const char* UFMethod = "Bayes"; //Bayes or Svd
			const int nIter = 3;

			//Xic0 eff
			const char* Cut3     = "eRec";//"Svd";
			const char* Cut3Flag = "stand";//"stand3";

			//Systematic error study: overrides existing cut variable
			if (b!=0) { Suffix = Form("%s_%s_Syst%i", TRIG, PERC, b); }
			if (b==1) { Cut1Flag = "vtight"; }
			if (b==2) { UFMethod = "Svd"; }

			//+++++++++++++++++++++++++++++++++++++++

			TH1D* H1a = (TH1D*)F_data->Get(Form("hRawPt_%s_%s", Cut1, Cut1Flag));
			TH1D* H1b = ApplyPreFilterEff(H1a, F_data, Cut1, Cut1Flag, Suffix, Show);
			ApplyBottomCorr(H1b, F_MC, BRFrac, NormF, V0xSec, Cut1, Cut1Flag, Show);

			TH1D* H1c = ApplyUnfolding(H1b, F_MC, Cut2, Cut2Flag, UFMethod, nIter, Suffix, Show);
			TH1D* H1d = ApplyXic0Eff(H1c, F_MC, Weighted, Cut3, Cut3Flag, Suffix, Show);
			TH1D* H1e = GetXS(H1d, NormF, V0xSec, BRXic0eXi, Suffix);
			//ApplyPromptFrac(H1e, F_MC, BRXic0eXi);//, true); return; //Show);

			TH1D* H1f = Cleanup(H1e);//, !(strcmp(TRIG, "MB"))?1.:2., 16.);

			//Special: extract pT weighting factors - require debugging (Nov. 15)
			//if (Weighted==false) GetWeightingFactor(H1f, Suffix);//, 1., 16.);

			//+++++++++++++++++++++++++++++++++++++++

			if (b==0) //Xsec
			{
				H1_xSec[a]    = (TH1D*)H1f->Clone();
				H1_syst[a][b] = (TH1D*)H1f->Clone(Form("%s_SystA", H1f->GetName()));
				H1_syst[a][b]->Reset();
			}
			else //Systematic
			{
				H1_syst[a][b] = (TH1D*)H1f->Clone();
				H1_syst[a][b]->Reset();
				for (int c=0; c<H1f->GetNbinsX(); c++)
				{
					if (H1_xSec[a]->GetBinContent(c+1) == 0) continue;
					const double systErr = fabs(H1_xSec[a]->GetBinContent(c+1) - H1f->GetBinContent(c+1));
					H1_syst[a][b]->SetBinContent(c+1, systErr);
				}//c
				H1_syst[a][0]->Add(H1_syst[a][b]);
			}

			//++++++++++++++++++++++++++++++++++++++++

			H1a->Delete();
			H1b->Delete();
			H1c->Delete();
			H1d->Delete();
			H1e->Delete();
			H1f->Delete();
			cout <<endl;
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
			//H1->GetXaxis()->SetRangeUser(0.0, 16.5);
			H1->GetYaxis()->SetRangeUser(2.E-3, 4.E3);
			H1->DrawCopy("9");
		}

		H1_xSec[a]->SetLineColor(COLOR[a]);
		H1_xSec[a]->SetMarkerColor(COLOR[a]);
		H1_xSec[a]->SetMarkerStyle(MARKER[a]);
		H1_xSec[a]->SetMarkerSize(1.65);
		H1_xSec[a]->DrawCopy("pe 9 same text0");
		H1_xSec[a]->SetName(Form("Xic0_xSec_%s_%s", Config[a][0], Config[a][1]));

		//Draw systematic error box
		for (int b=0; b<H1_xSec[a]->GetNbinsX(); b++)
		{
			if (H1_xSec[a]->GetBinContent(b+1) == 0) continue;
			if (H1_syst[a][0]->GetBinContent(b+1) == 0) continue;
			const float BinQW = H1_xSec[a]->GetBinWidth(b+1)/4;
			const float x1 = H1_xSec[a]->GetBinCenter(b+1) - BinQW;
			const float x2 = H1_xSec[a]->GetBinCenter(b+1) + BinQW;
			const float y1 = H1_xSec[a]->GetBinContent(b+1) - H1_syst[a][0]->GetBinContent(b+1);
			const float y2 = H1_xSec[a]->GetBinContent(b+1) + H1_syst[a][0]->GetBinContent(b+1);
			TBox* SystBox = new TBox(x1, y1, x2, y2);
			SystBox->SetLineColor(COLOR[a]);
			SystBox->SetFillStyle(0);
			SystBox->Draw("9 same");
		}

		TString CONFIG = Form("%s, [%s]", Config[a][0], Config[a][1]);
		CONFIG.ReplaceAll("to", ", ");
		CONFIG.ReplaceAll("p1", ".1");
		L1->AddEntry(H1_xSec[a], (const char*)CONFIG, "lp");
		if (a==nConfig-1) L1->Draw("9 same");
	}
	c1->Print(Form("c1z_xSecFin%s.png", INELLgt0?"_INEL0":""));

	//TFile* Fout = new TFile("out_xSec.root", "recreate");
	//for (int a=0; a<nConfig; a++) H1_xSec[a]->Write();
	//c1->Write();
	//Fout->Close();
	#endif

	#if 0
	TCanvas* c2;
	for (int a=0; a<nConfig; a++)
	{
		const int COLORSYS[] = {1, 2, 210, 4, 93, 6, 7};
		c2 = new TCanvas(Form("SystErr_%s_%s", Config[a][0], Config[a][1]), "", 800*2, 600*2);
		c2->cd()->SetGrid();
		for (int b=0; b<nSystErr+1; b++)
		{
			H1_syst[a][b]->SetLineColor(COLORSYS[b]);
			if (b==0) H1_syst[a][b]->SetMaximum(H1_syst[a][b]->GetMaximum()*1.2);
			if (b==0) H1_syst[a][b]->SetTitle(Form("%s;pT", c2->GetName()));
			if (b!=0) H1_syst[a][b]->SetFillColor(COLORSYS[b]);
			if (b!=0) H1_syst[a][b]->SetFillStyle(3003+b);
			H1_syst[a][b]->DrawCopy(b==0?"hist e text90":"hist e same");
		}//b
		//c2->Print(Form("%s.png", c2->GetName()));
	}//a
	#endif

	return;
}//Main
