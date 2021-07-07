#ifndef XI0CANAFUNCTION
#define XI0CANAFUNCTION

#include "AliNormalizationCounter.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <RooUnfold.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldSvd.h>
#include <TSVDUnfold_local.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

//Get normalization factor via AliNormalizationCounter: require original train output
//-------------------------------------------------------------------------------------------------
double GetNormFactor(const char* inFile, const char* TRIG, const char* PERC, bool INELLgt0 = false)
{
	//Check if original train output file exists
	TFile* F = TFile::Open(inFile);
	if (!F || F->IsZombie()) { cout <<Form("Cannot find file %s! Stop.\n", inFile); assert(false); }

	//Check if the counter exists
	const char* SDIR = "PWG3_D2H_Xic02eXipp13TeV_HM";
	const char* ANC_name = Form("ANC%s_%s_%s", INELLgt0?"INEL0":"", TRIG, PERC);
	AliNormalizationCounter* ANC = (AliNormalizationCounter*)F->Get(Form("%s/%s", SDIR, ANC_name));
	if (!ANC) {	cout <<Form("Cannot find ANC object %s/%s! Stop.\n", SDIR, ANC_name); assert(false); }

	const double fNorm = ANC->GetNEventsForNorm();
	cout <<Form("Normalization factor for %s_%s: %10.9f x 1.e9\n", TRIG, PERC, fNorm/1.e9);

	F->Close();
	return fNorm;
}//GetNormFactorANC

//----------------------
TH1D* ApplyPreFilterEff(
		TH1D* H1, TFile* F_data,
		const char* Cut, const char* CutFlag,
		const char* Suffix, bool Show = false
		)
{
	//Get prefilter efficiency
	TH1D* preF_eff_n = (TH1D*)F_data->Get(Form("hpre_%s_%s_nu", Cut, CutFlag)); //Numerator
	TH1D* preF_eff_d = (TH1D*)F_data->Get(Form("hpre_%s_%s_de", Cut, CutFlag)); //Denominator
	TH1D* preF_eff = (TH1D*)preF_eff_n->Clone(Form("hPreFEff_%s_%s", Cut, CutFlag)); preF_eff->Reset();
	preF_eff->Divide(preF_eff_n, preF_eff_d, 1, 1, "b"); //1, 1, b: Weighting on n/d and compute binomial error

	//Process return histogram: apply prefilter efficiency on raw yields
	TH1D* H1R = (TH1D*)H1->Clone(Form("%s_%s_preFCorr", H1->GetName(), Suffix));
	H1R->Divide(preF_eff);

	//+++++++++++++++++++++++++++++++++++++++++++

	//Draw
	if (Show)
	{
		gStyle->SetOptStat(0);
	    gStyle->SetPaintTextFormat("4.3f");
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900/2);
		c1->SetName(Form("c1a_%s_%s_%s_preFCorr", Suffix, Cut, CutFlag));
		c1->Divide(2, 1);

		c1->cd(1);//->SetGrid();
		preF_eff->GetYaxis()->SetRangeUser(preF_eff->GetMinimum()-0.03, preF_eff->GetMaximum()+0.03);
		preF_eff->SetLineColor(1);
		preF_eff->SetMarkerColor(2);
		preF_eff->SetMarkerSize(1.4);
		preF_eff->SetTitle(Form("%s, prefilter eff;pT;#epsilon", Suffix));
		preF_eff->DrawCopy("hist e text45");

		c1->cd(2)->SetGridy();
		TLegend* L1 = new TLegend(0.50, 0.65, 0.80, 0.85);
		for (int a=0; a<2; a++)
		{
			TH1D* H1Temp;
			if (a==0) H1Temp = (TH1D*)H1 ->Clone(Form("%s_showA", H1 ->GetName()));
			if (a==1) H1Temp = (TH1D*)H1R->Clone(Form("%s_showA", H1R->GetName()));

			H1Temp->SetLineColor(a+1);
			H1Temp->SetMarkerColor(a+1);
			H1Temp->SetMarkerStyle(20 + a*4);
			if (a==0) H1Temp->SetTitle(Form("%s, preFilter eff correction;pT", Suffix));
			H1Temp->DrawCopy(a==0?"hist e":"hist e same");
			L1->AddEntry(H1Temp, a==0?"Raw":"Raw / preFilter eff", "lp");
		}
		L1->Draw();

		c1->Print(Form("%s.png", c1->GetName()));
	}//Show

	preF_eff_n->Delete();
	preF_eff_d->Delete();
	preF_eff->Delete();
	return H1R;
}//ApplyPreFilterEff

//-------------------
TH1D* ApplyUnfolding(
		TH1D* H1, TFile* F_MC,
		const char* Cut, const char* CutFlag, const char* UFMethod, const int nIter,
		const char* Suffix, bool Show = false
		)
{
	if (Show) cout <<Form("Perform unfolding by using '%s' + '%s'... ", F_MC->GetName(), UFMethod);

	//Prepare input for unfolding
	TH1D* heXiPair   = (TH1D*)F_MC->Get(Form("hMCRecoLevPair_%s_%s", Cut, CutFlag));
	TH1D* hXic0      = (TH1D*)F_MC->Get(Form("hMCRecoLevXic0_%s_%s", Cut, CutFlag));
	TH2D* hRM_unfold = (TH2D*)F_MC->Get(Form("hRPM_%s_%s_un",        Cut, CutFlag));

	//Perform unfolding
    RooUnfoldResponse UFResponse(heXiPair, hXic0, hRM_unfold);
	TH1D* hReco = (TH1D*)hXic0->Clone(Form("hReco_%s_%s_%s", Cut, CutFlag, Suffix));
	hReco->Reset();

	if (!strcmp(UFMethod, "Bayes"))
	{
		TH1D* H1Temp = (TH1D*)H1->Clone(Form("%s_temp", H1->GetName()));
		RooUnfoldBayes UFBayes(&UFResponse, H1Temp, nIter);
		hReco = (TH1D*)UFBayes.Hreco();
		H1Temp->Delete();
	}
	else if (!strcmp(UFMethod, "Svd"))
	{
		cout <<endl;
		TH1D* H1Temp = (TH1D*)H1->Clone(Form("%s_temp", H1->GetName()));
		RooUnfoldSvd UFSvd(&UFResponse, H1Temp, nIter);
		hReco = (TH1D*)UFSvd.Hreco();
		H1Temp->Delete();
	
		if (Show)
		{
			TSVDUnfold_local* UFSvdLoc = (TSVDUnfold_local*)UFSvd.Impl();
			TH1D* fDHist = UFSvdLoc->GetD();
			for (int a=0; a<fDHist->GetNbinsX(); a++) cout <<Form("Bin %i: %f\n", a+1, fDHist->GetBinContent(a+1));
			//fDHist->Delete(); //Deletion causes segfault
		}
	}
	else { cout <<"Unknown unfolding method! Stop.\n"; assert(false); }
	
	//Process return histogram: match its binning to the input
	const int pTBinN = H1->GetNbinsX();	double pTBin[pTBinN];
	for (int a=0; a<pTBinN+1; a++) pTBin[a] = H1->GetXaxis()->GetBinLowEdge(a+1);
	TH1D* H1R = (TH1D*)hReco->Rebin(pTBinN, Form("%s_unfolded", H1->GetName()), pTBin);

	//+++++++++++++++++++++++++++++++++++++++++++
	
	if (Show)
	{
		gStyle->SetOptStat(0);
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900/2);
		c1->SetName(Form("c1b_%s_%s_%s_unfold", Suffix, Cut, CutFlag));
		c1->Divide(2, 1);

		c1->cd(1)->SetGrid();
		hRM_unfold->GetZaxis()->SetLabelSize(0.03);
		hRM_unfold->SetMarkerColor(2);
		hRM_unfold->SetMarkerSize(1.2);
		hRM_unfold->SetStats(false);
		hRM_unfold->SetTitle(Form("%s, RM (%s);pT (eXi);pT (Xic0)", Suffix, UFMethod));
		hRM_unfold->DrawCopy("colz");// text45");

		c1->cd(2)->SetGridy();
		TLegend* L1 = new TLegend(0.5, 0.65, 0.75, 0.85);
		for (int a=0; a<2; a++)
		{
			TH1D* H1Temp;
			if (a==0) H1Temp = (TH1D*)H1 ->Clone(Form("%s_showB", H1 ->GetName()));
			if (a==1) H1Temp = (TH1D*)H1R->Clone(Form("%s_showB", H1R->GetName()));

			H1Temp->SetLineColor(2 * (a+1));
			H1Temp->SetMarkerColor(2 * (a+1));
			H1Temp->SetMarkerStyle(24 - a*4);
			if (a==0) H1Temp->SetTitle(Form("%s, unfolding;pT", Suffix));
			H1Temp->DrawCopy(a==0?"hist e":"hist e same");
			L1->AddEntry(H1Temp, a==0?"Raw":"Unfolded", "lp");
		}//a
		L1->Draw("same");

		c1->Print(Form("%s.png", c1->GetName()));
	}//Show

	heXiPair->Delete();
	hXic0->Delete();
	hRM_unfold->Delete();
	hReco->Delete();
	return H1R;
}//ApplyUnfolding

//-----------------
TH1D* ApplyXic0Eff(
		TH1D* H1, TFile* F_MC, bool IsWeighted,
		const char* Cut, const char* CutFlag,
		const char* Suffix, bool Show = false
		)
{
	//Get Xic0 efficiency
	TH1D* Xic0_eff_n = (TH1D*)F_MC->Get(Form("hMCRecoLevXic0_%s_%s", Cut, CutFlag));
	TH1D* Xic0_eff_d = (TH1D*)F_MC->Get(Form("hMCGenInclusiveXic0_%s", IsWeighted?"W":"woW"));
	TH1D* Xic0_eff = (TH1D*)Xic0_eff_n->Clone(Form("hXic0Eff_%s_%s", Cut, CutFlag)); Xic0_eff->Reset();
	Xic0_eff->Divide(Xic0_eff_n, Xic0_eff_d, 1, 1, "b"); //1, 1, b: Weighting on n/d and compute binomial error

	//Process return histogram: apply Xic0 efficiency on input
	TH1D* H1R = (TH1D*)H1->Clone(Form("%s_Xic0EffCorr", H1->GetName()));
	H1R->Divide(Xic0_eff);

	//+++++++++++++++++++++++++++++++++++++++++++

	//Draw
	if (Show)
	{
		gStyle->SetOptStat(0);
	    gStyle->SetPaintTextFormat("4.3f");
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900/2);
		c1->SetName(Form("c1c_%s_%s_%s_Xic0Eff", Suffix, Cut, CutFlag));
		c1->Divide(2, 1);

		c1->cd(1);//->SetGrid();
		Xic0_eff->GetYaxis()->SetRangeUser(-0.01, Xic0_eff->GetMaximum()*1.25);
		Xic0_eff->SetLineColor(1);
		Xic0_eff->SetMarkerColor(210);
		Xic0_eff->SetMarkerSize(1.4);
		Xic0_eff->SetTitle(Form("%s, Xic0 eff;pT;#epsilon", Suffix));
		Xic0_eff->DrawCopy("hist e text45");

		c1->cd(2)->SetGridy();
		TLegend* L1 = new TLegend(0.50, 0.65, 0.85, 0.85);
		for (int a=0; a<2; a++)
		{
			TH1D* H1Temp;
			if (a==0) H1Temp = (TH1D*)H1 ->Clone(Form("%s_showC", H1 ->GetName()));
			if (a==1) H1Temp = (TH1D*)H1R->Clone(Form("%s_showC", H1R->GetName()));
			H1Temp->SetLineColor(a==0?4:210);
			H1Temp->SetMarkerColor(a==0?4:210);
			H1Temp->SetMarkerStyle(20 + a*4);

			const float scaleF = 0.01;
			if (a==0) H1Temp->GetYaxis()->SetRangeUser(-100, H1R->GetMaximum()*1.5*scaleF);
			if (a==0) H1Temp->SetTitle(Form("%s, Xic0 eff correction;pT", Suffix));
			if (a==1) H1Temp->Scale(scaleF);

			H1Temp->DrawCopy(a==0?"hist e":"hist e same");
			L1->AddEntry(H1Temp, a==0?"Raw":Form("(Raw / Xic0 eff) * %3.2f", scaleF), "lp");
		}
		L1->Draw();

		c1->Print(Form("%s.png", c1->GetName()));
	}//Show

	Xic0_eff_n->Delete();
	Xic0_eff_d->Delete();
	Xic0_eff->Delete();
	return H1R;
}//ApplyXic0Eff

//----------
TH1D* GetXS(
		TH1D* H1, double normF, double V0xSec, double Xic0SemiLBR,
		const char* Suffix, float xMin = -999., float xMax = 999.
		)
{
	TH1D* H1R = (TH1D*)H1->Clone(Form("Xic0XS_%s", Suffix));
	H1R->SetTitle("");
	H1R->Reset();

	for (int a=0; a<H1->GetNbinsX(); a++)
	{
		if (H1->GetBinCenter(a+1) < xMin || H1->GetBinCenter(a+1) > xMax) continue;

		const double Val = H1->GetBinContent(a+1);
		const double Err = H1->GetBinError(a+1);

		const double BinW   = H1->GetBinWidth(a+1);
		const double CConj  = 2.; //e-Xi pair charge conjugate
		const double DeltaY = 1.;
		const double IntL   = normF/V0xSec;
		const double xSecF  = BinW * CConj * DeltaY * IntL * Xic0SemiLBR;

		H1R->SetBinContent(a+1, Val/xSecF);
		H1R->SetBinError  (a+1, Err/xSecF);
	}

	return H1R;
}//GetXS

#if 0
//Get weighting factor
//---------------------------------------------------------------------------------------------------------------
void GetWeightingFactor(TFile* F_MC, TH1D* H1XS_data, vector<float> pTvec, const char* Suffix, bool Show = false)
{
	if (strcmp("out_MC_raw.root", F_MC->GetName())) { cout <<"You're using improper MC file! Stop.\n"; return; }

	TH1D* H1XS_MC = (TH1D*)F_MC->Get("hMCGenInclusiveXic0_woW");
	if (!H1XS_MC) { cout <<"Cannot find 'hMCGenInclusiveXic0_woW'! Stop.\n"; return; }
	else cout <<Form("Open '%s' and '%s' for pT matching...\n", H1XS_data->GetName(), H1XS_MC->GetName());

	if (H1XS_data->GetNbinsX() != H1XS_MC->GetNbinsX()) { cout <<"Binnings NOT match! Stop.\n"; return; }

	//-------------------------------------------

	enum {data, MC};
	TH1D* H1XS[2];
	H1XS[data] = (TH1D*)H1XS_data->Clone(Form("pTW_data_%s", Suffix));
	H1XS[MC] = (TH1D*)H1XS_MC->Clone(Form("pTW_MC_%s", Suffix));

	for (int a=0; a<2; a++) //Data and MC
	{
		//Empty the bin "0 < pT < 1" if it still filled
		const float Bin1Cnt = H1XS[a]->GetBinCenter(1);
		const float Bin1Val = H1XS[a]->GetBinContent(1);
		if (Bin1Cnt>0 && Bin1Cnt<1 && Bin1Val!=0) { H1XS[a]->SetBinContent(1, 0); H1XS[a]->SetBinError(1, 0); }

		//Normalize by bin width
		for (int b=0; b<H1XS[a]->GetNbinsX(); b++)
		{
			const float BinC = H1XS[a]->GetBinContent(b+1);
			const float BinW = H1XS[a]->GetBinWidth(b+1);
			if (BinC!=0 && fabs(BinW-1)>1.e-3) H1XS[a]->SetBinContent(b+1, BinC/BinW);
		}

		//Normalize by self integral
		const double INTEGRAL = H1XS[a]->Integral();
		H1XS[a]->Scale(1/INTEGRAL);
	}//a, data or MC

	//-------------------------------------------

	enum {CT, UD, DU};
	TH1D* H1Ratio[3];
	H1Ratio[CT] = (TH1D*)H1XS[data]->Clone(Form("pTW_ratio_%s", Suffix));
	H1Ratio[CT]->Divide(H1XS[MC]);
	H1Ratio[UD] = (TH1D*)H1Ratio[CT]->Clone(Form("%s_updown", H1Ratio[CT]->GetName()));
	H1Ratio[DU] = (TH1D*)H1Ratio[CT]->Clone(Form("%s_downup", H1Ratio[CT]->GetName()));
	H1Ratio[UD]->Reset(); H1Ratio[UD]->Sumw2(false);
	H1Ratio[DU]->Reset(); H1Ratio[DU]->Sumw2(false);

	for (int a=0; a<H1Ratio[CT]->GetNbinsX(); a++)
	{
		const float BinC  = H1Ratio[CT]->GetBinContent(a+1); if (BinC == 0) continue;
		const float BinEL = H1Ratio[CT]->GetBinErrorLow(a+1);
		const float BinEU = H1Ratio[CT]->GetBinErrorUp(a+1);
		const float pTMean = H1Ratio[CT]->GetBinCenter(a+1);
		if (pTMean < 4.0)
		{
			H1Ratio[UD]->SetBinContent(a+1, BinC + BinEU); 
			H1Ratio[DU]->SetBinContent(a+1, BinC - BinEL);
		}
		else if (pTMean < 5.0) // 4 < pT < 5
		{
			H1Ratio[UD]->SetBinContent(a+1, BinC);
			H1Ratio[DU]->SetBinContent(a+1, BinC);
		}
		else // pT > 5
		{
			H1Ratio[UD]->SetBinContent(a+1, BinC - BinEL);
			H1Ratio[DU]->SetBinContent(a+1, BinC + BinEU);
		}
	}//a, 3 ratio

	TF1* F1Ratio[3];
	for (int a=0; a<3; a++)
	{
		//F1Ratio[a] = new TF1(Form("F1Ratio_%i", a), "TMath::Exp([0] + [1]*x)", 1, pTvec.back());
		//H1Ratio[a]->Fit(F1Ratio[a]->GetName(), "EQR0", "", 1, pTvec.back());
		F1Ratio[a] = new TF1(Form("F1Ratio_%i", a), "expo", 1, pTvec.back());
		H1Ratio[a]->Fit(F1Ratio[a]->GetName(), "EQR0");
	}//a, 3 ratio

	//-------------------------------------------

	if (Show)
	{
		TCanvas* c1 = new TCanvas("c1", "", 1600, 1200/2);
		c1->SetName(Form("c1f_pTWeight_%s", Suffix));
		c1->Divide(2, 1);

		c1->cd(1)->SetLogy();
		TLegend* L1 = new TLegend(0.65, 0.7, 0.85, 0.85);
		L1->SetMargin(0.5);
		for (int a=0; a<2; a++)
		{
			if (a==0) H1XS[a]->GetYaxis()->SetRangeUser(1.e-7, 2.0);
			if (a==0) H1XS[a]->SetTitle(Form("%s;pT;dN/NdPT", Suffix));
			H1XS[a]->SetStats(false);
			H1XS[a]->SetLineColor(a+1);
			H1XS[a]->SetMarkerColor(a+1);
			H1XS[a]->SetMarkerSize(1.2);
			H1XS[a]->SetMarkerStyle(a*4 + 20);
			H1XS[a]->DrawCopy(a==0?"hist ep":"hist ep same");
			L1->AddEntry(H1XS[a], a==0?"data":"MC", "lp");
		}
		L1->Draw("same");

		c1->cd(2);
		TLegend* L2 = new TLegend(0.45, 0.5, 0.85, 0.85);
		L2->SetHeader("Anchor point: 4 < pT < 5", "C");
		L2->SetMargin(0.3);
		const int RCOLOR[3] = {1, 2, 4};
		for (int a=0; a<3; a++)
		{
			if (a==0) H1Ratio[a]->GetYaxis()->SetRangeUser(-0.5, 3.);
			if (a==0) H1Ratio[a]->SetTitle(Form("%s;pT;Weighting factor", Suffix));
			H1Ratio[a]->SetStats(false);
			H1Ratio[a]->SetLineColor(RCOLOR[a]);
			H1Ratio[a]->SetMarkerColor(RCOLOR[a]);
			H1Ratio[a]->SetMarkerStyle(a==0?20:24);
			H1Ratio[a]->SetMarkerSize(1.2);
			H1Ratio[a]->DrawCopy(a==0?"pe":"p same");

			F1Ratio[a]->SetLineColor(RCOLOR[a]);
			//F1Ratio[a]->SetLineStyle(2);
			F1Ratio[a]->Draw("l 9 same");

			const float p0Val = F1Ratio[a]->GetParameter(0);
			const float p1Val = F1Ratio[a]->GetParameter(1);
			const float p0Err = F1Ratio[a]->GetParError(0);
			const float p1Err = F1Ratio[a]->GetParError(1);
			L2->AddEntry(F1Ratio[a], Form("p0: %7.6f #pm %7.6f", p0Val, p0Err), "lp");
			L2->AddEntry((TObject*)0, Form("p1: %7.6f #pm %7.6f", p1Val, p1Err), "");
		}
		L2->Draw();

		c1->Print(Form("%s.png", c1->GetName()));
	}//Show

	for (int a=0; a<2; a++) H1XS[a]->Delete();
	for (int a=0; a<3; a++) H1Ratio[a]->Delete();
	//for (int a=0; a<3; a++) F1Ratio[a]->Delete();

	return;
}//GetWeightingFactor
#endif

#endif //XI0CANAFUNCTION
