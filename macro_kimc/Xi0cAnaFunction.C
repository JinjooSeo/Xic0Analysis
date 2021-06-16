#ifndef XI0CANAFUNCTION
#define XI0CANAFUNCTION

#include "AliNormalizationCounter.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
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

//-------------------------------------------------
TH1D* GetCMSLbSpectrum(TFile* F, bool Show = false)
{
	TDirectoryFile* sDIR = (TDirectoryFile*)F->Get("Table 2");
	TH1D* hCMSLb_val = (TH1D*)sDIR->Get("Hist1D_y1");
	TH1D* hCMSLb_err = (TH1D*)sDIR->Get("Hist1D_y1_e1");
	hCMSLb_val->Scale(0.001/4);
	hCMSLb_err->Scale(0.001/4);

	TH1D* hCMSLb = (TH1D*)hCMSLb_val->Clone("CMSLambda_b");
	for (int a=0; a<hCMSLb->GetNbinsX(); a++) hCMSLb->SetBinError(a+1, hCMSLb_err->GetBinContent(a+1));	

	//Draw
	if (Show)
	{
		gStyle->SetPaintTextFormat("4.3f");
		TCanvas* c1 = new TCanvas("c1y_CMSLb", "", 1600/2, 900/2); c1->cd()->SetGrid();
		TH1D* hCMSLb_temp = (TH1D*)hCMSLb->Clone();
		hCMSLb_temp->SetLineColor(1);
		hCMSLb_temp->SetMarkerSize(1.4);
		hCMSLb_temp->GetYaxis()->SetRangeUser(-0.01, hCMSLb_temp->GetMaximum()*1.2);
		hCMSLb_temp->SetTitle("CMS lambda_b");
		hCMSLb_temp->DrawCopy("hist e text0");
		c1->Print(Form("%s.png", c1->GetName())); delete c1;
	}

	hCMSLb_val->Delete();
	hCMSLb_err->Delete();
	return hCMSLb;
}//GetCMSLbSpectrum

//-----------------------------
TH1D* GetPreFilterCorrSpectrum(
		TFile* F_data, const char* Cut, const char* CutFlag, const char* Suffix = "", bool Show = false
		)
{
	//Get prefilter efficiency
	TH1D* preF_eff_n = (TH1D*)F_data->Get(Form("hpre_%s_%s_nu", Cut, CutFlag)); //Numerator
	TH1D* preF_eff_d = (TH1D*)F_data->Get(Form("hpre_%s_%s_de", Cut, CutFlag)); //Denominator

	TH1D* preF_eff = (TH1D*)preF_eff_d->Clone(Form("hpreff_%s_%s", Cut, CutFlag)); preF_eff->Reset();
	preF_eff->Divide(preF_eff_n, preF_eff_d, 1, 1, "b"); //1, 1, b: Weighting on n/d and compute binomial error

	//Return histogram: apply prefilter efficiency on measured yields
	TH1D* hMeas = (TH1D*)F_data->Get(Form("hRawPt_%s_%s", Cut, CutFlag))->Clone();
	hMeas->SetName(Form("%s_%s_preFCorr", Suffix, hMeas->GetName()));
	hMeas->Divide(preF_eff);

	if (Show)
	{
		gStyle->SetOptStat(0);
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900/2);
		c1->SetName(Form("c1a_preFCorr_%s_%s_%s", Cut, CutFlag, Suffix));
		c1->Divide(2, 1);

		c1->cd(1)->SetGrid();
		preF_eff->SetLineColor(1);
		preF_eff->SetMarkerColor(4);
		preF_eff->SetMarkerSize(1.4);
		preF_eff->GetYaxis()->SetRangeUser(0.9, 1.1);
		preF_eff->SetTitle(Form("%s, prefilter eff;pT;#epsilon", Suffix));
		preF_eff->DrawCopy("hist e text45");

		c1->cd(2)->SetGrid();
		TH1D* hMeas_orig = (TH1D*)F_data->Get(Form("hRawPt_%s_%s", Cut, CutFlag))->Clone();
		hMeas_orig->SetLineColor(1);
		hMeas_orig->SetMarkerColor(1);
		hMeas_orig->SetMarkerStyle(20);
		hMeas_orig->SetTitle(Form("%s, Raw yields;pT", Suffix));
		hMeas_orig->DrawCopy("hist e");
		TH1D* hMeas_temp = (TH1D*)hMeas->Clone(Form("%s_copy", hMeas->GetName()));
		hMeas_temp->SetLineColor(2);
		hMeas_temp->SetMarkerColor(2);
		hMeas_temp->SetMarkerStyle(24);
		hMeas_temp->DrawCopy("hist e same");
		TLegend* Leg = new TLegend(0.50, 0.65, 0.75, 0.85);
		Leg->AddEntry(hMeas_orig, "w/o preFilter corr", "lp");
		Leg->AddEntry(hMeas_temp, "w/ preFilter corr", "lp");
		Leg->Draw("same");

		c1->Print(Form("%s.png", c1->GetName())); delete c1;
		hMeas_orig->Delete();
		hMeas_temp->Delete();
	}

	preF_eff_n->Delete();
	preF_eff_d->Delete();
	preF_eff->Delete();
	return hMeas;
}//GetPreFilterCorrSpectrum

#if 0
//--------------------------
TH1D* GetBottomCorrSpectrum(
		TFile* F_MC, TH1D* hCMSLb, TH1D* hMeas_orig, vector<float>pTvec,
		const char* Cut, const char* CutFlag, const char* Suffix = "", bool Show = false
		)
{
	//Constants
	const double BC_ScaleFactor[9] =
	{
		1.53313, 1.69604, 1.806626, 1.887637, 1.950308, 2.018669, 2.121922, 2.249672, 2.439034
	};
	const double BE_ScaleFactor[9] =
	{
		0.002688548, 0.003208365, 0.004097259, 0.005324936, 0.006890631,
		0.009749647, 0.01760265,  0.03540223,  0.06652629
	};
	const double BRFraction = (3.9 * 10e-4) / (5.8 * 10e-5);
	const double Luminosity = 1.88554e+09/(57.8*1000000); //pp 13 TeV integrated luminosity

	const int npTbins = pTvec.size();
	double pTbins[npTbins];	for (int i=0; i<npTbins; i++) pTbins[i] = pTvec[i];

	//+++++++++++++++++++++++++++++++++++++++++++

	//1) Fit 7TeV CMS Lb spectrum with Tasllis function
	cout <<"\n GetBottomCorrSpectrum:: perform fit..." <<endl;
	TF1 *fTsallis = new TF1("Lambdab", "[0]*x*(pow(1+(sqrt(pow(x,2)+pow(5.619,2))-5.619)/(7.6*1.1),-7.6))", 0, 50);
	hCMSLb->Fit("Lambdab", "0");
	fTsallis->SetParameters(fTsallis->GetParameters());
	cout <<endl;

	//2) Multiply scale factor, to convert 7 TeV Lb to 13TeV Lb
	TH1D* hScaleFactor = new TH1D("hScaleFactor", "", 9, pTbins); //B meson ratio -> B(13TeV)/B(7TeV)
	for (int i=0; i<9; i++)
	{
		hScaleFactor->SetBinContent(i+1, BC_ScaleFactor[i]);
		hScaleFactor->SetBinError  (i+1, BE_ScaleFactor[i]);
	}
	TH1D* h13TeVLb = new TH1D("h13TeVLb", "", 9, pTbins);
	for (int i=0; i<9; i++)
	{
		const double tempLbVal = fTsallis->Eval( h13TeVLb->GetBinCenter(i+1) ) * hScaleFactor->GetBinContent(i+1);
		h13TeVLb->SetBinContent(i+1, tempLbVal); //Error is not assigned since we don't know the error at 1 to 10 pT
	}

	//3) Multiply Branching ratio fraction to convert to Xib
	TH1D* h13TeVXib = (TH1D*)h13TeVLb->Clone("h13TeVXib");
	h13TeVXib->Scale(BRFraction);

	//4) Calculate Xib yield
	TH1D* h13TeVXibRaw = new TH1D("h13TeVXibRaw","", 9, pTbins);
	for( int i=0; i<9; i++)
	{
		const double pTbinW = pTbins[i+1] - pTbins[i];
		const double tempXibVal = h13TeVXib->GetBinContent(i+1) * pTbinW*2 * Luminosity; //x2? //kimc
		h13TeVXibRaw->SetBinContent(i+1, tempXibVal);
	}

	//+++++++++++++++++++++++++++++++++++++++++++

	TH1D* hGenXib  = (TH1D*)F_MC->Get("XibGen05"); //Number of Xib in generation level
	TH1D* hRecoXib = (TH1D*)F_MC->Get(Form("hMCRecoLevXib_%s_%s", Cut, CutFlag)); //Number of Xib in reco lv
	TH1D* hXibEff  = (TH1D*)hRecoXib->Clone(Form("hXibEff_%s_%s", Cut, CutFlag)); hXibEff->Reset();
	hXibEff->Divide(hRecoXib, hGenXib, 1, 1, "b"); //Xib efficiency
	h13TeVXibRaw->Multiply(hXibEff);

	#if 1
	TH1D* hRecoeXi = (TH1D*)F_MC->Get(Form("hMCRecoLevPairXib_%s_%s", Cut, CutFlag)); //Number of eXi from Xib
	TH2D* hRM_Xib  = (TH2D*)F_MC->Get(Form("hRPM_%s_%s_Xib", Cut, CutFlag)); //Response matrix of Xib and eXi

	//5) Convert Xib spectrum to eXi spectrum
	RooUnfoldResponse Response(hRecoXib, hRecoeXi, hRM_Xib);
	RooUnfoldBinByBin Unfolding(&Response, h13TeVXibRaw);
	TH1D* heXiFromXib = (TH1D*)Unfolding.Hreco();

	//6) Return histogram: add eXi from Xib to eXi pair(RS - WS)
	TH1D* hMeas = (TH1D*)hMeas_orig->Clone(Form("%s_BCorr", hMeas_orig->GetName()));
	hMeas->Add(heXiFromXib);

	//+++++++++++++++++++++++++++++++++++++++++++
	
	if (Show)
	{
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900);
		c1->SetName(Form("c1b_XibCorr_%s_%s_%s", Cut, CutFlag, Suffix));
		c1->Divide(2, 2);

		c1->cd(1)->SetGrid();
		hXibEff->SetLineColor(1);
		hXibEff->SetMarkerColor(4);
		hXibEff->SetMarkerSize(1.3);
		hXibEff->SetMinimum(0.);
		hXibEff->SetTitle(Form("%s, Xib eff;pT", Suffix));
		hXibEff->DrawCopy("hist e text0");

		c1->cd(2)->SetGrid();
		hRM_Xib->SetStats(false);
		hRM_Xib->SetTitle(Form("%s, eXi-Xib RM;pT (Xib);pT (eXi)", Suffix));
		hRM_Xib->SetMarkerColor(2);
		hRM_Xib->SetMarkerSize(1.4);
		hRM_Xib->DrawCopy("colz text45");

		c1->cd(3)->SetGrid();
		heXiFromXib->SetStats(false);
		heXiFromXib->SetLineColor(2);
		heXiFromXib->SetMarkerColor(2);
		heXiFromXib->SetMarkerStyle(24);
		heXiFromXib->SetTitle(Form("%s, eXi from Xib;pT", Suffix));
		heXiFromXib->DrawCopy("hist e");
		h13TeVXibRaw->SetLineColor(1);
		h13TeVXibRaw->SetMarkerColor(1);
		h13TeVXibRaw->SetMarkerStyle(20);
		h13TeVXibRaw->DrawCopy("hist e same");
		TLegend* Leg1 = new TLegend(0.5, 0.65, 0.75, 0.85);
		Leg1->AddEntry(heXiFromXib, "eXi from Xib", "lp");
		Leg1->AddEntry(h13TeVXibRaw, "Xib raw", "lp");
		Leg1->Draw("same");

		c1->cd(4)->SetGrid();
		TH1D* hMeas_orig_temp = (TH1D*)hMeas_orig->Clone(Form("%s_temp", hMeas_orig->GetName()));
		hMeas_orig_temp->SetTitle(Form("%s, eXi yields w/ or w/o Xib corr;pT", Suffix));
		hMeas_orig_temp->SetLineColor(1);
		hMeas_orig_temp->SetMarkerColor(1);
		hMeas_orig_temp->SetMarkerStyle(20);
		hMeas_orig_temp->DrawCopy("hist e");
		TH1D* hMeas_temp = (TH1D*)hMeas->Clone(Form("%s_temp", hMeas->GetName()));
		hMeas_temp->SetLineColor(2);
		hMeas_temp->SetMarkerColor(2);
		hMeas_temp->SetMarkerStyle(24);
		hMeas_temp->DrawCopy("hist e same");
		TLegend* Leg2 = new TLegend(0.50, 0.65, 0.75, 0.85);
		Leg2->AddEntry(hMeas_temp, "Corrected", "lp");
		Leg2->AddEntry(hMeas_orig_temp, "NOT corrected", "lp");
		Leg2->Draw("same");

		c1->Print(Form("%s.png", c1->GetName())); delete c1;
		hMeas_orig_temp->Delete();
		hMeas_temp->Delete();
	}
	else
	{
		fTsallis    ->Delete();
		hScaleFactor->Delete();
		h13TeVLb    ->Delete();
		h13TeVXib   ->Delete();
		h13TeVXibRaw->Delete();
		hGenXib     ->Delete();
		hRecoXib    ->Delete();
		hXibEff     ->Delete();
		hRecoeXi    ->Delete();
		hRM_Xib     ->Delete();
		heXiFromXib ->Delete();
	}

	return hMeas;
	#endif
}//GetBottomBayronCorrectedSpectrum
#endif

//kimc: get either weighted or unweighted by providing relevant ROOT file
//-----------------------------------------------------------------------
TH1D* GetUnfoldedSpectrum(
		TFile* F_MC, TH1D* hMeas, vector<float>pTvec,
		const char* Cut, const char* CutFlag, const char* Method, const int nIter,
		const char* Suffix, bool Show = false
		)
{
	cout <<Form("GetUnfoldedSpectrum:: by using %s and method %s...", F_MC->GetName(), Method) <<endl;
	TH1D* hMeas_temp = (TH1D*)hMeas->Clone(Form("%s_temp", hMeas->GetName()));

	const int npTbins = pTvec.size();
	double pTbins[npTbins];	for (int i=0; i<npTbins; i++) pTbins[i] = pTvec[i];

	TH1D* heXiPair   = (TH1D*)F_MC->Get(Form("hMCRecoLevPair_%s_%s", Cut, CutFlag));
	TH1D* hXic0      = (TH1D*)F_MC->Get(Form("hMCRecoLevXic0_%s_%s", Cut, CutFlag));
	TH2D* hRM_unfold = (TH2D*)F_MC->Get(Form("hRPM_%s_%s_un",        Cut, CutFlag));

	RooUnfoldResponse Unfold(heXiPair, hXic0, hRM_unfold);
	TH1D* hReco = new TH1D(Form("hReco_%s_%s_%s", Cut, CutFlag, Suffix), "", npTbins-1, pTbins);
	if ( !strcmp(Method, "Bayes") )
	{
		RooUnfoldBayes UFBayes(&Unfold, hMeas_temp, nIter);
		hReco = (TH1D*)UFBayes.Hreco();
	}
	else if ( !strcmp(Method, "Svd") )
	{
		RooUnfoldSvd UFSvd(&Unfold, hMeas_temp, nIter);
		hReco = (TH1D*)UFSvd.Hreco();

		TSVDUnfold_local* UFSvdLoc = (TSVDUnfold_local*)UFSvd.Impl();
		TH1D* fDHist = UFSvdLoc->GetD();
		for (int i=0; i<fDHist->GetNbinsX(); i++) cout <<Form("%2i: %f", i, fDHist->GetBinContent(i+1)) <<endl;
	}
	else
	{
		cout <<"Unknown method! Stop.\n";
		hReco->Reset();
		return hReco;
	}

	//Return histogram
	TH1D* hUnfolded = new TH1D("hUnfolded_%s_%s", "", npTbins-1, pTbins);
	hUnfolded->SetName(Form("%s_unfolded", hMeas->GetName()));
	for (int i=0; i<npTbins; i++)
	{
		hUnfolded->SetBinContent(i+1, hReco->GetBinContent(i+1));
		hUnfolded->SetBinError  (i+1, hReco->GetBinError  (i+1));
	}

	if (Show)
	{
		TCanvas* c1 = new TCanvas("c1", "", 1600, 900/2);
		c1->SetName(Form("c1c_unfold_%s_%s_%s", Cut, CutFlag, Suffix));
		c1->Divide(2, 1);
		
		c1->cd(1)->SetGrid();
		hRM_unfold->SetTitle(Form("%s, eXi-Xic0 RM;pT (eXi);pT (Xic0)", Suffix));
		hRM_unfold->SetMarkerColor(2);
		hRM_unfold->SetMarkerSize(1.2);
		hRM_unfold->SetStats(false);
		hRM_unfold->DrawCopy("colz text45");

		c1->cd(2)->SetGrid();
		hMeas_temp->SetLineColor(1);
		hMeas_temp->SetMarkerColor(1);
		hMeas_temp->SetMarkerStyle(20);
		hMeas_temp->SetTitle(Form("%s, before/after unfolding;pT", hMeas->GetName()));
		hMeas_temp->DrawCopy("hist e");
		TH1D* hUnfold_temp = (TH1D*)hUnfolded->Clone(Form("%s_temp", hUnfolded->GetName()));
		hUnfold_temp->SetLineColor(2);
		hUnfold_temp->SetMarkerColor(2);
		hUnfold_temp->SetMarkerStyle(24);
		hUnfold_temp->DrawCopy("hist e same");
		TLegend* Leg = new TLegend(0.5, 0.65, 0.75, 0.85);
		Leg->AddEntry(hMeas_temp, "Before (original)", "lp");
		Leg->AddEntry(hUnfold_temp, "After (unfolded)", "lp");
		Leg->Draw("same");
		
		c1->Print(Form("%s.png", c1->GetName())); delete c1;
		hUnfold_temp->Delete();
	}

	hMeas_temp->Delete();
	heXiPair  ->Delete();
	hXic0     ->Delete();
	hRM_unfold->Delete();
	hReco     ->Delete();

	return hUnfolded;
}//GetUnfoldedSpectrum

//kimc: get either weighted or unweighted by providing relevant ROOT file
//-----------------------------------------------------------------------
TH1D* GetEfficiency(
		TFile* F_MC, vector<float>pTvec,
		const char* Cut, const char* CutFlag, bool IsWeighted,
		const char* Suffix, bool Show = false
		)
{
	TH1D* hGen  = (TH1D*)F_MC->Get(Form("hMCGenInclusiveXic0_%s", IsWeighted?"W":"woW"));
	TH1D* hReco = (TH1D*)F_MC->Get(Form("hMCRecoLevXic0_%s_%s", Cut, CutFlag));

	/*
	const int npTbins = pTvec.size();
	double pTbins[npTbins];	for (int i=0; i<npTbins; i++) pTbins[i] = pTvec[i];
	TH1D* hEff_d = new TH1D(Form("hEffd_%s_%s", Cut, CutFlag), "", npTbins-1, pTbins);
	TH1D* hEff_n = new TH1D(Form("hEffn_%s_%s", Cut, CutFlag), "", npTbins-1, pTbins);
	for (int i=0; i<npTbins; i++)
	{
		if (i>hGen->GetNbinsX() || i>hReco->GetNbinsX()) continue;
		hEff_d->SetBinContent(i+1, hGen ->GetBinContent(i+1));
		hEff_n->SetBinContent(i+1, hReco->GetBinContent(i+1));
		hEff_d->SetBinError(i+1, hGen ->GetBinError(i+1));
		hEff_n->SetBinError(i+1, hReco->GetBinError(i+1));
	}

	TH1D* hEff = new TH1D("hEff", "", npTbins-1, pTbins); hEff->Sumw2();
	hEff->SetName(Form("%s_hEff_%s_%s", Suffix, Cut, CutFlag));
	hEff->Divide(hEff_n, hEff_d, 1, 1, "b");
	*/

	TH1D* hEff = (TH1D*)hGen->Clone("hEff"); hEff->Reset();
	hEff->SetName(Form("%s_hEff_%s_%s", Suffix, Cut, CutFlag));
	hEff->Divide(hReco, hGen, 1, 1, "b");

	if (Show)
	{
		TCanvas* c1 = new TCanvas("c1", "", 1600/2, 900/2);
		c1->SetName(Form("c1d_eff_%s_%s_%s", Cut, CutFlag, Suffix)); c1->cd()->SetGrid();

		TH1D* hEffTemp = (TH1D*)hEff->Clone(Form("%s_temp", hEff->GetName()));
		if (hEffTemp->GetMaximum()*1.1 < 0.1) hEffTemp->GetYaxis()->SetRangeUser(0, hEffTemp->GetMaximum()*1.1);
		else hEffTemp->GetYaxis()->SetRangeUser(0, 0.15);
		hEffTemp->SetStats(false);
		hEffTemp->SetTitle(Form("Xic0 efficiency, %s;pT", Suffix));
		hEffTemp->SetLineColor(1);
		hEffTemp->SetMarkerColor(4);
		hEffTemp->SetMarkerSize(1.4);
		hEffTemp->DrawCopy("hist e text60");

		c1->Print(Form("%s.png", c1->GetName())); delete c1;
		hEffTemp->Delete();
	}

	hGen ->Delete();
	hReco->Delete();
	return hEff;
}//GetEfficiency

//---------------------------------------------------------------------
double GetNormFactor(TFile* F_data, const char* TRIG, const char* PERC)
{
	TH2D* hNorm = (TH2D*)F_data->Get("hNorm_multV0");

	//Normalization histogram's y bins: all (0), kINT7 (1), kHMV0 (2), kHMSPD (3), and 'kHMV0 || kHMSPD' (4)
	TH1D* hNorm1D = new TH1D(); //Extract multiplicity info of specific trigger
	if      (!strcmp(TRIG, "MB"))    hNorm1D = (TH1D*)hNorm->ProjectionX(Form("hNorm1D_%s", TRIG), 2, 2);
	else if (!strcmp(TRIG, "HMV0"))  hNorm1D = (TH1D*)hNorm->ProjectionX(Form("hNorm1D_%s", TRIG), 3, 3);
	else if (!strcmp(TRIG, "HMSPD")) hNorm1D = (TH1D*)hNorm->ProjectionX(Form("hNorm1D_%s", TRIG), 4, 4);

	//Specify multiplicity range
	float mult[2] = {0};

	if      (!strcmp(PERC, "0to0p1"))  { mult[0] =  0.0; mult[1] =   0.1; }
	else if (!strcmp(PERC, "0to100"))  { mult[0] =  0.0; mult[1] = 100.0; }
	else if (!strcmp(PERC, "0p1to30")) { mult[0] =  0.1; mult[1] =  30.0; }
	else if (!strcmp(PERC, "30to100")) { mult[0] = 30.0; mult[1] = 100.0; }
	else { cout <<"Unknown PERC setup: returning 0 nomalization factor...\n"; return 0.; }

	const int multBin1 = hNorm1D->GetXaxis()->FindBin(mult[0] + 1.e-4);
	const int multBin2 = hNorm1D->GetXaxis()->FindBin(mult[1] - 1.e-4);
	const double fNorm = hNorm1D->Integral(multBin1, multBin2);

	cout <<Form("\nNomalization factor for %s_%s: [%2.1f, %2.1f] (Bins [%i, %i]): %10.9f x 1.e9",
			TRIG, PERC, mult[0], mult[1], multBin1, multBin2, fNorm/1.e9) <<endl;

	hNorm1D->Delete();
	return fNorm;
}//GetNormFactor

//-------------------------------------------------------------
double GetNormFactorFromANC(const char* TRIG, const char* PERC)
{
	const char* F_orig = "./AnalysisResults_data.root";
	const char* SDIR   = "PWG3_D2H_Xic02eXipp13TeV_HM";
	cout <<Form("\nGetNormFactorFromANC - by using %s/%s/ANC_%s_%s", F_orig, SDIR, TRIG, PERC) <<endl;

	//Check if original train output file exists
	TFile* F = TFile::Open(F_orig);
	if (!F || F->IsZombie()) { cout <<"Cannot find the original train output with ANC! Stop.\n"; return 0.; }

	//Check if the counter exists
	AliNormalizationCounter* ANC = (AliNormalizationCounter*)F->Get(Form("%s/ANC_%s_%s", SDIR, TRIG, PERC));
	if (!ANC) {	cout <<Form("Cannot find the %s/ANC_%s_%s! Stop.\n", SDIR, TRIG, PERC);	return 0.; }

	const double fNorm = ANC->GetNEventsForNorm();
	cout <<Form("Normalization factor for %s_%s: %10.9f x 1.e9", TRIG, PERC, fNorm/1.e9) <<endl;

	F->Close();
	return fNorm;
}//GetNormFactorANC

//----------------
TH1D* GetXSection(
		TH1D* hRawXic0, TH1D* hEff, vector<float> pTvec, const double normF,
		const char* Suffix,	bool Show = false
		)
{
	const double BR      = 0.018; //Xic0 -> eXi: 1.8 +- 1.2 (%)
	const double DeltaY  = 1.0;
	const double V0ANDcs = 57.8 * 1000; //V0AND xSec @ 13 TeV: 57.8 mb +- 2.2 (%)
	const double IntL    = normF/V0ANDcs;

	const int npTbins = pTvec.size();
	double pTbins[npTbins];	for (int i=0; i<npTbins; i++) pTbins[i] = pTvec[i];

	//+++++++++++++++++++++++++++++++++++++++++++
	
	//Efficiency corrected total yields
	TH1D* hXic0Yield = (TH1D*)hRawXic0->Clone(Form("%s_yield", hRawXic0->GetName())); hXic0Yield->Reset();
	hXic0Yield->Divide(hRawXic0, hEff);

	//Return histogram
	TH1D* hXic0Xsec = new TH1D(Form("%s_hXsec", Suffix), ";pT", npTbins-1, pTbins);
	for (int i=0; i<npTbins-1; i++)
	{
		if (i==0) continue; //Jun 7 (2021), skip 1st bin (too small statistics)

		const double val = hXic0Yield->GetBinContent(i+1);
		const double err = hXic0Yield->GetBinError(i+1); 

		const double pTbinW = pTbins[i+1] - pTbins[i];
		const double xSecF  = 2 * pTbinW * DeltaY * IntL * BR;

		hXic0Xsec->SetBinContent(i+1, val/xSecF);
		hXic0Xsec->SetBinError  (i+1, err/xSecF);
		//cout <<i <<" " <<val <<" " <<err <<" " <<xSecF <<" " <<val/xSecF <<" " <<err/xSecF <<endl;
	}

	if (Show)
	{
		gStyle->SetOptStat(0);
		gStyle->SetPaintTextFormat("4.3f");

		TCanvas* c1 = new TCanvas("c1", "", 1600/2, 900/2);
		c1->SetName(Form("c1e_xSec_%s", Suffix));
		gPad->SetGrid();
		gPad->SetLogy();

		TH1D* hXsecTemp = (TH1D*)hXic0Xsec->Clone();
		hXsecTemp->SetTitle(Form("Xsec_%s;pT", Suffix));
		hXsecTemp->GetYaxis()->SetRangeUser(5.e-3, 5.e2);
		hXsecTemp->DrawCopy("hist e text0");

		c1->Print(Form("%s.png", c1->GetName())); delete c1;
		hXsecTemp->Delete();
	}

	return hXic0Xsec;
}//GetXSection

#endif //XI0CANAFUNCTION
