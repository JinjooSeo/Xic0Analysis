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
#include <TH1D.h>
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
#include <fstream>
#include <vector>
#include <map>
using namespace std;

//----------------------------------------------------------------------------
double GetNormFac(const char* inFile, const char* ANC_name, bool Show = false)
{
	//Check if original train output file exists
	TFile* F = TFile::Open(inFile);
	if (!F || F->IsZombie()) { cout <<Form("Cannot find file %s! Stop.\n", inFile); return -999.; }

	//Check if the counter exists
	const char* SDIR = "PWG3_D2H_Xic02eXipp13TeV_HM";
	AliNormalizationCounter* ANC = (AliNormalizationCounter*)F->Get(Form("%s/%s", SDIR, ANC_name));
	if (!ANC) {	cout <<Form("Cannot find ANC object %s/%s! Stop.\n", SDIR, ANC_name); return -999.; }

	const double fNorm = ANC->GetNEventsForNorm();
	if (Show) cout <<Form("Normalization factor from %s: %10.9f x 1.e9\n", ANC_name, fNorm/1.e9);

	F->Close();
	return fNorm;
}//GetNormFac (from ANC)

//-------------------------------------------------
double GetV0xSec(int Year, const char* Type = "pp")
{
	double xSec = -999;

	// https://twiki.cern.ch/twiki/bin/viewauth/ALICE/EventNormalization#proton_proton_at_sqrt_s_13_T_AN1
	// Last checked Sep. 27, 2021

	// V0 xSec @ pp 13 TeV
	// pp 2016: 58.44 (mb) +- 1.0/1.6/1.9 (%, correlated, uncorrelated, and total)
	// pp 2017: 58.10      +- 1.2/2.4/2.7
	// pp 2018: 57.52      +- 1.0/1.8/2.1

	// V0 @ p-Pb 5.02 TeV (2013) (* no vdM in 2016 - use 2013 value)
	// p-Pb: 2.09 (b) +- 3.7 (%)
	// Pb-p: 2.12     +- 3.4

	if (Year==2016 && !strcmp(Type, "pp")==true) xSec = 58.44;
	if (Year==2017 && !strcmp(Type, "pp")==true) xSec = 58.10;
	if (Year==2018 && !strcmp(Type, "pp")==true) xSec = 57.52;

	if (Year==2016 && !strcmp(Type, "pPb")==true) xSec = 2.09;
	if (Year==2016 && !strcmp(Type, "Pbp")==true) xSec = 2.12;

	return xSec;
}//GetV0xSec

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

//----------------------------------------------------
TH1D* ReadCMSLb(const char* inFile, bool Show = false)
{
	TFile* F = TFile::Open(inFile);
	if (!F || F->IsZombie()) { cout <<Form("Cannot open the file %s: stop.\n", inFile); assert(false); }
	TDirectoryFile* D = (TDirectoryFile*)F->Get("Table 2");

	TH1D* H1_CMSLb_val = (TH1D*)D->Get("Hist1D_y1");
	TH1D* H1_CMSLb_err = (TH1D*)D->Get("Hist1D_y1_e1");
	H1_CMSLb_val->Scale(0.001/4);
	H1_CMSLb_err->Scale(0.001/4);

	const double pTBins[] = {10, 13, 15, 18, 22, 28, 50};
	const int pTBinsN = sizeof(pTBins)/sizeof(pTBins[0]) - 1;
	TH1D* H1_CMSLb = new TH1D("CMSLb", "", pTBinsN, pTBins);
	for (int a=0; a<pTBinsN; a++)
	{
		H1_CMSLb->SetBinContent(a+1, H1_CMSLb_val->GetBinContent(a+1));
		H1_CMSLb->SetBinError  (a+1, H1_CMSLb_err->GetBinContent(a+1));
	}

	if (Show)
	{
		TCanvas* c1 = new TCanvas("CMSLb", "", 800, 600);
		c1->cd()->SetGrid();
		H1_CMSLb->DrawCopy("hist e");
		c1->Print(Form("%s.png", c1->GetName()));
		delete c1;
	}

	H1_CMSLb_val->Delete();
	H1_CMSLb_err->Delete();
	return H1_CMSLb;
}//ReadCMSLb

//--------------------------------------------------------------
TH1D* ReadFONLL(const char* inFile, TH1D* H1, bool Show = false)
{
	TH1D* H1FONLL = new TH1D(Form("H1_%s", inFile), ";pT", 210, 0, 21); //Finer binning to avoid overriding

	//Read txt file
	int Count = 0;
	ifstream in(inFile);
	if (!in.is_open()) { cout <<Form("Cannot find the file %s: stop.\n", inFile); assert(false); }
	while(in.peek() != EOF)
	{
		if (Count == 0) //Read 1st line and discard it
		{
			string firstLine;
			std::getline(in, firstLine);
			if (Show) cout <<Form("\nReading the file %s...\n", inFile);
		}
		else
		{
			float bPT, bXS;
			in >> bPT >> bXS;
			if (!in.good()) break;

			//Search corresponding bin and assign the value, make sure it does NOT override previous one
			const int tBin = H1FONLL->GetXaxis()->FindBin(bPT);
			H1FONLL->SetBinContent(tBin, bXS);
			H1FONLL->SetBinError(tBin, 0);
			if (Show) cout <<Form("Reading... %3i, %7.4f (-> bin %3i), %7.4f x 1.E6\n", Count, bPT, tBin, bXS/1.E6);
		}
	
		Count++;
	}//While loop: read txt and fill info to TH1

	//Get original binning
	vector<float> pTVec;
	for (int a=0; a<H1->GetNbinsX(); a++)
	{
		pTVec.push_back( H1->GetXaxis()->GetBinLowEdge(a+1) );
		if (a+1 == H1->GetNbinsX()) pTVec.push_back( H1->GetXaxis()->GetBinUpEdge(a+1) );
	}

	//Rebin to match the original binning
	double pTBins[pTVec.size()];
	for (unsigned int a=0; a<pTVec.size(); a++) pTBins[a] = pTVec[a];
	const double* pTBinsFix = pTBins;
	TH1D* H1FONLL_rb = (TH1D*)H1FONLL->Rebin(pTVec.size()-1, Form("%s_rb", H1FONLL->GetName()), pTBinsFix);

	//Normalize values with bin width
	TH1D* H1FONLL_rb_norm = (TH1D*)H1FONLL_rb->Clone(); H1FONLL_rb_norm->Reset();
	for (int a=0; a<H1FONLL_rb->GetNbinsX(); a++)
	{
		const double val0 = H1FONLL_rb->GetBinContent(a+1);
		const double val1 = val0/H1FONLL_rb->GetBinWidth(a+1);
		H1FONLL_rb_norm->SetBinContent(a+1, val1);
	}

	if (Show)
	{
		cout <<Form("\n%s after rebin:\n", inFile);
		for (int a=0; a<H1FONLL_rb->GetNbinsX(); a++)
		{
			const float pTLo = H1FONLL_rb->GetXaxis()->GetBinLowEdge(a+1);
			const float pTUp = H1FONLL_rb->GetXaxis()->GetBinUpEdge(a+1);
			const float Val  = H1FONLL_rb->GetBinContent(a+1);
			cout <<Form("Bin %2i, pT [%2.0f, %2.0f], %8.5f x 1.E6 \n", a+1, pTLo, pTUp, Val/1.E6);
		}

		TString inFileName = inFile;
		inFileName.ReplaceAll(".txt", "");
		TCanvas* c1 = new TCanvas(inFileName, "", 1600, 600); c1->Divide(2, 1); 
		c1->cd(1)->SetGrid(); H1FONLL->DrawCopy("hist e");
		c1->cd(2)->SetGrid(); H1FONLL_rb_norm->DrawCopy("hist e");
		c1->Print(Form("%s.png", c1->GetName()));
		delete c1;
	}

	H1FONLL->Delete();
	H1FONLL_rb->Delete();
	return H1FONLL_rb_norm;
}//ReadFONLL

//-------------------
void ApplyBottomCorr(
		TH1D* H1, TFile* F_MC,
		const double BRFrac, const double NormF, const double V0xs,
		const char* Cut, const char* CutFlag, bool Show = false
		)
{
	//Get CMS 7 TeV Lambda_b and fit it w/ Tsallis
	TH1D* H1_CMSLb = ReadCMSLb("./input/HEPData-ins1113442-v1-Table_2.root");
	TF1* fTsallis = new TF1("Tsallis", "[0]*x*(pow(1+(sqrt(pow(x,2)+pow(5.619,2))-5.619)/(7.6*1.1),-7.6))", 0,50);
	H1_CMSLb->Fit("Tsallis", "0");
	fTsallis->SetParameters(fTsallis->GetParameters());

	//Get ratio of Lambda_b btw 13/7 TeV
	TH1D* H1_FONLL_7 = ReadFONLL("./input/FONLL-Bmeson-dsdpt-sqrts7000-20GeV.txt", H1);
	TH1D* H1_FONLL_13 = ReadFONLL("./input/FONLL-Bmeson-dsdpt-sqrts13000-20GeV.txt", H1);
	TH1D* H1_FONLL_13to7 = (TH1D*)H1_FONLL_13->Clone("FONLL_13to7");
	H1_FONLL_13to7->Divide(H1_FONLL_7);

	//Get 13 TeV Lambda_b by scale up 7 TeV CMSLb
	TH1D* H1_Lb_13 = (TH1D*)H1->Clone("Lb_13"); H1_Lb_13->Reset();
	for (int a=0; a<H1_Lb_13->GetNbinsX(); a++)
	{
		const double Val = fTsallis->Eval(H1_Lb_13->GetBinCenter(a+1)) * H1_FONLL_13to7->GetBinContent(a+1);
		H1_Lb_13->SetBinContent(a+1, Val);
	}

	//Get 13 TeV Xib by multiplying Xic0 -> eXi BR
	TH1D* H1_Lb_13_BRFscaled = (TH1D*)H1_Lb_13->Clone("Lb_13_BRFscaled");
	H1_Lb_13_BRFscaled->Scale(BRFrac);

	//Estimate Xib yield
	TH1D* H1_Xib_13 = (TH1D*)H1_Lb_13_BRFscaled->Clone("Xib_13"); H1_Xib_13->Reset();
	for (int a=0; a<H1_Lb_13_BRFscaled->GetNbinsX(); a++)
	{
		const double BinC = H1_Lb_13_BRFscaled->GetBinContent(a+1);
		const double BinW = H1_Lb_13_BRFscaled->GetBinWidth(a+1);
		const double Lumi = (NormF/V0xs)/1.E3; //1.E3: scale matching
		//const double Lumi = 1.86437e+09 / (57.8*1000000); //!?
		H1_Xib_13->SetBinContent(a+1, BinC * BinW * Lumi * 2); //2: charge conjugate
	}

	//Apply efficiencty from MC
	TH1D* H1_Xib_gen  = (TH1D*)F_MC->Get("XibGen05"); //# of generated Xib
	TH1D* H1_Xib_reco = (TH1D*)F_MC->Get(Form("hMCRecoLevXib_%s_%s", Cut, CutFlag)); //# of reco Xib
	TH1D* H1_Xib_eff  = (TH1D*)H1_Xib_reco->Clone(Form("Xib_eff_%s_%s", Cut, CutFlag)); H1_Xib_eff->Reset();
	H1_Xib_eff->Divide(H1_Xib_reco, H1_Xib_gen, 1, 1, "b");
	H1_Xib_13->Multiply(H1_Xib_eff);

	//Convert Xib -> eXi
	TH1D* H1_eXi_reco = (TH1D*)F_MC->Get(Form("hMCRecoLevPairXib_%s_%s", Cut, CutFlag)); //# of eXi from Xib
	TH2D* H2_Xib_eXi_RM = (TH2D*)F_MC->Get(Form("hRPM_%s_%s_Xib", Cut, CutFlag)); //Xib <-> eXi response matrix
	RooUnfoldResponse Response(H1_Xib_reco, H1_eXi_reco, H2_Xib_eXi_RM);
	RooUnfoldBinByBin Unfolding(&Response, H1_Xib_13);
	TH1D* H1_eXiFromXib = (TH1D*)Unfolding.Hreco();

	//Add eXi from Xib to given eXi pair (RS-WS)
	H1->Add(H1_eXiFromXib);

	if (Show)
	{
		gStyle->SetOptStat(0);
		TCanvas* c1 = new TCanvas("c1b_Bcorr", "", 1600, 900);
		c1->Divide(2, 2);

		c1->cd(1)->SetGrid();
		H1_CMSLb->SetTitle("CMSLb;pT");
		H1_CMSLb->DrawCopy("hist e text0");
		fTsallis->SetLineColor(2);
		fTsallis->SetLineStyle(2);
		fTsallis->Draw("same");

		c1->cd(2)->SetGrid();
		H1_FONLL_13to7->SetTitle("FONLL_13/7;pT");
		H1_FONLL_13to7->DrawCopy("hist e text0");

		c1->cd(3)->SetGrid();
		H1_Xib_13->SetTitle("Xib_13TeV (BRFrac scaled + Eff corr);pT");
		H1_Xib_13->DrawCopy("hist e text0");

		c1->cd(4)->SetGrid();
		H1_eXiFromXib->SetTitle("eXi from Xib (unfolded);pT");
		H1_eXiFromXib->DrawCopy("hist e text0");

		c1->Print(Form("%s.png", c1->GetName()));
		delete c1;
	}

	fTsallis->Delete();
	H1_CMSLb->Delete();
	H1_FONLL_7->Delete();
	H1_FONLL_13->Delete();
	H1_FONLL_13to7->Delete();

	H1_Lb_13->Delete();
	H1_Lb_13_BRFscaled->Delete();
	H1_Xib_13->Delete();

	H1_Xib_gen->Delete();
	H1_Xib_reco->Delete();
	H1_Xib_eff->Delete();
	H1_eXi_reco->Delete();
	H2_Xib_eXi_RM->Delete();

	return;
}//ApplyBottomCorr

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

//--------------------------------------------------------------------------------------
TH1D* GetXS(TH1D* H1, double normF, double V0xSec, double BRXic0eXi, const char* Suffix)
{
	TH1D* H1R = (TH1D*)H1->Clone(Form("Xic0XS_%s", Suffix));
	H1R->SetTitle("");
	H1R->Reset();

	for (int a=0; a<H1->GetNbinsX(); a++)
	{
		const double Val = H1->GetBinContent(a+1);
		const double Err = H1->GetBinError(a+1);

		const double BinW   = H1->GetBinWidth(a+1);
		const double CConj  = 2.; //e-Xi pair charge conjugate
		const double DeltaY = 1.;
		const double IntL   = normF/V0xSec;
		const double xSecF  = BinW * CConj * DeltaY * IntL * BRXic0eXi;

		H1R->SetBinContent(a+1, Val/xSecF);
		H1R->SetBinError  (a+1, Err/xSecF);
	}

	return H1R;
}//GetXS

//----------------------------------------------------------------------------------
void ApplyPromptFrac(TH1D* H1Xsec, TFile* F_MC, double BRXic0eXi, bool Show = false)
{
	enum {AVG, MAX, MIN};
	enum {FDDOWN, PROMPT};
	const float LcFddownScaleF = 1.E-6/(0.0628 * 20); //pb to ub and divide by Lc2pKpi BR (PDG2020) + binning

	/*
	   Procedures

	   1. Preparation:

	   		a. Two Xic0 efficiency: feed-down and prompt
			b. Three Lc distributions from external source: central (avg), max, and min
			c. Get Xic0 - Lc ratio by using pol1 fit on Xic0-Lc raito from external source
			d. Xic0 feed-down distributions by using 'Lc' x 'Xic0-Lc ratio'
				d-1) 3 variations by Lc
				d-2) 3 variations by Xic0-Lc ratio

		2. Get prompt fraction:

			a. Get ' inclusive Xic0 x inclusive eff ', from previous step
			b. Get ' feed-dwon Xic0 x feed-down eff '
			c. Get ratio of feed-down/inclusive by using 2-a. and 2-b.
			d. Get prompt fraction: " 1.0 - c., pT bin by pT bin "
	*/

	//Get original binning from xSec histogram
	vector<float> pTVec;
	for (int a=0; a<H1Xsec->GetNbinsX(); a++)
	{
		pTVec.push_back( H1Xsec->GetXaxis()->GetBinLowEdge(a+1) );
		if (a+1 == H1Xsec->GetNbinsX()) pTVec.push_back( H1Xsec->GetXaxis()->GetBinUpEdge(a+1) );
	}
	double pTBins[pTVec.size()];
	for (unsigned int a=0; a<pTVec.size(); a++) pTBins[a] = pTVec[a];
	const double* pTBinsFix = pTBins;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//1 - a. Xic0 efficiency by MC, separated by feeddown or prompt 
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Originally uses MC unweighted - why? (Tested: MCraw = MC_MB_0to100 = MC_HMV0_0to0p1)
	//TFile* F_MC = TFile::Open("out_MC_raw.root");
	TH1D* H1_Xic0EffSep[2]; //Jinjoo: hFeeddownEff, hPromptEff
	for (int a=0; a<2; a++)
	{
		const char* Histo_d = (a==FDDOWN)?"hMCGenFeeddowmXic0_woW":"hMCGenPromptXic0_woW";
		const char* Histo_n = (a==FDDOWN)?"hnonprompt":"hprompt";
		TH1D* H1_Xic0Eff_d = (TH1D*)F_MC->Get(Histo_d);
		TH1D* H1_Xic0Eff_n = (TH1D*)F_MC->Get(Histo_n);

		H1_Xic0EffSep[a] = (TH1D*)H1_Xic0Eff_n->Clone(Form("PF_Xic0Eff_%i", a)); H1_Xic0EffSep[a]->Reset();
		H1_Xic0EffSep[a]->SetTitle(Form("Xic0Eff_%s", a==FDDOWN?"feed-down":"prompt"));
		H1_Xic0EffSep[a]->Divide(H1_Xic0Eff_n, H1_Xic0Eff_d, 1, 1, "b");

		H1_Xic0Eff_d->Delete();
		H1_Xic0Eff_n->Delete();
	}//a

	//1 - b. Feeddown Lc
	//+++++++++++++++++++++++++++++++++++++++++++

	TFile* F_MC_LcFddown = TFile::Open("./input/DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020.root");
	vector<const char*> v_LcFddown = 
	{
		"hLcpkpifromBpred_central_corr",
		"hLcpkpifromBpred_max_corr", //hFeeddownLcFONLL + FFMax
		"hLcpkpifromBpred_min_corr", //hFeeddownLcFONLL + FFMin
	};
	const int n_LcFddown = v_LcFddown.size();

	TH1D* H1_LcFddown[n_LcFddown]; //Jinjoo: hFeeddownLc
	for (int a=0; a<n_LcFddown; a++)
	{
		TH1D* H1_LcFddown_orig = (TH1D*)F_MC_LcFddown->Get(v_LcFddown[a]);
		TH1D* H1_LcFddown_rb = (TH1D*)H1_LcFddown_orig->Rebin(pTVec.size()-1, Form("PF_LcFddown_%i",a), pTBinsFix);
		H1_LcFddown_rb->Scale(LcFddownScaleF);

		H1_LcFddown[a] = (TH1D*)H1_LcFddown_rb->Clone(); H1_LcFddown[a]->Reset();
		for (unsigned int b=0; b<pTVec.size(); b++)
		{
			const double val = H1_LcFddown_rb->GetBinContent(b+1) / H1_LcFddown_rb->GetBinWidth(b+1);
			const double err = H1_LcFddown_rb->GetBinError  (b+1) / H1_LcFddown_rb->GetBinWidth(b+1);
			H1_LcFddown[a]->SetBinContent(b+1, val);
			H1_LcFddown[a]->SetBinError  (b+1, err);
		}//b, pTBins

		H1_LcFddown_orig->Delete();
		H1_LcFddown_rb  ->Delete();
	}//a, LcFddown

	//1 - c. Xic0 to Lc ratio, get pol1 parameters vs. pT
	//+++++++++++++++++++++++++++++++++++++++++++++++++++

	TFile* F_MC_Xic0ToLc = TFile::Open("./input/Xic0toLc_pp13TeV_new.root");
	TH1D* H1_Xic0ToLc = (TH1D*)F_MC_Xic0ToLc->Get("hRatio_Xic0toLc"); //Jinjoo: hRatioXicLc

	const float FitX0 = H1_Xic0ToLc->GetXaxis()->GetBinLowEdge(1);
	const float FitX1 = H1_Xic0ToLc->GetXaxis()->GetBinUpEdge(H1_Xic0ToLc->GetNbinsX());

	TF1* F1_Xic0ToLc = new TF1("F1_Xic0ToLc", "pol1");//, FitX0, FitX1); //Jinjoo: fFitFunction
	TF1* F1_Xic0ToLcMax = new TF1("F1_Xic0ToLcMax", "pol1");
	TF1* F1_Xic0ToLcMin = new TF1("F1_Xic0ToLcMin", "pol1");

	H1_Xic0ToLc->Fit(F1_Xic0ToLc->GetName(), "QL0", "kBlack", FitX0, FitX1);
	F1_Xic0ToLcMax->SetParameters(F1_Xic0ToLc->GetParameter(0)*2., F1_Xic0ToLc->GetParameter(1)*2.);
	F1_Xic0ToLcMin->SetParameters(F1_Xic0ToLc->GetParameter(0)/2., F1_Xic0ToLc->GetParameter(1)/2.);

	//1-d-1) Feeddown Xic0, variation by Lc + fixed Xic0ToLc ratio
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	vector<const char*> v_Xic0FddownVarLc = 
	{
		"hFeeddownXic0_LcAvg",
		"hFeeddownXic0_LcMax",
		"hFeeddownXic0_LcMin"
	};
	const int n_Xic0FddownVarLc = v_Xic0FddownVarLc.size();

	TH1D* H1_Xic0FddownVarLc[n_Xic0FddownVarLc]; //Jinjoo: hFeeddownXicVarLc
	for (int a=0; a<n_Xic0FddownVarLc; a++)
	{
		H1_Xic0FddownVarLc[a] = (TH1D*)H1_LcFddown[a]->Clone(v_Xic0FddownVarLc[a]);
		H1_Xic0FddownVarLc[a]->Reset();

		for (int b=0; b<H1_Xic0FddownVarLc[a]->GetNbinsX(); b++)
		{
			const float LcBinVal    = H1_LcFddown[a]->GetBinContent(b+1);
			const float LcBinErr    = H1_LcFddown[a]->GetBinError(b+1);
			const float LcBinCenter = H1_LcFddown[a]->GetBinCenter(b+1);
			const float Xic0ToLcR   = F1_Xic0ToLc->Eval(LcBinCenter);
			H1_Xic0FddownVarLc[a]->SetBinContent(b+1, LcBinVal * Xic0ToLcR);
			H1_Xic0FddownVarLc[a]->SetBinError  (b+1, LcBinErr * Xic0ToLcR);
		}//b
	}//a, Xic0FddownVarLc

	//1-d-2) Feeddown Xic0, fixed Lc + variation by Xic0ToLc ratio
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	vector<const char*> v_Xic0FddownVarR = 
	{
		"hFeeddownXic0_RatioAvg",
		"hFeeddownXic0_RatioMax",
		"hFeeddownXic0_RatioMin"
	};
	const int n_Xic0FddownVarR = v_Xic0FddownVarR.size();

	TH1D* H1_Xic0FddownVarR[n_Xic0FddownVarR]; //Jinjoo: hFeeddownXicVarRatio
	for (int a=0; a<n_Xic0FddownVarR; a++)
	{
		H1_Xic0FddownVarR[a] = (TH1D*)H1_LcFddown[AVG]->Clone(v_Xic0FddownVarR[a]);
		H1_Xic0FddownVarR[a]->Reset();

		for (int b=0; b<H1_Xic0FddownVarR[a]->GetNbinsX(); b++)
		{
			const float LcBinVal    = H1_LcFddown[AVG]->GetBinContent(b+1);
			const float LcBinErr    = H1_LcFddown[AVG]->GetBinError(b+1);
			const float LcBinCenter = H1_LcFddown[AVG]->GetBinCenter(b+1);
			float Xic0ToLcR = 9999;
			if      (a==AVG) Xic0ToLcR = F1_Xic0ToLc->Eval(LcBinCenter);
			else if (a==MAX) Xic0ToLcR = F1_Xic0ToLcMax->Eval(LcBinCenter);
			else if (a==MIN) Xic0ToLcR = F1_Xic0ToLcMin->Eval(LcBinCenter);
			H1_Xic0FddownVarR[a]->SetBinContent(b+1, LcBinVal * Xic0ToLcR);
			H1_Xic0FddownVarR[a]->SetBinError  (b+1, LcBinErr * Xic0ToLcR);
		}
	}//

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//Xic0 prompt fraction, variation by Lc
	//+++++++++++++++++++++++++++++++++++++++++++

	TH1D* H1_Xic0PromptVarLc[n_Xic0FddownVarLc];
	for (int a=0; a<n_Xic0FddownVarLc; a++)
	{
		H1_Xic0PromptVarLc[a] = (TH1D*)H1_LcFddown[a]->Clone(v_Xic0FddownVarLc[a]);
		H1_Xic0PromptVarLc[a]->SetTitle(Form("Xic0 prompt frac, variation by base Lc, %i", a));
		H1_Xic0PromptVarLc[a]->Reset();

		for (int b=0; b<H1_Xic0PromptVarLc[a]->GetNbinsX(); b++)
		{
			H1_Xic0PromptVarLc[a]->SetBinContent(b+1, 1);
			H1_Xic0PromptVarLc[a]->SetBinError  (b+1, 0);
		}//b

		TH1D* H1_xSecEff_inc = (TH1D*)H1Xsec->Clone(Form("%s_temp%i", H1Xsec->GetName(), a)); //xSec_inc * eff_inc
		TH1D* H1_xSecEff_fdn = (TH1D*)H1_Xic0FddownVarLc[a]->Clone(Form("%s_temp%i", v_Xic0FddownVarLc[a], a));
		TH1D* H1_xSecEff_fdnFrac = (TH1D*)H1_xSecEff_fdn->Clone(Form("xSecEff_varLc_fdnFrac%i", a));
		H1_xSecEff_fdnFrac->Reset();
		H1_xSecEff_fdnFrac->Divide(H1_xSecEff_fdn, H1_xSecEff_inc, 1, 1, "b");

		H1_Xic0PromptVarLc[a]->Add(H1_xSecEff_fdnFrac, -1); //Jinjoo: 1.73 or 1, NEED TO CONFIRM

		H1_xSecEff_inc->Delete();
		H1_xSecEff_fdn->Delete();
		H1_xSecEff_fdnFrac->Delete();
	}//a

	//Xic0 prompt fraction, variation by Xic0ToLc ratio (systematic?)
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	TH1D* H1_Xic0PromptVarR[n_Xic0FddownVarR];
	for (int a=0; a<n_Xic0FddownVarR; a++)
	{
		H1_Xic0PromptVarR[a] = (TH1D*)H1_LcFddown[a]->Clone(v_Xic0FddownVarR[a]);
		H1_Xic0PromptVarR[a]->SetTitle(Form("Xic0 prompt frac, variation by Xic0/Lc ratio, %i", a));
		H1_Xic0PromptVarR[a]->Reset();

		for (int b=0; b<H1_Xic0PromptVarR[a]->GetNbinsX(); b++)
		{
			H1_Xic0PromptVarR[a]->SetBinContent(b+1, 1);
			H1_Xic0PromptVarR[a]->SetBinError  (b+1, 0);
		}//b

		TH1D* H1_xSecEff_inc = (TH1D*)H1Xsec->Clone(Form("%s_temp%i", H1Xsec->GetName(), a)); //xSec_inc * eff_inc
		TH1D* H1_xSecEff_fdn = (TH1D*)H1_Xic0FddownVarR[a]->Clone(Form("%s_temp%i", v_Xic0FddownVarR[a], a));
		TH1D* H1_xSecEff_fdnFrac = (TH1D*)H1_xSecEff_fdn->Clone(Form("xSecEff_varLc_fdnFrac%i", a));
		H1_xSecEff_fdnFrac->Reset();
		H1_xSecEff_fdnFrac->Divide(H1_xSecEff_fdn, H1_xSecEff_inc, 1, 1, "b");

		H1_Xic0PromptVarR[a]->Add(H1_xSecEff_fdnFrac, -1); //Jinjoo: 1.73 or 1, NEED TO CONFIRM

		H1_xSecEff_inc->Delete();
		H1_xSecEff_fdn->Delete();
		H1_xSecEff_fdnFrac->Delete();
	}//a

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	if (Show)
	{
		TCanvas* c1 = new TCanvas("c1d_0", "", 1600, 1000); c1->Divide(3, 2);
		c1->cd(1); H1_Xic0EffSep[0]->DrawCopy("hist e text90");
		c1->cd(2); H1_Xic0EffSep[1]->DrawCopy("hist e text90");
		c1->cd(3); H1_Xic0ToLc->DrawCopy("hist e text90");
		c1->cd(4); H1_LcFddown[0]->DrawCopy("hist e text90");
		c1->cd(5); H1_LcFddown[1]->DrawCopy("hist e text90");
		c1->cd(6); H1_LcFddown[2]->DrawCopy("hist e text90");

		TCanvas* c2 = new TCanvas("c1d_1", "", 1600, 1000); c2->Divide(3, 2);
		c2->cd(1); H1_Xic0FddownVarLc[0]->DrawCopy("hist e");
		c2->cd(2); H1_Xic0FddownVarLc[1]->DrawCopy("hist e");
		c2->cd(3); H1_Xic0FddownVarLc[2]->DrawCopy("hist e");
		c2->cd(4); H1_Xic0FddownVarR[0]->DrawCopy("hist e");
		c2->cd(5); H1_Xic0FddownVarR[1]->DrawCopy("hist e");
		c2->cd(6); H1_Xic0FddownVarR[2]->DrawCopy("hist e");

		TCanvas* c3 = new TCanvas("c1d_2", "", 1600, 1000); c3->Divide(3, 2);
		c3->cd(1); H1_Xic0PromptVarLc[0]->DrawCopy("hist e text90");
		c3->cd(2); H1_Xic0PromptVarLc[1]->DrawCopy("hist e text90");
		c3->cd(3); H1_Xic0PromptVarLc[2]->DrawCopy("hist e text90");
		c3->cd(4); H1_Xic0PromptVarR[0]->DrawCopy("hist e text90");
		c3->cd(5); H1_Xic0PromptVarR[1]->DrawCopy("hist e text90");
		c3->cd(6); H1_Xic0PromptVarR[2]->DrawCopy("hist e text90");

		//c1->Print(Form("%s.png", c1->GetName()));
		//c2->Print(Form("%s.png", c2->GetName()));
		//c3->Print(Form("%s.png", c3->GetName()));
	}

	//APPLY PROMPT FRACTION TO XSEC
	cout <<"Applying prompt fraction...\n";
	H1Xsec->Multiply(H1_Xic0PromptVarLc[AVG]);

	//Cleanup
	F_MC_LcFddown->Close();
	F_MC_Xic0ToLc->Close();
	return;
}//ApplyPromptFrac

//----------------------
void GetWeightingFactor(
		TH1D* H1XS_data, const char* Suffix,
		float xMin = 0., float xMax = 20.,
		bool Show = true
		)
{
	TFile* F = TFile::Open("./out_MC_raw.root", "READ");
	if (!F || F->IsZombie()) { cout <<Form("Cannot open the file 'out_MC_raw.root'! Stop.\n"); return; }

	TH1D* H1XS_mc = (TH1D*)F->Get("hMCGenInclusiveXic0_woW");
	if (!H1XS_mc) { cout <<"Cannot find 'hMCGenInclusiveXic0_woW'! Stop.\n"; return; }

	const int maxBin = H1XS_data->GetNbinsX();
	if ( (H1XS_data->GetNbinsX() != H1XS_mc->GetNbinsX()) ||
		 (H1XS_data->GetXaxis()->GetBinLowEdge(1) != H1XS_mc->GetXaxis()->GetBinLowEdge(1)) ||
		 (H1XS_data->GetXaxis()->GetBinUpEdge(maxBin) != H1XS_mc->GetXaxis()->GetBinUpEdge(maxBin)) )
	{
		cout <<"Histograms NOT match! Stop.\n";
		return;
	}

	//+++++++++++++++++++++++++++++++++++++++++++

	enum {DT, MC};
	TH1D* H1xs[2];
	for (int a=0; a<2; a++) //data or MC
	{
		if (a==DT) H1xs[a] = (TH1D*)H1XS_data->Clone("pTW_data"); //Already normalized by bin width in xs routine
		else
		{
			H1xs[a] = (TH1D*)H1XS_mc->Clone("pTW_mc");

			//Normalization by bin width
			for (int b=0; b<H1xs[a]->GetNbinsX(); b++)
			{
				const double BinC = H1xs[a]->GetBinContent(b+1);
				const double BinE = H1xs[a]->GetBinError(b+1);
				const double BinW = H1xs[a]->GetBinWidth(b+1);
				H1xs[a]->SetBinContent(b+1, BinC/BinW);
				H1xs[a]->SetBinError  (b+1, BinE/BinW);
			}
		}//a==MC

		for (int b=0; b<H1xs[a]->GetNbinsX(); b++)
		{
			const float xLo = H1xs[a]->GetXaxis()->GetBinLowEdge(b+1);
			const float xUp = H1xs[a]->GetXaxis()->GetBinUpEdge(b+1);
			if (xLo<xMin || xUp>xMax)
			{
				H1xs[a]->SetBinContent(b+1, 0);
				H1xs[a]->SetBinError  (b+1, 0);
			}
		}//b


		/*
		//Set range
		if (xMin>0. || xMax<20.)
		{
			const int xMinB = H1xs[a]->GetXaxis()->FindBin(xMin + 1.E-3);
			const int xMaxB = H1xs[a]->GetXaxis()->FindBin(xMax - 1.E-3);
			const float xMinV = H1xs[a]->GetXaxis()->GetBinLowEdge(xMinB);
			const float xMaxV = H1xs[a]->GetXaxis()->GetBinUpEdge(xMaxB);
			H1xs[a]->GetXaxis()->SetRangeUser(xMinV, xMaxV);
		}
		*/

		//Normalization by self integral: checked the value changes w/ range
		const double INTEGRAL = H1xs[a]->Integral();
		H1xs[a]->Scale(1./INTEGRAL);
	}//a

	//Ratio
	enum {CT, UD, DU}; //Center, Updown, and Downup
	TH1D* H1Ratio[3];
	TF1*  F1Ratio[3];
	for (int a=0; a<3; a++)
	{
		if (a==CT)
		{
			H1Ratio[a] = (TH1D*)H1xs[DT]->Clone(Form("pTW_ratio%i", a));
			H1Ratio[a]->Divide(H1xs[MC]);
		}
		else //UD, DU
		{
			const int pTBinN = H1Ratio[CT]->GetNbinsX(); double pTBin[pTBinN];
			for (int x=0; x<pTBinN+1; x++) pTBin[x] = H1Ratio[CT]->GetXaxis()->GetBinLowEdge(x+1);

			H1Ratio[a] = new TH1D();
			H1Ratio[a]->SetBins(pTBinN, pTBin);
			H1Ratio[a]->SetName(Form("pTW_ratio%i", a));

			for (int b=0; b<pTBinN; b++)
			{
				const float pT  = H1Ratio[CT]->GetBinCenter (b+1);
				const float Val = H1Ratio[CT]->GetBinContent(b+1);
				const float Err = H1Ratio[CT]->GetBinError  (b+1);

				if (pT < 4.0)
				{
					if (a==UD) H1Ratio[a]->SetBinContent(b+1, Val + Err);
					if (a==DU) H1Ratio[a]->SetBinContent(b+1, Val - Err);
				}
				else if (pT < 5.0) //Anchor point: 4 < pT < 5
				{
					H1Ratio[a]->SetBinContent(b+1, Val);
				}
				else
				{
					if (a==UD) H1Ratio[a]->SetBinContent(b+1, Val - Err);
					if (a==DU) H1Ratio[a]->SetBinContent(b+1, Val + Err);
				}
			}//b
		}//UD, DU

		//Fit
		F1Ratio[a] = new TF1(Form("F1Ratio_%i", a), "expo", xMin, xMax);
		H1Ratio[a]->Fit(F1Ratio[a]->GetName(), "EQR0");
	}//a

	//+++++++++++++++++++++++++++++++++++++++++++

#if 1
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
            if (a==0) H1xs[a]->GetXaxis()->SetRangeUser(xMin, xMax);
            if (a==0) H1xs[a]->GetYaxis()->SetRangeUser(1.e-5, 3.0);
            if (a==0) H1xs[a]->SetTitle(Form("%s;pT;dN/NdPT", Suffix));
            H1xs[a]->SetStats(false);
            H1xs[a]->SetLineColor(a+1);
            H1xs[a]->SetMarkerColor(a+1);
            H1xs[a]->SetMarkerSize(1.2);
            H1xs[a]->SetMarkerStyle(a*4 + 20);
            H1xs[a]->DrawCopy(a==0?"hist ep":"hist ep same");
            L1->AddEntry(H1xs[a], a==0?"data":"MC", "lp");
        }
        L1->Draw("same");

        c1->cd(2);
        const int RCOLOR[3] = {1, 2, 4};
        TLegend* L2 = new TLegend(0.45, 0.5, 0.85, 0.85);
        L2->SetHeader("Anchor point: 4 < pT < 5", "C");
        L2->SetMargin(0.3);
        for (int a=0; a<3; a++)
        {
            if (a==0) H1Ratio[a]->GetXaxis()->SetRangeUser(xMin, xMax);
            if (a==0) H1Ratio[a]->GetYaxis()->SetRangeUser(-0.5, 3.);
            if (a==0) H1Ratio[a]->SetTitle(Form("%s;pT;Weighting factor", Suffix));
            H1Ratio[a]->SetStats(false);
            H1Ratio[a]->SetLineColor(RCOLOR[a]);
            H1Ratio[a]->SetMarkerColor(RCOLOR[a]);
            H1Ratio[a]->SetMarkerStyle(a==0?20:24);
            H1Ratio[a]->SetMarkerSize(1.2);
            H1Ratio[a]->DrawCopy(a==0?"pe":"p same");
            F1Ratio[a]->SetLineColor(RCOLOR[a]);
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
#endif

	for (int a=0; a<3; a++) H1Ratio[a]->Delete();
	for (int a=0; a<2; a++) H1xs[a]->Delete();
	H1XS_mc->Delete();
	F->Close();
	return;
}//GetWeightingFactor

//--------------------------------------------------------------------
TH1D* Cleanup(TH1D* H1, const float xMin = 0., const float xMax = 20.)
{
	TH1D* H1R = (TH1D*)H1->Clone(Form("%s_CLEANED", H1->GetName())); H1R->Reset();

	for (int a=0; a<H1->GetNbinsX(); a++)
	{
		const float xPosLo = H1->GetXaxis()->GetBinLowEdge(a+1);
		const float xPosUp = H1->GetXaxis()->GetBinUpEdge(a+1);
		if (xPosLo<xMin || xPosUp>xMax) continue;

		const double xVal = H1->GetBinContent(a+1);
		const double xErr = H1->GetBinError(a+1);
		H1R->SetBinContent(a+1, xVal);
		H1R->SetBinError  (a+1, xErr);
	}

	return H1R;
}//Cleanup

#endif //XI0CANAFUNCTION
