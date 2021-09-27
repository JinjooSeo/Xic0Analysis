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

//---------------------------------------------------------
double GetNormFac(const char* inFile, const char* ANC_name)
{
	//Check if original train output file exists
	TFile* F = TFile::Open(inFile);
	if (!F || F->IsZombie()) { cout <<Form("Cannot find file %s! Stop.\n", inFile); return -999.; }

	//Check if the counter exists
	const char* SDIR = "PWG3_D2H_Xic02eXipp13TeV_HM";
	AliNormalizationCounter* ANC = (AliNormalizationCounter*)F->Get(Form("%s/%s", SDIR, ANC_name));
	if (!ANC) {	cout <<Form("Cannot find ANC object %s/%s! Stop.\n", SDIR, ANC_name); return -999.; }

	const double fNorm = ANC->GetNEventsForNorm();
	cout <<Form("Normalization factor from %s: %10.9f x 1.e9\n", ANC_name, fNorm/1.e9);

	F->Close();
	return fNorm;
}//GetNormFac (from ANC)

//----------------------------------------------------------------
double GetRunFrac(TH1F* H1, int RunS, int RunE, bool Show = false) //S: start, E: end
{
	const int LimLo = H1->GetXaxis()->GetBinLowEdge(1);
	const int LimUp = H1->GetXaxis()->GetBinUpEdge(H1->GetNbinsX());
	if (RunS<LimLo || RunE>LimUp) { cout <<Form("Out of range [%i, %i]: stop.\n", LimLo, LimUp); return -999.; }

	//Get fraction of given run range in the given H1
	const int BinS = H1->GetXaxis()->FindBin(RunS + 0.1);
	const int BinE = H1->GetXaxis()->FindBin(RunE - 0.1);
	const double Yield_part = H1->Integral(BinS, BinE);
	const double Yield_full = H1->Integral();

	double Frac = -999;
	if (Yield_part>=0 && Yield_full!=0) Frac = Yield_part/Yield_full;
	if (Show) cout <<Form("Fraction of [%i, %i] in the [%i, %i]: %4.3f\n", RunS, RunE, LimLo, LimUp, Frac);

	return Frac;
}//GetRunFrac

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
	TH1D* H1FONLL = new TH1D(Form("H1_%s", inFile), ";pT", 210, 0, 21); //Finer binning to avoid override

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
		TFile* F_MC, TH1D* H1,
		const double BRFrac, const double NormF, const double V0xs,
		const char* Cut, const char* CutFlag, bool Show = false
		)
{
	//Get CMS 7 TeV Lambda_b and fit it w/ Tsallis
	TH1D* H1_CMSLb = ReadCMSLb("CMSLb.root");//, Show);
	TF1* fTsallis = new TF1("Tsallis", "[0]*x*(pow(1+(sqrt(pow(x,2)+pow(5.619,2))-5.619)/(7.6*1.1),-7.6))", 0,50);
	H1_CMSLb->Fit("Tsallis", "0");
	fTsallis->SetParameters(fTsallis->GetParameters());

	//Get ratio of Lambda_b btw 13/7 TeV
	TH1D* H1_FONLL_7 = ReadFONLL("FONLL_7.txt", H1);//, Show);
	TH1D* H1_FONLL_13 = ReadFONLL("FONLL_13.txt", H1);//, Show);
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

//-------------------
void ApplyPromptFrac()
{

	return;
}//ApplyPromptFrac

//----------------------
void GetWeightingFactor(
		TH1D* H1XS_data, const char* Suffix,
		float xMin = -999., float xMax = 999.,
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

	//Data/MC histograms
	enum {DATA, MC};
	TH1D* H1xs[2];
	for (int a=0; a<2; a++) //data or MC
	{
		if (a==DATA) H1xs[a] = (TH1D*)H1XS_data->Clone("pTW_data");
		else         H1xs[a] = (TH1D*)H1XS_mc->Clone("pTW_mc");

		//Apply x range + normalization by bin width
		for (int b=0; b<H1xs[a]->GetNbinsX(); b++)
		{
			if (H1xs[a]->GetBinCenter(b+1)<xMin || H1xs[a]->GetBinCenter(b+1)>xMax) 
			{
				H1xs[a]->SetBinContent(b+1, 0);
				H1xs[a]->SetBinError(b+1, 0);
			}
			if (fabs(H1xs[a]->GetBinWidth(b+1) - 1) > 1.E-5)
			{
				const double BinC = H1xs[a]->GetBinContent(b+1);
				const double BinE = H1xs[a]->GetBinError(b+1);
				const double BinW = H1xs[a]->GetBinWidth(b+1);
				H1xs[a]->SetBinContent(b+1, BinC/BinW);
				H1xs[a]->SetBinError(b+1, BinE/BinW);
			}
		}//b, loop over bins

		//Normalization by self integral
		const double INTEGRAL = H1xs[a]->Integral();
		H1xs[a]->Scale(1./INTEGRAL);
	}//a

	//Ratio
	enum {CT, UD, DU}; //Center, Updown, and Downup
	TH1D* H1Ratio[3];
	TF1* F1Ratio[3];
	for (int a=0; a<3; a++)
	{
		if (a==CT)
		{
			H1Ratio[a] = (TH1D*)H1xs[DATA]->Clone(Form("pTW_ratio%i", a));
			H1Ratio[a]->Divide(H1xs[MC]);
		}
		else
		{
			H1Ratio[a] = (TH1D*)H1Ratio[CT]->Clone(Form("pTW_ratio%i", a));
			H1Ratio[a]->Sumw2(false);
			H1Ratio[a]->Reset();

			for (int b=0; b<H1Ratio[CT]->GetNbinsX(); b++)
			{
				const float pT  = H1Ratio[CT]->GetBinCenter(b+1);
				const float Val = H1Ratio[CT]->GetBinContent(b+1);
				const float Err = H1Ratio[CT]->GetBinError(b+1); //There's no reason the error to be asymmetric
				if (pT < 4.0)
				{
					if (a==UD) H1Ratio[a]->SetBinContent(b+1, Val + Err);
					if (a==DU) H1Ratio[a]->SetBinContent(b+1, Val - Err);
				}
				else if (pT < 5.0) H1Ratio[a]->SetBinContent(b+1, Val);
				else
				{
					if (a==UD) H1Ratio[a]->SetBinContent(b+1, Val - Err);
					if (a==DU) H1Ratio[a]->SetBinContent(b+1, Val + Err);
				}
			}//b
		}//UD and DU

		F1Ratio[a] = new TF1(Form("F1Ratio_%i", a), "expo", xMin, xMax);
		H1Ratio[a]->Fit(F1Ratio[a]->GetName(), "EQR0");
	}//a

	//+++++++++++++++++++++++++++++++++++++++++++

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

	for (int a=0; a<3; a++) H1Ratio[a]->Delete();
	for (int a=0; a<2; a++) H1xs[a]->Delete();
	H1XS_mc->Delete();
	F->Close();
	return;
}//GetWeightingFactor

#endif //XI0CANAFUNCTION
