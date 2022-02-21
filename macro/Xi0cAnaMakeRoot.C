#ifndef XI0CANAMAKEROOT
#define XI0CANAMAKEROOT

/*
	Make 2nd ROOT files for Xi0c analysis, originally written by J. Seo

	Chong Kim
	Inha Univ. / Pusan National Univ.
	kimc@cern.ch
*/

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
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <iterator>
using namespace std;

TH1F* MakeTH1(const Char_t *name, Int_t NumOfBin, Double_t* binning);
TH2F* MakeTH2(const Char_t *name, Int_t NumOfBin, Double_t* binning);

void nSigmaPlot(TFile* F, const char* SDIR); //nSigmaTOF and TPC plot
void XiMassvsPt(TFile* F, const char* SDIR); //Xi mass distribution for various pT bin
void XiCutDistribution(TFile *F, const char* SDIR); //look Xi cut value distribution from prompt an feeddown Xic0
void eXiPairTree(TFile* F, bool IsMC, const char* SDIR, const char* TRIG, double* MultPerc, double* WFitPar);

//----------------------------------------------------------------------------------
void Xi0cAnaMakeRoot(const char* inFile = "/Users/jinjoo/Desktop/Xic0/data/Download/train_0201/AnalysisResultsData.root", bool IsMC = false)
{
    double MultPerc[2] = {  0, 100}; //Multiplicity percentile
    //double WFitPar[2]  = {0.889618,  -0.329188}; //default
    double WFitPar[2]  = {1.0, 1.0};
    //double WFitPar[2]  = { 1.43224e+00, -4.44314e-01}; //var1
    //double WFitPar[2]  = { 1.99904e-01, -1.98923e-01}; //var2

    const char* TRIG = "MB"; //MB, HMV0, HMSPD, or HMOR
	const char* SDIR = "PWG3_D2H_Xic02eXipp13TeV_HM"; //Sub directory in the file, !

    const char* PERC = Form("%ito%i", (int)MultPerc[0], (int)MultPerc[1]);
    const char* TYPE = (IsMC)?"MC":"data";
	if (IsMC) TYPE = Form("%s%s", TYPE, (fabs(WFitPar[0]-1.0)<1.e-5 && fabs(WFitPar[1]-1.0)<1.e-5)?"raw":"wgt");

	cout <<"\nGenerating ROOT file in following setup:" <<endl;
	cout <<Form("- Type: %s", TYPE) <<endl;
	cout <<Form("- Trigger: %s", TRIG) <<endl;
	cout <<Form("- Multiplicity percentile: [%2.1f, %2.1f]", MultPerc[0], MultPerc[1]) <<endl;
	cout <<Form("- Weight fit parameters: [%5.4f, %5.4f] \n", WFitPar[0], WFitPar[1]) <<endl;

	//+++++++++++++++++++++++++++++++++

    TFile* F = TFile::Open(inFile);
    if (!F || F->IsZombie()) { cout <<Form("Cannot open %s!\n", inFile); return; }
    TFile* G = new TFile(Form("out_%s_%s_%s.root", TYPE, TRIG, PERC), "recreate");

    //nSigmaPlot(F, SDIR);
    //XiMassvsPt(F, SDIR);
    //XiCutDistribution(F, SDIR);
    eXiPairTree(F, IsMC, SDIR, TRIG, MultPerc, WFitPar);

    G->Write();
    G->Close();
    F->Close();

    return;
}//Main

//------------------------------------------------------------------
TH1F* MakeTH1(const Char_t *name, Int_t NumOfBin, Double_t* binning)
{
	TH1F* H1 = new TH1F(name, "", NumOfBin, binning);
	H1->Sumw2();
	return H1;
}

//------------------------------------------------------------------
TH2F* MakeTH2(const Char_t *name, Int_t NumOfBin, Double_t* binning)
{
	TH2F* H2 = new TH2F(name, "", NumOfBin, binning, NumOfBin, binning);
	H2->Sumw2();
	return H2;
}

//-----------------------------------------
void nSigmaPlot(TFile* F, const char* SDIR)
{
	//Link hist object
	TObject* hist = F->Get(Form("%s/histogram", SDIR));

    TH2D* hTPC = (TH2D*) hist->FindObject("nSigmaTPCvsPt")->Clone("nSigmaTPCvsPt");
    TH2D* hTOF = (TH2D*) hist->FindObject("nSigmaTOFvsPt")->Clone("nSigmaTOFvsPt");

    TH1D * hTPC05_1 = hTPC->ProjectionY("hTPC05_1", 1 ,2);
    TH1D * hTPC1_2  = hTPC->ProjectionY("hTPC1_2",  3, 4);
    TH1D * hTPC2_3  = hTPC->ProjectionY("hTPC2_3",  5, 6);
    TH1D * hTPC3_4  = hTPC->ProjectionY("hTPC3_4",  7, 8);
    TH1D * hTPC4_5  = hTPC->ProjectionY("hTPC4_5", 9, 10);
    return;
}

//-----------------------------------------
void XiMassvsPt(TFile* F, const char* SDIR)
{
	//Link hist object
	TObject* hist = F->Get(Form("%s/histogram", SDIR));

    TH2D* XiMassvsPt = (TH2D*) hist->FindObject("hXimassvsPt")->Clone("massvspt");

    TH1D *hXi0_1  = XiMassvsPt->ProjectionX("hXiMass0_1",  0,  1);
    TH1D *hXi1_2  = XiMassvsPt->ProjectionX("hXiMass1_2",  1,  2);
    TH1D *hXi2_3  = XiMassvsPt->ProjectionX("hXiMass2_3",  2,  3);
    TH1D *hXi3_4  = XiMassvsPt->ProjectionX("hXiMass3_4",  3,  4);
    TH1D *hXi4_5  = XiMassvsPt->ProjectionX("hXiMass4_5",  4,  5);
    TH1D *hXi5_6  = XiMassvsPt->ProjectionX("hXiMass5_6",  5,  6);
    TH1D *hXi6_7  = XiMassvsPt->ProjectionX("hXiMass6_7",  6,  7);
    TH1D *hXi7_8  = XiMassvsPt->ProjectionX("hXiMass7_8",  7,  8);
    TH1D *hXi8_9  = XiMassvsPt->ProjectionX("hXiMass8_9",  8,  9);
    TH1D *hXi9_10 = XiMassvsPt->ProjectionX("hXiMass9_10", 9, 10);

    return;
}

//------------------------------------------------
void XiCutDistribution(TFile *F, const char* SDIR)
{
	//Link hist object
	TObject* hist = F->Get(Form("%s/histogram", SDIR));

    TH1D* hC  = (TH1D*) hist->FindObject("C_flag")->Clone("C_flag");
    TH1D* hB  = (TH1D*) hist->FindObject("B_flag")->Clone("B_flag");
    TH1D* hBc = (TH1D*) hist->FindObject("Bcut_flag")->Clone("Bcut_flag");
    TH1D* hCc = (TH1D*) hist->FindObject("Ccut_flag")->Clone("Ccut_flag");
    TH1D* hec = (TH1D*) hist->FindObject("e_c_flag")->Clone("hec");
    TH1D* heb = (TH1D*) hist->FindObject("e_b_flag")->Clone("heb");

	//kimc: removed suffix " Min "
    TH1D* b1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertex_b")->Clone("hDCAV0PrToPrimVertex_b");
    TH1D* b2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertex_b")->Clone("hDCAV0PiToPrimVertex_b");
    TH1D* b3 = (TH1D*) hist->FindObject("hDCABachToPrimVertex_b")->Clone("hDCABachToPrimVertex_b");
    TH1D* b4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertex_b")->Clone("hDCAV0ToPrimVertex_b");
    TH1D* b5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXi_b")->Clone("hV0CosineOfPoiningAngleXi_b");
    TH1D* b6 = (TH1D*) hist->FindObject("hCascDecayLength_b")->Clone("hCascDecayLength_b");
    TH1D* b7 = (TH1D*) hist->FindObject("hDecayLengthV0_b")->Clone("hDecayLengthV0_b");

    TH1D* c1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertex_c")->Clone("hDCAV0PrToPrimVertex_c");
    TH1D* c2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertex_c")->Clone("hDCAV0PiToPrimVertex_c");
    TH1D* c3 = (TH1D*) hist->FindObject("hDCABachToPrimVertex_c")->Clone("hDCABachToPrimVertex_c");
    TH1D* c4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertex_c")->Clone("hDCAV0ToPrimVertex_c");
    TH1D* c5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXi_c")->Clone("hV0CosineOfPoiningAngleXi_c");
    TH1D* c6 = (TH1D*) hist->FindObject("hCascDecayLength_c")->Clone("hCascDecayLength_c");
    TH1D* c7 = (TH1D*) hist->FindObject("hDecayLengthV0_c")->Clone("hDecayLengthV0_c");

    TH1D* d1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertex")->Clone("hDCAV0PrToPrimVertex");
    TH1D* d2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertex")->Clone("hDCAV0PiToPrimVertex");
    TH1D* d3 = (TH1D*) hist->FindObject("hDCABachToPrimVertex")->Clone("hDCABachToPrimVertex");
    TH1D* d4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertex")->Clone("hDCAV0ToPrimVertex");
    TH1D* d5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXi")->Clone("hV0CosineOfPoiningAngleXi");
    TH1D* d6 = (TH1D*) hist->FindObject("hCascDecayLength")->Clone("hCascDecayLength");
    TH1D* d7 = (TH1D*) hist->FindObject("hDecayLengthV0")->Clone("hDecayLengthV0");

    return;
}//XiCutDistribution

//----------------------------------------------------------------------------------------------------------
void eXiPairTree(TFile* F, bool IsMC, const char* SDIR, const char* TRIG, double* MultPerc, double* WFitPar)
{
	//Trigger bits
	if ( strcmp(TRIG, "MB") && strcmp(TRIG, "HMV0") && strcmp(TRIG, "HMSPD") && strcmp(TRIG, "HMOR") )
	{
		cout <<"No valid TRIG provided! Stop.\n";
		return;
	}

	UInt_t TrigMB    = 2;
	UInt_t TrigHMV0  = 65536;
	UInt_t TrigHMSPD = 8;

	//MC weight fit parameters
	TF1* fWeightFit = new TF1("fWeightFit", "1", 0, 20);
	//fWeightFit->SetParameter(0, WFitPar[0]);
	//fWeightFit->SetParameter(1, WFitPar[1]);

	Double_t binning[] = {0.,1.,2.,3.,4.,5.,6.,8.,12.,16.,20}; //11
  Double_t binning2[8] = {1.,2.,3.,4.,5.,6.,8.,12.};

	//Link tree
	//*******************************************

	//!
	//Histograms container
	TObject* hist = F->Get(Form("%s/histogram", SDIR));

	//!
    //EventTree
    TTree* EventTree = (TTree*)F->Get(Form("%s/EventTree", SDIR));
	if (!EventTree) { cout <<"Cannot link the EventTree!" <<endl; return; }
    Float_t fRunNumber, fCentrality, fCentralSPD, fNSPDTracklets, fNeXiPair, fZvtx;
    UInt_t fTrigBit;
    EventTree->SetBranchAddress("fRunNumber",     &fRunNumber);
    EventTree->SetBranchAddress("fCentrality",    &fCentrality);
    EventTree->SetBranchAddress("fCentralSPD",    &fCentralSPD);
    EventTree->SetBranchAddress("fNSPDTracklets", &fNSPDTracklets);
    EventTree->SetBranchAddress("fNeXiPair",      &fNeXiPair);
    EventTree->SetBranchAddress("fVtxZ",          &fZvtx);
    EventTree->SetBranchAddress("fTrigBit",       &fTrigBit);

    TTree* Pair = (TTree*)F->Get(Form("%s/eXiTree", SDIR)); //!
	if (!Pair) { cout <<"Cannot link the eXiTree!" <<endl; return; }
	Float_t pTe, echarge, pTv, vcharge, Massv;
	Float_t cosoa, In_Mass, Pt, nSigmaTOF, nSigmaTPC;
	Float_t TPCCluster, ITSCluster, TPCPIDCluster, CascDecayLength, V0DecayLength;
	Float_t DCABachToPrimVertex, V0CosineOfPoiningAngleXi, DCAV0ToPrimVertex, DCAPosToPrimVertex, DCANegToPrimVertex;
	Float_t e_minmass, e_minmass_ss, phi, erap, Xirap;
	Float_t pionTPCCluster, protonTPCCluster, b_pionTPCCluster, e_crossedratio, e_findable;
	Float_t pion_crossedratio, pion_findable, proton_crossedratio, proton_findable, bpion_crossedratio;
	Float_t bpion_findable, XiCosineOfPoiningAngle, pTpion, pTproton, pTbach;
	Float_t MassLambda, MassAntiLambda;
    Pair->SetBranchAddress("pTe", &pTe);
    Pair->SetBranchAddress("echarge", &echarge);
    Pair->SetBranchAddress("TOFnSigma", &nSigmaTOF);
    Pair->SetBranchAddress("TPCnSigma", &nSigmaTPC);
    //Pair->SetBranchAddress("TPC", &TPCCluster);
    Pair->SetBranchAddress("TPCPID", &TPCPIDCluster);
    Pair->SetBranchAddress("ITS", &ITSCluster);
    Pair->SetBranchAddress("e_crossedrows", &e_crossedratio);
    Pair->SetBranchAddress("e_findable", &e_findable);
    Pair->SetBranchAddress("phi", &phi);
    Pair->SetBranchAddress("erap", &erap);
    Pair->SetBranchAddress("e_minmass", &e_minmass);
    Pair->SetBranchAddress("e_minmass_ss", &e_minmass_ss);
    Pair->SetBranchAddress("pTv", &pTv);
    Pair->SetBranchAddress("vcharge", &vcharge);
    Pair->SetBranchAddress("Massv", &Massv);
    Pair->SetBranchAddress("MassLambda", &MassLambda);
    Pair->SetBranchAddress("MassAntiLambda", &MassAntiLambda);  //mod
    Pair->SetBranchAddress("Xirap", &Xirap);  //mod
    Pair->SetBranchAddress("V0DecayLength", &V0DecayLength);
    Pair->SetBranchAddress("CascDecayLength", &CascDecayLength);
    Pair->SetBranchAddress("DCABachToPrimVertex", &DCABachToPrimVertex);
    Pair->SetBranchAddress("DCAV0NegToPrimVertex", &DCANegToPrimVertex);
    Pair->SetBranchAddress("DCAV0PosToPrimVertex", &DCAPosToPrimVertex);
    Pair->SetBranchAddress("V0CosineOfPoiningAngleXi", &V0CosineOfPoiningAngleXi);
    Pair->SetBranchAddress("XiCosineOfPoiningAngle", &XiCosineOfPoiningAngle);  //modify
    Pair->SetBranchAddress("DCAV0ToPrimVertex", &DCAV0ToPrimVertex);   ///new
    Pair->SetBranchAddress("cosoa", &cosoa);  //mod
    Pair->SetBranchAddress("pion_crossedrows", &pion_crossedratio);
    Pair->SetBranchAddress("pion_findable", &pion_findable);
    Pair->SetBranchAddress("proton_crossedrows", &proton_crossedratio);
    Pair->SetBranchAddress("proton_findable", &proton_findable);
    Pair->SetBranchAddress("bpion_crossedratio", &bpion_crossedratio);
    Pair->SetBranchAddress("bpion_findable", &bpion_findable);
    Pair->SetBranchAddress("pTpion", &pTpion);
    Pair->SetBranchAddress("pTproton", &pTproton);
    Pair->SetBranchAddress("pTbach", &pTbach);
    Pair->SetBranchAddress("In_Mass", &In_Mass);
    Pair->SetBranchAddress("eXiPt", &Pt);

	//!
  TTree* MCTree = (TTree*)F->Get(Form("%s/MCTree", SDIR));
	if (!MCTree) { cout <<"Cannot link the MCTree!" <<endl; return; }
	Float_t mcpTe, mcecharge, mcpTv, mcvcharge, mcpteXi, mcptXic0, mcc_flag, mcb_flag, mcXib, XibeXi, mcXibMass;
	MCTree->SetBranchAddress("mcpTe", &mcpTe);
	MCTree->SetBranchAddress("mcecharge", &mcecharge);
	MCTree->SetBranchAddress("mcpTv", &mcpTv);
	MCTree->SetBranchAddress("mcvcharge", &mcvcharge);
	MCTree->SetBranchAddress("mcpTeXi", &mcpteXi);
	MCTree->SetBranchAddress("mcpTXic0", &mcptXic0);
	MCTree->SetBranchAddress("c_flag", &mcc_flag);
	MCTree->SetBranchAddress("b_flag", &mcb_flag);
	MCTree->SetBranchAddress("mcpTXib", &mcXib);
	MCTree->SetBranchAddress("mceXipTb", &XibeXi);
	MCTree->SetBranchAddress("mcXibMass", &mcXibMass);

  TTree* MCXicTree = (TTree*)F->Get(Form("%s/MCXicTree", SDIR)); //!
  Float_t mcGenXic0pT, mcGenepT, mcGenXipT, mcGenXic0rap, mcGenerap, mcGenXirap, mcGencflag, mcGenbflag;
  if(IsMC){
    if (!MCXicTree) { cout <<"Cannot link the MCXicTree!" <<endl; return; }

    MCXicTree->SetBranchAddress("Xic0_pT", &mcGenXic0pT);
  	MCXicTree->SetBranchAddress("e_pT", &mcGenepT);
  	MCXicTree->SetBranchAddress("Xi_pT", &mcGenXipT);
  	MCXicTree->SetBranchAddress("Xic0_rap", &mcGenXic0rap);
  	MCXicTree->SetBranchAddress("e_rap", &mcGenerap);
  	MCXicTree->SetBranchAddress("Xi_rap", &mcGenXirap);
  	MCXicTree->SetBranchAddress("c_flag", &mcGencflag);
  	MCXicTree->SetBranchAddress("b_flag", &mcGenbflag);
  }

	//Histograms, raw yield
	//*******************************************

    TH1F *hPtRS_TePID = MakeTH1("hPtRS_TePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TePID = MakeTH1("hPtWS_TePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SePID = MakeTH1("hPtRS_SePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SePID = MakeTH1("hPtWS_SePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_LePID = MakeTH1("hPtRS_LePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LePID = MakeTH1("hPtWS_LePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTePID = MakeTH1("hPtRS_VTePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTePID = MakeTH1("hPtWS_VTePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLePID = MakeTH1("hPtRS_VLePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLePID = MakeTH1("hPtWS_VLePID",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_TXiPID = MakeTH1("hPtRS_TXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TXiPID = MakeTH1("hPtWS_TXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiPID = MakeTH1("hPtRS_SXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPID = MakeTH1("hPtWS_SXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_LXiPID = MakeTH1("hPtRS_LXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LXiPID = MakeTH1("hPtWS_LXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTXiPID = MakeTH1("hPtRS_VTXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTXiPID = MakeTH1("hPtWS_VTXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLXiPID = MakeTH1("hPtRS_VLXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLXiPID = MakeTH1("hPtWS_VLXiPID",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_TXiPIDRec = MakeTH1("hPtRS_TXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TXiPIDRec = MakeTH1("hPtWS_TXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiPIDRec = MakeTH1("hPtRS_SXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPIDRec = MakeTH1("hPtWS_SXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_LXiPIDRec = MakeTH1("hPtRS_LXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LXiPIDRec = MakeTH1("hPtWS_LXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTXiPIDRec = MakeTH1("hPtRS_VTXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTXiPIDRec = MakeTH1("hPtWS_VTXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLXiPIDRec = MakeTH1("hPtRS_VLXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLXiPIDRec = MakeTH1("hPtWS_VLXiPIDRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_TXiPIDLRec = MakeTH1("hPtRS_TXiPIDLRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TXiPIDLRec = MakeTH1("hPtWS_TXiPIDLRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiPIDSRec = MakeTH1("hPtRS_SXiPIDSRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPIDSRec = MakeTH1("hPtWS_SXiPIDSRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_LXiPIDTRec = MakeTH1("hPtRS_LXiPIDTRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LXiPIDTRec = MakeTH1("hPtWS_LXiPIDTRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTXiPIDVLRec = MakeTH1("hPtRS_VTXiPIDVLRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTXiPIDVLRec = MakeTH1("hPtWS_VTXiPIDVLRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLXiPIDVTRec = MakeTH1("hPtRS_VLXiPIDVTRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLXiPIDVTRec = MakeTH1("hPtWS_VLXiPIDVTRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_SXiPID_V0 = MakeTH1("hPtRS_SXiPID_V0",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPID_V0 = MakeTH1("hPtWS_SXiPID_V0",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiPID_Xi = MakeTH1("hPtRS_SXiPID_Xi",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPID_Xi = MakeTH1("hPtWS_SXiPID_Xi",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiPID_DCAb = MakeTH1("hPtRS_SXiPID_DCAb",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiPID_DCAb = MakeTH1("hPtWS_SXiPID_DCAb",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_LeRec = MakeTH1("hPtRS_LeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LeRec = MakeTH1("hPtWS_LeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SeRec = MakeTH1("hPtRS_SeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SeRec = MakeTH1("hPtWS_SeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_TeRec = MakeTH1("hPtRS_TeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TeRec = MakeTH1("hPtWS_TeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTeRec = MakeTH1("hPtRS_VTeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTeRec = MakeTH1("hPtWS_VTeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLeRec = MakeTH1("hPtRS_VLeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLeRec = MakeTH1("hPtWS_VLeRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_LXiRec = MakeTH1("hPtRS_LXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LXiRec = MakeTH1("hPtWS_LXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SXiRec = MakeTH1("hPtRS_SXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SXiRec = MakeTH1("hPtWS_SXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_TXiRec = MakeTH1("hPtRS_TXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TXiRec = MakeTH1("hPtWS_TXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VTXiRec = MakeTH1("hPtRS_VTXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VTXiRec = MakeTH1("hPtWS_VTXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_VLXiRec = MakeTH1("hPtRS_VLXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_VLXiRec = MakeTH1("hPtWS_VLXiRec",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_LeIM = MakeTH1("hPtRS_LeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LeIM = MakeTH1("hPtWS_LeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SeIM = MakeTH1("hPtRS_SeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SeIM = MakeTH1("hPtWS_SeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_TeIM = MakeTH1("hPtRS_TeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TeIM = MakeTH1("hPtWS_TeIM",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hPtRS_LOA = MakeTH1("hPtRS_LOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_LOA = MakeTH1("hPtWS_LOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_TOA = MakeTH1("hPtRS_TOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_TOA = MakeTH1("hPtWS_TOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hPtRS_SOA = MakeTH1("hPtRS_SOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hPtWS_SOA = MakeTH1("hPtWS_SOA",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2D *hMassPtRS = new TH2D("hMassPtRS","",10,1.3,3.3,60,0,12); hMassPtRS->Sumw2();
    TH2D *hMassPtWS = new TH2D("hMassPtWS","",10,1.3,3.3,60,0,12); hMassPtWS->Sumw2();

	//Histograms, prefilter
	//*******************************************

    TH1F *hpre_eRec_loose_de = MakeTH1("hpre_eRec_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_eRec_loose_nu = MakeTH1("hpre_eRec_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_eRec_stand_de = MakeTH1("hpre_eRec_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_eRec_stand_nu = MakeTH1("hpre_eRec_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_eRec_tight_de = MakeTH1("hpre_eRec_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_eRec_tight_nu = MakeTH1("hpre_eRec_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_eRec_vloose_de = MakeTH1("hpre_eRec_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_eRec_vloose_nu = MakeTH1("hpre_eRec_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_eRec_vtight_de = MakeTH1("hpre_eRec_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_eRec_vtight_nu = MakeTH1("hpre_eRec_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_XiRec_loose_de = MakeTH1("hpre_XiRec_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiRec_loose_nu = MakeTH1("hpre_XiRec_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiRec_stand_de = MakeTH1("hpre_XiRec_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiRec_stand_nu = MakeTH1("hpre_XiRec_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiRec_tight_de = MakeTH1("hpre_XiRec_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiRec_tight_nu = MakeTH1("hpre_XiRec_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiRec_vloose_de = MakeTH1("hpre_XiRec_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiRec_vloose_nu = MakeTH1("hpre_XiRec_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiRec_vtight_de = MakeTH1("hpre_XiRec_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiRec_vtight_nu = MakeTH1("hpre_XiRec_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_ePID_loose_de = MakeTH1("hpre_ePID_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_ePID_loose_nu = MakeTH1("hpre_ePID_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_ePID_stand_de = MakeTH1("hpre_ePID_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_ePID_stand_nu = MakeTH1("hpre_ePID_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_ePID_tight_de = MakeTH1("hpre_ePID_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_ePID_tight_nu = MakeTH1("hpre_ePID_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_ePID_vloose_de = MakeTH1("hpre_ePID_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_ePID_vloose_nu = MakeTH1("hpre_ePID_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_ePID_vtight_de = MakeTH1("hpre_ePID_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_ePID_vtight_nu = MakeTH1("hpre_ePID_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_XiPID_loose_de = MakeTH1("hpre_XiPID_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPID_loose_nu = MakeTH1("hpre_XiPID_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPID_stand_de = MakeTH1("hpre_XiPID_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPID_stand_nu = MakeTH1("hpre_XiPID_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPID_tight_de = MakeTH1("hpre_XiPID_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPID_tight_nu = MakeTH1("hpre_XiPID_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPID_vloose_de = MakeTH1("hpre_XiPID_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPID_vloose_nu = MakeTH1("hpre_XiPID_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPID_vtight_de = MakeTH1("hpre_XiPID_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPID_vtight_nu = MakeTH1("hpre_XiPID_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_XiPIDRec_loose_de = MakeTH1("hpre_XiPIDRec_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDRec_loose_nu = MakeTH1("hpre_XiPIDRec_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDRec_stand_de = MakeTH1("hpre_XiPIDRec_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDRec_stand_nu = MakeTH1("hpre_XiPIDRec_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDRec_tight_de = MakeTH1("hpre_XiPIDRec_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDRec_tight_nu = MakeTH1("hpre_XiPIDRec_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDRec_vloose_de = MakeTH1("hpre_XiPIDRec_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDRec_vloose_nu = MakeTH1("hpre_XiPIDRec_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDRec_vtight_de = MakeTH1("hpre_XiPIDRec_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDRec_vtight_nu = MakeTH1("hpre_XiPIDRec_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_XiPIDinvRec_loose_de = MakeTH1("hpre_XiPIDinvRec_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDinvRec_loose_nu = MakeTH1("hpre_XiPIDinvRec_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDinvRec_stand_de = MakeTH1("hpre_XiPIDinvRec_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDinvRec_stand_nu = MakeTH1("hpre_XiPIDinvRec_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDinvRec_tight_de = MakeTH1("hpre_XiPIDinvRec_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDinvRec_tight_nu = MakeTH1("hpre_XiPIDinvRec_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDinvRec_vloose_de = MakeTH1("hpre_XiPIDinvRec_vloose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDinvRec_vloose_nu = MakeTH1("hpre_XiPIDinvRec_vloose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_XiPIDinvRec_vtight_de = MakeTH1("hpre_XiPIDinvRec_vtight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_XiPIDinvRec_vtight_nu = MakeTH1("hpre_XiPIDinvRec_vtight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_OA_loose_de = MakeTH1("hpre_OA_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_OA_loose_nu = MakeTH1("hpre_OA_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_OA_stand_de = MakeTH1("hpre_OA_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_OA_stand_nu = MakeTH1("hpre_OA_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_OA_tight_de = MakeTH1("hpre_OA_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_OA_tight_nu = MakeTH1("hpre_OA_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_IM_loose_de = MakeTH1("hpre_IM_loose_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_IM_loose_nu = MakeTH1("hpre_IM_loose_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_IM_stand_de = MakeTH1("hpre_IM_stand_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_IM_stand_nu = MakeTH1("hpre_IM_stand_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_IM_tight_de = MakeTH1("hpre_IM_tight_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_IM_tight_nu = MakeTH1("hpre_IM_tight_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_Bayes_stand2_de = MakeTH1("hpre_Bayes_stand2_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand2_nu = MakeTH1("hpre_Bayes_stand2_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Bayes_stand3_de = MakeTH1("hpre_Bayes_stand3_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand3_nu = MakeTH1("hpre_Bayes_stand3_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Bayes_stand4_de = MakeTH1("hpre_Bayes_stand4_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand4_nu = MakeTH1("hpre_Bayes_stand4_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Bayes_stand5_de = MakeTH1("hpre_Bayes_stand5_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand5_nu = MakeTH1("hpre_Bayes_stand5_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Bayes_stand6_de = MakeTH1("hpre_Bayes_stand6_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand6_nu = MakeTH1("hpre_Bayes_stand6_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Bayes_stand7_de = MakeTH1("hpre_Bayes_stand7_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Bayes_stand7_nu = MakeTH1("hpre_Bayes_stand7_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hpre_Svd_stand3_de = MakeTH1("hpre_Svd_stand3_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Svd_stand3_nu = MakeTH1("hpre_Svd_stand3_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Svd_stand4_de = MakeTH1("hpre_Svd_stand4_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Svd_stand4_nu = MakeTH1("hpre_Svd_stand4_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F *hpre_Svd_stand5_de = MakeTH1("hpre_Svd_stand5_de",(sizeof(binning) / sizeof(binning[0]))-1,binning);
	TH1F *hpre_Svd_stand5_nu = MakeTH1("hpre_Svd_stand5_nu",(sizeof(binning) / sizeof(binning[0]))-1,binning);

	//Histograms, MC, usingunfolding and efficiency calculation
	//*********************************************************

    TH2F* hRPM_eRec_loose_un        = MakeTH2("hRPM_eRec_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning); //response matrix using unfolding
    TH2F* hRPM_eRec_loose           = MakeTH2("hRPM_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);    //response matrix using refolding
    TH1F* hMCRecoLevXic0_eRec_loose = MakeTH1("hMCRecoLevXic0_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning); //Xic0 in generation level
    TH1F* hMCRecoLevPair_eRec_loose = MakeTH1("hMCRecoLevPair_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning); //eXi pair decay from Xic0 in generation level
    TH2F* hRPM_eRec_stand_un        = MakeTH2("hRPM_eRec_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_stand           = MakeTH2("hRPM_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_eRec_stand = MakeTH1("hMCRecoLevXic0_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_eRec_stand = MakeTH1("hMCRecoLevPair_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_tight_un        = MakeTH2("hRPM_eRec_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_tight           = MakeTH2("hRPM_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_eRec_tight = MakeTH1("hMCRecoLevXic0_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_eRec_tight = MakeTH1("hMCRecoLevPair_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vloose_un        = MakeTH2("hRPM_eRec_vloose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vloose           = MakeTH2("hRPM_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_eRec_vloose = MakeTH1("hMCRecoLevXic0_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_eRec_vloose = MakeTH1("hMCRecoLevPair_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vtight_un        = MakeTH2("hRPM_eRec_vtight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vtight           = MakeTH2("hRPM_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_eRec_vtight = MakeTH1("hMCRecoLevXic0_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_eRec_vtight = MakeTH1("hMCRecoLevPair_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_ePID_loose_un        = MakeTH2("hRPM_ePID_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_loose           = MakeTH2("hRPM_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_ePID_loose = MakeTH1("hMCRecoLevXic0_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_ePID_loose = MakeTH1("hMCRecoLevPair_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_stand_un        = MakeTH2("hRPM_ePID_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_stand           = MakeTH2("hRPM_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_ePID_stand = MakeTH1("hMCRecoLevXic0_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_ePID_stand = MakeTH1("hMCRecoLevPair_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_tight_un        = MakeTH2("hRPM_ePID_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_tight           = MakeTH2("hRPM_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_ePID_tight = MakeTH1("hMCRecoLevXic0_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_ePID_tight = MakeTH1("hMCRecoLevPair_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vloose_un        = MakeTH2("hRPM_ePID_vloose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vloose           = MakeTH2("hRPM_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_ePID_vloose = MakeTH1("hMCRecoLevXic0_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_ePID_vloose = MakeTH1("hMCRecoLevPair_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vtight_un        = MakeTH2("hRPM_ePID_vtight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vtight           = MakeTH2("hRPM_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_ePID_vtight = MakeTH1("hMCRecoLevXic0_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_ePID_vtight = MakeTH1("hMCRecoLevPair_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_XiRec_loose_un        = MakeTH2("hRPM_XiRec_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_loose           = MakeTH2("hRPM_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiRec_loose = MakeTH1("hMCRecoLevXic0_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiRec_loose = MakeTH1("hMCRecoLevPair_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_stand_un        = MakeTH2("hRPM_XiRec_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_stand           = MakeTH2("hRPM_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiRec_stand = MakeTH1("hMCRecoLevXic0_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiRec_stand = MakeTH1("hMCRecoLevPair_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_tight_un        = MakeTH2("hRPM_XiRec_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_tight           = MakeTH2("hRPM_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiRec_tight = MakeTH1("hMCRecoLevXic0_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiRec_tight = MakeTH1("hMCRecoLevPair_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vloose_un        = MakeTH2("hRPM_XiRec_vloose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vloose           = MakeTH2("hRPM_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiRec_vloose = MakeTH1("hMCRecoLevXic0_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiRec_vloose = MakeTH1("hMCRecoLevPair_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vtight_un        = MakeTH2("hRPM_XiRec_vtight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vtight           = MakeTH2("hRPM_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiRec_vtight = MakeTH1("hMCRecoLevXic0_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiRec_vtight = MakeTH1("hMCRecoLevPair_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

	TH2F* hRPM_XiPID_loose_un        = MakeTH2("hRPM_XiPID_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_loose           = MakeTH2("hRPM_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiPID_loose = MakeTH1("hMCRecoLevXic0_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiPID_loose = MakeTH1("hMCRecoLevPair_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_stand_un        = MakeTH2("hRPM_XiPID_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_stand           = MakeTH2("hRPM_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiPID_stand = MakeTH1("hMCRecoLevXic0_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiPID_stand = MakeTH1("hMCRecoLevPair_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_tight_un        = MakeTH2("hRPM_XiPID_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_tight           = MakeTH2("hRPM_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiPID_tight = MakeTH1("hMCRecoLevXic0_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiPID_tight = MakeTH1("hMCRecoLevPair_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vloose_un        = MakeTH2("hRPM_XiPID_vloose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vloose           = MakeTH2("hRPM_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiPID_vloose = MakeTH1("hMCRecoLevXic0_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiPID_vloose = MakeTH1("hMCRecoLevPair_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vtight_un        = MakeTH2("hRPM_XiPID_vtight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vtight           = MakeTH2("hRPM_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_XiPID_vtight = MakeTH1("hMCRecoLevXic0_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_XiPID_vtight = MakeTH1("hMCRecoLevPair_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_op_loose_un        = MakeTH2("hRPM_OA_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_loose           = MakeTH2("hRPM_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_op_loose = MakeTH1("hMCRecoLevXic0_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_op_loose = MakeTH1("hMCRecoLevPair_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_stand_un        = MakeTH2("hRPM_OA_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_stand           = MakeTH2("hRPM_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_op_stand = MakeTH1("hMCRecoLevXic0_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_op_stand = MakeTH1("hMCRecoLevPair_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_tight_un        = MakeTH2("hRPM_OA_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_tight           = MakeTH2("hRPM_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_op_tight = MakeTH1("hMCRecoLevXic0_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_op_tight = MakeTH1("hMCRecoLevPair_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_im_loose_un        = MakeTH2("hRPM_IM_loose_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_loose           = MakeTH2("hRPM_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_im_loose = MakeTH1("hMCRecoLevXic0_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_im_loose = MakeTH1("hMCRecoLevPair_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_stand_un        = MakeTH2("hRPM_IM_stand_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_stand           = MakeTH2("hRPM_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_im_stand = MakeTH1("hMCRecoLevXic0_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_im_stand = MakeTH1("hMCRecoLevPair_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_tight_un        = MakeTH2("hRPM_IM_tight_un",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_tight           = MakeTH2("hRPM_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXic0_im_tight = MakeTH1("hMCRecoLevXic0_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPair_im_tight = MakeTH1("hMCRecoLevPair_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_eRec_loose_Xib          = MakeTH2("hRPM_eRec_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_eRec_loose     = MakeTH1("hMCRecoLevXib_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_eRec_loose = MakeTH1("hMCRecoLevPairXib_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_stand_Xib          = MakeTH2("hRPM_eRec_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_eRec_stand     = MakeTH1("hMCRecoLevXib_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_eRec_stand = MakeTH1("hMCRecoLevPairXib_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_tight_Xib          = MakeTH2("hRPM_eRec_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_eRec_tight     = MakeTH1("hMCRecoLevXib_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_eRec_tight = MakeTH1("hMCRecoLevPairXib_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vloose_Xib          = MakeTH2("hRPM_eRec_vloose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_eRec_vloose     = MakeTH1("hMCRecoLevXib_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_eRec_vloose = MakeTH1("hMCRecoLevPairXib_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_eRec_vtight_Xib          = MakeTH2("hRPM_eRec_vtight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_eRec_vtight     = MakeTH1("hMCRecoLevXib_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_eRec_vtight = MakeTH1("hMCRecoLevPairXib_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_ePID_loose_Xib          = MakeTH2("hRPM_ePID_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_ePID_loose     = MakeTH1("hMCRecoLevXib_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_ePID_loose = MakeTH1("hMCRecoLevPairXib_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_stand_Xib          = MakeTH2("hRPM_ePID_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_ePID_stand     = MakeTH1("hMCRecoLevXib_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_ePID_stand = MakeTH1("hMCRecoLevPairXib_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_tight_Xib          = MakeTH2("hRPM_ePID_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_ePID_tight     = MakeTH1("hMCRecoLevXib_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_ePID_tight = MakeTH1("hMCRecoLevPairXib_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vloose_Xib          = MakeTH2("hRPM_ePID_vloose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_ePID_vloose     = MakeTH1("hMCRecoLevXib_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_ePID_vloose = MakeTH1("hMCRecoLevPairXib_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_ePID_vtight_Xib          = MakeTH2("hRPM_ePID_vtight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_ePID_vtight     = MakeTH1("hMCRecoLevXib_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_ePID_vtight = MakeTH1("hMCRecoLevPairXib_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_XiPID_loose_Xib          = MakeTH2("hRPM_XiPID_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiPID_loose     = MakeTH1("hMCRecoLevXib_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiPID_loose = MakeTH1("hMCRecoLevPairXib_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_stand_Xib          = MakeTH2("hRPM_XiPID_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiPID_stand     = MakeTH1("hMCRecoLevXib_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiPID_stand = MakeTH1("hMCRecoLevPairXib_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_tight_Xib          = MakeTH2("hRPM_XiPID_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiPID_tight     = MakeTH1("hMCRecoLevXib_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiPID_tight = MakeTH1("hMCRecoLevPairXib_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vloose_Xib          = MakeTH2("hRPM_XiPID_vloose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiPID_vloose     = MakeTH1("hMCRecoLevXib_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiPID_vloose = MakeTH1("hMCRecoLevPairXib_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiPID_vtight_Xib          = MakeTH2("hRPM_XiPID_vtight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiPID_vtight     = MakeTH1("hMCRecoLevXib_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiPID_vtight = MakeTH1("hMCRecoLevPairXib_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_XiRec_loose_Xib          = MakeTH2("hRPM_XiRec_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiRec_loose     = MakeTH1("hMCRecoLevXib_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiRec_loose = MakeTH1("hMCRecoLevPairXib_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_stand_Xib          = MakeTH2("hRPM_XiRec_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiRec_stand     = MakeTH1("hMCRecoLevXib_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiRec_stand = MakeTH1("hMCRecoLevPairXib_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_tight_Xib          = MakeTH2("hRPM_XiRec_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiRec_tight     = MakeTH1("hMCRecoLevXib_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiRec_tight = MakeTH1("hMCRecoLevPairXib_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vloose_Xib          = MakeTH2("hRPM_XiRec_vloose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiRec_vloose     = MakeTH1("hMCRecoLevXib_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiRec_vloose = MakeTH1("hMCRecoLevPairXib_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_XiRec_vtight_Xib          = MakeTH2("hRPM_XiRec_vtight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_XiRec_vtight     = MakeTH1("hMCRecoLevXib_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_XiRec_vtight = MakeTH1("hMCRecoLevPairXib_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_im_loose_Xib          = MakeTH2("hRPM_IM_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_im_loose     = MakeTH1("hMCRecoLevXib_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_im_loose = MakeTH1("hMCRecoLevPairXib_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_stand_Xib          = MakeTH2("hRPM_IM_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_im_stand     = MakeTH1("hMCRecoLevXib_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_im_stand = MakeTH1("hMCRecoLevPairXib_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_im_tight_Xib          = MakeTH2("hRPM_IM_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_im_tight     = MakeTH1("hMCRecoLevXib_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_im_tight = MakeTH1("hMCRecoLevPairXib_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH2F* hRPM_op_loose_Xib          = MakeTH2("hRPM_OA_loose_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_op_loose     = MakeTH1("hMCRecoLevXib_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_op_loose = MakeTH1("hMCRecoLevPairXib_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_stand_Xib          = MakeTH2("hRPM_OA_stand_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_op_stand     = MakeTH1("hMCRecoLevXib_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_op_stand = MakeTH1("hMCRecoLevPairXib_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* hRPM_op_tight_Xib          = MakeTH2("hRPM_OA_tight_Xib",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevXib_op_tight     = MakeTH1("hMCRecoLevXib_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPairXib_op_tight = MakeTH1("hMCRecoLevPairXib_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_eRec_loose = MakeTH1("hMCRecoLevPromptXic0_eRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_eRec_stand = MakeTH1("hMCRecoLevPromptXic0_eRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_eRec_tight = MakeTH1("hMCRecoLevPromptXic0_eRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_eRec_vloose = MakeTH1("hMCRecoLevPromptXic0_eRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_eRec_vtight = MakeTH1("hMCRecoLevPromptXic0_eRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_XiRec_loose = MakeTH1("hMCRecoLevPromptXic0_XiRec_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiRec_stand = MakeTH1("hMCRecoLevPromptXic0_XiRec_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiRec_tight = MakeTH1("hMCRecoLevPromptXic0_XiRec_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiRec_vloose = MakeTH1("hMCRecoLevPromptXic0_XiRec_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiRec_vtight = MakeTH1("hMCRecoLevPromptXic0_XiRec_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_ePID_loose = MakeTH1("hMCRecoLevPromptXic0_ePID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_ePID_stand = MakeTH1("hMCRecoLevPromptXic0_ePID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_ePID_tight = MakeTH1("hMCRecoLevPromptXic0_ePID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_ePID_vloose = MakeTH1("hMCRecoLevPromptXic0_ePID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_ePID_vtight = MakeTH1("hMCRecoLevPromptXic0_ePID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_XiPID_loose = MakeTH1("hMCRecoLevPromptXic0_XiPID_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiPID_stand = MakeTH1("hMCRecoLevPromptXic0_XiPID_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiPID_tight = MakeTH1("hMCRecoLevPromptXic0_XiPID_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiPID_vloose = MakeTH1("hMCRecoLevPromptXic0_XiPID_vloose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_XiPID_vtight = MakeTH1("hMCRecoLevPromptXic0_XiPID_vtight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_IM_loose = MakeTH1("hMCRecoLevPromptXic0_IM_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_IM_stand = MakeTH1("hMCRecoLevPromptXic0_IM_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_IM_tight = MakeTH1("hMCRecoLevPromptXic0_IM_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCRecoLevPromptXic0_OA_loose = MakeTH1("hMCRecoLevPromptXic0_OA_loose",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_OA_stand = MakeTH1("hMCRecoLevPromptXic0_OA_stand",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH1F* hMCRecoLevPromptXic0_OA_tight = MakeTH1("hMCRecoLevPromptXic0_OA_tight",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F* hMCGenInclusiveXic0_woW = MakeTH1("hMCGenInclusiveXic0_woW",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenPromptXic0_woW = MakeTH1("hMCGenPromptXic0_woW",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenFeeddowmXic0_woW = MakeTH1("hMCGenFeeddowmXic0_woW",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);

    TH1F* hMCGenInclusiveXic0_W_rap08 = MakeTH1("hMCGenInclusiveXic0_W_rap08",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenPromptXic0_W_rap08 = MakeTH1("hMCGenPromptXic0_W_rap08",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenFeeddowmXic0_W_rap08 = MakeTH1("hMCGenFeeddowmXic0_W_rap08",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);

    TH1F* hMCGenInclusiveXic0_W = MakeTH1("hMCGenInclusiveXic0_W",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenPromptXic0_W = MakeTH1("hMCGenPromptXic0_W",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);
    TH1F* hMCGenFeeddowmXic0_W = MakeTH1("hMCGenFeeddowmXic0_W",(sizeof(binning2) / sizeof(binning2[0]))-1,binning2);

    TH2F* hevseXi = MakeTH2("hevseXi",(sizeof(binning) / sizeof(binning[0]))-1,binning);
    TH2F* heXivsXic0 = MakeTH2("heXivsXic0",(sizeof(binning) / sizeof(binning[0]))-1,binning);

    TH1F *hprompt    = MakeTH1("hprompt",   7, binning2);
	  TH1F *hnonprompt = MakeTH1("hnonprompt" , 7, binning2);
	  TH1F *hinclu     = MakeTH1("hinclu",     7, binning2);
    TH2F* hRPM_preliminary = new TH2F("hWRPM_preliminary","",60,1,20,60,1,20); hRPM_preliminary->Sumw2();

    TH1F* hOA_Data1 = new TH1F("oa_data1","",50,0,180); hOA_Data1->Sumw2();
     TH1F* hOA_Data2 = new TH1F("oa_data2","",50,0,180); hOA_Data2->Sumw2();
    TH1F* hOA_MC = new TH1F("oa_mc","",50,0,180); hOA_MC->Sumw2();

	#if 1
	//Loop over tracks
	//*****************************************************

	//!
	const Int_t nTracks = Pair->GetEntriesFast();
	const Int_t nEvents = EventTree->GetEntriesFast();
	if (nEvents != nTracks) { cout <<"EventTree and eXiTree does NOT match! Stop.\n"; return; }

		Int_t test = 0;

    if(IsMC){
      for(int i =0; i< MCXicTree->GetEntriesFast(); i++){
        MCXicTree->GetEntry(i);

        if(fabs(mcGenXic0pT)!=9999 && fabs(mcGenXic0rap)<0.5 ){
          hMCGenInclusiveXic0_woW->Fill(mcGenXic0pT);
          hMCGenInclusiveXic0_W->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
          if(mcGencflag == 1)hMCGenPromptXic0_woW->Fill(mcGenXic0pT);
          if(mcGencflag == 1)hMCGenPromptXic0_W->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
          if(mcGenbflag == 1)hMCGenFeeddowmXic0_woW->Fill(mcGenXic0pT);
          if(mcGenbflag == 1)hMCGenFeeddowmXic0_W->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
        }
        if(fabs(mcGenXic0pT)!=9999 && fabs(mcGenXic0rap)<0.8 ){
          hMCGenInclusiveXic0_W_rap08->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
          if(mcGencflag == 1)hMCGenPromptXic0_W_rap08->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
          if(mcGenbflag == 1)hMCGenFeeddowmXic0_W_rap08->Fill(mcGenXic0pT,fWeightFit->Eval(mcGenXic0pT));
        }
      }
    }

	for (Int_t i=0; i<nTracks; i++)
	{
		//!
		//Eventwise cut, kimc
		//+++++++++++++++++++++++++++++++++++++++

		EventTree->GetEntry(i);

		//Sanity check by run number
		if ( (fRunNumber < 252000) || (fRunNumber > 295000) )
		{
			cout <<Form("Invalid RunNumber %f detected in track %i: skip.", fRunNumber, i) <<endl;
			continue;
		}

		//Trigger, apply only to the data
		Bool_t TrigFired = false;
		if      ( !strcmp(TRIG, "MB")    && (fTrigBit & TrigMB)    ) TrigFired = true;
		else if ( !strcmp(TRIG, "HMV0")  && (fTrigBit & TrigHMV0)  ) TrigFired = true;
		else if ( !strcmp(TRIG, "HMSPD") && (fTrigBit & TrigHMSPD) ) TrigFired = true;
		else if ( !strcmp(TRIG, "HMOR") )
		{
			if ( (fTrigBit & TrigHMV0) || (fTrigBit & TrigHMSPD) ) TrigFired = true;
		}
		if (IsMC==false && TrigFired==false) continue;

		//Multiplicity percentile
		Float_t tempMP = 0;
		if      ( !strcmp(TRIG, "MB")    ) tempMP = fCentrality;
		else if ( !strcmp(TRIG, "HMV0")  ) tempMP = fCentrality;
		else if ( !strcmp(TRIG, "HMSPD") ) tempMP = fCentralSPD;
		else if ( !strcmp(TRIG, "HMOR")  )
		{
			if ( (fTrigBit & TrigHMV0) && (fTrigBit & TrigHMSPD) ) tempMP = fCentralSPD;
			else if ( fTrigBit & TrigHMSPD ) tempMP = fCentralSPD;
			else if ( fTrigBit & TrigHMV0 ) tempMP = fCentrality;
		}
		if ( (tempMP < MultPerc[0]) || (tempMP > MultPerc[1]) ) continue;

		test++;
		//eXi pair - flags
		//+++++++++++++++++++++++++++++++++++++++

		Pair->GetEntry(i);

		if (fabs(Massv - 1.32171) > 0.008) continue; //Xi mass tolerance
		if (In_Mass < 1.3)                 continue; //pair mass low limit
		if (fabs(pTe == 999))              continue; //dummy tree reject

		//Bool_t isparticle = (vcharge>0)?kFALSE:kTRUE; //kimc: what's the purpose? Masking it
		Float_t VL_e_nsigma_cut = -4.3 + (1.17*pTe) - (0.094*pTe*pTe);
		Float_t  L_e_nsigma_cut = -4.1 + (1.17*pTe) - (0.094*pTe*pTe); ///need to modify
		Float_t  S_e_nsigma_cut = -3.9 + (1.17*pTe) - (0.094*pTe*pTe);
		Float_t  T_e_nsigma_cut = -3.7 + (1.17*pTe) - (0.094*pTe*pTe); ///need to modify
		Float_t VT_e_nsigma_cut = -3.5 + (1.15*pTe) - (0.09 *pTe*pTe);
		if (pTe >= 5)
		{
			VL_e_nsigma_cut = -4.3 + (1.17*5) - (0.094*25);
			 L_e_nsigma_cut = -4.1 + (1.17*5) - (0.094*25);
			 S_e_nsigma_cut = -3.9 + (1.17*5) - (0.094*25);
			 T_e_nsigma_cut = -3.7 + (1.17*5) - (0.094*25);
			VT_e_nsigma_cut = -3.5 + (1.15*5) - (0.09 *25);
		}

		Bool_t Xi_Topology_VLoose_flag = kFALSE;
		Bool_t Xi_Topology_Loose_flag = kFALSE;
		Bool_t Xi_Topology_Stand_flag = kFALSE;
		Bool_t Xi_Topology_Tight_flag = kFALSE;
		Bool_t Xi_Topology_VTight_flag = kFALSE;
		if ( V0DecayLength>1.1 && CascDecayLength>0.5 && DCABachToPrimVertex>0.05 && DCANegToPrimVertex>0.05 &&
			 DCAPosToPrimVertex>0.05 && V0CosineOfPoiningAngleXi>0.98 && XiCosineOfPoiningAngle>0.98 &&
			 DCAV0ToPrimVertex>0.05) Xi_Topology_VLoose_flag = kTRUE;
		if ( V0DecayLength>1.55 && CascDecayLength>0.5 && DCABachToPrimVertex>0.05 && DCANegToPrimVertex>0.07 &&
			 DCAPosToPrimVertex>0.07 && V0CosineOfPoiningAngleXi>0.981 && XiCosineOfPoiningAngle>0.981 &&
			 DCAV0ToPrimVertex>0.07 ) Xi_Topology_Loose_flag = kTRUE;
		if ( V0DecayLength>2.67 && CascDecayLength>0.5 && DCABachToPrimVertex>0.07 && DCANegToPrimVertex>0.095 &&
			 DCAPosToPrimVertex>0.095 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&
			 DCAV0ToPrimVertex>0.09 ) Xi_Topology_Stand_flag = kTRUE;
		if ( V0DecayLength>3.6 && CascDecayLength>0.6 && DCABachToPrimVertex>0.09 && DCANegToPrimVertex>0.11 &&
			 DCAPosToPrimVertex>0.11 && V0CosineOfPoiningAngleXi>0.9839 && XiCosineOfPoiningAngle>0.9839 &&
			 DCAV0ToPrimVertex>0.10 ) Xi_Topology_Tight_flag = kTRUE;
		if ( V0DecayLength>4.39 && CascDecayLength>0.7 && DCABachToPrimVertex>0.11 && DCANegToPrimVertex>0.13 &&
			 DCAPosToPrimVertex>0.13 && V0CosineOfPoiningAngleXi>0.985 && XiCosineOfPoiningAngle>0.985 &&
			 DCAV0ToPrimVertex>0.12) Xi_Topology_VTight_flag = kTRUE;

		Bool_t Xi_Recon_VLoose_flag = kFALSE;
		Bool_t Xi_Recon_Loose_flag = kFALSE;
		Bool_t Xi_Recon_Stand_flag = kFALSE;
		Bool_t Xi_Recon_Tight_flag = kFALSE;
		Bool_t Xi_Recon_VTight_flag = kFALSE;
		if ( pion_findable>0 && proton_findable>0 && bpion_findable>0) //Update: Require all denominators > 0, kimc
		{
      if ( pion_crossedratio/pion_findable>0.74 && proton_crossedratio/proton_findable>0.74 &&
				 bpion_crossedratio/bpion_findable>0.74 && pion_crossedratio>66 &&
				 proton_crossedratio>66 && bpion_crossedratio>66 ) Xi_Recon_VLoose_flag = kTRUE;
			if ( pion_crossedratio/pion_findable>0.77 && proton_crossedratio/proton_findable>0.77 &&
				 bpion_crossedratio/bpion_findable>0.77 && pion_crossedratio>66 &&
				 proton_crossedratio>66 && bpion_crossedratio>66 ) Xi_Recon_Loose_flag = kTRUE;
			if ( pion_crossedratio/pion_findable>0.77 && proton_crossedratio/proton_findable>0.77 &&
				 bpion_crossedratio/bpion_findable>0.77 && pion_crossedratio>70 &&
				 proton_crossedratio>70 && bpion_crossedratio>70 ) Xi_Recon_Stand_flag = kTRUE;
			if ( pion_crossedratio/pion_findable>0.79 && proton_crossedratio/proton_findable>0.79 &&
				 bpion_crossedratio/bpion_findable>0.79 && pion_crossedratio>70 &&
				 proton_crossedratio>70 && bpion_crossedratio>70 ) Xi_Recon_Tight_flag = kTRUE;
			if ( pion_crossedratio/pion_findable>0.79 && proton_crossedratio/proton_findable>0.79 &&
				 bpion_crossedratio/bpion_findable>0.79 && pion_crossedratio>74 &&
				 proton_crossedratio>74 && bpion_crossedratio>74 ) Xi_Recon_VTight_flag = kTRUE;
		}

		Bool_t e_Recon_VLoose_flag = kFALSE;
		Bool_t e_Recon_Loose_flag = kFALSE;
		Bool_t e_Recon_Stand_flag = kFALSE;
		Bool_t e_Recon_Tight_flag = kFALSE;
		Bool_t e_Recon_VTight_flag = kFALSE;
		if (e_findable > 0) //Update: Require all denominators > 0, kimc
		{
			if ( e_crossedratio>65 && TPCPIDCluster>40 && e_crossedratio/e_findable>0.7 &&
				 ITSCluster>=3) e_Recon_VLoose_flag = kTRUE;
			if ( e_crossedratio>65 && TPCPIDCluster>45 && e_crossedratio/e_findable>0.75 &&
				 ITSCluster>=3) e_Recon_Loose_flag = kTRUE;
			if ( e_crossedratio>70 && TPCPIDCluster>50 && e_crossedratio/e_findable>0.8 &&
				 ITSCluster>=3) e_Recon_Stand_flag = kTRUE;
			if ( e_crossedratio>75 && TPCPIDCluster>55 && e_crossedratio/e_findable>0.85 &&
				 ITSCluster>=3) e_Recon_Tight_flag = kTRUE;
			if ( e_crossedratio>85 && TPCPIDCluster>60 && e_crossedratio/e_findable>0.9 &&
				 ITSCluster>=3) e_Recon_VTight_flag = kTRUE;
		}

		Bool_t e_PID_VLoose_flag = kFALSE;
		Bool_t e_PID_Loose_flag = kFALSE;
		Bool_t e_PID_Stand_flag = kFALSE;
		Bool_t e_PID_Tight_flag = kFALSE;
		Bool_t e_PID_VTight_flag = kFALSE;
		if (fabs(nSigmaTOF)<=3 && nSigmaTPC>=VL_e_nsigma_cut && nSigmaTPC<=3) e_PID_VLoose_flag = kTRUE;
		if (fabs(nSigmaTOF)<=3 && nSigmaTPC>=L_e_nsigma_cut && nSigmaTPC<=3) e_PID_Loose_flag = kTRUE;
		if (fabs(nSigmaTOF)<=3 && nSigmaTPC>=S_e_nsigma_cut && nSigmaTPC<=3) e_PID_Stand_flag = kTRUE;
		if (fabs(nSigmaTOF)<=3 && nSigmaTPC>=T_e_nsigma_cut && nSigmaTPC<=3) e_PID_Tight_flag = kTRUE;
		if (fabs(nSigmaTOF)<=3 && nSigmaTPC>=VT_e_nsigma_cut && nSigmaTPC<=3) e_PID_VTight_flag = kTRUE;

		Bool_t OPAngle_Loose_flag = kFALSE;
		Bool_t OPAngle_Stand_flag = kFALSE;
		Bool_t OPAngle_Tight_flag = kFALSE;
		if (cosoa>cos( 90*(3.141592/180))) OPAngle_Loose_flag = kTRUE; //NO LOOSE CUT FOR OPENGING ANGLE
		if (cosoa>cos( 90*(3.141592/180))) OPAngle_Stand_flag = kTRUE;
		if (cosoa>cos( 70*(3.141592/180))) OPAngle_Tight_flag = kTRUE;

		Bool_t PairMass_Loose_flag = kFALSE;
		Bool_t PairMass_Stand_flag = kFALSE;
		Bool_t PairMass_Tight_flag = kFALSE;
		if (In_Mass<2.7) PairMass_Loose_flag = kTRUE;
		if (In_Mass<2.5) PairMass_Stand_flag = kTRUE;
		if (In_Mass<2.3) PairMass_Tight_flag = kTRUE;

		//Fill
		//+++++++++++++++++++++++++++++++++++++++

		//Fill1
		if (echarge*vcharge<0 && e_minmass>0.05) //Right Sign
        {
			if (Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && e_Recon_Stand_flag) hMassPtRS->Fill(In_Mass,Pt);

            if (Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (e_Recon_Tight_flag) hPtRS_TeRec->Fill(Pt);
                if (e_Recon_Stand_flag) hPtRS_SeRec->Fill(Pt);
                if (e_Recon_Loose_flag) hPtRS_LeRec->Fill(Pt);
                if (e_Recon_VTight_flag) hPtRS_VTeRec->Fill(Pt);
                if (e_Recon_VLoose_flag) hPtRS_VLeRec->Fill(Pt);
            }
            if (e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (Xi_Recon_Tight_flag) hPtRS_TXiRec->Fill(Pt);
                if (Xi_Recon_Stand_flag) hPtRS_SXiRec->Fill(Pt);
                if (Xi_Recon_Loose_flag) hPtRS_LXiRec->Fill(Pt);
                if (Xi_Recon_VTight_flag) hPtRS_VTXiRec->Fill(Pt);
                if (Xi_Recon_VLoose_flag) hPtRS_VLXiRec->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (e_PID_Tight_flag) hPtRS_TePID->Fill(Pt);
                if (e_PID_Stand_flag) hPtRS_SePID->Fill(Pt);
                if (e_PID_Loose_flag) hPtRS_LePID->Fill(Pt);
                if (e_PID_VTight_flag) hPtRS_VTePID->Fill(Pt);
                if (e_PID_VLoose_flag) hPtRS_VLePID->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (Xi_Topology_Tight_flag) hPtRS_TXiPID->Fill(Pt);
                if (Xi_Topology_Stand_flag) hPtRS_SXiPID->Fill(Pt);
                if (Xi_Topology_Loose_flag) hPtRS_LXiPID->Fill(Pt);
                if (Xi_Topology_VTight_flag) hPtRS_VTXiPID->Fill(Pt);
                if (Xi_Topology_VLoose_flag) hPtRS_VLXiPID->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				Xi_Topology_Stand_flag && PairMass_Stand_flag)
			{
                if (OPAngle_Tight_flag) hPtRS_TOA->Fill(Pt);
                if (OPAngle_Stand_flag) hPtRS_SOA->Fill(Pt);
                if (OPAngle_Loose_flag) hPtRS_LOA->Fill(Pt);
                hOA_Data1->Fill((acos(cosoa)/3.141592)*180);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				Xi_Topology_Stand_flag && OPAngle_Stand_flag)
			{
                if (PairMass_Tight_flag) hPtRS_TeIM->Fill(Pt);
                if (PairMass_Stand_flag) hPtRS_SeIM->Fill(Pt);
                if (PairMass_Loose_flag) hPtRS_LeIM->Fill(Pt);
            }
        }//RS, Fill1

		//Fill2
        if (echarge*vcharge>0 && e_minmass>0.05) //Wrong Sign
        {
            if (Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && e_Recon_Stand_flag) hMassPtWS->Fill(In_Mass,Pt);

            if (Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (e_Recon_Tight_flag) hPtWS_TeRec->Fill(Pt);
                if (e_Recon_Stand_flag) hPtWS_SeRec->Fill(Pt);
                if (e_Recon_Loose_flag) hPtWS_LeRec->Fill(Pt);
                if (e_Recon_VTight_flag) hPtWS_VTeRec->Fill(Pt);
                if (e_Recon_VLoose_flag) hPtWS_VLeRec->Fill(Pt);
            }
            if (e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (Xi_Recon_Tight_flag) hPtWS_TXiRec->Fill(Pt);
                if (Xi_Recon_Stand_flag) hPtWS_SXiRec->Fill(Pt);
                if (Xi_Recon_Loose_flag) hPtWS_LXiRec->Fill(Pt);
                if (Xi_Recon_VTight_flag) hPtWS_VTXiRec->Fill(Pt);
                if (Xi_Recon_VLoose_flag) hPtWS_VLXiRec->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (e_PID_Tight_flag) hPtWS_TePID->Fill(Pt);
                if (e_PID_Stand_flag) hPtWS_SePID->Fill(Pt);
                if (e_PID_Loose_flag) hPtWS_LePID->Fill(Pt);
                if (e_PID_VTight_flag) hPtWS_VTePID->Fill(Pt);
                if (e_PID_VLoose_flag) hPtWS_VLePID->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				OPAngle_Stand_flag && PairMass_Stand_flag)
			{
                if (Xi_Topology_Tight_flag) hPtWS_TXiPID->Fill(Pt);
                if (Xi_Topology_Stand_flag) hPtWS_SXiPID->Fill(Pt);
                if (Xi_Topology_Loose_flag) hPtWS_LXiPID->Fill(Pt);
                if (Xi_Topology_VTight_flag) hPtWS_VTXiPID->Fill(Pt);
                if (Xi_Topology_VLoose_flag) hPtWS_VLXiPID->Fill(Pt);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				Xi_Topology_Stand_flag && PairMass_Stand_flag)
			{
                if (OPAngle_Tight_flag) hPtWS_TOA->Fill(Pt);
                if (OPAngle_Stand_flag) hPtWS_SOA->Fill(Pt);
                if (OPAngle_Loose_flag) hPtWS_LOA->Fill(Pt);
				hOA_Data2->Fill((acos(cosoa)/3.141592)*180);
            }
            if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				Xi_Topology_Stand_flag && OPAngle_Stand_flag)
			{
                if (PairMass_Tight_flag) hPtWS_TeIM->Fill(Pt);
                if (PairMass_Stand_flag) hPtWS_SeIM->Fill(Pt);
                if (PairMass_Loose_flag) hPtWS_LeIM->Fill(Pt);
            }
        }//WS, Fill2

		//Fill3
        if (echarge*vcharge>0 && PairMass_Stand_flag && Xi_Recon_Stand_flag &&
			e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag)
		{
            if (e_Recon_Tight_flag) hpre_eRec_tight_de->Fill(Pt);
            if (e_Recon_Stand_flag) hpre_eRec_stand_de->Fill(Pt);
            if (e_Recon_Loose_flag) hpre_eRec_loose_de->Fill(Pt);
            if (e_Recon_VTight_flag) hpre_eRec_vtight_de->Fill(Pt);
            if (e_Recon_VLoose_flag) hpre_eRec_vloose_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (e_Recon_Tight_flag) hpre_eRec_tight_nu->Fill(Pt);
                if (e_Recon_Stand_flag) hpre_eRec_stand_nu->Fill(Pt);
                if (e_Recon_Loose_flag) hpre_eRec_loose_nu->Fill(Pt);
                if (e_Recon_VTight_flag) hpre_eRec_vtight_nu->Fill(Pt);
                if (e_Recon_VLoose_flag) hpre_eRec_vloose_nu->Fill(Pt);
            }
        }//Fill3

		//Fill4
        if (echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag &&
			e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag)
		{
            if (Xi_Recon_Loose_flag) hpre_XiRec_loose_de->Fill(Pt);
            if (Xi_Recon_Stand_flag) hpre_XiRec_stand_de->Fill(Pt);
            if (Xi_Recon_Tight_flag) hpre_XiRec_tight_de->Fill(Pt);
            if (Xi_Recon_VLoose_flag) hpre_XiRec_vloose_de->Fill(Pt);
            if (Xi_Recon_VTight_flag) hpre_XiRec_vtight_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (Xi_Recon_Loose_flag) hpre_XiRec_loose_nu->Fill(Pt);
                if (Xi_Recon_Stand_flag) hpre_XiRec_stand_nu->Fill(Pt);
                if (Xi_Recon_Tight_flag) hpre_XiRec_tight_nu->Fill(Pt);
                if (Xi_Recon_VLoose_flag) hpre_XiRec_vloose_nu->Fill(Pt);
                if (Xi_Recon_VTight_flag) hpre_XiRec_vtight_nu->Fill(Pt);
            }
        }//Fill4

		//Fill5
        if (echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag &&
			Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag)
		{
            if (e_PID_Loose_flag) hpre_ePID_loose_de->Fill(Pt);
            if (e_PID_Stand_flag) hpre_ePID_stand_de->Fill(Pt);
            if (e_PID_Tight_flag) hpre_ePID_tight_de->Fill(Pt);
            if (e_PID_VLoose_flag) hpre_ePID_vloose_de->Fill(Pt);
            if (e_PID_VTight_flag) hpre_ePID_vtight_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (e_PID_Loose_flag) hpre_ePID_loose_nu->Fill(Pt);
                if (e_PID_Stand_flag) hpre_ePID_stand_nu->Fill(Pt);
                if (e_PID_Tight_flag) hpre_ePID_tight_nu->Fill(Pt);
                if (e_PID_VLoose_flag) hpre_ePID_vloose_nu->Fill(Pt);
                if (e_PID_VTight_flag) hpre_ePID_vtight_nu->Fill(Pt);
            }
        }//Fill5

		//Fill6
        if (echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag &&
			Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag)
		{
            if (Xi_Topology_Loose_flag) hpre_XiPID_loose_de->Fill(Pt);
            if (Xi_Topology_Stand_flag) hpre_XiPID_stand_de->Fill(Pt);
            if (Xi_Topology_Tight_flag) hpre_XiPID_tight_de->Fill(Pt);
            if (Xi_Topology_VLoose_flag) hpre_XiPID_vloose_de->Fill(Pt);
            if (Xi_Topology_VTight_flag) hpre_XiPID_vtight_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (Xi_Topology_Loose_flag) hpre_XiPID_loose_nu->Fill(Pt);
                if (Xi_Topology_Stand_flag) hpre_XiPID_stand_nu->Fill(Pt);
                if (Xi_Topology_Tight_flag) hpre_XiPID_tight_nu->Fill(Pt);
                if (Xi_Topology_VLoose_flag) hpre_XiPID_vloose_nu->Fill(Pt);
                if (Xi_Topology_VTight_flag) hpre_XiPID_vtight_nu->Fill(Pt);
            }
        }//Fill6

		//Fill7
        if (echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag &&
			Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag)
		{
            if (OPAngle_Loose_flag) hpre_OA_loose_de->Fill(Pt);
            if (OPAngle_Stand_flag) hpre_OA_stand_de->Fill(Pt);
            if (OPAngle_Tight_flag) hpre_OA_tight_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (OPAngle_Loose_flag) hpre_OA_loose_nu->Fill(Pt);
                if (OPAngle_Stand_flag) hpre_OA_stand_nu->Fill(Pt);
                if (OPAngle_Tight_flag) hpre_OA_tight_nu->Fill(Pt);
            }
        }//Fill7

		//Fill8
        if (echarge*vcharge>0 && OPAngle_Stand_flag && e_Recon_Stand_flag &&
			Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag)
		{
            if (PairMass_Loose_flag) hpre_IM_loose_de->Fill(Pt);
            if (PairMass_Stand_flag) hpre_IM_stand_de->Fill(Pt);
            if (PairMass_Tight_flag) hpre_IM_tight_de->Fill(Pt);
            if (e_minmass_ss>0.05)
			{
                if (PairMass_Loose_flag) hpre_IM_loose_nu->Fill(Pt);
                if (PairMass_Stand_flag) hpre_IM_stand_nu->Fill(Pt);
                if (PairMass_Tight_flag) hpre_IM_tight_nu->Fill(Pt);
            }
        }//WS, Fill8

		//MC selection
		//+++++++++++++++++++++++++++++++++++++++

		if (IsMC)
		{
			const float nullV = 9999;
			MCTree->GetEntry(i);

			//Fill b
			if (Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag &&
				PairMass_Stand_flag)
			{
				if (e_Recon_VTight_flag)
				{
					if ((fabs(mcptXic0)!=nullV) && mcptXic0>0) //kimc: mcptXic0>0 to avoid floating point exception
					{
						hRPM_eRec_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_eRec_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_eRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_eRec_vtight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_eRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_eRec_vtight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_eRec_vtight->Fill(mcXib);
						hMCRecoLevPairXib_eRec_vtight->Fill(XibeXi);
					}
				}
				if (e_Recon_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_eRec_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_eRec_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_eRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_eRec_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_eRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_eRec_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_eRec_tight->Fill(mcXib);
						hMCRecoLevPairXib_eRec_tight->Fill(XibeXi);
					}
				}
				if (e_Recon_Stand_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_eRec_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_eRec_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_eRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_eRec_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_eRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));

            hevseXi->Fill(pTe,Pt);
            heXivsXic0->Fill(Pt,mcptXic0);

					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_eRec_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_eRec_stand->Fill(mcXib);
						hMCRecoLevPairXib_eRec_stand->Fill(XibeXi);
					}
				}
				if (e_Recon_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_eRec_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_eRec_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_eRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_eRec_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_eRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_eRec_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_eRec_loose->Fill(mcXib);
						hMCRecoLevPairXib_eRec_loose->Fill(XibeXi);
					}
				}
				if (e_Recon_VLoose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_eRec_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_eRec_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_eRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_eRec_vloose->Fill(Pt,fWeightFit->Eval(mcptXic0));
					    if (mcc_flag==1) hMCRecoLevPromptXic0_eRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_eRec_vloose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_eRec_vloose->Fill(mcXib);
						hMCRecoLevPairXib_eRec_vloose->Fill(XibeXi);
					}
				}
			}//Fill b

			//Fill c
			if (e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag &&
				PairMass_Stand_flag)
			{
				if (Xi_Recon_VTight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiRec_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiRec_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiRec_vtight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPromptXic0_XiRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiRec_vtight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiRec_vtight->Fill(mcXib);
						hMCRecoLevPairXib_XiRec_vtight->Fill(XibeXi);
					}
				}
				if (Xi_Recon_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiRec_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiRec_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiRec_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiRec_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiRec_tight->Fill(mcXib);
						hMCRecoLevPairXib_XiRec_tight->Fill(XibeXi);
					}
				}
				if (Xi_Recon_Stand_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiRec_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiRec_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiRec_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiRec_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiRec_stand->Fill(mcXib);
						hMCRecoLevPairXib_XiRec_stand->Fill(XibeXi);
					}
				}
				if (Xi_Recon_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiRec_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiRec_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiRec_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiRec_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiRec_loose->Fill(mcXib);
						hMCRecoLevPairXib_XiRec_loose->Fill(XibeXi);
					}
				}
				if (Xi_Recon_VLoose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiRec_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiRec_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiRec_vloose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiRec_vloose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiRec_vloose->Fill(mcXib);
						hMCRecoLevPairXib_XiRec_vloose->Fill(XibeXi);
					}
				}
			}//Fill c

			//Fill d
			if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag &&
				PairMass_Stand_flag)
			{
				if (e_PID_VTight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_ePID_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_ePID_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_ePID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_ePID_vtight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_ePID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_ePID_vtight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_ePID_vtight->Fill(mcXib);
						hMCRecoLevPairXib_ePID_vtight->Fill(XibeXi);
					}
				}
				if (e_PID_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_ePID_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_ePID_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_ePID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_ePID_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_ePID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_ePID_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_ePID_tight->Fill(mcXib);
						hMCRecoLevPairXib_ePID_tight->Fill(XibeXi);
					}
				}
				if (e_PID_Stand_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_ePID_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_ePID_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_ePID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_ePID_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_ePID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_ePID_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_ePID_stand->Fill(mcXib);
						hMCRecoLevPairXib_ePID_stand->Fill(XibeXi);
					}
				}
				if (e_PID_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_ePID_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_ePID_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_ePID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_ePID_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_ePID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_ePID_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_ePID_loose->Fill(mcXib);
						hMCRecoLevPairXib_ePID_loose->Fill(XibeXi);
					}
				}
				if (e_PID_VLoose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_ePID_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_ePID_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_ePID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_ePID_vloose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_ePID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_ePID_vloose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_ePID_vloose->Fill(mcXib);
						hMCRecoLevPairXib_ePID_vloose->Fill(XibeXi);
					}
				}
			}//Fill d

			//Fill e
			if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag &&
				PairMass_Stand_flag)
			{
				if (Xi_Topology_VTight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiPID_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiPID_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiPID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiPID_vtight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiPID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiPID_vtight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiPID_vtight->Fill(mcXib);
						hMCRecoLevPairXib_XiPID_vtight->Fill(XibeXi);
					}
				}
				if (Xi_Topology_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiPID_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiPID_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiPID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiPID_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiPID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiPID_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiPID_tight->Fill(mcXib);
						hMCRecoLevPairXib_XiPID_tight->Fill(XibeXi);
					}
				}
				if (Xi_Topology_Stand_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiPID_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiPID_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiPID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiPID_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiPID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiPID_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiPID_stand->Fill(mcXib);
						hMCRecoLevPairXib_XiPID_stand->Fill(XibeXi);
					}
				}
				if (Xi_Topology_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiPID_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiPID_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiPID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiPID_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiPID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiPID_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiPID_loose->Fill(mcXib);
						hMCRecoLevPairXib_XiPID_loose->Fill(XibeXi);
					}
				}
				if (Xi_Topology_VLoose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_XiPID_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_XiPID_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_XiPID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_XiPID_vloose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_XiPID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_XiPID_vloose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_XiPID_vloose->Fill(mcXib);
						hMCRecoLevPairXib_XiPID_vloose->Fill(XibeXi);
					}
				}
			}//Fill e

			//Fill f
			if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				PairMass_Stand_flag)
			{
				if (OPAngle_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_op_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_op_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_op_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_op_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_OA_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_op_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_op_tight->Fill(mcXib);
						hMCRecoLevPairXib_op_tight->Fill(XibeXi);
					}
				}
				if (OPAngle_Stand_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_op_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_op_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_op_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_op_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_OA_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_op_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_op_stand->Fill(mcXib);
						hMCRecoLevPairXib_op_stand->Fill(XibeXi);
					}
				}
				if (OPAngle_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_op_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_op_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_op_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_op_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_OA_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_op_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_op_loose->Fill(mcXib);
						hMCRecoLevPairXib_op_loose->Fill(XibeXi);
					}
				}

        if(fabs(mcptXic0)!=nullV && mcptXic0>0) hOA_MC->Fill((acos(cosoa)/3.141592)*180);
			}//Fill f

			//Fill g
			if (e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag &&
				OPAngle_Stand_flag)
			{
				if (PairMass_Tight_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_im_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_im_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_im_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_im_tight->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_IM_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_im_tight_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_im_tight->Fill(mcXib);
						hMCRecoLevPairXib_im_tight->Fill(XibeXi);
					}
				}
				if (PairMass_Stand_flag && mcptXic0>0)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_im_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_im_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_im_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_im_stand->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_IM_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=nullV)
					{
						hRPM_im_stand_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_im_stand->Fill(mcXib);
						hMCRecoLevPairXib_im_stand->Fill(XibeXi);
					}
				}
				if (PairMass_Loose_flag)
				{
					if (fabs(mcptXic0)!=nullV && mcptXic0>0)
					{
						hRPM_im_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
						hRPM_im_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevXic0_im_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
						hMCRecoLevPair_im_loose->Fill(Pt,fWeightFit->Eval(mcptXic0));
						if (mcc_flag==1) hMCRecoLevPromptXic0_IM_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
					}
					if (echarge*vcharge>0 && fabs(mcXib)!=9999)
					{
						hRPM_im_loose_Xib->Fill(mcXib,XibeXi);
						hMCRecoLevXib_im_loose->Fill(mcXib);
						hMCRecoLevPairXib_im_loose->Fill(XibeXi);
					}
				}
			}//Fill g

			if (echarge*vcharge<0 && e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag &&
				Xi_Topology_Stand_flag)
			{
				if (mcc_flag == 1) { if (OPAngle_Stand_flag && PairMass_Stand_flag) hprompt->Fill(mcptXic0); }
				if (mcb_flag == 1) { if (OPAngle_Stand_flag && PairMass_Stand_flag) hnonprompt->Fill(mcptXic0); }

				//kimc: added mcptXic0 != nullV, otherwise it causes floating exception error during eval
				if (fabs(mcptXic0)!=nullV && mcptXic0>0 && OPAngle_Stand_flag && PairMass_Stand_flag)
				{
					hinclu->Fill(mcptXic0, fWeightFit->Eval(mcptXic0));
				}
			}
		}//IsMC
	}//i, loop over tracks
	#endif

	//*****************************************************

    TH1F* hMCTrueXic0 = new TH1F;
    TH1F* hMCGenLevXic0_inc = new TH1F;
    TH1F* hMCGenLevXic0_incW = new TH1F;
    TH1F* hMCGenLevXic0_inc8 = new TH1F;
    TH1F* hMCGenLevXic0_incW8 = new TH1F;
    TH1F* hMCGenLevXic01_p = new TH1F;
    TH1F* hMCGenLevXic02_p = new TH1F;
    TH1F* hMCGenLevXic_p = new TH1F;
    TH1F* hMCGenLevXic01_np = new TH1F;
    TH1F* hMCGenLevXic02_np = new TH1F;
    TH1F* hMCRecoLevXic0_inc = new TH1F;
    TH1F* hMCRecoLevXic0_np = new TH1F;
    TH1F* hMCRecoLevXic0_p = new TH1F;
    TH1F* eff_inc = new TH1F;
    TH1F* eff_inc_2 = new TH1F;
    TH1F* eff_p = new TH1F;
    TH1F* eff_np = new TH1F;
    TH1F* NonPromptXicRap = new TH1F;
    TH1F* PromptXicRap = new TH1F;
    TH1F* XicRap = new TH1F;
    TH2F* XicRap2 = new TH2F;
    TH1F* XibGen = new TH1F;
    TH1F* XibGen05 = new TH1F;

    if (IsMC)
    {
		//generation level for 0.5 rap., ptbinning1
        hMCGenLevXic0_inc   = (TH1F*)hist->FindObject("hTrueXic0")->Clone("hMCGenLevXic0_inc");
        //hMCGenLevXic0_incW  = (TH1F*)hist->FindObject("hTrueXic0W")->Clone("hMCGenLevXic0_incW");
        //hMCGenLevXic0_inc8  = (TH1F*)hist->FindObject("hTrueXic0_rap8")->Clone("hMCGenLevXic0_inc8");
        //hMCGenLevXic0_incW8 = (TH1F*)hist->FindObject("hTrueXic0W_rap8")->Clone("hMCGenLevXic0_incW8");

		//generation level for 0.5 rap., ptbinning1
        hMCGenLevXic01_p  = (TH1F*)hist->FindObject("hXic0PtFromCharm1")->Clone("hMCGenLevXic01_p");
        hMCGenLevXic02_p  = (TH1F*)hist->FindObject("hXic0PtFromCharm2")->Clone("hMCGenLevXic02_p");
        hMCGenLevXic01_np = (TH1F*)hist->FindObject("hXic0PtFromBottom1")->Clone("hMCGenLevXic01_np");
        hMCGenLevXic02_np = (TH1F*)hist->FindObject("hXic0PtFromBottom2")->Clone("hMCGenLevXic02_np");

		//reconstruction level, ptbinning1
        hMCRecoLevXic0_inc = (TH1F*)hist->FindObject("hGenXic0Pt")->Clone("hMCRecoLevXic0_inc");

        hMCRecoLevXic0_np = (TH1F*)hist->FindObject("hGenXic0PtFromXib")->Clone("hMCRecoLevXic0_np");
        hMCRecoLevXic0_p = (TH1F*)hist->FindObject("hGenXic0PtFromXic")->Clone("hMCRecoLevXic0_p");
		//hMCGenLevXic_p = (TH1F*)hist->FindObject("hXic0PtFromCharmW")->Clone("hMCGenLevXic_pW"); //segfault - kimc

        XibGen   = (TH1F*)hist->FindObject("XibGen")->Clone("XibGen");
        XibGen05 = (TH1F*)hist->FindObject("XibGen05")->Clone("XibGen05");

		//! - not exist in old code?
        /*eff_inc = (TH1F*)hinclu->Clone("eff_inc");
        eff_inc->Divide(eff_inc,hMCGenLevXic0_inc,1,1,"b");

        eff_inc_2 = (TH1F*)hinclu->Clone("eff_inc_2");
        eff_inc_2->Divide(eff_inc_2,hMCGenLevXic0_incW,1,1,"b");

        eff_p = (TH1F*)hMCGenLevXic01_p->Clone("eff_p");
        eff_p->Add(hMCGenLevXic02_p,1);
        eff_p->Divide(hprompt,eff_p,1,1,"b");

        eff_np = (TH1F*)hMCGenLevXic01_np->Clone("eff_np");
        eff_np->Add(hMCGenLevXic02_np,1);
        eff_np->Divide(hnonprompt,eff_np,1,1,"b");*/

		//kimc - what's the purpose?
        TH1F* hPromptXic0 = (TH1F*)hMCGenLevXic01_p->Clone("hMCGenLevXic_p");
        hPromptXic0->Add(hMCGenLevXic02_p,2);

        NonPromptXicRap = (TH1F*)hist->FindObject("hNonPromptXicRap")->Clone("hNonPromptXicRap");
        PromptXicRap    = (TH1F*)hist->FindObject("hPromptXicRap")->Clone("hPromptXicRap");

        XicRap  = (TH1F*)hist->FindObject("hXicRap")->Clone("hXicRap");
        XicRap2 = (TH2F*)hist->FindObject("hXicPtRap")->Clone("hXicPtRap"); //! - not exist in old code
    }

	//*****************************************************

	#if 1
    TH1F *hXic0_Bayes_stand2 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand2");
    TH1F *hXic0_Bayes_stand3 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand3");
    TH1F *hXic0_Bayes_stand4 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand4");
    TH1F *hXic0_Bayes_stand5 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand5");
    TH1F *hXic0_Bayes_stand6 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand6");
    TH1F *hXic0_Bayes_stand7 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand7");
    TH1F *hXic0_Svd_stand3 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand3");
    TH1F *hXic0_Svd_stand4 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand4");
    TH1F *hXic0_Svd_stand5 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand5");

      TH1F *hXic0_BottomBaryon_default = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_BottomBaryon_default");
      TH1F *hXic0_BottomBaryon_variation1 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_BottomBaryon_variation1");
      TH1F *hXic0_BottomBaryon_variation2 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_BottomBaryon_variation2");
      TH1F *hXic0_Weighting_default = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Weighting_default");
      TH1F *hXic0_Weighting_variation1 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Weighting_variation1");
      TH1F *hXic0_Weighting_variation2 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Weighting_variation2");

    TH2F* hRPM_Unfold_stand2 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand2");
    TH2F* hRPM_Unfold_stand3 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand3");
    TH2F* hRPM_Unfold_stand4 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand4");
    TH2F* hRPM_Unfold_stand5 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand5");
    TH2F* hRPM_Unfold_stand6 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand6");
    TH2F* hRPM_Unfold_stand7 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand7");
    TH2F* hRPM_Unfold_stand2_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand2_un");
    TH2F* hRPM_Unfold_stand3_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand3_un");
    TH2F* hRPM_Unfold_stand4_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand4_un");
    TH2F* hRPM_Unfold_stand5_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand5_un");
    TH2F* hRPM_Unfold_stand6_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand6_un");
    TH2F* hRPM_Unfold_stand7_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Bayes_stand7_un");
    TH1F* hMCRecoLevXic0_Unfold_stand2 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand2");
    TH1F* hMCRecoLevPair_Unfold_stand2 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand2");
    TH1F* hMCRecoLevXic0_Unfold_stand3 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand3");
    TH1F* hMCRecoLevPair_Unfold_stand3 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand3");
    TH1F* hMCRecoLevXic0_Unfold_stand4 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand4");
    TH1F* hMCRecoLevPair_Unfold_stand4 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand4");
    TH1F* hMCRecoLevXic0_Unfold_stand5 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand5");
    TH1F* hMCRecoLevPair_Unfold_stand5 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand5");
    TH1F* hMCRecoLevXic0_Unfold_stand6 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand6");
    TH1F* hMCRecoLevPair_Unfold_stand6 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand6");
    TH1F* hMCRecoLevXic0_Unfold_stand7 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand7");
    TH1F* hMCRecoLevPair_Unfold_stand7 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand7");

      TH2F* hRPM_BottomBaryon_default = (TH2F*) hRPM_im_stand->Clone("hRPM_BottomBaryon_default");
      TH2F* hRPM_BottomBaryon_variation1 = (TH2F*) hRPM_im_stand->Clone("hRPM_BottomBaryon_variation1");
      TH2F* hRPM_BottomBaryon_variation2 = (TH2F*) hRPM_im_stand->Clone("hRPM_BottomBaryon_variation2");
      TH2F* hRPM_Weighting_default = (TH2F*) hRPM_im_stand->Clone("hRPM_Weighting_default");
      TH2F* hRPM_Weighting_variation1 = (TH2F*) hRPM_im_stand->Clone("hRPM_Weighting_variation1");
      TH2F* hRPM_Weighting_variation2 = (TH2F*) hRPM_im_stand->Clone("hRPM_Weighting_variation2");
      TH2F* hRPM_BottomBaryon_default_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_BottomBaryon_default_un");
      TH2F* hRPM_BottomBaryon_variation1_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_BottomBaryon_variation1_un");
      TH2F* hRPM_BottomBaryon_variation2_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_BottomBaryon_variation2_un");
      TH2F* hRPM_Weighting_default_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Weighting_default_un");
      TH2F* hRPM_Weighting_variation1_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Weighting_variation1_un");
      TH2F* hRPM_Weighting_variation2_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Weighting_variation2_un");
      TH1F* hMCRecoLevXic0_BottomBaryon_default = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_BottomBaryon_default");
      TH1F* hMCRecoLevPair_BottomBaryon_default = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_BottomBaryon_default");
      TH1F* hMCRecoLevXic0_BottomBaryon_variation1 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_BottomBaryon_variation1");
      TH1F* hMCRecoLevPair_BottomBaryon_variation1 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_BottomBaryon_variation1");
      TH1F* hMCRecoLevXic0_BottomBaryon_variation2 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_BottomBaryon_variation2");
      TH1F* hMCRecoLevPair_BottomBaryon_variation2 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_BottomBaryon_variation2");
      TH1F* hMCRecoLevXic0_Weighting_default = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Weighting_default");
      TH1F* hMCRecoLevPair_Weighting_default = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Weighting_default");
      TH1F* hMCRecoLevXic0_Weighting_variation1 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Weighting_variation1");
      TH1F* hMCRecoLevPair_Weighting_variation1 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Weighting_variation1");
      TH1F* hMCRecoLevXic0_Weighting_variation2 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Weighting_variation2");
      TH1F* hMCRecoLevPair_Weighting_variation2 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Weighting_variation2");

    TH1F* hRPM_Bayes_stand2_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand2_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand2     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand2");
    TH1F* hMCRecoLevPairXib_Bayes_stand2 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand2");
    TH1F* hRPM_Bayes_stand3_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand3_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand3     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand3");
    TH1F* hMCRecoLevPairXib_Bayes_stand3 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand3");
    TH1F* hRPM_Bayes_stand4_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand4_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand4     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand4");
    TH1F* hMCRecoLevPairXib_Bayes_stand4 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand4");
    TH1F* hRPM_Bayes_stand5_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand5_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand5     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand5");
    TH1F* hMCRecoLevPairXib_Bayes_stand5 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand5");
    TH1F* hRPM_Bayes_stand6_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand6_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand6     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand6");
    TH1F* hMCRecoLevPairXib_Bayes_stand6 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand6");
    TH1F* hRPM_Bayes_stand7_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand7_Xib");
    TH1F* hMCRecoLevXib_Bayes_stand7     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand7");
    TH1F* hMCRecoLevPairXib_Bayes_stand7 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand7");

    TH1F* hRPM_Svd_stand3_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Svd_stand3_Xib");
    TH1F* hMCRecoLevXib_Svd_stand3     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand3");
    TH1F* hMCRecoLevPairXib_Svd_stand3 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand3");
    TH1F* hRPM_Svd_stand4_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Svd_stand4_Xib");
    TH1F* hMCRecoLevXib_Svd_stand4     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand4");
    TH1F* hMCRecoLevPairXib_Svd_stand4 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand4");
    TH1F* hRPM_Svd_stand5_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Svd_stand5_Xib");
    TH1F* hMCRecoLevXib_Svd_stand5     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand5");
    TH1F* hMCRecoLevPairXib_Svd_stand5 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand5");

      TH1F* hRPM_BottomBaryon_default_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_BottomBaryon_default_Xib");
      TH1F* hMCRecoLevXib_BottomBaryon_default     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_BottomBaryon_default");
      TH1F* hMCRecoLevPairXib_BottomBaryon_default = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_BottomBaryon_default");
      TH1F* hRPM_BottomBaryon_variation1_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_BottomBaryon_variation1_Xib");
      TH1F* hMCRecoLevXib_BottomBaryon_variation1     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_BottomBaryon_variation1");
      TH1F* hMCRecoLevPairXib_BottomBaryon_variation1 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_BottomBaryon_variation1");
      TH1F* hRPM_BottomBaryon_variation2_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_BottomBaryon_variation2_Xib");
      TH1F* hMCRecoLevXib_BottomBaryon_variation2     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_BottomBaryon_variation2");
      TH1F* hMCRecoLevPairXib_BottomBaryon_variation2 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_BottomBaryon_variation2");
      TH1F* hRPM_Weighting_default_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Weighting_default_Xib");
      TH1F* hMCRecoLevXib_Weighting_default     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Weighting_default");
      TH1F* hMCRecoLevPairXib_Weighting_default = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Weighting_default");
      TH1F* hRPM_Weighting_variation1_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Weighting_variation1_Xib");
      TH1F* hMCRecoLevXib_Weighting_variation1     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Weighting_variation1");
      TH1F* hMCRecoLevPairXib_Weighting_variation1 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Weighting_variation1");
      TH1F* hRPM_Weighting_variation2_Xib          = (TH1F*)hRPM_op_stand_Xib->Clone("hRPM_Weighting_variation2_Xib");
      TH1F* hMCRecoLevXib_Weighting_variation2     = (TH1F*)hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Weighting_variation2");
      TH1F* hMCRecoLevPairXib_Weighting_variation2 = (TH1F*)hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Weighting_variation2");

    TH2F* hRPM_Unfold_stand3_s_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Svd_stand3_un");
    TH2F* hRPM_Unfold_stand4_s_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Svd_stand4_un");
    TH2F* hRPM_Unfold_stand5_s_un = (TH2F*) hRPM_im_stand_un->Clone("hRPM_Svd_stand5_un");
    TH2F* hRPM_Unfold_stand3_s = (TH2F*) hRPM_im_stand->Clone("hRPM_Svd_stand3");
    TH2F* hRPM_Unfold_stand4_s = (TH2F*) hRPM_im_stand->Clone("hRPM_Svd_stand4");
    TH2F* hRPM_Unfold_stand5_s = (TH2F*) hRPM_im_stand->Clone("hRPM_Svd_stand5");
    TH1F* hMCRecoLevXic0_Unfold_stand3_s = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Svd_stand3");
    TH1F* hMCRecoLevPair_Unfold_stand3_s = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Svd_stand3");
    TH1F* hMCRecoLevXic0_Unfold_stand4_s = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Svd_stand4");
    TH1F* hMCRecoLevPair_Unfold_stand4_s = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Svd_stand4");
    TH1F* hMCRecoLevXic0_Unfold_stand5_s = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Svd_stand5");
    TH1F* hMCRecoLevPair_Unfold_stand5_s = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Svd_stand5");

    TH1F *hRawPt_VTePID = (TH1F*) hPtRS_VTePID->Clone("hRawPt_ePID_vtight"); hRawPt_VTePID->Add(hPtWS_VTePID,-1);
    TH1F *hRawPt_TePID  = (TH1F*) hPtRS_TePID->Clone("hRawPt_ePID_tight");   hRawPt_TePID->Add(hPtWS_TePID,-1);
    TH1F *hRawPt_SePID  = (TH1F*) hPtRS_SePID->Clone("hRawPt_ePID_stand");   hRawPt_SePID->Add(hPtWS_SePID,-1);
    TH1F *hRawPt_LePID  = (TH1F*) hPtRS_LePID->Clone("hRawPt_ePID_loose");   hRawPt_LePID->Add(hPtWS_LePID,-1);
    TH1F *hRawPt_VLePID = (TH1F*) hPtRS_VLePID->Clone("hRawPt_ePID_vloose"); hRawPt_VLePID->Add(hPtWS_VLePID,-1);

    TH1F *hRawPt_VTXiPID = (TH1F*)hPtRS_VTXiPID->Clone("hRawPt_XiPID_vtight"); hRawPt_VTXiPID->Add(hPtWS_VTXiPID,-1);
    TH1F *hRawPt_TXiPID  = (TH1F*)hPtRS_TXiPID->Clone("hRawPt_XiPID_tight");   hRawPt_TXiPID->Add(hPtWS_TXiPID,-1);
    TH1F *hRawPt_SXiPID  = (TH1F*)hPtRS_SXiPID->Clone("hRawPt_XiPID_stand");   hRawPt_SXiPID->Add(hPtWS_SXiPID,-1);
    TH1F *hRawPt_LXiPID  = (TH1F*)hPtRS_LXiPID->Clone("hRawPt_XiPID_loose");   hRawPt_LXiPID->Add(hPtWS_LXiPID,-1);
    TH1F *hRawPt_VLXiPID = (TH1F*)hPtRS_VLXiPID->Clone("hRawPt_XiPID_vloose"); hRawPt_VLXiPID->Add(hPtWS_VLXiPID,-1);

    TH1F *hRawPt_VTXiPIDRec = (TH1F*)hPtRS_VTXiPIDRec->Clone("hRawPt_XiPIDRec_vtight");
	hRawPt_VTXiPIDRec->Add(hPtWS_VTXiPIDRec,-1);
    TH1F *hRawPt_TXiPIDRec = (TH1F*)hPtRS_TXiPIDRec->Clone("hRawPt_XiPIDRec_tight");
	hRawPt_TXiPIDRec->Add(hPtWS_TXiPIDRec,-1);
    TH1F *hRawPt_SXiPIDRec = (TH1F*) hPtRS_SXiPIDRec->Clone("hRawPt_XiPIDRec_stand");
	hRawPt_SXiPIDRec->Add(hPtWS_SXiPIDRec,-1);
    TH1F *hRawPt_LXiPIDRec = (TH1F*) hPtRS_LXiPIDRec->Clone("hRawPt_XiPIDRec_loose");
	hRawPt_LXiPIDRec->Add(hPtWS_LXiPIDRec,-1);
    TH1F *hRawPt_VLXiPIDRec = (TH1F*) hPtRS_VLXiPIDRec->Clone("hRawPt_XiPIDRec_vloose");
	hRawPt_VLXiPIDRec->Add(hPtWS_VLXiPIDRec,-1);

    TH1F *hRawPt_VTXiPIDVLRec = (TH1F*) hPtRS_VTXiPIDVLRec->Clone("hRawPt_XiPIDinvRec_vtight");
	hRawPt_VTXiPIDVLRec->Add(hPtWS_VTXiPIDVLRec,-1);
    TH1F *hRawPt_TXiPIDLRec = (TH1F*) hPtRS_TXiPIDLRec->Clone("hRawPt_XiPIDinvRec_tight");
	hRawPt_TXiPIDLRec->Add(hPtWS_TXiPIDLRec,-1);
    TH1F *hRawPt_SXiPIDSRec = (TH1F*) hPtRS_SXiPIDSRec->Clone("hRawPt_XiPIDinvRec_stand");
	hRawPt_SXiPIDSRec->Add(hPtWS_SXiPIDSRec,-1);
    TH1F *hRawPt_LXiPIDTRec = (TH1F*) hPtRS_LXiPIDTRec->Clone("hRawPt_XiPIDinvRec_loose");
	hRawPt_LXiPIDTRec->Add(hPtWS_LXiPIDTRec,-1);
    TH1F *hRawPt_VLXiPIDVTRec = (TH1F*) hPtRS_VLXiPIDVTRec->Clone("hRawPt_XiPIDinvRec_vloose");
	hRawPt_VLXiPIDVTRec->Add(hPtWS_VLXiPIDVTRec,-1);

    TH1F *hRawPt_SXiPID_V0 = (TH1F*) hPtRS_SXiPID_V0->Clone("hRawPt_XiPID_stand_V0");
	hRawPt_SXiPID_V0->Add(hPtWS_SXiPID_V0,-1);
    TH1F *hRawPt_SXiPID_Xi = (TH1F*) hPtRS_SXiPID_Xi->Clone("hRawPt_XiPID_stand_Xi");
	hRawPt_SXiPID_Xi->Add(hPtWS_SXiPID_Xi,-1);
    TH1F *hRawPt_SXiPID_DCAb = (TH1F*) hPtRS_SXiPID_DCAb->Clone("hRawPt_XiPID_stand_DCAb");
	hRawPt_SXiPID_DCAb->Add(hPtWS_SXiPID_DCAb,-1);

    TH1F *hRawPt_VTXiRec = (TH1F*)hPtRS_VTXiRec->Clone("hRawPt_XiRec_vtight"); hRawPt_VTXiRec->Add(hPtWS_VTXiRec,-1);
    TH1F *hRawPt_TXiRec  = (TH1F*)hPtRS_TXiRec->Clone("hRawPt_XiRec_tight");   hRawPt_TXiRec->Add(hPtWS_TXiRec,-1);
    TH1F *hRawPt_SXiRec  = (TH1F*)hPtRS_SXiRec->Clone("hRawPt_XiRec_stand");   hRawPt_SXiRec->Add(hPtWS_SXiRec,-1);
    TH1F *hRawPt_LXiRec  = (TH1F*)hPtRS_LXiRec->Clone("hRawPt_XiRec_loose");   hRawPt_LXiRec->Add(hPtWS_LXiRec,-1);
    TH1F *hRawPt_VLXiRec = (TH1F*)hPtRS_VLXiRec->Clone("hRawPt_XiRec_vloose"); hRawPt_VLXiRec->Add(hPtWS_VLXiRec,-1);

    TH1F *hRawPt_VLeRec = (TH1F*) hPtRS_VLeRec->Clone("hRawPt_eRec_vloose"); hRawPt_VLeRec->Add(hPtWS_VLeRec,-1);
    TH1F *hRawPt_LeRec  = (TH1F*) hPtRS_LeRec->Clone("hRawPt_eRec_loose");   hRawPt_LeRec->Add(hPtWS_LeRec,-1);
    TH1F *hRawPt_SeRec  = (TH1F*) hPtRS_SeRec->Clone("hRawPt_eRec_stand");   hRawPt_SeRec->Add(hPtWS_SeRec,-1);
    TH1F *hRawPt_TeRec  = (TH1F*) hPtRS_TeRec->Clone("hRawPt_eRec_tight");   hRawPt_TeRec->Add(hPtWS_TeRec,-1);
    TH1F *hRawPt_VTeRec = (TH1F*) hPtRS_VTeRec->Clone("hRawPt_eRec_vtight"); hRawPt_VTeRec->Add(hPtWS_VTeRec,-1);

    TH1F *hRawPt_LOA = (TH1F*) hPtRS_LOA->Clone("hRawPt_OA_loose"); hRawPt_LOA->Add(hPtWS_LOA,-1);
    TH1F *hRawPt_TOA = (TH1F*) hPtRS_TOA->Clone("hRawPt_OA_tight"); hRawPt_TOA->Add(hPtWS_TOA,-1);
    TH1F *hRawPt_SOA = (TH1F*) hPtRS_SOA->Clone("hRawPt_OA_stand"); hRawPt_SOA->Add(hPtWS_SOA,-1);
    TH1F *hRawPt_LeIM = (TH1F*) hPtRS_LeIM->Clone("hRawPt_IM_loose"); hRawPt_LeIM->Add(hPtWS_LeIM,-1);
    TH1F *hRawPt_TeIM = (TH1F*) hPtRS_TeIM->Clone("hRawPt_IM_tight"); hRawPt_TeIM->Add(hPtWS_TeIM,-1);
    TH1F *hRawPt_SeIM = (TH1F*) hPtRS_SeIM->Clone("hRawPt_IM_stand"); hRawPt_SeIM->Add(hPtWS_SeIM,-1);

    TH1F *hRawPt_Bayes_stand2 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand2");
    TH1F *hRawPt_Bayes_stand3 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand3");
    TH1F *hRawPt_Bayes_stand4 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand4");
    TH1F *hRawPt_Bayes_stand5 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand5");
    TH1F *hRawPt_Bayes_stand6 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand6");
    TH1F *hRawPt_Bayes_stand7 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand7");
    TH1F *hRawPt_Svd_stand3 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Svd_stand3");
    TH1F *hRawPt_Svd_stand4 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Svd_stand4");
    TH1F *hRawPt_Svd_stand5 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Svd_stand5");

    TH1F *hRawPt_BottomBaryon_default = (TH1F*) hRawPt_SeIM->Clone("hRawPt_BottomBaryon_default");
    TH1F *hRawPt_BottomBaryon_variation1 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_BottomBaryon_variation1");
    TH1F *hRawPt_BottomBaryon_variation2 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_BottomBaryon_variation2");
    TH1F *hRawPt_Weighting_default = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Weighting_default");
    TH1F *hRawPt_Weighting_variation1 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Weighting_variation1");
    TH1F *hRawPt_Weighting_variation2 = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Weighting_variation2");

	//*****************************************************

    TH1F* NumOfEvtperRun = (TH1F*) hist->FindObject("NumOfEvtperRun")->Clone("NumOfEvtperRun");
    TH1F* NumOfe         = (TH1F*) hist->FindObject("NumOfe")->Clone("NumOfe");
    TH1F* NumOfXi        = (TH1F*) hist->FindObject("NumOfXi")->Clone("NumOfXi");

    TH1D *hDS = (TH1D*) hist->FindObject("DSElectronPair")->Clone("hDS");
    TH1D* hSS = (TH1D*) hist->FindObject("SSElectronPair")->Clone("hSS");

    TH1F* NumOfXiperEvt = (TH1F*)NumOfXi->Clone("NumOfXiperEvt"); NumOfXiperEvt->Divide(NumOfEvtperRun);
    TH1F* NumOfeperEvt  = (TH1F*)NumOfe->Clone("NumOfeperEvt");   NumOfeperEvt->Divide(NumOfEvtperRun);

    TH1F* hpreff_eRec_loose = (TH1F*) hpre_eRec_loose_nu->Clone("hpreff_eRec_loose");
	hpreff_eRec_loose->Divide(hpreff_eRec_loose,hpre_eRec_loose_de,1,1,"b");
    TH1F* hpreff_eRec_stand = (TH1F*) hpre_eRec_stand_nu->Clone("hpreff_eRec_stand");
	hpreff_eRec_stand->Divide(hpreff_eRec_stand,hpre_eRec_stand_de,1,1,"b");
    TH1F* hpreff_eRec_tight = (TH1F*) hpre_eRec_tight_nu->Clone("hpreff_eRec_tight");
	hpreff_eRec_tight->Divide(hpreff_eRec_tight,hpre_eRec_tight_de,1,1,"b");

    TH1F* hpreff_XiRec_loose = (TH1F*) hpre_XiRec_loose_nu->Clone("hpreff_XiRec_loose");
	hpreff_XiRec_loose->Divide(hpreff_XiRec_loose,hpre_XiRec_loose_de,1,1,"b");
    TH1F* hpreff_XiRec_stand = (TH1F*) hpre_XiRec_stand_nu->Clone("hpreff_XiRec_stand");
	hpreff_XiRec_stand->Divide(hpreff_XiRec_stand,hpre_XiRec_stand_de,1,1,"b");
    TH1F* hpreff_XiRec_tight = (TH1F*) hpre_XiRec_tight_nu->Clone("hpreff_XiRec_tight");
	hpreff_XiRec_tight->Divide(hpreff_XiRec_tight,hpre_XiRec_tight_de,1,1,"b");

    TH1F* hpreff_ePID_loose = (TH1F*) hpre_ePID_loose_nu->Clone("hpreff_ePID_loose");
	hpreff_ePID_loose->Divide(hpreff_ePID_loose,hpre_ePID_loose_de,1,1,"b");
    TH1F* hpreff_ePID_stand = (TH1F*) hpre_ePID_stand_nu->Clone("hpreff_ePID_stand");
	hpreff_ePID_stand->Divide(hpreff_ePID_stand,hpre_ePID_stand_de,1,1,"b");
    TH1F* hpreff_ePID_tight = (TH1F*) hpre_ePID_tight_nu->Clone("hpreff_ePID_tight");
	hpreff_ePID_tight->Divide(hpreff_ePID_tight,hpre_ePID_tight_de,1,1,"b");

    TH1F* hpreff_XiPID_loose = (TH1F*) hpre_XiPID_loose_nu->Clone("hpreff_XiPID_loose");
	hpreff_XiPID_loose->Divide(hpreff_XiPID_loose,hpre_XiPID_loose_de,1,1,"b");
    TH1F* hpreff_XiPID_stand = (TH1F*) hpre_XiPID_stand_nu->Clone("hpreff_XiPID_stand");
	hpreff_XiPID_stand->Divide(hpreff_XiPID_stand,hpre_XiPID_stand_de,1,1,"b");
    TH1F* hpreff_XiPID_tight = (TH1F*) hpre_XiPID_tight_nu->Clone("hpreff_XiPID_tight");
	hpreff_XiPID_tight->Divide(hpreff_XiPID_tight,hpre_XiPID_tight_de,1,1,"b");

    TH1F* hpreff_OA_loose = (TH1F*) hpre_OA_loose_nu->Clone("hpreff_OA_loose");
	hpreff_OA_loose->Divide(hpreff_OA_loose,hpre_OA_loose_de,1,1,"b");
    TH1F* hpreff_OA_stand = (TH1F*) hpre_OA_stand_nu->Clone("hpreff_OA_stand");
	hpreff_OA_stand->Divide(hpreff_OA_stand,hpre_OA_stand_de,1,1,"b");
    TH1F* hpreff_OA_tight = (TH1F*) hpre_OA_tight_nu->Clone("hpreff_OA_tight");
	hpreff_OA_tight->Divide(hpreff_OA_tight,hpre_OA_tight_de,1,1,"b");

    TH1F* hpreff_IM_loose = (TH1F*) hpre_IM_loose_nu->Clone("hpreff_IM_loose");
	hpreff_IM_loose->Divide(hpreff_IM_loose,hpre_IM_loose_de,1,1,"b");
    TH1F* hpreff_IM_stand = (TH1F*) hpre_IM_stand_nu->Clone("hpreff_IM_stand");
	hpreff_IM_stand->Divide(hpreff_IM_stand,hpre_IM_stand_de,1,1,"b");
    TH1F* hpreff_IM_tight = (TH1F*) hpre_IM_tight_nu->Clone("hpreff_IM_tight");
	hpreff_IM_tight->Divide(hpreff_IM_tight,hpre_IM_tight_de,1,1,"b");

    hpre_Bayes_stand2_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand2_nu");
	hpre_Bayes_stand2_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand2_de");
    hpre_Bayes_stand3_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand3_nu");
	hpre_Bayes_stand3_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand3_de");
    hpre_Bayes_stand4_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand4_nu");
	hpre_Bayes_stand4_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand4_de");
    hpre_Bayes_stand5_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand5_nu");
	hpre_Bayes_stand5_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand5_de");
    hpre_Bayes_stand6_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand6_nu");
	hpre_Bayes_stand6_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand6_de");
    hpre_Bayes_stand7_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand7_nu");
	hpre_Bayes_stand7_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand7_de");

    hpre_Svd_stand3_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand3_nu");
	hpre_Svd_stand3_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand3_de");
    hpre_Svd_stand4_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand4_nu");
	hpre_Svd_stand4_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand4_de");
    hpre_Svd_stand5_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand5_nu");
	hpre_Svd_stand5_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand5_de");

    TH1F* hpre_Bayes_stand_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand_nu");
	TH1F* hpre_Bayes_stand_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand_de");

    TH1F* hpre_BottomBaryon_default_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_BottomBaryon_default_nu");
    TH1F* hpre_BottomBaryon_default_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_BottomBaryon_default_de");
    TH1F* hpre_BottomBaryon_variation1_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_BottomBaryon_variation1_nu");
    TH1F* hpre_BottomBaryon_variation1_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_BottomBaryon_variation1_de");
    TH1F* hpre_BottomBaryon_variation2_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_BottomBaryon_variation2_nu");
    TH1F* hpre_BottomBaryon_variation2_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_BottomBaryon_variation2_de");
    TH1F* hpre_Weighting_default_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Weighting_default_nu");
    TH1F* hpre_Weighting_default_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Weighting_default_de");
    TH1F* hpre_Weighting_variation1_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Weighting_variation1_nu");
    TH1F* hpre_Weighting_variation1_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Weighting_variation1_de");
    TH1F* hpre_Weighting_variation2_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Weighting_variation2_nu");
    TH1F* hpre_Weighting_variation2_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Weighting_variation2_de");

    TH1F* hRawPt_Svd_stand = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand");
    TH2F* hRPM_Unfold_stand = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand");
    TH2F* hRPM_Unfold_stand_pt2 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand_pt2");
    TH1F* hMCRecoLevXic0_Unfold_stand = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand");
    TH1F* hMCRecoLevPair_Unfold_stand = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand");
    TH1F* hMCRecoLevXic0_Unfold_stand_pt2 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand_pt2");
    TH1F* hMCRecoLevPair_Unfold_stand_pt2 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand_pt2");
    TH1F* hXic0_Svd_stand = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hXic0_Bayes_stand");
    TH1D* hEventNumbers = (TH1D*) hist->FindObject("hEventNumbers")->Clone("hEventNumbers");

	//kimc: new normalization histograms: trigger vs. multiplicity
	if (!IsMC)
	{
		TH2D* hNorm_multV0  = (TH2D*)hist->FindObject("hNorm_multV0")->Clone();
		TH2D* hNorm_multSPD = (TH2D*)hist->FindObject("hNorm_multSPD")->Clone();
	}
	#endif

 cout << test << endl;
	return;
}//MakeAnalysisHistogram

#endif //XI0CANAMAKEROOT
