#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TTree.h>
#include <TEfficiency.h>

using namespace std;

void XiCutDistribution(TObject *hist);  // look Xi cut value distribution from prompt an feeddown Xic0
void eXiPairTree(TFile* input, Bool_t IsMC); //generate main histogram which used for systematic study
void MakeAnalysisHistogram(TFile* input, Bool_t IsMC); //choose which one do you want to check
void nSigmaPlot(TObject* hist); //nSigmaTOF and TPC plot
void XiMassvsPt(TObject* hist); //Xi mass distribution for various pT bin

TH2F* MakeTH2(const Char_t *name, Int_t NumOfBin, Double_t* binning);
TH1F* MakeTH1(const Char_t *name, Int_t NumOfBin, Double_t* binning);

void MakeROOTResultXic0()
{
  Char_t *data;
  TFile *input;
  TFile *GridInput;
  Bool_t IsMC;
  Bool_t IsTrainOutput;
  TFile* TrainInput;

  cout << "IS MC ? (1 or 0) ";
  cin >> IsMC;
  cout << "IS Train Output ? (1 or 0) ";
  cin >> IsTrainOutput;
  cout << endl;
  if(IsMC) input = TFile::Open("/Users/jinjoo/Desktop/Xic0/data/Download/MC0815/AnalysisResultsMC.root");
  if(!IsMC) input = TFile::Open("/Users/jinjoo/Desktop/Xic0/data/Download/train_data/AnalysisResults_data.root");
  if(IsTrainOutput) TrainInput = (TFile*) input->Get("PWG3_D2H_Xic02eXipp13TeV");
  if(!IsTrainOutput) GridInput = (TFile*) input;

   //MakeAnalysisHistogram(inputdata, IsMC);
  if(IsTrainOutput) MakeAnalysisHistogram(TrainInput, IsMC);
  if(!IsTrainOutput) MakeAnalysisHistogram(GridInput, IsMC);
}

void MakeAnalysisHistogram(TFile* input, Bool_t IsMC){
  TObject* hist = input->Get("histogram");
  TFile* file;

  if(IsMC){
    cout << "---------------------MC reconstructed histogram---------------------" << endl;
    file = new TFile("MCWeightedAnalysisHistogram.root","recreate");
  }
  if(!IsMC){
    cout << "---------------------Data reconstructed histogram---------------------" << endl;
    file = new TFile("DataAnalysisHistogram.root","recreate");
  }
  //nSigmaPlot(hist);                  //choose
  //XiMassvsPt(hist);                  //choose
  //XiCutDistribution(hist)            //choose
  eXiPairTree(input,IsMC);             //choose       ---------------- difne the weighting function in the fisrt line in eXiPairTree
  cout << "---------------------Reconstructed histogram is done---------------------" << endl;

  file->Write();
  file->Close();

  return;
}

void XiCutDistribution(TObject* hist){
  TH1D* hC = (TH1D*) hist->FindObject("C_flag"); hC = (TH1D*) hC->Clone("C_flag");
  TH1D* hB = (TH1D*) hist->FindObject("B_flag"); hB = (TH1D*) hB->Clone("B_flag");
  TH1D* hBc = (TH1D*) hist->FindObject("Bcut_flag"); hBc = (TH1D*) hBc->Clone("Bcut_flag");
  TH1D* hCc = (TH1D*) hist->FindObject("Ccut_flag"); hCc = (TH1D*) hCc->Clone("Ccut_flag");
  TH1D* hec = (TH1D*) hist->FindObject("e_c_flag"); hec = (TH1D*) hec->Clone("hec");
  TH1D* heb = (TH1D*) hist->FindObject("e_b_flag"); heb = (TH1D*) heb->Clone("heb");

  TH1D* b1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertexMin_b");   b1 = (TH1D*) b1->Clone("hDCAV0PrToPrimVertexMin_b");
  TH1D* b2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertexMin_b");   b2 = (TH1D*) b2->Clone("hDCAV0PiToPrimVertexMin_b");
  TH1D* b3 = (TH1D*) hist->FindObject("hDCABachToPrimVertexMin_b");   b3 = (TH1D*) b3->Clone("hDCABachToPrimVertexMin_b");
  TH1D* b4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertexMin_b");   b4 = (TH1D*) b4->Clone("hDCAV0ToPrimVertexMin_b");
  TH1D* b5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXiMin_b");   b5 = (TH1D*) b5->Clone("hV0CosineOfPoiningAngleXiMin_b");
  TH1D* b6 = (TH1D*) hist->FindObject("hCascDecayLengthMin_b");   b6 = (TH1D*) b6->Clone("hCascDecayLengthMin_b");
  TH1D* b7 = (TH1D*) hist->FindObject("hDecayLengthV0_b");   b7 = (TH1D*) b7->Clone("hDecayLengthV0_b");

  TH1D* c1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertexMin_c"); c1 = (TH1D*) c1->Clone("hDCAV0PrToPrimVertexMin_c");
  TH1D* c2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertexMin_c"); c2 = (TH1D*) c2->Clone("hDCAV0PiToPrimVertexMin_c");
  TH1D* c3 = (TH1D*) hist->FindObject("hDCABachToPrimVertexMin_c"); c3 = (TH1D*) c3->Clone("hDCABachToPrimVertexMin_c");
  TH1D* c4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertexMin_c"); c4 = (TH1D*) c4->Clone("hDCAV0ToPrimVertexMin_c");
  TH1D* c5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXiMin_c"); c5 = (TH1D*) c5->Clone("hV0CosineOfPoiningAngleXiMin_c");
  TH1D* c6 = (TH1D*) hist->FindObject("hCascDecayLengthMin_c"); c6 = (TH1D*) c6->Clone("hCascDecayLengthMin_c");
  TH1D* c7 = (TH1D*) hist->FindObject("hDecayLengthV0_c"); c7 = (TH1D*) c7->Clone("hDecayLengthV0_c");

  TH1D* d1 = (TH1D*) hist->FindObject("hDCAV0PrToPrimVertexMin"); d1 = (TH1D*) d1->Clone("hDCAV0PrToPrimVertexMin");
  TH1D* d2 = (TH1D*) hist->FindObject("hDCAV0PiToPrimVertexMin"); d2 = (TH1D*) d2->Clone("hDCAV0PiToPrimVertexMin");
  TH1D* d3 = (TH1D*) hist->FindObject("hDCABachToPrimVertexMin"); d3 = (TH1D*) d3->Clone("hDCABachToPrimVertexMin");
  TH1D* d4 = (TH1D*) hist->FindObject("hDCAV0ToPrimVertexMin"); d4 = (TH1D*) d4->Clone("hDCAV0ToPrimVertexMin");
  TH1D* d5 = (TH1D*) hist->FindObject("hV0CosineOfPoiningAngleXiMin"); d5 = (TH1D*) d5->Clone("hV0CosineOfPoiningAngleXiMin");
  TH1D* d6 = (TH1D*) hist->FindObject("hCascDecayLengthMin"); d6 = (TH1D*) d6->Clone("hCascDecayLengthMin");
  TH1D* d7 = (TH1D*) hist->FindObject("hDecayLengthV0"); d7 = (TH1D*) d7->Clone("hDecayLengthV0");

  return;
}

void nSigmaPlot(TObject* hist){
  TH2D* hTPC = (TH2D*) hist->FindObject("nSigmaTPCvsPt"); hTPC = (TH2D*) hTPC->Clone("nSigmaTPCvsPt");
  TH2D* hTOF = (TH2D*) hist->FindObject("nSigmaTOFvsPt"); hTOF = (TH2D*) hTOF->Clone("nSigmaTOFvsPt");

  TH1D * hTPC05_1 = hTPC->ProjectionY("hTPC05_1",1,2);
  TH1D * hTPC1_2 = hTPC->ProjectionY("hTPC1_2",3,4);
  TH1D * hTPC2_3 = hTPC->ProjectionY("hTPC2_3",5,6);
  TH1D * hTPC3_4 = hTPC->ProjectionY("hTPC3_4",7,8);
  TH1D * hTPC4_5 = hTPC->ProjectionY("hTPC4_5",9,10);

  return;
}

void XiMassvsPt(TObject* hist){
  TH2D* XiMassvsPt = (TH2D*) hist->FindObject("hXimassvsPt"); XiMassvsPt = (TH2D*) XiMassvsPt->Clone("massvspt");

  TH1D *hXi0_1 = XiMassvsPt->ProjectionX("hXiMass0_1",0,1);
  TH1D *hXi1_2 = XiMassvsPt->ProjectionX("hXiMass1_2",1,2);
  TH1D *hXi2_3 = XiMassvsPt->ProjectionX("hXiMass2_3",2,3);
  TH1D *hXi3_4 = XiMassvsPt->ProjectionX("hXiMass3_4",3,4);
  TH1D *hXi4_5 = XiMassvsPt->ProjectionX("hXiMass4_5",4,5);
  TH1D *hXi5_6 = XiMassvsPt->ProjectionX("hXiMass5_6",5,6);
  TH1D *hXi6_7 = XiMassvsPt->ProjectionX("hXiMass6_7",6,7);
  TH1D *hXi7_8 = XiMassvsPt->ProjectionX("hXiMass7_8",7,8);
  TH1D *hXi8_9 = XiMassvsPt->ProjectionX("hXiMass8_9",8,9);
  TH1D *hXi9_10 = XiMassvsPt->ProjectionX("hXiMass9_10",9,10);

  return;
}

TH1F* MakeTH1(const Char_t *name, Int_t NumOfBin, Double_t* binning){
  TH1F* hist = new TH1F(name,"",NumOfBin,binning);
  hist->Sumw2();
  return hist;
}

TH2F* MakeTH2(const Char_t *name, Int_t NumOfBin, Double_t* binning){
  TH2F* hist = new TH2F(name,"",NumOfBin,binning,NumOfBin,binning);
  hist->Sumw2();
  return hist;
}

void eXiPairTree(TFile* input, Bool_t IsMC){
  TF1 * fWeightFit = new TF1("fWeightFit","expo",0,12);
  fWeightFit -> SetParameter(0,7.48584e-01);
  fWeightFit -> SetParameter(1,-1.47142e-01);

 // Double_t binning1[8] = {0,1.,2.,3.2,4.4,6.,8.,12.};
  //Double_t binning1[8] = {1.,2.,3.,4.,5.,6.,8.,12.};
  Double_t binning1[10] = {1,2.,3.,4.,5.,6.,8.,12.,16.,20};
  Double_t binning2[10] = {1,2.,3.,4.,5.,6.,8.,12.,16.,20};
  Double_t binning3[8] = {1.,2.,3.,4.,5.,6.,8.,12.};

  TTree *Pair = (TTree *) input->Get("eXiTree");
  Float_t pTe; Float_t echarge;
  Float_t pTv; Float_t vcharge;
  Float_t Massv; Float_t cosoa; Float_t In_Mass; Float_t Pt;
  Float_t nSigmaTOF; Float_t nSigmaTPC; Float_t TPCCluster; Float_t ITSCluster; Float_t TPCPIDCluster;
  Float_t CascDecayLength; Float_t V0DecayLength; Float_t DCABachToPrimVertex;
  Float_t V0CosineOfPoiningAngleXi; Float_t DCAV0ToPrimVertex; Float_t DCAPosToPrimVertex; Float_t DCANegToPrimVertex;
  Float_t e_minmass; Float_t e_minmass_ss; Float_t phi; Float_t erap; Float_t Xirap;
  Float_t pionTPCCluster; Float_t protonTPCCluster; Float_t b_pionTPCCluster;
  Float_t e_crossedratio; Float_t e_findable; Float_t pion_crossedratio; Float_t pion_findable;
  Float_t proton_crossedratio; Float_t proton_findable; Float_t bpion_crossedratio; Float_t bpion_findable;
  Float_t XiCosineOfPoiningAngle; Float_t pTpion; Float_t pTproton; Float_t pTbach;
  Float_t MassLambda; Float_t MassAntiLambda;

  Pair->SetBranchAddress("pTe",&pTe);
  Pair->SetBranchAddress("echarge",&echarge);
  Pair->SetBranchAddress("TOFnSigma",&nSigmaTOF);
  Pair->SetBranchAddress("TPCnSigma",&nSigmaTPC);
  //Pair->SetBranchAddress("TPC",&TPCCluster);
  Pair->SetBranchAddress("TPCPID",&TPCPIDCluster);
  Pair->SetBranchAddress("ITS",&ITSCluster);
  Pair->SetBranchAddress("e_crossedrows",&e_crossedratio);
  Pair->SetBranchAddress("e_findable",&e_findable);
  Pair->SetBranchAddress("phi",&phi);
  Pair->SetBranchAddress("erap",&erap);
  Pair->SetBranchAddress("e_minmass",&e_minmass);
  Pair->SetBranchAddress("e_minmass_ss",&e_minmass_ss);
  Pair->SetBranchAddress("pTv",&pTv);
  Pair->SetBranchAddress("vcharge",&vcharge);
  Pair->SetBranchAddress("Massv",&Massv);
  Pair->SetBranchAddress("MassLambda",&MassLambda);
  Pair->SetBranchAddress("MassAntiLambda",&MassAntiLambda);  //mod
  Pair->SetBranchAddress("Xirap",&Xirap);  //mod
  Pair->SetBranchAddress("V0DecayLength",&V0DecayLength);
  Pair->SetBranchAddress("CascDecayLength",&CascDecayLength);
  Pair->SetBranchAddress("DCABachToPrimVertex",&DCABachToPrimVertex);
  Pair->SetBranchAddress("DCAV0NegToPrimVertex",&DCANegToPrimVertex);
  Pair->SetBranchAddress("DCAV0PosToPrimVertex",&DCAPosToPrimVertex);
  Pair->SetBranchAddress("V0CosineOfPoiningAngleXi",&V0CosineOfPoiningAngleXi);
  Pair->SetBranchAddress("XiCosineOfPoiningAngle",&XiCosineOfPoiningAngle);  //modify
  Pair->SetBranchAddress("DCAV0ToPrimVertex",&DCAV0ToPrimVertex);   ///new
  Pair->SetBranchAddress("cosoa",&cosoa);  //mod
  Pair->SetBranchAddress("pion_crossedrows",&pion_crossedratio);
  Pair->SetBranchAddress("pion_findable",&pion_findable);
  Pair->SetBranchAddress("proton_crossedrows",&proton_crossedratio);
  Pair->SetBranchAddress("proton_findable",&proton_findable);
  Pair->SetBranchAddress("bpion_crossedratio",&bpion_crossedratio);
  Pair->SetBranchAddress("bpion_findable",&bpion_findable);
  Pair->SetBranchAddress("pTpion",&pTpion);
  Pair->SetBranchAddress("pTproton",&pTproton);
  Pair->SetBranchAddress("pTbach",&pTbach);
  Pair->SetBranchAddress("In_Mass",&In_Mass);
  Pair->SetBranchAddress("eXiPt",&Pt);

  TTree *MCTree = (TTree *) input->Get("MCTree");

  Float_t mcpTe; Float_t mcecharge;
  Float_t mcpTv; Float_t mcvcharge;
  Float_t mcpteXi; Float_t mcptXic0; Float_t mcc_flag; Float_t mcb_flag; Float_t mcXib; Float_t XibeXi; Float_t mcXibMass;

  MCTree->SetBranchAddress("mcpTe",&mcpTe);
  MCTree->SetBranchAddress("mcecharge",&mcecharge);
  MCTree->SetBranchAddress("mcpTv",&mcpTv);
  MCTree->SetBranchAddress("mcvcharge",&mcvcharge);
  MCTree->SetBranchAddress("mcpTeXi",&mcpteXi);
  MCTree->SetBranchAddress("mcpTXic0",&mcptXic0);
  MCTree->SetBranchAddress("c_flag",&mcc_flag);
  MCTree->SetBranchAddress("b_flag",&mcb_flag);
  MCTree->SetBranchAddress("mcpTXib",&mcXib);
  MCTree->SetBranchAddress("mceXipTb",&XibeXi);
  MCTree->SetBranchAddress("mcXibMass",&mcXibMass);

//-------------------------- raw yield --------------------------//
  TH1F *hPtRS_TePID = MakeTH1("hPtRS_TePID",9,binning1); TH1F *hPtWS_TePID = MakeTH1("hPtWS_TePID",9,binning1);
  TH1F *hPtRS_SePID = MakeTH1("hPtRS_SePID",9,binning1); TH1F *hPtWS_SePID = MakeTH1("hPtWS_SePID",9,binning1);
  TH1F *hPtRS_LePID = MakeTH1("hPtRS_LePID",9,binning1); TH1F *hPtWS_LePID = MakeTH1("hPtWS_LePID",9,binning1);
  TH1F *hPtRS_VTePID = MakeTH1("hPtRS_VTePID",9,binning1); TH1F *hPtWS_VTePID = MakeTH1("hPtWS_VTePID",9,binning1);
  TH1F *hPtRS_VLePID = MakeTH1("hPtRS_VLePID",9,binning1); TH1F *hPtWS_VLePID = MakeTH1("hPtWS_VLePID",9,binning1);

  TH1F *hPtRS_TXiPID = MakeTH1("hPtRS_TXiPID",9,binning1); TH1F *hPtWS_TXiPID = MakeTH1("hPtWS_TXiPID",9,binning1);
  TH1F *hPtRS_SXiPID = MakeTH1("hPtRS_SXiPID",9,binning1); TH1F *hPtWS_SXiPID = MakeTH1("hPtWS_SXiPID",9,binning1);
  TH1F *hPtRS_LXiPID = MakeTH1("hPtRS_LXiPID",9,binning1); TH1F *hPtWS_LXiPID = MakeTH1("hPtWS_LXiPID",9,binning1);
  TH1F *hPtRS_VTXiPID = MakeTH1("hPtRS_VTXiPID",9,binning1); TH1F *hPtWS_VTXiPID = MakeTH1("hPtWS_VTXiPID",9,binning1);
  TH1F *hPtRS_VLXiPID = MakeTH1("hPtRS_VLXiPID",9,binning1); TH1F *hPtWS_VLXiPID = MakeTH1("hPtWS_VLXiPID",9,binning1);

  TH1F *hPtRS_TXiPIDRec = MakeTH1("hPtRS_TXiPIDRec",9,binning1); TH1F *hPtWS_TXiPIDRec = MakeTH1("hPtWS_TXiPIDRec",9,binning1);
  TH1F *hPtRS_SXiPIDRec = MakeTH1("hPtRS_SXiPIDRec",9,binning1); TH1F *hPtWS_SXiPIDRec = MakeTH1("hPtWS_SXiPIDRec",9,binning1);
  TH1F *hPtRS_LXiPIDRec = MakeTH1("hPtRS_LXiPIDRec",9,binning1); TH1F *hPtWS_LXiPIDRec = MakeTH1("hPtWS_LXiPIDRec",9,binning1);
  TH1F *hPtRS_VTXiPIDRec = MakeTH1("hPtRS_VTXiPIDRec",9,binning1); TH1F *hPtWS_VTXiPIDRec = MakeTH1("hPtWS_VTXiPIDRec",9,binning1);
  TH1F *hPtRS_VLXiPIDRec = MakeTH1("hPtRS_VLXiPIDRec",9,binning1); TH1F *hPtWS_VLXiPIDRec = MakeTH1("hPtWS_VLXiPIDRec",9,binning1);
  TH1F *hPtRS_TXiPIDLRec = MakeTH1("hPtRS_TXiPIDLRec",9,binning1); TH1F *hPtWS_TXiPIDLRec = MakeTH1("hPtWS_TXiPIDLRec",9,binning1);
  TH1F *hPtRS_SXiPIDSRec = MakeTH1("hPtRS_SXiPIDSRec",9,binning1); TH1F *hPtWS_SXiPIDSRec = MakeTH1("hPtWS_SXiPIDSRec",9,binning1);
  TH1F *hPtRS_LXiPIDTRec = MakeTH1("hPtRS_LXiPIDTRec",9,binning1); TH1F *hPtWS_LXiPIDTRec = MakeTH1("hPtWS_LXiPIDTRec",9,binning1);
  TH1F *hPtRS_VTXiPIDVLRec = MakeTH1("hPtRS_VTXiPIDVLRec",9,binning1); TH1F *hPtWS_VTXiPIDVLRec = MakeTH1("hPtWS_VTXiPIDVLRec",9,binning1);
  TH1F *hPtRS_VLXiPIDVTRec = MakeTH1("hPtRS_VLXiPIDVTRec",9,binning1); TH1F *hPtWS_VLXiPIDVTRec = MakeTH1("hPtWS_VLXiPIDVTRec",9,binning1);

  TH1F *hPtRS_SXiPID_V0 = MakeTH1("hPtRS_SXiPID_V0",9,binning1); TH1F *hPtWS_SXiPID_V0 = MakeTH1("hPtWS_SXiPID_V0",9,binning1);
  TH1F *hPtRS_SXiPID_Xi = MakeTH1("hPtRS_SXiPID_Xi",9,binning1); TH1F *hPtWS_SXiPID_Xi = MakeTH1("hPtWS_SXiPID_Xi",9,binning1);
  TH1F *hPtRS_SXiPID_DCAb = MakeTH1("hPtRS_SXiPID_DCAb",9,binning1); TH1F *hPtWS_SXiPID_DCAb = MakeTH1("hPtWS_SXiPID_DCAb",9,binning1);

  TH1F *hPtRS_LeRec = MakeTH1("hPtRS_LeRec",9,binning1); TH1F *hPtWS_LeRec = MakeTH1("hPtWS_LeRec",9,binning1);
  TH1F *hPtRS_SeRec = MakeTH1("hPtRS_SeRec",9,binning1); TH1F *hPtWS_SeRec = MakeTH1("hPtWS_SeRec",9,binning1);
  TH1F *hPtRS_TeRec = MakeTH1("hPtRS_TeRec",9,binning1); TH1F *hPtWS_TeRec = MakeTH1("hPtWS_TeRec",9,binning1);
  TH1F *hPtRS_VTeRec = MakeTH1("hPtRS_VTeRec",9,binning1); TH1F *hPtWS_VTeRec = MakeTH1("hPtWS_VTeRec",9,binning1);
  TH1F *hPtRS_VLeRec = MakeTH1("hPtRS_VLeRec",9,binning1); TH1F *hPtWS_VLeRec = MakeTH1("hPtWS_VLeRec",9,binning1);

  TH1F *hPtRS_LXiRec = MakeTH1("hPtRS_LXiRec",9,binning1); TH1F *hPtWS_LXiRec = MakeTH1("hPtWS_LXiRec",9,binning1);
  TH1F *hPtRS_SXiRec = MakeTH1("hPtRS_SXiRec",9,binning1); TH1F *hPtWS_SXiRec = MakeTH1("hPtWS_SXiRec",9,binning1);
  TH1F *hPtRS_TXiRec = MakeTH1("hPtRS_TXiRec",9,binning1); TH1F *hPtWS_TXiRec = MakeTH1("hPtWS_TXiRec",9,binning1);
  TH1F *hPtRS_VTXiRec = MakeTH1("hPtRS_VTXiRec",9,binning1); TH1F *hPtWS_VTXiRec = MakeTH1("hPtWS_VTXiRec",9,binning1);
  TH1F *hPtRS_VLXiRec = MakeTH1("hPtRS_VLXiRec",9,binning1); TH1F *hPtWS_VLXiRec = MakeTH1("hPtWS_VLXiRec",9,binning1);

  TH1F *hPtRS_LeIM = MakeTH1("hPtRS_LeIM",9,binning1); TH1F *hPtWS_LeIM = MakeTH1("hPtWS_LeIM",9,binning1);
  TH1F *hPtRS_SeIM = MakeTH1("hPtRS_SeIM",9,binning1); TH1F *hPtWS_SeIM = MakeTH1("hPtWS_SeIM",9,binning1);
  TH1F *hPtRS_TeIM = MakeTH1("hPtRS_TeIM",9,binning1); TH1F *hPtWS_TeIM = MakeTH1("hPtWS_TeIM",9,binning1);

  TH1F *hPtRS_LOA = MakeTH1("hPtRS_LOA",9,binning1); TH1F *hPtWS_LOA = MakeTH1("hPtWS_LOA",9,binning1);
  TH1F *hPtRS_TOA = MakeTH1("hPtRS_TOA",9,binning1); TH1F *hPtWS_TOA = MakeTH1("hPtWS_TOA",9,binning1);
  TH1F *hPtRS_SOA = MakeTH1("hPtRS_SOA",9,binning1); TH1F *hPtWS_SOA = MakeTH1("hPtWS_SOA",9,binning1);

  TH2D *hMassPtRS = new TH2D("hMassPtRS","",10,1.3,3.3,60,0,12); hMassPtRS->Sumw2();
  TH2D *hMassPtWS = new TH2D("hMassPtWS","",10,1.3,3.3,60,0,12); hMassPtWS->Sumw2();

  //-------------------------- prefilter --------------------------//

  TH1F *hpre_eRec_loose_de = MakeTH1("hpre_eRec_loose_de",9,binning1);   TH1F *hpre_eRec_loose_nu = MakeTH1("hpre_eRec_loose_nu",9,binning1);
  TH1F *hpre_eRec_stand_de = MakeTH1("hpre_eRec_stand_de",9,binning1);   TH1F *hpre_eRec_stand_nu = MakeTH1("hpre_eRec_stand_nu",9,binning1);
  TH1F *hpre_eRec_tight_de = MakeTH1("hpre_eRec_tight_de",9,binning1);   TH1F *hpre_eRec_tight_nu = MakeTH1("hpre_eRec_tight_nu",9,binning1);
  TH1F *hpre_eRec_vloose_de = MakeTH1("hpre_eRec_vloose_de",9,binning1);   TH1F *hpre_eRec_vloose_nu = MakeTH1("hpre_eRec_vloose_nu",9,binning1);
  TH1F *hpre_eRec_vtight_de = MakeTH1("hpre_eRec_vtight_de",9,binning1);   TH1F *hpre_eRec_vtight_nu = MakeTH1("hpre_eRec_vtight_nu",9,binning1);

  TH1F *hpre_XiRec_loose_de = MakeTH1("hpre_XiRec_loose_de",9,binning1);   TH1F *hpre_XiRec_loose_nu = MakeTH1("hpre_XiRec_loose_nu",9,binning1);
  TH1F *hpre_XiRec_stand_de = MakeTH1("hpre_XiRec_stand_de",9,binning1);   TH1F *hpre_XiRec_stand_nu = MakeTH1("hpre_XiRec_stand_nu",9,binning1);
  TH1F *hpre_XiRec_tight_de = MakeTH1("hpre_XiRec_tight_de",9,binning1);   TH1F *hpre_XiRec_tight_nu = MakeTH1("hpre_XiRec_tight_nu",9,binning1);
  TH1F *hpre_XiRec_vloose_de = MakeTH1("hpre_XiRec_vloose_de",9,binning1);   TH1F *hpre_XiRec_vloose_nu = MakeTH1("hpre_XiRec_vloose_nu",9,binning1);
  TH1F *hpre_XiRec_vtight_de = MakeTH1("hpre_XiRec_vtight_de",9,binning1);   TH1F *hpre_XiRec_vtight_nu = MakeTH1("hpre_XiRec_vtight_nu",9,binning1);

  TH1F *hpre_ePID_loose_de = MakeTH1("hpre_ePID_loose_de",9,binning1);   TH1F *hpre_ePID_loose_nu = MakeTH1("hpre_ePID_loose_nu",9,binning1);
  TH1F *hpre_ePID_stand_de = MakeTH1("hpre_ePID_stand_de",9,binning1);   TH1F *hpre_ePID_stand_nu = MakeTH1("hpre_ePID_stand_nu",9,binning1);
  TH1F *hpre_ePID_tight_de = MakeTH1("hpre_ePID_tight_de",9,binning1);   TH1F *hpre_ePID_tight_nu = MakeTH1("hpre_ePID_tight_nu",9,binning1);
  TH1F *hpre_ePID_vloose_de = MakeTH1("hpre_ePID_vloose_de",9,binning1);   TH1F *hpre_ePID_vloose_nu = MakeTH1("hpre_ePID_vloose_nu",9,binning1);
  TH1F *hpre_ePID_vtight_de = MakeTH1("hpre_ePID_vtight_de",9,binning1);   TH1F *hpre_ePID_vtight_nu = MakeTH1("hpre_ePID_vtight_nu",9,binning1);

  TH1F *hpre_XiPID_loose_de = MakeTH1("hpre_XiPID_loose_de",9,binning1);   TH1F *hpre_XiPID_loose_nu = MakeTH1("hpre_XiPID_loose_nu",9,binning1);
  TH1F *hpre_XiPID_stand_de = MakeTH1("hpre_XiPID_stand_de",9,binning1);   TH1F *hpre_XiPID_stand_nu = MakeTH1("hpre_XiPID_stand_nu",9,binning1);
  TH1F *hpre_XiPID_tight_de = MakeTH1("hpre_XiPID_tight_de",9,binning1);   TH1F *hpre_XiPID_tight_nu = MakeTH1("hpre_XiPID_tight_nu",9,binning1);
  TH1F *hpre_XiPID_vloose_de = MakeTH1("hpre_XiPID_vloose_de",9,binning1);   TH1F *hpre_XiPID_vloose_nu = MakeTH1("hpre_XiPID_vloose_nu",9,binning1);
  TH1F *hpre_XiPID_vtight_de = MakeTH1("hpre_XiPID_vtight_de",9,binning1);   TH1F *hpre_XiPID_vtight_nu = MakeTH1("hpre_XiPID_vtight_nu",9,binning1);

  TH1F *hpre_XiPIDRec_loose_de = MakeTH1("hpre_XiPIDRec_loose_de",9,binning1);   TH1F *hpre_XiPIDRec_loose_nu = MakeTH1("hpre_XiPIDRec_loose_nu",9,binning1);
  TH1F *hpre_XiPIDRec_stand_de = MakeTH1("hpre_XiPIDRec_stand_de",9,binning1);   TH1F *hpre_XiPIDRec_stand_nu = MakeTH1("hpre_XiPIDRec_stand_nu",9,binning1);
  TH1F *hpre_XiPIDRec_tight_de = MakeTH1("hpre_XiPIDRec_tight_de",9,binning1);   TH1F *hpre_XiPIDRec_tight_nu = MakeTH1("hpre_XiPIDRec_tight_nu",9,binning1);
  TH1F *hpre_XiPIDRec_vloose_de = MakeTH1("hpre_XiPIDRec_vloose_de",9,binning1);   TH1F *hpre_XiPIDRec_vloose_nu = MakeTH1("hpre_XiPIDRec_vloose_nu",9,binning1);
  TH1F *hpre_XiPIDRec_vtight_de = MakeTH1("hpre_XiPIDRec_vtight_de",9,binning1);   TH1F *hpre_XiPIDRec_vtight_nu = MakeTH1("hpre_XiPIDRec_vtight_nu",9,binning1);

  TH1F *hpre_XiPIDinvRec_loose_de = MakeTH1("hpre_XiPIDinvRec_loose_de",9,binning1);   TH1F *hpre_XiPIDinvRec_loose_nu = MakeTH1("hpre_XiPIDinvRec_loose_nu",9,binning1);
  TH1F *hpre_XiPIDinvRec_stand_de = MakeTH1("hpre_XiPIDinvRec_stand_de",9,binning1);   TH1F *hpre_XiPIDinvRec_stand_nu = MakeTH1("hpre_XiPIDinvRec_stand_nu",9,binning1);
  TH1F *hpre_XiPIDinvRec_tight_de = MakeTH1("hpre_XiPIDinvRec_tight_de",9,binning1);   TH1F *hpre_XiPIDinvRec_tight_nu = MakeTH1("hpre_XiPIDinvRec_tight_nu",9,binning1);
  TH1F *hpre_XiPIDinvRec_vloose_de = MakeTH1("hpre_XiPIDinvRec_vloose_de",9,binning1);   TH1F *hpre_XiPIDinvRec_vloose_nu = MakeTH1("hpre_XiPIDinvRec_vloose_nu",9,binning1);
  TH1F *hpre_XiPIDinvRec_vtight_de = MakeTH1("hpre_XiPIDinvRec_vtight_de",9,binning1);   TH1F *hpre_XiPIDinvRec_vtight_nu = MakeTH1("hpre_XiPIDinvRec_vtight_nu",9,binning1);

  TH1F *hpre_OA_loose_de = MakeTH1("hpre_OA_loose_de",9,binning1);   TH1F *hpre_OA_loose_nu = MakeTH1("hpre_OA_loose_nu",9,binning1);
  TH1F *hpre_OA_stand_de = MakeTH1("hpre_OA_stand_de",9,binning1);   TH1F *hpre_OA_stand_nu = MakeTH1("hpre_OA_stand_nu",9,binning1);
  TH1F *hpre_OA_tight_de = MakeTH1("hpre_OA_tight_de",9,binning1);   TH1F *hpre_OA_tight_nu = MakeTH1("hpre_OA_tight_nu",9,binning1);

  TH1F *hpre_IM_loose_de = MakeTH1("hpre_IM_loose_de",9,binning1);   TH1F *hpre_IM_loose_nu = MakeTH1("hpre_IM_loose_nu",9,binning1);
  TH1F *hpre_IM_stand_de = MakeTH1("hpre_IM_stand_de",9,binning1);   TH1F *hpre_IM_stand_nu = MakeTH1("hpre_IM_stand_nu",9,binning1);
  TH1F *hpre_IM_tight_de = MakeTH1("hpre_IM_tight_de",9,binning1);   TH1F *hpre_IM_tight_nu = MakeTH1("hpre_IM_tight_nu",9,binning1);

  TH1F *hpre_Bayes_stand2_de = MakeTH1("hpre_Bayes_stand2_de",9,binning1);   TH1F *hpre_Bayes_stand2_nu = MakeTH1("hpre_Bayes_stand2_nu",9,binning1);
  TH1F *hpre_Bayes_stand3_de = MakeTH1("hpre_Bayes_stand3_de",9,binning1);   TH1F *hpre_Bayes_stand3_nu = MakeTH1("hpre_Bayes_stand3_nu",9,binning1);
  TH1F *hpre_Bayes_stand4_de = MakeTH1("hpre_Bayes_stand4_de",9,binning1);   TH1F *hpre_Bayes_stand4_nu = MakeTH1("hpre_Bayes_stand4_nu",9,binning1);
  TH1F *hpre_Bayes_stand5_de = MakeTH1("hpre_Bayes_stand5_de",9,binning1);   TH1F *hpre_Bayes_stand5_nu = MakeTH1("hpre_Bayes_stand5_nu",9,binning1);
  TH1F *hpre_Bayes_stand6_de = MakeTH1("hpre_Bayes_stand6_de",9,binning1);   TH1F *hpre_Bayes_stand6_nu = MakeTH1("hpre_Bayes_stand6_nu",9,binning1);
  TH1F *hpre_Bayes_stand7_de = MakeTH1("hpre_Bayes_stand7_de",9,binning1);   TH1F *hpre_Bayes_stand7_nu = MakeTH1("hpre_Bayes_stand7_nu",9,binning1);

  TH1F *hpre_Svd_stand3_de = MakeTH1("hpre_Svd_stand3_de",9,binning1);   TH1F *hpre_Svd_stand3_nu = MakeTH1("hpre_Svd_stand3_nu",9,binning1);
  TH1F *hpre_Svd_stand4_de = MakeTH1("hpre_Svd_stand4_de",9,binning1);   TH1F *hpre_Svd_stand4_nu = MakeTH1("hpre_Svd_stand4_nu",9,binning1);
  TH1F *hpre_Svd_stand5_de = MakeTH1("hpre_Svd_stand5_de",9,binning1);   TH1F *hpre_Svd_stand5_nu = MakeTH1("hpre_Svd_stand5_nu",9,binning1);

  //-------------------------- MC historgram using unfolding and efficiency calculation --------------------------//

  TH2F* hRPM_eRec_loose_un = MakeTH2("hRPM_eRec_loose_un",9,binning2);  //response matrix using unfolding
  TH2F* hRPM_eRec_loose = MakeTH2("hRPM_eRec_loose",9,binning2);  //response matrix using refolding
  TH1F* hMCRecoLevXic0_eRec_loose = MakeTH1("hMCRecoLevXic0_eRec_loose",9,binning2);  //Xic0 in generation level
  TH1F* hMCRecoLevPair_eRec_loose = MakeTH1("hMCRecoLevPair_eRec_loose",9,binning2);  //eXi pair decay from Xic0 in generation level
  TH2F* hRPM_eRec_stand_un = MakeTH2("hRPM_eRec_stand_un",9,binning2);
  TH2F* hRPM_eRec_stand = MakeTH2("hRPM_eRec_stand",9,binning2);
  TH1F* hMCRecoLevXic0_eRec_stand = MakeTH1("hMCRecoLevXic0_eRec_stand",9,binning2);
  TH1F* hMCRecoLevPair_eRec_stand = MakeTH1("hMCRecoLevPair_eRec_stand",9,binning2);
  TH2F* hRPM_eRec_tight_un = MakeTH2("hRPM_eRec_tight_un",9,binning2);
  TH2F* hRPM_eRec_tight = MakeTH2("hRPM_eRec_tight",9,binning2);
  TH1F* hMCRecoLevXic0_eRec_tight = MakeTH1("hMCRecoLevXic0_eRec_tight",9,binning2);
  TH1F* hMCRecoLevPair_eRec_tight = MakeTH1("hMCRecoLevPair_eRec_tight",9,binning2);
  TH2F* hRPM_eRec_vloose_un = MakeTH2("hRPM_eRec_vloose_un",9,binning2);
  TH2F* hRPM_eRec_vloose = MakeTH2("hRPM_eRec_vloose",9,binning2);
  TH1F* hMCRecoLevXic0_eRec_vloose = MakeTH1("hMCRecoLevXic0_eRec_vloose",9,binning2);
  TH1F* hMCRecoLevPair_eRec_vloose = MakeTH1("hMCRecoLevPair_eRec_vloose",9,binning2);
  TH2F* hRPM_eRec_vtight_un = MakeTH2("hRPM_eRec_vtight_un",9,binning2);
  TH2F* hRPM_eRec_vtight = MakeTH2("hRPM_eRec_vtight",9,binning2);
  TH1F* hMCRecoLevXic0_eRec_vtight = MakeTH1("hMCRecoLevXic0_eRec_vtight",9,binning2);
  TH1F* hMCRecoLevPair_eRec_vtight = MakeTH1("hMCRecoLevPair_eRec_vtight",9,binning2);

  TH2F* hRPM_ePID_loose_un = MakeTH2("hRPM_ePID_loose_un",9,binning2);
  TH2F* hRPM_ePID_loose = MakeTH2("hRPM_ePID_loose",9,binning2);
  TH1F* hMCRecoLevXic0_ePID_loose = MakeTH1("hMCRecoLevXic0_ePID_loose",9,binning2);
  TH1F* hMCRecoLevPair_ePID_loose = MakeTH1("hMCRecoLevPair_ePID_loose",9,binning2);
  TH2F* hRPM_ePID_stand_un = MakeTH2("hRPM_ePID_stand_un",9,binning2);
  TH2F* hRPM_ePID_stand = MakeTH2("hRPM_ePID_stand",9,binning2);
  TH1F* hMCRecoLevXic0_ePID_stand = MakeTH1("hMCRecoLevXic0_ePID_stand",9,binning2);
  TH1F* hMCRecoLevPair_ePID_stand = MakeTH1("hMCRecoLevPair_ePID_stand",9,binning2);
  TH2F* hRPM_ePID_tight_un = MakeTH2("hRPM_ePID_tight_un",9,binning2);
  TH2F* hRPM_ePID_tight = MakeTH2("hRPM_ePID_tight",9,binning2);
  TH1F* hMCRecoLevXic0_ePID_tight = MakeTH1("hMCRecoLevXic0_ePID_tight",9,binning2);
  TH1F* hMCRecoLevPair_ePID_tight = MakeTH1("hMCRecoLevPair_ePID_tight",9,binning2);
  TH2F* hRPM_ePID_vloose_un = MakeTH2("hRPM_ePID_vloose_un",9,binning2);
  TH2F* hRPM_ePID_vloose = MakeTH2("hRPM_ePID_vloose",9,binning2);
  TH1F* hMCRecoLevXic0_ePID_vloose = MakeTH1("hMCRecoLevXic0_ePID_vloose",9,binning2);
  TH1F* hMCRecoLevPair_ePID_vloose = MakeTH1("hMCRecoLevPair_ePID_vloose",9,binning2);
  TH2F* hRPM_ePID_vtight_un = MakeTH2("hRPM_ePID_vtight_un",9,binning2);
  TH2F* hRPM_ePID_vtight = MakeTH2("hRPM_ePID_vtight",9,binning2);
  TH1F* hMCRecoLevXic0_ePID_vtight = MakeTH1("hMCRecoLevXic0_ePID_vtight",9,binning2);
  TH1F* hMCRecoLevPair_ePID_vtight = MakeTH1("hMCRecoLevPair_ePID_vtight",9,binning2);

  TH2F* hRPM_XiRec_loose_un = MakeTH2("hRPM_XiRec_loose_un",9,binning2);
  TH2F* hRPM_XiRec_loose = MakeTH2("hRPM_XiRec_loose",9,binning2);
  TH1F* hMCRecoLevXic0_XiRec_loose = MakeTH1("hMCRecoLevXic0_XiRec_loose",9,binning2);
  TH1F* hMCRecoLevPair_XiRec_loose = MakeTH1("hMCRecoLevPair_XiRec_loose",9,binning2);
  TH2F* hRPM_XiRec_stand_un = MakeTH2("hRPM_XiRec_stand_un",9,binning2);
  TH2F* hRPM_XiRec_stand = MakeTH2("hRPM_XiRec_stand",9,binning2);
  TH1F* hMCRecoLevXic0_XiRec_stand = MakeTH1("hMCRecoLevXic0_XiRec_stand",9,binning2);
  TH1F* hMCRecoLevPair_XiRec_stand = MakeTH1("hMCRecoLevPair_XiRec_stand",9,binning2);
  TH2F* hRPM_XiRec_tight_un = MakeTH2("hRPM_XiRec_tight_un",9,binning2);
  TH2F* hRPM_XiRec_tight = MakeTH2("hRPM_XiRec_tight",9,binning2);
  TH1F* hMCRecoLevXic0_XiRec_tight = MakeTH1("hMCRecoLevXic0_XiRec_tight",9,binning2);
  TH1F* hMCRecoLevPair_XiRec_tight = MakeTH1("hMCRecoLevPair_XiRec_tight",9,binning2);
  TH2F* hRPM_XiRec_vloose_un = MakeTH2("hRPM_XiRec_vloose_un",9,binning2);
  TH2F* hRPM_XiRec_vloose = MakeTH2("hRPM_XiRec_vloose",9,binning2);
  TH1F* hMCRecoLevXic0_XiRec_vloose = MakeTH1("hMCRecoLevXic0_XiRec_vloose",9,binning2);
  TH1F* hMCRecoLevPair_XiRec_vloose = MakeTH1("hMCRecoLevPair_XiRec_vloose",9,binning2);
  TH2F* hRPM_XiRec_vtight_un = MakeTH2("hRPM_XiRec_vtight_un",9,binning2);
  TH2F* hRPM_XiRec_vtight = MakeTH2("hRPM_XiRec_vtight",9,binning2);
  TH1F* hMCRecoLevXic0_XiRec_vtight = MakeTH1("hMCRecoLevXic0_XiRec_vtight",9,binning2);
  TH1F* hMCRecoLevPair_XiRec_vtight = MakeTH1("hMCRecoLevPair_XiRec_vtight",9,binning2);

  TH2F* hRPM_XiPID_loose_un = MakeTH2("hRPM_XiPID_loose_un",9,binning2);
  TH2F* hRPM_XiPID_loose = MakeTH2("hRPM_XiPID_loose",9,binning2);
  TH1F* hMCRecoLevXic0_XiPID_loose = MakeTH1("hMCRecoLevXic0_XiPID_loose",9,binning2);
  TH1F* hMCRecoLevPair_XiPID_loose = MakeTH1("hMCRecoLevPair_XiPID_loose",9,binning2);
  TH2F* hRPM_XiPID_stand_un = MakeTH2("hRPM_XiPID_stand_un",9,binning2);
  TH2F* hRPM_XiPID_stand = MakeTH2("hRPM_XiPID_stand",9,binning2);
  TH1F* hMCRecoLevXic0_XiPID_stand = MakeTH1("hMCRecoLevXic0_XiPID_stand",9,binning2);
  TH1F* hMCRecoLevPair_XiPID_stand = MakeTH1("hMCRecoLevPair_XiPID_stand",9,binning2);
  TH2F* hRPM_XiPID_tight_un = MakeTH2("hRPM_XiPID_tight_un",9,binning2);
  TH2F* hRPM_XiPID_tight = MakeTH2("hRPM_XiPID_tight",9,binning2);
  TH1F* hMCRecoLevXic0_XiPID_tight = MakeTH1("hMCRecoLevXic0_XiPID_tight",9,binning2);
  TH1F* hMCRecoLevPair_XiPID_tight = MakeTH1("hMCRecoLevPair_XiPID_tight",9,binning2);
  TH2F* hRPM_XiPID_vloose_un = MakeTH2("hRPM_XiPID_vloose_un",9,binning2);
  TH2F* hRPM_XiPID_vloose = MakeTH2("hRPM_XiPID_vloose",9,binning2);
  TH1F* hMCRecoLevXic0_XiPID_vloose = MakeTH1("hMCRecoLevXic0_XiPID_vloose",9,binning2);
  TH1F* hMCRecoLevPair_XiPID_vloose = MakeTH1("hMCRecoLevPair_XiPID_vloose",9,binning2);
  TH2F* hRPM_XiPID_vtight_un = MakeTH2("hRPM_XiPID_vtight_un",9,binning2);
  TH2F* hRPM_XiPID_vtight = MakeTH2("hRPM_XiPID_vtight",9,binning2);
  TH1F* hMCRecoLevXic0_XiPID_vtight = MakeTH1("hMCRecoLevXic0_XiPID_vtight",9,binning2);
  TH1F* hMCRecoLevPair_XiPID_vtight = MakeTH1("hMCRecoLevPair_XiPID_vtight",9,binning2);

  TH2F* hRPM_op_loose_un = MakeTH2("hRPM_OA_loose_un",9,binning2);
  TH2F* hRPM_op_loose = MakeTH2("hRPM_OA_loose",9,binning2);
  TH1F* hMCRecoLevXic0_op_loose = MakeTH1("hMCRecoLevXic0_OA_loose",9,binning2);
  TH1F* hMCRecoLevPair_op_loose = MakeTH1("hMCRecoLevPair_OA_loose",9,binning2);
  TH2F* hRPM_op_stand_un = MakeTH2("hRPM_OA_stand_un",9,binning2);
  TH2F* hRPM_op_stand = MakeTH2("hRPM_OA_stand",9,binning2);
  TH1F* hMCRecoLevXic0_op_stand = MakeTH1("hMCRecoLevXic0_OA_stand",9,binning2);
  TH1F* hMCRecoLevPair_op_stand = MakeTH1("hMCRecoLevPair_OA_stand",9,binning2);
  TH2F* hRPM_op_tight_un = MakeTH2("hRPM_OA_tight_un",9,binning2);
  TH2F* hRPM_op_tight = MakeTH2("hRPM_OA_tight",9,binning2);
  TH1F* hMCRecoLevXic0_op_tight = MakeTH1("hMCRecoLevXic0_OA_tight",9,binning2);
  TH1F* hMCRecoLevPair_op_tight = MakeTH1("hMCRecoLevPair_OA_tight",9,binning2);

  TH2F* hRPM_im_loose_un = MakeTH2("hRPM_IM_loose_un",9,binning2);
  TH2F* hRPM_im_loose = MakeTH2("hRPM_IM_loose",9,binning2);
  TH1F* hMCRecoLevXic0_im_loose = MakeTH1("hMCRecoLevXic0_IM_loose",9,binning2);
  TH1F* hMCRecoLevPair_im_loose = MakeTH1("hMCRecoLevPair_IM_loose",9,binning2);
  TH2F* hRPM_im_stand_un = MakeTH2("hRPM_IM_stand_un",9,binning2);
  TH2F* hRPM_im_stand = MakeTH2("hRPM_IM_stand",9,binning2);
  TH1F* hMCRecoLevXic0_im_stand = MakeTH1("hMCRecoLevXic0_IM_stand",9,binning2);
  TH1F* hMCRecoLevPair_im_stand = MakeTH1("hMCRecoLevPair_IM_stand",9,binning2);
  TH2F* hRPM_im_tight_un = MakeTH2("hRPM_IM_tight_un",9,binning2);
  TH2F* hRPM_im_tight = MakeTH2("hRPM_IM_tight",9,binning2);
  TH1F* hMCRecoLevXic0_im_tight = MakeTH1("hMCRecoLevXic0_IM_tight",9,binning2);
  TH1F* hMCRecoLevPair_im_tight = MakeTH1("hMCRecoLevPair_IM_tight",9,binning2);

  TH2F* hRPM_eRec_loose_Xib = MakeTH2("hRPM_eRec_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_eRec_loose = MakeTH1("hMCRecoLevXib_eRec_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_eRec_loose = MakeTH1("hMCRecoLevPairXib_eRec_loose",9,binning2);
  TH2F* hRPM_eRec_vloose_Xib = MakeTH2("hRPM_eRec_vloose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_eRec_vloose = MakeTH1("hMCRecoLevXib_eRec_vloose",9,binning2);
  TH1F* hMCRecoLevPairXib_eRec_vloose = MakeTH1("hMCRecoLevPairXib_eRec_vloose",9,binning2);
  TH2F* hRPM_eRec_stand_Xib = MakeTH2("hRPM_eRec_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_eRec_stand = MakeTH1("hMCRecoLevXib_eRec_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_eRec_stand = MakeTH1("hMCRecoLevPairXib_eRec_stand",9,binning2);
  TH2F* hRPM_eRec_tight_Xib = MakeTH2("hRPM_eRec_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_eRec_tight = MakeTH1("hMCRecoLevXib_eRec_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_eRec_tight = MakeTH1("hMCRecoLevPairXib_eRec_tight",9,binning2);
  TH2F* hRPM_eRec_vtight_Xib = MakeTH2("hRPM_eRec_vtight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_eRec_vtight = MakeTH1("hMCRecoLevXib_eRec_vtight",9,binning2);
  TH1F* hMCRecoLevPairXib_eRec_vtight = MakeTH1("hMCRecoLevPairXib_eRec_vtight",9,binning2);

  TH2F* hRPM_ePID_loose_Xib = MakeTH2("hRPM_ePID_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_ePID_loose = MakeTH1("hMCRecoLevXib_ePID_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_ePID_loose = MakeTH1("hMCRecoLevPairXib_ePID_loose",9,binning2);
  TH2F* hRPM_ePID_vloose_Xib = MakeTH2("hRPM_ePID_vloose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_ePID_vloose = MakeTH1("hMCRecoLevXib_ePID_vloose",9,binning2);
  TH1F* hMCRecoLevPairXib_ePID_vloose = MakeTH1("hMCRecoLevPairXib_ePID_vloose",9,binning2);
  TH2F* hRPM_ePID_stand_Xib = MakeTH2("hRPM_ePID_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_ePID_stand = MakeTH1("hMCRecoLevXib_ePID_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_ePID_stand = MakeTH1("hMCRecoLevPairXib_ePID_stand",9,binning2);
  TH2F* hRPM_ePID_tight_Xib = MakeTH2("hRPM_ePID_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_ePID_tight = MakeTH1("hMCRecoLevXib_ePID_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_ePID_tight = MakeTH1("hMCRecoLevPairXib_ePID_tight",9,binning2);
  TH2F* hRPM_ePID_vtight_Xib = MakeTH2("hRPM_ePID_vtight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_ePID_vtight = MakeTH1("hMCRecoLevXib_ePID_vtight",9,binning2);
  TH1F* hMCRecoLevPairXib_ePID_vtight = MakeTH1("hMCRecoLevPairXib_ePID_vtight",9,binning2);

  TH2F* hRPM_XiPID_loose_Xib = MakeTH2("hRPM_XiPID_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiPID_loose = MakeTH1("hMCRecoLevXib_XiPID_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_XiPID_loose = MakeTH1("hMCRecoLevPairXib_XiPID_loose",9,binning2);
  TH2F* hRPM_XiPID_vloose_Xib = MakeTH2("hRPM_XiPID_vloose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiPID_vloose = MakeTH1("hMCRecoLevXib_XiPID_vloose",9,binning2);
  TH1F* hMCRecoLevPairXib_XiPID_vloose = MakeTH1("hMCRecoLevPairXib_XiPID_vloose",9,binning2);
  TH2F* hRPM_XiPID_stand_Xib = MakeTH2("hRPM_XiPID_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiPID_stand = MakeTH1("hMCRecoLevXib_XiPID_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_XiPID_stand = MakeTH1("hMCRecoLevPairXib_XiPID_stand",9,binning2);
  TH2F* hRPM_XiPID_tight_Xib = MakeTH2("hRPM_XiPID_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiPID_tight = MakeTH1("hMCRecoLevXib_XiPID_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_XiPID_tight = MakeTH1("hMCRecoLevPairXib_XiPID_tight",9,binning2);
  TH2F* hRPM_XiPID_vtight_Xib = MakeTH2("hRPM_XiPID_vtight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiPID_vtight = MakeTH1("hMCRecoLevXib_XiPID_vtight",9,binning2);
  TH1F* hMCRecoLevPairXib_XiPID_vtight = MakeTH1("hMCRecoLevPairXib_XiPID_vtight",9,binning2);

  TH2F* hRPM_XiRec_loose_Xib = MakeTH2("hRPM_XiRec_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiRec_loose = MakeTH1("hMCRecoLevXib_XiRec_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_XiRec_loose = MakeTH1("hMCRecoLevPairXib_XiRec_loose",9,binning2);
  TH2F* hRPM_XiRec_vloose_Xib = MakeTH2("hRPM_XiRec_vloose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiRec_vloose = MakeTH1("hMCRecoLevXib_XiRec_vloose",9,binning2);
  TH1F* hMCRecoLevPairXib_XiRec_vloose = MakeTH1("hMCRecoLevPairXib_XiRec_vloose",9,binning2);
  TH2F* hRPM_XiRec_stand_Xib = MakeTH2("hRPM_XiRec_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiRec_stand = MakeTH1("hMCRecoLevXib_XiRec_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_XiRec_stand = MakeTH1("hMCRecoLevPairXib_XiRec_stand",9,binning2);
  TH2F* hRPM_XiRec_tight_Xib = MakeTH2("hRPM_XiRec_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiRec_tight = MakeTH1("hMCRecoLevXib_XiRec_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_XiRec_tight = MakeTH1("hMCRecoLevPairXib_XiRec_tight",9,binning2);
  TH2F* hRPM_XiRec_vtight_Xib = MakeTH2("hRPM_XiRec_vtight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_XiRec_vtight = MakeTH1("hMCRecoLevXib_XiRec_vtight",9,binning2);
  TH1F* hMCRecoLevPairXib_XiRec_vtight = MakeTH1("hMCRecoLevPairXib_XiRec_vtight",9,binning2);

  TH2F* hRPM_im_loose_Xib = MakeTH2("hRPM_IM_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_im_loose = MakeTH1("hMCRecoLevXib_IM_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_im_loose = MakeTH1("hMCRecoLevPairXib_IM_loose",9,binning2);
  TH2F* hRPM_im_stand_Xib = MakeTH2("hRPM_IM_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_im_stand = MakeTH1("hMCRecoLevXib_IM_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_im_stand = MakeTH1("hMCRecoLevPairXib_IM_stand",9,binning2);
  TH2F* hRPM_im_tight_Xib = MakeTH2("hRPM_IM_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_im_tight = MakeTH1("hMCRecoLevXib_IM_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_im_tight = MakeTH1("hMCRecoLevPairXib_IM_tight",9,binning2);

  TH2F* hRPM_op_loose_Xib = MakeTH2("hRPM_OA_loose_Xib",9,binning2);
  TH1F* hMCRecoLevXib_op_loose = MakeTH1("hMCRecoLevXib_OA_loose",9,binning2);
  TH1F* hMCRecoLevPairXib_op_loose = MakeTH1("hMCRecoLevPairXib_OA_loose",9,binning2);
  TH2F* hRPM_op_stand_Xib = MakeTH2("hRPM_OA_stand_Xib",9,binning2);
  TH1F* hMCRecoLevXib_op_stand = MakeTH1("hMCRecoLevXib_OA_stand",9,binning2);
  TH1F* hMCRecoLevPairXib_op_stand = MakeTH1("hMCRecoLevPairXib_OA_stand",9,binning2);
  TH2F* hRPM_op_tight_Xib = MakeTH2("hRPM_OA_tight_Xib",9,binning2);
  TH1F* hMCRecoLevXib_op_tight = MakeTH1("hMCRecoLevXib_OA_tight",9,binning2);
  TH1F* hMCRecoLevPairXib_op_tight = MakeTH1("hMCRecoLevPairXib_OA_tight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_eRec_stand = MakeTH1("hMCRecoLevPromptXic0_eRec_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_eRec_vloose = MakeTH1("hMCRecoLevPromptXic0_eRec_vloose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_eRec_loose = MakeTH1("hMCRecoLevPromptXic0_eRec_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_eRec_tight = MakeTH1("hMCRecoLevPromptXic0_eRec_tight",9,binning2);
  TH1F* hMCRecoLevPromptXic0_eRec_vtight = MakeTH1("hMCRecoLevPromptXic0_eRec_vtight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_XiRec_stand = MakeTH1("hMCRecoLevPromptXic0_XiRec_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiRec_vloose = MakeTH1("hMCRecoLevPromptXic0_XiRec_vloose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiRec_loose = MakeTH1("hMCRecoLevPromptXic0_XiRec_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiRec_tight = MakeTH1("hMCRecoLevPromptXic0_XiRec_tight",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiRec_vtight = MakeTH1("hMCRecoLevPromptXic0_XiRec_vtight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_ePID_stand = MakeTH1("hMCRecoLevPromptXic0_ePID_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_ePID_vloose = MakeTH1("hMCRecoLevPromptXic0_ePID_vloose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_ePID_loose = MakeTH1("hMCRecoLevPromptXic0_ePID_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_ePID_tight = MakeTH1("hMCRecoLevPromptXic0_ePID_tight",9,binning2);
  TH1F* hMCRecoLevPromptXic0_ePID_vtight = MakeTH1("hMCRecoLevPromptXic0_ePID_vtight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_XiPID_stand = MakeTH1("hMCRecoLevPromptXic0_XiPID_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiPID_vloose = MakeTH1("hMCRecoLevPromptXic0_XiPID_vloose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiPID_loose = MakeTH1("hMCRecoLevPromptXic0_XiPID_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiPID_tight = MakeTH1("hMCRecoLevPromptXic0_XiPID_tight",9,binning2);
  TH1F* hMCRecoLevPromptXic0_XiPID_vtight = MakeTH1("hMCRecoLevPromptXic0_XiPID_vtight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_IM_stand = MakeTH1("hMCRecoLevPromptXic0_IM_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_IM_loose = MakeTH1("hMCRecoLevPromptXic0_IM_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_IM_tight = MakeTH1("hMCRecoLevPromptXic0_IM_tight",9,binning2);

  TH1F* hMCRecoLevPromptXic0_OA_stand = MakeTH1("hMCRecoLevPromptXic0_OA_stand",9,binning2);
  TH1F* hMCRecoLevPromptXic0_OA_loose = MakeTH1("hMCRecoLevPromptXic0_OA_loose",9,binning2);
  TH1F* hMCRecoLevPromptXic0_OA_tight = MakeTH1("hMCRecoLevPromptXic0_OA_tight",9,binning2);

  TH1F *hprompt = MakeTH1("hprompt",7,binning1); TH1F *hnonprompt = MakeTH1("hnonprompt",7,binning1); TH1F *hinclu = MakeTH1("hinclu",7,binning1);
  TH2F* hRPM_preliminary = new TH2F("hWRPM_preliminary","",60,1,20,60,1,20); hRPM_preliminary->Sumw2();

  TH1F* htest1 = MakeTH1("test1",9,binning2);
  TH1F* htest2 = MakeTH1("test2",9,binning2);
  TH2F* htest3 = MakeTH2("test3",9,binning2);
  TH1F* htest4 = new TH1F("invariant mass","",30,1,6); htest4->Sumw2();
  TH1F* htest5 = new TH1F("opening angle","",30,-1,1); htest5->Sumw2();

        Int_t nTrks = Pair->GetEntriesFast();

        for (Int_t i=0; i<nTrks; i++) //loose cut
        {
          Pair->GetEntry(i);
          MCTree->GetEntry(i);

          if(fabs(Massv-1.32171)>0.008) continue;  //Xi mass tolerance
          if(In_Mass<1.3) continue;  //pair mass low limit
          if(fabs(pTe==999)) continue;  //dummy tree reject

          Bool_t isparticle = kTRUE; if(vcharge>0) isparticle = kFALSE;
          Float_t VT_e_nsigma_cut = -3.5+(1.15*pTe)-(0.09*pTe*pTe);
          Float_t T_e_nsigma_cut = -3.7+(1.17*pTe)-(0.094*pTe*pTe);  ///need to modify
          Float_t S_e_nsigma_cut = -3.9+(1.17*pTe)-(0.094*pTe*pTe);
          Float_t L_e_nsigma_cut = -4.1+(1.17*pTe)-(0.094*pTe*pTe); ///need to modify
          Float_t VL_e_nsigma_cut = -4.3+(1.17*pTe)-(0.094*pTe*pTe);

          if(pTe>=5){
            T_e_nsigma_cut = -3.7+(1.17*5)-(0.094*25);
            VT_e_nsigma_cut = -3.5+(1.15*5)-(0.09*25);
            S_e_nsigma_cut = -3.9+(1.17*5)-(0.094*25);
          VL_e_nsigma_cut = -4.3+(1.17*5)-(0.094*25);
            L_e_nsigma_cut = -4.1+(1.17*5)-(0.094*25);
          }

          Bool_t Xi_Topology_VLoose_flag = kFALSE; if(V0DecayLength>0.02 && CascDecayLength>0.02 && DCABachToPrimVertex>0.01
            && DCANegToPrimVertex>0.05 && DCAPosToPrimVertex>0.05 && V0CosineOfPoiningAngleXi>0.98 && XiCosineOfPoiningAngle>0.98 &&DCAV0ToPrimVertex>0.01) Xi_Topology_VLoose_flag = kTRUE;
          Bool_t Xi_Topology_Loose_flag = kFALSE; if(V0DecayLength>1.55 && CascDecayLength>0.29 && DCABachToPrimVertex>0.0146
            && DCANegToPrimVertex>0.061 && DCAPosToPrimVertex>0.061 && V0CosineOfPoiningAngleXi>0.981 && XiCosineOfPoiningAngle>0.981 &&DCAV0ToPrimVertex>0.02) Xi_Topology_Loose_flag = kTRUE;
          Bool_t Xi_Topology_Stand_flag = kFALSE; if(V0DecayLength>2.67 && CascDecayLength>0.38 && DCABachToPrimVertex>0.0204
            && DCANegToPrimVertex>0.073 && DCAPosToPrimVertex>0.073 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&DCAV0ToPrimVertex>0.03) Xi_Topology_Stand_flag = kTRUE;
          Bool_t Xi_Topology_Tight_flag = kFALSE; if(V0DecayLength>3.6 && CascDecayLength>0.53 && DCABachToPrimVertex>0.037
            && DCANegToPrimVertex>0.088 && DCAPosToPrimVertex>0.088 && V0CosineOfPoiningAngleXi>0.9839 && XiCosineOfPoiningAngle>0.9839 &&DCAV0ToPrimVertex>0.04) Xi_Topology_Tight_flag = kTRUE;
          Bool_t Xi_Topology_VTight_flag = kFALSE; if(V0DecayLength>4.39 && CascDecayLength>0.72 && DCABachToPrimVertex>0.037
            && DCANegToPrimVertex>0.102 && DCAPosToPrimVertex>0.102 && V0CosineOfPoiningAngleXi>0.985 && XiCosineOfPoiningAngle>0.985 &&DCAV0ToPrimVertex>0.06) Xi_Topology_VTight_flag = kTRUE;

          Bool_t Xi_Recon_VTight_flag = kFALSE;
          Bool_t Xi_Recon_Tight_flag = kFALSE;
          Bool_t Xi_Recon_Stand_flag = kFALSE;
          Bool_t Xi_Recon_Loose_flag = kFALSE;
          Bool_t Xi_Recon_VLoose_flag = kFALSE;
          if(pion_crossedratio/pion_findable>0.81 && proton_crossedratio/proton_findable>0.81 && bpion_crossedratio/bpion_findable>0.81 && pion_crossedratio>80 && proton_crossedratio>80 && bpion_crossedratio>80) Xi_Recon_VTight_flag = kTRUE;
          if(pion_crossedratio/pion_findable>0.79 && proton_crossedratio/proton_findable>0.79 && bpion_crossedratio/bpion_findable>0.79 && pion_crossedratio>75 && proton_crossedratio>75 && bpion_crossedratio>75) Xi_Recon_Tight_flag = kTRUE;
          if(pion_crossedratio/pion_findable>0.77 && proton_crossedratio/proton_findable>0.77 && bpion_crossedratio/bpion_findable>0.77 && pion_crossedratio>70 && proton_crossedratio>70 && bpion_crossedratio>70) Xi_Recon_Stand_flag = kTRUE;
          if(pion_crossedratio/pion_findable>0.75 && proton_crossedratio/proton_findable>0.75 && bpion_crossedratio/bpion_findable>0.75 && pion_crossedratio>65 && proton_crossedratio>65 && bpion_crossedratio>65) Xi_Recon_Loose_flag = kTRUE;
          if(pion_crossedratio/pion_findable>0.70 && proton_crossedratio/proton_findable>0.70 && bpion_crossedratio/bpion_findable>0.70 && pion_crossedratio>65 && proton_crossedratio>65 && bpion_crossedratio>65) Xi_Recon_VLoose_flag = kTRUE;

          Bool_t e_Recon_VTight_flag = kFALSE;  if(e_crossedratio>85 && TPCPIDCluster>60 && e_crossedratio/e_findable>0.9 && ITSCluster>=3) e_Recon_VTight_flag = kTRUE;
          Bool_t e_Recon_Tight_flag = kFALSE;  if(e_crossedratio>75 && TPCPIDCluster>55 && e_crossedratio/e_findable>0.85 && ITSCluster>=3) e_Recon_Tight_flag = kTRUE;
          Bool_t e_Recon_Stand_flag = kFALSE; if(e_crossedratio>70 && TPCPIDCluster>50 && e_crossedratio/e_findable>0.8 && ITSCluster>=3) e_Recon_Stand_flag = kTRUE;
          Bool_t e_Recon_Loose_flag = kFALSE;  if(e_crossedratio>65 && TPCPIDCluster>45 && e_crossedratio/e_findable>0.8 && ITSCluster>=3) e_Recon_Loose_flag = kTRUE;
          Bool_t e_Recon_VLoose_flag = kFALSE;  if(e_crossedratio>65 && TPCPIDCluster>40 && e_crossedratio/e_findable>0.75 && ITSCluster>=3) e_Recon_VLoose_flag = kTRUE;

          Bool_t e_PID_VTight_flag = kFALSE; if(fabs(nSigmaTOF)<=2 && nSigmaTPC>=VT_e_nsigma_cut && nSigmaTPC<=3) e_PID_VTight_flag = kTRUE;
          Bool_t e_PID_Tight_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=T_e_nsigma_cut && nSigmaTPC<=3) e_PID_Tight_flag = kTRUE;
          Bool_t e_PID_Stand_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=S_e_nsigma_cut && nSigmaTPC<=3) e_PID_Stand_flag = kTRUE;
          Bool_t e_PID_Loose_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=L_e_nsigma_cut && nSigmaTPC<=3) e_PID_Loose_flag = kTRUE;
          Bool_t e_PID_VLoose_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=VL_e_nsigma_cut && nSigmaTPC<=3) e_PID_VLoose_flag = kTRUE;

          Bool_t OPAngle_Tight_flag = kFALSE; if(cosoa>cos(70*(3.141592/180))) OPAngle_Tight_flag = kTRUE;
          Bool_t OPAngle_Stand_flag = kFALSE; if(cosoa>cos(90*(3.141592/180))) OPAngle_Stand_flag = kTRUE;
          Bool_t OPAngle_Loose_flag = kFALSE; if(cosoa>cos(110*(3.141592/180)))OPAngle_Loose_flag = kTRUE;

          Bool_t PairMass_Tight_flag = kFALSE; if(In_Mass<2.3) PairMass_Tight_flag = kTRUE;
          Bool_t PairMass_Stand_flag = kFALSE; if(In_Mass<2.5) PairMass_Stand_flag = kTRUE;
          Bool_t PairMass_Loose_flag = kFALSE; if(In_Mass<2.7) PairMass_Loose_flag = kTRUE;

          if(echarge*vcharge<0 && e_minmass>0.05) //Right Sign
          {
            if(Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && e_Recon_Stand_flag) hMassPtRS->Fill(In_Mass,Pt);
            if(Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_Recon_VTight_flag) hPtRS_VTeRec->Fill(Pt);
              if(e_Recon_Tight_flag) hPtRS_TeRec->Fill(Pt);
              if(e_Recon_Stand_flag) hPtRS_SeRec->Fill(Pt);
              if(e_Recon_Loose_flag) hPtRS_LeRec->Fill(Pt);
              if(e_Recon_VLoose_flag) hPtRS_VLeRec->Fill(Pt);
            }
            if(e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Recon_VTight_flag) hPtRS_VTXiRec->Fill(Pt);
              if(Xi_Recon_Tight_flag) hPtRS_TXiRec->Fill(Pt);
              if(Xi_Recon_Stand_flag) hPtRS_SXiRec->Fill(Pt);
              if(Xi_Recon_Loose_flag) hPtRS_LXiRec->Fill(Pt);
              if(Xi_Recon_VLoose_flag) hPtRS_VLXiRec->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_PID_VTight_flag) hPtRS_VTePID->Fill(Pt);
              if(e_PID_Tight_flag) hPtRS_TePID->Fill(Pt);
              if(e_PID_Stand_flag) hPtRS_SePID->Fill(Pt);
              if(e_PID_Loose_flag) hPtRS_LePID->Fill(Pt);
              if(e_PID_VLoose_flag) hPtRS_VLePID->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Topology_VTight_flag) hPtRS_VTXiPID->Fill(Pt);
              if(Xi_Topology_Tight_flag) hPtRS_TXiPID->Fill(Pt);
              if(Xi_Topology_Stand_flag) hPtRS_SXiPID->Fill(Pt);
              if(Xi_Topology_Loose_flag) hPtRS_LXiPID->Fill(Pt);
              if(Xi_Topology_VLoose_flag) hPtRS_VLXiPID->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && PairMass_Stand_flag){
              if(OPAngle_Tight_flag) hPtRS_TOA->Fill(Pt);
              if(OPAngle_Stand_flag) hPtRS_SOA->Fill(Pt);
              if(OPAngle_Loose_flag) hPtRS_LOA->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
              if(PairMass_Tight_flag) hPtRS_TeIM->Fill(Pt);
              if(PairMass_Stand_flag) hPtRS_SeIM->Fill(Pt);
              if(PairMass_Loose_flag) hPtRS_LeIM->Fill(Pt);
            }
          } //RS

          if(echarge*vcharge>0 && e_minmass>0.05) //Wrong Sign
          {
            if(Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && e_Recon_Stand_flag) hMassPtWS->Fill(In_Mass,Pt);
            if(Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_Recon_VTight_flag) hPtWS_VTeRec->Fill(Pt);
              if(e_Recon_Tight_flag) hPtWS_TeRec->Fill(Pt);
              if(e_Recon_Stand_flag) hPtWS_SeRec->Fill(Pt);
              if(e_Recon_Loose_flag) hPtWS_LeRec->Fill(Pt);
              if(e_Recon_VLoose_flag) hPtWS_VLeRec->Fill(Pt);
            }
            if(e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Recon_VTight_flag) hPtWS_VTXiRec->Fill(Pt);
              if(Xi_Recon_Tight_flag) hPtWS_TXiRec->Fill(Pt);
              if(Xi_Recon_Stand_flag) hPtWS_SXiRec->Fill(Pt);
              if(Xi_Recon_Loose_flag) hPtWS_LXiRec->Fill(Pt);
              if(Xi_Recon_VLoose_flag) hPtWS_VLXiRec->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_PID_VTight_flag) hPtWS_VTePID->Fill(Pt);
              if(e_PID_Tight_flag) hPtWS_TePID->Fill(Pt);
              if(e_PID_Stand_flag) hPtWS_SePID->Fill(Pt);
              if(e_PID_Loose_flag) hPtWS_LePID->Fill(Pt);
              if(e_PID_VLoose_flag) hPtWS_VLePID->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Topology_VTight_flag) hPtWS_VTXiPID->Fill(Pt);
              if(Xi_Topology_Tight_flag) hPtWS_TXiPID->Fill(Pt);
              if(Xi_Topology_Stand_flag) hPtWS_SXiPID->Fill(Pt);
              if(Xi_Topology_Loose_flag) hPtWS_LXiPID->Fill(Pt);
              if(Xi_Topology_VLoose_flag) hPtWS_VLXiPID->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && PairMass_Stand_flag){
              if(OPAngle_Tight_flag) hPtWS_TOA->Fill(Pt);
              if(OPAngle_Stand_flag) hPtWS_SOA->Fill(Pt);
              if(OPAngle_Loose_flag) hPtWS_LOA->Fill(Pt);
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
              if(PairMass_Tight_flag) hPtWS_TeIM->Fill(Pt);
              if(PairMass_Stand_flag) hPtWS_SeIM->Fill(Pt);
              if(PairMass_Loose_flag) hPtWS_LeIM->Fill(Pt);
            }
          } //WS

          if(echarge*vcharge>0 && PairMass_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
            if(e_Recon_VTight_flag) hpre_eRec_vtight_de->Fill(Pt);
            if(e_Recon_Tight_flag) hpre_eRec_tight_de->Fill(Pt);
            if(e_Recon_Stand_flag) hpre_eRec_stand_de->Fill(Pt);
            if(e_Recon_Loose_flag) hpre_eRec_loose_de->Fill(Pt);
            if(e_Recon_VLoose_flag) hpre_eRec_vloose_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(e_Recon_VTight_flag) hpre_eRec_vtight_nu->Fill(Pt);
              if(e_Recon_Tight_flag) hpre_eRec_tight_nu->Fill(Pt);
              if(e_Recon_Stand_flag) hpre_eRec_stand_nu->Fill(Pt);
              if(e_Recon_Loose_flag) hpre_eRec_loose_nu->Fill(Pt);
              if(e_Recon_VLoose_flag) hpre_eRec_vloose_nu->Fill(Pt);
            }
          }
          if(echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
            if(Xi_Recon_VLoose_flag) hpre_XiRec_vloose_de->Fill(Pt);
            if(Xi_Recon_Loose_flag) hpre_XiRec_loose_de->Fill(Pt);
            if(Xi_Recon_Stand_flag) hpre_XiRec_stand_de->Fill(Pt);
            if(Xi_Recon_Tight_flag) hpre_XiRec_tight_de->Fill(Pt);
            if(Xi_Recon_VTight_flag) hpre_XiRec_vtight_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(Xi_Recon_VLoose_flag) hpre_XiRec_vloose_nu->Fill(Pt);
              if(Xi_Recon_Loose_flag) hpre_XiRec_loose_nu->Fill(Pt);
              if(Xi_Recon_Stand_flag) hpre_XiRec_stand_nu->Fill(Pt);
              if(Xi_Recon_Tight_flag) hpre_XiRec_tight_nu->Fill(Pt);
              if(Xi_Recon_VTight_flag) hpre_XiRec_vtight_nu->Fill(Pt);
            }
          }
          if(echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
            if(e_PID_VLoose_flag) hpre_ePID_vloose_de->Fill(Pt);
            if(e_PID_Loose_flag) hpre_ePID_loose_de->Fill(Pt);
            if(e_PID_Stand_flag) hpre_ePID_stand_de->Fill(Pt);
            if(e_PID_Tight_flag) hpre_ePID_tight_de->Fill(Pt);
            if(e_PID_VTight_flag) hpre_ePID_vtight_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(e_PID_VLoose_flag) hpre_ePID_vloose_nu->Fill(Pt);
              if(e_PID_Loose_flag) hpre_ePID_loose_nu->Fill(Pt);
              if(e_PID_Stand_flag) hpre_ePID_stand_nu->Fill(Pt);
              if(e_PID_Tight_flag) hpre_ePID_tight_nu->Fill(Pt);
              if(e_PID_VTight_flag) hpre_ePID_vtight_nu->Fill(Pt);
            }
          }
          if(echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag){
            if(Xi_Topology_VLoose_flag) hpre_XiPID_vloose_de->Fill(Pt);
            if(Xi_Topology_Loose_flag) hpre_XiPID_loose_de->Fill(Pt);
            if(Xi_Topology_Stand_flag) hpre_XiPID_stand_de->Fill(Pt);
            if(Xi_Topology_Tight_flag) hpre_XiPID_tight_de->Fill(Pt);
            if(Xi_Topology_VTight_flag) hpre_XiPID_vtight_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(Xi_Topology_VLoose_flag) hpre_XiPID_vloose_nu->Fill(Pt);
              if(Xi_Topology_Loose_flag) hpre_XiPID_loose_nu->Fill(Pt);
              if(Xi_Topology_Stand_flag) hpre_XiPID_stand_nu->Fill(Pt);
              if(Xi_Topology_Tight_flag) hpre_XiPID_tight_nu->Fill(Pt);
              if(Xi_Topology_VTight_flag) hpre_XiPID_vtight_nu->Fill(Pt);
            }
          }
          if(echarge*vcharge>0 && PairMass_Stand_flag && e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag){
            if(OPAngle_Loose_flag) hpre_OA_loose_de->Fill(Pt);
            if(OPAngle_Stand_flag) hpre_OA_stand_de->Fill(Pt);
            if(OPAngle_Tight_flag) hpre_OA_tight_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(OPAngle_Loose_flag) hpre_OA_loose_nu->Fill(Pt);
              if(OPAngle_Stand_flag) hpre_OA_stand_nu->Fill(Pt);
              if(OPAngle_Tight_flag) hpre_OA_tight_nu->Fill(Pt);
            }
          }
          if(echarge*vcharge>0 && OPAngle_Stand_flag && e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag){
            if(PairMass_Loose_flag) hpre_IM_loose_de->Fill(Pt);
            if(PairMass_Stand_flag) hpre_IM_stand_de->Fill(Pt);
            if(PairMass_Tight_flag) hpre_IM_tight_de->Fill(Pt);
            if(e_minmass_ss>0.05){
              if(PairMass_Loose_flag) hpre_IM_loose_nu->Fill(Pt);
              if(PairMass_Stand_flag) hpre_IM_stand_nu->Fill(Pt);
              if(PairMass_Tight_flag) hpre_IM_tight_nu->Fill(Pt);
            }
          } //WS

          if(IsMC){
            if(echarge*vcharge>0 && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && PairMass_Stand_flag && OPAngle_Stand_flag && e_Recon_Stand_flag){
              if(fabs(mcXib) == 9999.) continue;
              htest3->Fill(mcXib,XibeXi);
              htest1->Fill(XibeXi);
              htest2->Fill(mcXib);
            }
            if(fabs(mcXib) != 9999  && echarge*vcharge>0 && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && e_Recon_Stand_flag) htest4->Fill(In_Mass);
			if(fabs(mcXib) != 9999 && echarge*vcharge>0 && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && PairMass_Stand_flag && e_Recon_Stand_flag) htest5->Fill(cosoa);
            //if(fabs(mcptXic0) == 9999.) continue;

            if(Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_Recon_VTight_flag){
                if(fabs(mcptXic0) != 9999) hRPM_eRec_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_eRec_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_eRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_eRec_vtight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_eRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_eRec_vtight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_eRec_vtight->Fill(mcXib);
                  hMCRecoLevPairXib_eRec_vtight->Fill(XibeXi);
                }
              }
              if(e_Recon_Tight_flag){
                if(fabs(mcptXic0) != 9999) hRPM_eRec_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_eRec_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_eRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_eRec_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_eRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_eRec_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_eRec_tight->Fill(mcXib);
                  hMCRecoLevPairXib_eRec_tight->Fill(XibeXi);
                }
              }
              if(e_Recon_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_eRec_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_eRec_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_eRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_eRec_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_eRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_eRec_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_eRec_stand->Fill(mcXib);
                  hMCRecoLevPairXib_eRec_stand->Fill(XibeXi);
                }
              }
              if(e_Recon_Loose_flag){
                if(fabs(mcptXic0) != 9999) hRPM_eRec_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_eRec_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_eRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_eRec_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_eRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_eRec_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_eRec_loose->Fill(mcXib);
                  hMCRecoLevPairXib_eRec_loose->Fill(XibeXi);
                }
              }
              if(e_Recon_VLoose_flag){
                if(fabs(mcptXic0) != 9999) hRPM_eRec_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_eRec_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_eRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_eRec_vloose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_eRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_eRec_vloose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_eRec_vloose->Fill(mcXib);
                  hMCRecoLevPairXib_eRec_vloose->Fill(XibeXi);
                }
              }
            }
            if(e_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Recon_VTight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiRec_vtight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiRec_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiRec_vtight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiRec_vtight->Fill(mcXib);
                  hMCRecoLevPairXib_XiRec_vtight->Fill(XibeXi);
                }
              }
              if(Xi_Recon_Tight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiRec_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiRec_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiRec_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiRec_tight->Fill(mcXib);
                  hMCRecoLevPairXib_XiRec_tight->Fill(XibeXi);
                }
              }
              if(Xi_Recon_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiRec_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiRec_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiRec_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiRec_stand->Fill(mcXib);
                  hMCRecoLevPairXib_XiRec_stand->Fill(XibeXi);
                }
              }
              if(Xi_Recon_Loose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiRec_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiRec_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiRec_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiRec_loose->Fill(mcXib);
                  hMCRecoLevPairXib_XiRec_loose->Fill(XibeXi);
                }
              }
              if(Xi_Recon_VLoose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiRec_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiRec_vloose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiRec_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiRec_vloose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiRec_vloose->Fill(mcXib);
                  hMCRecoLevPairXib_XiRec_vloose->Fill(XibeXi);
                }
              }
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(e_PID_VTight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_ePID_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_ePID_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_ePID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_ePID_vtight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_ePID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_ePID_vtight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_ePID_vtight->Fill(mcXib);
                  hMCRecoLevPairXib_ePID_vtight->Fill(XibeXi);
                }
              }
              if(e_PID_Tight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_ePID_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_ePID_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_ePID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_ePID_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_ePID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_ePID_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_ePID_tight->Fill(mcXib);
                  hMCRecoLevPairXib_ePID_tight->Fill(XibeXi);
                }
              }
              if(e_PID_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_ePID_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_ePID_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_ePID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_ePID_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_ePID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_ePID_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_ePID_stand->Fill(mcXib);
                  hMCRecoLevPairXib_ePID_stand->Fill(XibeXi);
                }
              }
              if(e_PID_Loose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_ePID_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_ePID_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_ePID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_ePID_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_ePID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_ePID_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_ePID_loose->Fill(mcXib);
                  hMCRecoLevPairXib_ePID_loose->Fill(XibeXi);
                }
              }
              if(e_PID_VLoose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_ePID_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_ePID_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_ePID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_ePID_vloose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_ePID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_ePID_vloose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_ePID_vloose->Fill(mcXib);
                  hMCRecoLevPairXib_ePID_vloose->Fill(XibeXi);
                }
              }
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && OPAngle_Stand_flag && PairMass_Stand_flag){
              if(Xi_Topology_VTight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_vtight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_vtight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiPID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiPID_vtight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiPID_vtight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiPID_vtight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiPID_vtight->Fill(mcXib);
                  hMCRecoLevPairXib_XiPID_vtight->Fill(XibeXi);
                }
              }
              if(Xi_Topology_Tight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiPID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiPID_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiPID_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiPID_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiPID_tight->Fill(mcXib);
                  hMCRecoLevPairXib_XiPID_tight->Fill(XibeXi);
                }
              }
              if(Xi_Topology_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiPID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiPID_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiPID_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiPID_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiPID_stand->Fill(mcXib);
                  hMCRecoLevPairXib_XiPID_stand->Fill(XibeXi);
                }
              }
              if(Xi_Topology_Loose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiPID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiPID_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiPID_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiPID_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiPID_loose->Fill(mcXib);
                  hMCRecoLevPairXib_XiPID_loose->Fill(XibeXi);
                }
              }
              if(Xi_Topology_VLoose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_vloose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_XiPID_vloose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_XiPID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_XiPID_vloose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_XiPID_vloose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_XiPID_vloose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_XiPID_vloose->Fill(mcXib);
                  hMCRecoLevPairXib_XiPID_vloose->Fill(XibeXi);
                }
              }
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && PairMass_Stand_flag){
              if(OPAngle_Tight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_op_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_op_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_op_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_op_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_OA_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_op_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_op_tight->Fill(mcXib);
                  hMCRecoLevPairXib_op_tight->Fill(XibeXi);
                }
              }
              if(OPAngle_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_op_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_op_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_op_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_op_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_OA_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_op_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_op_stand->Fill(mcXib);
                  hMCRecoLevPairXib_op_stand->Fill(XibeXi);
                }
              }
              if(OPAngle_Loose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_op_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_op_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_op_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_op_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_OA_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_op_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_op_loose->Fill(mcXib);
                  hMCRecoLevPairXib_op_loose->Fill(XibeXi);
                }
              }
            }
            if(e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag && OPAngle_Stand_flag){
              if(PairMass_Tight_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_im_tight->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_im_tight_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_im_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_im_tight->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_IM_tight->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_im_tight_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_im_tight->Fill(mcXib);
                  hMCRecoLevPairXib_im_tight->Fill(XibeXi);
                }
              }
              if(PairMass_Stand_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_im_stand->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_im_stand_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_im_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_im_stand->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_IM_stand->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_im_stand_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_im_stand->Fill(mcXib);
                  hMCRecoLevPairXib_im_stand->Fill(XibeXi);
                }
              }
              if(PairMass_Loose_flag) {
                if(fabs(mcptXic0) != 9999) hRPM_im_loose->Fill(mcptXic0,Pt,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hRPM_im_loose_un->Fill(Pt,mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevXic0_im_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(fabs(mcptXic0) != 9999) hMCRecoLevPair_im_loose->Fill(Pt);
                if(fabs(mcptXic0) != 9999 &&  mcc_flag == 1) hMCRecoLevPromptXic0_IM_loose->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
                if(echarge*vcharge>0 && fabs(mcXib) != 9999){
                  hRPM_im_loose_Xib->Fill(mcXib,XibeXi);
                  hMCRecoLevXib_im_loose->Fill(mcXib);
                  hMCRecoLevPairXib_im_loose->Fill(XibeXi);
                }
              }
            }
            if(echarge*vcharge<0 && e_Recon_Stand_flag && Xi_Recon_Stand_flag && e_PID_Stand_flag && Xi_Topology_Stand_flag){
              if(mcc_flag == 1){ if(OPAngle_Stand_flag && PairMass_Stand_flag) hprompt->Fill(mcptXic0); }
              if(mcb_flag == 1){ if(OPAngle_Stand_flag && PairMass_Stand_flag) hnonprompt->Fill(mcptXic0); }
              if(OPAngle_Stand_flag && PairMass_Stand_flag) hinclu->Fill(mcptXic0,fWeightFit->Eval(mcptXic0));
            }
          } //IsMC
        } //nTreeValue

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

        if(IsMC)
        {
          TObject* hist = input->Get("histogram");

          hMCGenLevXic0_inc = (TH1F*) hist->FindObject("hTrueXic0");
          hMCGenLevXic0_inc = (TH1F*) hMCGenLevXic0_inc->Clone("hMCGenLevXic0_inc");  //generation level for 0.5 rap. , ptbinning1
          hMCGenLevXic0_incW = (TH1F*) hist->FindObject("hTrueXic0W");
          hMCGenLevXic0_incW = (TH1F*) hMCGenLevXic0_incW->Clone("hMCGenLevXic0_incW");
          hMCGenLevXic0_inc8 = (TH1F*) hist->FindObject("hTrueXic0_rap8");
          hMCGenLevXic0_inc8 = (TH1F*) hMCGenLevXic0_inc8->Clone("hMCGenLevXic0_inc8");  //generation level for 0.5 rap. , ptbinning1
          hMCGenLevXic0_incW8 = (TH1F*) hist->FindObject("hTrueXic0W_rap8");
          hMCGenLevXic0_incW8 = (TH1F*) hMCGenLevXic0_incW8->Clone("hMCGenLevXic0_incW8");
          hMCGenLevXic01_p = (TH1F*) hist->FindObject("hXic0PtFromCharm1");
          hMCGenLevXic01_p = (TH1F*) hMCGenLevXic01_p->Clone("hMCGenLevXic01_p");  //generation level for 0.5 rap. , ptbinning1
          hMCGenLevXic02_p = (TH1F*) hist->FindObject("hXic0PtFromCharm2");
          hMCGenLevXic02_p = (TH1F*) hMCGenLevXic02_p->Clone("hMCGenLevXic02_p");  //generation level for 0.5 rap. , ptbinning1
          hMCGenLevXic01_np = (TH1F*) hist->FindObject("hXic0PtFromBottom1");
          hMCGenLevXic01_np = (TH1F*) hMCGenLevXic01_np->Clone("hMCGenLevXic01_np");  //generation level for 0.5 rap. , ptbinning1
          hMCGenLevXic02_np = (TH1F*) hist->FindObject("hXic0PtFromBottom2");
          hMCGenLevXic02_np = (TH1F*) hMCGenLevXic02_np->Clone("hMCGenLevXic02_np");  //generation level for 0.5 rap. , ptbinning1
          hMCRecoLevXic0_inc = (TH1F*) hist->FindObject("hGenXic0Pt");
          hMCRecoLevXic0_inc = (TH1F*) hMCRecoLevXic0_inc->Clone("hMCRecoLevXic0_inc");  //reconstruction level, ptbinning1
          hMCRecoLevXic0_np = (TH1F*) hist->FindObject("hGenXic0PtFromXib");
          hMCRecoLevXic0_np = (TH1F*) hMCRecoLevXic0_np->Clone("hMCRecoLevXic0_np");
          hMCRecoLevXic0_p = (TH1F*) hist->FindObject("hGenXic0PtFromXic");
          hMCRecoLevXic0_p = (TH1F*) hMCRecoLevXic0_p->Clone("hMCRecoLevXic0_p");
          hMCGenLevXic_p = (TH1F*) hist->FindObject("hXic0PtFromCharmW");
          hMCGenLevXic_p = (TH1F*) hMCGenLevXic_p->Clone("hMCGenLevXic_pW");
          XibGen = (TH1F*) hist->FindObject("XibGen");
          XibGen = (TH1F*) XibGen->Clone("XibGen");
          XibGen05 = (TH1F*) hist->FindObject("XibGen05");
          XibGen05 = (TH1F*) XibGen05->Clone("XibGen05");
          eff_inc = (TH1F*) hinclu->Clone("eff_inc");
          eff_inc->Divide(eff_inc,hMCGenLevXic0_inc,1,1,"b");
          eff_p = (TH1F*) hMCGenLevXic01_p->Clone("eff_p");
          eff_p->Add(hMCGenLevXic02_p,1);
          eff_p->Divide(hprompt,eff_p,1,1,"b"); //eff_p->Scale(1.4);
          eff_np = (TH1F*) hMCGenLevXic01_np->Clone("eff_np");
          eff_np->Add(hMCGenLevXic02_np,1);
          eff_np->Divide(hnonprompt,eff_np,1,1,"b"); //eff_np->Scale(1.4);
          eff_inc_2 = (TH1F*) hinclu->Clone("eff_inc_2");
          eff_inc_2->Divide(eff_inc_2,hMCGenLevXic0_incW,1,1,"b");
          TH1F* hPromptXic0 = (TH1F*) hMCGenLevXic01_p->Clone("hMCGenLevXic_p");
          hPromptXic0->Add(hMCGenLevXic02_p,2);
          NonPromptXicRap = (TH1F*) hist->FindObject("hNonPromptXicRap");
          NonPromptXicRap = (TH1F*) NonPromptXicRap->Clone("hNonPromptXicRap");
          PromptXicRap = (TH1F*) hist->FindObject("hPromptXicRap");
          PromptXicRap = (TH1F*) PromptXicRap->Clone("hPromptXicRap");
          XicRap = (TH1F*) hist->FindObject("hXicRap");
          XicRap = (TH1F*) XicRap->Clone("hXicRap");
          XicRap2 = (TH2F*) hist->FindObject("hXicPtRap");
          XicRap2 = (TH2F*) XicRap2->Clone("hXicPtRap");
        }

        TH1F *hXic0_Bayes_stand2 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand2");
        TH1F *hXic0_Bayes_stand3 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand3");
        TH1F *hXic0_Bayes_stand4 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand4");
        TH1F *hXic0_Bayes_stand5 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand5");
        TH1F *hXic0_Bayes_stand6 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand6");
        TH1F *hXic0_Bayes_stand7 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Bayes_stand7");
        TH1F *hXic0_Svd_stand3 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand3");
        TH1F *hXic0_Svd_stand4 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand4");
        TH1F *hXic0_Svd_stand5 = (TH1F*) hMCRecoLevXic0_XiPID_stand->Clone("hXic0_Svd_stand5");

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

        TH1F* hRPM_Bayes_stand2_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand2_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand2 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand2");
        TH1F* hMCRecoLevPairXib_Bayes_stand2 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand2");
        TH1F* hRPM_Bayes_stand3_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand3_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand3 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand3");
        TH1F* hMCRecoLevPairXib_Bayes_stand3 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand3");
        TH1F* hRPM_Bayes_stand4_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand4_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand4 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand4");
        TH1F* hMCRecoLevPairXib_Bayes_stan4 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand4");
        TH1F* hRPM_Bayes_stand5_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand5_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand5 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand5");
        TH1F* hMCRecoLevPairXib_Bayes_stand5 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand5");
        TH1F* hRPM_Bayes_stand6_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand6_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand6 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand6");
        TH1F* hMCRecoLevPairXib_Bayes_stand6 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand6");
        TH1F* hRPM_Bayes_stand7_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Bayes_stand7_Xib");
        TH1F* hMCRecoLevXib_Bayes_stand7 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Bayes_stand7");
        TH1F* hMCRecoLevPairXib_Bayes_stand7 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Bayes_stand7");

        TH1F* hRPM_Svd_stand3_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Svd_stand3_Xib");
        TH1F* hMCRecoLevXib_Svd_stand3 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand3");
        TH1F* hMCRecoLevPairXib_Svd_stand3 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand3");
        TH1F* hRPM_Svd_stand4_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Svd_stand4_Xib");
        TH1F* hMCRecoLevXib_Svd_stand4 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand4");
        TH1F* hMCRecoLevPairXib_Svd_stand4 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand4");
        TH1F* hRPM_Svd_stand5_Xib = (TH1F*) hRPM_op_stand_Xib->Clone("hRPM_Svd_stand5_Xib");
        TH1F* hMCRecoLevXib_Svd_stand5 = (TH1F*) hMCRecoLevXib_op_stand->Clone("hMCRecoLevXib_Svd_stand5");
        TH1F* hMCRecoLevPairXib_Svd_stand5 = (TH1F*) hMCRecoLevPairXib_op_stand->Clone("hMCRecoLevPairXib_Svd_stand5");




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
        TH1F *hRawPt_TePID = (TH1F*) hPtRS_TePID->Clone("hRawPt_ePID_tight"); hRawPt_TePID->Add(hPtWS_TePID,-1);
        TH1F *hRawPt_SePID = (TH1F*) hPtRS_SePID->Clone("hRawPt_ePID_stand"); hRawPt_SePID->Add(hPtWS_SePID,-1);
        TH1F *hRawPt_LePID = (TH1F*) hPtRS_LePID->Clone("hRawPt_ePID_loose"); hRawPt_LePID->Add(hPtWS_LePID,-1);
        TH1F *hRawPt_VLePID = (TH1F*) hPtRS_VLePID->Clone("hRawPt_ePID_vloose"); hRawPt_VLePID->Add(hPtWS_VLePID,-1);

        TH1F *hRawPt_VTXiPID = (TH1F*) hPtRS_VTXiPID->Clone("hRawPt_XiPID_vtight"); hRawPt_VTXiPID->Add(hPtWS_VTXiPID,-1);
        TH1F *hRawPt_TXiPID = (TH1F*) hPtRS_TXiPID->Clone("hRawPt_XiPID_tight"); hRawPt_TXiPID->Add(hPtWS_TXiPID,-1);
        TH1F *hRawPt_SXiPID = (TH1F*) hPtRS_SXiPID->Clone("hRawPt_XiPID_stand"); hRawPt_SXiPID->Add(hPtWS_SXiPID,-1);
        TH1F *hRawPt_LXiPID = (TH1F*) hPtRS_LXiPID->Clone("hRawPt_XiPID_loose"); hRawPt_LXiPID->Add(hPtWS_LXiPID,-1);
        TH1F *hRawPt_VLXiPID = (TH1F*) hPtRS_VLXiPID->Clone("hRawPt_XiPID_vloose"); hRawPt_VLXiPID->Add(hPtWS_VLXiPID,-1);

        TH1F *hRawPt_VTXiPIDRec = (TH1F*) hPtRS_VTXiPIDRec->Clone("hRawPt_XiPIDRec_vtight"); hRawPt_VTXiPIDRec->Add(hPtWS_VTXiPIDRec,-1);
        TH1F *hRawPt_TXiPIDRec = (TH1F*) hPtRS_TXiPIDRec->Clone("hRawPt_XiPIDRec_tight"); hRawPt_TXiPIDRec->Add(hPtWS_TXiPIDRec,-1);
        TH1F *hRawPt_SXiPIDRec = (TH1F*) hPtRS_SXiPIDRec->Clone("hRawPt_XiPIDRec_stand"); hRawPt_SXiPIDRec->Add(hPtWS_SXiPIDRec,-1);
        TH1F *hRawPt_LXiPIDRec = (TH1F*) hPtRS_LXiPIDRec->Clone("hRawPt_XiPIDRec_loose"); hRawPt_LXiPIDRec->Add(hPtWS_LXiPIDRec,-1);
        TH1F *hRawPt_VLXiPIDRec = (TH1F*) hPtRS_VLXiPIDRec->Clone("hRawPt_XiPIDRec_vloose"); hRawPt_VLXiPIDRec->Add(hPtWS_VLXiPIDRec,-1);

        TH1F *hRawPt_VTXiPIDVLRec = (TH1F*) hPtRS_VTXiPIDVLRec->Clone("hRawPt_XiPIDinvRec_vtight"); hRawPt_VTXiPIDVLRec->Add(hPtWS_VTXiPIDVLRec,-1);
        TH1F *hRawPt_TXiPIDLRec = (TH1F*) hPtRS_TXiPIDLRec->Clone("hRawPt_XiPIDinvRec_tight"); hRawPt_TXiPIDLRec->Add(hPtWS_TXiPIDLRec,-1);
        TH1F *hRawPt_SXiPIDSRec = (TH1F*) hPtRS_SXiPIDSRec->Clone("hRawPt_XiPIDinvRec_stand"); hRawPt_SXiPIDSRec->Add(hPtWS_SXiPIDSRec,-1);
        TH1F *hRawPt_LXiPIDTRec = (TH1F*) hPtRS_LXiPIDTRec->Clone("hRawPt_XiPIDinvRec_loose"); hRawPt_LXiPIDTRec->Add(hPtWS_LXiPIDTRec,-1);
        TH1F *hRawPt_VLXiPIDVTRec = (TH1F*) hPtRS_VLXiPIDVTRec->Clone("hRawPt_XiPIDinvRec_vloose"); hRawPt_VLXiPIDVTRec->Add(hPtWS_VLXiPIDVTRec,-1);

        TH1F *hRawPt_SXiPID_V0 = (TH1F*) hPtRS_SXiPID_V0->Clone("hRawPt_XiPID_stand_V0"); hRawPt_SXiPID_V0->Add(hPtWS_SXiPID_V0,-1);
        TH1F *hRawPt_SXiPID_Xi = (TH1F*) hPtRS_SXiPID_Xi->Clone("hRawPt_XiPID_stand_Xi"); hRawPt_SXiPID_Xi->Add(hPtWS_SXiPID_Xi,-1);
        TH1F *hRawPt_SXiPID_DCAb = (TH1F*) hPtRS_SXiPID_DCAb->Clone("hRawPt_XiPID_stand_DCAb"); hRawPt_SXiPID_DCAb->Add(hPtWS_SXiPID_DCAb,-1);

        TH1F *hRawPt_VTXiRec = (TH1F*) hPtRS_VTXiRec->Clone("hRawPt_XiRec_vtight"); hRawPt_VTXiRec->Add(hPtWS_VTXiRec,-1);
        TH1F *hRawPt_TXiRec = (TH1F*) hPtRS_TXiRec->Clone("hRawPt_XiRec_tight"); hRawPt_TXiRec->Add(hPtWS_TXiRec,-1);
        TH1F *hRawPt_SXiRec = (TH1F*) hPtRS_SXiRec->Clone("hRawPt_XiRec_stand"); hRawPt_SXiRec->Add(hPtWS_SXiRec,-1);
        TH1F *hRawPt_LXiRec = (TH1F*) hPtRS_LXiRec->Clone("hRawPt_XiRec_loose"); hRawPt_LXiRec->Add(hPtWS_LXiRec,-1);
        TH1F *hRawPt_VLXiRec = (TH1F*) hPtRS_VLXiRec->Clone("hRawPt_XiRec_vloose"); hRawPt_VLXiRec->Add(hPtWS_VLXiRec,-1);

        TH1F *hRawPt_VLeRec = (TH1F*) hPtRS_VLeRec->Clone("hRawPt_eRec_vloose"); hRawPt_VLeRec->Add(hPtWS_VLeRec,-1);
        TH1F *hRawPt_LeRec = (TH1F*) hPtRS_LeRec->Clone("hRawPt_eRec_loose"); hRawPt_LeRec->Add(hPtWS_LeRec,-1);
        TH1F *hRawPt_SeRec = (TH1F*) hPtRS_SeRec->Clone("hRawPt_eRec_stand"); hRawPt_SeRec->Add(hPtWS_SeRec,-1);
        TH1F *hRawPt_TeRec = (TH1F*) hPtRS_TeRec->Clone("hRawPt_eRec_tight"); hRawPt_TeRec->Add(hPtWS_TeRec,-1);
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

        TObject* hist = input->Get("histogram");

        TH1F* NumOfEvtperRun = (TH1F*) hist->FindObject("NumOfEvtperRun");
        NumOfEvtperRun = (TH1F*) NumOfEvtperRun->Clone("NumOfEvtperRun");
        TH1F* NumOfe = (TH1F*) hist->FindObject("NumOfe");
        NumOfe = (TH1F*) NumOfe->Clone("NumOfe");
        TH1F* NumOfXi = (TH1F*) hist->FindObject("NumOfXi");
        NumOfXi = (TH1F*) NumOfXi->Clone("NumOfXi");

        TH1D *hDS = (TH1D*) hist->FindObject("DSElectronPair");
        TH1D* hSS = (TH1D*) hist->FindObject("SSElectronPair");
        hDS = (TH1D*) hDS->Clone("hDS");
        hSS = (TH1D*) hSS->Clone("hSS");

        TH1F* NumOfXiperEvt = new TH1F;
        NumOfXiperEvt = (TH1F*) NumOfXi->Clone("NumOfXiperEvt"); NumOfXiperEvt->Divide(NumOfEvtperRun);
        TH1F* NumOfeperEvt = new TH1F;
        NumOfeperEvt = (TH1F*) NumOfe->Clone("NumOfeperEvt"); NumOfeperEvt->Divide(NumOfEvtperRun);

        TH1F* hpreff_eRec_loose = (TH1F*) hpre_eRec_loose_nu->Clone("hpreff_eRec_loose"); hpreff_eRec_loose->Divide(hpreff_eRec_loose,hpre_eRec_loose_de,1,1,"b");
        TH1F* hpreff_eRec_stand = (TH1F*) hpre_eRec_stand_nu->Clone("hpreff_eRec_stand"); hpreff_eRec_stand->Divide(hpreff_eRec_stand,hpre_eRec_stand_de,1,1,"b");
        TH1F* hpreff_eRec_tight = (TH1F*) hpre_eRec_tight_nu->Clone("hpreff_eRec_tight"); hpreff_eRec_tight->Divide(hpreff_eRec_tight,hpre_eRec_tight_de,1,1,"b");
        TH1F* hpreff_XiRec_loose = (TH1F*) hpre_XiRec_loose_nu->Clone("hpreff_XiRec_loose"); hpreff_XiRec_loose->Divide(hpreff_XiRec_loose,hpre_XiRec_loose_de,1,1,"b");
        TH1F* hpreff_XiRec_stand = (TH1F*) hpre_XiRec_stand_nu->Clone("hpreff_XiRec_stand"); hpreff_XiRec_stand->Divide(hpreff_XiRec_stand,hpre_XiRec_stand_de,1,1,"b");
        TH1F* hpreff_XiRec_tight = (TH1F*) hpre_XiRec_tight_nu->Clone("hpreff_XiRec_tight"); hpreff_XiRec_tight->Divide(hpreff_XiRec_tight,hpre_XiRec_tight_de,1,1,"b");
        TH1F* hpreff_ePID_loose = (TH1F*) hpre_ePID_loose_nu->Clone("hpreff_ePID_loose"); hpreff_ePID_loose->Divide(hpreff_ePID_loose,hpre_ePID_loose_de,1,1,"b");
        TH1F* hpreff_ePID_stand = (TH1F*) hpre_ePID_stand_nu->Clone("hpreff_ePID_stand"); hpreff_ePID_stand->Divide(hpreff_ePID_stand,hpre_ePID_stand_de,1,1,"b");
        TH1F* hpreff_ePID_tight = (TH1F*) hpre_ePID_tight_nu->Clone("hpreff_ePID_tight"); hpreff_ePID_tight->Divide(hpreff_ePID_tight,hpre_ePID_tight_de,1,1,"b");
        TH1F* hpreff_XiPID_loose = (TH1F*) hpre_XiPID_loose_nu->Clone("hpreff_XiPID_loose"); hpreff_XiPID_loose->Divide(hpreff_XiPID_loose,hpre_XiPID_loose_de,1,1,"b");
        TH1F* hpreff_XiPID_stand = (TH1F*) hpre_XiPID_stand_nu->Clone("hpreff_XiPID_stand"); hpreff_XiPID_stand->Divide(hpreff_XiPID_stand,hpre_XiPID_stand_de,1,1,"b");
        TH1F* hpreff_XiPID_tight = (TH1F*) hpre_XiPID_tight_nu->Clone("hpreff_XiPID_tight"); hpreff_XiPID_tight->Divide(hpreff_XiPID_tight,hpre_XiPID_tight_de,1,1,"b");
        TH1F* hpreff_OA_loose = (TH1F*) hpre_OA_loose_nu->Clone("hpreff_OA_loose"); hpreff_OA_loose->Divide(hpreff_OA_loose,hpre_OA_loose_de,1,1,"b");
        TH1F* hpreff_OA_stand = (TH1F*) hpre_OA_stand_nu->Clone("hpreff_OA_stand"); hpreff_OA_stand->Divide(hpreff_OA_stand,hpre_OA_stand_de,1,1,"b");
        TH1F* hpreff_OA_tight = (TH1F*) hpre_OA_tight_nu->Clone("hpreff_OA_tight"); hpreff_OA_tight->Divide(hpreff_OA_tight,hpre_OA_tight_de,1,1,"b");
        TH1F* hpreff_IM_loose = (TH1F*) hpre_IM_loose_nu->Clone("hpreff_IM_loose"); hpreff_IM_loose->Divide(hpreff_IM_loose,hpre_IM_loose_de,1,1,"b");
        TH1F* hpreff_IM_stand = (TH1F*) hpre_IM_stand_nu->Clone("hpreff_IM_stand"); hpreff_IM_stand->Divide(hpreff_IM_stand,hpre_IM_stand_de,1,1,"b");
        TH1F* hpreff_IM_tight = (TH1F*) hpre_IM_tight_nu->Clone("hpreff_IM_tight"); hpreff_IM_tight->Divide(hpreff_IM_tight,hpre_IM_tight_de,1,1,"b");

        hpre_Bayes_stand2_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand2_nu"); hpre_Bayes_stand2_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand2_de");
        hpre_Bayes_stand3_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand3_nu"); hpre_Bayes_stand3_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand3_de");
        hpre_Bayes_stand4_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand4_nu"); hpre_Bayes_stand4_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand4_de");
        hpre_Bayes_stand5_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand5_nu"); hpre_Bayes_stand5_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand5_de");
        hpre_Bayes_stand6_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand6_nu"); hpre_Bayes_stand6_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand6_de");
        hpre_Bayes_stand7_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand7_nu"); hpre_Bayes_stand7_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand7_de");
        hpre_Svd_stand3_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand3_nu"); hpre_Svd_stand3_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand3_de");
        hpre_Svd_stand4_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand4_nu"); hpre_Svd_stand4_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand4_de");
        hpre_Svd_stand5_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Svd_stand5_nu"); hpre_Svd_stand5_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Svd_stand5_de");

        TH1F* hpre_Bayes_stand_nu = (TH1F*) hpre_IM_stand_nu->Clone("hpre_Bayes_stand_nu"); TH1F* hpre_Bayes_stand_de = (TH1F*) hpre_IM_stand_de->Clone("hpre_Bayes_stand_de");

        TH1F *hRawPt_Svd_stand = (TH1F*) hRawPt_SeIM->Clone("hRawPt_Bayes_stand");
        TH2F* hRPM_Unfold_stand = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand");
		    TH2F* hRPM_Unfold_stand_pt2 = (TH2F*) hRPM_im_stand->Clone("hRPM_Bayes_stand_pt2");
        TH1F* hMCRecoLevXic0_Unfold_stand = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand");
        TH1F* hMCRecoLevPair_Unfold_stand = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand");
        TH1F* hMCRecoLevXic0_Unfold_stand_pt2 = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hMCRecoLevXic0_Bayes_stand_pt2");
        TH1F* hMCRecoLevPair_Unfold_stand_pt2 = (TH1F*) hMCRecoLevPair_im_stand->Clone("hMCRecoLevPair_Bayes_stand_pt2");
        TH1F *hXic0_Svd_stand = (TH1F*) hMCRecoLevXic0_im_stand->Clone("hXic0_Bayes_stand");


        TH1D* hEventNumbers = (TH1D*) hist->FindObject("hEventNumbers"); hEventNumbers = (TH1D*) hEventNumbers->Clone("hEventNumbers");

        return;
}
