#include "DrawTool.C"
#include "GetPreFilterCorrectedSpectrum.C"
#include "GetBottomBayronCorrectedSpectrum.C"
#include "GetUnfoldedSpectrum.C"
#include "GetEfficiency.C"
#include "GetPromptEfficiency.C"
#include "GetCrossSection.C"
#include "GetPromptFraction.C"
#include "GetWeightingFactor.C"
#include "GetITSTPCMatchingEfficiency.C"

TFile *DataROOTFile = TFile::Open("out_data_MB_0to100.root");  //hist
TFile *MCROOTFile = TFile::Open("out_MCraw_MB_0to100.root");   ///////AnalysisHistogram of mc   //MChist
TFile* WeightedMCROOTFile = TFile::Open("out_MCwgt_MB_0to100.root");
TFile* WeightedMCROOTFileVar1 = TFile::Open("out_MCwgt_MB_0to100_var1.root");
TFile* WeightedMCROOTFileVar2 = TFile::Open("out_MCwgt_MB_0to100_var2.root");

//----------------------------Setting parameter//----------------------------//
Bool_t IsWeighted = kTRUE;
const Char_t* method = "Bayes";
Int_t NumOfIteration = 3;
TH1D* hCMSLb;
Double_t BRFraction  = (3.9*10e-4)/(5.8*10e-5);
Int_t NumberOfBin = 7;
Double_t ptBinning[8] = {1.,2.,3.,4.,5.,6.,8.,12.};
//----------------------------Setting parameter//----------------------------//

TH2D* hevseXi = (TH2D*) MCROOTFile->Get("hevseXi");
TH2D* heXivsXic0 = (TH2D*) MCROOTFile->Get("heXivsXic0");

TH1D* GetSystematicUncertaintyOfUnfolding(const Char_t* method, vector <TString> iteration, vector <Int_t> NumOfIterationUnfolding, Bool_t IsPtDependent, Bool_t DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfCuts(const Char_t* cut, vector <TString> CutFlag, Int_t iteration, Bool_t IsPtDependent, Bool_t  DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfBottomBayron(Double_t LbUncertainty, Double_t BRUncertainty, Bool_t DrawOption);
TH1D* GetSystematicUncertaintyOfWeighting(Double_t LbUncertainty, Double_t BRUncertainty);
//TH1D* GetSystematicUncertaintyFromRatio(TH1D* hnu, TH1D hde, Bool_t DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfWeighting(TFile *MCROOTFileVar1, TFile *MCROOTFileVar2, Bool_t DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfITSTPCMatching(TH2D* hevseXi, TH2D* heXivsXic0);
TH1D* GetCMSLbSpectrum();
Double_t GetBarlowValue(TH1D* def, TH1D* var, int bin);

void GetSystematicError(){
  hCMSLb = GetCMSLbSpectrum();
  vector <TString> CutList1 = {"eRec","ePID","XiRec","XiPID","OA","IM","Bayes","Svd"};
  vector <TString> CutFlag1 = {"stand","vloose","loose","tight","vtight"};
  vector <TString> CutFlag2 = {"stand","loose","tight"};
  vector <TString> CutFlag5 = {"stand","tight"};
  vector <TString> CutFlag3 = {"stand3","stand2","stand4","stand5","stand6","stand7"};
  vector <Int_t> NoI1 = {3,2,4,5,6,7};
  vector <TString> CutFlag4 = {"stand3","stand4","stand5"};
  vector <Int_t> NoI2 = {3,4,5};

  TH1D** hSysUn = new TH1D*[16];

 /* TCanvas* can = new TCanvas;
  TH1D* hOpeningAngle_d1 = (TH1D*) DataROOTFile->Get("oa_data1");
TH1D* hOpeningAngle_d2 = (TH1D*) DataROOTFile->Get("oa_data2");
TH1D* hOpeningAngle_d = (TH1D*) hOpeningAngle_d1->Clone("hOpeningAngle_d");
hOpeningAngle_d->Add(hOpeningAngle_d2,-1);
		hOpeningAngle_d->Scale(1/hOpeningAngle_d->Integral());
hOpeningAngle_d1->Scale(1/hOpeningAngle_d1->Integral());
hOpeningAngle_d2->Scale(1/hOpeningAngle_d2->Integral());
  TH1D* hOpeningAngle_m = (TH1D*) MCROOTFile->Get("oa_mc"); hOpeningAngle_m->Scale(1/hOpeningAngle_m->Integral());
  TH1D* hOpeningAngle_md = (TH1D*) MCROOTFile->Get("oa_data"); hOpeningAngle_md->Scale(1/hOpeningAngle_md->Integral());
  hOpeningAngle_d->Draw(); HistSty(hOpeningAngle_d,kBlack,kDot);
  hOpeningAngle_m->Draw("SAME"); HistSty(hOpeningAngle_m,kRed,kDot);
  hOpeningAngle_md->Draw("SAME"); HistSty(hOpeningAngle_md,kBlue,kDot);
  hOpeningAngle_d1->Draw("SAME"); HistSty(hOpeningAngle_d1,kBlack,kFullCircle);
		hOpeningAngle_d2->Draw("SAME"); HistSty(hOpeningAngle_d2,kBlack,kCircle);
  cout << "Data -> O to 90 : " << hOpeningAngle_d->Integral(0,25) << endl;
  cout << "Data -> 0 to 70 : " << hOpeningAngle_d->Integral(0,20) << endl;
  cout << "Data -> 0 to 110 : " << hOpeningAngle_d->Integral(0,30) << endl;
  cout << "MC -> 0 to 90 : " << hOpeningAngle_m->Integral(0,25) << endl;
  cout << "MC -> 0 to 70 : " << hOpeningAngle_m->Integral(0,20) << endl;
  cout << "MC -> 0 to 110 : " << hOpeningAngle_m->Integral(0,30) << endl;*/
/*  TFile* pF = TFile::Open("OA_stand.root");
  TH1D* hstand = (TH1D*) pF->Get("25.935484319.000000_ratio");
  TFile* pE = TFile::Open("OA_tight.root");
  TH1D* htight = (TH1D*) pE->Get("25.93548441.000000_ratio");
  TH1D* hRatio = (TH1D*) GetRatio(htight,hstand,"b");
  hRatio->Draw();
  SetAxis(hRatio,"#it{p}_{T} GeV/#it{c})","Ratio");
*/



  TH1D* htest1 = GetPreFilterCorrectedSpectrum(DataROOTFile, "OA", "stand", kFALSE);//
  TH1D* htest2 = GetBottomBayronCorrectedSpectrum(WeightedMCROOTFile, hCMSLb, BRFraction, "Svd", "stand3", htest1, kFALSE);
  TH1D* htest3 = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, htest2, "Svd", "stand3", "Bayes", 3, kTRUE, kFALSE);
  TH1D* htest4 = GetEfficiency(WeightedMCROOTFile, MCROOTFile, "Svd", "stand3", kTRUE, kFALSE);
  TH1D* htest5 = GetCrossSection(htest3, htest4, "Svd", "stand3", kFALSE);
  //GetWeightingFactor(MCROOTFile,htest5,kFALSE);
  TH1D* htest6 = GetPromptFraction(MCROOTFile,WeightedMCROOTFile,htest5,htest4,kFALSE,kTRUE);
/*
  TFile *file = new TFile("To_Tiantian.root","recreate");

  TH1D* hunf = (TH1D*) htest3->Clone("WeightedUnfoldedYield");
  TH1D* hefff = (TH1D*) htest4->Clone("WeightedInclusiveEfficiency");
  //hunf->Write();
  //hefff->Write();

  file->Write();
  file->Close();
*/
/*
   hSysUn[0] = GetSystematicUncertaintyOfCuts("eRec",CutFlag1,NumOfIteration,kFALSE,kFALSE);
   hSysUn[1] = GetSystematicUncertaintyOfCuts("ePID",CutFlag1,NumOfIteration,kTRUE,kFALSE);
   hSysUn[2] = GetSystematicUncertaintyOfCuts("XiRec",CutFlag1,NumOfIteration,kTRUE,kFALSE);
   hSysUn[3] = GetSystematicUncertaintyOfCuts("XiPID",CutFlag1,NumOfIteration,kTRUE,kFALSE);
   hSysUn[4] = GetSystematicUncertaintyOfCuts("OA",CutFlag5,NumOfIteration,kTRUE,kFALSE);  //discuss
   hSysUn[5] = GetSystematicUncertaintyOfCuts("IM",CutFlag2,NumOfIteration,kTRUE,kFALSE);
   hSysUn[6] = GetSystematicUncertaintyOfUnfolding("Bayes",CutFlag3,NoI1,kFALSE,kFALSE);
   hSysUn[7] = GetSystematicUncertaintyOfUnfolding("Svd",CutFlag4,NoI2,kTRUE,kFALSE);
   hSysUn[8] = GetSystematicUncertaintyOfWeighting(WeightedMCROOTFileVar1,WeightedMCROOTFileVar2,kFALSE);
   hSysUn[9] = GetSystematicUncertaintyOfBottomBayron(0.50,0.33,kFALSE);
   hSysUn[10] = GetSystematicUncertaintyOfITSTPCMatching(hevseXi,heXivsXic0);
*/
/*
 Double_t sRapidity[7] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
 Double_t sWeighting[7] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.03};
 Double_t sPromptFractionU_Lcvar[7] = {0.0087, 0.0102, 0.0099, 0.0106, 0.0130, 0.0144, 0.0119};
 Double_t sPromptFractionD_Lcvar[7] = {0.0125, 0.0140, 0.0132, 0.0142, 0.0179, 0.0214, 0.0174};
 Double_t sPromptFractionU_Ratiovar[7] = {0.0170, 0.0220, 0.0238, 0.0278, 0.0368, 0.0450, 0.0404};
 Double_t sPromptFractionD_Ratiovar[7] = {0.0162, 0.0210, 0.0223, 0.0264, 0.0349, 0.0427, 0.0384};
 hSysUn[8] = new TH1D("hWeighting","",7,ptBinning);
 hSysUn[11] = new TH1D("hRapidity","",7,ptBinning);
 hSysUn[12] = new TH1D("hPromptFractionUp_Lcvar","",7,ptBinning);
 hSysUn[13] = new TH1D("hPromptFractionDown_Lcvar","",7,ptBinning);
 hSysUn[14] = new TH1D("hPromptFractionUp_Ratiovar","",7,ptBinning);
 hSysUn[15] = new TH1D("hPromptFractionDown_Ratiovar","",7,ptBinning);
for(int i=1; i<8; i++) {
  hSysUn[8]->SetBinContent(i,sWeighting[i-1]);
  hSysUn[11]->SetBinContent(i,sRapidity[i-1]);
  hSysUn[12]->SetBinContent(i,sPromptFractionU_Lcvar[i-1]);
  hSysUn[13]->SetBinContent(i,sPromptFractionD_Lcvar[i-1]);
  hSysUn[14]->SetBinContent(i,sPromptFractionU_Ratiovar[i-1]);
  hSysUn[15]->SetBinContent(i,sPromptFractionD_Ratiovar[i-1]);
}

  TFile* pF = TFile::Open("PromptInput_middle.root");
  TH1D* hPromptFraction = (TH1D*) pF->Get("hPromptFractionVarLc_0");
  TFile* pE = TFile::Open("PromptEfficiency.root");
  TH1D* hPromptEfficiency = (TH1D*) pE->Get("efficiency_Svd_stand3;1");
  TH1D* hInclusiveEfficiency = (TH1D*) pF->Get("efficiency_Svd_stand3;1");
  TH1D* hPromptXic0_tmp = (TH1D*) pF->Get("hXic0CrossSection_Svd_stand3;1");
  TH1D* hPromptXic0 = (TH1D*) hPromptXic0_tmp->Clone("PromptCrossSection");
  hPromptXic0->Multiply(hInclusiveEfficiency);
  hPromptXic0->Divide(hPromptXic0,hPromptEfficiency);
  hPromptXic0->Multiply(hPromptFraction);


 vector <TString> SysList = {"eReconstruction","ePID","XiRec","XiPID","OpeningAngle","InvariantMassOfPair","UnfoldingOfBayes","UnfoldingOfSVD",
 "Weighting","BottomBaryon","ITSTPCMatching","RapidityRange","PromptFractionUp_Lcvar","PromptFractionDown_Lcvar","PromptFractionUp_Ratiovar","PromptFractionDown_Ratiovar"};

 TH1D* hTotalSysUnc = new TH1D("Total Systematic Uncertainty","",NumberOfBin,ptBinning);
for(int i=1; i<8; i++){
 Double_t StdDev = sqrt(pow(hSysUn[0]->GetBinContent(i),2)+pow(hSysUn[1]->GetBinContent(i),2)+pow(hSysUn[2]->GetBinContent(i),2)+pow(hSysUn[3]->GetBinContent(i),2)
 +pow(hSysUn[4]->GetBinContent(i),2)+pow(hSysUn[6]->GetBinContent(i),2)+pow(hSysUn[5]->GetBinContent(i),2)+pow(hSysUn[7]->GetBinContent(i),2)
 +pow(hSysUn[8]->GetBinContent(i),2)+pow(hSysUn[9]->GetBinContent(i),2)+pow(hSysUn[10]->GetBinContent(i),2)+pow(hSysUn[11]->GetBinContent(i),2));
 //+pow(hSysUn[13]->GetBinContent(i),2));
 hTotalSysUnc->SetBinContent(i,StdDev);
}

TBox *h13TeVBox[7];
TCanvas* can = new TCanvas;
SetStyle();
hPromptXic0->Draw(); HistSty(hPromptXic0,kRed,kFullCircle);
SetAxis(hPromptXic0,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");

for(int i=1; i<8; i++){
  h13TeVBox[i] = new TBox(hPromptXic0->GetBinCenter(i)-hTotalSysUnc->GetBinWidth(i)/4,
  hPromptXic0->GetBinContent(i)-hPromptXic0->GetBinContent(i)*hTotalSysUnc->GetBinContent(i),
  hPromptXic0->GetBinCenter(i)+hTotalSysUnc->GetBinWidth(i)/4,
  hPromptXic0->GetBinContent(i)+hPromptXic0->GetBinContent(i)*hTotalSysUnc->GetBinContent(i));

  h13TeVBox[i]->SetFillStyle(0);
  h13TeVBox[i]->SetLineColor(kRed);
  h13TeVBox[i]->Draw("SAME");
}



TFile *file = new TFile("Xic0SemileptonicResults_New.root","recreate");
for(int i = 0; i<16; i++) {
  hSysUn[i]->SetTitle(Form("%s",SysList[i].Data()));
  hSysUn[i]->Write();
}
hPromptXic0->Write();
file->Write();
file->Close();

*/
/*  SetStyle();
  TCanvas* can = new TCanvas;
    TH1D* hTotalSysUnc = new TH1D("Total Systematic Uncertainty","",NumberOfBin,ptBinning);
  TLegend *leg_sys = new TLegend(0.55,0.65,0.76,0.82);
for(int i=1; i<8; i++){
    Double_t StdDev = sqrt(pow(hSysUn[0]->GetBinContent(i),2)+pow(hSysUn[1]->GetBinContent(i),2)+pow(hSysUn[2]->GetBinContent(i),2)+pow(hSysUn[3]->GetBinContent(i),2)
    +pow(hSysUn[4]->GetBinContent(i),2)+pow(hSysUn[6]->GetBinContent(i),2)+pow(hSysUn[5]->GetBinContent(i),2)+pow(hSysUn[7]->GetBinContent(i),2)) ;

//    Double_t StdDev = sqrt(pow(hSysUn[6]->GetBinContent(i),2)+pow(hSysUn[5]->GetBinContent(i),2)) ;

   cout << "i : "<< i << " StdDev : " <<StdDev << endl;

    hTotalSysUnc->SetBinContent(i,StdDev);
  }
        Color_t HistColorSet[9] = {kBlack,kRed,kGreen,kYellow,kViolet,kBlue,kYellow+2,kViolet+2,kBlue+2};

        for (int i = 0; i<8; i++) {
        hSysUn[i]->Draw("HIST SAME");
        HistSty(hSysUn[i],HistColorSet[i],kDot);
        leg_sys->AddEntry(hSysUn[i],CutList1[i].Data());
}
        hTotalSysUnc->Draw("HIST SAME");
        HistSty(hTotalSysUnc,HistColorSet[8],kDot);
        leg_sys->AddEntry(hTotalSysUnc,"TotalSystematicError");
        LegSty(leg_sys);
        leg_sys->Draw();
  hSysUn[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hSysUn[0]->GetYaxis()->SetTitle("Systematic Uncertainty");
  hSysUn[0]->GetXaxis()->SetRangeUser(1,12);
  hSysUn[0]->GetYaxis()->SetRangeUser(0,0.7);

*/

  /*TH1D* hpromptEff = (TH1D*) MCROOTFile->Get("eff_p");
  TH1D* hpromptXS = (TH1D*) htest5->Clone("PromptXic0");
  hpromptXS->Multiply(htest4);
  hpromptXS->Divide(hpromptEff);
  hpromptXS->Multiply(htest6);
  hpromptXS->Draw(); HistSty(hpromptXS,kRed,kCircle);
  htest5->Draw("SAME"); HistSty(htest5,kBlack,kCircle);*/

/*  TFile* mf = TFile::Open("MergedCrossSection_Xic0_pp13TeV.root");
  TH1D* hXSold = (TH1D*) mf->Get("hSigma_stat");
  TFile* old = TFile::Open("XSold.root");
  TH1D* hXSold2 = (TH1D*) old->Get("hXic0CrossSection_eRec_stand");
  htest5->Draw(""); HistSty(htest5,kRed,kCircle);
  hXSold->Draw("SAME"); HistSty(hXSold,kBlue,kCircle);
  hXSold2->Draw("SAME"); HistSty(hXSold2,kBlack,kCircle);

  SetAxis(htest5,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
  TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
  leg1->AddEntry(htest5,"new");
  leg1->AddEntry(hXSold2,"old");
  leg1->AddEntry(hXSold,"merged");
  LegSty(leg1);
  leg1->Draw();*/

  /*TFile* rap08 = TFile::Open("wXSrap08.root");
  TH1D* hXSrap08 = (TH1D*) rap08->Get("hXic0CrossSection_Svd_stand3");
  TFile* rap05 = TFile::Open("wXSrap05.root");
  TH1D* hXSrap05 = (TH1D*) rap05->Get("hXic0CrossSection_Svd_stand3");
  TH1D* hRatio = (TH1D*) GetRatio(hXSrap05,hXSrap08,"b");
  hRatio->Draw();*/
  /*hXSrap05->Draw(""); HistSty(hXSrap05,kRed,kCircle);
  hXSrap08->Draw("SAME"); HistSty(hXSrap08,kBlue,kCircle);

    SetAxis(htest5,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(hXSrap05,"|y|<0.5");
    leg1->AddEntry(hXSrap08,"|y|<0.8");
    LegSty(leg1);
    leg1->Draw();*/

/*  TFile *file = new TFile("Xic0Semileptonic.root","update");
		//file->Delete("PromptXic0");
		//file->Delete("");
		//file->Delete();
 SetAxis(hpromptXS,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(hpromptXS,"prompt");
    leg1->AddEntry(htest5,"inclusive");
    LegSty(leg1);
    leg1->Draw();
//hpromptXS->Write();
//htest5->Write();

file->Write();
file->Close();*/
  return;
}

TH1D* GetSystematicUncertaintyOfCuts(const Char_t* cut, vector <TString> CutFlag, Int_t iteration, Bool_t IsPtDependent, Bool_t  DrawOption = kFALSE){
  TH1D** heff = new TH1D*[CutFlag.size()];
  TH1D** hMeasuredSpectrum = new TH1D*[CutFlag.size()];
  TH1D** hRawSpectrum = new TH1D*[CutFlag.size()];
  TH1D** hXSection = new TH1D*[CutFlag.size()];
  TH1D** hXRatio = new TH1D*[CutFlag.size()];
  TH1D* hRMS = new TH1D(Form("RMS_%s",cut),"",NumberOfBin,ptBinning);
  TH1D* hRMS_average = new TH1D(Form("RMS_%s_average",cut),"",NumberOfBin,ptBinning);
  Double_t StdDev[NumberOfBin]; for (int i=0; i<NumberOfBin; i++) StdDev[i] = 0.0;
  Double_t sum = 0.0;

  for(int i=0; i<CutFlag.size(); i++) {
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile, cut, CutFlag[i].Data());
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction, cut, CutFlag[i].Data(), hMeasuredSpectrum[i]);
    hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], cut, CutFlag[i].Data(), method, NumOfIteration, IsWeighted);
    heff[i] = GetEfficiency(WeightedMCROOTFile, MCROOTFile, cut, CutFlag[i].Data(), IsWeighted);
    hXSection[i] = GetCrossSection(hRawSpectrum[i], heff[i], cut, CutFlag[i].Data());
    hXRatio[i] = (TH1D*) hXSection[i]->Clone(Form("%s/%s_%s",CutFlag[i].Data(),CutFlag[0].Data(),cut));
    hXRatio[i]->Divide(hXRatio[i],hXSection[0],1,1,"b");

    for (int j=0; j<NumberOfBin; j++) StdDev[j] += pow(1-fabs(hXRatio[i]->GetBinContent(j+1)),2)/(CutFlag.size()-1);
  }
  for(int i=1; i<NumberOfBin+1;i++) {
    hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
    sum += sqrt(StdDev[i-1]);
  }
  for(int i=1; i<NumberOfBin+1; i++) hRMS_average->SetBinContent(i,sum/NumberOfBin);

  TH1D* Sys_Unc = new TH1D;
  if(IsPtDependent) Sys_Unc = (TH1D*) hRMS->Clone("Systematic Uncertainty");
  else Sys_Unc = (TH1D*) hRMS_average->Clone("Systematic Uncertainty");

  if(DrawOption){
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("can%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    Color_t HistColor[9] = {kBlack,kGreen,kBlue,kRed,kViolet,kGreen+2,kBlue+2,kRed+2,kViolet+2};
    TLegend *leg1 = new TLegend();
    for(int i=0; i<CutFlag.size(); i++){
      hXSection[i]->Draw("SAME");
      HistSty(hXSection[i],HistColor[i],kCircle);
      leg1->AddEntry(hXSection[i],CutFlag[i].Data());
    }
    SetAxis(hXSection[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    hXSection[0]->GetXaxis()->SetRangeUser(1,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<CutFlag.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",CutFlag[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(1,12);
    hXRatio[1]->GetYaxis()->SetRangeUser(0.5,1.5);
    leg2->SetBorderSize(0); LegSty(leg2); leg2->Draw();

    can[2]->cd();
    TLegend *leg3 = new TLegend();
    for(int i=0; i<CutFlag.size(); i++){
      heff[i]->Draw("SAME");
      HistSty(heff[i],HistColor[i],kCircle);
      leg3->AddEntry(heff[i],CutFlag[i].Data());
    }
    SetAxis(heff[0],"#it{p}_{T} (GeV/#it{c})","Efficiency*Acceptance");
    heff[0]->GetXaxis()->SetRangeUser(1,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(1,12);
    hRMS->GetYaxis()->SetRangeUser(0,0.5);
    hRMS->Draw();
    hRMS_average->SetLineStyle(10);
    hRMS_average->SetLineWidth(3);
    hRMS_average->Draw("SAME");

    can[4]->cd();
    TLegend *leg4 = new TLegend();
    TH1D** hBarlow = new TH1D*[CutFlag.size()];
    for(int i=0; i<CutFlag.size(); i++) hBarlow[i] = (TH1D*) hXRatio[i]->Clone();
    for (int bin = 1; bin < hBarlow[1]->GetNbinsX()+1; bin++) {
      for(int j=1; j<CutFlag.size(); j++)
      hBarlow[j]->SetBinContent(bin,GetBarlowValue(hXRatio[0],hXRatio[j],bin));
    }
    for(int i=1; i<CutFlag.size(); i++){
      hBarlow[i]->Draw("HIST SAME");
      HistSty(hBarlow[i],HistColor[i],kCircle);
      leg4->AddEntry(hBarlow[i],CutFlag[i].Data());
    }
    SetAxis(hBarlow[1],"#it{p}_{T} (GeV/#it{c})","B");
    hBarlow[1]->GetXaxis()->SetRangeUser(1,12);
    hBarlow[1]->GetYaxis()->SetRangeUser(-1,4);
    leg4->SetBorderSize(0); LegSty(leg4); leg4->Draw();

    can[0]->SaveAs(Form("CrossSection_%s.pdf",cut));
    can[1]->SaveAs(Form("CrossSectionRatio_%s.pdf",cut));
    can[2]->SaveAs(Form("Efficiency_%s.pdf",cut));
    can[3]->SaveAs(Form("RMS_%s.pdf",cut));
    can[4]->SaveAs(Form("BarlowCheck_%s.pdf",cut));
  }

  if(!DrawOption){
    delete[] heff;
    delete[] hMeasuredSpectrum;
    delete[] hRawSpectrum;
    delete[] hXSection;
    delete[] hXRatio;
  }

  return Sys_Unc;
}

TH1D* GetSystematicUncertaintyOfBottomBayron(Double_t LbUncertainty, Double_t BRUncertainty, Bool_t DrawOption = kFALSE){
  TH1D** heff = new TH1D*[3];
  TH1D** hMeasuredSpectrum = new TH1D*[3];
  TH1D** hRawSpectrum = new TH1D*[3];
  TH1D** hXSection = new TH1D*[3];
  TH1D** hXRatio = new TH1D*[3];
  TH1D* hRMS = new TH1D(Form("RMS_%s","BottomBaryon"),"",NumberOfBin,ptBinning);
  TH1D* hRMS_average = new TH1D(Form("RMS_%s_average","BottomBaryon"),"",NumberOfBin,ptBinning);
  Double_t StdDev[NumberOfBin]; for (int i=0; i<NumberOfBin; i++) StdDev[i] = 0.0;
  Double_t sum = 0.0;

  vector <TString> Variation = {"default","variation1","variation2"};

  for(int i=0; i<3; i++) {
    if(i==1) hCMSLb->Scale(1+LbUncertainty);
    if(i==2) {
      hCMSLb->Scale(1/1+LbUncertainty); //origianl Lb
      BRFraction = BRFraction+BRUncertainty;
    }
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile, "BottomBaryon", Variation[i].Data());
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction, "BottomBaryon", Variation[i].Data(), hMeasuredSpectrum[i]);
    hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], "BottomBaryon", Variation[i].Data(), method, NumOfIteration, IsWeighted);
    heff[i] = GetEfficiency(WeightedMCROOTFile, MCROOTFile, "BottomBaryon", Variation[i].Data(), IsWeighted);
    hXSection[i] = GetCrossSection(hRawSpectrum[i], heff[i], "BottomBaryon", Variation[i].Data());
    hXRatio[i] = (TH1D*) hXSection[i]->Clone(Form("%d/%d_%s",i,0,"BottomBaryon"));
    hXRatio[i]->Divide(hXRatio[i],hXSection[0],1,1,"b");

    for (int j=0; j<NumberOfBin; j++) StdDev[j] += pow(1-fabs(hXRatio[i]->GetBinContent(j+1)),2);
  }

  for(int i=1; i<NumberOfBin+1;i++) {
    hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
    sum += sqrt(StdDev[i-1]);
  }
  for(int i=1; i<NumberOfBin+1; i++) hRMS_average->SetBinContent(i,sum/NumberOfBin);

  for(int i=1; i<NumberOfBin+1;i++) hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
  TH1D* Sys_Unc = (TH1D*) hRMS->Clone("Systematic Uncertainty");

  if(DrawOption){
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("can%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    Color_t HistColor[9] = {kBlack,kGreen,kBlue,kRed,kViolet,kGreen+2,kBlue+2,kRed+2,kViolet+2};
    TLegend *leg1 = new TLegend();
    for(int i=0; i<Variation.size(); i++){
      hXSection[i]->Draw("SAME");
      HistSty(hXSection[i],HistColor[i],kCircle);
      leg1->AddEntry(hXSection[i],Variation[i].Data());
    }
    SetAxis(hXSection[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    hXSection[0]->GetXaxis()->SetRangeUser(1,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<Variation.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",Variation[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(1,12);
    hXRatio[1]->GetYaxis()->SetRangeUser(0.5,1.5);
    leg2->SetBorderSize(0); LegSty(leg2); leg2->Draw();

    can[2]->cd();
    TLegend *leg3 = new TLegend();
    for(int i=0; i<Variation.size(); i++){
      heff[i]->Draw("SAME");
      HistSty(heff[i],HistColor[i],kCircle);
      leg3->AddEntry(heff[i],Variation[i].Data());
    }
    SetAxis(heff[0],"#it{p}_{T} (GeV/#it{c})","Efficiency*Acceptance");
    heff[0]->GetXaxis()->SetRangeUser(1,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(1,12);
    hRMS->GetYaxis()->SetRangeUser(0,0.5);
    hRMS->Draw();
    hRMS_average->SetLineStyle(10);
    hRMS_average->SetLineWidth(3);
    hRMS_average->Draw("SAME");

    can[4]->cd();
    TLegend *leg4 = new TLegend();
    TH1D** hBarlow = new TH1D*[Variation.size()];
    for(int i=0; i<Variation.size(); i++) hBarlow[i] = (TH1D*) hXRatio[i]->Clone();
    for (int bin = 1; bin < hBarlow[1]->GetNbinsX()+1; bin++) {
      for(int j=1; j<Variation.size(); j++)
      hBarlow[j]->SetBinContent(bin,GetBarlowValue(hXRatio[0],hXRatio[j],bin));
    }
    for(int i=1; i<Variation.size(); i++){
      hBarlow[i]->Draw("HIST SAME");
      HistSty(hBarlow[i],HistColor[i],kCircle);
      leg4->AddEntry(hBarlow[i],Variation[i].Data());
    }
    SetAxis(hBarlow[1],"#it{p}_{T} (GeV/#it{c})","B");
    hBarlow[1]->GetXaxis()->SetRangeUser(1,12);
    hBarlow[1]->GetYaxis()->SetRangeUser(-1,4);
    leg4->SetBorderSize(0); LegSty(leg4); leg4->Draw();

    can[0]->SaveAs(Form("CrossSection_%s.pdf","BottomBaryon"));
    can[1]->SaveAs(Form("CrossSectionRatio_%s.pdf","BottomBaryon"));
    can[2]->SaveAs(Form("Efficiency_%s.pdf","BottomBaryon"));
    can[3]->SaveAs(Form("RMS_%s.pdf","BottomBaryon"));
    can[4]->SaveAs(Form("BarlowCheck_%s.pdf","BottomBaryon"));
  }
  if(!DrawOption){
    delete[] heff;
    delete[] hMeasuredSpectrum;
    delete[] hRawSpectrum;
    delete[] hXSection;
    delete[] hXRatio;
  }

  return Sys_Unc;
}

TH1D* GetCMSLbSpectrum(){
  TFile* ROOTCMS;
  ROOTCMS = TFile::Open("HEPData-ins1113442-v1-Table_2.root");
  TDirectoryFile* dCMSLamb = (TDirectoryFile*) ROOTCMS->Get("Table 2");
  TH1D* hCMSLamb_tmp = (TH1D*) dCMSLamb->Get("Hist1D_y1");
  TH1D* hCMSLambe_tmp = (TH1D*) dCMSLamb->Get("Hist1D_y1_e1");
  double cmsbin[7] = {10,13,15,18,22,28,50};
  hCMSLamb_tmp->Scale(0.001/4);
  hCMSLambe_tmp->Scale(0.001/4);
  hCMSLb = new TH1D("lambdab","",6,cmsbin);
  for (int i=1; i<7; i++){
    hCMSLb->SetBinContent(i,hCMSLamb_tmp->GetBinContent(i));
    hCMSLb->SetBinError(i,hCMSLambe_tmp->GetBinContent(i));
  }

  delete hCMSLamb_tmp;
  delete hCMSLambe_tmp;

  return hCMSLb;
}

TH1D* GetSystematicUncertaintyOfWeighting(TFile *MCROOTFileVar1, TFile *MCROOTFileVar2, Bool_t DrawOption = kFALSE){
  TH1D** heff = new TH1D*[3];
  TH1D** hMeasuredSpectrum = new TH1D*[3];
  TH1D** hRawSpectrum = new TH1D*[3];
  TH1D** hXSection = new TH1D*[3];
  TH1D** hXRatio = new TH1D*[3];
  TH1D* hRMS = new TH1D(Form("RMS_%s","Weighting"),"",NumberOfBin,ptBinning);
  TH1D* hRMS_average = new TH1D(Form("RMS_%s_average","Weighting"),"",NumberOfBin,ptBinning);
  Double_t StdDev[NumberOfBin]; for (int i=0; i<NumberOfBin; i++) StdDev[i] = 0.0;
  Double_t sum = 0.0;

  vector <TString> Variation = {"default","variation1","variation2"};

  for(int i=0; i<3; i++) {
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile,"Weighting", Variation[i].Data());
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction,"Weighting",Variation[i].Data(), hMeasuredSpectrum[i]);
    if(i==0){
      hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i],"Weighting",Variation[i].Data(),method, NumOfIteration, IsWeighted);
      heff[i] = GetEfficiency(WeightedMCROOTFile, MCROOTFile, "Weighting",Variation[i].Data(), IsWeighted);
    }
    if(i==1){
      hRawSpectrum[i] = GetUnfoldedSpectrum(MCROOTFileVar1, MCROOTFile, hMeasuredSpectrum[i],"Weighting",Variation[i].Data(), method, NumOfIteration, IsWeighted);
      heff[i] = GetEfficiency(MCROOTFileVar1, MCROOTFile,"Weighting",Variation[i].Data(), IsWeighted);
    }
    if(i==2){
      hRawSpectrum[i] = GetUnfoldedSpectrum(MCROOTFileVar2, MCROOTFile, hMeasuredSpectrum[i],"Weighting",Variation[i].Data(), method, NumOfIteration, IsWeighted);
      heff[i] = GetEfficiency(MCROOTFileVar2, MCROOTFile,"Weighting",Variation[i].Data(), IsWeighted);
    }
    hXSection[i] = GetCrossSection(hRawSpectrum[i], heff[i],"Weighting",Variation[i].Data());
    hXRatio[i] = (TH1D*) hXSection[i]->Clone(Form("%d/%d_%s",i,0,"Weighting"));
    hXRatio[i]->Divide(hXRatio[i],hXSection[0],1,1,"b");

    for (int j=0; j<NumberOfBin; j++) StdDev[j] += pow(1-fabs(hXRatio[i]->GetBinContent(j+1)),2);
  }

  for(int i=1; i<NumberOfBin+1;i++) {
    hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
    sum += sqrt(StdDev[i-1]);
  }
  for(int i=1; i<NumberOfBin+1; i++) hRMS_average->SetBinContent(i,sum/NumberOfBin);

  for(int i=1; i<NumberOfBin+1;i++) hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
  TH1D* Sys_Unc = (TH1D*) hRMS->Clone("Systematic Uncertainty");

  if(DrawOption){
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("can%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    Color_t HistColor[9] = {kBlack,kGreen,kBlue,kRed,kViolet,kGreen+2,kBlue+2,kRed+2,kViolet+2};
    TLegend *leg1 = new TLegend();
    for(int i=0; i<Variation.size(); i++){
      hXSection[i]->Draw("SAME");
      HistSty(hXSection[i],HistColor[i],kCircle);
      leg1->AddEntry(hXSection[i],Variation[i].Data());
    }
    SetAxis(hXSection[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    hXSection[0]->GetXaxis()->SetRangeUser(1,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<Variation.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",Variation[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(1,12);
    hXRatio[1]->GetYaxis()->SetRangeUser(0.5,1.5);
    leg2->SetBorderSize(0); LegSty(leg2); leg2->Draw();

    can[2]->cd();
    TLegend *leg3 = new TLegend();
    for(int i=0; i<Variation.size(); i++){
      heff[i]->Draw("SAME");
      HistSty(heff[i],HistColor[i],kCircle);
      leg3->AddEntry(heff[i],Variation[i].Data());
    }
    SetAxis(heff[0],"#it{p}_{T} (GeV/#it{c})","Efficiency*Acceptance");
    heff[0]->GetXaxis()->SetRangeUser(1,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(1,12);
    hRMS->GetYaxis()->SetRangeUser(0,0.5);
    hRMS->Draw();
    hRMS_average->SetLineStyle(10);
    hRMS_average->SetLineWidth(3);
    hRMS_average->Draw("SAME");

    can[4]->cd();
    TLegend *leg4 = new TLegend();
    TH1D** hBarlow = new TH1D*[Variation.size()];
    for(int i=0; i<Variation.size(); i++) hBarlow[i] = (TH1D*) hXRatio[i]->Clone();
    for (int bin = 1; bin < hBarlow[1]->GetNbinsX()+1; bin++) {
      for(int j=1; j<Variation.size(); j++)
      hBarlow[j]->SetBinContent(bin,GetBarlowValue(hXRatio[0],hXRatio[j],bin));
    }
    for(int i=1; i<Variation.size(); i++){
      hBarlow[i]->Draw("HIST SAME");
      HistSty(hBarlow[i],HistColor[i],kCircle);
      leg4->AddEntry(hBarlow[i],Variation[i].Data());
    }
    SetAxis(hBarlow[1],"#it{p}_{T} (GeV/#it{c})","B");
    hBarlow[1]->GetXaxis()->SetRangeUser(1,12);
    hBarlow[1]->GetYaxis()->SetRangeUser(-1,4);
    leg4->SetBorderSize(0); LegSty(leg4); leg4->Draw();

    can[0]->SaveAs(Form("CrossSection_%s.pdf","Weighting"));
    can[1]->SaveAs(Form("CrossSectionRatio_%s.pdf","Weighting"));
    can[2]->SaveAs(Form("Efficiency_%s.pdf","Weighting"));
    can[3]->SaveAs(Form("RMS_%s.pdf","Weighting"));
    can[4]->SaveAs(Form("BarlowCheck_%s.pdf","Weighting"));
  }

  if(!DrawOption){
    delete[] heff;
    delete[] hMeasuredSpectrum;
    delete[] hRawSpectrum;
    delete[] hXSection;
    delete[] hXRatio;
  }

  return Sys_Unc;
}


TH1D* GetSystematicUncertaintyOfUnfolding(const Char_t* method, vector <TString> iteration, vector <Int_t> NumOfIterationUnfolding, Bool_t IsPtDependent, Bool_t DrawOption = kFALSE){
  TH1D** heff = new TH1D*[iteration.size()];
  TH1D** hMeasuredSpectrum = new TH1D*[iteration.size()];
  TH1D** hRawSpectrum = new TH1D*[iteration.size()];
  TH1D** hXSection = new TH1D*[iteration.size()];
  TH1D** hXRatio = new TH1D*[iteration.size()];
  TH1D* hRMS = new TH1D(Form("RMS_%s",method),"",NumberOfBin,ptBinning);
  TH1D* hRMS_average = new TH1D(Form("RMS_%s_average",method),"",NumberOfBin,ptBinning);
  Double_t StdDev[NumberOfBin]; for (int i=0; i<NumberOfBin; i++) StdDev[i] = 0.0;
  Double_t sum = 0.0;

const Char_t* Bayesian = "Bayes";

  for(int i=0; i<iteration.size(); i++) {
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile, method, iteration[i].Data());
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction, method, iteration[i].Data(), hMeasuredSpectrum[i]);
    if(i == 0){
      hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), "Bayes", NumOfIterationUnfolding[i], IsWeighted);
    }
    else {
      if(method == Bayesian){
        hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), method, NumOfIterationUnfolding[i], IsWeighted);
      }
      else hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), method, NumOfIterationUnfolding[i], IsWeighted);
    }

    heff[i] = GetEfficiency(WeightedMCROOTFile, MCROOTFile, method, iteration[i].Data(), IsWeighted);
    hXSection[i] = GetCrossSection(hRawSpectrum[i], heff[i], method, iteration[i].Data());
    hXRatio[i] = (TH1D*) hXSection[i]->Clone(Form("%s/%s_%s",iteration[i].Data(),iteration[0].Data(),method));
    hXRatio[i]->Divide(hXRatio[i],hXSection[0],1,1,"b");

    for (int j=0; j<NumberOfBin; j++) StdDev[j] += pow(1-fabs(hXRatio[i]->GetBinContent(j+1)),2)/(iteration.size()-1);
  }

  for(int i=1; i<NumberOfBin+1;i++) {
    hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
    sum += sqrt(StdDev[i-1]);
  }
  for(int i=1; i<NumberOfBin+1; i++) hRMS_average->SetBinContent(i,sum/NumberOfBin);

  TH1D* Sys_Unc = new TH1D;
  if(IsPtDependent) Sys_Unc = (TH1D*) hRMS->Clone("Systematic Uncertainty");
  else Sys_Unc = (TH1D*) hRMS_average->Clone("Systematic Uncertainty");

  if(DrawOption){
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("can%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    Color_t HistColor[9] = {kBlack,kGreen,kBlue,kRed,kViolet,kGreen+2,kBlue+2,kRed+2,kViolet+2};
    TLegend *leg1 = new TLegend();
    for(int i=0; i<iteration.size(); i++){
      hXSection[i]->Draw("SAME");
      HistSty(hXSection[i],HistColor[i],kCircle);
      leg1->AddEntry(hXSection[i],iteration[i].Data());
    }
    SetAxis(hXSection[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    hXSection[0]->GetXaxis()->SetRangeUser(1,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<iteration.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",iteration[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(1,12);
    hXRatio[1]->GetYaxis()->SetRangeUser(0.5,1.5);
    leg2->SetBorderSize(0); LegSty(leg2); leg2->Draw();

    can[2]->cd();
    TLegend *leg3 = new TLegend();
    for(int i=0; i<iteration.size(); i++){
      heff[i]->Draw("SAME");
      HistSty(heff[i],HistColor[i],kCircle);
      leg3->AddEntry(heff[i],iteration[i].Data());
    }
    SetAxis(heff[0],"#it{p}_{T} (GeV/#it{c})","Efficiency*Acceptance");
    heff[0]->GetXaxis()->SetRangeUser(1,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(1,12);
    hRMS->GetYaxis()->SetRangeUser(0,0.5);
    hRMS->Draw();
    hRMS_average->SetLineStyle(10);
    hRMS_average->SetLineWidth(3);
    hRMS_average->Draw("SAME");

    can[4]->cd();
    TLegend *leg4 = new TLegend();
    TH1D** hBarlow = new TH1D*[iteration.size()];
    for(int i=0; i<iteration.size(); i++) hBarlow[i] = (TH1D*) hXRatio[i]->Clone();
    for (int bin = 1; bin < hBarlow[1]->GetNbinsX()+1; bin++) {
      for(int j=1; j<iteration.size(); j++)
      hBarlow[j]->SetBinContent(bin,GetBarlowValue(hXRatio[0],hXRatio[j],bin));
    }
    for(int i=1; i<iteration.size(); i++){
      hBarlow[i]->Draw("HIST SAME");
      HistSty(hBarlow[i],HistColor[i],kCircle);
      leg4->AddEntry(hBarlow[i],iteration[i].Data());
    }
    SetAxis(hBarlow[1],"#it{p}_{T} (GeV/#it{c})","B");
    hBarlow[1]->GetXaxis()->SetRangeUser(1,12);
    hBarlow[1]->GetYaxis()->SetRangeUser(-1,4);
    leg4->SetBorderSize(0); LegSty(leg4); leg4->Draw();

    can[0]->SaveAs(Form("CrossSection_%s.pdf",method));
    can[1]->SaveAs(Form("CrossSectionRatio_%s.pdf",method));
    can[2]->SaveAs(Form("Efficiency_%s.pdf",method));
    can[3]->SaveAs(Form("RMS_%s.pdf",method));
    can[4]->SaveAs(Form("BarlowCheck_%s.pdf",method));
  }

  if(!DrawOption){
    delete[] heff;
    delete[] hMeasuredSpectrum;
    delete[] hRawSpectrum;
    delete[] hXSection;
    delete[] hXRatio;
  }

  return Sys_Unc;
}

TH1D* GetSystematicUncertaintyOfITSTPCMatching(TH2D* hevseXi, TH2D* heXivsXic0){
  Double_t SysUnc_DPG[8] = {0.011, 0.018, 0.027, 0.027, 0.023, 0.023, 0.024, 0.03};
  Double_t sys1[8] ={0,0,0,0,0,0,0,0};
  Double_t sys2[8] ={0,0,0,0,0,0,0,0};
  double sum1 = 0.0;
  double sum2 = 0.0;

  for (int i=1; i<9; i++){
    sum1 = 0.0;
    for(int j=1; j<9; j++){
      if(!hevseXi->GetBinContent(i,j)) continue;
      sys1[i-1] += hevseXi->GetBinContent(i,j)*SysUnc_DPG[j-1];
      sum1 += hevseXi->GetBinContent(i,j);
    }
    sys1[i-1] = sys1[i-1]/sum1;
cout << sys1[i-1] << endl;  
}
cout << "--" << endl;

  for (int i=1; i<9; i++){
    sum2 = 0.0;
    for(int j=1; j<9; j++){
      if(!heXivsXic0->GetBinContent(i,j)) continue;
      sys2[i-1] += heXivsXic0->GetBinContent(i,j)*sys1[j-1];
      sum2 += heXivsXic0->GetBinContent(i,j);
    }
    sys2[i-1] = sys2[i-1]/sum2;
		cout << sys2[i-1] << endl; 
  }

  TH1D* hSysUnc_ITSTPC = new TH1D("hSysUnc_ITSTPC","",7,ptBinning);
  for(int i=1; i<8; i++){
    hSysUnc_ITSTPC->SetBinContent(i,sys2[i]);
  }

  return hSysUnc_ITSTPC;
}

Double_t GetBarlowValue(TH1D* def, TH1D* var, int bin){
  double error = fabs(var->GetBinContent(bin)-def->GetBinContent(bin)); // abs(Variation - Default) -> "Error"
  double errorFraction = error/var->GetBinContent(bin); // abs(Variation - Default)/Deafult  ????????
  double deltasigma = sqrt(abs(pow(var->GetBinError(bin), 2)-pow(def->GetBinError(bin), 2))); //sqrt(abs(variation_stat.error^2-default_stat.error^2)) -> "delta sigma"
  double deltasigmafraction = deltasigma/var->GetBinContent(bin);  //deltasigma/variation
  //if (errorFraction_L / deltasigmafraction_L > 1) hDivided_L->SetBinContent(bin, TMath::Abs(hDivided_L->GetBinContent(bin)-1));
  //else hDivided_L->SetBinContent(bin,0);
  if (deltasigma == 0) return 0;
  return errorFraction/deltasigmafraction;
}

//TH1D* GetSystematicUncertaintyOfRapRange()
