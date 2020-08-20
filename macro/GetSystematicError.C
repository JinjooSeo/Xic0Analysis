#include "DrawTool.C"
#include "GetPreFilterCorrectedSpectrum.C"
#include "GetBottomBayronCorrectedSpectrum.C"
#include "GetUnfoldedSpectrum.C"
#include "GetEfficiency.C"
#include "GetPromptEfficiency.C"
#include "GetCrossSection.C"
#include "GetPromptFraction.C"

TFile *DataROOTFile = new TFile("DataAnalysisHistogram.root");  //hist
TFile *MCROOTFile = TFile::Open("MCAnalysisHistogram.root");   ///////AnalysisHistogram of mc   //MChist
TFile* WeightedMCROOTFile = TFile::Open("MCWeightedAnalysisHistogram.root");

TFile* XiSystROOTFile = TFile::Open("XiCutSystematics.root");
TFile* WeightedROOTFile = TFile::Open("WeightedCS.root");
TFile* UnweightedROOTFile = TFile::Open("UnweightedCS.root");
TFile* WeightedROOTFile8 = TFile::Open("WeightedCS8.root");
TFile* PromptEff = TFile::Open("PromptEff.root");
TFile* PromptFraction = TFile::Open("PromptFraction.root");

//----------------------------Setting parameter//----------------------------//
Bool_t IsWeighted = kTRUE;
const Char_t* method = "Bayes";
Int_t NumOfIteration = 3;
TH1D* hCMSLb;
Double_t BRFraction  = (3.9*10e-4)/(5.8*10e-5);
Int_t NumberOfBin = 6;
Double_t ptBinning[7] = {2.,3.,4.,5.,6.,8.,12.};
//----------------------------Setting parameter//----------------------------//

TH1D* GetSystematicUncertaintyOfUnfolding(const Char_t* method, vector <TString> iteration, vector <Int_t> NumOfIterationUnfolding, Bool_t IsPtDependent, Bool_t DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfCuts(const Char_t* cut, vector <TString> CutFlag, Int_t iteration, Bool_t IsPtDependent, Bool_t  DrawOption = kFALSE);
TH1D* GetSystematicUncertaintyOfBottomBayron(Double_t LbUncertainty, Double_t BRUncertainty);
TH1D* GetSystematicUncertaintyFromRatio(TH1D* hnu, TH1D hde, Bool_t DrawOption = kFALSE);
TH1D* GetCMSLbSpectrum();
Double_t GetBarlowValue(TH1D* def, TH1D* var, int bin);

void GetSystematicError(){
  hCMSLb = GetCMSLbSpectrum();
  vector <TString> CutList1 = {"eRec","ePID","IM","OA","Bayes","Svd"};
  vector <TString> CutFlag1 = {"stand","vloose","loose","tight","vtight"};
  vector <TString> CutFlag2 = {"stand","loose","tight"};
  vector <TString> CutFlag3 = {"stand3","stand2","stand4","stand5","stand6","stand7"};
  vector <Int_t> NoI1 = {3,2,4,5,6,7};
  vector <TString> CutFlag4 = {"stand3","stand4","stand5"};
  vector <Int_t> NoI2 = {3,6,7};

  TH1D** hSysUn = new TH1D*[SysList.size()];

  //TH1D* htest1 = GetPreFilterCorrectedSpectrum(DataROOTFile, "eRec", "stand", kFALSE);
  //TH1D* htest2 = GetBottomBayronCorrectedSpectrum(WeightedMCROOTFile, hCMSLb, BRFraction, "eRec", "stand", htest1, kFALSE);
  //TH1D* htest3 = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, htest2, "Bayes", "stand3", "Bayes", 3, kTRUE, kFALSE);
  //TH1D* htest4 = GetEfficiency(WeightedMCROOTFile, MCROOTFile, "Bayes", "stand3", kTRUE, kFALSE);
  //TH1D* htest5 = GetCrossSection(htest3, htest4, "Bayes", "stand3", kFALSE);
  //TH1D* test = GetSystematicUncertaintyOfUnfolding("Bayes", CutFlag3, NoI1, kTRUE, kFALSE);
  //TH1D* htest6 = GetPromptFraction(htest5,htest4,kFALSE,kTRUE);

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
    hXSection[0]->GetXaxis()->SetRangeUser(2,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<CutFlag.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",CutFlag[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(2,12);
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
    heff[0]->GetXaxis()->SetRangeUser(2,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(2,12);
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
    hBarlow[1]->GetXaxis()->SetRangeUser(2,12);
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

TH1D* GetSystematicUncertaintyOfBottomBayron(Double_t LbUncertainty, Double_t BRUncertainty){
  TH1D** heff = new TH1D*[3];
  TH1D** hMeasuredSpectrum = new TH1D*[3];
  TH1D** hRawSpectrum = new TH1D*[3];
  TH1D** hXSection = new TH1D*[3];
  TH1D** hXRatio = new TH1D*[3];
  TH1D* hRMS = new TH1D(Form("RMS_%s","BottomBaryon"),"",NumberOfBin,ptBinning);
  TH1D* hRMS_average = new TH1D(Form("RMS_%s_average","BottomBaryon"),"",NumberOfBin,ptBinning);
  Double_t StdDev[NumberOfBin]; for (int i=0; i<NumberOfBin; i++) StdDev[i] = 0.0;
  Double_t sum = 0.0;

  for(int i=0; i<3; i++) {
    if(i==1) hCMSLb->Scale(1+LbUncertainty);
    if(i==2) {
      hCMSLb->Scale(1/1+LbUncertainty); //origianl Lb
      BRFraction = BRFraction+BRUncertainty;
    }
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile, "XiRec", "stand");
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction, "XiRec", "stand", hMeasuredSpectrum[i]);
    hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], "XiRec", "stand", "Bayes", NumOfIteration, IsWeighted);
    heff[i] = GetEfficiency(WeightedMCROOTFile, MCROOTFile, "XiRec", "stand", IsWeighted);
    hXSection[i] = GetCrossSection(hRawSpectrum[i], heff[i], "XiRec", "stand");
    hXRatio[i] = (TH1D*) hXSection[i]->Clone(Form("%d/%d_%s",i,0,"BottomBaryon"));
    hXRatio[i]->Divide(hXRatio[i],hXSection[0],1,1,"b");

    for (int j=0; j<NumberOfBin; j++) StdDev[j] += pow(1-fabs(hXRatio[i]->GetBinContent(j+1)),2);
  }
  BRFraction = BRFraction-BRUncertainty; //origianl BR Fraction

  for(int i=1; i<NumberOfBin+1;i++) hRMS->SetBinContent(i,sqrt(StdDev[i-1]));
  TH1D* Sys_Unc = (TH1D*) hRMS->Clone("Systematic Uncertainty");

  delete[] heff;
  delete[] hMeasuredSpectrum;
  delete[] hRawSpectrum;
  delete[] hXSection;
  delete[] hXRatio;

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

const Char_t* Bayess = "Bayes";

  for(int i=0; i<iteration.size(); i++) {
    hMeasuredSpectrum[i] = GetPreFilterCorrectedSpectrum(DataROOTFile, method, iteration[i].Data());
    hMeasuredSpectrum[i] = GetBottomBayronCorrectedSpectrum(MCROOTFile, hCMSLb, BRFraction, method, iteration[i].Data(), hMeasuredSpectrum[i]);
    //if(i == 0){
    //  hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), "Bayes", NumOfIterationUnfolding[i], IsWeighted);
  //  }
    if(method == Bayess){
    hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), method, NumOfIterationUnfolding[i], IsWeighted);
  }
  else hRawSpectrum[i] = GetUnfoldedSpectrum(WeightedMCROOTFile, MCROOTFile, hMeasuredSpectrum[i], method, iteration[i].Data(), method, NumOfIterationUnfolding[i], kFALSE);
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
    hXSection[0]->GetXaxis()->SetRangeUser(2,12);
    leg1->SetBorderSize(0); LegSty(leg1); leg1->Draw();

    can[1]->cd();
    TLegend *leg2 = new TLegend();
    for(int i=1; i<iteration.size(); i++){
      hXRatio[i]->Draw("HIST SAME");
      HistSty(hXRatio[i],HistColor[i],kCircle);
      leg2->AddEntry(hXRatio[i],Form("%s/stand",iteration[i].Data()));
    }
    SetAxis(hXRatio[1],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hXRatio[1]->GetXaxis()->SetRangeUser(2,12);
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
    heff[0]->GetXaxis()->SetRangeUser(2,12);
    heff[0]->GetYaxis()->SetRangeUser(0,0.14);
    leg3->SetBorderSize(0); LegSty(leg3); leg3->Draw();

    can[3]->cd();
    SetAxis(hRMS,"#it{p}_{T} (GeV/#it{c})","RMS");
    HistSty(hRMS,kBlack,kCircle);
    HistSty(hRMS_average,kBlue,kCircle);
    hRMS->GetXaxis()->SetRangeUser(2,12);
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
    hBarlow[1]->GetXaxis()->SetRangeUser(2,12);
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

Double_t GetBarlowValue(TH1D* def, TH1D* var, int bin){
  double error = fabs(var->GetBinContent(bin)-def->GetBinContent(bin)); // abs(Variation - Default) -> "Error"
  double errorFraction = error/var->GetBinContent(bin); // abs(Variation - Default)/Deafult  ????????
  double deltasigma = sqrt(abs(pow(var->GetBinError(bin), 2)-pow(def->GetBinError(bin), 2))); //sqrt(abs(variation_stat.error^2-default_stat.error^2)) -> "delta sigma"
  double deltasigmafraction = deltasigma/var->GetBinContent(bin);  //deltasigma/variation
  //if (errorFraction_L / deltasigmafraction_L > 1) hDivided_L->SetBinContent(bin, TMath::Abs(hDivided_L->GetBinContent(bin)-1));
  //else hDivided_L->SetBinContent(bin,0);
  return errorFraction/deltasigmafraction;
}

//TH1D* GetSystematicUncertaintyOfWeighting()
//TH1D* GetSystematicUncertaintyOfRapRange()
