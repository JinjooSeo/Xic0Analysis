#include "DrawTool.C"
#include "GetPreFilterCorrectedSpectrum.C"

void GetUnfoldingLimit(){
  #ifdef __CINT__
    gSystem->Load("libRooUnfold");
  #endif

  Bool_t DrawOption=kTRUE;

  TFile *DataROOTFile = TFile::Open("out_data_MB_0to100.root");  //default
  TFile *MCROOTFile = TFile::Open("out_MCraw_MB_0to100.root");   //default

  TFile *DataVar1 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar1.root");  //0.1 1 2 3...
  TFile *MCVar1 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar1.root");   //0.1 1 2 3...

  TFile *DataVar2 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar2.root");  //0.3 1 2 3...
  TFile *MCVar2 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar2.root");   //0.3 1 2 3...

  TFile *DataVar3 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar3.root");  //0.5 1 2 3...
  TFile *MCVar3 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar3.root");   //0.5 1 2 3...

  TFile *DataVar4 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar4.root");  //0 0.1 1 2 3...
  TFile *MCVar4 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar4.root");   //0 0.1 1 2 3...

  TFile *DataVar5 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar5.root");  //0 0.3 1 2 3...
  TFile *MCVar5 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar5.root");   //0 0.3 1 2 3...

  TFile *DataVar6 = TFile::Open("CheckSystematics/out_data_MB_0to100_binVar6.root");  //0 0.5 1 2 3...
  TFile *MCVar6 = TFile::Open("CheckSystematics/out_MCraw_MB_0to100_binVar6.root");   //0 0.5 1 2 3...

  TH1D* hRawYield_default = GetPreFilterCorrectedSpectrum(DataROOTFile, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_default = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_default = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_refold_default = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_default = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_default(hUnweigthedeXiPair_default,hUnweigthedXic0_default,Mat_unfold_default);
  double Unfoldbinning[11] = {0,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_default = new TH1D(Form("hUnweightedReco_%s_%s","eRec","stand"),"",10,Unfoldbinning);
  RooUnfoldBayes unfolding_unweighted_default (&unfold_unweighted_default, hRawYield_default, 3); // OR
  hUnweightedReco_default = (TH1D*) unfolding_unweighted_default.Hreco();

  TH1D* hRawYield_Var1 = GetPreFilterCorrectedSpectrum(DataVar1, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var1 = (TH1D*) MCVar1->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var1 = (TH1D*) MCVar1->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var1 = (TH2D*) MCVar1->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_Var1(hUnweigthedeXiPair_Var1,hUnweigthedXic0_Var1,Mat_unfold_Var1);
  double Unfoldbinning_Var1[11] = {0.1,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var1 = new TH1D(Form("hUnweightedReco_%s_%s_Var1","eRec","stand"),"",10,Unfoldbinning_Var1);
  RooUnfoldBayes unfolding_unweighted_Var1 (&unfold_unweighted_Var1, hRawYield_Var1, 3); // OR
  hUnweightedReco_Var1 = (TH1D*) unfolding_unweighted_Var1.Hreco();

  TH1D* hRawYield_Var2 = GetPreFilterCorrectedSpectrum(DataVar2, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var2 = (TH1D*) MCVar2->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var2 = (TH1D*) MCVar2->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var2 = (TH2D*) MCVar2->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_Var2(hUnweigthedeXiPair_Var2,hUnweigthedXic0_Var2,Mat_unfold_Var2);
  double Unfoldbinning_Var2[11] = {0.3,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var2 = new TH1D(Form("hUnweightedReco_%s_%s_Var2","eRec","stand"),"",10,Unfoldbinning_Var2);
  RooUnfoldBayes unfolding_unweighted_Var2 (&unfold_unweighted_Var2, hRawYield_Var2, 3); // OR
  hUnweightedReco_Var2 = (TH1D*) unfolding_unweighted_Var2.Hreco();

  TH1D* hRawYield_Var3 = GetPreFilterCorrectedSpectrum(DataVar3, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var3 = (TH1D*) MCVar3->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var3 = (TH1D*) MCVar3->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var3 = (TH2D*) MCVar3->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_Var3(hUnweigthedeXiPair_Var3,hUnweigthedXic0_Var3,Mat_unfold_Var3);
  double Unfoldbinning_Var3[11] = {0.5,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var3 = new TH1D(Form("hUnweightedReco_%s_%s_Var3","eRec","stand"),"",10,Unfoldbinning_Var3);
  RooUnfoldBayes unfolding_unweighted_Var3 (&unfold_unweighted_Var3, hRawYield_Var3, 3); // OR
  hUnweightedReco_Var3 = (TH1D*) unfolding_unweighted_Var3.Hreco();

  TH1D* hRawYield_Var4 = GetPreFilterCorrectedSpectrum(DataVar4, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var4 = (TH1D*) MCVar4->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var4 = (TH1D*) MCVar4->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var4 = (TH2D*) MCVar4->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_Var4(hUnweigthedeXiPair_Var4,hUnweigthedXic0_Var4,Mat_unfold_Var4);
  double Unfoldbinning_Var4[12] = {0,0.1,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var4 = new TH1D(Form("hUnweightedReco_%s_%s_Var4","eRec","stand"),"",11,Unfoldbinning_Var4);
  RooUnfoldBayes unfolding_unweighted_Var4 (&unfold_unweighted_Var4, hRawYield_Var4, 3); // OR
  hUnweightedReco_Var4 = (TH1D*) unfolding_unweighted_Var4.Hreco();

  TH1D* hRawYield_Var5 = GetPreFilterCorrectedSpectrum(DataVar5, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var5 = (TH1D*) MCVar5->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var5 = (TH1D*) MCVar5->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var5 = (TH2D*) MCVar5->Get(Form("hRPM_%s_%s_un","eRec","stand"));
  RooUnfoldResponse unfold_unweighted_Var5(hUnweigthedeXiPair_Var5,hUnweigthedXic0_Var5,Mat_unfold_Var5);
  double Unfoldbinning_Var5[12] = {0,0.3,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var5 = new TH1D(Form("hUnweightedReco_%s_%s_Var5","eRec","stand"),"",11,Unfoldbinning_Var5);
  RooUnfoldBayes unfolding_unweighted_Var5 (&unfold_unweighted_Var5, hRawYield_Var5, 3); // OR
  hUnweightedReco_Var5 = (TH1D*) unfolding_unweighted_Var5.Hreco();

  TH1D* hRawYield_Var6 = GetPreFilterCorrectedSpectrum(DataVar6, "eRec", "stand", kFALSE);
  TH1D* hUnweigthedXic0_Var6 = (TH1D*) MCVar6->Get(Form("hMCRecoLevXic0_%s_%s","eRec","stand"));
  TH1D* hUnweigthedeXiPair_Var6 = (TH1D*) MCVar6->Get(Form("hMCRecoLevPair_%s_%s","eRec","stand"));
  TH2D* Mat_unfold_Var6 = (TH2D*) MCVar6->Get(Form("hRPM_%s_%s_un","ePID","stand"));
  RooUnfoldResponse unfold_unweighted_Var6(hUnweigthedeXiPair_Var6,hUnweigthedXic0_Var6,Mat_unfold_Var6);
  double Unfoldbinning_Var6[12] = {0,0.5,1,2,3,4,5,6,8,12,16,20};
  TH1D *hUnweightedReco_Var6 = new TH1D(Form("hUnweightedReco_%s_%s_Var6","eRec","stand"),"",11,Unfoldbinning_Var6);
  RooUnfoldBayes unfolding_unweighted_Var6 (&unfold_unweighted_Var6, hRawYield_Var6, 3); // OR
  hUnweightedReco_Var6 = (TH1D*) unfolding_unweighted_Var6.Hreco();

  double ptbinning[8] = {1.,2., 3., 4., 5., 6., 8., 12.};
  TH1D* hdef = (TH1D*) hUnweightedReco_default->Rebin(7,"hdef",ptbinning);
  TH1D* hvar1 = (TH1D*) hUnweightedReco_Var1->Rebin(7,"hvar1",ptbinning);
  TH1D* hvar2 = (TH1D*) hUnweightedReco_Var2->Rebin(7,"hvar2",ptbinning);
  TH1D* hvar3 = (TH1D*) hUnweightedReco_Var3->Rebin(7,"hvar3",ptbinning);
  TH1D* hvar4 = (TH1D*) hUnweightedReco_Var4->Rebin(7,"hvar4",ptbinning);
  TH1D* hvar5 = (TH1D*) hUnweightedReco_Var5->Rebin(7,"hvar5",ptbinning);
  TH1D* hvar6 = (TH1D*) hUnweightedReco_Var6->Rebin(7,"hvar6",ptbinning);

  TH1D* hvar1Ratio = (TH1D*) GetRatio(hvar1,hdef,"b");
  TH1D* hvar2Ratio = (TH1D*) GetRatio(hvar2,hdef,"b");
  TH1D* hvar3Ratio = (TH1D*) GetRatio(hvar3,hdef,"b");
  TH1D* hvar4Ratio = (TH1D*) GetRatio(hvar4,hdef,"b");
  TH1D* hvar5Ratio = (TH1D*) GetRatio(hvar5,hdef,"b");
  TH1D* hvar6Ratio = (TH1D*) GetRatio(hvar6,hdef,"b");


  if(DrawOption){
    SetStyle();
    int nCan = 2;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("XsectionCan%d",i),"",650,500);

    can[0]->cd();
    HistSty(hUnweightedReco_default,kBlack,kFullCircle); hUnweightedReco_default->Draw();
    HistSty(hUnweightedReco_Var1,kRed,kFullCircle); hUnweightedReco_Var1->Draw("SAME");
    HistSty(hUnweightedReco_Var2,kBlue,kFullCircle); hUnweightedReco_Var2->Draw("SAME");
    HistSty(hUnweightedReco_Var3,kGreen,kFullCircle); hUnweightedReco_Var3->Draw("SAME");
    HistSty(hUnweightedReco_Var4,kRed+2,kCircle); hUnweightedReco_Var4->Draw("SAME");
    HistSty(hUnweightedReco_Var5,kBlue+2,kCircle); hUnweightedReco_Var5->Draw("SAME");
    HistSty(hUnweightedReco_Var6,kGreen+2,kCircle); hUnweightedReco_Var6->Draw("SAME");

    SetAxis(hUnweightedReco_default,"#it{p}_{T} (GeV/#it{c})","Unfolded yield");
    hUnweightedReco_default->GetXaxis()->SetRangeUser(0,12);
    TLegend *leg0 = new TLegend(0.55,0.65,0.76,0.82);
    leg0->AddEntry(hUnweightedReco_default,"Default");
    leg0->AddEntry(hUnweightedReco_Var1,"Var1");
    leg0->AddEntry(hUnweightedReco_Var2,"Var2");
    leg0->AddEntry(hUnweightedReco_Var3,"Var3");
    leg0->AddEntry(hUnweightedReco_Var4,"Var4");
    leg0->AddEntry(hUnweightedReco_Var5,"Var5");
    leg0->AddEntry(hUnweightedReco_Var6,"Var6");
    LegSty(leg0);
    leg0->Draw();

    can[1]->cd();
    HistSty(hvar1Ratio,kRed,kFullCircle); hvar1Ratio->Draw("HIST"); hvar1Ratio->SetLineStyle(kSolid);
    HistSty(hvar2Ratio,kBlue,kFullCircle); hvar2Ratio->Draw("HIST SAME"); hvar1Ratio->SetLineStyle(kSolid);
    HistSty(hvar3Ratio,kGreen,kFullCircle); hvar3Ratio->Draw("HIST SAME"); hvar1Ratio->SetLineStyle(kSolid);
    HistSty(hvar4Ratio,kRed+2,kCircle); hvar4Ratio->Draw("HIST SAME"); hvar1Ratio->SetLineStyle(kDashDotted);
    HistSty(hvar5Ratio,kBlue+2,kCircle); hvar5Ratio->Draw("HIST SAME"); hvar1Ratio->SetLineStyle(kDashDotted);
    HistSty(hvar6Ratio,kGreen+2,kCircle); hvar6Ratio->Draw("HIST SAME"); hvar1Ratio->SetLineStyle(kDashDotted);
    SetAxis(hvar1Ratio,"#it{p}_{T} (GeV/#it{c})","Variation/Default");
    hvar1Ratio->GetYaxis()->SetRangeUser(0.95,1.05);
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(hvar1Ratio,"Var1");
    leg1->AddEntry(hvar2Ratio,"Var2");
    leg1->AddEntry(hvar3Ratio,"Var3");
    leg1->AddEntry(hvar4Ratio,"Var4");
    leg1->AddEntry(hvar5Ratio,"Var5");
    leg1->AddEntry(hvar6Ratio,"Var6");
    LegSty(leg1);
    leg1->Draw();

    can[0]->SaveAs("ComparisonOfUnfoldingVarBin.pdf");
    can[1]->SaveAs("RatioOfUnfoldingVarBin.pdf");
  }

  return;
}
