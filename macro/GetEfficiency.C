//#include "DrawTool.C"

TH1D* GetEfficiency(TFile* WeightedROOTFile, TFile* MCROOTFile, const Char_t* Cut, const Char_t* CutFlag, Bool_t IsWeighted, Bool_t DrawOption = kFALSE){
  TH1D* hGenXic0_woW = new TH1D;
  TH1D* hRecoXic0_woW = new TH1D;
  TH1D* hefficiency_woW = new TH1D;
  TH1D* hGenXic0_W_rap08 = new TH1D;
  TH1D* hRecoXic0_W_rap08 = new TH1D;
  TH1D* hefficiency_W_rap08 = new TH1D;
  TH1D* hGenXic0_W = new TH1D;
  TH1D* hRecoXic0_W = new TH1D;
  TH1D* hefficiency_W = new TH1D;
  Double_t bin[8] = {1.,2.,3.,4.,5.,6.,8.,12.};

  hGenXic0_W = (TH1D*) WeightedROOTFile->Get("hMCGenInclusiveXic0_W");
  hRecoXic0_W = (TH1D*) WeightedROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
  hefficiency_W = (TH1D*) hRecoXic0_W->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
  hefficiency_W->Divide(hefficiency_W,hGenXic0_W,1,1,"b");

  hGenXic0_woW = (TH1D*) MCROOTFile->Get("hMCGenInclusiveXic0_woW");
  hRecoXic0_woW = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
  hefficiency_woW = (TH1D*) hRecoXic0_woW->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
  hefficiency_woW->Divide(hefficiency_woW,hGenXic0_woW,1,1,"b");

  /*hGenXic0_W_rap08 = (TH1D*) WeightedROOTFile->Get("hMCGenInclusiveXic0_W_rap08");
  hRecoXic0_W_rap08 = (TH1D*) WeightedROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
  hefficiency_W_rap08 = (TH1D*) hRecoXic0_W_rap08->Rebin(7,Form("efficiency_%s_%s_rap08",Cut,CutFlag),bin);
  hefficiency_W_rap08->Divide(hefficiency_W_rap08,hGenXic0_W_rap08,1,1,"b");*/

  if(DrawOption){
    SetStyle();
    int nCan = 3;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("EfficiencyCan%d",i),"",650,500);

    can[0]->cd();
    hefficiency_W->Draw(); HistSty(hefficiency_W,kRed,kFullCircle);
    SetAxis(hefficiency_W,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","Efficiency");

    can[1]->cd();
    hefficiency_W->Draw();
    hefficiency_woW->Draw("SAME"); HistSty(hefficiency_woW,kBlack,kCircle);
    //hefficiency_W_rap08->Draw("SAME"); HistSty(hefficiency_W_rap08,kBlue,kCircle);
    SetAxis(hefficiency_W,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","Efficiency");
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(hefficiency_W,"weighted");
    leg1->AddEntry(hefficiency_woW,"unweighted");
    leg1->AddEntry(hefficiency_W_rap08,"weighted, |y|<0.8");
    LegSty(leg1);
    leg1->Draw();

    can[2]->cd();
    TH1D* hratio1 = GetRatio(hefficiency_W, hefficiency_woW, "b");
    //TH1D* hratio2 = GetRatio(hefficiency_W_rap08, hefficiency_W, "b");
    hratio1->Draw(); HistSty(hratio1,kBlack,kFullCircle);
    //hratio2->Draw(); HistSty(hratio2,kBlue,kFullCircle);
    SetAxis(hratio1,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","Weighted/UnWeighted");
    /*TLegend *leg2 = new TLegend(0.55,0.65,0.76,0.82);
    leg2->AddEntry(hratio1,"unweighted");
    leg2->AddEntry(hratio2,"weighted, |y|<0.8");
    LegSty(leg2);
    leg2->Draw();*/

    can[0]->SaveAs(Form("WeightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));
    can[1]->SaveAs(Form("WeightedUnweightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));
    can[2]->SaveAs(Form("RatioOfWeightedUnweightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));

    delete[] can;
  }

  if(IsWeighted) return hefficiency_W;
  if(!IsWeighted) return hefficiency_woW;
  return hefficiency_W;
}
