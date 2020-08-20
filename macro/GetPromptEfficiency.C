//#include "DrawTool.C"

TH1D* GetPromptEfficiency(TFile* WeightedROOTFile, TFile* MCROOTFile, const Char_t* Cut, const Char_t* CutFlag, Bool_t IsWeighted, Bool_t DrawOption = kFALSE){
  TH1D* hGenXic0 = new TH1D;
  TH1D* hRecoXic0 = new TH1D;
  TH1D* hefficiency = new TH1D;
  Double_t bin[8] = {1.,2.,3.,4.,5.,6.,8.,12.};

  if(IsWeighted){
    hGenXic0 = (TH1D*) MCROOTFile->Get("hMCGenLevXic_p");
    hRecoXic0 = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPromptXic0_%s_%s",Cut,CutFlag));
    hefficiency = (TH1D*) hRecoXic0->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
    hefficiency->Divide(hefficiency,hGenXic0,1,1,"b");
  }
  else{
    hGenXic0 = (TH1D*) MCROOTFile->Get("hMCGenLevXic_p");
    hRecoXic0 = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPromptXic0_%s_%s",Cut,CutFlag));
    hefficiency = (TH1D*) hRecoXic0->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
    hefficiency->Divide(hefficiency,hGenXic0,1,1,"b");
  }

  if(DrawOption){
    SetStyle();
    int nCan = 3;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("EfficiencyCan%d",i),"",650,500);

    TH1D* heff1 = new TH1D;
    TH1D* heff2 = new TH1D;
    TH1D* hGenXic02 = new TH1D;
    TH1D* hRecoXic02 = new TH1D;

    if(IsWeighted){
      heff1 = (TH1D*) hefficiency->Clone(Form("efficiency_%s_%s_w",Cut,CutFlag));
      hGenXic02 = (TH1D*) MCROOTFile->Get("hMCGenLevXic0_inc");
      hRecoXic02 = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
      heff2 = (TH1D*) hRecoXic02->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
      heff2->Divide(heff2,hGenXic02,1,1,"b");
    }
    if(!IsWeighted){
      heff2 = (TH1D*) hefficiency->Clone(Form("efficiency_%s_%s_wo",Cut,CutFlag));
      hGenXic02 = (TH1D*) WeightedROOTFile->Get("hMCGenLevXic0_incW");
      hRecoXic02 = (TH1D*) WeightedROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
      heff1 = (TH1D*) hRecoXic02->Rebin(7,Form("efficiency_%s_%s",Cut,CutFlag),bin);
      heff1->Divide(heff1,hGenXic02,1,1,"b");
    }

    can[0]->cd();
    heff1->Draw(); HistSty(heff1,kRed,kFullCircle);
    SetAxis(heff1,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","Efficiency");

    can[1]->cd();
    heff1->Draw();
    heff2->Draw("SAME"); HistSty(heff2,kBlack,kCircle);
    SetAxis(heff1,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","Efficiency");
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(heff1,"weighted");
    leg1->AddEntry(heff2,"unweighted");
    LegSty(leg1);
    leg1->Draw();

    can[2]->cd();
    TH1D* hratio = GetRatio(heff1, heff2, "b");
    hratio->Draw(); HistSty(hratio,kBlack,kFullCircle);
    SetAxis(hratio,"#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})","weighted/unweighted");


    can[0]->SaveAs(Form("PromptWeightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));
    can[1]->SaveAs(Form("PromptWeightedUnweightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));
    can[2]->SaveAs(Form("PromptRatioOfWeightedUnweightedEfficiencyOfXic_%s_%s.pdf",Cut,CutFlag));

    delete[] can;
  }

  delete hRecoXic0;

  return hefficiency;
}
