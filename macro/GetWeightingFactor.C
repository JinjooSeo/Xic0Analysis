//#include "DrawTool.C"

TF1* GetWeightingFactor(TFile* MCROOTFile,TH1D* hXS_Xic0_data, Bool_t DrawOption=kFALSE){
  TH1D* hXS_Xic0_mc = (TH1D*) MCROOTFile->Get("hMCGenInclusiveXic0_woW");
  hXS_Xic0_mc->SetBinContent(7,hXS_Xic0_mc->GetBinContent(7)/4);
  hXS_Xic0_mc->SetBinContent(6,hXS_Xic0_mc->GetBinContent(6)/2);

  TF1 *fXS_fit_mc = new TF1("fXS_fit_mc","expo",1,12);
  TF1 *fXS_fit_data = new TF1("fXS_fit_data","expo",1,12);

  //normalize histograms
  hXS_Xic0_data->Scale(1/hXS_Xic0_data->Integral());
  hXS_Xic0_mc->Scale(1/hXS_Xic0_mc->Integral());

  TH1D* hWeightingFactor = (TH1D*) GetRatio(hXS_Xic0_data,hXS_Xic0_mc,"s");

  Double_t ptBin[8] = {1.,2.,3.,4.,5.,6.,8.,12.};
  TH1D* hWF_center = new TH1D("hWF_center","",7,ptBin);
  TH1D* hWF_updown = new TH1D("hWF_updown","",7,ptBin);
  TH1D* hWF_downup = new TH1D("hWF_downup","",7,ptBin);
  for(int i=1; i<8; i++){
    hWF_center->SetBinContent(i,hWeightingFactor->GetBinContent(i));
    hWF_center->SetBinError(i,hWeightingFactor->GetBinError(i));
    if(i<4){
      hWF_updown->SetBinContent(i,hWeightingFactor->GetBinContent(i)+hWeightingFactor->GetBinErrorUp(i));
      hWF_downup->SetBinContent(i,hWeightingFactor->GetBinContent(i)-hWeightingFactor->GetBinErrorLow(i));
    }
    else if(i==4){
      hWF_updown->SetBinContent(i,hWeightingFactor->GetBinContent(i));
      hWF_downup->SetBinContent(i,hWeightingFactor->GetBinContent(i));
    }
    else{
      hWF_updown->SetBinContent(i,hWeightingFactor->GetBinContent(i)-hWeightingFactor->GetBinErrorLow(i));
      hWF_downup->SetBinContent(i,hWeightingFactor->GetBinContent(i)+hWeightingFactor->GetBinErrorUp(i));
    }
  }

  TF1 *fWeightingFactor_center = new TF1("fWeightingFactor_center","expo",1,12);
  hWF_center->Fit(fWeightingFactor_center,"L 0");
  TF1 *fWeightingFactor_updown = new TF1("fWeightingFactor_updown","expo",1,12);
  hWF_updown->Fit(fWeightingFactor_updown,"L 0");
  TF1 *fWeightingFactor_downup = new TF1("fWeightingFactor_downup","expo",1,12);
  hWF_downup->Fit(fWeightingFactor_downup,"L 0");

  cout << "------------------Weighing Function Parameter------------------"<< endl;
  cout << "- Fit Function : Expo" << endl;
  cout << "- Constant : " << fWeightingFactor_center->GetParameter(0) << endl;
  cout << "- Slope : " << fWeightingFactor_center->GetParameter(1) << endl;
  cout << "---------------------------------------------------------------"<< endl;

  if(DrawOption){
    SetStyle();
    int nCan = 2;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("XsectionCan%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    HistSty(hXS_Xic0_data,kBlue,kCircle); hXS_Xic0_data->Draw();
    HistSty(hXS_Xic0_mc,kRed,kCircle); hXS_Xic0_mc->Draw("SAME");
    SetAxis(hXS_Xic0_data,"#it{p}_{T} (GeV/#it{c})","dN/Nd#it{p}_{T}");
    TLegend *leg0 = new TLegend(0.55,0.65,0.76,0.82);
    leg0->AddEntry(hXS_Xic0_data,"Data");
    leg0->AddEntry(hXS_Xic0_mc,"MC");
    LegSty(leg0);
    leg0->Draw();

    can[1]->cd();
    HistSty(hWeightingFactor,kBlack,kCircle); hWeightingFactor->Draw();
    HistSty(hWF_center,kBlack,kFullCircle); hWF_center->Draw("P SAME");
    HistSty(hWF_updown,kRed,kFullCircle); hWF_updown->Draw("P SAME");
    HistSty(hWF_downup,kBlue,kFullCircle); hWF_downup->Draw("P SAME");
    fWeightingFactor_center->Draw("SAME"); fWeightingFactor_center->SetLineColor(kBlack);
    fWeightingFactor_updown->Draw("SAME"); fWeightingFactor_updown->SetLineColor(kRed);
    fWeightingFactor_downup->Draw("SAME"); fWeightingFactor_downup->SetLineColor(kBlue);
    SetAxis(hWeightingFactor,"#it{p}_{T} (GeV/#it{c})","Weighting Factor");
    hWeightingFactor->GetYaxis()->SetRangeUser(0,3);

    can[0]->SaveAs("ComparisonofXicSpectrum_DataMC.pdf");
    can[1]->SaveAs("WeightingFactor.pdf");
  }

  return fWeightingFactor_center;
}
