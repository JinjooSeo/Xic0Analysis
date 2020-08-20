//#include "DrawTool.C"

TH1D* GetPreFilterCorrectedSpectrum(TFile* DataROOTFile, const Char_t* Cut, const Char_t* CutFlag, Bool_t DrawOption = kFALSE){
  TH1D* hMeas_tmp = (TH1D*) DataROOTFile->Get(Form("hRawPt_%s_%s",Cut,CutFlag));
  TH1D* hMeas = (TH1D*) hMeas_tmp->Clone(Form("hRawPt_%s_%s_use",Cut,CutFlag));
  TH1D* prefilter_eff_nu = (TH1D*) DataROOTFile->Get(Form("hpre_%s_%s_nu",Cut,CutFlag));  //prefilter
  TH1D* prefilter_eff_de = (TH1D*) DataROOTFile->Get(Form("hpre_%s_%s_de",Cut,CutFlag));  //prefilter
  TH1D* prefilter_eff = (TH1D*) prefilter_eff_nu->Clone(Form("hpreff_%s_%s",Cut,CutFlag));  //prefilter
  prefilter_eff->Divide(prefilter_eff,prefilter_eff_de,1,1,"b");
  hMeas->Divide(prefilter_eff);

  if(DrawOption){  //----------------------------------------------Draw histogram on canvas
    SetStyle();
    int nCan = 2;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("PrefilterCan%d",i),"",650,500);

    can[0]->cd();
    TH1D* hMeas_origin = (TH1D*) hMeas_tmp->Clone(Form("hMeas_origin_%s_%s",Cut,CutFlag));
    hMeas_origin->Draw(); HistSty(hMeas_origin,kBlack,kFullCircle);
    hMeas->Draw("SAME"); HistSty(hMeas,kRed,kFullCircle);
    SetAxis(hMeas_origin,"#it{p}_{T}(e#Xi) (GeV/#it{c})","Entries");
    TLegend *leg0 = new TLegend(0.55,0.65,0.76,0.82);
    leg0->AddEntry(hMeas_origin,"not corrected");
    leg0->AddEntry(hMeas,"corrected");
    LegSty(leg0);
    leg0->Draw();

    can[1]->cd();
    prefilter_eff->Draw(); HistSty(prefilter_eff,kBlack,kFullCircle);
    prefilter_eff->GetYaxis()->SetRangeUser(0.9,1.1);
    SetAxis(prefilter_eff,"#it{p}_{T}(e#Xi) (GeV/#it{c})","PreFilter Efficiency");

    can[0]->SaveAs(Form("PreFilterCorrectedSpectrum_%s_%s.pdf",Cut,CutFlag));
    can[1]->SaveAs(Form("PreFilterEfficiency_%s_%s.pdf",Cut,CutFlag));

    delete[] can;
  }

    delete hMeas_tmp;
    delete prefilter_eff_nu;
    delete prefilter_eff_de;

    return hMeas;
}
