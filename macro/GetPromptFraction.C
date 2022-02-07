//#include "DrawTool.C"

TH1D* GetPromptFraction(TFile* UnweightedMCROOT, TFile* WeightedMCROOT, TH1D* hXic0CrossSection, TH1D* hInclusiveEff, Bool_t GetSysError, Bool_t DrawOption = kFALSE){
  TFile* PythiaROOT;    PythiaROOT = TFile::Open("HistogramXic.root");
  TFile* FeeddownLcROOT;    FeeddownLcROOT = TFile::Open("DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020.root");
  TFile* FeeddownLcROOTFONLLcent;    FeeddownLcROOTFONLLcent = TFile::Open("DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020.root");
  TFile* FeeddownLcROOTFFcent;    FeeddownLcROOTFFcent = TFile::Open("DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020.root");
  TFile* RatioXicLcROOT;    RatioXicLcROOT = TFile::Open("Xic0toLc_pp13TeV_new.root");

  TH1D** hFeeddownLc = new TH1D*[7];
  TH1D** hFeeddownXicVarLc = new TH1D*[7];
  TH1D** hFeeddownXicVarRatio = new TH1D*[3];
  TH1D** hPromptFractionVarLc = new TH1D*[7];
  TH1D** hPromptFractionVarLc_tmp = new TH1D*[7];
  TH1D** hPromptFractionVarRatio = new TH1D*[3];
  TH1D** hPromptFractionVarRatio_tmp = new TH1D*[3];
  TH1D** hFeeddownXicVarLc_tmp = new TH1D*[7];
  TH1D** hFeeddownXicVarRatio_tmp = new TH1D*[3];
  TString NameOfFeeddownLc[7]={"hFeeddownLc","hFeeddownLcMax","hFeeddownLcMin","LhFeeddownLcFFMax","hFeeddownLcFFMin","hFeeddownLcFONLLMax","hFeeddownLcFONLLMin"};
  TString NameOfFeeddownXic1[7]={"hFeeddownXic","hFeeddownXicMax","hFeeddownXicMin","hFeeddownXicFFMax","hFeeddownXicFFMin","hFeeddownXicFONLLMax","hFeeddownXicFONLLMin"};
  TString NameOfFeeddownXic2[3]={"hFeeddownXic","hFeeddownXicRatioMax","hFeeddownXicRatioMin"};

  TH1D* hMCFeedDownXic0 = (TH1D*) PythiaROOT->Get("hXicZeroFromBSpectrumY05_px__3;1");
  TH1D* hMCPromptXic0 = (TH1D*) PythiaROOT->Get("hXicZeroSpectrumY05_px__2;1");
  TH1D* hInclusiveLc = (TH1D*) FeeddownLcROOT->Get("hLcpkpipred_central");
  hFeeddownLc[0] = (TH1D*) FeeddownLcROOT->Get("hLcpkpifromBpred_central_corr");
  hFeeddownLc[1] = (TH1D*) FeeddownLcROOT->Get("hLcpkpifromBpred_max_corr");  //hFeeddownLcFONLL+FFMax
  hFeeddownLc[2] = (TH1D*) FeeddownLcROOT->Get("hLcpkpifromBpred_min_corr");  //hFeeddownLcFONLL+FFMin
  hFeeddownLc[3] = (TH1D*) FeeddownLcROOTFONLLcent->Get("hLcpkpifromBpred_max_corr");  //hFeeddownLcFFMax
  hFeeddownLc[4] = (TH1D*) FeeddownLcROOTFONLLcent->Get("hLcpkpifromBpred_min_corr");  //hFeeddownLcFFMin
  hFeeddownLc[5] = (TH1D*) FeeddownLcROOTFFcent->Get("hLcpkpifromBpred_max_corr");  //hFeeddownLcFONLLMax
  hFeeddownLc[6] = (TH1D*) FeeddownLcROOTFFcent->Get("hLcpkpifromBpred_min_corr");  //hFeeddownLcFONLLMin
  TH1D* hRatioXicLc = (TH1D*) RatioXicLcROOT->Get("hRatio_Xic0toLc");

  TH1D* hPromptEff_de = (TH1D*) UnweightedMCROOT->Get("hMCGenPromptXic0_woW");
  TH1D* hPromptEff = (TH1D*) UnweightedMCROOT->Get("hprompt");
  hPromptEff->Divide(hPromptEff,hPromptEff_de,1,1,"b");
  TH1D* hFeeddownEff_de = (TH1D*) UnweightedMCROOT->Get("hMCGenFeeddowmXic0_woW");
  TH1D* hFeeddownEff = (TH1D*) UnweightedMCROOT->Get("hnonprompt");
  hFeeddownEff->Divide(hFeeddownEff,hFeeddownEff_de,1,1,"b");

//hFeeddownEff->Draw();
//hInclusiveEff->Draw();


  //1) Feeddown Lc spectrum--------------------------------------------------------------
  double ptbinning[8] = {1., 2., 3, 4, 5, 6, 8., 12};
  //double ptbinning2[7] = {2., 3, 4, 5, 6, 8., 12};
  hMCFeedDownXic0 = (TH1D*) hMCFeedDownXic0->Rebin(7,"hMCFeedDownXic0",ptbinning);
  hMCPromptXic0 = (TH1D*) hMCPromptXic0->Rebin(7,"hMCPromptXic0",ptbinning);
  for(int i=0; i<7; i++) {
    hFeeddownLc[i] = (TH1D*) hFeeddownLc[i]->Rebin(7,NameOfFeeddownLc[i].Data(),ptbinning);
    hFeeddownLc[i]->Scale(1e-6/(0.0628*20)); //pb to ub and divide Lc2pKpi branching ratio (PDG2020) + binning
    hFeeddownLc[i]->Scale(0.616);
    for(int j=1; j<8; j++){
      hFeeddownLc[i]->SetBinContent(j,hFeeddownLc[i]->GetBinContent(j)/(hFeeddownLc[i]->GetBinWidth(j)));
      hFeeddownLc[i]->SetBinError(j,hFeeddownLc[i]->GetBinError(j)/(hFeeddownLc[i]->GetBinWidth(j)));
    } //normalization
  }
  hMCFeedDownXic0->Scale(1e-6/0.018); //pb to ub and divide Xic branching ratio -for systematic !!NOT USE!!
  hMCPromptXic0->Scale(1e-6/0.018); //pb to ub and divide Xic branching ratio -for systematic !!NOT USE!!
  hMCPromptXic0->SetBinContent(6,hMCPromptXic0->GetBinContent(6)/2);
  hMCPromptXic0->SetBinContent(7,hMCPromptXic0->GetBinContent(7)/4);
  hMCFeedDownXic0->SetBinContent(6,hMCFeedDownXic0->GetBinContent(6)/2);
  hMCFeedDownXic0->SetBinContent(7,hMCFeedDownXic0->GetBinContent(7)/4);

  //2) Ratio of Xic and Lc--------------------------------------------------------------
  TF1 *fFitFunction = new TF1("RatioFit","pol1",1,12);
  hRatioXicLc->Fit(fFitFunction,"L 0");
  TF1 *fFitFunctionMax = new TF1("RatioFitMax","pol1",1,12); // for systematic
  fFitFunctionMax->SetParameter(0,fFitFunction->GetParameter(0)*2.0);
  fFitFunctionMax->SetParameter(1,fFitFunction->GetParameter(1)*2.0);
  TF1 *fFitFunctionMin = new TF1("RatioFitMin","pol1",1,12); // for systematic
  fFitFunctionMin->SetParameter(0,fFitFunction->GetParameter(0)*0.05);
  fFitFunctionMin->SetParameter(1,fFitFunction->GetParameter(1)*0.05);


  //3) Feeddown Xic spectrum--------------------------------------------------------------
  for(int i=0; i<7; i++) {
    hFeeddownXicVarLc[i] = (TH1D*) hFeeddownLc[i]->Clone(NameOfFeeddownXic1[i].Data());
    for(int j=1; j<8; j++){
      hFeeddownXicVarLc[i]->SetBinContent(j,hFeeddownXicVarLc[i]->GetBinContent(j)*fFitFunction->Eval(hFeeddownXicVarLc[i]->GetBinCenter(j)));
      hFeeddownXicVarLc[i]->SetBinError(j,hFeeddownXicVarLc[i]->GetBinError(j)*fFitFunction->Eval(hFeeddownXicVarLc[i]->GetBinCenter(j)));
    }
  }

  for(int i=0; i<3; i++) {
    hFeeddownXicVarRatio[i] = (TH1D*) hFeeddownLc[0]->Clone(NameOfFeeddownXic2[i].Data());
    for(int j=1; j<8; j++){
      if(i == 0){
        hFeeddownXicVarRatio[i]->SetBinContent(j,hFeeddownXicVarRatio[i]->GetBinContent(j)*fFitFunction->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
        hFeeddownXicVarRatio[i]->SetBinError(j,hFeeddownXicVarRatio[i]->GetBinError(j)*fFitFunction->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
      }
      if(i == 1){
        hFeeddownXicVarRatio[i]->SetBinContent(j,hFeeddownXicVarRatio[i]->GetBinContent(j)*fFitFunctionMax->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
        hFeeddownXicVarRatio[i]->SetBinError(j,hFeeddownXicVarRatio[i]->GetBinError(j)*fFitFunctionMax->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
      }
      if(i == 2){
        hFeeddownXicVarRatio[i]->SetBinContent(j,hFeeddownXicVarRatio[i]->GetBinContent(j)*fFitFunctionMin->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
        hFeeddownXicVarRatio[i]->SetBinError(j,hFeeddownXicVarRatio[i]->GetBinError(j)*fFitFunctionMin->Eval(hFeeddownXicVarRatio[i]->GetBinCenter(j)));
      }
    }
  }  //need to modify

  //4) Prompt Fraction of Xic--------------------------------------------------------------
  TH1D* hInclusiveEff_tmp = new TH1D("hInclusiveEff_tmp","",7,ptbinning);
  TH1D* hFeeddownEff_tmp = new TH1D("hFeeddownEff_tmp","",7,ptbinning);

  for(int i=0; i<7; i++) {
    hPromptFractionVarLc_tmp[i] = (TH1D*) hXic0CrossSection->Clone(Form("hPromptFractionVarLc_tmp%d",i));
    hPromptFractionVarLc[i] = new TH1D(Form("hPromptFractionVarLc_%d",i),"",7,ptbinning);
    hFeeddownXicVarLc_tmp[i] = new TH1D(Form("hFeeddownXicVarLc_tmp%d",i),"",7,ptbinning);
    for(int j=1; j<8; j++){
      hPromptFractionVarLc[i]->SetBinContent(j,1);
      hPromptFractionVarLc[i]->SetBinError(j,0);
      hFeeddownXicVarLc_tmp[i]->SetBinContent(j,hFeeddownXicVarLc[i]->GetBinContent(j));
      hFeeddownXicVarLc_tmp[i]->SetBinError(j,hFeeddownXicVarLc[i]->GetBinError(j));
    }
  }
  for(int i=0; i<3; i++) {
      hPromptFractionVarRatio_tmp[i] = (TH1D*) hXic0CrossSection->Clone(Form("hPromptFractionVarRatio_tmp%d",i));
      hPromptFractionVarRatio[i] = new TH1D(Form("hPromptFractionVarRatio%d",i),"",7,ptbinning);
      hFeeddownXicVarRatio_tmp[i] = new TH1D(Form("hFeeddownXicVarRatio_tmp%d",i),"",7,ptbinning);
    for(int j=1; j<8; j++){
      hPromptFractionVarRatio[i]->SetBinContent(j,1);
      hPromptFractionVarRatio[i]->SetBinError(j,0);
      hFeeddownXicVarRatio_tmp[i]->SetBinContent(j,hFeeddownXicVarRatio[i]->GetBinContent(j));
      hFeeddownXicVarRatio_tmp[i]->SetBinError(j,hFeeddownXicVarRatio[i]->GetBinError(j));
    }
  }
//hInclusiveEff->Draw();

  for(int i=1; i<8; i++){
    hInclusiveEff_tmp->SetBinContent(i,hInclusiveEff->GetBinContent(i));
    hInclusiveEff_tmp->SetBinError(i,hInclusiveEff->GetBinError(i));
    hFeeddownEff_tmp->SetBinContent(i,hFeeddownEff->GetBinContent(i));
    hFeeddownEff_tmp->SetBinError(i,hFeeddownEff->GetBinError(i));
  }

  for(int i=0; i<7; i++){
    hPromptFractionVarLc_tmp[i]->Multiply(hInclusiveEff_tmp); //XS_inclusive * e_inclusive
    hFeeddownXicVarLc_tmp[i]->Multiply(hFeeddownEff_tmp); //XS_feeddown * e_feeddown
    hPromptFractionVarLc_tmp[i]->Divide(hFeeddownXicVarLc_tmp[i],hPromptFractionVarLc_tmp[i],1,1,"b"); //N_feeddown / N_inclusive
    hPromptFractionVarLc[i]->Add(hPromptFractionVarLc_tmp[i],-1);  //1.73 or 1 NEED TO CONFIRM
  }
  for(int i=0; i<3; i++){
    hPromptFractionVarRatio_tmp[i]->Multiply(hInclusiveEff_tmp);
    hFeeddownXicVarRatio_tmp[i]->Multiply(hFeeddownEff_tmp);
    hPromptFractionVarRatio_tmp[i]->Divide(hFeeddownXicVarRatio_tmp[i],hPromptFractionVarRatio_tmp[i],1,1,"b");
    hPromptFractionVarRatio[i]->Add(hPromptFractionVarRatio_tmp[i],-1); //1.73 or 1 NEED TO CONFIRM
  }

  TH1D* hfeedDumy = (TH1D*) hMCFeedDownXic0->Clone("hfeedDumy"); hfeedDumy->Multiply(hFeeddownEff_tmp);  //7
  TH1D* hPrDumy = (TH1D*) hMCPromptXic0->Clone("hPrDumy"); hPrDumy->Multiply(hPromptEff); //7
  TH1D* hIncDumy = (TH1D*) hXic0CrossSection->Clone("hIncDumy"); hIncDumy->Multiply(hInclusiveEff_tmp); //6

  TH1D* hNbtmp = new TH1D("hNbtmp","",7,ptbinning);
  for (int i=1; i<8; i++){ hNbtmp->SetBinError(i,hfeedDumy->GetBinError(i)); hNbtmp->SetBinContent(i,hfeedDumy->GetBinContent(i)); }
  hNbtmp->Divide(hNbtmp,hIncDumy,1,1,"b");
  TH1D* hfc = (TH1D*) hPrDumy->Clone("hfc"); hfc->Add(hfeedDumy,1); hfc->Divide(hPrDumy,hfc,1,1,"b");

  TH1D* hNb = new TH1D("hNb","",7,ptbinning);
  for (int i=1; i<8; i++){ hNb->SetBinError(i,0); hNb->SetBinContent(i,1); }
  hNb->Add(hNbtmp,-1);


  TFile *file = new TFile("PromptInput_middle.root","recreate");
  hFeeddownLc[0]->Write();
  hFeeddownXicVarLc[0]->Write();
  hFeeddownEff->Write();
  hPromptFractionVarLc[0]->Write();
  hXic0CrossSection->Write();
  hInclusiveEff->Write();
  file->Write();
  file->Close();


  if(DrawOption){
    SetStyle();
    int nCan = 13;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("XsectionCan%d",i),"",650,500);

    can[0]->cd();
    gPad->SetLogy();
    hMCPromptXic0->Draw(""); HistSty(hMCPromptXic0,kBlue,kFullCircle);
    hMCFeedDownXic0->Draw("SAME"); HistSty(hMCFeedDownXic0,kRed,kFullCircle);
    SetAxis(hMCPromptXic0,"#it{p}_{T} (GeV/#it{c})","Entries");
    hMCPromptXic0->GetYaxis()->SetRangeUser(0.001,10);
    TLegend *leg0 = new TLegend();
    leg0->AddEntry(hMCPromptXic0,"Prompt Xic0");
    leg0->AddEntry(hMCFeedDownXic0,"Feeddown Xic0");
    LegSty(leg0); leg0->Draw();

    can[1]->cd();
    hRatioXicLc->Draw(); hRatioXicLc->GetYaxis()->SetRangeUser(-0.2,2.0);
    fFitFunction->Draw("SAME"); fFitFunction->SetLineColor(kRed);
    fFitFunctionMax->Draw("SAME"); fFitFunctionMax->SetLineColor(kRed);
    fFitFunctionMin->Draw("SAME"); fFitFunctionMin->SetLineColor(kRed);
    SetAxis(hRatioXicLc,"#it{p}_{T} (GeV/#it{c})","#Xi_{c}/#Lambda_{c}");

    can[2]->cd();
    gPad->SetLogy();
    Color_t HistColor[9] = {kBlack,kBlue,kRed,kBlue+2,kRed+2,kYellow+2,kGreen};
    TLegend *leg2 = new TLegend();
    for(int i=0; i<7; i++){
      hFeeddownLc[i]->Draw("SAME");
      HistSty(hFeeddownLc[i],HistColor[i],kCircle);
      leg2->AddEntry(hFeeddownLc[i],NameOfFeeddownLc[i].Data());
    }
    HistSty(hFeeddownLc[0],HistColor[0],kFullCircle);
    HistSty(hFeeddownLc[1],HistColor[1],kFullCircle);
    HistSty(hFeeddownLc[2],HistColor[2],kFullCircle);
    SetAxis(hFeeddownLc[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    LegSty(leg2); leg2->Draw();

    can[3]->cd();
    gPad->SetLogy();
    TLegend *leg3 = new TLegend();
    for(int i=0; i<7; i++){
      hFeeddownXicVarLc[i]->Draw("SAME");
      HistSty(hFeeddownXicVarLc[i],HistColor[i],kCircle);
      leg3->AddEntry(hFeeddownXicVarLc[i],NameOfFeeddownXic1[i].Data());
    }
    HistSty(hFeeddownXicVarLc[0],HistColor[0],kFullCircle);
    HistSty(hFeeddownXicVarLc[1],HistColor[1],kFullCircle);
    HistSty(hFeeddownXicVarLc[2],HistColor[2],kFullCircle);
    SetAxis(hFeeddownXicVarLc[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    LegSty(leg3); leg3->Draw();

    can[4]->cd();
    gPad->SetLogy();
    TLegend *leg4 = new TLegend();
    for(int i=0; i<3; i++){
      hFeeddownXicVarRatio[i]->Draw("SAME");
      HistSty(hFeeddownXicVarRatio[i],HistColor[i+3],kCircle);
      leg4->AddEntry(hFeeddownXicVarRatio[i],NameOfFeeddownXic2[i].Data());
    }
    HistSty(hFeeddownXicVarRatio[0],HistColor[0],kFullCircle);
    SetAxis(hFeeddownXicVarRatio[0],"#it{p}_{T} (GeV/#it{c})","d^{2}#sigma/d#it{p}_{T}dy (#mub(GeV/#it{c})^{-1})");
    LegSty(leg4); leg4->Draw();

    can[5]->cd();
    TLegend *leg5 = new TLegend();
    for(int i=0; i<7; i++){
      hPromptFractionVarLc[i]->Draw("SAME");
      HistSty(hPromptFractionVarLc[i],HistColor[i],kCircle);
      leg5->AddEntry(hPromptFractionVarLc[i],NameOfFeeddownXic1[i].Data());
    }
    hPromptFractionVarLc[0]->GetYaxis()->SetRangeUser(0.8,1.1);
    HistSty(hPromptFractionVarLc[0],HistColor[0],kFullCircle);
    SetAxis(hPromptFractionVarLc[0],"#it{p}_{T} (GeV/#it{c})","Prompt fraction");
    LegSty(leg5); leg5->Draw();

    can[6]->cd();
    TLegend *leg6 = new TLegend();
    for(int i=0; i<3; i++){
      hPromptFractionVarRatio[i]->Draw("SAME");
      HistSty(hPromptFractionVarRatio[i],HistColor[i+3],kCircle);
      leg6->AddEntry(hPromptFractionVarRatio[i],NameOfFeeddownXic2[i].Data());
    }
    hPromptFractionVarRatio[0]->GetYaxis()->SetRangeUser(0.8,1.1);
    HistSty(hPromptFractionVarRatio[0],HistColor[0],kFullCircle);
    SetAxis(hPromptFractionVarRatio[0],"#it{p}_{T} (GeV/#it{c})","Prompt fraction");
    LegSty(leg6); leg6->Draw();

    can[7]->cd();
    TLegend *leg7 = new TLegend();
    hPromptFractionVarLc[0]->Draw();
    hPromptFractionVarLc[1]->Draw("SAME");
    hPromptFractionVarLc[2]->Draw("SAME");
    hPromptFractionVarRatio[1]->Draw("SAME");
    hPromptFractionVarRatio[2]->Draw("SAME");
    leg7->AddEntry(hPromptFractionVarLc[0],"cent");
    leg7->AddEntry(hPromptFractionVarLc[1],"FONLL+FF max");
    leg7->AddEntry(hPromptFractionVarLc[2],"FONLL+FF min");
    leg7->AddEntry(hPromptFractionVarRatio[1],"Ratio max");
    leg7->AddEntry(hPromptFractionVarRatio[2],"Ratio min");
    LegSty(leg7); leg7->Draw();

    can[8]->cd();
    TH1D* hRatioLcMax = GetRatio(hPromptFractionVarLc[1], hPromptFractionVarLc[0], "b");
    TH1D* hRatioLcMin = GetRatio(hPromptFractionVarLc[2], hPromptFractionVarLc[0], "b");
    TH1D* hRatioRatioMax = GetRatio(hPromptFractionVarRatio[1], hPromptFractionVarLc[0], "b");
    TH1D* hRatioRatioMin = GetRatio(hPromptFractionVarRatio[2], hPromptFractionVarLc[0], "b");
    hRatioLcMax->Draw(); hRatioLcMin->Draw("SAME"); hRatioRatioMax->Draw("SAME"); hRatioRatioMin->Draw("SAME");
    TLegend *leg8 = new TLegend();
    leg8->AddEntry(hRatioLcMax,"FONLL+FF max");
    leg8->AddEntry(hRatioLcMin,"FONLL+FF min");
    leg8->AddEntry(hRatioRatioMax,"Ratio max");
    leg8->AddEntry(hRatioRatioMin,"Ratio min");
    LegSty(leg8); leg8->Draw();
    SetAxis(hRatioLcMax,"#it{p}_{T} (GeV/#it{c})","Ratio");
    hRatioLcMax->GetYaxis()->SetRangeUser(0.8,1.2);

    can[9]->cd();
    hPromptFractionVarLc[0]->Draw();
    hNb->Draw("SAME"); HistSty(hNb,HistColor[1],kFullCircle);
    hfc->Draw("SAME"); HistSty(hfc,HistColor[2],kFullCircle);
    TLegend *leg9 = new TLegend();
    leg9->AddEntry(hPromptFractionVarLc[0],"Scaling #Lambda_{c}");
    leg9->AddEntry(hNb,"Nb method");
    leg9->AddEntry(hfc,"fc method");
    LegSty(leg9); leg9->Draw();

    can[10]->cd();
    TH1D** hRatioVarLc = new TH1D*[6];
    TLegend *leg10 = new TLegend();
    for(int i=0; i<6; i++){
      hRatioVarLc[i] = GetRatio(hPromptFractionVarLc[i+1], hPromptFractionVarLc[0], "b");
      hRatioVarLc[i]->Draw("SAME");
      HistSty(hRatioVarLc[i],HistColor[i],kCircle);
      leg10->AddEntry(hRatioVarLc[i],NameOfFeeddownXic1[i+1].Data());
    }
    LegSty(leg10); leg10->Draw();
    SetAxis(hRatioVarLc[0],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hRatioVarLc[0]->GetYaxis()->SetRangeUser(0.8,1.2);

    can[11]->cd();
    TH1D** hRatioVarLc2 = new TH1D*[2];
    TLegend *leg11 = new TLegend();
    for(int i=0; i<2; i++){
      hRatioVarLc2[i] = GetRatio(hPromptFractionVarRatio[i+1], hPromptFractionVarLc[0], "b");
      hRatioVarLc2[i]->Draw("SAME");
      HistSty(hRatioVarLc2[i],HistColor[i+3],kCircle);
      leg11->AddEntry(hRatioVarLc2[i],NameOfFeeddownXic2[i+1].Data());
    }
    LegSty(leg11); leg11->Draw();
    SetAxis(hRatioVarLc2[0],"#it{p}_{T} (GeV/#it{c})","Ratio");
    hRatioVarLc2[0]->GetYaxis()->SetRangeUser(0.8,1.2);

    can[12]->cd();
    TH1D* hRatioNb = GetRatio(hNb, hPromptFractionVarLc[0], "b");
    TH1D* hRatiofc = GetRatio(hfc, hPromptFractionVarLc[0], "b");
    hRatioLcMax->Draw(); hRatioLcMin->Draw("SAME"); hRatioRatioMax->Draw("SAME"); hRatioRatioMin->Draw("SAME"); hRatiofc->Draw("SAME"); hRatioNb->Draw("SAME");
    TLegend *leg12 = new TLegend();
    leg12->AddEntry(hRatioLcMax,"FONLL+FF max");
    leg12->AddEntry(hRatioLcMin,"FONLL+FF min");
    leg12->AddEntry(hRatioRatioMax,"Ratio max");
    leg12->AddEntry(hRatioRatioMin,"Ratio min");
    leg12->AddEntry(hRatioNb,"Nb method");
    leg12->AddEntry(hRatiofc,"fc method");
    LegSty(leg12); leg12->Draw();
    SetAxis(hRatioLcMax,"#it{p}_{T} (GeV/#it{c})","Ratio");
    hRatioLcMax->GetYaxis()->SetRangeUser(0.8,1.2);

    can[0]->SaveAs("Xic0GeneratedByPythia8.pdf");
    can[1]->SaveAs("RatioOfXicAndLc.pdf");
    can[2]->SaveAs("FeedDownLc.pdf");
    can[3]->SaveAs("FeedDownXicVarLc.pdf");
    can[4]->SaveAs("FeedDownXicVarRatio.pdf");
    can[5]->SaveAs("PromptFractionVarLc.pdf");
    can[6]->SaveAs("PromptFractionVarRatio.pdf");
    can[7]->SaveAs("PromptFractionForSystematic.pdf");
    can[8]->SaveAs("RatioOfPromptFractionForSystematic.pdf");
    can[9]->SaveAs("PromptFractionOfNbFcScaling.pdf");
    can[10]->SaveAs("RatioOfPromptFractionVarLc.pdf");
    can[11]->SaveAs("RatioOfPromptFractionVarRatio.pdf");
    can[12]->SaveAs("RatioOfPromptFractionVarRatioSys.pdf");
  }


  return hPromptFractionVarLc[0];
}
