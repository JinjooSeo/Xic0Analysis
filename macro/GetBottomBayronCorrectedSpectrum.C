//#include "DrawTool.C"

TH1D* GetBottomBayronCorrectedSpectrum(TFile* MCROOTFile, TH1D* hCMSLb, Double_t BRFraction, const Char_t* Cut, const Char_t* CutFlag, TH1D* hMeas_tmp, Bool_t DrawOption = kFALSE){
  #ifdef __CINT__
    gSystem->Load("libRooUnfold");
  #endif

  //1) Fit 7TeV CMS Lb spectrum with Tasllis function---------------------------------------------------------------------------------
  TF1 *fTsallis = new TF1("Lambdab","[0]*x*(pow(1+(sqrt(pow(x,2)+pow(5.619,2))-5.619)/(7.6*1.1),-7.6))",0,50);
  hCMSLb->Fit("Lambdab");
  fTsallis->SetParameters(fTsallis->GetParameters());

  //2) Multiply scale factor to convert to 13TeV Lb---------------------------------------------------------------------------------
  double binning[10] = {1,2,3,4,5,6,8,12,16,20};
  TH1D *hScaleFactor = new TH1D("hScaleFactor","",9, binning);  //B meson ratio -> B(13TeV)/B(7TeV)
  Double_t BC_ScaleFactor[9] = {1.53313, 1.69604, 1.806626, 1.887637, 1.950308, 2.018669, 2.121922, 2.249672, 2.439034};
  Double_t BE_ScaleFactor[9] = {0.002688548, 0.003208365, 0.004097259, 0.005324936, 0.006890631, 0.009749647, 0.01760265, 0.03540223, 0.06652629};
  for(int i = 0; i<9; i++){
    hScaleFactor->SetBinContent(i+1,BC_ScaleFactor[i]);
    hScaleFactor->SetBinError(i+1,BE_ScaleFactor[i]);
  }

  TH1D* h13TeVLb = new TH1D("h13TeVLb","",9,binning);
  for(int i=1; i<10; i++) h13TeVLb->SetBinContent(i,fTsallis->Eval(h13TeVLb->GetBinCenter(i))*hScaleFactor->GetBinContent(i)); //Error is not assigned since we don't know the error at 1 to 10 pT

  //3) Multiply Branching ratio fraction to conver to Xib spectrum---------------------------------------------------------------------------------
  TH1D* h13TeVXib = (TH1D*) h13TeVLb->Clone("h13TeVXib");
  h13TeVXib->Scale(BRFraction);

  //4) Calculate Xib yield---------------------------------------------------------------------------------
  double ptbin[9] = {1,1,1,1,1,2,4,4,4};
  TH1D* h13TeVXibRaw = new TH1D("h13TeVXibRaw","",9,binning);
  Double_t Luminosity=1.88554e+09/(57.8*1000000);  //pp 13TeV integrated luminosity
  for(int i=1; i<10; i++) h13TeVXibRaw->SetBinContent(i,h13TeVXib->GetBinContent(i)*ptbin[i-1]*2*Luminosity);

  TH1D* hGenXib = (TH1D*) MCROOTFile->Get("XibGen05");  //Number of Xib in generation level
  TH1D* hRecoXib = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXib_%s_%s",Cut,CutFlag));  //Number of Xib in reconstruction level
  TH1D* hRecoeXi = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPairXib_%s_%s",Cut,CutFlag));  //Number of eXi from Xib
  TH2D* hRM_Xib = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s_Xib",Cut,CutFlag));  //Response matrix of Xib and eXi
  TH1D* hXibEff = (TH1D*) hRecoXib->Clone(Form("hXibEff_%s_%s",Cut,CutFlag));
  hXibEff->Divide(hXibEff,hGenXib,1,1,"b");  //Xib efficiency
  h13TeVXibRaw->Multiply(hXibEff);

  //5) Convert to eXi spectrum from Xib spectrum---------------------------------------------------------------------------------
  RooUnfoldResponse response(hRecoXib,hRecoeXi,hRM_Xib);
  RooUnfoldBinByBin unfolding(&response, h13TeVXibRaw);  //BinByBin refolding
  TH1D* heXiFromXib = (TH1D*) unfolding.Hreco();

  //6) Add eXi from Xib to eXi pair(RS-WS)---------------------------------------------------------------------------------
  TH1D* hMeas = (TH1D*) hMeas_tmp->Clone(Form("hRawPt_%s_%s_use",Cut,CutFlag));
  hMeas->Add(heXiFromXib);

  if(DrawOption){
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("BottomBayronCan%d",i),"",650,500);

    can[0]->cd();
    hRM_Xib->Draw("COLZ");
    SetAxis(hRM_Xib,"#it{p}_{T}(#Xi_{b}) (GeV/#it{c})","#it{p}_{T}(e#Xi) (GeV/#it{c})");

    can[1]->cd();
    hXibEff->Draw(); HistSty(hXibEff,kBlack,kFullCircle);
    SetAxis(hXibEff,"#it{p}_{T}(#Xi_{b}) (GeV/#it{c})","efficiency");

    can[2]->cd();
    heXiFromXib->Draw(); HistSty(heXiFromXib,kRed,kFullCircle);
    h13TeVXibRaw->Draw("SAME"); HistSty(h13TeVXibRaw,kBlack,kFullCircle);
    SetAxis(heXiFromXib,"#it{p}_{T} (GeV/#it{c})","Entries");
    TLegend *leg3 = new TLegend(0.55,0.65,0.76,0.82);
    leg3->AddEntry(heXiFromXib,"eXi");
    leg3->AddEntry(h13TeVXibRaw,"Xib");
    LegSty(leg3);
    leg3->Draw();

    can[3]->cd();
    hMeas->Draw(); HistSty(hMeas,kRed,kFullCircle);
    hMeas_tmp->Draw("SAME"); HistSty(hMeas_tmp,kBlack,kFullCircle);
    SetAxis(hMeas,"#it{p}_{T}(e#Xi) (GeV/#it{c})","Entries");
    TLegend *leg4 = new TLegend(0.55,0.65,0.76,0.82);
    leg4->AddEntry(hMeas,"corrected");
    leg4->AddEntry(hMeas_tmp,"not corrected");
    LegSty(leg4);
    leg4->Draw();

    can[4]->cd();
    TH1D* hratio = GetRatio(hMeas, hMeas_tmp, "b");
    hratio->Draw(); HistSty(hratio,kBlack,kFullCircle);
    SetAxis(hratio,"#it{p}_{T}(e#Xi) (GeV/#it{c})","corrected/not corrected");

    can[0]->SaveAs(Form("ResponseMatrixOfXib_%s_%s.pdf",Cut,CutFlag));
    can[1]->SaveAs(Form("EfficiencyOfXib_%s_%s.pdf",Cut,CutFlag));
    can[2]->SaveAs(Form("UnfoldRefoldOfXib_%s_%s.pdf",Cut,CutFlag));
    can[3]->SaveAs(Form("CorrectedUncorrectedOfXib_%s_%s.pdf",Cut,CutFlag));
    can[4]->SaveAs(Form("RatioOfCorrectedUncorrectedOfXib_%s_%s.pdf",Cut,CutFlag));

    delete[] can;
  }

  if(!DrawOption){
    delete h13TeVXibRaw;
    delete hGenXib;
    delete hRecoXib;
    delete hRecoeXi;
    delete fTsallis;
    delete hScaleFactor;
    delete h13TeVLb;
    delete h13TeVXib;
    delete hRM_Xib;
  }

  return hMeas;
}
