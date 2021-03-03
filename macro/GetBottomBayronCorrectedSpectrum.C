//#include "DrawTool.C"

TH1D* ReadFONLL(TString filename, Int_t nPtBins, Double_t* binning);
TH1D* GetBottomBayronCorrectedSpectrum(TFile* MCROOTFile, TH1D* hCMSLb, Double_t BRFraction, const Char_t* Cut, const Char_t* CutFlag, TH1D* hMeas_tmp, Bool_t DrawOption = kFALSE){
  #ifdef __CINT__
    gSystem->Load("libRooUnfold");
  #endif

  TString dFONLLB13TeV = "../input/FONLL-Bmeson-dsdpt-sqrts13000-20GeV.txt";
  TString dFONLLB5TeV = "../input/FONLL-Bmeson-dsdpt-sqrts7000-20GeV.txt";

  //1) Fit 7TeV CMS Lb spectrum with Tasllis function---------------------------------------------------------------------------------
  TF1 *fTsallis = new TF1("Lambdab","[0]*x*(pow(1+(sqrt(pow(x,2)+pow(5.619,2))-5.619)/(7.6*1.1),-7.6))",0,50);
  hCMSLb->Fit("Lambdab","0");
  fTsallis->SetParameters(fTsallis->GetParameters());

  //2) Multiply scale factor to convert to 13TeV Lb---------------------------------------------------------------------------------
  Double_t binning[11] = {0,1,2,3,4,5,6,8,12,16,20};
  TH1D* hBmeson13TeV = ReadFONLL(dFONLLB13TeV,10,binning);
  TH1D* hBmeson5TeV = ReadFONLL(dFONLLB5TeV,10,binning);
  TH1D *hScaleFactor = GetRatio(hBmeson13TeV,hBmeson5TeV,"s");  //B meson ratio -> B(13TeV)/B(7TeV)

  TH1D* h13TeVLb = new TH1D("h13TeVLb","",10,binning);
  for(int i=1; i<11; i++) h13TeVLb->SetBinContent(i,fTsallis->Eval(h13TeVLb->GetBinCenter(i))*hScaleFactor->GetBinContent(i)); //Error is not assigned since we don't know the error at 1 to 10 pT

  //3) Multiply Branching ratio fraction to conver to Xib spectrum---------------------------------------------------------------------------------
  TH1D* h13TeVXib = (TH1D*) h13TeVLb->Clone("h13TeVXib");
  h13TeVXib->Scale(BRFraction);

  //4) Calculate Xib yield---------------------------------------------------------------------------------
  double ptbin[10] = {1,1,1,1,1,1,2,4,4,4};
  TH1D* h13TeVXibRaw = new TH1D("h13TeVXibRaw","",10,binning);
  Double_t Luminosity=1.86437e+09/(57.8*1000000);  //pp 13TeV integrated luminosity
  for(int i=1; i<11; i++) h13TeVXibRaw->SetBinContent(i,h13TeVXib->GetBinContent(i)*ptbin[i-1]*2*Luminosity);

  TH1D* hGenXib = (TH1D*) MCROOTFile->Get("XibGen05");  //Number of Xib in generation level
  TH1D* hRecoXib = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXib_%s_%s",Cut,CutFlag));  //Number of Xib in reconstruction level
  TH1D* hRecoeXi = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPairXib_%s_%s",Cut,CutFlag));  //Number of eXi from Xib
  TH2D* hRM_Xib = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s_Xib",Cut,CutFlag));  //Response matrix of Xib and eXi
  TH1D* hXibEff_tmp = (TH1D*) hRecoXib->Clone(Form("hXibEff_%s_%s",Cut,CutFlag));
  double binning2[10] = {1,2,3,4,5,6,8,12,16,20};
  TH1D* hXibEff = (TH1D*) hXibEff_tmp->Rebin(9,Form("hXibEff_%s_%s_rebin",Cut,CutFlag),binning2);
  hXibEff->Divide(hXibEff,hGenXib,1,1,"b");  //Xib efficiency
  TH1D* hXibEff2 = new TH1D(Form("hXibEff_%s_%s_rebin2",Cut,CutFlag),"",10,binning);
  for(int i=2; i<11; i++){
    hXibEff2->SetBinContent(i,hXibEff->GetBinContent(i-1));
    hXibEff2->SetBinError(i,hXibEff->GetBinError(i-1));
  }
  hXibEff2->SetBinContent(1,0);
  h13TeVXibRaw->Multiply(hXibEff2);

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
    delete hBmeson5TeV;
    delete hBmeson13TeV;
  }

  return hMeas;
}

TH1D* ReadFONLL(TString filename, Int_t nPtBins, Double_t* binning){
  FILE* infil=fopen(filename.Data(),"r");
  Char_t line[101];
  Char_t* rc;
  for(Int_t il=0; il<16; il++){
    rc=fgets(line,101,infil);
    if(strstr(line,"central")) break;
  }
  Float_t pt,csc;
  Int_t iPt=0;
  Double_t x[100],y[100];
  Bool_t ok;
  while(!feof(infil)){
    ok=fscanf(infil,"%f %f",&pt,&csc);
    if(feof(infil)) break;
    x[iPt]=pt;
    y[iPt]=csc;
    iPt++;
  }
  fclose(infil);

  TH1D* hfonll=new TH1D(filename.Data(),"",100,0,20);
  for(Int_t iBin=0; iBin<iPt; iBin++){
    hfonll->SetBinContent(iBin+1,y[iBin]);
    hfonll->SetBinError(iBin+1,0.);
  }
  TH1D* hfonll_rebin = (TH1D*) hfonll->Rebin(nPtBins,filename.Data(),binning);
  delete hfonll;

  return hfonll_rebin;
}
