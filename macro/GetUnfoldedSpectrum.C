//#include "DrawTool.C"

TH1D* GetUnfoldedSpectrum(TFile* WeightedROOTFile, TFile* MCROOTFile, TH1D* hMeas, const Char_t* Cut, const Char_t* CutFlag, const Char_t* method, Int_t iteration, Bool_t IsWeighted, Bool_t DrawOption = kFALSE){
  #ifdef __CINT__
    gSystem->Load("libRooUnfold");
  #endif

  double Unfoldbinning[10] = {1,2,3,4,5,6,8,12,16,20};
  double ptbinning[7] = {2., 3., 4., 5., 6., 8., 12.};

  //1) Prepare Histogram---------------------------------------------------------------------------------
  TH1D* hUnweigthedXic0 = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
  TH1D* hUnweigthedeXiPair = (TH1D*) MCROOTFile->Get(Form("hMCRecoLevPair_%s_%s",Cut,CutFlag));
  TH2D* Mat_refold = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s",Cut,CutFlag));
  TH2D* Mat_unfold = (TH2D*) MCROOTFile->Get(Form("hRPM_%s_%s_un",Cut,CutFlag));

  TH1D* hWeigthedXic0 = (TH1D*) WeightedROOTFile->Get(Form("hMCRecoLevXic0_%s_%s",Cut,CutFlag));
  TH1D* hWeigthedeXiPair = (TH1D*) WeightedROOTFile->Get(Form("hMCRecoLevPair_%s_%s",Cut,CutFlag));
  TH2D* WeigthedMat_refold = (TH2D*) WeightedROOTFile->Get(Form("hRPM_%s_%s",Cut,CutFlag));
  TH2D* WeigthedMat_unfold = (TH2D*) WeightedROOTFile->Get(Form("hRPM_%s_%s_un",Cut,CutFlag));

  //2) Unfold step of unfolding---------------------------------------------------------------------------------
  RooUnfoldResponse unfold_unweighted(hUnweigthedeXiPair,hUnweigthedXic0,Mat_unfold);
  RooUnfoldResponse refold_unweighted(hUnweigthedXic0,hUnweigthedeXiPair,Mat_refold);
  RooUnfoldResponse unfold_weighted(hWeigthedeXiPair,hWeigthedXic0,WeigthedMat_unfold);
  RooUnfoldResponse refold_weighted(hWeigthedXic0,hWeigthedeXiPair,WeigthedMat_refold);

  TH1D *hUnweightedReco = new TH1D(Form("hUnweightedReco_%s_%s",Cut,CutFlag),"",9,Unfoldbinning);
  TH1D *hWeightedReco = new TH1D(Form("hWeightedReco_%s_%s",Cut,CutFlag),"",9,Unfoldbinning);

  const Char_t* Bayes = "Bayes";
  const Char_t* Svd = "Svd";
  TString Bayes2 = "Bayes";
  TString Svd2 = "Svd";
  if(method == Bayes ||method == Bayes2){
    RooUnfoldBayes unfolding_unweighted (&unfold_unweighted, hMeas, iteration); // OR
    hUnweightedReco = (TH1D*) unfolding_unweighted.Hreco();
    RooUnfoldBayes unfolding_weighted (&unfold_weighted, hMeas, iteration); // OR
    hWeightedReco = (TH1D*) unfolding_weighted.Hreco();
  }
  if(method == Svd || method == Svd2 ){
    RooUnfoldSvd     unfold_unweighted_svd (&unfold_unweighted, hMeas, iteration);   // OR
    hUnweightedReco = (TH1D*) unfold_unweighted_svd.Hreco();
    RooUnfoldSvd unfolding_weighted_svd (&unfold_weighted, hMeas, iteration); // OR
    hWeightedReco = (TH1D*) unfolding_weighted_svd.Hreco();

    TSVDUnfold_local* SVDUnfold = (TSVDUnfold_local*) unfold_unweighted_svd.Impl();
    TH1D* fDHist = SVDUnfold->GetD();
    for (int i = 0; i<10; i++) cout <<i << " : "<< fDHist->GetBinContent(i) << endl;
  }

  cout << "method=---------------------------------" << method << endl;

  //3) Output---------------------------------------------------------------------------------
  TH1D* hUnweightedUnfolded = new TH1D(Form("hWoUnfolded_%s_%s",Cut,CutFlag),"",6,ptbinning);
  for(int i=1; i<7; i++){
    hUnweightedUnfolded->SetBinContent(i,hUnweightedReco->GetBinContent(i+1));
    hUnweightedUnfolded->SetBinError(i,hUnweightedReco->GetBinError(i+1));
  }
  TH1D* hWeightedUnfolded = new TH1D(Form("hWUnfolded_%s_%s",Cut,CutFlag),"",6,ptbinning);
  for(int i=1; i<7; i++){
    hWeightedUnfolded->SetBinContent(i,hWeightedReco->GetBinContent(i+1));
    hWeightedUnfolded->SetBinError(i,hWeightedReco->GetBinError(i+1));
  }

  if(DrawOption){
    SetStyle();
    int nCan = 2;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("UnfoldCan%d",i),"",650,500);

    can[0]->cd();
    WeigthedMat_unfold->Draw("COLZ");
    SetAxis(WeigthedMat_unfold,"#it{p}_{T}(e#Xi) (GeV/#it{c})","#it{p}_{T}(#Xi_{c}^{0}) (GeV/#it{c})");

    RooUnfoldBayes refolding_unweighted(&refold_unweighted, hUnweightedReco, 1);
    TH1D *hUnweightedBack = (TH1D*) refolding_unweighted.Hreco();
    RooUnfoldBayes refolding_weighted(&refold_weighted, hWeightedReco, 1);
    TH1D *hWeightedBack = (TH1D*) refolding_weighted.Hreco();

    can[1]->cd();
    hUnweightedReco->Draw(); HistSty(hUnweightedReco,kBlue,kFullCircle);
    hWeightedReco->Draw("SAME"); HistSty(hWeightedReco,kGreen,kFullCircle);
    hMeas->Draw("SAME"); HistSty(hMeas,kBlack,kFullCircle);
    hUnweightedBack->Draw("SAME"); HistSty(hUnweightedBack,kRed,kCircle);
    hWeightedBack->Draw("SAME"); HistSty(hWeightedBack,kYellow+2,kCircle);
    hUnweightedReco->GetYaxis()->SetRangeUser(0,1200);
    SetAxis(hUnweightedReco,"#it{p}_{T}(e#Xi) (GeV/#it{c})","Entries");
    TLegend *leg1 = new TLegend(0.55,0.65,0.76,0.82);
    leg1->AddEntry(hMeas,"measured");
    leg1->AddEntry(hUnweightedReco,"unfolded");
    leg1->AddEntry(hWeightedReco,"weighted unfolded");
    leg1->AddEntry(hUnweightedBack,"refolded");
    leg1->AddEntry(hWeightedBack,"weighted refolded");
    LegSty(leg1);
    leg1->Draw();

    can[0]->SaveAs(Form("ResponseMatrixOfXic_%s_%s.pdf",Cut,CutFlag));
    can[1]->SaveAs(Form("UnfoldRefoldOfXic_%s_%s.pdf",Cut,CutFlag));

    delete[] can;
  }

  if(!DrawOption){
    delete hUnweigthedXic0;
    delete hUnweigthedeXiPair;
    delete Mat_refold;
    delete Mat_unfold;
    delete hWeigthedXic0;
    delete hWeigthedeXiPair;
    delete WeigthedMat_refold;
    delete WeigthedMat_unfold;
    delete hUnweightedReco;
    delete hWeightedReco;
  }

  if(!IsWeighted) return hUnweightedUnfolded;
  if(IsWeighted) return hWeightedUnfolded;

  return hWeightedUnfolded;    //default
}
