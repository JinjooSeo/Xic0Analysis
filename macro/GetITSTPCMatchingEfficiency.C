//#include "DrawTool.C"

TH1D* GetITSTPCMatchingEfficiency(TH2D* hevseXi, TH2D* heXivsXic0){
  Double_t SysUnc_DPG[8] = {0.011, 0.018, 0.027, 0.027, 0.023, 0.023, 0.024, 0.03}; //
  Double_t sys1[8];
  Double_t sys2[8];

  for (int i=1; i<9; i++){
    double sum1 = 0.0;
    for(int j=1; j<9; j++){
      sys1[i-1] = hevseXi->GetBinContent(j,i)*SysUnc_DPG[j-1];
      sum1 += hevseXi->GetBinContent(j,i);
    }
    sys1[i-1] = sys1[i-1]/sum1;
		cout << sys1[i-1] << endl;
  }
		cout << "---------------------" << endl; 

  for (int i=1; i<9; i++){
    double sum2 = 0.0;
    for(int j=1; j<9; j++){
      sys2[i-1] = heXivsXic0->GetBinContent(j,i)*sys1[j-1];
      sum2 += heXivsXic0->GetBinContent(j,i);
    }
    sys2[i-1] = sys2[i-1]/sum2;
		cout << sys2[i-1] << endl;
  }

  Double_t ptBinning[8] = {1.,2.,3.,4.,5.,6.,8.,12.};
  TH1D* hSysUnc_ITSTPC = new TH1D("","",7,ptBinning);
  for(int i=1; i<8; i++){
    hSysUnc_ITSTPC->SetBinContent(i,sys2[i]);
  }

  return hSysUnc_ITSTPC;
}
