#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "TSVDUnfold_local.h"
#endif

void SetStyle(Bool_t graypalette=kFALSE);
void LegSty(TLegend *leg);
void HistSty(TH1* hist, Color_t col, Style_t sty);
void SetAxis(TH1D* hist, const Char_t* Xaxis, const Char_t* Yaxis);
void SetAxis(TH2D* hist, const Char_t* Xaxis, const Char_t* Yaxis);

void DrawTool(){
  #ifdef __CINT__
    gSystem->Load("libRooUnfold");
  #endif
}

TH1D* GetRatio(TH1D* nu, TH1D* de, const Char_t* b){
  TH1D* hratio = (TH1D*) nu->Clone(Form("%f%f_ratio",nu->GetBinContent(1),de->GetBinContent(1)));
  const Char_t* bino = "b";
  if(bino == b) hratio->Divide(hratio,de,1,1,"b");
  else hratio->Divide(de);

  return hratio;
}

void LegSty(TLegend *leg){
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.045);
	leg->SetEntrySeparation(0.1);
}

void HistSty(TH1* hist, Color_t col, Style_t sty){
  hist->SetLineColor(col);
  hist->SetMarkerColor(col);
  hist->SetMarkerStyle(sty);
}

void SetStyle(Bool_t graypalette) {
  //cout << "Setting style!" << endl;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.15,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");  gStyle->SetTickLength(0.03,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}

void SetAxis(TH1D* hist, const Char_t* Xaxis, const Char_t* Yaxis){
  hist->GetXaxis()->SetTitle(Xaxis);
  hist->GetYaxis()->SetTitle(Yaxis);
}

void SetAxis(TH2D* hist, const Char_t* Xaxis, const Char_t* Yaxis){
  hist->GetXaxis()->SetTitle(Xaxis);
  hist->GetYaxis()->SetTitle(Yaxis);
}
