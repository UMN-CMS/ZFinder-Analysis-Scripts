#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "TMatrixDSparse.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TArray.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldSvd.h"
#include <sstream>
#include <string>

#include <iostream>
#include <fstream>
using namespace RooFit;
//using namespace RooStats;

//const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
RooBinning phistarbin(nphistar,phistarBins);

void CreateRatio(TH1D* h_ratio, TH1D* graph, TH1D* graphmc){
  double y,error,ymc;
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    ymc= graphmc->GetBinContent(iphistar+1);
    y=   graph->GetBinContent(iphistar+1);
    error=graph->GetBinError(iphistar+1);
    h_ratio->SetBinContent(iphistar+1,y/ymc);
    h_ratio->SetBinError(iphistar+1,error/ymc);
  }
  return;
}

void TestMuonElec(){
  TH1D* h_Muon_Abs=new TH1D("h_Muon_Abs","h_Muon_Abs",nphistar,phistarBins);
  TH1D* h_Muon_Norm=new TH1D("h_Muon_Norm","h_Muon_Norm",nphistar,phistarBins);
  TH1D* h_Elec_Abs=new TH1D("h_Elec_Abs","h_Elec_Abs",nphistar,phistarBins);
  TH1D* h_Elec_Norm=new TH1D("h_Elec_Norm","h_Elec_Norm",nphistar,phistarBins);
  TH1D* r_Muon_Abs=new TH1D("h_Muon_Abs","h_Muon_Abs",nphistar,phistarBins);
  TH1D* r_Muon_Norm=new TH1D("h_Muon_Norm","h_Muon_Norm",nphistar,phistarBins);
  TH1D* r_Elec_Abs=new TH1D("h_Elec_Abs","h_Elec_Abs",nphistar,phistarBins);
  TH1D* r_Elec_Norm=new TH1D("h_Elec_Norm","h_Elec_Norm",nphistar,phistarBins);
  TH1D* Pull_Abs=new TH1D("Pull","Pull",20,-5,5);
  TH1D* Pull_Norm=new TH1D("Pull","Pull",20,-5,5);
  h_Muon_Abs->Sumw2();
  h_Muon_Norm->Sumw2();
  h_Elec_Abs->Sumw2();
  h_Elec_Norm->Sumw2();
  r_Muon_Abs->Sumw2();
  r_Muon_Norm->Sumw2();
  r_Elec_Abs->Sumw2();
  r_Elec_Norm->Sumw2();
  Pull_Abs->Sumw2();
  Pull_Norm->Sumw2();

  double Total_Muon_Abs=0;
  double Total_Elec_Abs=0;

  TFile fm("../Madgraph_Gen_PhiStar.root");
  TH1D *h_MT_Abs  = (TH1D*)fm.Get("PtZGen_Absolute");
  TH1D *h_MT_Norm = (TH1D*)fm.Get("PtZGen_Normalized");

  TFile f2("Final_Hist_Abs_MG_Dressed.root");
  TH1D *h_T_Abs = (TH1D*)f2.Get("h_mc");
  for (uint i=0; i<nphistar; i++){
    double val=h_T_Abs->GetBinContent(i+1);
    double error=h_T_Abs->GetBinError(i+1);
    cout<<"Abs, Bin:"<<i+1<<" val:"<<val<<"+-"<<error<<endl;
    h_Muon_Abs->SetBinContent(i+1,h_MT_Abs->GetBinContent(i+1)/1000.);
    h_Muon_Abs->SetBinError(i+1,h_MT_Abs->GetBinError(i+1)/1000.);
    h_Elec_Abs->SetBinContent(i+1,val);
    h_Elec_Abs->SetBinError(i+1,error);
    double pull_abs=h_Muon_Abs->GetBinContent(i+1)-h_Elec_Abs->GetBinContent(i+1);
    double error_pull_abs=sqrt(h_Elec_Abs->GetBinError(i+1)*h_Elec_Abs->GetBinError(i+1)+h_Muon_Abs->GetBinError(i+1)*h_Muon_Abs->GetBinError(i+1));
    pull_abs=pull_abs/error_pull_abs;
    Pull_Abs->Fill(pull_abs);
    double dphistar=phistarBins[i+1]-phistarBins[i];
    Total_Muon_Abs+=dphistar*h_MT_Abs->GetBinContent(i+1)/1000.;
    Total_Elec_Abs+=dphistar*val;
  }
  TFile f1("Final_Hist_Norm_MG_Dressed.root");
  TH1D *h_T_Norm = (TH1D*)f1.Get("h_mc");

  h_T_Abs->Divide(h_T_Norm);
  for (int i=0; i<nphistar; i++){
    cout<<"Bin "<<i<<" "<<h_T_Abs->GetBinContent(i+1)<<endl;
  }
  h_T_Abs->Draw();
  return;

 
  cout<<"Total Muon = "<<Total_Muon_Abs<<endl;
  cout<<"Total Elec = "<<Total_Elec_Abs<<endl;

  double Total_Muon=0;
  double Total_Elec=0;

  for (uint i=0; i<nphistar; i++){
    double val=h_T_Norm->GetBinContent(i+1);
    double error=h_T_Norm->GetBinError(i+1);
    cout<<"Norm, Bin:"<<i+1<<" val:"<<val<<"+-"<<error<<endl;
    h_Muon_Norm->SetBinContent(i+1,h_MT_Norm->GetBinContent(i+1));
    h_Muon_Norm->SetBinError(i+1,h_MT_Norm->GetBinError(i+1));
    h_Elec_Norm->SetBinContent(i+1,val);
    h_Elec_Norm->SetBinError(i+1,error);
    double pull_norm=h_Muon_Norm->GetBinContent(i+1)-h_Elec_Norm->GetBinContent(i+1);
    double error_pull_norm=sqrt(h_Elec_Norm->GetBinError(i+1)*h_Elec_Norm->GetBinError(i+1)+h_Muon_Norm->GetBinError(i+1)*h_Muon_Norm->GetBinError(i+1));
    pull_norm=pull_norm/error_pull_norm;
    Pull_Norm->Fill(pull_norm);
    double dphistar=phistarBins[i+1]-phistarBins[i];
    Total_Muon+=dphistar*h_MT_Norm->GetBinContent(i+1);
    Total_Elec+=dphistar*val;
  }

  cout<<"Total Muon = "<<Total_Muon<<endl;
  cout<<"Total Elec = "<<Total_Elec<<endl;
  
  TCanvas* ME_Abs = new TCanvas("ME_Abs","ME_Abs",800,900);
  ME_Abs->cd();
  ME_Abs->SetLogx();
  h_Elec_Abs->GetXaxis()->SetRangeUser(0.001,3.2);
  h_Elec_Abs->GetXaxis()->SetTitle("#phi^{*}");
  h_Elec_Abs->GetXaxis()->SetTitleOffset(0.8);
  h_Elec_Abs->GetXaxis()->SetLabelOffset(-0.01);
  h_Elec_Abs->GetXaxis()->SetTitleSize(0.04);
  h_Elec_Abs->GetXaxis()->SetLabelSize(0.04);
  h_Elec_Abs->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi_{#eta}^{*} (pb)");
  h_Elec_Abs->SetTitle("");
  h_Elec_Abs->GetYaxis()->SetTitleOffset(1.2);
  h_Elec_Abs->GetYaxis()->SetTitleSize(0.04);
  h_Elec_Abs->GetYaxis()->SetLabelSize(0.03);
  h_Elec_Abs->SetStats(0);
  h_Elec_Abs->SetBit( TH1::kNoTitle, true );
  h_Elec_Abs->SetLineColor(1);
  h_Elec_Abs->SetMarkerColor(1);
  h_Elec_Abs->SetMarkerStyle(20);
  h_Elec_Abs->Draw();
  h_Muon_Abs->SetLineColor(2);
  h_Muon_Abs->SetMarkerColor(2);
  h_Muon_Abs->SetMarkerStyle(21);
  h_Muon_Abs->Draw("same");
  //  h_Elec_Abs->Draw("same");
  TLegend* leg_Abs = new TLegend(0.45,0.77,0.85,0.91);
  leg_Abs->SetFillStyle(0);
  leg_Abs->SetBorderSize(0);
  leg_Abs->SetLineWidth(1);
  leg_Abs->SetNColumns(1);
  leg_Abs->SetTextFont(42);
  leg_Abs->SetTextSize(0.04);
  leg_Abs->AddEntry(h_Elec_Abs,"Electron","PL");
  leg_Abs->AddEntry(h_Muon_Abs,"Muon","PL");
  leg_Abs->Draw(); 

  CreateRatio(r_Elec_Abs,h_Elec_Abs,h_Elec_Abs);
  CreateRatio(r_Muon_Abs,h_Muon_Abs,h_Elec_Abs);

  TCanvas* r_ME_Abs = new TCanvas("r_ME_Abs","r_ME_Abs",800,900);
  r_ME_Abs->cd();
  r_ME_Abs->SetLogx();
  r_Elec_Abs->GetXaxis()->SetRangeUser(0.001,3.2);
  r_Elec_Abs->GetXaxis()->SetTitle("#phi^{*}");
  r_Elec_Abs->GetXaxis()->SetTitleOffset(0.8);
  r_Elec_Abs->GetXaxis()->SetLabelOffset(-0.01);
  r_Elec_Abs->GetXaxis()->SetTitleSize(0.04);
  r_Elec_Abs->GetXaxis()->SetLabelSize(0.04);
  r_Elec_Abs->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi_{#eta}^{*} (pb) / Electron distribution");
  r_Elec_Abs->SetTitle("");
  r_Elec_Abs->GetYaxis()->SetTitleOffset(1.2);
  r_Elec_Abs->GetYaxis()->SetTitleSize(0.04);
  r_Elec_Abs->GetYaxis()->SetLabelSize(0.03);
  r_Elec_Abs->SetStats(0);
  r_Elec_Abs->SetBit( TH1::kNoTitle, true );
  r_Elec_Abs->SetLineColor(1);
  r_Elec_Abs->SetMarkerColor(1);
  r_Elec_Abs->SetMarkerStyle(20);
  r_Muon_Abs->SetLineColor(2);
  r_Muon_Abs->SetMarkerColor(2);
  r_Muon_Abs->SetMarkerStyle(21);
  r_Elec_Abs->Draw();
  r_Muon_Abs->Draw("same");
  //  r_Elec_Abs->Draw("same");
  leg_Abs->Draw(); 

  TCanvas* ME_Norm = new TCanvas("ME_Norm","ME_Norm",800,900);
  ME_Norm->cd();
  ME_Norm->SetLogx();
  h_Elec_Norm->GetXaxis()->SetRangeUser(0.001,3.2);
  h_Elec_Norm->GetXaxis()->SetTitle("#phi^{*}");
  h_Elec_Norm->GetXaxis()->SetTitleOffset(0.8);
  h_Elec_Norm->GetXaxis()->SetLabelOffset(-0.01);
  h_Elec_Norm->GetXaxis()->SetTitleSize(0.04);
  h_Elec_Norm->GetXaxis()->SetLabelSize(0.04);
  h_Elec_Norm->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi_{#eta}^{*}");
  h_Elec_Norm->SetTitle("");
  h_Elec_Norm->GetYaxis()->SetTitleOffset(1.2);
  h_Elec_Norm->GetYaxis()->SetTitleSize(0.04);
  h_Elec_Norm->GetYaxis()->SetLabelSize(0.03);
  h_Elec_Norm->SetStats(0);
  h_Elec_Norm->SetBit( TH1::kNoTitle, true );
  h_Elec_Norm->SetLineColor(1);
  h_Elec_Norm->SetMarkerColor(1);
  h_Elec_Norm->SetMarkerStyle(20);
  // h_Elec_Norm->Draw();
  h_Muon_Norm->SetLineColor(2);
  h_Muon_Norm->SetMarkerColor(2);
  h_Muon_Norm->SetMarkerStyle(21);
  h_Elec_Norm->Draw();
  h_Muon_Norm->Draw("same");
  leg_Abs->Draw(); 

  CreateRatio(r_Elec_Norm,h_Elec_Norm,h_Elec_Abs);
  CreateRatio(r_Muon_Norm,h_Muon_Norm,h_Muon_Abs);

  TCanvas* r_ME_Norm = new TCanvas("r_ME_Norm","r_ME_Norm",800,900);
  r_ME_Norm->cd();
  r_ME_Norm->SetLogx();
  r_Elec_Norm->GetXaxis()->SetRangeUser(0.001,3.2);
  r_Elec_Norm->GetXaxis()->SetTitle("#phi^{*}");
  r_Elec_Norm->GetXaxis()->SetTitleOffset(0.8);
  r_Elec_Norm->GetXaxis()->SetLabelOffset(-0.01);
  r_Elec_Norm->GetXaxis()->SetTitleSize(0.04);
  r_Elec_Norm->GetXaxis()->SetLabelSize(0.04);
  r_Elec_Norm->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi_{#eta}^{*} / Electron distribution");
  r_Elec_Norm->SetTitle("");
  r_Elec_Norm->GetYaxis()->SetTitleOffset(1.2);
  r_Elec_Norm->GetYaxis()->SetTitleSize(0.04);
  r_Elec_Norm->GetYaxis()->SetLabelSize(0.03);
  r_Elec_Norm->SetStats(0);
  r_Elec_Norm->SetBit( TH1::kNoTitle, true );
  r_Elec_Norm->SetLineColor(1);
  r_Elec_Norm->SetMarkerColor(1);
  r_Elec_Norm->SetMarkerStyle(20);
  r_Muon_Norm->SetLineColor(2);
  r_Muon_Norm->SetMarkerColor(2);
  r_Muon_Norm->SetMarkerStyle(21);
  r_Elec_Norm->Draw();
  r_Muon_Norm->Draw("same");
  //  r_Elec_Norm->Draw("same");
  leg_Abs->Draw(); 

  TCanvas* CPull_Abs = new TCanvas("CPull_Abs","CPull_Abs",800,900);
  CPull_Abs->cd();
  Pull_Abs->GetXaxis()->SetTitle("(Muon-Electron)/#sqrt(2)#sigma(Electron_noLumi)");
  Pull_Abs->GetXaxis()->SetTitleOffset(1.0);
  Pull_Abs->GetXaxis()->SetLabelOffset(0.0);
  Pull_Abs->GetXaxis()->SetTitleSize(0.04);
  Pull_Abs->GetXaxis()->SetLabelSize(0.04);
  Pull_Abs->SetBit( TH1::kNoTitle, true );
  Pull_Abs->SetStats(0);
  Pull_Abs->SetLineColor(1);
  Pull_Abs->SetMarkerColor(1);
  Pull_Abs->GetYaxis()->SetTitle("#phi^{*} bins");
  Pull_Abs->GetYaxis()->SetTitleOffset(1.0);
  Pull_Abs->GetYaxis()->SetTitleSize(0.04);
  Pull_Abs->GetYaxis()->SetLabelSize(0.04);
  Pull_Abs->Draw();

  TCanvas* CPull_Norm = new TCanvas("CPull_Norm","CPull_Norm",800,900);
  CPull_Norm->cd();
  Pull_Norm->GetXaxis()->SetTitle("(Muon-Electron)/#sqrt(2)#sigma(Electron)");
  Pull_Norm->GetXaxis()->SetTitleOffset(1.0);
  Pull_Norm->GetXaxis()->SetLabelOffset(0.0);
  Pull_Norm->GetXaxis()->SetTitleSize(0.04);
  Pull_Norm->GetXaxis()->SetLabelSize(0.04);
  Pull_Norm->SetBit( TH1::kNoTitle, true );
  Pull_Norm->SetStats(0);
  Pull_Norm->SetLineColor(1);
  Pull_Norm->SetMarkerColor(1);
  Pull_Norm->GetYaxis()->SetTitle("#phi^{*} bins");
  Pull_Norm->GetYaxis()->SetTitleOffset(1.0);
  Pull_Norm->GetYaxis()->SetTitleSize(0.04);
  Pull_Norm->GetYaxis()->SetLabelSize(0.04);
  Pull_Norm->Draw();
}
