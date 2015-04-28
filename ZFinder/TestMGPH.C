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

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
RooBinning phistarbin(nphistar,phistarBins);

std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/Signal_reco.root";

std::string reco_name="Combined Single Reco";

double vec_MadGraph[34][2]={{3844,0.31}, 
			    {3847,0.37},
			    {3756,0.37},
			    {3698,0.39}, 
			    {3538,0.37}, 
			    {3421,0.38}, 
			    {3246,0.32}, 
			    {3044,0.34}, 
			    {2849,0.34}, 
			    {2665,0.31},
			    {2438,0.31},
			    {2258,0.39},
			    {2092,0.33},
			    {1883,0.33},
			    {1705,0.31},
			    {1521,0.31},
			    {1345,0.32},
			    {1184,0.33},
			    {1034,0.32},
			    {884,0.30},
			    {743,0.31},
			    {611,0.31},
			    {483,0.31},
			    {374,0.31},
			    {274,0.30},
			    {184,0.30},
			    {109,0.30},
			    {57,0.36},
			    {28.8,0.45},
			    {14.8,0.63},
			    {7.8,0.68},
			    {3.6,0.81},
			    {1.81,1.07},
			    {0.92,1.31}};

double vec_Powheg[34][2]={{3895,0.74},
			  {3855,0.86},
			  {3751,0.89},
			  {3667,0.90},
			  {3587,0.83},
			  {3358,0.85},
			  {3258,0.77},
			  {3074,0.77},
			  {2862,0.83},
			  {2671,0.70},
			  {2451,0.66},
			  {2252,0.87},
			  {2057,0.73},
			  {1891,0.69},
			  {1694,0.68},
			  {1502,0.66},
			  {1341,0.68},
			  {1195,0.70},
			  {1021,0.69},
			  {887,0.68},
			  {743,0.71},
			  {603,0.68},
			  {483,0.68},
			  {370,0.67},
			  {269,0.66},
			  {182,0.66},
			  {106,0.82},
			  {56,0.81},
			  {28.5,1.02},
			  {15.0,1.42},
			  {7.7,1.64},
			  {3.5,2.02},
			  {1.79,2.64},
			  {0.88,3.01}};

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc){
  double x,y,errorl,errorh,xmc,ymc,errorlmc,errorhmc;
  TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graphmc->GetPoint(iphistar,xmc,ymc);
    graph->GetPoint(iphistar,x,y);
    errorl=graph->GetErrorYlow(iphistar);
    errorh=graph->GetErrorYhigh(iphistar);
    g_ratio->SetPoint(iphistar,x,y/ymc);
    g_ratio->SetPointError(iphistar, 0, 0, errorl/ymc, errorh/ymc);
  }
  return g_ratio;
}

void TestMGPH(){
  TGraphAsymmErrors* g_MG=new TGraphAsymmErrors(nphistar);
  TGraphAsymmErrors* g_PH=new TGraphAsymmErrors(nphistar);
  TH1D* Pull=new TH1D("Pull","Pull",20,-5,5);
  Pull->Sumw2();

  for (int i=0; i<nphistar; i++){
    g_MG->SetPoint(i,(phistarBins[i]+phistarBins[i+1])/2.,vec_MadGraph[i][0]);
    g_MG->SetPointError(i,0,0,vec_MadGraph[i][0]*vec_MadGraph[i][1]/100.,vec_MadGraph[i][0]*vec_MadGraph[i][1]/100.);
    g_PH->SetPoint(i,(phistarBins[i]+phistarBins[i+1])/2.,vec_Powheg[i][0]);
    g_PH->SetPointError(i,0,0,vec_Powheg[i][0]*vec_Powheg[i][1]/100.,vec_Powheg[i][0]*vec_Powheg[i][1]/100.);
    double pull=vec_Powheg[i][0]-vec_MadGraph[i][0];
    double error_pull=sqrt((vec_Powheg[i][0]*vec_Powheg[i][1]/100.)*(vec_Powheg[i][0]*vec_Powheg[i][1]/100.)+(vec_MadGraph[i][0]*vec_MadGraph[i][1]/100.)*(vec_MadGraph[i][0]*vec_MadGraph[i][1]/100.));
    pull=pull/error_pull;
    Pull->Fill(pull);
  }
  
  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);
  
  TCanvas* MGPH = new TCanvas("MGPH","MGPH",800,900);
  MGPH->cd();
  MGPH->SetLogx();
  g_MG->GetXaxis()->SetRangeUser(0.001,3.2);
  g_MG->GetXaxis()->SetTitle("#phi^{*}");
  g_MG->GetXaxis()->SetTitleOffset(0.8);
  g_MG->GetXaxis()->SetLabelOffset(-0.01);
  g_MG->GetXaxis()->SetTitleSize(0.04);
  g_MG->GetXaxis()->SetLabelSize(0.04);
  g_MG->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi_{#eta}^{*} (pb)");
  g_MG->SetTitle("");
  g_MG->GetYaxis()->SetTitleOffset(1.2);
  g_MG->GetYaxis()->SetTitleSize(0.04);
  g_MG->GetYaxis()->SetLabelSize(0.03);
    //    g_MG->SetStats(0);
  g_MG->SetBit( TH1::kNoTitle, true );
  g_MG->SetLineColor(1);
  g_MG->SetMarkerColor(1);
  g_MG->SetMarkerStyle(20);
  g_MG->Draw("PA");
  g_PH->SetLineColor(2);
  g_PH->SetMarkerColor(2);
  g_PH->SetMarkerStyle(21);
  g_PH->Draw("Psame");
  g_MG->Draw("Psame");
  TLegend* leg_PM = new TLegend(0.45,0.77,0.85,0.91);
  leg_PM->SetFillStyle(0);
  leg_PM->SetBorderSize(0);
  leg_PM->SetLineWidth(1);
  leg_PM->SetNColumns(1);
  leg_PM->SetTextFont(42);
  leg_PM->SetTextSize(0.04);
  leg_PM->AddEntry(g_MG,"Data unfolded with MadGraph","PL");
  leg_PM->AddEntry(g_PH,"Data unfolded with Powheg","PL");
  leg_PM->Draw(); 

  TGraphAsymmErrors* r_MG=CreateRatio(g_MG,g_MG);
  TGraphAsymmErrors* r_PH=CreateRatio(g_PH,g_MG);

  TCanvas* r_MGPH = new TCanvas("r_MGPH","r_MGPH",800,900);
  r_MGPH->cd();
  r_MGPH->SetLogx();
  r_MG->GetXaxis()->SetRangeUser(0.001,3.2);
  r_MG->GetXaxis()->SetTitle("#phi^{*}");
  r_MG->GetXaxis()->SetTitleOffset(0.8);
  r_MG->GetXaxis()->SetLabelOffset(-0.01);
  r_MG->GetXaxis()->SetTitleSize(0.04);
  r_MG->GetXaxis()->SetLabelSize(0.04);
  r_MG->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi_{#eta}^{*} (pb)");
  r_MG->SetTitle("");
  r_MG->GetYaxis()->SetTitleOffset(1.2);
  r_MG->GetYaxis()->SetTitleSize(0.04);
  r_MG->GetYaxis()->SetLabelSize(0.03);
    //    r_MG->SetStats(0);
  r_MG->SetBit( TH1::kNoTitle, true );
  r_MG->SetLineColor(1);
  r_MG->SetMarkerColor(1);
  r_MG->SetMarkerStyle(20);
  r_MG->Draw("PA");
  r_PH->SetLineColor(2);
  r_PH->SetMarkerColor(2);
  r_PH->SetMarkerStyle(21);
  r_PH->Draw("Psame");
  r_MG->Draw("Psame");
  leg_PM->Draw(); 

  TCanvas* CPull = new TCanvas("CPull","CPull",800,900);
  CPull->cd();
  Pull->GetXaxis()->SetTitle("(Powheg-MadGraph)/{#sigma}(Powheg-MadGraph)");
  Pull->GetXaxis()->SetTitleOffset(1.0);
  Pull->GetXaxis()->SetLabelOffset(0.0);
  Pull->GetXaxis()->SetTitleSize(0.04);
  Pull->GetXaxis()->SetLabelSize(0.04);
  Pull->SetBit( TH1::kNoTitle, true );
  Pull->SetStats(0);
  Pull->SetLineColor(1);
  Pull->SetMarkerColor(1);
  Pull->GetYaxis()->SetTitle("#phi^{*} bins");
  Pull->GetYaxis()->SetTitleOffset(1.0);
  Pull->GetYaxis()->SetTitleSize(0.04);
  Pull->GetYaxis()->SetLabelSize(0.04);
  Pull->Draw();
}
