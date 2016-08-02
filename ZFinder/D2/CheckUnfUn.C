#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TRandom.h"
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
#include <sstream>
#include <string>

#include <iomanip>

#include <iostream>
#include <fstream>

const bool doNorm=1;
const int elec=1;
const std::string Tag="";

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,5,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
const double yBins[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
size_t ny=(sizeof(yBins)/sizeof(yBins[0]))-1;
size_t nbins=nphistar*ny;

vector<TGraphAsymmErrors*> SplitGraph(TGraphAsymmErrors* graph){
  vector<TGraphAsymmErrors*> v;
  for (uint i=0; i<ny; i++){
    TGraphAsymmErrors* g=new TGraphAsymmErrors(nphistar);
    for (uint j=0; j<nphistar; j++){
      int bin=i*nphistar+j;
      double x, y;
      graph->GetPoint(bin,x,y);
      g->SetPoint(j,(phistarBins[j] + phistarBins[j + 1]) / 2.,y);
      g->SetPointError(j,0,0,graph->GetErrorYlow(bin),graph->GetErrorYhigh(bin));
    }
    v.push_back(g);
  }
  return v;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* ph, TGraphAsymmErrors* mg){
  double x,y,errorl,errorh,ymg,errorlmg,errorhmg;
  TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbins);
  for (size_t ibin=0; ibin<nbins; ibin++){
    ph->GetPoint(ibin,x,y);
    mg->GetPoint(ibin,x,ymg);
    errorlmg=mg->GetErrorYlow(ibin)/ymg;
    errorhmg=mg->GetErrorYhigh(ibin)/ymg;
    errorl=ph->GetErrorYlow(ibin)/y;
    errorh=ph->GetErrorYhigh(ibin)/y;
    g_ratio->SetPoint(ibin,x,(ymg/y)-1);
    g_ratio->SetPointError(ibin, 0, 0, sqrt((errorl*errorl)+(errorlmg*errorlmg)), sqrt((errorh*errorh)+(errorhmg*errorhmg)));
    cout<<ibin<<" "<<y<<" "<<ymg<<" "<<ymg/y-1<<" "<<errorl*100<<" "<<errorh*100<<" "<<errorlmg*100<<" "<<errorhmg*100<<" "<<sqrt((errorl*errorl)+(errorlmg*errorlmg))*100<<" "<<sqrt((errorh*errorh)+(errorhmg*errorhmg))*100<<endl;
  }
  return g_ratio;
}


TGraphAsymmErrors* CreateCopy(TGraphAsymmErrors* graph, double scale=1.0){
  double x,y,errorl,errorh;
  TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nbins);
  for (size_t ibin=0; ibin<nbins; ibin++){
    graph->GetPoint(ibin,x,y);
    errorl=graph->GetErrorYlow(ibin);
    errorh=graph->GetErrorYhigh(ibin);
    g_copy->SetPoint(ibin,x,y*scale);
    g_copy->SetPointError(ibin, 0, 0, errorl*scale, errorh*scale);
  }
  return g_copy;
}

void PlotFinal(){

  TGraphAsymmErrors* g_ph;
  TGraphAsymmErrors* g_mg;
  TGraphAsymmErrors* g_r;

  std::string textn="Output/Data_Graph_MCStat_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  cout<<textn.c_str()<<endl;
  TFile p(textn.c_str());
  TGraphAsymmErrors* g_ph_temp= (TGraphAsymmErrors*)p.Get("Graph");
  g_ph=CreateCopy(g_ph_temp);

  textn="Output/Data_Graph_MCStat_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="MG_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  cout<<textn.c_str()<<endl;
  TFile m(textn.c_str());
  TGraphAsymmErrors* g_mg_temp= (TGraphAsymmErrors*)m.Get("Graph");
  g_mg=CreateCopy(g_mg_temp);

  g_r=CreateRatio(g_ph, g_mg);
  vector<TGraphAsymmErrors*> ratio=SplitGraph(g_r);

  for (uint i=0; i<ny; i++){ 
  //  for (uint i=0; i<5; i++){ 
    std::ostringstream strs;
    strs << i;
    std::string Canvasname="UnfUn_Bin"+strs.str();
    TCanvas* FinalPhiTot = new TCanvas(Canvasname.c_str(),Canvasname.c_str(),800,900);
    gPad->SetPad("p1","p1",0,0, 1,1,kWhite,0,0);
    // gPad->SetBottomMargin(0.01);
    // gPad->SetTopMargin(0.06);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);
    ratio[i]->SetLineWidth(2);
    ratio[i]->SetMarkerStyle(20);
    ratio[i]->GetXaxis()->SetRangeUser(0.001,10);
    ratio[i]->GetXaxis()->SetTitle("#phi*");
    ratio[i]->SetTitle(0);      
    ratio[i]->GetYaxis()->SetTitle("MG/PH-1");
    ratio[i]->GetXaxis()->CenterTitle();
    ratio[i]->GetYaxis()->CenterTitle();
    ratio[i]->Draw("APE");
    TLatex mark;
    mark.SetTextSize(0.043);
    mark.SetTextFont(42);
    //mark.SetTextSize(0.035);
    mark.SetNDC(true);
    if (i==0) mark.DrawLatex(0.2,0.88,"0.0 < |y_{ee}| < 0.4");
    if (i==1) mark.DrawLatex(0.2,0.88,"0.4 < |y_{ee}| < 0.8");
    if (i==2) mark.DrawLatex(0.2,0.88,"0.8 < |y_{ee}| < 1.2");
    if (i==3) mark.DrawLatex(0.2,0.88,"1.2 < |y_{ee}| < 1.6");
    if (i==4) mark.DrawLatex(0.2,0.88,"1.6 < |y_{ee}| < 2.0");
    if (i==5) mark.DrawLatex(0.2,0.88,"2.0 < |y_{ee}| < 2.4");

    std::string plotname="Plots/Unf_Un_";
    plotname+="Bin"+strs.str()+"_";
    plotname+=Tag;
    if (doNorm) plotname+="Norm_";
    else        plotname+="Abs_";
    if (elec==0)plotname+="Dressed.";
    if (elec==1)plotname+="Born.";
    if (elec==2)plotname+="Naked.";
    FinalPhiTot->SaveAs((plotname+"pdf").c_str());
  }
}

