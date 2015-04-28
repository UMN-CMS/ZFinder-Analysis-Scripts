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
const int elec=0;
const int doMG=1;
const std::string Tag="";
const std::string Type="combined"; //elec, muon or combined

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

void PrintFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, bool dataonly=0){
  ofstream outputfile;
  std::string textname="Table_Values_";
  textname+=Tag;
  if (!dataonly){
    textname+=Type;
    textname+="_";
  }
  if (dataonly) textname+="MuEl_";
  if (Type=="elec" && !doMG) textname+="PH_";
  if (doNorm) textname+="Norm_";
  else        textname+="Abs_";
  if (elec==0)textname+="Dressed.txt";
  if (elec==1)textname+="Born.txt";
  if (elec==2)textname+="Naked.txt";
  outputfile.open(textname.c_str());
  std::string tableheader="";
  if (!dataonly){
    if (!doNorm) tableheader="$\\phi^*$ range & \\multicolumn{2}{c}{Data \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{MadGraph \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Powheg \\frac{d\\sigma^{fid}}{d\\phi*}(pb)}\\\\ \\hline";
    else         tableheader="$\\phi^*$ range & \\multicolumn{2}{c}{Data \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{MadGraph \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Powheg \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}}\\\\ \\hline";
  }
  else {
    if (!doNorm) tableheader="$\\phi^*$ range & \\multicolumn{2}{c}{Data: Z \\rightarrow ll \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Data: Z \\rightarrow ee \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Data: Z \\rightarrow \\mu\\mu \\frac{d\\sigma^{fid}}{d\\phi*}(pb)}\\\\ \\hline";
    else         tableheader="$\\phi^*$ range & \\multicolumn{2}{c}{Data: Z \\rightarrow ll \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Data: Z \\rightarrow ee \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Data: Z \\rightarrow \\mu\\mu \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}}\\\\ \\hline";
  }
  outputfile << tableheader << "\n";
  std::cout<< tableheader <<endl;
  for (size_t i=0; i<nphistar;i++){
    double x,y,xmg,ymg,xph,yph;
    g_data_final->GetPoint(i,x,y);
    g_mg_final->GetPoint(i,xmg,ymg);
    g_ph_final->GetPoint(i,xph,yph);
    double temp_d=g_data_final->GetErrorYhigh(i);
    double temp_m=g_mg_final->GetErrorYhigh(i);
    double temp_p=g_ph_final->GetErrorYhigh(i);
    int n_d=0;
    int n_m=0;
    int n_p=0;
    while (temp_d<1){
      temp_d=temp_d*10.;
      n_d++;
    }
    while (temp_m<1){
      temp_m=temp_m*10.;
      n_m++;
    }
    while (temp_p<1){
      temp_p=temp_p*10.;
      n_p++;
    }
    outputfile <<  std::fixed <<std::setprecision(3)<< phistarBins[i]<<"-"<<phistarBins[i+1]<<" & "<<std::setprecision(n_d)<<y<<" & "<<g_data_final->GetErrorYhigh(i)<<" & "<<std::setprecision(n_m)<<ymg<<" & "<<g_mg_final->GetErrorYhigh(i)<<" & "<<std::setprecision(n_p)<<yph<<" & "<<g_ph_final->GetErrorYhigh(i)<<" \\\\ \\hline"<< "\n";
  }
  outputfile.close();
}

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_dummy_phistar, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2=0){
  TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot","FinalPhiTot",800,900);
  FinalPhiTot->Divide(1,2);
  FinalPhiTot->cd(1);
  gPad->SetPad("p1","p1",0,2.5/9.0, 1,1,kWhite,0,0);
  gPad->SetBottomMargin(0.01);
  gPad->SetTopMargin(0.06);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  g_dummy_phistar->GetXaxis()->SetRangeUser(0.001,3.0);
  if (doNorm) g_dummy_phistar->GetYaxis()->SetRangeUser(0.0012,100.0);
  else g_dummy_phistar->GetYaxis()->SetRangeUser(0.8,40000.0);
  g_dummy_phistar->Draw("A2");
  g_mg_final->SetMarkerColor(kBlue-7);
  g_mg_final->SetLineColor(kBlue-7);
  g_mg_final->SetMarkerSize(1);
  g_mg_final->SetLineWidth(2);
  g_mg_final->SetMarkerStyle(21);
  g_mg_final->Draw("PEsame");
  g_ph_final->SetMarkerColor(kRed);
  g_ph_final->SetLineColor(kRed);
  g_ph_final->SetMarkerSize(1);
  g_ph_final->SetLineWidth(2);
  g_ph_final->SetMarkerStyle(22);
  g_ph_final->Draw("PEsame");
  g_data_final->SetFillColor(kYellow);
  g_data_final->SetMarkerSize(1);
  g_data_final->SetLineWidth(2);
  g_data_final->SetMarkerStyle(20);
  g_data_final->Draw("PEsame");
  g_data_final->SetFillColor(kYellow);
  
  TLegend* leg = new TLegend(0.53,0.72,0.85,0.91);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineWidth(1);
  leg->SetNColumns(1);
  leg->SetTextFont(42);
  
  if (!isPlot2){
    leg->AddEntry(g_data_final,"2012 data","PEF");
    if (Type=="elec"){
      leg->AddEntry(g_mg_final,"Z #rightarrow ee MadGraph","P");
      leg->AddEntry(g_ph_final,"Z #rightarrow ee Powheg","P");
    }
    if (Type=="muon"){
      leg->AddEntry(g_mg_final,"Z #rightarrow #mu#mu MadGraph","P");
      leg->AddEntry(g_ph_final,"Z #rightarrow #mu#mu Powheg","P");
    }
    if (Type=="combined"){
      leg->AddEntry(g_mg_final,"Z #rightarrow ll MadGraph","P");
      leg->AddEntry(g_ph_final,"Z #rightarrow ll Powheg","P");
    }
  }
  else {
   leg->AddEntry(g_mg_final,"2012 data: Z #rightarrow ee","P");
   leg->AddEntry(g_ph_final,"2012 data: Z #rightarrow #mu#mu","P");
   leg->AddEntry(g_data_final,"MadGraph: Z #rightarrow ll","PEF");
  }
  leg->Draw();

  TLatex mark;
  // mark.SetTextSize(0.05);
  mark.SetTextSize(0.035);
  mark.SetNDC(true);
  mark.DrawLatex(0.745,0.95,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.19,0.895,"CMS Preliminary");
  if (Type=="elec" && !isPlot2){
    mark.DrawLatex(0.19,0.20,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
    mark.DrawLatex(0.19,0.13,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
    mark.DrawLatex(0.19,0.06,"60 GeV < M_{ee} < 120 GeV");
  }
  if (Type=="muon"&& !isPlot2){
    mark.DrawLatex(0.19,0.20,"|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
    mark.DrawLatex(0.19,0.13,"p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
    mark.DrawLatex(0.19,0.06,"60 GeV < M_{#mu#mu} < 120 GeV");
  }
  if (Type=="combined" || isPlot2){
    mark.DrawLatex(0.19,0.20,"|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
    mark.DrawLatex(0.19,0.13,"p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
    mark.DrawLatex(0.19,0.06,"60 GeV < M_{ll} < 120 GeV");
  }
  FinalPhiTot->cd(2);
  gPad->SetPad("p2","p2",0,0,1,2.5/9.0,kWhite,0,0);
  gPad->SetBottomMargin(0.37);
  gPad->SetTopMargin(0.01);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetLogx(1);

  if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC  ");
  else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data   ");
  g_ratio_phistar->GetXaxis()->SetRangeUser(0.001,3.0);
  g_ratio_phistar->GetYaxis()->SetRangeUser(0.76,1.24);
  if (isPlot2) g_ratio_phistar->GetYaxis()->SetRangeUser(0.88,1.12);
  if (isPlot2 && (!doNorm)) g_ratio_phistar->GetYaxis()->SetRangeUser(0.88,1.22);
  g_ratio_phistar->GetYaxis()->SetTitleOffset(0.45);
  g_ratio_phistar->SetFillColor(kYellow);
  g_ratio_phistar->Draw("AE2");
  g_ratio_mg_phistar->SetMarkerSize(1);
  g_ratio_mg_phistar->SetLineWidth(2);
  g_ratio_mg_phistar->SetMarkerStyle(21);
  g_ratio_mg_phistar->SetMarkerColor(kBlue-7);
  g_ratio_mg_phistar->SetLineColor(kBlue-7);
  g_ratio_mg_phistar->Draw("PEsame");
  g_ratio_ph_phistar->SetMarkerSize(1);
  g_ratio_ph_phistar->SetLineWidth(2);
  g_ratio_ph_phistar->SetMarkerStyle(22);
  g_ratio_ph_phistar->SetMarkerColor(kRed);
  g_ratio_ph_phistar->SetLineColor(kRed);
  g_ratio_ph_phistar->Draw("PEsame");
  std::string plotname="ZShape_";
  plotname+=Tag;
  if (isPlot2)plotname+="MuEl";
  else plotname+=Type;
  plotname+="_";
  if (Type=="elec" && !doMG) plotname+="PH_";
  if (doNorm) plotname+="Norm_";
  else        plotname+="Abs_";
  if (elec==0)plotname+="Dressed.";
  if (elec==1)plotname+="Born.";
  if (elec==2)plotname+="Naked.";
  FinalPhiTot->SaveAs((plotname+"pdf").c_str());
  FinalPhiTot->SaveAs((plotname+"png").c_str());
  FinalPhiTot->SaveAs((plotname+"C").c_str());

  TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio","FinalPhiRatio",800,900);
  FinalPhiRatio->cd();
  FinalPhiRatio->SetLogx();
  if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC");
  else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data");
  g_ratio_phistar->SetFillColor(kYellow);
  g_ratio_phistar->GetYaxis()->SetRangeUser(0.7,1.3);
  g_ratio_phistar->GetYaxis()->SetTitleOffset(1.2);
  g_ratio_phistar->GetYaxis()->SetTitleSize(0.04);
  g_ratio_phistar->GetYaxis()->SetLabelSize(0.04);
  g_ratio_phistar->GetXaxis()->SetTitleSize(0.04);
  g_ratio_phistar->GetXaxis()->SetLabelSize(0.04);
  g_ratio_phistar->GetXaxis()->SetLabelOffset(-0.01);
  g_ratio_phistar->Draw("AE2");
  g_ratio_mg_phistar->Draw("PEsame");
  g_ratio_ph_phistar->Draw("PEsame");
  mark.SetTextSize(0.03);
  mark.DrawLatex(0.7,0.907,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.19,0.87,"CMS Preliminary");
  if (Type=="elec" && !isPlot2){
    mark.DrawLatex(0.15,0.25,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
    mark.DrawLatex(0.15,0.20,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
    mark.DrawLatex(0.15,0.15,"60 GeV < M_{ee} < 120 GeV");
  }
  if (Type=="muon"&& !isPlot2){
    mark.DrawLatex(0.15,0.25,"|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
    mark.DrawLatex(0.15,0.20,"p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
    mark.DrawLatex(0.15,0.15,"60 GeV < M_{#mu#mu} < 120 GeV");
  }
  if (Type=="combined" || isPlot2){
    mark.DrawLatex(0.15,0.25,"|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
    mark.DrawLatex(0.15,0.20,"p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
    mark.DrawLatex(0.15,0.15,"60 GeV < M_{ll} < 120 GeV");
  }
  TLegend* leg2 = new TLegend(0.53,0.72,0.85,0.88);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetLineWidth(1);
  leg2->SetNColumns(1);
  leg2->SetTextFont(42);
  
  if (!isPlot2){
    leg2->AddEntry(g_ratio_phistar,"2012 data","F");
    if (Type=="elec"){
      leg2->AddEntry(g_mg_final,"Z #rightarrow ee MadGraph","P");
      leg2->AddEntry(g_ph_final,"Z #rightarrow ee Powheg","P");
    }
    if (Type=="muon"){
      leg2->AddEntry(g_mg_final,"Z #rightarrow #mu#mu MadGraph","P");
      leg2->AddEntry(g_ph_final,"Z #rightarrow #mu#mu Powheg","P");
    }
    if (Type=="combined"){
      leg2->AddEntry(g_mg_final,"Z #rightarrow ll MadGraph","P");
      leg2->AddEntry(g_ph_final,"Z #rightarrow ll Powheg","P");
    }
  }
  else {
   leg2->AddEntry(g_mg_final,"2012 data: Z #rightarrow ee","P");
   leg2->AddEntry(g_ph_final,"2012 data: Z #rightarrow #mu#mu","P");
   leg2->AddEntry(g_data_final,"MadGraph: Z #rightarrow ll","PEF");
  }
  leg2->Draw();
  plotname="ZShape_Ratio";
  plotname+=Tag;
  if (isPlot2)plotname+="MuEl";
  else plotname+=Type;
  plotname+="_";
  if (Type=="elec" && !doMG) plotname+="PH_";
  if (doNorm) plotname+="Norm_";
  else        plotname+="Abs_";
  if (elec==0)plotname+="Dressed.";
  if (elec==1)plotname+="Born.";
  if (elec==2)plotname+="Naked.";
  FinalPhiRatio->SaveAs((plotname+"pdf").c_str());
  FinalPhiRatio->SaveAs((plotname+"png").c_str());
  FinalPhiRatio->SaveAs((plotname+"C").c_str());
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph){
  vector<TGraphAsymmErrors*> g_dummyvec;
  // double x,y,errorl,errorh;
  double x,y;
  TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graph->GetPoint(iphistar,x,y);
    // errorl=graph->GetErrorYlow(iphistar);
    // errorh=graph->GetErrorYhigh(iphistar);
    g_dummy->SetPoint(iphistar,x,y);
    g_dummy->SetPointError(iphistar, x-phistarBins[iphistar], x-phistarBins[iphistar], 0, 0);
  }
  g_dummy->GetXaxis()->SetRangeUser(0.001,3.2);
  g_dummy->GetXaxis()->SetTitleOffset(1.05); 
  g_dummy->GetXaxis()->SetTitle(0);
  g_dummy->GetXaxis()->SetTitleSize(0.05);
  g_dummy->GetXaxis()->SetLabelSize(0.05);
  g_dummy->GetYaxis()->SetTitleOffset(1.05); 
  //  if (!doNorm)g_dummy->GetYaxis()->SetTitleOffset(1.05); 
  g_dummy->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi* ");
  if (!doNorm) g_dummy->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi* (pb)");
  g_dummy->GetYaxis()->SetTitleSize(0.05);
  g_dummy->GetYaxis()->SetLabelSize(0.05);
  g_dummy->SetFillColor(kWhite);
  g_dummy->SetTitle(0);
  return g_dummy;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData){
  double x,y,errorl,errorh,xmc,ymc,errorlmc,errorhmc;
  TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graph->GetPoint(iphistar,x,y);
    errorl=graph->GetErrorYlow(iphistar);
    errorh=graph->GetErrorYhigh(iphistar);
    if (!isData){
      graphmc->GetPoint(iphistar,xmc,ymc);
      errorlmc=graph->GetErrorYlow(iphistar);
      errorhmc=graph->GetErrorYhigh(iphistar);
      g_ratio->SetPoint(iphistar,x,ymc/y);
      g_ratio->SetPointError(iphistar, 0, 0, errorlmc/y, errorhmc/y);
    }
    else{
      g_ratio->SetPoint(iphistar,x,1);
      g_ratio->SetPointError(iphistar, x-phistarBins[iphistar], x-phistarBins[iphistar], errorl/y, errorh/y);
    }
  }
  if (isData){
    g_ratio->SetLineWidth(2);
    g_ratio->GetXaxis()->SetRangeUser(0.001,3.2);
    g_ratio->GetXaxis()->SetTitle("#phi*");
    g_ratio->GetXaxis()->SetTitleOffset(1.05);
    g_ratio->GetXaxis()->SetTitleSize(0.12);
    g_ratio->GetXaxis()->SetLabelSize(0.12);
    g_ratio->GetYaxis()->SetTitle("MC/Data  ");
    g_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    g_ratio->GetYaxis()->SetTitleOffset(0.32);
    g_ratio->GetYaxis()->SetTitleSize(0.12);
    g_ratio->GetYaxis()->SetLabelSize(0.12);
    g_ratio->GetYaxis()->SetNdivisions(3,5,0);
    g_ratio->SetTitle(0);      
  }
  return g_ratio;
}

TGraphAsymmErrors* ConvertToTGraph(TH1D* h){
  TGraphAsymmErrors* g=new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    g->SetPoint(iphistar, (phistarBins[iphistar]+phistarBins[iphistar+1])/2., h->GetBinContent(iphistar+1));
    g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar+1), h->GetBinError(iphistar+1));
  }
  return g;
}

TGraphAsymmErrors* CreateCopy(TGraphAsymmErrors* graph, double scale=1.0){
  double x,y,errorl,errorh;
  TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graph->GetPoint(iphistar,x,y);
    errorl=graph->GetErrorYlow(iphistar);
    errorh=graph->GetErrorYhigh(iphistar);
    g_copy->SetPoint(iphistar,x,y*scale);
    g_copy->SetPointError(iphistar, 0, 0, errorl*scale, errorh*scale);
  }
  return g_copy;
}

void MakeFinalPlots(){

  TGraphAsymmErrors* g_data_final;
  TGraphAsymmErrors* g_mg_final;
  TGraphAsymmErrors* g_ph_final;

  if (Type=="elec"){
    std::string textn="Data_Graph_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    if (doMG) textn+="MG_";
    else textn+="PH_";
    if (elec==0)textn+="Dressed.root";
    if (elec==1)textn+="Born.root";
    if (elec==2)textn+="Naked.root";
    cout<<textn.c_str()<<endl;
    TFile dg(textn.c_str());
    TGraphAsymmErrors* g_data_temp= (TGraphAsymmErrors*)dg.Get("Graph");
    g_data_final=CreateCopy(g_data_temp);
    textn="Data_Graph_MC_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    textn+="MG_";
    if (elec==0)textn+="Dressed.root";
    if (elec==1)textn+="Born.root";
    if (elec==2)textn+="Naked.root";
    cout<<textn.c_str()<<endl;
    TFile mg(textn.c_str());  
    TGraphAsymmErrors* g_mg_temp= (TGraphAsymmErrors*)mg.Get("Graph");
    g_mg_final=CreateCopy(g_mg_temp);
    textn="Data_Graph_MC_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    textn+="PH_";
    if (elec==0)textn+="Dressed.root";
    if (elec==1)textn+="Born.root";
    if (elec==2)textn+="Naked.root";
    cout<<textn.c_str()<<endl;
    TFile fph(textn.c_str());  
    TGraphAsymmErrors* g_ph_temp= (TGraphAsymmErrors*)fph.Get("Graph");
    g_ph_final=CreateCopy(g_ph_temp);
  }

  if (Type=="combined"){
    std::string textn="Comb_Hist_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    textn+="MG_";
    textn+="Dressed.root";
    TFile tr(textn.c_str());
    TH1D* h_data_temp= (TH1D*)tr.Get("h_Comb");
    TGraphAsymmErrors* g_data_temp=ConvertToTGraph(h_data_temp);
    g_data_final=CreateCopy(g_data_temp);
    textn="Comb_Hist_MC_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    textn+="MG_";
    textn+="Dressed.root";
    TFile tr2(textn.c_str());
    TH1D* h_mg_temp= (TH1D*)tr2.Get("h_Comb");
    TGraphAsymmErrors* g_mg_temp=ConvertToTGraph(h_mg_temp);
    g_mg_final=CreateCopy(g_mg_temp);
    textn="Comb_Hist_MC_";
    textn+=Tag;
    if (doNorm) textn+="Norm_";
    else        textn+="Abs_";
    textn+="PH_";
    textn+="Dressed.root";
    TFile tr3(textn.c_str());
    TH1D* h_ph_temp= (TH1D*)tr3.Get("h_Comb");
    TGraphAsymmErrors* g_ph_temp=ConvertToTGraph(h_ph_temp);
    g_ph_final=CreateCopy(g_ph_temp);
  }

  cout<<"going to make ratio plots"<<endl;
  TGraphAsymmErrors* g_dummy_phistar      = CreateDummy(g_data_final);
  TGraphAsymmErrors* g_ratio_phistar      = CreateRatio(g_data_final,g_data_final,1);
  TGraphAsymmErrors* g_ratio_mg_phistar   = CreateRatio(g_data_final,g_mg_final,0);
  TGraphAsymmErrors* g_ratio_ph_phistar   = CreateRatio(g_data_final,g_ph_final,0);
 
  PrintFinal(g_data_final,g_mg_final,g_ph_final);
  PlotFinal(g_data_final,g_mg_final,g_ph_final,g_dummy_phistar,g_ratio_phistar,g_ratio_mg_phistar,g_ratio_ph_phistar);
}

TH1D* GetHistMuon(std::string name){
  TH1D* h_Muon=new TH1D("h_Comb","h_Comb",nphistar,phistarBins);
  h_Muon->Sumw2();
  TFile fm(name.c_str());
  TH1F *h_temp = (TH1F*)fm.Get("hRecoClone1");
  for (uint i=0; i<nphistar; i++){
    h_Muon->SetBinContent(i+1,h_temp->GetBinContent(i+1));
    h_Muon->SetBinError(i+1,h_temp->GetBinError(i+1));
  }
  return h_Muon;
} 

void MakeFinalPlots2(){

  TGraphAsymmErrors* g_mg_final;
  TGraphAsymmErrors* g_data_muon;
  TGraphAsymmErrors* g_data_elec;

  std::string textn="Data_Graph_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="MG_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  cout<<textn.c_str()<<endl;
  TFile dg(textn.c_str());
  TGraphAsymmErrors* g_data_temp= (TGraphAsymmErrors*)dg.Get("Graph");
  g_data_elec=CreateCopy(g_data_temp);
  textn="Data_Graph_MC_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn="Comb_Hist_MC_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="MG_";
  textn+="Dressed.root";
  TFile tr2(textn.c_str());
  TH1D* h_mg_temp= (TH1D*)tr2.Get("h_Comb");
  TGraphAsymmErrors* g_mg_temp=ConvertToTGraph(h_mg_temp);
  g_mg_final=CreateCopy(g_mg_temp);  
  //   textn+="MG_";
  // if (elec==0)textn+="Dressed.root";
  // if (elec==1)textn+="Born.root";
  // if (elec==2)textn+="Naked.root";
  // cout<<textn.c_str()<<endl;
  // TFile mg(textn.c_str());  
  // TGraphAsymmErrors* g_mg_temp= (TGraphAsymmErrors*)mg.Get("Graph");
  // g_mg_final=CreateCopy(g_mg_temp);
  std::string muon_name="Muon_CrossSection_For_Combination/Madgraph_Dressed_";
  if (doNorm) muon_name+="Normal_";
  else        muon_name+="Absolute_";
  muon_name+="Central_Full_Errors.root";
  TH1D* h_Muon_temp=GetHistMuon(muon_name);
  if (!doNorm){
    h_Muon_temp->Scale(1./1000.);
  }
  g_data_muon=ConvertToTGraph(h_Muon_temp);
 
  cout<<"going to make ratio plots"<<endl;
  TGraphAsymmErrors* g_dummy_phistar      = CreateDummy(g_mg_final);
  TGraphAsymmErrors* g_ratio_elec   = CreateRatio(g_mg_final,g_data_elec,0);
  TGraphAsymmErrors* g_ratio_mg     = CreateRatio(g_mg_final,g_mg_final,1);
  TGraphAsymmErrors* g_ratio_muon   = CreateRatio(g_mg_final,g_data_muon,0);
  //
  //  PrintFinal(g_data_final,g_mg_final,g_ph_final);
  PlotFinal(g_mg_final,g_data_elec,g_data_muon,g_dummy_phistar,g_ratio_mg,g_ratio_elec,g_ratio_muon,1);
  textn="Comb_Hist_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="MG_";
  textn+="Dressed.root";
  TFile tr_comb(textn.c_str());
  TH1D* h_comb_temp= (TH1D*)tr_comb.Get("h_Comb");
  TGraphAsymmErrors* g_comb_temp=ConvertToTGraph(h_comb_temp);
  TGraphAsymmErrors* g_comb_final=CreateCopy(g_comb_temp);

  PrintFinal(g_comb_final,g_data_elec,g_data_muon,1);
}
