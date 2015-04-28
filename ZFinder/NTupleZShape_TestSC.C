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
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldSvd.h"
#include <sstream>
#include <string>

#include <iomanip>
#include <iostream>
#include <fstream>

const bool doNorm=0;
const int Ntoys=500;
const int doMG=1;

const double nMC=-1;

const double TEfficiencyEtaBins[7] = {-2.1, -1.478, -0.8, 0.0, 0.8, 1.478, 2.1};
const int TEfficiencyETBins[5] = {30, 40, 50, 200, 2000};
const double TEfficiencyData[6][4][2] = {{{0.915,  0.001}, {0.878,  0.001}, {0.850,  0.001}, {0.895,  0.021}},
					 {{0.912,  0.001}, {0.907,  0.001}, {0.894,  0.001}, {0.920,  0.010}},
					 {{0.899,  0.001}, {0.931,  0.001}, {0.919,  0.001}, {0.962,  0.006}},
					 {{0.902,  0.001}, {0.931,  0.001}, {0.920,  0.001}, {0.939,  0.007}},
					 {{0.917,  0.001}, {0.912,  0.001}, {0.898,  0.001}, {0.923,  0.010}},
					 {{0.918,  0.001}, {0.886,  0.001}, {0.856,  0.001}, {0.896,  0.019}}};
const double TEfficiencyMC[6][4][2]   = {{{0.939,  0.001}, {0.892,  0.001}, {0.866,  0.002}, {0.863,  0.036}},
					 {{0.930,  0.001}, {0.920,  0.001}, {0.913,  0.001}, {0.916,  0.016}},
					 {{0.918,  0.001}, {0.950,  0.001}, {0.943,  0.001}, {0.980,  0.007}},
					 {{0.920,  0.001}, {0.951,  0.001}, {0.945,  0.001}, {0.979,  0.007}},
					 {{0.932,  0.001}, {0.922,  0.001}, {0.915,  0.001}, {0.917,  0.015}},
					 {{0.937,  0.001}, {0.894,  0.001}, {0.871,  0.002}, {0.956,  0.022}}};

const double EfficiencyEtaBins[6] = {0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};
const int EfficiencyETBins[7] = {10, 15, 20, 30, 40, 50, 1000000};
const double EfficiencySF[5][6][3] = {{{0.977,0.024,0.054}, {0.997,0.009,0.031}, {0.982,0.003,0.012}, {0.988,0.001,0.008}, {0.990,0.001,0.004}, {0.990,0.001,0.004}},
                                      {{0.977,0.024,0.054}, {0.997,0.009,0.031}, {0.993,0.002,0.012}, {0.993,0.001,0.008}, {0.993,0.001,0.004}, {0.991,0.001,0.004}},
                                      {{1.076,0.152,0.095}, {0.952,0.025,0.070}, {1.016,0.012,0.020}, {0.985,0.004,0.009}, {0.987,0.004,0.004}, {0.974,0.009,0.006}},
                                      {{1.096,0.036,0.060}, {1.008,0.010,0.031}, {0.988,0.003,0.012}, {0.993,0.002,0.008}, {0.992,0.001,0.004}, {0.990,0.003,0.004}},
                                      {{1.096,0.036,0.060}, {1.008,0.010,0.031}, {1.002,0.004,0.012}, {1.004,0.002,0.008}, {1.005,0.002,0.004}, {0.998,0.004,0.004}}};

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

vector<std::string> Get_File_Data(){
  vector<std::string> vec;
  vec.push_back("/afs/cern.ch/work/r/ruckstuh/public/Data.root");
  return vec;
}
std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco.root";
std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Powheg_reco="/afs/cern.ch/work/r/ruckstuh/public/Powheg_reco.root";
std::string File_Powheg_gen= "/afs/cern.ch/work/r/ruckstuh/public/Powheg_gen.root";
std::string File_tt=         "/afs/cern.ch/work/r/ruckstuh/public/TTbar.root";
std::string File_tautau=     "/afs/cern.ch/work/r/ruckstuh/public/TauTau.root";
std::string File_tbarw=      "/afs/cern.ch/work/r/ruckstuh/public/TbarW.root";
std::string File_tw=         "/afs/cern.ch/work/r/ruckstuh/public/TW.root";
std::string File_ww=         "/afs/cern.ch/work/r/ruckstuh/public/WW.root";
std::string File_wz=         "/afs/cern.ch/work/r/ruckstuh/public/WZ.root";
std::string File_zz=         "/afs/cern.ch/work/r/ruckstuh/public/ZZ.root";

std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

void NormalizeGraph(vector<TGraphAsymmErrors*> &graph){
  for (size_t i=0; i<graph.size();i++){
    double xstot=0;
    double x,y,errorl,errorh;
    for (size_t iphistar=0; iphistar<nphistar; iphistar++){
      graph[i]->GetPoint(iphistar,x,y);
      xstot+=y;
    }
    //cout<<i<<":Normalisation:"<<xstot<<endl;
    for (size_t iphistar=0; iphistar<nphistar; iphistar++){
      double dphistar=phistarBins[iphistar+1]-phistarBins[iphistar];
      double norm=dphistar*xstot;
      if (!doNorm) norm=dphistar;
      graph[i]->GetPoint(iphistar,x,y);
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      graph[i]->SetPoint(iphistar,x,y/norm);
      graph[i]->SetPointError(iphistar, 0, 0, errorl/norm, errorh/norm);
    }
  }
}

TGraphAsymmErrors* CalcTotalSysU(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal, bool same=0, bool cteq=0){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y,errorl,errorh;
    graph_nominal->GetPoint(iphistar,x,y);
    g_un->SetPoint(iphistar,x,y);
    errorl=0;
    errorh=0;
    int last=0;
    int now=0;
    double diflast=0;
    for (size_t i=same; i<graph.size();i++){
      now=0;
      double xtemp,ytemp;
      graph[i]->GetPoint(iphistar,xtemp,ytemp);
      if (xtemp!=x) cout<<"This is really weird and wrong"<<endl;
      double dif=ytemp-y;
      if (dif/y>0.1)cout<<"WOW very large error: "<<dif/y<<" "<<y<<"  "<<ytemp<<endl;
      if (dif>0){errorh=sqrt(errorh*errorh+dif*dif);now=1.0;}
      if (dif<0){errorl=sqrt(errorl*errorl+dif*dif);now=-1.0;}
      if (i>0 && i%2!=same){
	if (last!=0 && now!=0 && last==now && fabs(diflast/y)>0.0001 && fabs(dif/y)>0.0001) cout<<"both change go in the same dirrection: "<< i<<"  "<<iphistar<<"  "<<dif/y<<" "<<diflast/y<<endl;
      }
      last=now;
      diflast=dif;
    }
    if (cteq){
      errorl=errorl/1.645;
      errorh=errorh/1.645;
    }
    g_un->SetPointError(iphistar, 0, 0, errorl, errorh);
  }
  return g_un;
}

TGraphAsymmErrors* CalcTotalSysU_effstat(TGraphAsymmErrors* graph, TGraphAsymmErrors* graph_nominal, TGraphAsymmErrors* graph_mctoy){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y,errorl,errorh;
    double x2,y2,errorl2,errorh2;
    double errorl3,errorh3;
    graph_nominal->GetPoint(iphistar,x,y);
    graph->GetPoint(iphistar,x2,y2);
    if (x!=x2 || y !=y2) cout<<"eff stat: nominal values don't agree: "<<x<<" "<<x2<<" : "<<y<<" "<<y2<<endl;
    g_un->SetPoint(iphistar,x,y);
    errorl=graph_nominal->GetErrorYlow(iphistar);
    errorh=graph_nominal->GetErrorYhigh(iphistar);
    errorl2=graph->GetErrorYlow(iphistar);
    errorh2=graph->GetErrorYhigh(iphistar);
    errorl3=graph_mctoy->GetErrorYlow(iphistar);
    errorh3=graph_mctoy->GetErrorYhigh(iphistar);
    if (errorl2<errorl || errorh2<errorh) std::string File_Signal_reco_cteq="/afs/cern.ch/work/r/ruckstuh/MG_cteq_reco.root";
std::string File_Signal_gen_cteq= "/afs/cern.ch/work/r/ruckstuh/MG_cteq_gen.root";
cout<<"eff stat: uncertainty size don't agree: "<<errorl<<" "<<errorl2<<" : "<<errorh<<" "<<errorh2<<endl;
    g_un->SetPointError(iphistar, 0, 0, sqrt((errorl2*errorl2)-(errorl*errorl)+(errorl3*errorl3)), sqrt((errorh2*errorh2)-(errorh*errorh)+(errorh3*errorh3)));
  }
  return g_un;
}

TGraphAsymmErrors* CalcTotalSysU_nnpdf(vector<TGraphAsymmErrors*> graph, bool useMedian=0){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y;
    vector<double> phis;
    for (size_t i=0; i<graph.size();i++){
      graph[i]->GetPoint(iphistar,x,y);
      phis.push_back(y);
    }
    g_un->SetPoint(iphistar,x,TMath::Mean(phis.begin(),phis.end()));
    double rms=TMath::RMS(phis.begin(),phis.end());
    g_un->SetPointError(iphistar, 0, 0, rms, rms);
    if (useMedian){
     std::sort(phis.begin(),phis.end());
     double med=phis[Ntoys/2];
     int idx68=0.1587*(Ntoys+1)-0.5;
     double min68=phis[idx68];
     double max68=phis[(Ntoys+1)-idx68];
     g_un->SetPoint(iphistar,x,med);
     g_un->SetPointError(iphistar, 0, 0, med-min68, max68-med);
    }
  }
  return g_un;
}

TGraphAsymmErrors* CalcTotalSysU_toymc(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal){
  TGraphAsymmErrors* g_toy=CalcTotalSysU_nnpdf(graph,1);
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y,xmc,ymc;
    graph_nominal->GetPoint(iphistar,x,y);
    g_un->SetPoint(iphistar,x,y);
    g_toy->GetPoint(iphistar,xmc,ymc);
    double eh=g_toy->GetErrorYhigh(iphistar);
    double el=g_toy->GetErrorYhigh(iphistar);
    eh=ymc-y+eh;
    el=ymc-y+el;
    if (eh<0) eh=0;
    if (el<0) el=0;
    g_un->SetPointError(iphistar, 0, 0, el, eh);
  }
  return g_un;
}

TGraphAsymmErrors* CalcTotalSysU_fsr(TGraphAsymmErrors* graph, TGraphAsymmErrors* graph_nominal){
 TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y,xtemp,ytemp,error;
    graph_nominal->GetPoint(iphistar,x,y);
    // g_un->SetPoint(iphistar,x,y);
    g_un->SetPoint(iphistar,x,y);
    graph->GetPoint(iphistar,xtemp,ytemp);
    if (xtemp!=x) cout<<"This is really weird and wrong"<<endl;
    error=fabs(ytemp-y);
    g_un->SetPointError(iphistar, 0, 0, error, error);
  }
  return g_un;
}

TGraphAsymmErrors* CalcTotalSysU_pileup(TGraphAsymmErrors* graph_1, TGraphAsymmErrors* graph_2, TGraphAsymmErrors* graph_nominal){
 TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y,xtemp1,ytemp1,xtemp2,ytemp2;
    double errorh =0;
    double errorl =0;
    graph_nominal->GetPoint(iphistar,x,y);
    // g_un->SetPoint(iphistar,x,y);
    g_un->SetPoint(iphistar,x,y);
    graph_1->GetPoint(iphistar,xtemp1,ytemp1);
    graph_2->GetPoint(iphistar,xtemp2,ytemp2);
    if (xtemp1!=x) cout<<"This is really weird and wrong"<<endl;
    if (xtemp2!=x) cout<<"This is really weird and wrong"<<endl;
    if (y-ytemp1>0) errorl=fabs(ytemp1-y);
    else errorh=fabs(ytemp1-y);
    if (y-ytemp2>0) {
      if (errorl==0) errorl=fabs(ytemp2-y);
      else{ 
	if (errorl<fabs(ytemp2-y)) errorl=fabs(ytemp2-y);
	cout<<"pile-up: both errors on same side:"<<y-ytemp1<<" "<<y-ytemp2<<endl;
      }
    }
    else  {
      if (errorh==0) errorh=fabs(ytemp2-y);
      else{ 
	if (errorh<fabs(ytemp2-y)) errorh=fabs(ytemp2-y);
	cout<<"pile-up: both errors on same side:"<<y-ytemp1<<" "<<y-ytemp2<<endl;
      }
    }
    g_un->SetPointError(iphistar, 0, 0, errorl, errorh);
  }
  return g_un;
}

TGraphAsymmErrors* GetDataFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  uint nun=graph.size();
  ofstream outputfile;
  std::string textname="Table_Un_Norm.txt";
  if (!doNorm) textname="Table_Un_Abs.txt";
  if (!doMG){
    textname="Table_Un_PH_Norm.txt";
    if (!doNorm) textname="Table_Un_PH_Abs.txt";
  }
  outputfile.open(textname.c_str());
  outputfile << "$\\phi^*$ range & Total & Unfolding";
  if (!doNorm) outputfile << "& Lumi.";
  outputfile << "& MC stat. & Bkg. & Eff. & PDF & Pile-up & FSR & Ang. res.\\\\ \\hline"<< "\n";
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double h_un=0;
    double h_pdf=0;
    double h_eff=0;
    double h_bg=0;
    double h_fsr=0;
    double h_pu=0; 
    double h_sc=0;
    double h_mc=0;
    double l_eff=0;
    double l_bg=0;
    double l_fsr=0;
    double l_pu=0; 
    double l_sc=0;
    double l_mc=0;
    double x,y,errorl,errorh;
    graph[0]->GetPoint(iphistar,x,y);
    vector<double> pdf_max,pdf_min;
    for (uint i=nun-3; i<nun; i++){
      if (slist[i]!="cteq" && slist[i]!="mstw" && slist[i]!="nnpdf"){std::cout<<"This should be a pdf uncertainty:"<<slist[i]<<std::endl;}
      graph[i]->GetPoint(iphistar,x,y);
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      pdf_max.push_back(y+errorh);
      pdf_min.push_back(y-errorl);
    }
    double max=TMath::Max(pdf_max[0],pdf_max[1]);
    max=TMath::Max(pdf_max[2],max);
    double min=TMath::Min(pdf_min[0],pdf_min[1]);
    min=TMath::Min(pdf_min[2],min);
    double mean=(max+min)/2.0;
    double rms_pdf=(max-min)/2.0;
    h_pdf=rms_pdf/mean;
    double error_sys_max=h_pdf*h_pdf;
    double error_sys_min=h_pdf*h_pdf;
    for (uint i=0; i<nun-3; i++){
      if (slist[i]=="cteq" || slist[i]=="mstw" || slist[i]=="nnpdf"){std::cout<<"This shouldn't be a pdf uncertainty:"<<slist[i]<<std::endl;}
      graph[i]->GetPoint(iphistar,x,y);
      // if (y<min || y>max) {std::cout<<"Originial phi* outside pdfrange:"<<y<<"  Range: "<<min<<"   "<<max<<std::endl; }
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      double error_max=errorh/y;
      double error_min=errorl/y;
      error_sys_max+=error_max*error_max;
      error_sys_min+=error_min*error_min;
      if (slist[i]=="unfolding") {h_un =error_max;}
      if (slist[i]=="eff")       {h_eff=error_max;l_eff=error_min;}
      if (slist[i]=="bg")        {h_bg =error_max;l_bg =error_min;}
      if (slist[i]=="fsr")       {h_fsr=error_max;l_fsr=error_min;}
      if (slist[i]=="pileup")    {h_pu=error_max;l_pu=error_min;}
      if (slist[i]=="sc")        {h_sc=error_max;l_sc=error_min;}
      if (slist[i]=="mcstat")    {h_mc=error_max;l_mc=error_min;}
      //cout<<error_sys_min<<endl;
    }
    if (!doNorm) error_sys_max+=0.026*0.026;
    if (!doNorm) error_sys_min+=0.026*0.026;
    //cout<<error_sys_min<<endl;
    g_un->SetPoint(iphistar,x,mean);
    g_un->SetPointError(iphistar,0,0,mean*sqrt(error_sys_min),mean*sqrt(error_sys_max));
    if (doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f : %9.2f : %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n",int(iphistar),mean,sqrt(error_sys_max)*mean, sqrt(error_sys_min)*mean,sqrt(error_sys_max)*100.,h_un*100., h_mc*100., h_pdf*100., h_eff*100., h_bg*100., h_fsr*100., h_pu*100., h_sc*100.);
    if (!doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f: %9.2f: %9.2f: %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n",int(iphistar),mean,sqrt(error_sys_max)*mean, sqrt(error_sys_min)*mean,sqrt(error_sys_max)*100.,h_un*100., 2.6, h_mc*100., h_pdf*100., h_eff*100., h_bg*100., h_fsr*100., h_pu*100., h_sc*100.);
    
    outputfile << std::fixed;
    outputfile << std::setprecision(3) <<phistarBins[iphistar]<<"-"<<phistarBins[iphistar+1]<<" & "<< std::setprecision(2)<<sqrt(error_sys_max)*100.<<"\\% & "<<h_un*100.<<"\\% & ";
    if (!doNorm) outputfile << "2.6\\% &";
    outputfile <<h_mc*100.<<"\\% & "<<h_bg*100.<<"\\% & "<<h_eff*100.<<"\\% & "<<h_pdf*100.<<"\\% & "<<h_pu*100.<<"\\% & "<<h_fsr*100.<<"\\% & "<<h_sc*100.<<"\\%  \\\\ \\hline"<< "\n";
  }
  outputfile.close();
  return g_un;
}

TGraphAsymmErrors* GetMCFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  uint nun=graph.size();
  ofstream outputfile;
  std::string textname="TableMC_Un_Norm.txt";
  if (!doNorm) textname="TableMC_Un_Abs.txt";
  if (!doMG) {
    textname="TableMC_Un_PH_Norm.txt";
    if (!doNorm) textname="TableMC_Un_PH_Abs.txt";
  }
  outputfile.open(textname.c_str());
  outputfile << "$\\phi^*$ range & Total & Stat.";
  outputfile << "& PDF & Pile-up & FSR \\\\ \\hline"<< "\n";
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double e_un=0;
    double e_pdf=0;
    double e_fsr=0;
    double e_pu=0;
    double x,y,errorl,errorh;
    graph[0]->GetPoint(iphistar,x,y);
    vector<double> pdf_max,pdf_min;
    for (uint i=nun-3; i<nun; i++){
      if (slist[i]!="cteq" && slist[i]!="mstw" && slist[i]!="nnpdf"){std::cout<<"This should be a pdf uncertainty:"<<slist[i]<<std::endl;}
      graph[i]->GetPoint(iphistar,x,y);
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      pdf_max.push_back(y+errorh);
      pdf_min.push_back(y-errorl);
    }
    double max=TMath::Max(pdf_max[0],pdf_max[1]);
    max=TMath::Max(pdf_max[2],max);
    double min=TMath::Min(pdf_min[0],pdf_min[1]);
    min=TMath::Min(pdf_min[2],min);
    double mean=(max+min)/2.0;
    double rms_pdf=(max-min)/2.0;
    e_pdf=rms_pdf/mean;
    double error_sys_max=e_pdf*e_pdf;
    double error_sys_min=e_pdf*e_pdf;
    for (uint i=0; i<nun-3; i++){
     if (slist[i]=="cteq" || slist[i]=="mstw" || slist[i]=="nnpdf"){std::cout<<"This shouldn't be a pdf uncertainty:"<<slist[i]<<std::endl;}
      graph[i]->GetPoint(iphistar,x,y);
      if (y<min || y>max) {std::cout<<"Originial phi* outside pdfrange:"<<y<<"  Range: "<<min<<"   "<<max<<std::endl; }
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      double error_max=errorh/y;
      double error_min=errorl/y;
      error_sys_max+=error_max*error_max;
      error_sys_min+=error_min*error_min;
      if (slist[i]=="unfolding") {e_un =error_max;}
      if (slist[i]=="fsr")       {e_fsr=error_max;}
      if (slist[i]=="pileup")    {e_pu=error_max;}
       //cout<<error_sys_min<<endl;
    }
    //cout<<error_sys_min<<endl;
    g_un->SetPoint(iphistar,x,mean);
    g_un->SetPointError(iphistar,0,0,mean*sqrt(error_sys_min),mean*sqrt(error_sys_max));
    printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f : %9.2f%% : %9.2f : %9.2f%% : %9.2f%% \n",int(iphistar),mean,sqrt(error_sys_max)*mean, sqrt(error_sys_min)*mean,sqrt(error_sys_max)*100.,e_un*100., e_pdf*100., e_fsr*100., e_pu*100.);
    
    outputfile << std::fixed;
    outputfile << std::setprecision(3) <<phistarBins[iphistar]<<"-"<<phistarBins[iphistar+1]<<" & "<< std::setprecision(2)<<sqrt(error_sys_max)*100.<<"\\% & "<<e_un*100.<<"\\% & ";
    outputfile <<e_pdf*100.<<"\\% & "<<e_pu*100.<<"\\% & "<<e_fsr*100.<<"\\%  \\\\ \\hline"<< "\n";
  }
  outputfile.close();
  cout<<"GetMCFinal"<<endl;
  return g_un;
}

void PrintFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mc_final){
  ofstream outputfile;
  std::string textname="Table_Values_Norm.txt";
  if (!doNorm) textname="Table_Values_Abs.txt";
  if (!doMG){
    textname="Table_Values_PH_Norm.txt";
    if (!doNorm) textname="Table_Values_PH_Abs.txt";
  }
  outputfile.open(textname.c_str());
  if (doNorm) outputfile << "$\\phi^*$ range & \\multicolumn{2}{c}{Data (pb)} & \\multicolumn{2}{c}{MadGraph (pb)}\\\\ \\hline"<< "\n";
  else outputfile << "$\\phi^*$ range & \\multicolumn{2}{c}{Data} & \\multicolumn{2}{c}{MadGraph}\\\\ \\hline"<< "\n";
  std::cout<< "$\\phi^*$ range & Data & MadGraph \\\\ \\hline"<<endl;
  for (size_t i=0; i<nphistar;i++){
    double x,y,xmc,ymc;
    g_data_final->GetPoint(i,x,y);
    g_mc_final->GetPoint(i,xmc,ymc);
    double temp_d=g_data_final->GetErrorYhigh(i);
    double temp_m=g_mc_final->GetErrorYhigh(i);
    int n_d=0;
    int n_m=0;
    while (temp_d<1){
      temp_d=temp_d*10.;
      n_d++;
    }
    while (temp_m<1){
      temp_m=temp_m*10.;
      n_m++;
    }
    outputfile <<  std::fixed <<std::setprecision(3)<< phistarBins[i]<<"-"<<phistarBins[i+1]<<" & "<<std::setprecision(n_d)<<y<<" & "<<g_data_final->GetErrorYhigh(i)<<" & "<<std::setprecision(n_m)<<ymc<<" & "<<g_mc_final->GetErrorYhigh(i)<<" \\\\ \\hline"<< "\n";
  }
  outputfile.close();
}

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mc_final, TGraphAsymmErrors* g_dummy_phistar, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratiomc_phistar){
  TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot","FinalPhiTot",800,900);
  FinalPhiTot->Divide(1,2);
  FinalPhiTot->cd(1);
  gPad->SetPad("p1","p1",0,2.5/9.0, 1,1,kWhite,0,0);
  gPad->SetBottomMargin(0.01);
  gPad->SetTopMargin(0.06);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetLogx(1);
  g_dummy_phistar->GetXaxis()->SetRangeUser(0.001,3.0);
  g_dummy_phistar->Draw("A2");
  g_mc_final->SetMarkerColor(kBlue-7);
  g_mc_final->SetLineColor(kBlue-7);
  g_mc_final->Draw("PEsame");
  g_data_final->Draw("PEsame");
  
  TLegend* leg = new TLegend(0.53,0.77,0.90,0.91);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineWidth(1);
  leg->SetNColumns(1);
  leg->SetTextFont(42);
  
  leg->AddEntry(g_mc_final,"2012 data","P");
  leg->AddEntry(g_data_final,"SC","P");
  leg->Draw();
  
  TLatex mark;
  mark.SetTextSize(0.05);
  mark.SetTextFont(42);
  mark.SetNDC(true);
  mark.DrawLatex(0.19,0.30,"#sqrt{s} = 8 TeV");
  mark.DrawLatex(0.19,0.23,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
  mark.DrawLatex(0.19,0.16,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
  mark.DrawLatex(0.19,0.09,"60 GeV < M_{ee} < 120 GeV");
  mark.DrawLatex(0.6,0.6,"CMS Preliminary");
  
  FinalPhiTot->cd(2);
  gPad->SetPad("p2","p2",0,0,1,2.5/9.0,kWhite,0,0);
  gPad->SetBottomMargin(0.37);
  gPad->SetTopMargin(0.01);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetLogx(1);
  
  g_ratiomc_phistar->GetXaxis()->SetRangeUser(0.001,3.0);
  g_ratiomc_phistar->GetYaxis()->SetRangeUser(0.8,1.2);
  g_ratiomc_phistar->Draw("AE2");
  g_ratio_phistar->Draw("PEsame");
  if (doMG){
    FinalPhiTot->Print("ZShape_Phistar.png");
    FinalPhiTot->Print("ZShape_PhiStar.C");
  }
  else{
    FinalPhiTot->Print("ZShape_Powheg_Phistar.png");
    FinalPhiTot->Print("ZShape_Powheg_PhiStar.C");
  }
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph){
  vector<TGraphAsymmErrors*> g_dummyvec;
  double x,y,errorl,errorh;
  TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graph->GetPoint(iphistar,x,y);
    errorl=graph->GetErrorYlow(iphistar);
    errorh=graph->GetErrorYhigh(iphistar);
    g_dummy->SetPoint(iphistar,x,y);
    g_dummy->SetPointError(iphistar, x-phistarBins[iphistar], x-phistarBins[iphistar], errorl, errorh);
  }
  g_dummy->GetXaxis()->SetRangeUser(0.001,3.2);
  g_dummy->GetXaxis()->SetTitleOffset(1.05); 
  g_dummy->GetXaxis()->SetTitle(0);
  g_dummy->GetXaxis()->SetTitleSize(0.20*2/6.5);
  g_dummy->GetXaxis()->SetLabelSize(0.20*2/6.5);
  g_dummy->GetYaxis()->SetTitleOffset(1.05); 
  g_dummy->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi_{#eta}^{*} ");
  if (!doNorm) g_dummy->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi_{#eta}^{*} (pb)");
  g_dummy->GetYaxis()->SetTitleSize(0.20*2/6.5);
  g_dummy->GetYaxis()->SetLabelSize(0.20*2/6.5);
  g_dummy->SetFillColor(kWhite);
  g_dummy->SetTitle(0);
  return g_dummy;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData){
  double x,y,errorl,errorh,xmc,ymc,errorlmc,errorhmc;
  TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    graphmc->GetPoint(iphistar,xmc,ymc);
    errorlmc=graphmc->GetErrorYlow(iphistar);
    errorhmc=graphmc->GetErrorYhigh(iphistar);
    if (isData){
      graph->GetPoint(iphistar,x,y);
      errorl=graph->GetErrorYlow(iphistar);
      errorh=graph->GetErrorYhigh(iphistar);
      g_ratio->SetPoint(iphistar,x,y/ymc);
      g_ratio->SetPointError(iphistar, 0, 0, errorl/ymc, errorh/ymc);
    }
    else{
      g_ratio->SetPoint(iphistar,xmc,1);
      g_ratio->SetPointError(iphistar, xmc-phistarBins[iphistar], xmc-phistarBins[iphistar], errorlmc/ymc, errorhmc/ymc);
      g_ratio->SetLineColor(kBlue-7);
    }
  }
  if (!isData){
    g_ratio->SetLineWidth(2);
    g_ratio->SetMarkerColor(kBlue-7);
    g_ratio->SetFillColor(kBlue-7);
    g_ratio->GetXaxis()->SetRangeUser(0.001,3.2);
    g_ratio->GetXaxis()->SetTitle("#phi^{*}");
    g_ratio->GetXaxis()->SetTitleOffset(1.05);
    g_ratio->GetXaxis()->SetTitleSize(0.20 * 2/2.5);
    g_ratio->GetXaxis()->SetLabelSize(0.20 * 2/2.5);
    g_ratio->GetYaxis()->SetTitle("Data/MC");
    g_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    g_ratio->GetYaxis()->SetTitleOffset(0.32);
    g_ratio->GetYaxis()->SetTitleSize(0.20 * 2/2.5);
    g_ratio->GetYaxis()->SetLabelSize(0.20 * 2/2.5);
    g_ratio->GetYaxis()->SetNdivisions(3,5,0);
    g_ratio->SetTitle(0);      
  }
  return g_ratio;
}

TGraphAsymmErrors* CreateRatio2(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc){
  double x,y,errorl,errorh,xmc,ymc;
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

void PreparePlots(TGraphAsymmErrors &graph, bool isData){
  graph.SetMarkerSize(1);
  graph.SetLineWidth(2);
  if (isData){
    graph.SetMarkerStyle(20);
  }
  else{
    graph.SetMarkerColor(kBlue-7);
    graph.SetLineColor(kBlue-7);
    graph.SetMarkerStyle(21);
  }
}

void PlotEff(TH1D* h_eff){
  TCanvas* Efficiency = new TCanvas("Efficiency","Efficiency",800,900);
  Efficiency->cd();
  Efficiency->SetLogx();
  h_eff->GetXaxis()->SetRangeUser(0.001,3.2);
  h_eff->GetXaxis()->SetTitle("#phi^{*}_{generated}");
  h_eff->GetXaxis()->SetTitleOffset(0.8);
  h_eff->GetXaxis()->SetTitleSize(0.04);
  h_eff->GetXaxis()->SetLabelOffset(-0.01);
  h_eff->GetXaxis()->SetLabelSize(0.04);
  h_eff->GetYaxis()->SetTitle("N_{reconstructed}/N_{generated}");
  h_eff->GetYaxis()->SetTitleOffset(1.2);
  h_eff->GetYaxis()->SetTitleSize(0.04);
  h_eff->GetYaxis()->SetLabelSize(0.04);
  h_eff->SetStats(0);
  h_eff->SetBit( TH1::kNoTitle, true );
  h_eff->SetLineColor(1);
  h_eff->Draw();
  TLegend* leg_eff = new TLegend(0.45,0.77,0.85,0.91);
  leg_eff->SetFillStyle(0);
  leg_eff->SetBorderSize(0);
  leg_eff->SetLineWidth(1);
  leg_eff->SetNColumns(1);
  leg_eff->SetTextFont(42);
  leg_eff->SetTextSize(0.04);
  leg_eff->AddEntry(h_eff,"MadGraph","PL");
  leg_eff->Draw();
  return;
}

void PrintBG(TH1D* Data,TH1D* h_tt,TH1D* h_tautau,TH1D* h_tbarw,TH1D* h_tw,TH1D* h_ww,TH1D* h_wz,TH1D* h_zz){
  double data_sel=Data->GetSumOfWeights();
  double tt_sel=h_tt->GetSumOfWeights();
  double tautau_sel=h_tautau->GetSumOfWeights();
  double tbarw_sel=h_tbarw->GetSumOfWeights();
  double tw_sel=h_tw->GetSumOfWeights();
  double ww_sel=h_ww->GetSumOfWeights();
  double wz_sel=h_wz->GetSumOfWeights();
  double zz_sel=h_zz->GetSumOfWeights();

  double t_bg=tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel;

  cout<<"Weighted number of events:"<<endl;
  cout<<"Data: "<<data_sel<<"  tt: "<<tt_sel<<"  tautau: "<<tautau_sel<<"  tbarw: "<<tbarw_sel<<"  tw: "<<tw_sel<<" singletop: " <<tbarw_sel+tw_sel<<"  ww: "<<ww_sel<<"  wz: "<<wz_sel<<"  zz: "<<zz_sel<< endl;
  cout<<"ratio:"          <<"  tt: "<<tt_sel*100./data_sel<<"  tautau: "<<tautau_sel*100./data_sel<<"  tbarw: "<<tbarw_sel*100./data_sel<<"  tw: "<<tw_sel*100./data_sel<<" singletop: " <<(tbarw_sel+tw_sel)*100./data_sel<< "  ww: "<<ww_sel*100./data_sel<<"  wz: "<<wz_sel*100./data_sel<<"  zz: "<<zz_sel*100./data_sel<< "data: "<< (data_sel-t_bg)*100./data_sel<<endl;
  cout<<"ratio of bg:"    <<"  tt: "<<tt_sel*100./t_bg    <<"  tautau: "<<tautau_sel*100./t_bg    <<"  tbarw: "<<tbarw_sel*100./t_bg    <<"  tw: "<<tw_sel*100./t_bg    <<" singletop: " <<(tbarw_sel+tw_sel)*100./t_bg    << "  ww: "<<ww_sel*100./t_bg    <<"  wz: "<<wz_sel*100./t_bg    <<"  zz: "<<zz_sel*100./t_bg  <<endl;
  cout<<"total: "<<t_bg<<" "<<(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel)<<"  persentage: "<<(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel)*100./data_sel<<endl;
  cout<<"wz+zz: "<<(wz_sel+zz_sel)<<"  persentage: "<<(wz_sel+zz_sel)*100./data_sel<<" of background: "<<(wz_sel+zz_sel)*100./(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel)<<endl;
  cout<<"e+mu: "<<(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel)<<"  persentage: "<<(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel)*100./data_sel<<" of background: "<<(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel)*100./(tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel)<<endl;
}

void GetToyResponse(vector<RooUnfoldResponse*> &BinM, vector<TH2D*> BinMigration){
  for (int i=0; i<Ntoys+1; i++){
    RooUnfoldResponse* BinM_temp  =new RooUnfoldResponse (0,0,BinMigration[i]);
    BinM.push_back(BinM_temp);
  }
  return;
}

vector<TH1D*> RemoveBG(TH1D* Data,TH1D* bg_sf,vector<TH1D*> h_tt,vector<TH1D*> h_tautau,vector<TH1D*> h_tbarw,vector<TH1D*> h_tw,vector<TH1D*> h_ww,vector<TH1D*> h_wz,vector<TH1D*> h_zz){
  vector<TH1D*> h_data;
  for (uint i=0; i<h_tt.size();i++){
    TH1D* bgtemp=(TH1D*)h_tt[i]->Clone();
    bgtemp->Add(h_tautau[i],1.0);
    bgtemp->Add(h_tbarw[i],1.0);
    bgtemp->Add(h_tw[i],1.0);
    bgtemp->Add(h_ww[i],1.0);
    bgtemp->Multiply(bg_sf);
    TH1D* datatemp  =(TH1D*)Data->Clone();
    datatemp->Add(bgtemp,-1.0);
    datatemp->Add(h_wz[i],-1.0);
    datatemp->Add(h_zz[i],-1.0);
    h_data.push_back(datatemp);
  }
  return h_data;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, vector<TH1D*> mc_truegen){
  vector<TH1D*> h_eff;
  for (uint i=0; i<mc_truereco.size(); i++){
    TH1D* h_efftemp=new TH1D("h_eff","h_eff",nphistar,phistarBins);
    TH1D* gentemp  =(TH1D*)mc_truegen[i]->Clone();
    h_efftemp->Divide(mc_truereco[i],gentemp,1.,1.,"B");
    h_eff.push_back(h_efftemp);
  }
  return h_eff;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, TH1D* mc_truegen){
  vector<TH1D*> h_eff;
  for (uint i=0; i<mc_truereco.size(); i++){
    TH1D* h_efftemp=new TH1D("h_eff","h_eff",nphistar,phistarBins);
    TH1D* gentemp  =(TH1D*)mc_truegen->Clone();
    h_efftemp->Divide(mc_truereco[i],gentemp,1.,1.,"B");
    h_eff.push_back(h_efftemp);
  }
  return h_eff;
}
TGraphAsymmErrors* ConvertToTGraph(TH1D* h){
  TGraphAsymmErrors* g=new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    g->SetPoint(iphistar, (phistarBins[iphistar]+phistarBins[iphistar+1])/2., h->GetBinContent(iphistar+1));
    g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar+1), h->GetBinError(iphistar+1));
  }
  return g;
}

TH1D* ConvertToHist(TGraphAsymmErrors* g, bool data){
  TH1D* h_temp;
  if (data) h_temp=new TH1D("h_data","h_data",nphistar,phistarBins);
  if (!data) h_temp=new TH1D("h_mc","h_mc",nphistar,phistarBins);
  h_temp->Sumw2();
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y;
    g->GetPoint(iphistar,x,y);
    double error =g->GetErrorYhigh(iphistar);
    h_temp->SetBinContent(iphistar+1,y);
    h_temp->SetBinError(iphistar+1,error);
  }
  return h_temp;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(vector<RooUnfoldResponse*> BinM, TH1D* h_data, TH1D* h_eff){
  vector<TGraphAsymmErrors *> g_data;
  size_t n= BinM.size();
  for (size_t u=0; u<n; u++){
    // cout<<u<<" "<<b<<" "<<e<<" "<<d<<endl;
    RooUnfoldBayes  unfoldBay_data(BinM[u],h_data, 4);
    TH1D* h_BinMBay_data= (TH1D*) unfoldBay_data.Hreco();  
    TH1D* h_Unfolded_temp= (TH1D*) unfoldBay_data.Hreco();
    h_Unfolded_temp->Divide(h_BinMBay_data,h_eff);
    
    TGraphAsymmErrors* g_data_temp=ConvertToTGraph(h_Unfolded_temp);    
    g_data.push_back(g_data_temp);
  }
  return g_data;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(RooUnfoldResponse* BinM, vector<TH1D*> h_data, TH1D* h_eff){
  vector<TGraphAsymmErrors *> g_data;
  size_t n= h_data.size();
  for (size_t u=0; u<n; u++){
    RooUnfoldBayes  unfoldBay_data(BinM,h_data[u], 4);
    TH1D* h_BinMBay_data= (TH1D*) unfoldBay_data.Hreco();  
    TH1D* h_Unfolded_temp= (TH1D*) unfoldBay_data.Hreco();
    h_Unfolded_temp->Divide(h_BinMBay_data,h_eff);
    
    TGraphAsymmErrors* g_data_temp=ConvertToTGraph(h_Unfolded_temp);    
    g_data.push_back(g_data_temp);
  }
  return g_data;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(RooUnfoldResponse* BinM, TH1D* h_data, TH1D* h_eff){
  vector<TGraphAsymmErrors *> g_data;
  RooUnfoldBayes  unfoldBay_data(BinM,h_data, 4);
  TH1D* h_BinMBay_data= (TH1D*) unfoldBay_data.Hreco();  
  TH1D* h_Unfolded_temp= (TH1D*) unfoldBay_data.Hreco();
  h_Unfolded_temp->Divide(h_BinMBay_data,h_eff);
    
  TGraphAsymmErrors* g_data_temp=ConvertToTGraph(h_Unfolded_temp);    
  g_data.push_back(g_data_temp);

  return g_data;
}

vector<TH1D*> EmptyHVec(int n){  
  vector<TH1D*> h;
  for (int i=0; i<n; i++){
    TH1D *phistartemp= new TH1D("phistar","phistar",nphistar,phistarBins);
    phistartemp->Sumw2();
    h.push_back(phistartemp);
  }
  return h;
}
vector<RooUnfoldResponse*> EmptyBinMVec(int n){  
  vector<RooUnfoldResponse*> BinM;
  for (int i=0; i<n; i++){
    TH1D *phistartemp= new TH1D("phistar","phistar",nphistar,phistarBins);
    phistartemp->Sumw2();
    RooUnfoldResponse* responsetemp=new RooUnfoldResponse (phistartemp,phistartemp);
    BinM.push_back(responsetemp);
  }
  return BinM;
}

void FillRecoEffFluc(int &idx, vector<TH1D*> &h_eff, double weight, double phistar, double E0_pt, double E0_eta, double E1_pt, double E1_eta, bool doBinM, vector<RooUnfoldResponse*> BinM_eff, double phistar_true=0){
  for (int e=0; e<5; e++){
    for (int p=2; p<6; p++){
      double weight_temp_h=weight;
      double weight_temp_l=weight;
      if (E0_pt>EfficiencyETBins[p] && E0_pt<EfficiencyETBins[p+1] && fabs(E0_eta)>EfficiencyEtaBins[e] && fabs(E0_eta)<EfficiencyEtaBins[e+1]){
	double corr= sqrt(pow(EfficiencySF[e][p][1]/EfficiencySF[e][p][0],2.)+pow(EfficiencySF[e][p][2]/EfficiencySF[e][p][0],2.));
	weight_temp_h=weight_temp_h*(1+corr);
	weight_temp_l=weight_temp_l*(1-corr);
      }
      if (E1_pt>EfficiencyETBins[p] && E1_pt<EfficiencyETBins[p+1] && fabs(E1_eta)>EfficiencyEtaBins[e] && fabs(E1_eta)<EfficiencyEtaBins[e+1]){
	double corr= sqrt(pow(EfficiencySF[e][p][1]/EfficiencySF[e][p][0],2.)+pow(EfficiencySF[e][p][2]/EfficiencySF[e][p][0],2.));
	weight_temp_h=weight_temp_h*(1+corr);
	weight_temp_l=weight_temp_l*(1-corr);
      }
      if (!doBinM){
	h_eff[idx]->Fill(phistar,weight_temp_h);
	h_eff[idx+1]->Fill(phistar,weight_temp_l);
      }
      else {
	h_eff[idx]->Fill(phistar_true,weight_temp_h);
	h_eff[idx+1]->Fill(phistar_true,weight_temp_l);
	BinM_eff[idx]         ->Fill(phistar,phistar_true,weight_temp_h); 
	BinM_eff[idx+1]       ->Fill(phistar,phistar_true,weight_temp_l);
      }
      idx=idx+2;
    }
  }
}

void FillTrigEffFluc(int &idx, vector<TH1D*> &h_eff, double weight, double phistar, int binP0, int binE0, int binP1, int binE1, bool doBinM, vector<RooUnfoldResponse*> BinM_eff, double phistar_true=0){
  for (int e=0; e<6; e++){
    for (int p=0; p<4; p++){
      double weight_temp_h_d=weight;
      double weight_temp_l_d=weight;
      double weight_temp_h_m=weight;
      double weight_temp_l_m=weight;
      bool bin_ele0=0;
      bool bin_ele1=0;
      if (binE0==e&& binP0==p) bin_ele0=1;
      if (binE1==e&& binP1==p) bin_ele1=1;
      if ((bin_ele0 && (binE1<0 || binP1<0)) || (bin_ele1 && (binE0<0 || binP0<0))){
	double eff_old = TEfficiencyData[e][p][0]/TEfficiencyMC[e][p][0];
	double eff_d_h = (TEfficiencyData[e][p][0]+TEfficiencyData[e][p][1])/TEfficiencyMC[e][p][0];
	double eff_d_l = (TEfficiencyData[e][p][0]-TEfficiencyData[e][p][1])/TEfficiencyMC[e][p][0];
	double eff_m_h = TEfficiencyData[e][p][0]/(TEfficiencyMC[e][p][0]+TEfficiencyMC[e][p][1]);
	double eff_m_l = TEfficiencyData[e][p][0]/(TEfficiencyMC[e][p][0]-TEfficiencyMC[e][p][1]);
	weight_temp_h_d=weight_temp_h_d*eff_d_h/eff_old;
	weight_temp_l_d=weight_temp_l_d*eff_d_l/eff_old;
	weight_temp_h_m=weight_temp_h_m*eff_m_h/eff_old;
	weight_temp_l_m=weight_temp_l_m*eff_m_l/eff_old;
      }
      else if (bin_ele0 && bin_ele1){
	double eff_old=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_h=(1 - (1-TEfficiencyData[binE0][binP0][0]-TEfficiencyData[binE0][binP0][1]) * (1-TEfficiencyData[binE1][binP1][0]-TEfficiencyData[binE1][binP1][1])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_l=(1 - (1-TEfficiencyData[binE0][binP0][0]+TEfficiencyData[binE0][binP0][1]) * (1-TEfficiencyData[binE1][binP1][0]+TEfficiencyData[binE1][binP1][1])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_m_h=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]-TEfficiencyMC[binE0][binP0][1]) * (1-TEfficiencyMC[binE1][binP1][0]-TEfficiencyMC[binE1][binP1][1]));
	double eff_m_l=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]+TEfficiencyMC[binE0][binP0][1]) * (1-TEfficiencyMC[binE1][binP1][0]+TEfficiencyMC[binE1][binP1][1]));
	weight_temp_h_d=weight_temp_h_d*eff_d_h/eff_old;
	weight_temp_l_d=weight_temp_l_d*eff_d_l/eff_old;
	weight_temp_h_m=weight_temp_h_m*eff_m_h/eff_old;
	weight_temp_l_m=weight_temp_l_m*eff_m_l/eff_old;
      }
      else if (bin_ele0){
	double eff_old=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_h=(1 - (1-TEfficiencyData[binE0][binP0][0]-TEfficiencyData[binE0][binP0][1]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_l=(1 - (1-TEfficiencyData[binE0][binP0][0]+TEfficiencyData[binE0][binP0][1]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_m_h=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]-TEfficiencyMC[binE0][binP0][1]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_m_l=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]+TEfficiencyMC[binE0][binP0][1]) * (1-TEfficiencyMC[binE1][binP1][0]));
	weight_temp_h_d=weight_temp_h_d*eff_d_h/eff_old;
	weight_temp_l_d=weight_temp_l_d*eff_d_l/eff_old;
	weight_temp_h_m=weight_temp_h_m*eff_m_h/eff_old;
	weight_temp_l_m=weight_temp_l_m*eff_m_l/eff_old;
      }
      else if (bin_ele1){
	double eff_old=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_h=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0]-TEfficiencyData[binE1][binP1][1])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_d_l=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0]+TEfficiencyData[binE1][binP1][1])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]));
	double eff_m_h=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]-TEfficiencyMC[binE1][binP1][1]));
	double eff_m_l=(1 - (1-TEfficiencyData[binE0][binP0][0]) * (1-TEfficiencyData[binE1][binP1][0])) / (1 - (1-TEfficiencyMC[binE0][binP0][0]) * (1-TEfficiencyMC[binE1][binP1][0]+TEfficiencyMC[binE1][binP1][1]));
	weight_temp_h_d=weight_temp_h_d*eff_d_h/eff_old;
	weight_temp_l_d=weight_temp_l_d*eff_d_l/eff_old;
	weight_temp_h_m=weight_temp_h_m*eff_m_h/eff_old;
	weight_temp_l_m=weight_temp_l_m*eff_m_l/eff_old;
      }
      if (!doBinM){
	h_eff[idx]->Fill(phistar,weight_temp_h_d);
	h_eff[idx+1]->Fill(phistar,weight_temp_l_d);
	h_eff[idx+2]->Fill(phistar,weight_temp_h_m);
	h_eff[idx+3]->Fill(phistar,weight_temp_l_m);
      }
      else{
	h_eff[idx]->Fill(phistar_true,weight_temp_h_d);
	h_eff[idx+1]->Fill(phistar_true,weight_temp_l_d);
	h_eff[idx+2]->Fill(phistar_true,weight_temp_h_m);
	h_eff[idx+3]->Fill(phistar_true,weight_temp_l_m);
	BinM_eff[idx]         ->Fill(phistar,phistar_true,weight_temp_h_d); 
	BinM_eff[idx+1]       ->Fill(phistar,phistar_true,weight_temp_l_d);
	BinM_eff[idx+2]       ->Fill(phistar,phistar_true,weight_temp_h_m); 
	BinM_eff[idx+3]       ->Fill(phistar,phistar_true,weight_temp_l_m);
      }
      idx=idx+4;
    }
  }
}
void GetBGPhiStar(std::string FileName, double sampleweight, vector<TH1D*> &h_eff, vector<TH1D*> &h_cteq, vector<TH1D*> &h_mstw, vector<TH1D*> &h_nnpdf, vector<TH1D*> &h_fsr_pileup){
  gErrorIgnoreLevel = kError;
  cout<<"reading data for "<<FileName<<endl;
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  nfiles=t->Add(FileName.c_str());

  t->SetBranchStatus("event_info",0); //to disable all branches
  t->SetBranchStatus("truth",0); //to disable all branches
  TBranch *b_reco=t->GetBranch("reco");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_sc=b_reco->GetLeaf("z_phistar_sc");
  TLeaf *l_e0_pt  =b_reco->GetLeaf("e_pt0");
  TLeaf *l_e0_eta =b_reco->GetLeaf("e_eta0");
  TLeaf *l_e1_pt  =b_reco->GetLeaf("e_pt1");
  TLeaf *l_e1_eta =b_reco->GetLeaf("e_eta1");
  int nweights;
  int nwcteq;
  int nwmstw;
  int nwnnpdf;
  t->SetBranchAddress("weight_size",&nweights);
  t->SetBranchAddress("weight_cteq_size",&nwcteq);
  t->SetBranchAddress("weight_mstw_size",&nwmstw);
  t->SetBranchAddress("weight_nnpdf_size",&nwnnpdf);
  t->GetEntry(0);
  double weights[nweights];
  int weightid[nweights];
  double weights_cteq[nwcteq];
  double weights_mstw[nwmstw];
  double weights_nnpdf[nwnnpdf];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);
  t->SetBranchAddress("weights_cteq",&weights_cteq);
  t->SetBranchAddress("weights_mstw",&weights_mstw);
  t->SetBranchAddress("weights_nnpdf",&weights_nnpdf);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr",&weight_fsr);
  // float weight_cteq6;
  // t->SetBranchAddress("weight_cteq6",&weight_cteq6);

  int n_effvariations=5*4*2+1+4*6*2*2;
  h_eff=EmptyHVec(n_effvariations);
  h_cteq=EmptyHVec(nwcteq);
  h_mstw=EmptyHVec(nwmstw);
  h_nnpdf=EmptyHVec(nwnnpdf);
  h_fsr_pileup=EmptyHVec(4);

  cout<<"Entries: "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    double E0_pt  =l_e0_pt  ->GetValue();
    double E0_eta =l_e0_eta ->GetValue();
    double E1_pt  =l_e1_pt  ->GetValue();
    double E1_eta =l_e1_eta ->GetValue();
    double phistar=l_phistar->GetValue();
    double phistar_sc=l_phistar_sc->GetValue();
    double weight =sampleweight;
    //    double pdfnorm=weight_cteq6;
    double pdfnorm=weights_cteq[0];
    double weightpu_0=0;
    double weightpu_p=0;
    double weightpu_m=0;

    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
      if (weightid[w]==2) weightpu_0=weights[w];
      if (weightid[w]==3) weightpu_p=weights[w];
      if (weightid[w]==4) weightpu_m=weights[w];
    }
    if (weightpu_0==0 || weightpu_p==0 || weightpu_m==0) cout<<"pile-up weights not there"<<endl;

    h_eff[0]->Fill(phistar,weight);
    for (int w=0; w<nwcteq; w++) {h_cteq[w] ->Fill(phistar,weight*weights_cteq[w]/pdfnorm);}
    for (int w=0; w<nwmstw; w++) {h_mstw[w] ->Fill(phistar,weight*weights_mstw[w]/pdfnorm);}
    for (int w=0; w<nwnnpdf; w++){h_nnpdf[w]->Fill(phistar,weight*weights_nnpdf[w]/pdfnorm);}
    h_fsr_pileup[0]->Fill(phistar,weight*weight_fsr);
    h_fsr_pileup[1]->Fill(phistar_sc,weight);
    if (weightpu_0!=0){
      h_fsr_pileup[2]->Fill(phistar,weight*weightpu_p/weightpu_0);
      h_fsr_pileup[3]->Fill(phistar,weight*weightpu_m/weightpu_0);
    }
    int idx=1;
    vector<RooUnfoldResponse*> dummy;
    FillRecoEffFluc(idx,h_eff,weight,phistar,E0_pt,E0_eta,E1_pt,E1_eta,0,dummy);
    if (idx==1) cout<<"idx did not increase when it should have"<<endl;
    int binE0=-1;
    int binP0=-1;
    int binE1=-1;
    int binP1=-1;
    for (int e=0; e<6; e++){
      if (E0_eta>TEfficiencyEtaBins[e] && E0_eta<TEfficiencyEtaBins[e+1]) binE0=e;
      if (E1_eta>TEfficiencyEtaBins[e] && E1_eta<TEfficiencyEtaBins[e+1]) binE1=e;
    }
    for (int p=0; p<4; p++){
      if (E0_pt>TEfficiencyETBins[p] && E0_pt<TEfficiencyETBins[p+1]) binP0=p;
      if (E1_pt>TEfficiencyETBins[p] && E1_pt<TEfficiencyETBins[p+1]) binP1=p;
    }
    FillTrigEffFluc(idx,h_eff,weight,phistar,binP0,binE0,binP1,binE1,0,dummy);
  }
  cout<<"done reading data for "<<FileName<<endl;
}

vector<TH1D*> GetDataPhiStar(double sampleweight, vector<TH1D*> &h_phistar_sc){
  cout<<"reading data"<<endl;
  TH2D* BinMigration_sc=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
  BinMigration_sc->Sumw2();
  vector<std::string> File_Data=Get_File_Data();
  vector<TH1D*> h_phistar;
  for (int nt=0; nt<Ntoys+1; nt++){
    TH1D *h_phistar_temp= new TH1D("phistar","phistar",nphistar,phistarBins);
    h_phistar_temp->Sumw2();
    h_phistar.push_back(h_phistar_temp);
    h_phistar_sc.push_back(h_phistar_temp);
  }
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  for (uint i=0; i<File_Data.size(); i++){
    nfiles=t->Add(File_Data[i].c_str());
  }
  TBranch *b_reco=t->GetBranch("reco");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_sc=b_reco->GetLeaf("z_phistar_sc");
  cout<<"Entries: "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    h_phistar[0]->Fill(l_phistar->GetValue(),sampleweight);
    h_phistar_sc[0]->Fill(l_phistar_sc->GetValue(),sampleweight);
    BinMigration_sc->Fill(l_phistar->GetValue(),l_phistar_sc->GetValue(),1);
    for (int nt=0; nt<Ntoys; nt++){
      double x=gRandom->Poisson(1);
      h_phistar[nt+1]->Fill(l_phistar->GetValue(),sampleweight*x);
      h_phistar_sc[nt+1]->Fill(l_phistar_sc->GetValue(),sampleweight*x);

    } 
  }
  cout<<"filled data phistar histogram"<<endl;
  BinMigration_sc->GetXaxis()->SetRangeUser(0.001,3.2);
  BinMigration_sc->GetXaxis()->SetTitle("#phi*_{SC} (reconstructed)");
  BinMigration_sc->GetXaxis()->SetTitleOffset(0.8);
  BinMigration_sc->GetXaxis()->SetTitleSize(0.04);
  BinMigration_sc->GetXaxis()->SetLabelOffset(-0.01);
  BinMigration_sc->GetXaxis()->SetLabelSize(0.04);
  BinMigration_sc->GetYaxis()->SetTitleOffset(1.05);
  BinMigration_sc->GetYaxis()->SetTitleSize(0.04);
  BinMigration_sc->GetYaxis()->SetLabelSize(0.04);
  BinMigration_sc->GetYaxis()->SetRangeUser(0.001,3.2);
  BinMigration_sc->GetYaxis()->SetTitle("#phi* (generated)");
  BinMigration_sc->GetZaxis()->SetTitleOffset(-0.004);
  BinMigration_sc->SetStats(0);
  BinMigration_sc->SetBit( TH2::kNoTitle, true );
  TCanvas* BinM_SC_GSF = new TCanvas("BinM_SC_GSF","BinM_SC_GSF",800,900);
  BinM_SC_GSF->cd();
  BinM_SC_GSF->SetLogx();
  BinM_SC_GSF->SetLogy();
  BinM_SC_GSF->SetLogz();
  BinMigration_sc->Draw("COLZ");
  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);
  mark.DrawLatex(0.19,0.80,"MadGraph");
  return h_phistar;
}

void GetGenPhiStar(double sampleweight, TH1D* &h_phistar, vector<TH1D*> &h_cteq, vector<TH1D*> &h_mstw, vector<TH1D*> &h_nnpdf, vector<TH1D*> &h_fsr_pileup){
  gErrorIgnoreLevel = kError;
  cout<<"reading signal gen"<<endl;
  h_phistar= new TH1D("phistar","phistar",nphistar,phistarBins);
  h_phistar->Sumw2();

  TChain* t = new TChain(gen_name.c_str(),gen_name.c_str());
  int nfiles;
  if (doMG) nfiles=t->Add(File_Signal_gen.c_str());
  else nfiles=t->Add(File_Powheg_gen.c_str());
  t->SetBranchStatus("event_info",0); //to disable all branches
  t->SetBranchStatus("reco",0); //to disable all branches
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar=b_truth->GetLeaf("z_phistar_dressed");
  int nweights;
  int nwcteq;
  int nwmstw;
  int nwnnpdf;
  t->SetBranchAddress("weight_size",&nweights);
  t->SetBranchAddress("weight_cteq_size",&nwcteq);
  t->SetBranchAddress("weight_mstw_size",&nwmstw);
  t->SetBranchAddress("weight_nnpdf_size",&nwnnpdf);
  t->GetEntry(0);

  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr",&weight_fsr);

  cout<<"Entries: "<<t->GetEntries()<<endl;

  h_cteq=EmptyHVec(nwcteq);
  h_mstw=EmptyHVec(nwmstw);
  h_nnpdf=EmptyHVec(nwnnpdf);
  h_fsr_pileup=EmptyHVec(4);
 
  for (int i=0; i<t->GetEntries();i++){
  //  for (int i=0; i<50000;i++){
    if (nMC!=-1 && i>nMC) continue;
    t->GetEntry(i);
    double weight=sampleweight;
     double weightpu_0=0;
    double weightpu_p=0;
    double weightpu_m=0;
    //cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;

     for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2) {weight=weight*weights[w];}
      if (weightid[w]==2) weightpu_0=weights[w];
      if (weightid[w]==3) weightpu_p=weights[w];
      if (weightid[w]==4) weightpu_m=weights[w];
    }
    if (weightpu_0==0 || weightpu_p==0 || weightpu_m==0) cout<<"pile-up weights not there"<<endl;
    h_phistar->Fill(l_phistar->GetValue(),weight);
    h_fsr_pileup[0]->Fill(l_phistar->GetValue(),weight*weight_fsr);
    h_fsr_pileup[1]->Fill(l_phistar->GetValue(),weight);
    if (weightpu_0!=0){
      h_fsr_pileup[2]->Fill(l_phistar->GetValue(),weight*weightpu_p/weightpu_0);
      h_fsr_pileup[3]->Fill(l_phistar->GetValue(),weight*weightpu_m/weightpu_0);
    }
  }
  cout<<"done reading signal gen"<<endl;
  return;
}

void GetBinM(double sampleweight, vector<RooUnfoldResponse*> &BinM_eff, vector<RooUnfoldResponse*> &BinM_mcstat, vector<RooUnfoldResponse*> &BinM_sc, vector<RooUnfoldResponse*> &BinM_cteq, vector<RooUnfoldResponse*> &BinM_mstw, vector<RooUnfoldResponse*> &BinM_nnpdf, vector<RooUnfoldResponse*> &BinM_fsr_pileup, vector<TH1D*> &h_eff, vector<TH1D*> &h_cteq, vector<TH1D*> &h_mstw, vector<TH1D*> &h_nnpdf, vector<TH1D*> &h_fsr_pileup){
  gErrorIgnoreLevel = kError;
  cout<<"reading data for "<<File_Signal_reco<<"  "<<reco_name<<endl;
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  if (doMG) nfiles=t->Add(File_Signal_reco.c_str());
  else nfiles=t->Add(File_Powheg_reco.c_str());

  t->SetBranchStatus("event_info",0); //to disable all branches
  TBranch *b_reco=t->GetBranch("reco");
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_sc=b_reco->GetLeaf("z_phistar_sc");
  TLeaf *l_e0_pt  =b_reco->GetLeaf("e_pt0");
  TLeaf *l_e0_eta =b_reco->GetLeaf("e_eta0");
  TLeaf *l_e1_pt  =b_reco->GetLeaf("e_pt1");
  TLeaf *l_e1_eta =b_reco->GetLeaf("e_eta1");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  int nweights;
  int nwcteq;
  int nwmstw;
  int nwnnpdf;
  t->SetBranchAddress("weight_size",&nweights);
  t->SetBranchAddress("weight_cteq_size",&nwcteq);
  t->SetBranchAddress("weight_mstw_size",&nwmstw);
  t->SetBranchAddress("weight_nnpdf_size",&nwnnpdf);
  t->GetEntry(0);

  double weights[nweights];
  int weightid[nweights];
   t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr",&weight_fsr);

  cout<<"Entries: "<<t->GetEntries()<<endl;

  int n_effvariations=5*4*2+1+4*6*2*2;
  h_eff=EmptyHVec(n_effvariations);
  h_cteq=EmptyHVec(nwcteq);
  h_mstw=EmptyHVec(nwmstw);
  h_nnpdf=EmptyHVec(nwnnpdf);
  h_fsr_pileup=EmptyHVec(4);
  BinM_eff=EmptyBinMVec(n_effvariations);
  BinM_cteq=EmptyBinMVec(nwcteq);
  BinM_mstw=EmptyBinMVec(nwmstw);
  BinM_nnpdf=EmptyBinMVec(nwnnpdf);
  BinM_fsr_pileup=EmptyBinMVec(4);

  TH1D *h_rec= new TH1D("phistar","phistar",nphistar,phistarBins);
  h_rec->Sumw2();
  TH1D *h_rec_sc= new TH1D("phistar","phistar",nphistar,phistarBins);
  h_rec_sc->Sumw2();
  vector<TH2D*> BinMigration, BinMigration_sc;
  for (int i=0; i<Ntoys+1; i++){
    TH2D* BinMigration_temp=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
    BinMigration_temp->Sumw2();
    TH2D* BinMigration_sc_temp=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
    BinMigration_sc_temp->Sumw2();
    BinMigration.push_back(BinMigration_temp);
    BinMigration_sc.push_back(BinMigration_sc_temp);
  }

  for (int i=0; i<t->GetEntries();i++){
   if (nMC!=-1 && i>nMC) continue;
  //for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double E0_pt  =l_e0_pt  ->GetValue();
    double E0_eta =l_e0_eta ->GetValue();
    double E1_pt  =l_e1_pt  ->GetValue();
    double E1_eta =l_e1_eta ->GetValue();
    double phistar=l_phistar->GetValue();
    double phistar_sc=l_phistar_sc->GetValue();
    double phistar_true=l_phistar_true->GetValue();
    double pdfnorm=0;
    double weight =sampleweight;
    double weightpu_0=0;
    double weightpu_p=0;
    double weightpu_m=0;
    //    cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;
   
     for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
      if (weightid[w]==2) weightpu_0=weights[w];
      if (weightid[w]==3) weightpu_p=weights[w];
      if (weightid[w]==4) weightpu_m=weights[w];
    }
    if (weightpu_0==0 || weightpu_p==0 || weightpu_m==0) cout<<"pile-up weights not there"<<endl;
    if (weightpu_0==0) cout<<"pile-up nominal weight 0"<<endl;
    h_eff[0]      ->Fill(phistar_true,weight);
    BinM_eff[0]   ->Fill(phistar,phistar_true,weight); 
    h_rec         ->Fill(phistar,weight);
    h_rec_sc          ->Fill(phistar_sc,weight);
    BinMigration[0]  ->Fill(phistar,phistar_true,weight);
    BinMigration_sc[0]   ->Fill(phistar_sc,phistar_true,weight);
    double x;
    for (int nt=0; nt<Ntoys; nt++){
      x=gRandom->Poisson(1);
      BinMigration[nt+1]      ->Fill(phistar,phistar_true,x*weight);
      BinMigration_sc[nt+1]   ->Fill(phistar_sc,phistar_true,x*weight);
    } 
    h_fsr_pileup[0]   ->Fill(phistar_true,weight*weight_fsr);
    BinM_fsr_pileup[0]->Fill(phistar,phistar_true,weight*weight_fsr);
    h_fsr_pileup[1]   ->Fill(phistar_true,weight);
    BinM_fsr_pileup[1]->Fill(phistar_sc,phistar_true,weight);
    if (weightpu_0!=0){
      h_fsr_pileup[2]   ->Fill(phistar_true,weight*weightpu_p/weightpu_0);
      BinM_fsr_pileup[2]->Fill(phistar,phistar_true,weight*weightpu_p/weightpu_0);
      h_fsr_pileup[3]   ->Fill(phistar_true,weight*weightpu_m/weightpu_0);
      BinM_fsr_pileup[3]->Fill(phistar,phistar_true,weight*weightpu_m/weightpu_0);
    }
    int idx=1;
    FillRecoEffFluc(idx,h_eff,weight,phistar,E0_pt,E0_eta,E1_pt,E1_eta,1,BinM_eff,phistar_true);
    int binE0=-1;
    int binP0=-1;
    int binE1=-1;
    int binP1=-1;
    for (int e=0; e<6; e++){
      if (E0_eta>TEfficiencyEtaBins[e] && E0_eta<TEfficiencyEtaBins[e+1]) binE0=e;
      if (E1_eta>TEfficiencyEtaBins[e] && E1_eta<TEfficiencyEtaBins[e+1]) binE1=e;
    }
    for (int p=0; p<4; p++){
      if (E0_pt>TEfficiencyETBins[p] && E0_pt<TEfficiencyETBins[p+1]) binP0=p;
      if (E1_pt>TEfficiencyETBins[p] && E1_pt<TEfficiencyETBins[p+1]) binP1=p;
    }
    FillTrigEffFluc(idx,h_eff,weight,phistar,binP0,binE0,binP1,binE1,1,BinM_eff,phistar_true);
  }

  BinMigration_sc[0]->GetXaxis()->SetRangeUser(0.001,3.2);
  BinMigration_sc[0]->GetXaxis()->SetTitle("#phi*_{SC} (reconstructed)");
  BinMigration_sc[0]->GetXaxis()->SetTitleOffset(0.8);
  BinMigration_sc[0]->GetXaxis()->SetTitleSize(0.04);
  BinMigration_sc[0]->GetXaxis()->SetLabelOffset(-0.01);
  BinMigration_sc[0]->GetXaxis()->SetLabelSize(0.04);
  BinMigration_sc[0]->GetYaxis()->SetTitleOffset(1.05);
  BinMigration_sc[0]->GetYaxis()->SetTitleSize(0.04);
  BinMigration_sc[0]->GetYaxis()->SetLabelSize(0.04);
  BinMigration_sc[0]->GetYaxis()->SetRangeUser(0.001,3.2);
  BinMigration_sc[0]->GetYaxis()->SetTitle("#phi* (generated)");
  BinMigration_sc[0]->GetZaxis()->SetTitleOffset(-0.004);
  BinMigration_sc[0]->SetStats(0);
  BinMigration_sc[0]->SetBit( TH2::kNoTitle, true );
  TCanvas* BinM_SC = new TCanvas("BinM_SC","BinM_SC",800,900);
  BinM_SC->cd();
  BinM_SC->SetLogx();
  BinM_SC->SetLogy();
  BinM_SC->SetLogz();
  BinMigration_sc[0]->Draw("COLZ");
  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);
  mark.DrawLatex(0.19,0.80,"MadGraph");
 
  GetToyResponse(BinM_mcstat, BinMigration);
  GetToyResponse(BinM_sc, BinMigration_sc);
  cout<<"done reading data for "<<File_Signal_reco<<"  "<<reco_name<<endl;
  return;
}

void NTupleZShape(){
  double Lumi = 19700.;
  double ttbar_weight  = 1./(4246440./23.64);
  double tautau_weight = 1./(47271600./1966.7);
  double tbarw_weight  = 1./(493460./11.1);
  double tw_weight     = 1./(497658./11.1);
  double ww_weight     = 1./(10000430./54.84);
  double wz_weight     = 1./(10000280./33.21);
  double zz_weight     = 1./(9799908./17.0);
  double signal_weight = 1./(30459500./3531.89);
  if (!doMG) signal_weight = 1./(3297045./1966.7);

  vector<TH1D*> Data_sc;
  vector<TH1D*> Data=GetDataPhiStar(1./Lumi, Data_sc);
  vector<TH1D*> h_tt_eff,     h_tt_cteq,     h_tt_mstw,     h_tt_nnpdf,     h_tt_fsr_pileup;
  vector<TH1D*> h_tautau_eff, h_tautau_cteq, h_tautau_mstw, h_tautau_nnpdf, h_tautau_fsr_pileup;
  vector<TH1D*> h_tbarw_eff,  h_tbarw_cteq,  h_tbarw_mstw,  h_tbarw_nnpdf,  h_tbarw_fsr_pileup;
  vector<TH1D*> h_tw_eff,     h_tw_cteq,     h_tw_mstw,     h_tw_nnpdf,     h_tw_fsr_pileup;
  vector<TH1D*> h_ww_eff,     h_ww_cteq,     h_ww_mstw,     h_ww_nnpdf,     h_ww_fsr_pileup;
  vector<TH1D*> h_wz_eff,     h_wz_cteq,     h_wz_mstw,     h_wz_nnpdf,     h_wz_fsr_pileup;
  vector<TH1D*> h_zz_eff,     h_zz_cteq,     h_zz_mstw,     h_zz_nnpdf,     h_zz_fsr_pileup;

  GetBGPhiStar(File_tt,    ttbar_weight, h_tt_eff,     h_tt_cteq,     h_tt_mstw,     h_tt_nnpdf,     h_tt_fsr_pileup); 
  GetBGPhiStar(File_tautau,tautau_weight,h_tautau_eff, h_tautau_cteq, h_tautau_mstw, h_tautau_nnpdf, h_tautau_fsr_pileup);
  GetBGPhiStar(File_tbarw, tbarw_weight, h_tbarw_eff,  h_tbarw_cteq,  h_tbarw_mstw,  h_tbarw_nnpdf,  h_tbarw_fsr_pileup);
  GetBGPhiStar(File_tw,    tw_weight,    h_tw_eff,     h_tw_cteq,     h_tw_mstw,     h_tw_nnpdf,     h_tw_fsr_pileup);
  GetBGPhiStar(File_ww,    ww_weight,    h_ww_eff,     h_ww_cteq,     h_ww_mstw,     h_ww_nnpdf,     h_ww_fsr_pileup);
  GetBGPhiStar(File_wz,    wz_weight,    h_wz_eff,     h_wz_cteq,     h_wz_mstw,     h_wz_nnpdf,     h_wz_fsr_pileup);
  GetBGPhiStar(File_zz,    zz_weight,    h_zz_eff,     h_zz_cteq,     h_zz_mstw,     h_zz_nnpdf,     h_zz_fsr_pileup);
  
  // PrintBG(Data,h_tt_eff[0]  ,h_tautau_eff[0]  ,h_tbarw_eff[0]  ,h_tw_eff[0]  ,h_ww_eff[0]  ,h_wz_eff[0]  ,h_zz_eff[0])

  // TFile f_bg("data_mc_emu_ratio.root");
  // TCanvas *c_bg = (TCanvas*)f_bg.Get("Canvas_1");
  // TH1F *bg_sf_full = (TH1F*)c_bg->FindObject("hPull");
  // bg_sf_full->Sumw2();  

  TH1D* bg_sf=new TH1D("bg","bg",nphistar,phistarBins);
  for (uint i=0; i<nphistar; i++){
    bg_sf->SetBinContent(i+1,1.);
    bg_sf->SetBinError(i+1,0.);
  }

  vector<TH1D*> h_data_eff  =RemoveBG(Data[0],bg_sf,h_tt_eff  ,h_tautau_eff  ,h_tbarw_eff  ,h_tw_eff  ,h_ww_eff  ,h_wz_eff  ,h_zz_eff);
  vector<TH1D*> h_data_fsr_pileup; 
  vector<TH1D*> h_data_stat; 
  vector<TH1D*> h_data_sc_stat; 
 
  for (uint idx=0; idx<h_tt_fsr_pileup.size();idx++){
    TH1D* bgtemp=(TH1D*)h_tt_fsr_pileup[idx]->Clone();
    bgtemp->Add(h_tautau_fsr_pileup[idx],1.0);
    bgtemp->Add(h_tbarw_fsr_pileup[idx],1.0);
    bgtemp->Add(h_tw_fsr_pileup[idx],1.0);
    bgtemp->Add(h_ww_fsr_pileup[idx],1.0);
    bgtemp->Multiply(bg_sf);
    TH1D* datatemp  =(TH1D*)Data[0]->Clone();
    if (idx==1) datatemp  =(TH1D*)Data_sc[0]->Clone();
    datatemp->Add(bgtemp,-1.0);
    datatemp->Add(h_wz_fsr_pileup[idx],-1.0);
    datatemp->Add(h_zz_fsr_pileup[idx],-1.0);
    h_data_fsr_pileup.push_back(datatemp);
  }

  for (uint idx=0; idx<Ntoys+1;idx++){
    TH1D* bgtemp=(TH1D*)h_tt_eff[0]->Clone();
    bgtemp->Add(h_tautau_eff[0],1.0);
    bgtemp->Add(h_tbarw_eff[0],1.0);
    bgtemp->Add(h_tw_eff[0],1.0);
    bgtemp->Add(h_ww_eff[0],1.0);
    bgtemp->Multiply(bg_sf);
    TH1D* datatemp  =(TH1D*)Data[idx]->Clone();
    datatemp->Add(bgtemp,-1.0);
    datatemp->Add(h_wz_eff[0],-1.0);
    datatemp->Add(h_zz_eff[0],-1.0);
    h_data_stat.push_back(datatemp);
  }

  for (uint idx=0; idx<Ntoys+1;idx++){
    TH1D* bgtemp=(TH1D*)h_tt_eff[0]->Clone();
    bgtemp->Add(h_tautau_eff[0],1.0);
    bgtemp->Add(h_tbarw_eff[0],1.0);
    bgtemp->Add(h_tw_eff[0],1.0);
    bgtemp->Add(h_ww_eff[0],1.0);
    bgtemp->Multiply(bg_sf);
    TH1D* datatemp  =(TH1D*)Data_sc[idx]->Clone();
    datatemp->Add(bgtemp,-1.0);
    datatemp->Add(h_wz_eff[0],-1.0);
    datatemp->Add(h_zz_eff[0],-1.0);
    h_data_sc_stat.push_back(datatemp);
  }

  vector<RooUnfoldResponse*> BinM_eff, BinM_mcstat, BinM_sc, BinM_cteq, BinM_mstw, BinM_nnpdf, BinM_fsr_pileup;
  vector<TH1D*> mc_truereco_eff, mc_truereco_cteq, mc_truereco_mstw, mc_truereco_nnpdf, mc_truereco_fsr_pileup;
  GetBinM(signal_weight, BinM_eff, BinM_mcstat, BinM_sc, BinM_cteq, BinM_mstw, BinM_nnpdf, BinM_fsr_pileup, mc_truereco_eff, mc_truereco_cteq, mc_truereco_mstw, mc_truereco_nnpdf, mc_truereco_fsr_pileup);
  
  TH1D* mc_truegen;
  vector<TH1D*> mc_truegen_cteq, mc_truegen_mstw, mc_truegen_nnpdf, mc_truegen_fsr_pileup;
  GetGenPhiStar(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_mstw, mc_truegen_nnpdf, mc_truegen_fsr_pileup);
  
  vector<TH1D*> h_eff_eff  = GetEff(mc_truereco_eff  ,mc_truegen);
  
  vector<TGraphAsymmErrors *> g_data_phistar_mcstat =GetUnfoldedData(BinM_mcstat,h_data_eff[0],h_eff_eff[0]);
  vector<TGraphAsymmErrors *> g_data_phistar_sc     =GetUnfoldedData(BinM_sc    ,h_data_fsr_pileup[1],h_eff_eff[0]);
  vector<TGraphAsymmErrors *> g_data_phistar_stat   =GetUnfoldedData(BinM_mcstat[0],h_data_stat,h_eff_eff[0]);
  vector<TGraphAsymmErrors *> g_data_phistar_sc_stat=GetUnfoldedData(BinM_sc    [0],h_data_sc_stat,h_eff_eff[0]);

  NormalizeGraph(g_data_phistar_mcstat);
  NormalizeGraph(g_data_phistar_sc);
  NormalizeGraph(g_data_phistar_stat);
  NormalizeGraph(g_data_phistar_sc_stat);
 
  vector<TGraphAsymmErrors *> g_ratio_phistar_mcstat;
  vector<TGraphAsymmErrors *> g_ratio_phistar_scstat;
  for (int t=0; t<Ntoys+1; t++){
    TGraphAsymmErrors* g_ratiosc_phistar      = CreateRatio(g_data_phistar_sc[t],g_data_phistar_mcstat[t],1);
    TGraphAsymmErrors* g_ratiosc_phistar_stat = CreateRatio(g_data_phistar_sc_stat[t],g_data_phistar_stat[t],1);
    g_ratio_phistar_scstat.push_back(g_ratiosc_phistar_stat);
    g_ratio_phistar_mcstat.push_back(g_ratiosc_phistar);
  }

  TGraphAsymmErrors* g_syst_phistar_sc     = CalcTotalSysU_nnpdf(g_ratio_phistar_mcstat,1);
  TGraphAsymmErrors* g_syst_phistar_stat   = CalcTotalSysU_nnpdf(g_ratio_phistar_scstat,1);

  for (uint i=0; i<nphistar-1; i++){
    double errorh=g_syst_phistar_sc->GetErrorYhigh(i);
    double errorh2=g_syst_phistar_stat->GetErrorYhigh(i);
    double errorl=g_syst_phistar_sc->GetErrorYlow(i);
    double errorl2=g_syst_phistar_stat->GetErrorYlow(i);
    cout<<errorh<<"  "<<errorh2<<"  "<<errorl<<"  "<<errorl2<<endl;
    errorh=sqrt(errorh*errorh+errorh2*errorh2);
    errorl=sqrt(errorl*errorl+errorl2*errorl2);
    g_syst_phistar_sc->SetPointError(i,0,0,errorl,errorh);
  }

  TH1D* Pull2=new TH1D("Pull2","Pull2",20,-5,5);
  Pull2->Sumw2();
  for (uint i=0; i<nphistar-1; i++){
    double x,y;
    g_syst_phistar_sc->GetPoint(i,x,y);
    double error;
    if (y>1) error=g_syst_phistar_sc->GetErrorYlow(i);
    else error=g_syst_phistar_sc->GetErrorYhigh(i);
    double pull=(y-1)/error;
    Pull2->Fill(pull, 1);
  }

  TCanvas* CPull1 = new TCanvas("CPull1","CPull1",800,900);
  CPull1->cd();
  Pull2->GetXaxis()->SetTitle("(SC/GSF-1)/#sigma");
  Pull2->GetXaxis()->SetTitleOffset(1.0);
  Pull2->GetXaxis()->SetLabelOffset(0.0);
  Pull2->GetXaxis()->SetTitleSize(0.04);
  Pull2->GetXaxis()->SetLabelSize(0.04);
  Pull2->SetBit( TH1::kNoTitle, true );
  Pull2->SetStats(0);
  Pull2->SetLineColor(1);
  Pull2->SetMarkerColor(1);
  Pull2->GetYaxis()->SetTitle("#phi* bins");
  Pull2->GetYaxis()->SetTitleOffset(1.0);
  Pull2->GetYaxis()->SetTitleSize(0.04);
  Pull2->GetYaxis()->SetLabelSize(0.04);
  Pull2->Draw();
  TLegend* leg = new TLegend(0.53,0.80,0.90,0.86);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineWidth(1);
  leg->SetNColumns(1);
  leg->SetTextFont(42);
  leg->AddEntry(Pull2,"2012 data","PL");
  leg->Draw();

  TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot","FinalPhiTot",800,900);
  FinalPhiTot->cd();
  FinalPhiTot->SetLogx();
  g_syst_phistar_sc->GetXaxis()->SetRangeUser(0.001,3.0);
  g_syst_phistar_sc->GetYaxis()->SetRangeUser(0.8,1.2);
  g_syst_phistar_sc->GetXaxis()->SetTitle("#phi*");
  g_syst_phistar_sc->GetXaxis()->SetTitleOffset(0.8);
  g_syst_phistar_sc->GetXaxis()->SetTitleSize(0.04);
  g_syst_phistar_sc->GetXaxis()->SetLabelOffset(-0.01);
  g_syst_phistar_sc->GetXaxis()->SetLabelSize(0.04);
  g_syst_phistar_sc->GetYaxis()->SetTitle("SC/GSF");
  g_syst_phistar_sc->GetYaxis()->SetTitleOffset(1.2);
  g_syst_phistar_sc->GetYaxis()->SetTitleSize(0.04);
  g_syst_phistar_sc->GetYaxis()->SetLabelSize(0.03);
  g_syst_phistar_sc->GetYaxis()->SetRangeUser(0.95,1.05);
  g_syst_phistar_sc->SetMarkerStyle(20);
  g_syst_phistar_sc->Draw("APE");
  g_syst_phistar_stat->SetMarkerStyle(22);
  g_syst_phistar_stat->SetMarkerColor(3);
  g_syst_phistar_stat->SetLineColor(3);
  g_syst_phistar_stat->Draw("PEsame");
  g_ratio_phistar_mcstat[0]->SetMarkerStyle(21);
  g_ratio_phistar_mcstat[0]->SetMarkerColor(2);
  g_ratio_phistar_mcstat[0]->SetLineColor(2);
  g_ratio_phistar_mcstat[0]->Draw("PEsame");
  TLegend* leg2 = new TLegend(0.53,0.80,0.90,0.86);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetLineWidth(1);
  leg2->SetNColumns(1);
  leg2->SetTextFont(42);
  leg2->AddEntry(g_syst_phistar_sc,"2012 data","PL");
  leg2->Draw();
}
