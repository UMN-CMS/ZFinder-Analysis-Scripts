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
using namespace RooFit;
//using namespace RooStats;

const bool doNorm=0;
const int Ntoys=500;
const int elec=0;
const int doMG=1;

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/MadGraph_gen.root";
std::string File_Powheg_gen= "/afs/cern.ch/work/r/ruckstuh/Powheg_gen.root";

std::string gen_name ="Combined Gen Cuts Reco";

void NormalizeGraph(vector<TGraphAsymmErrors*> &graph){
  for (size_t i=0; i<graph.size();i++){
    double xstot=0;
    double x,y,errorl,errorh;
    for (size_t iphistar=0; iphistar<nphistar-1; iphistar++){
      graph[i]->GetPoint(iphistar,x,y);
      xstot+=y;
    }
    cout<<i<<":Normalisation:"<<xstot<<endl;
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


TGraphAsymmErrors* GetMCFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  uint nun=graph.size();
  ofstream outputfile;
  std::string textname="TableMC_UnTestGen_";
  if (doNorm) textname+="Norm_";
  else        textname+="Abs_";
  if (doMG)   textname+="MG_";
  else        textname+="PH_";
  if (elec==0)textname+="Dressed.txt";
  if (elec==1)textname+="Born.txt";
  if (elec==2)textname+="Naked.txt";
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
  std::string textname="Table_Values_TestGen_";
  if (doNorm) textname+="Norm_";
  else        textname+="Abs_";
  if (doMG)   textname+="MG_";
  else        textname+="PH_";
  if (elec==0)textname+="Dressed.txt";
  if (elec==1)textname+="Born.txt";
  if (elec==2)textname+="Naked.txt";
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

vector<TH1D*> EmptyHVec(int n){  
  vector<TH1D*> h;
  for (int i=0; i<n; i++){
    TH1D *phistartemp= new TH1D("phistar","phistar",nphistar,phistarBins);
    phistartemp->Sumw2();
    h.push_back(phistartemp);
  }
  return h;
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
  if (elec==1) l_phistar=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar=b_truth->GetLeaf("z_phistar_naked");
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
  float weight_cteq6;
  if (doMG) t->SetBranchAddress("weight_cteq6",&weight_cteq6);

  cout<<"Entries: "<<t->GetEntries()<<endl;

  h_cteq=EmptyHVec(nwcteq);
  h_mstw=EmptyHVec(nwmstw);
  h_nnpdf=EmptyHVec(nwnnpdf);
  h_fsr_pileup=EmptyHVec(4);
 
  for (int i=0; i<t->GetEntries();i++){
  //  for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double weight=sampleweight;
    double pdfnorm=weights_cteq[0];
    if (doMG) pdfnorm=weight_cteq6;
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
    if (pdfnorm!=0){    
      for (int w=0; w<nwcteq; w++) {h_cteq[w] ->Fill(l_phistar->GetValue(),weight*weights_cteq[w]/pdfnorm);}
      for (int w=0; w<nwmstw; w++) {h_mstw[w] ->Fill(l_phistar->GetValue(),weight*weights_mstw[w]/pdfnorm);}
      for (int w=0; w<nwnnpdf; w++){h_nnpdf[w]->Fill(l_phistar->GetValue(),weight*weights_nnpdf[w]/pdfnorm);}
    }
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


void NTupleZShape(){
  double signal_weight = 3531.89/30459500.;//3504
  if (!doMG) signal_weight = 1966.7/3297045.;

  TH1D* mc_truegen;
  vector<TH1D*> mc_truegen_cteq, mc_truegen_mstw, mc_truegen_nnpdf, mc_truegen_fsr_pileup;
  GetGenPhiStar(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_mstw, mc_truegen_nnpdf, mc_truegen_fsr_pileup);
  
  cout<<"Get MC distribution"<<endl;
  vector<TGraphAsymmErrors*> g_MC_phistar;
  TGraphAsymmErrors* g_MC_phistar_temp1=ConvertToTGraph(mc_truegen);
  g_MC_phistar.push_back(g_MC_phistar_temp1);
  
  vector<TGraphAsymmErrors*> g_MC_phistar_cteq;
  for (size_t u=0; u<mc_truegen_cteq.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_cteq[u]);
    g_MC_phistar_cteq.push_back(g_MC_phistar_temp);
  }
  vector<TGraphAsymmErrors*> g_MC_phistar_mstw;
  for (size_t u=0; u<mc_truegen_mstw.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_mstw[u]);
    g_MC_phistar_mstw.push_back(g_MC_phistar_temp);
  }
  vector<TGraphAsymmErrors*> g_MC_phistar_nnpdf;
  for (size_t u=0; u<mc_truegen_nnpdf.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_nnpdf[u]);
    g_MC_phistar_nnpdf.push_back(g_MC_phistar_temp);
  }
  vector<TGraphAsymmErrors*> g_MC_phistar_fsr_pileup;
  for (size_t u=0; u<mc_truegen_fsr_pileup.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_fsr_pileup[u]);
    g_MC_phistar_fsr_pileup.push_back(g_MC_phistar_temp);
  }

  NormalizeGraph(g_MC_phistar);
  NormalizeGraph(g_MC_phistar_cteq);
  NormalizeGraph(g_MC_phistar_mstw);
  NormalizeGraph(g_MC_phistar_nnpdf);
  NormalizeGraph(g_MC_phistar_fsr_pileup);
 
  TGraphAsymmErrors* g_syst_mc_cteq = CalcTotalSysU(g_MC_phistar_cteq,g_MC_phistar_cteq[0],1,1);
  TGraphAsymmErrors* g_syst_mc_mstw = CalcTotalSysU(g_MC_phistar_mstw,g_MC_phistar_mstw[0],1);
  TGraphAsymmErrors* g_syst_mc_nnpdf= CalcTotalSysU_nnpdf(g_MC_phistar_nnpdf);
  TGraphAsymmErrors* g_syst_mc_fsr  = CalcTotalSysU_fsr(g_MC_phistar_fsr_pileup[0],g_MC_phistar[0]);
  TGraphAsymmErrors* g_syst_mc_pu   = CalcTotalSysU_pileup(g_MC_phistar_fsr_pileup[2],g_MC_phistar_fsr_pileup[3],g_MC_phistar[0]);
  
  vector<TGraphAsymmErrors *> g_mc_syst;
  g_mc_syst.push_back(g_MC_phistar[0]);
  g_mc_syst.push_back(g_syst_mc_fsr);
  g_mc_syst.push_back(g_syst_mc_pu);
  g_mc_syst.push_back(g_syst_mc_cteq);
  g_mc_syst.push_back(g_syst_mc_mstw);
  g_mc_syst.push_back(g_syst_mc_nnpdf);
  vector<std::string> syst_mc_list;
  syst_mc_list.push_back("unfolding");
  syst_mc_list.push_back("fsr");
  syst_mc_list.push_back("pileup");
  syst_mc_list.push_back("cteq");
  syst_mc_list.push_back("mstw");
  syst_mc_list.push_back("nnpdf");
  
  TGraphAsymmErrors* g_mc_final=GetMCFinal(g_mc_syst, syst_mc_list);
  
  TH1D* h_mc= ConvertToHist(g_mc_final,0);

  std::string textn="Final_Hist_TestGen_";
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  TFile tr(textn.c_str(),"RECREATE");
  h_mc->Write();

  PrintFinal(g_mc_final,g_mc_final);
}
