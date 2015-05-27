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
const int elec=0;
const int doMG=1;
const std::string Tag="MC_";

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Powheg_gen= "/afs/cern.ch/work/r/ruckstuh/public/Powheg_gen.root";
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


TGraphAsymmErrors* GetMCFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist, bool doLumi=1){
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
  uint nun=graph.size();
  ofstream outputfile;
  std::string textname="TableMC_Un_";
  textname+=Tag;
  if (doNorm) textname+="Norm_";
  else        textname+="Abs_";
  if (doMG)   textname+="MG_";
  else        textname+="PH_";
  if (elec==0)textname+="Dressed.txt";
  if (elec==1)textname+="Born.txt";
  if (elec==2)textname+="Naked.txt";
  outputfile.open(textname.c_str());
  outputfile << "$\\phi^*$ range & Total & Stat.";
  if (!doMG || doNorm) outputfile << " & PDF";
  outputfile <<" & Pile-up & FSR \\\\ \\hline"<< "\n";
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double e_un=0;
    double e_pdf=0;
    double e_fsr=0;
    double e_pu=0;
    double x,y,errorl,errorh;
    graph[0]->GetPoint(iphistar,x,y);
    double x_nom=x;
    double y_nom=y;;
    double error_sys_max=0;
    double error_sys_min=0;
    if (doMG && !doNorm && doLumi) error_sys_max=0.033*0.033;
    if (doMG && !doNorm && doLumi) error_sys_min=0.033*0.033;
    for (uint i=0; i<nun; i++){
      graph[i]->GetPoint(iphistar,x,y);
      errorl=graph[i]->GetErrorYlow(iphistar);
      errorh=graph[i]->GetErrorYhigh(iphistar);
      double error_max=errorh/y;
      double error_min=errorl/y;
      error_sys_max+=error_max*error_max;
      error_sys_min+=error_min*error_min;
      if (slist[i]=="unfolding") {e_un =error_max; x_nom=x; y_nom=y;}
      if (slist[i]=="fsr")       {e_fsr=error_max;}
      if (slist[i]=="pileup")    {e_pu =error_max;}
      if (slist[i]=="cteq")      {e_pdf=error_max;}
       //cout<<error_sys_min<<endl;
    }
    //cout<<error_sys_min<<endl;
    g_un->SetPoint(iphistar,x_nom,y_nom);
    g_un->SetPointError(iphistar,0,0,y_nom*sqrt(error_sys_min),y_nom*sqrt(error_sys_max));
    printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f : %9.2f%% : %9.2f : %9.2f%% : %9.2f%% \n",int(iphistar),y_nom,sqrt(error_sys_max)*y_nom, sqrt(error_sys_min)*y_nom,sqrt(error_sys_max)*100.,e_un*100., e_pdf*100., e_fsr*100., e_pu*100.);
    
    outputfile << std::fixed;
    outputfile << std::setprecision(3) <<phistarBins[iphistar]<<"-"<<phistarBins[iphistar+1]<<" & "<< std::setprecision(2)<<sqrt(error_sys_max)*100.<<"\\% & "<<e_un*100.<<"\\% & ";
    if (!doMG) outputfile <<e_pdf*100.<<"\\% & ";
    else if (!doNorm) outputfile <<3.3<<"\\% & ";
    outputfile <<e_pu*100.<<"\\% & "<<e_fsr*100.<<"\\%  \\\\ \\hline"<< "\n";
  }
  outputfile.close();
  cout<<"GetMCFinal"<<endl;
  return g_un;
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


void GetGenPhiStar(double sampleweight, TH1D* &h_phistar, vector<TH1D*> &h_cteq, vector<TH1D*> &h_fsr_pileup){
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
  t->SetBranchAddress("weight_size",&nweights);
  t->SetBranchAddress("weight_cteq_size",&nwcteq);
  t->GetEntry(0);

  double weights[nweights];
  int weightid[nweights];
  double weights_cteq[nwcteq];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);
  t->SetBranchAddress("weights_cteq",&weights_cteq);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr",&weight_fsr);

  cout<<"Entries: "<<t->GetEntries()<<endl;

  h_cteq=EmptyHVec(nwcteq);
  h_fsr_pileup=EmptyHVec(4);
 
  for (int i=0; i<t->GetEntries();i++){
  //  for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double weight=sampleweight;
    double pdfnorm=weights_cteq[0];
     //cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;

     for (int w=0; w<nweights;w++){
       //if (weightid[w]==1 || weightid[w]==2) {weight=weight*weights[w];}
       if (weightid[w]==1) {weight=weight*weights[w];}
     }
     h_phistar->Fill(l_phistar->GetValue(),weight);
    if (pdfnorm!=0){    
      for (int w=0; w<nwcteq; w++) {h_cteq[w] ->Fill(l_phistar->GetValue(),weight*weights_cteq[w]/pdfnorm);}
    }
    h_fsr_pileup[0]->Fill(l_phistar->GetValue(),weight*weight_fsr);
    h_fsr_pileup[1]->Fill(l_phistar->GetValue(),weight);
   }
  cout<<"done reading signal gen"<<endl;
  return;
}


void NTupleZShape_MC(){
  double signal_weight = 3531.89/30459500.;//3504
  if (!doMG) signal_weight = 1966.7/3297045.;

  TH1D* mc_truegen;
  vector<TH1D*> mc_truegen_cteq, mc_truegen_fsr_pileup;
  GetGenPhiStar(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_fsr_pileup);
  
  cout<<"Get MC distribution"<<endl;
  vector<TGraphAsymmErrors*> g_MC_phistar;
  TGraphAsymmErrors* g_MC_phistar_temp1=ConvertToTGraph(mc_truegen);
  g_MC_phistar.push_back(g_MC_phistar_temp1);
  
  vector<TGraphAsymmErrors*> g_MC_phistar_cteq;
  for (size_t u=0; u<mc_truegen_cteq.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_cteq[u]);
    g_MC_phistar_cteq.push_back(g_MC_phistar_temp);
  }
  vector<TGraphAsymmErrors*> g_MC_phistar_fsr_pileup;
  for (size_t u=0; u<mc_truegen_fsr_pileup.size(); u++){
    TGraphAsymmErrors* g_MC_phistar_temp=ConvertToTGraph(mc_truegen_fsr_pileup[u]);
    g_MC_phistar_fsr_pileup.push_back(g_MC_phistar_temp);
  }

  NormalizeGraph(g_MC_phistar);
  NormalizeGraph(g_MC_phistar_cteq);
  NormalizeGraph(g_MC_phistar_fsr_pileup);

  TGraphAsymmErrors* g_syst_mc_cteq = CalcTotalSysU(g_MC_phistar_cteq,g_MC_phistar_cteq[0],1,1);
  TGraphAsymmErrors* g_syst_mc_fsr  = CalcTotalSysU_fsr(g_MC_phistar_fsr_pileup[0],g_MC_phistar[0]);
  
  vector<TGraphAsymmErrors *> g_mc_syst_muon;
  g_mc_syst_muon.push_back(g_MC_phistar[0]);
  vector<std::string> syst_mc_list_muon;
  syst_mc_list_muon.push_back("unfolding");
  TGraphAsymmErrors* g_mc_final_muon=GetMCFinal(g_mc_syst_muon, syst_mc_list_muon,0);
   
  vector<TGraphAsymmErrors *> g_mc_syst;
  g_mc_syst.push_back(g_MC_phistar[0]);
  g_mc_syst.push_back(g_syst_mc_fsr);
  if (!doMG) g_mc_syst.push_back(g_syst_mc_cteq);
  vector<std::string> syst_mc_list;
  syst_mc_list.push_back("unfolding");
  syst_mc_list.push_back("fsr");
  if (!doMG) syst_mc_list.push_back("cteq");
   
  TGraphAsymmErrors* g_mc_final=GetMCFinal(g_mc_syst, syst_mc_list);
  
  PreparePlots(*g_mc_final,0);
  
  TH1D* h_mc= ConvertToHist(g_mc_final_muon,0);

  std::string textn="Final_Hist_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  TFile tr(textn.c_str(),"RECREATE");
  h_mc->Write();

  textn="Data_Graph_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  TFile tr2(textn.c_str(),"RECREATE");
  g_mc_final->Write();

}
