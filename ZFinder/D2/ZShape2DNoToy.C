#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TRandom.h"
#include "TMatrix.h"
#include "TMatrixD.h"
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
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
//#include <locale>
using namespace std;

const bool debug = false; // makes print statments and ends many of the for loops earlier so the code runs faster if true

// const int Ntoys = 500;
// const int Ntoys2 = 100;
const int Ntoys = 500 ;
const int Ntoys2 = 100 ;
const int elec = 1;
const int doMG = 0;

const std::string Tag = "";
const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 5, 10};
const double yBins[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
// const std::string Tag = "Debug_";
//const double phistarBins[] = {0.000, 0.016, 0.039, 0.072, 0.145, 0.391, 1.496, 5., 10.};
//const double yBins[] = {0.0,0.8,1.6,2.4};
// const std::string Tag = "1D_";
// const double yBins[] = {0.0,10000.0};
// const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
size_t ny=(sizeof(yBins)/sizeof(yBins[0]))-1;
size_t nbins=nphistar*ny;

const double TEfficiencyEtaBins[11] = {-2.1, -2.0, -1.556, -1.442, -0.8, 0.0, 0.8, 1.442, 1.556, 2.0, 2.1};
const double TEfficiencyETBins[5]   = {30, 40, 50, 70, 5000};
const double TEfficiencyData[10][4][2] = {{{0.741, 0.003},  {0.773, 0.003},  {0.780, 0.005},  {0.790, 0.010}},
					  {{0.734, 0.001},  {0.772, 0.001},  {0.786, 0.002},  {0.792, 0.005}},
					  {{0.725, 0.003},  {0.821, 0.002},  {0.809, 0.004},  {0.848, 0.010}},
					  {{0.893, 0.0005}, {0.9396, 0.0003},{0.9509, 0.0006},{0.966, 0.001}},
					  {{0.9213, 0.0004},{0.9528, 0.0002},{0.9601, 0.0004},{0.969, 0.001}},
					  {{0.9174, 0.0004},{0.9473, 0.0004},{0.9561, 0.0004},{0.963, 0.001}},
					  {{0.8964, 0.0005},{0.9424, 0.0003},{0.9533, 0.0006},{0.966, 0.001}},
					  {{0.714, 0.003},  {0.823, 0.002},  {0.827, 0.004},  {0.861, 0.010}},
					  {{0.758, 0.001},  {0.800, 0.001},  {0.811, 0.002},  {0.823, 0.005}},
					  {{0.764, 0.003},  {0.792, 0.002},  {0.797, 0.005},  {0.820, 0.010}}};
const double TEfficiencyMC[10][4][2]   = {{{0.734, 0.004},  {0.769, 0.003},  {0.771, 0.004},  {0.760, 0.020}},
				  	  {{0.736, 0.002},  {0.768, 0.004},  {0.779, 0.003},  {0.789, 0.008}},
					  {{0.791, 0.004},  {0.847, 0.003},  {0.850, 0.006},  {0.870, 0.020}},
					  {{0.9395, 0.0006},{0.9612, 0.0004},{0.9690, 0.0007},{0.980, 0.002}},
					  {{0.9469, 0.0005},{0.9670, 0.0003},{0.9745, 0.0005},{0.982, 0.001}},
					  {{0.9466, 0.0005},{0.9665, 0.0003},{0.9739, 0.0005},{0.982, 0.001}},
					  {{0.9364, 0.0007},{0.9597, 0.0004},{0.9668, 0.0008},{0.979, 0.002}},
					  {{0.779, 0.004},  {0.841, 0.003},  {0.842, 0.006},  {0.860, 0.020}},
					  {{0.749, 0.002},  {0.786, 0.002},  {0.798, 0.003},  {0.810, 0.008}},
					  {{0.737, 0.004},  {0.769, 0.004},  {0.779, 0.007},  {0.820, 0.020}}};
const double EfficiencyEtaBins[6] = {0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};
const double EfficiencyETBins[5]  = {18, 30, 40, 50, 1000000};
const double EfficiencySF[5][4][3]       = {{{0.982, 0.003, 0.012},{0.988, 0.001, 0.008},{0.990, 0.001, 0.004},{0.990, 0.001, 0.004}},
					    {{0.993, 0.002, 0.012},{0.993, 0.001, 0.008},{0.993, 0.001, 0.004},{0.991, 0.001, 0.004}},
					    {{1.016, 0.012, 0.020},{0.985, 0.004, 0.009},{0.987, 0.004, 0.004},{0.974, 0.009, 0.006}},
					    {{0.988, 0.003, 0.012},{0.993, 0.002, 0.008},{0.992, 0.001, 0.004},{0.990, 0.003, 0.004}},
					    {{1.002, 0.004, 0.012},{1.004, 0.002, 0.008},{1.005, 0.002, 0.004},{0.998, 0.004, 0.004}}};
const double EfficiencyMediumSF[5][4][3] = {{{0.986, 0.002, 0.001},{1.002, 0.001, 0.001},{1.005, 0.001, 0.001},{1.004, 0.001, 0.001}},
					    {{0.959, 0.003, 0.003},{0.980, 0.001, 0.001},{0.988, 0.001, 0.001},{0.988, 0.002, 0.002}},
					    {{0.967, 0.007, 0.013},{0.950, 0.006, 0.007},{0.958, 0.005, 0.005},{0.966, 0.009, 0.009}},
					    {{0.941, 0.005, 0.005},{0.967, 0.003, 0.003},{0.992, 0.002, 0.002},{1.000, 0.003, 0.003}},
					    {{1.020, 0.003, 0.003},{1.021, 0.003, 0.003},{1.019, 0.002, 0.002},{1.022, 0.004, 0.004}}};
const double EfficiencyTightSF[5][4][3]  = {{{0.960, 0.003, 0.003},{0.978, 0.001, 0.001},{0.981, 0.001, 0.001},{0.982, 0.002, 0.002}},
					    {{0.936, 0.004, 0.004},{0.958, 0.002, 0.002},{0.969, 0.001, 0.001},{0.969, 0.002, 0.002}},
					    {{0.933, 0.015, 0.017},{0.907, 0.008, 0.008},{0.904, 0.004, 0.004},{0.929, 0.011, 0.011}},
					    {{0.879, 0.007, 0.007},{0.909, 0.003, 0.003},{0.942, 0.002, 0.002},{0.957, 0.004, 0.004}},
					    {{0.974, 0.004, 0.004},{0.987, 0.004, 0.004},{0.991, 0.003, 0.003},{0.999, 0.005, 0.005}}};


//this makes it so if we don't want to look at YSeperation
std::string File_Data = "/afs/cern.ch/work/r/ruckstuh/public/Data_R9.root";
std::string File_Data_Pt_L = "/afs/cern.ch/work/r/ruckstuh/public/Data_Low_R9.root";
std::string File_Data_Pt_H = "/afs/cern.ch/work/r/ruckstuh/public/Data_High_R9.root";
std::string File_Signal_reco = "/afs/cern.ch/work/r/ruckstuh/public/MG_dressed_reco.root";
std::string File_Signal_reco_born = "/afs/cern.ch/work/r/ruckstuh/public/MG_born_reco.root";
std::string File_Signal_reco_bare = "/afs/cern.ch/work/r/ruckstuh/public/MG_bare_reco.root";
std::string File_Signal_gen = "/afs/cern.ch/work/r/ruckstuh/public/MG_dressed_gen.root";
std::string File_Signal_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/MG_born_gen.root";
std::string File_Signal_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/MG_bare_gen.root";
std::string File_Powheg_reco = "/afs/cern.ch/work/r/ruckstuh/public/PH_dressed_reco.root";
std::string File_Powheg_reco_born = "/afs/cern.ch/work/r/ruckstuh/public/PH_born_reco.root";
std::string File_Powheg_reco_bare = "/afs/cern.ch/work/r/ruckstuh/public/PH_bare_reco.root";
std::string File_Powheg_gen = "/afs/cern.ch/work/r/ruckstuh/public/PH_dressed_gen.root";
std::string File_Powheg_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/PH_born_gen.root";
std::string File_Powheg_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/PH_bare_gen.root";
std::string File_tt = "/afs/cern.ch/work/r/ruckstuh/public/TTbar.root";
std::string File_tautau = "/afs/cern.ch/work/r/ruckstuh/public/TauTau.root";
std::string File_tbarw = "/afs/cern.ch/work/r/ruckstuh/public/TbarW.root";
std::string File_tw = "/afs/cern.ch/work/r/ruckstuh/public/TW.root";
std::string File_ww = "/afs/cern.ch/work/r/ruckstuh/public/WW.root";
std::string File_wz = "/afs/cern.ch/work/r/ruckstuh/public/WZ.root";
std::string File_zz = "/afs/cern.ch/work/r/ruckstuh/public/ZZ.root";

std::string reco_name = "Combined Single Reco";
// std::string reco_name2 = "Combined Single 3020 Reco";
std::string reco_name_en_l = "Combined Single Lowered Threshold Reco";
std::string reco_name_en_h = "Combined Single Higher Threshold Reco";
//std::string reco_name_en="Combined Single Reco";
std::string gen_name = "Combined Gen Cuts Reco";

void NoNeg(TH1D* hist){
  for (uint i=0; i<nbins; i++){
    if (hist->GetBinContent(i+1)<0){
      double val =hist->GetBinContent(i+1);
      double error=hist->GetBinError(i+1);
      error = val+error;
      if (error<0) error=0;
      hist->SetBinContent(i+1,0);
      hist->SetBinError(i+1,error);
      cout<<"BIN was nagative. bin:"<<i<<endl;
    }
  }  
}

void NoNeg(vector<TH1D*> & hist){
  for (uint i=0; i<hist.size(); i++){
    NoNeg(hist[i]);
  }
}

void DeleteVec(vector<TH1D*> hist){
  for (uint i=0; i<hist.size();i++){
    delete hist[i];
  }
}

void DeleteVec(vector<RooUnfoldResponse*> hist){
  for (uint i=0; i<hist.size();i++){
    delete hist[i];
  }
}

void DeleteVec(vector<TGraphAsymmErrors*> hist){
  for (uint i=0; i<hist.size();i++){
    delete hist[i];
  }
}

int GetBin(double phistar, double y){
  TAxis* A_phistar=new TAxis(nphistar,phistarBins);
  TAxis* A_y=new TAxis(ny,yBins);
  int bin_phistar=A_phistar->FindBin(phistar)-1;
  int bin_y=A_y->FindBin(fabs(y))-1;
  int bin=bin_y*nphistar+bin_phistar;
  if (bin_phistar<0 || bin_y<0) bin=-1;
  if (bin_phistar>=int(nphistar) || bin_y>=int(ny)) bin=nbins;
  delete A_phistar;
  delete A_y;
  return bin;
}

TGraphAsymmErrors * GetUnfoldingSysUn(bool norm, TGraphAsymmErrors* NominalGraph){
  TGraphAsymmErrors* g = (TGraphAsymmErrors *) NominalGraph->Clone();
  TFile f("SysUnFunctionValues.root");
  TH1D *syserr=0;
  if(ny==1 && nphistar==35){
    if (norm) syserr= (TH1D*) f.Get("NormErrorValues");
    else syserr= (TH1D*) f.Get("AbsoluteErrorValues");
  }
  // else if(ny==6 && nphistar ==36){
  // }
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y;
    g->GetPoint(ibin, x, y);
    double er=0; ///TO BE CHANGED
    if (syserr){
      er=fabs(syserr->GetBinContent(ibin+1));
      // cout<<"syserr exist"<<endl;
    }
    g->SetPointError(ibin, 0, 0,er*y, er*y);  
  }
  //  delete syserr;
  return g;
}

TGraphAsymmErrors * ErrorSmoother(TGraphAsymmErrors* OrginalGraph) {
  TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins); //used for smoothing the errors out
  double x, y, y2;
  double lowerrFracD, centerrFracD, higherrFracD, lowerrFracH, centerrFracH, higherrFracH; //should be fraction of whole unless fucked up
  double newErrorhigh;
  double newErrorlow;
   
  OrginalGraph->GetPoint(0, x, y);
  g->SetPoint(0, x, y);
  centerrFracD = OrginalGraph->GetErrorYlow(0) / y;
  centerrFracH = OrginalGraph->GetErrorYhigh(0) / y;
  OrginalGraph->GetPoint(1, x, y2);
  higherrFracD = OrginalGraph->GetErrorYlow(1) / y2;
  higherrFracH = OrginalGraph->GetErrorYhigh(1) / y2;
  newErrorlow = (2 * centerrFracD + higherrFracD) / 3 * y;
  newErrorhigh = (2 * centerrFracH + higherrFracH) / 3 * y;
  g->SetPointError(0, 0, 0, newErrorlow, newErrorhigh);

  for (size_t ibin = 1; ibin < nbins - 2; ibin++) {
    OrginalGraph->GetPoint(ibin, x, y);
    g->SetPoint(ibin, x, y);
    centerrFracD = OrginalGraph->GetErrorYlow(ibin) / y;
    centerrFracH = OrginalGraph->GetErrorYhigh(ibin) / y;
    OrginalGraph->GetPoint(ibin - 1, x, y2);
    lowerrFracD = OrginalGraph->GetErrorYlow(ibin - 1) / y2;
    lowerrFracH = OrginalGraph->GetErrorYhigh(ibin - 1) / y2;
    OrginalGraph->GetPoint(ibin + 1, x, y2);
    higherrFracD = OrginalGraph->GetErrorYlow(ibin + 1) / y2;
    higherrFracH = OrginalGraph->GetErrorYhigh(ibin + 1) / y2;
    newErrorlow = (higherrFracD + lowerrFracD + centerrFracD) / 3;
    newErrorhigh = (higherrFracH + lowerrFracH + centerrFracH) / 3;
    g->SetPointError(ibin, 0, 0, newErrorlow*y, newErrorhigh * y);
  }
  OrginalGraph->GetPoint(nbins - 2, x, y);
  g->SetPoint(nbins - 2, x, y);
  centerrFracD = OrginalGraph->GetErrorYlow(nbins - 2) / y;
  centerrFracH = OrginalGraph->GetErrorYhigh(nbins - 2) / y;
  OrginalGraph->GetPoint(nbins - 3, x, y2);
  lowerrFracD = OrginalGraph->GetErrorYlow(nbins - 3) / y2;
  lowerrFracH = OrginalGraph->GetErrorYhigh(nbins - 3) / y2;
  newErrorlow = (2 * centerrFracD + lowerrFracD) / 3 * y;
  newErrorhigh = (2 * centerrFracH + lowerrFracH) / 3 * y;
  g->SetPointError(nbins - 2, 0, 0, newErrorlow, newErrorhigh);
  
  OrginalGraph->GetPoint(nbins - 1, x, y);
  g->SetPoint(nbins - 1, x, y);
  g->SetPointError(nbins - 1, 0, 0, OrginalGraph->GetErrorYlow(nbins - 1), OrginalGraph->GetErrorYhigh(nbins - 1));
  return g;
}

TMatrixD CalcCovM_lumi(TGraphAsymmErrors* graph){
  TMatrixD cov(nbins,nbins);
  double error=0.026;
  for(uint j = 0; j < nbins; j++) {
    double xj, yj, xk, yk;
    graph->GetPoint(j, xj, yj);
    for(uint k = 0; k < nbins; k++) {
      graph->GetPoint(k, xk, yk);
      double Covariance=(error*yj)*(error*yk);
      cov(j,k)=Covariance;
    }
  }
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"LUMI:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  return cov;
}

TMatrixD CalcCovM_unfsys(TGraphAsymmErrors* graph, bool norm){
  TFile f("SysUnFunctionValues.root");
  TH1D *syserr=0;
  if(ny==1 && nphistar==35){
    if (norm) syserr= (TH1D*) f.Get("NormErrorValues");
    else syserr= (TH1D*) f.Get("AbsoluteErrorValues");
  }
  // if(ny==2 && nphistar==36){
  // }
  TMatrixD cov(nbins,nbins);
  for(uint j = 0; j < nbins; j++) {
    double xj, yj, xk, yk;
    double errorj, errork;
    graph->GetPoint(j, xj, yj);
    errorj=0; ///TO BE CHANGED
    if (syserr) {
      errorj=syserr->GetBinContent(j+1);
      // cout<<"syserr exist"<<endl;
    }
    for(uint k = 0; k < nbins; k++) {
      graph->GetPoint(k, xk, yk);
      errork=0; ///TO BE CHANGED
      if (syserr) errork=syserr->GetBinContent(k+1);
      double Covariance=(errorj*yj)*(errork*yk);
      cov(j,k)=Covariance;
    }
  }
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"UNFSYS:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  //  delete syserr;
  return cov;
}

TMatrixD CalcCovM_down(TGraphAsymmErrors* graphd, TGraphAsymmErrors* graphn){
  TMatrixD cov(nbins,nbins);
  for(uint j = 0; j < nbins; j++) {
    double xjd, yjd, xkd, ykd;
    double xjn, yjn, xkn, ykn;
    graphd->GetPoint(j, xjd, yjd);
    graphn->GetPoint(j, xjn, yjn);
    for(uint k = 0; k < nbins; k++) {
      graphd->GetPoint(k, xkd, ykd);
      graphn->GetPoint(k, xkn, ykn);
      double Covariance=(yjd-yjn)*(ykd-ykn);
      cov(j,k)=Covariance;
    }
  }
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"DOWN:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  return cov;
}
TMatrixD CalcCovM_downup(TGraphAsymmErrors* graphd, TGraphAsymmErrors* graphu, TGraphAsymmErrors* graphn){
  TMatrixD cov(nbins,nbins);
  for(uint j = 0; j < nbins; j++) {
    double xjd, yjd, xkd, ykd;
    double xju, yju, xku, yku;
    double xjn, yjn, xkn, ykn;
    graphd->GetPoint(j, xjd, yjd);
    graphu->GetPoint(j, xju, yju);
    graphn->GetPoint(j, xjn, yjn);
    for(uint k = 0; k < nbins; k++) {
      graphd->GetPoint(k, xkd, ykd);
      graphu->GetPoint(k, xku, yku);
      graphn->GetPoint(k, xkn, ykn);
      double Covariance=((yjd-yjn)*(ykd-ykn)+(yju-yjn)*(yku-ykn))/2.;
      cov(j,k)=Covariance;
    }
  }
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"DOWNUP:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  return cov;
}
TMatrixD CalcCovM_cteq(vector<TGraphAsymmErrors*> graph){
  TMatrixD cov(nbins,nbins);
  int n=(graph.size()-1)/2;
  for(int i=0; i<n; i++){
    TMatrixD covi=CalcCovM_downup(graph[i*2+1], graph[i*2+2], graph[0]);
    if (i==0) cov=covi;
    else cov+=covi;
  }
  cov*=(1/(1.645*1.645));
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"CTEQ:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  return cov;
}

TMatrixD CalcCovM_uncor(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2){
  TMatrixD cov(nbins,nbins);
  for(uint j = 0; j < nbins; j++) {
    for(uint k = 0; k < nbins; k++) {
      if (j!=k) cov(j,k)=0;
      else{
	double x1, y1, errorl1, errorh1;
	double x2, y2, errorl2, errorh2;
	graph1->GetPoint(j, x1, y1);
	graph2->GetPoint(j, x2, y2);
	if (x1 != x2 || y1 != y2) cout << "eff stat: nominal values don't agree: " << x1 << " " << x2 << " : " << y1 << " " << y2 << endl;
	errorl1 = graph1->GetErrorYlow(j);
	errorh1 = graph1->GetErrorYhigh(j);
	errorl2 = graph2->GetErrorYlow(j);
	errorh2 = graph2->GetErrorYhigh(j);
	if (errorl2 < errorl1 || errorh2 < errorh1) cout << "eff stat: uncertainty size don't agree: " << errorl1 << " " << errorl2 << " : " << errorh1 << " " << errorh2 << endl;
	double errorl=(errorl2 * errorl2)-(errorl1 * errorl1);
	double errorh=(errorh2 * errorh2)-(errorh1 * errorh1);
	double Covariance=TMath::Max(errorl,errorh);
	cov(j,k)=Covariance;
      }
    }
  }
  return cov;
}

TMatrixD CalcCovM_toy(vector<TGraphAsymmErrors*> graphs, int part=1, bool eff=0){
  double Mean_PROD_SUM[nbins][nbins];
  double Mean_SUM[nbins]; //holds the sum of all the toys for the bin
  double xj, yj, xk, yk;
  int nt = Ntoys;
  if (eff) nt = Ntoys2;
  uint  gs = 1 + (nt * (part - 1));
  for(uint i = gs; i < gs+nt; i++) {
    for(uint j = 0; j < nbins; j++) {
      if (i == gs) {
        Mean_SUM[j] = 0;
      }
      graphs[i]->GetPoint(j, xj, yj);
      Mean_SUM[j] += yj;
      for(uint k = 0; k < nbins; k++) {
        if (i == gs) {
          Mean_PROD_SUM[j][k] = 0;
        }
	graphs[i]->GetPoint(k, xk, yk);
	Mean_PROD_SUM[j][k] += yj*yk;
      }
    }
  }
  TMatrixD cov(nbins,nbins);
  for(uint j = 0; j < nbins; j++) {
    for(uint k = 0; k < nbins; k++) {
      double Covariance = ((Mean_PROD_SUM[j][k]) - (Mean_SUM[j] * Mean_SUM[k]/double(nt)))/double(nt-1);
      cov(j,k)=Covariance;
    }
  }
  //test symmetry.
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (cov(j,k)!=cov(k,j)) cout<<"TOY:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<cov(j,k)<<" != "<<cov(k,j)<<endl;
    }
  }
  return cov;
}

TH1D * ConvertToHist(TGraphAsymmErrors* g, std::string name) {
  TH1D* h_temp;
  h_temp = new TH1D(name.c_str(), name.c_str(), nbins, 0, nbins);
  h_temp->Sumw2();
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y;
    g->GetPoint(ibin, x, y);
    double error = g->GetErrorYhigh(ibin);
    h_temp->SetBinContent(ibin + 1, y);
    h_temp->SetBinError(ibin + 1, error);
  }
  return h_temp;
}

vector<TH2D*> GetEffTMCToys(bool MC = 1) {
  gErrorIgnoreLevel = kError;
  vector<TH2D*> EffToys;
  for (int t = 0; t < Ntoys2 + 1; t++) {
    TH2D *h_eff = new TH2D("h_eff", "h_eff", 10, TEfficiencyEtaBins, 4, TEfficiencyETBins);
    EffToys.push_back(h_eff);
  }
  for (int e = 0; e < 10; e++) {
    for (int p = 0; p < 4; p++) {
      double mean = TEfficiencyMC[e][p][0];
      double error = TEfficiencyMC[e][p][1];
      if (!MC) {
        mean = TEfficiencyData[e][p][0];
        error = TEfficiencyData[e][p][1];
      }
      EffToys[0]->SetBinContent(e + 1, p + 1, mean);
      for (int t = 0; t < Ntoys2; t++) {
        // gRandom->SetSeed(2537+t*200+e*10+p);
        double x = gRandom->Gaus(mean, error);
        if (x > 1) cout << "Error efficiency is larger then 1" << endl;
        EffToys[t + 1]->SetBinContent(e + 1, p + 1, x);
      }
    }
  }
  return EffToys;
}

vector<TH2D*> GetEffSFToys(int type = 0) {
  gErrorIgnoreLevel = kError;
  vector<TH2D*> EffToys;
  for (int t = 0; t < Ntoys2 + 1; t++) {
    TH2D *h_eff = new TH2D("h_eff", "h_eff", 5, EfficiencyEtaBins, 4, EfficiencyETBins);
    EffToys.push_back(h_eff);
  }
  for (int e = 0; e < 5; e++) {
    for (int p = 0; p < 4; p++) {
      double mean = EfficiencySF[e][p][0];
      double error = sqrt(EfficiencySF[e][p][1] * EfficiencySF[e][p][1] + EfficiencySF[e][p][2] * EfficiencySF[e][p][2]);
      double errorl = sqrt(EfficiencySF[e][p][1] * EfficiencySF[e][p][1] + EfficiencySF[e][p][2] * EfficiencySF[e][p][2]);
      if (type == 1) {
        mean = EfficiencyMediumSF[e][p][0];
        error = EfficiencyMediumSF[e][p][1];
        errorl = EfficiencyMediumSF[e][p][2];
      }
      if (type == 2) {
        mean = EfficiencyTightSF[e][p][0];
        error = EfficiencyTightSF[e][p][1];
        errorl = EfficiencyTightSF[e][p][2];
      }
      EffToys[0]->SetBinContent(e + 1, p + 1, mean);
      for (int t = 0; t < Ntoys2; t++) {
        // gRandom->SetSeed(27349302+t*200+e*10+p);
        double x = gRandom->Gaus(mean, error);
        if (x < mean) {
          x = mean - ((mean - x) * errorl / error);
        }
        EffToys[t + 1]->SetBinContent(e + 1, p + 1, x);
      }
    }
  }
  //  cout<<"doei"<<endl;
 return EffToys;
}

vector<TH1D*> GetToyBg() {
  TH1F *bg_sf_full=0;
  TFile f_bg("ratio_data_mc_emu.root");
  if (ny==1 && nphistar==35){
    TCanvas *c_bg = (TCanvas*) f_bg.Get("Canvas_1");
    bg_sf_full = (TH1F*) c_bg->FindObject("hPull");
  }
  if (ny==6 && nphistar==36){
    bg_sf_full= new TH1F("bg_sf","bg_sf",nbins,0,nbins);
    for (uint i=0; i<ny-1; i++){
      std::string filename="ratio_data_mc_bin45.root";
      if (i==0) filename="ratio_data_mc_bin0.root";
      if (i==1) filename="ratio_data_mc_bin1.root";
      if (i==2) filename="ratio_data_mc_bin2.root";
      if (i==3) filename="ratio_data_mc_bin3.root";
      TFile f_bg_bin(filename.c_str());
      TH1F* temp=(TH1F*) f_bg_bin.Get("hPull");
      for (uint j=0; j<nphistar; j++){
    	int bin=i*nphistar+j;
	if (filename=="ratio_data_mc_bin45.root" && j==nphistar-1){
	  bg_sf_full->SetBinContent(bin+1, 1);
	  bg_sf_full->SetBinError(bin+1, 1);	
	}
	else{
	  bg_sf_full->SetBinContent(bin+1, temp->GetBinContent(j+1));
	  bg_sf_full->SetBinError(bin+1, temp->GetBinContent(j+1));	
	}
      }
    }
  }
  vector<TH1D*> bg_sf;
  for (int t = 0; t < Ntoys + 1; t++) {
    TH1D* bg_sf_temp = new TH1D("bg", "bg", nbins, 0, nbins);
    bg_sf_temp->Sumw2();
    for (uint i = 0; i < nbins; i++) {//NEED TO BE CHANGED
      double val =1;
      if (bg_sf_full) {
	val= bg_sf_full->GetBinContent(i + 1);
	// cout<<"bg_sf_full exists"<<endl;
      }
      if (t == 0) bg_sf_temp->SetBinContent(i + 1, val);//NEED TO BE CHANGED
      else {
	double error = 0;
	if (bg_sf_full) error = bg_sf_full->GetBinError(i + 1);
	double x = gRandom->Gaus(val, error);
	if (x<0) x=0;
	bg_sf_temp->SetBinContent(i + 1, x);
      }
      bg_sf_temp->SetBinError(i + 1, 0.);
    }
    bg_sf.push_back(bg_sf_temp);
  }
  f_bg.Close();
  return bg_sf;
}

TH1D* GetBgSS() {
  TH1D *bg = new TH1D("bg", "bg", nbins, 0, nbins);
  bg->Sumw2();
  if(ny==1 && nphistar==35){
    TFile f_ss("SS.root");
    TH1D* ss_full = (TH1D*) f_ss.Get("qcd_phistar");
    for (uint i = 0; i < nbins; i++) {
      double val= ss_full->GetBinContent(i + 1);
      if (i == nbins - 1) val=0;
      bg->SetBinContent(i + 1, val);
      bg->SetBinError(i + 1, 0.);
    }
    f_ss.Close();
  }
  else if (ny==6 && nphistar==36){
    for (uint j=0; j<ny; j++){
      std::ostringstream strs;
      strs << j;
      std::string filename="SSYBin"+strs.str()+".root";
      TFile f_ss(filename.c_str());
      TH1D* ss=(TH1D*) f_ss.Get("qcd/phistar");
      for (uint i = 0; i < nphistar; i++) {
	int bin=j*nphistar+i;
	bg->SetBinContent(bin+1,ss->GetBinContent(i+1)/19712.);
	bg->SetBinError(bin+1, 0.);
      }
      f_ss.Close();
    }
  }
  else {
    for (uint i = 0; i < nbins; i++) {//NEED TO BE CHANGED
      bg->SetBinContent(i + 1, 0.);
      bg->SetBinError(i + 1, 0.);
    }
  }
  return bg;
}

void NormalizeGraph(vector<TGraphAsymmErrors*> &graph, bool doNorm = 0) {
  for (size_t i = 0; i < graph.size(); i++) {
    TGraphAsymmErrors *graphtemp = (TGraphAsymmErrors *) graph[i]->Clone();
    double xstot = 0;
    double x, y, errorl, errorh;
    for (size_t ibin = 0; ibin < nbins; ibin++) {
      graphtemp->GetPoint(ibin, x, y);
      xstot += y;
    }
    // cout<<i<<":Normalisation:"<<xstot<<endl;
    for (size_t ibin = 0; ibin < nbins; ibin++) {
      double dphistar = phistarBins[(ibin % nphistar) + 1] - phistarBins[ibin % nphistar]; 
      double dy = yBins[ibin/nphistar + 1] - yBins[ibin/nphistar];
      if (ny==1) dy=1.;
      // cout<<ibin<<" "<<ibin/nphistar<<" "<<ibin % nphistar<<" "<<dphistar<<" "<<dy<<endl;
      double norm = dphistar*dy*xstot;
      double errorbin2_h = 0;
      double errorbin2_l = 0;
      graphtemp->GetPoint(ibin, x, y);
      if (doNorm) {
        for (uint j = 0; j < nbins; j++) {
          double errorbin_j_on_i_h = 0;
          double errorbin_j_on_i_l = 0;
          if (ibin == j) {
            errorbin_j_on_i_h = graphtemp->GetErrorYhigh(j) * ((1 / xstot)-(y / (xstot * xstot)));
            errorbin_j_on_i_l = graphtemp->GetErrorYlow(j) * ((1 / xstot)-(y / (xstot * xstot)));
          } else {
            errorbin_j_on_i_h = graphtemp->GetErrorYhigh(j)*(y / (xstot * xstot));
            errorbin_j_on_i_l = graphtemp->GetErrorYlow(j)*(y / (xstot * xstot));
          }
          errorbin2_h += errorbin_j_on_i_h*errorbin_j_on_i_h;
          errorbin2_l += errorbin_j_on_i_l*errorbin_j_on_i_l;
        }
      }
      if (!doNorm) norm = dphistar*dy;
      errorl = graphtemp->GetErrorYlow(ibin);
      errorh = graphtemp->GetErrorYhigh(ibin);
      graph[i]->SetPoint(ibin, x, y / norm);
      graph[i]->SetPointError(ibin, 0, 0, errorl / norm, errorh / norm);
      if (doNorm) graph[i]->SetPointError(ibin, 0, 0, sqrt(errorbin2_l) / (dphistar*dy), sqrt(errorbin2_h) / (dphistar*dy));
    }
    delete graphtemp;
  }
}
void NormalizeGraph(vector<TGraphAsymmErrors*> &graph, bool doNorm, TMatrixD &CovM) {
  for (size_t i = 0; i < graph.size(); i++) {
    TGraphAsymmErrors * graphtemp = (TGraphAsymmErrors *) graph[i]->Clone();
    double xstot = 0;
    double x, y, errorl, errorh;
    for (size_t ibin = 0; ibin < nbins; ibin++) {
      graphtemp->GetPoint(ibin, x, y);
      xstot += y;
    }
    // cout<<i<<":Normalisation:"<<xstot<<endl;
    for (uint ibin = 0; ibin < nbins; ibin++) {
      double dphistar = phistarBins[(ibin % nphistar) + 1] - phistarBins[ibin % nphistar];
      double dy = yBins[ibin/nphistar + 1] - yBins[ibin/nphistar];
      if (ny==1) dy=1.;
      double norm = dphistar*dy*xstot;
      double errorbin2_h = 0;
      double errorbin2_l = 0;
      graphtemp->GetPoint(ibin, x, y);
      if (doNorm) {
        for (uint j = 0; j < nbins; j++) {
          double errorbin_j_on_i_h = 0;
          double errorbin_j_on_i_l = 0;
          if (ibin == j) {
            errorbin_j_on_i_h = graphtemp->GetErrorYhigh(j) * ((1 / xstot)-(y / (xstot * xstot)));
            errorbin_j_on_i_l = graphtemp->GetErrorYlow(j) * ((1 / xstot)-(y / (xstot * xstot)));
          } else {
            errorbin_j_on_i_h = graphtemp->GetErrorYhigh(j)*(y / (xstot * xstot));
            errorbin_j_on_i_l = graphtemp->GetErrorYlow(j)*(y / (xstot * xstot));
          }
          errorbin2_h += errorbin_j_on_i_h*errorbin_j_on_i_h;
          errorbin2_l += errorbin_j_on_i_l*errorbin_j_on_i_l;
        }
      }
      if (!doNorm) norm = dphistar*dy;
      errorl = graphtemp->GetErrorYlow(ibin);
      errorh = graphtemp->GetErrorYhigh(ibin);
      graph[i]->SetPoint(ibin, x, y / norm);
      graph[i]->SetPointError(ibin, 0, 0, errorl / norm, errorh / norm);
      if (doNorm) graph[i]->SetPointError(ibin, 0, 0, sqrt(errorbin2_l) / (dphistar*dy), sqrt(errorbin2_h) / (dphistar*dy));
      if (i==0 && ibin!=nbins){
	for (uint j = 0; j < nbins; j++) {
 	  double dphistarj = phistarBins[(j % nphistar) + 1] - phistarBins[j % nphistar];
	  double dyj = yBins[j/nphistar + 1] - yBins[j/nphistar];
	  if (ny==1) dyj=1.;
 	  double normj = dphistarj*dyj*xstot;
	  if (!doNorm) normj = dphistarj*dyj;
	  CovM(ibin,j)=CovM(ibin,j)/(normj*norm);
	}
      }
    }
    delete graphtemp;
  }
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (fabs(CovM(j,k)/CovM(k,j)-1)>0.00001) cout<<"norm:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<CovM(j,k)<<" != "<<CovM(k,j)<<endl;
    }
  }
}

TGraphAsymmErrors * CalcTotalSysU_updown(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal, bool same = 0, bool cteq = 0) { // adds everything in quadrature (up down)
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, errorl, errorh;
    graph_nominal->GetPoint(ibin, x, y);
    g_un->SetPoint(ibin, x, y);
    errorl = 0;
    errorh = 0;
    // int last = 0;
    //    int now = 0;
    // double diflast = 0;
    for (size_t i = same; i < graph.size(); i++) {
      //now = 0;
      double xtemp, ytemp;
      graph[i]->GetPoint(ibin, xtemp, ytemp);
      if (xtemp != x) cout << "This is really weird and wrong" << endl;
      double dif = ytemp - y;
      if (dif / y > 0.1)cout << "WOW very large error: " << dif / y << " " << y << " " << ytemp << endl;
      if (dif > 0) {
        errorh = sqrt(errorh * errorh + dif * dif);
	// now = 1.0;
      }
      if (dif < 0) {
        errorl = sqrt(errorl * errorl + dif * dif);
	//   now = -1.0;
      }
      if (i > 0 && i % 2 != same) {
        // if (last != 0 && now != 0 && last == now && fabs(diflast / y) > 0.0001 && fabs(dif / y) > 0.0001) cout << "both change go in the same dirrection: " << i << " " << ibin << " " << dif / y << " " << diflast / y << endl;
      }
      // last = now;
      // diflast = dif;
    }
    if (cteq) {
      errorl = errorl / 1.645;
      errorh = errorh / 1.645;
    }
    g_un->SetPointError(ibin, 0, 0, errorl, errorh);
  }
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_comb3(TGraphAsymmErrors* graph, TGraphAsymmErrors* graph_nominal, TGraphAsymmErrors* graph_mctoy, bool add = 1) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, errorl, errorh;
    double x2, y2, errorl2, errorh2;
    double errorl3, errorh3;
    graph_nominal->GetPoint(ibin, x, y);
    graph->GetPoint(ibin, x2, y2);
    if (x != x2 || y != y2) cout << "eff stat: nominal values don't agree: " << x << " " << x2 << " : " << y << " " << y2 << endl;
    g_un->SetPoint(ibin, x, y);
    errorl = graph_nominal->GetErrorYlow(ibin);
    errorh = graph_nominal->GetErrorYhigh(ibin);
    errorl2 = graph->GetErrorYlow(ibin);
    errorh2 = graph->GetErrorYhigh(ibin);
    errorl3 = graph_mctoy->GetErrorYlow(ibin);
    errorh3 = graph_mctoy->GetErrorYhigh(ibin);
    if (!add) {
      if (errorl2 < errorl || errorh2 < errorh) cout << "eff stat: uncertainty size don't agree: " << errorl << " " << errorl2 << " : " << errorh << " " << errorh2 << endl;
      g_un->SetPointError(ibin, 0, 0, sqrt((errorl2 * errorl2)-(errorl * errorl)+(errorl3 * errorl3)), sqrt((errorh2 * errorh2)-(errorh * errorh)+(errorh3 * errorh3)));
    } else {
      g_un->SetPointError(ibin, 0, 0, sqrt((errorl2 * errorl2)+(errorl * errorl)+(errorl3 * errorl3)), sqrt((errorh2 * errorh2)+(errorh * errorh)+(errorh3 * errorh3)));
    }
  }
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_comb5(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2, TGraphAsymmErrors* graph3, TGraphAsymmErrors* graph4, TGraphAsymmErrors * graph5) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, errorl, errorh;
    double errorl2, errorh2;
    double errorl3, errorh3;
    double errorl4, errorh4;
    double errorl5, errorh5;
    graph1->GetPoint(ibin, x, y);
    g_un->SetPoint(ibin, x, y);
    errorl = graph1->GetErrorYlow(ibin);
    errorh = graph1->GetErrorYhigh(ibin);
    errorl2 = graph2->GetErrorYlow(ibin);
    errorh2 = graph2->GetErrorYhigh(ibin);
    errorl3 = graph3->GetErrorYlow(ibin);
    errorh3 = graph3->GetErrorYhigh(ibin);
    errorl4 = graph4->GetErrorYlow(ibin);
    errorh4 = graph4->GetErrorYhigh(ibin);
    errorl5 = graph5->GetErrorYlow(ibin);
    errorh5 = graph5->GetErrorYhigh(ibin);
    double errorltot = sqrt((errorl * errorl)+(errorl2 * errorl2)+(errorl3 * errorl3)+(errorl4 * errorl4)+(errorl5 * errorl5));
    double errorhtot = sqrt((errorh * errorh)+(errorh2 * errorh2)+(errorh3 * errorh3)+(errorh4 * errorh4)+(errorh5 * errorh5));
    // cout<<ibin<<" "<<errorl/y<<" "<<errorl2/y<<" "<<errorl3/y<<" "<<errorl4/y<<" "<<errorl5/y<<endl;
    // cout<<ibin<<" "<<errorh/y<<" "<<errorh2/y<<" "<<errorh3/y<<" "<<errorh4/y<<" "<<errorh5/y<<endl;
    g_un->SetPointError(ibin, 0, 0, errorltot, errorhtot);
  }
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_toyvariation(vector<TGraphAsymmErrors*> graph, bool useMedian = 0, int part = 1, bool eff = 0) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y;
    vector<double> phis;
    graph[0]->GetPoint(ibin, x, y);
    phis.push_back(y);
    int nt = Ntoys;
    if (eff) nt = Ntoys2;
    uint gs = 1 + (nt * (part - 1));
    for (size_t i = gs; i < gs + nt; i++) {
      graph[i]->GetPoint(ibin, x, y);
      phis.push_back(y);
    }
    g_un->SetPoint(ibin, x, TMath::Mean(phis.begin(), phis.end()));
    double rms = TMath::RMS(phis.begin(), phis.end());
    g_un->SetPointError(ibin, 0, 0, rms, rms);
    if (useMedian && !eff) {
      std::sort(phis.begin(), phis.end());
      double med = phis[(nt - 1) / 2];
      int idx68 = 0.1587 * (nt) - 0.5;
      double min68 = phis[idx68];
      double max68 = phis[(nt) - idx68];
      g_un->SetPoint(ibin, x, med);
      g_un->SetPointError(ibin, 0, 0, med - min68, max68 - med);
    }
  }
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_toymc(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal, bool isBG = 0, int part = 1, bool eff = 0) {
  TGraphAsymmErrors* g_toy = CalcTotalSysU_toyvariation(graph, 0, part, eff);
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, xmc, ymc;
    graph_nominal->GetPoint(ibin, x, y);
    g_un->SetPoint(ibin, x, y);
    g_toy->GetPoint(ibin, xmc, ymc);
    double eh = ymc + g_toy->GetErrorYhigh(ibin);
    double el = ymc - g_toy->GetErrorYlow(ibin);
    eh = eh - y;
    el = el - y;
    // cout << ibin << " " << ymc << " " << eh << " " << el << " " << y << endl;
    if (eh < 0) eh = 0;
    if (el > 0) el = 0;
    //cout << ibin << " " << ymc << " " << eh << " " << el << " " << y << endl;
    if (isBG) {
      double ymc1, ymc2;
      double ymc3, ymc4;
      graph[graph.size() - 1]->GetPoint(ibin, x, ymc1);
      graph[graph.size() - 2]->GetPoint(ibin, x, ymc2);
      graph[graph.size() - 3]->GetPoint(ibin, x, ymc3);
      graph[graph.size() - 4]->GetPoint(ibin, x, ymc4);
      if (ymc1 > y) eh = sqrt(pow(eh, 2) + pow((ymc1 - y), 2));
      if (ymc2 > y) eh = sqrt(pow(eh, 2) + pow((ymc2 - y), 2));
      if (ymc3 > y) eh = sqrt(pow(eh, 2) + pow((ymc3 - y), 2));
      if (ymc4 > y) eh = sqrt(pow(eh, 2) + pow((ymc4 - y), 2));
      if (ymc1 < y) el = sqrt(pow(el, 2) + pow((ymc1 - y), 2));
      if (ymc2 < y) el = sqrt(pow(el, 2) + pow((ymc2 - y), 2));
      if (ymc3 < y) el = sqrt(pow(el, 2) + pow((ymc3 - y), 2));
      if (ymc4 < y) el = sqrt(pow(el, 2) + pow((ymc4 - y), 2));
    }
    g_un->SetPointError(ibin, 0, 0, el, eh);
  }
  delete g_toy;
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_fsr(TGraphAsymmErrors* graph, TGraphAsymmErrors * graph_nominal) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, xtemp, ytemp, error;
    graph_nominal->GetPoint(ibin, x, y);
    // g_un->SetPoint(ibin,x,y);
    g_un->SetPoint(ibin, x, y);
    graph->GetPoint(ibin, xtemp, ytemp);
    if (xtemp != x) cout << "This is really weird and wrong" << endl;
    error = fabs(ytemp - y);
    g_un->SetPointError(ibin, 0, 0, error, error);
  }
  return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_pileup(TGraphAsymmErrors* graph_1, TGraphAsymmErrors* graph_2, TGraphAsymmErrors * graph_nominal) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, xtemp1, ytemp1, xtemp2, ytemp2;
    double errorh = 0;
    double errorl = 0;
    graph_nominal->GetPoint(ibin, x, y);
    // g_un->SetPoint(ibin,x,y);
    g_un->SetPoint(ibin, x, y);
    graph_1->GetPoint(ibin, xtemp1, ytemp1);
    graph_2->GetPoint(ibin, xtemp2, ytemp2);
    if (xtemp1 != x) cout << "This is really weird and wrong" << endl;
    if (xtemp2 != x) cout << "This is really weird and wrong" << endl;
    if (y - ytemp1 > 0) errorl = fabs(ytemp1 - y);
    else errorh = fabs(ytemp1 - y);
    if (y - ytemp2 > 0) {
      if (errorl == 0) errorl = fabs(ytemp2 - y);
      else {
        if (errorl < fabs(ytemp2 - y)) errorl = fabs(ytemp2 - y);
        // cout << "pile-up: both errors on same side:" << y - ytemp1 << " " << y - ytemp2 << endl;
      }
    } else {
      if (errorh == 0) errorh = fabs(ytemp2 - y);
      else {
        if (errorh < fabs(ytemp2 - y)) errorh = fabs(ytemp2 - y);
        // cout << "pile-up: both errors on same side:" << y - ytemp1 << " " << y - ytemp2 << endl;
      }
    }
    g_un->SetPointError(ibin, 0, 0, errorl, errorh);
  }
  return g_un;
}

TGraphAsymmErrors * GetDataFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist, bool doNorm = 0, bool doLumi = 1) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  uint nun = graph.size();
  ofstream outputfile;
  std::string textname = "Output/Table_Un_";
  textname += Tag;
  if (doNorm) textname += "Norm_";
  else textname += "Abs_";
  if (doMG) textname += "MG_";
  else textname += "PH_";
  if (elec == 0)textname += "Dressed.txt";
  if (elec == 1)textname += "Born.txt";
  if (elec == 2)textname += "Naked.txt";
  outputfile.open(textname.c_str());
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double h_un = 0;
    double h_pdf = 0;
    double h_pt = 0;
    double l_pt = 0;
    double h_eff = 0;
    double h_bg = 0;
    double h_pu = 0;
    double h_mc = 0;
    double h_uns = 0;
    double l_eff = 0;
    double l_bg = 0;
    double l_pu = 0;
    double l_mc = 0;
    double l_pdf = 0; //Set but not used, delete?
    double x, y, errorl, errorh;
    graph[0]->GetPoint(ibin, x, y);
    double x_nom = x;
    double y_nom = y;
    ;
    double error_sys_max = 0;
    double error_sys_min = 0;
    double error_sys_tot = 0;
    double error_sys2_max = 0;
    double error_sys2_min = 0;
    double error_sys2_tot = 0;
    for (uint i = 0; i < nun; i++) {
      graph[i]->GetPoint(ibin, x, y);
      errorl = graph[i]->GetErrorYlow(ibin);
      errorh = graph[i]->GetErrorYhigh(ibin);
      double error_max = errorh / y;
      double error_min = errorl / y;
      double error_tot = (error_max + error_min) / 2.;
      error_sys_max += error_max*error_max;
      error_sys_min += error_min*error_min;
      error_sys_tot += error_tot*error_tot;
      if (slist[i] != "unfolding") {
        error_sys2_max += error_max*error_max;
        error_sys2_min += error_min*error_min;
        error_sys2_tot += error_tot*error_tot;
      }
      if (slist[i] == "unfolding") {
        h_un = error_max;
      }
      if (slist[i] == "eff") {
        h_eff = error_max;
        l_eff = error_min;
        x_nom = x;
        y_nom = y;
      }
      if (slist[i] == "bg") {
        h_bg = error_max;
        l_bg = error_min;
      }
      if (slist[i] == "pileup") {
        h_pu = error_max;
        l_pu = error_min;
      }
      if (slist[i] == "mcstat") {
        h_mc = error_max;
        l_mc = error_min;
      }
      if (slist[i] == "cteq") {
        h_pdf = error_max;
        l_pdf = error_min;
      }
      if (slist[i] == "pt") {
        h_pt = error_max;
        l_pt = error_min;
      }
      if (slist[i] == "unsyst") {
        h_uns = error_max;
      }
      //cout<<error_sys_min<<endl;
    }
    if (!doNorm && doLumi) error_sys_max += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys_min += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys_tot += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys2_max += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys2_min += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys2_tot += 0.026 * 0.026;
    //cout<<error_sys_min<<endl;

    if (ibin%nphistar==0){
      outputfile << std::fixed;
      outputfile << std::setprecision(1);
      outputfile << "\n";
      outputfile <<  yBins[ibin/nphistar] << "|Z(y)|" << yBins[(ibin/nphistar) + 1]<<"\n";
      outputfile << "\n";
      outputfile << "$\\phi^*$ range & Unfolding Syst & MC stat. & Pile-up & Background & Energy scale & Efficiencies ";
      if (!doMG) outputfile << "& PDF ";
      outputfile << "& Total syst. & Stat. & Total \\\\ \\hline" << "\n";
    }

    g_un->SetPoint(ibin, x_nom, y_nom);
    g_un->SetPointError(ibin, 0, 0, y_nom * sqrt(error_sys_min), y_nom * sqrt(error_sys_max));
    outputfile << std::fixed;
    outputfile << std::setprecision(3) << phistarBins[ibin%nphistar] << "-" << phistarBins[(ibin%nphistar) + 1] << " & " << std::setprecision(2) <<  h_uns * 100. << " & " << (h_mc + l_mc) / 2. * 100. << " & " << (h_pu + l_pu) / 2. * 100. << " & " << (h_bg + l_bg) / 2. * 100. << " & " << (h_pt + l_pt) / 2. * 100. << " & " << (h_eff + l_eff) / 2. * 100.;
    if (!doMG) outputfile << " & " << (h_pdf + l_pdf) / 2. * 100.;
    outputfile << " & " << sqrt(error_sys2_tot)*100 << " & " << h_un * 100. << " & " << sqrt(error_sys_tot)*100 << " \\\\ \\hline" << "\n";
  }
  outputfile.close();
  return g_un;
}
TGraphAsymmErrors* GetDataFinal(vector<TMatrixD> cov, vector<std::string> slist, TGraphAsymmErrors* graph, bool doNorm = 0, bool doLumi = 1) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  uint nun = slist.size();
  ofstream outputfile;
  std::string textname = "Output/TableM_Un_";
  textname += Tag;
  if (doNorm) textname += "Norm_";
  else textname += "Abs_";
  if (doMG) textname += "MG_";
  else textname += "PH_";
  if (elec == 0)textname += "Dressed.txt";
  if (elec == 1)textname += "Born.txt";
  if (elec == 2)textname += "Naked.txt";
  outputfile.open(textname.c_str());
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double h_un = 0;
    double h_pdf = 0;
    double h_pt = 0;
    double h_eff = 0;
    double h_bg = 0;
    double h_pu = 0;
    double h_mc = 0;
    double h_uns = 0;
    double x, y;
    graph->GetPoint(ibin, x, y);
    double error_sys_tot = 0;
    double error_sys2_tot = 0;
    for (uint i = 0; i < nun; i++) {
      double error = cov[i](ibin,ibin)/(y*y);
      error_sys_tot += error;
      if (slist[i] != "unfolding") error_sys2_tot += error;
      if (slist[i] == "unfolding") h_un = sqrt(error);
      if (slist[i] == "eff")       h_eff= sqrt(error);
      if (slist[i] == "bg")        h_bg = sqrt(error);
      if (slist[i] == "pileup")    h_pu = sqrt(error);
      if (slist[i] == "mcstat")    h_mc = sqrt(error);
      if (slist[i] == "cteq")      h_pdf= sqrt(error);
      if (slist[i] == "pt")        h_pt = sqrt(error);
      if (slist[i] == "unsyst")    h_uns= sqrt(error);
    }
    if (!doNorm && doLumi) error_sys_tot += 0.026 * 0.026;
    if (!doNorm && doLumi) error_sys2_tot += 0.026 * 0.026;

    if (ibin%nphistar==0){
      outputfile << std::fixed;
      outputfile << std::setprecision(1);
      outputfile << "\n";
      outputfile <<  yBins[ibin/nphistar] << "|Z(y)|" << yBins[(ibin/nphistar) + 1]<<"\n";
      outputfile << "\n";
      outputfile << "$\\phi^*$ range & Unfolding Syst & MC stat. & Pile-up & Background & Energy scale & Efficiencies ";
      if (!doMG) outputfile << "& PDF ";
      outputfile << "& Total syst. & Stat. & Total \\\\ \\hline" << "\n";
    }

    g_un->SetPoint(ibin, x, y);
    g_un->SetPointError(ibin, 0, 0, y * sqrt(error_sys_tot), y* sqrt(error_sys_tot));
    outputfile << std::fixed;
    outputfile << std::setprecision(3) << phistarBins[ibin%nphistar] << "-" << phistarBins[(ibin%nphistar) + 1] << " & " << std::setprecision(2) <<  h_uns * 100. << " & " << h_mc * 100. << " & " << h_pu * 100. << " & " << h_bg * 100. << " & " << h_pt * 100. << " & " << h_eff * 100.;
    if (!doMG) outputfile << " & " << h_pdf * 100.;
    outputfile << " & " << sqrt(error_sys2_tot)*100 << " & " << h_un * 100. << " & " << sqrt(error_sys_tot)*100 << " \\\\ \\hline" << "\n";
  }
  outputfile.close();
  return g_un;
}


void PrintBG(TH1D* Data, TH1D* h_tt, TH1D* h_tautau, TH1D* h_tbarw, TH1D* h_tw, TH1D* h_ww, TH1D* h_wz, TH1D* h_zz, TH1D* h_QCD) {
  double data_sel = Data->GetSumOfWeights();
  double tt_sel = h_tt->GetSumOfWeights();
  double tautau_sel = h_tautau->GetSumOfWeights();
  double tbarw_sel = h_tbarw->GetSumOfWeights();
  double tw_sel = h_tw->GetSumOfWeights();
  double ww_sel = h_ww->GetSumOfWeights();
  double wz_sel = h_wz->GetSumOfWeights();
  double zz_sel = h_zz->GetSumOfWeights();
  double qcd_sel = h_QCD->GetSumOfWeights();
  double t_bg = tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel;
  cout << "Weighted number of events:" << endl;
  cout << "Data: " << data_sel << " tt: " << tt_sel << " tautau: " << tautau_sel << " tbarw: " << tbarw_sel << " tw: " << tw_sel << " singletop: " << tbarw_sel + tw_sel << " ww: " << ww_sel << " wz: " << wz_sel << " zz: " << zz_sel << " qcd: " << qcd_sel << endl;
  cout << "ratio:" << " tt: " << tt_sel * 100. / data_sel << " tautau: " << tautau_sel * 100. / data_sel << " tbarw: " << tbarw_sel * 100. / data_sel << " tw: " << tw_sel * 100. / data_sel << " singletop: " << (tbarw_sel + tw_sel)*100. / data_sel << " ww: " << ww_sel * 100. / data_sel << " wz: " << wz_sel * 100. / data_sel << " zz: " << zz_sel * 100. / data_sel << " qcd: " << qcd_sel * 100. / data_sel << "data: " << (data_sel - t_bg)*100. / data_sel << endl;
  cout << "ratio of bg:" << " tt: " << tt_sel * 100. / t_bg << " tautau: " << tautau_sel * 100. / t_bg << " tbarw: " << tbarw_sel * 100. / t_bg << " tw: " << tw_sel * 100. / t_bg << " singletop: " << (tbarw_sel + tw_sel)*100. / t_bg << " ww: " << ww_sel * 100. / t_bg << " wz: " << wz_sel * 100. / t_bg << " zz: " << zz_sel * 100. / t_bg << " qcd: " << qcd_sel * 100. / t_bg << endl;
  cout << "total: " << t_bg << " " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel) << " persentage: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel)*100. / data_sel << endl;
  cout << "wz+zz: " << (wz_sel + zz_sel) << " persentage: " << (wz_sel + zz_sel)*100. / data_sel << " of background: " << (wz_sel + zz_sel)*100. / (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel) << endl;
  cout << "e+mu: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel) << " persentage: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel)*100. / data_sel << " of background: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel)*100. / (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel) << endl;
}

void GetToyResponse(vector<RooUnfoldResponse*> &BinM, TH2D * BinMigration) {
  // RooUnfoldResponse* BinM1 =new RooUnfoldResponse (h_reco,h_gen,BinMigration);
  RooUnfoldResponse* BinM1 = new RooUnfoldResponse(0, 0, BinMigration);
  BinM.push_back(BinM1);
  double x;
  TH2D* BinMigrationtemp = new TH2D("BinMigration", "BinMigration", nbins, 0, nbins, nbins, 0, nbins);
  BinMigrationtemp->Sumw2();
  for (int i = 0; i < Ntoys; i++) {
    for (uint j = 0; j < nbins; j++) {
      for (uint k = 0; k < nbins; k++) {
        double mean = BinMigration->GetBinContent(j + 1, k + 1);
        if (mean == 0) continue;
        double error = BinMigration->GetBinError(j + 1, k + 1);
	x = gRandom->Gaus(mean, error);
	if (x<0) x=0;
        BinMigrationtemp->SetBinContent(j + 1, k + 1, x);
        BinMigrationtemp->SetBinError(j + 1, k + 1, error);
      }
    }
    RooUnfoldResponse* BinM1temp = new RooUnfoldResponse(0, 0, BinMigrationtemp);
    BinM.push_back(BinM1temp);
  }
  delete BinMigrationtemp;
  return;
}

vector<TH1D*> RemoveBG(TH1D* Data, TH1D* bg_sf, TH1D* bg_ss, vector<TH1D*> h_tt, vector<TH1D*> h_tautau, vector<TH1D*> h_tbarw, vector<TH1D*> h_tw, vector<TH1D*> h_ww, vector<TH1D*> h_wz, vector<TH1D*> h_zz) {
  vector<TH1D*> h_data;
  for (uint i = 0; i < h_tt.size(); i++) {
    TH1D* bgtemp = (TH1D*) h_tt[i]->Clone();
    bgtemp->Add(h_tautau[i], 1.0);
    bgtemp->Add(h_tbarw[i], 1.0);
    bgtemp->Add(h_tw[i], 1.0);
    bgtemp->Add(h_ww[i], 1.0);
    bgtemp->Multiply(bg_sf);
    TH1D* datatemp = (TH1D*) Data->Clone();
    datatemp->Add(bgtemp, -1.0);
    datatemp->Add(h_wz[i], -1.0);
    datatemp->Add(h_zz[i], -1.0);
    datatemp->Add(bg_ss, -1.0);
    h_data.push_back(datatemp);
    delete bgtemp;
  }
  return h_data;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, vector<TH1D*> mc_truegen) {
  vector<TH1D*> h_eff;
  for (uint i = 0; i < mc_truereco.size(); i++) {
    TH1D* h_efftemp = new TH1D("h_eff", "h_eff", nbins, 0, nbins);
    TH1D* gentemp = (TH1D*) mc_truegen[i]->Clone();
    h_efftemp->Divide(mc_truereco[i], gentemp, 1., 1., "B");
    h_eff.push_back(h_efftemp);
    delete gentemp;
  }
  return h_eff;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, TH1D * mc_truegen) {
  vector<TH1D*> h_eff;
  for (uint i = 0; i < mc_truereco.size(); i++) {
    TH1D* h_efftemp = new TH1D("h_eff", "h_eff", nbins, 0, nbins);
    TH1D* gentemp = (TH1D*) mc_truegen->Clone();
    h_efftemp->Divide(mc_truereco[i], gentemp, 1., 1., "B");
    h_eff.push_back(h_efftemp);
    delete gentemp;
  }
  return h_eff;
}

TGraphAsymmErrors * ConvertToTGraph(TH1D * h) {
  TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    g->SetPoint(ibin, ibin, h->GetBinContent(ibin + 1));
    g->SetPointError(ibin, 0, 0, h->GetBinError(ibin + 1), h->GetBinError(ibin + 1));
  }
  return g;
}

vector<TGraphAsymmErrors *> CreateCopy(vector<TGraphAsymmErrors *> graphvec) {
  vector<TGraphAsymmErrors *> newvec;
  for (int i = 0; i < ((int) graphvec.size()); i++) {//just getting rid of the warnings
    TGraphAsymmErrors *temp = (TGraphAsymmErrors *) graphvec[i]->Clone();
    newvec.push_back(temp);
  }
  return newvec;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(vector<RooUnfoldResponse*> BinM, vector<TH1D*> h_data, vector<TH1D*> h_eff, int getFirst = 0) {
  cout << "Unfolding:" << getFirst << endl;
  vector<TGraphAsymmErrors *> g_data;
  size_t n = BinM.size();
  if (getFirst == 1)n = h_data.size();
  for (size_t u = 0; u < n; u++) {
    if (u%100==0) cout<<"Unfolding "<<u<<" of "<<n<<endl;
    uint b = u;
    uint e = u;
    uint d = u;
    if (getFirst == 1) {
      b = 0;
      e = 0;
    }
    if (getFirst == 2) {
      d = 0;
      e = 0;
    }
    if (getFirst == 3) {
      d = 1;
      e = 1;
    }
    // cout<<u<<" "<<b<<" "<<e<<" "<<d<<endl;
    RooUnfoldBayes unfoldBay_data(BinM[b], h_data[d], 4);
    if (u!=0)unfoldBay_data.SetVerbose(0);
    TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
    NoNeg(h_BinMBay_data);
    TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
    h_Unfolded_temp->Divide(h_BinMBay_data, h_eff[e]);
    TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
    g_data.push_back(g_data_temp);
    delete h_BinMBay_data;
    delete h_Unfolded_temp;
  }
  return g_data;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(vector<RooUnfoldResponse*> BinM, TH1D* h_data, vector<TH1D*> h_eff) {
  vector<TGraphAsymmErrors *> g_data;
  size_t n = BinM.size();
  for (size_t u = 0; u < n; u++) {
    if (u%100==0) cout<<"Unfolding "<<u<<" of "<<n<<endl;
    RooUnfoldBayes unfoldBay_data(BinM[u], h_data, 4);
    unfoldBay_data.SetVerbose(0);
    TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
    NoNeg(h_BinMBay_data);
    TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
    h_Unfolded_temp->Divide(h_BinMBay_data, h_eff[u]);
    TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
    g_data.push_back(g_data_temp);
    delete h_BinMBay_data;
    delete h_Unfolded_temp;
  }
  return g_data;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(RooUnfoldResponse* BinM, TH1D* h_data, TH1D * h_eff, TMatrixD& CovM) {
  vector<TGraphAsymmErrors *> g_data;
  RooUnfoldBayes unfoldBay_data(BinM, h_data, 4);
  //  unfoldBay_data.SetVerbose(0);
  TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
  NoNeg(h_BinMBay_data);
  TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
  h_Unfolded_temp->Divide(h_BinMBay_data, h_eff);
  TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
  g_data.push_back(g_data_temp);
  TMatrixD CovM_temp= unfoldBay_data.Ereco();
  for(uint j = 0; j < nbins; j++) {
    for(uint k = 0; k < nbins; k++) {
      CovM(j,k)=CovM_temp(j,k)*(1/h_eff->GetBinContent(j+1))*(1/h_eff->GetBinContent(k+1));
    }
  }
  //test symmetry.
  for(uint j = 0; j < nbins; j++) {
    for(uint k = j+1; k < nbins; k++) {
      if (fabs(CovM(j,k)/CovM(k,j)-1)>0.00001) cout<<"Stat:Covarience matrix is not symmetric, "<<j<<" "<<k<<":"<<CovM(j,k)<<" != "<<CovM(k,j)<<endl;
    }
  }
  delete h_BinMBay_data;
  delete h_Unfolded_temp;
  return g_data;
}

vector<TH1D*> EmptyHVec(int n) {
  vector<TH1D*> h;
  for (int i = 0; i < n; i++) {
    TH1D *phistartemp = new TH1D("phistar", "phistar", nbins, 0, nbins);
    phistartemp->Sumw2();
    h.push_back(phistartemp);
  }
  return h;
}

vector<RooUnfoldResponse*> EmptyBinMVec(int n) {
  vector<RooUnfoldResponse*> BinM;
  TH1D *phistartemp = new TH1D("phistar", "phistar", nbins, 0, nbins);
  phistartemp->Sumw2();
  for (int i = 0; i < n; i++) {
    RooUnfoldResponse* responsetemp = new RooUnfoldResponse(phistartemp, phistartemp);
    BinM.push_back(responsetemp); 
  }
  delete phistartemp;
  return BinM;
}

/*ask about the blue code to Jeremey*/
void FillRecoEffFluc(int &idx, vector<TH2D*> ToyEff, vector<TH1D*> &h_eff, double weight, int bin, double E0_pt, double E0_eta, double E1_pt, double E1_eta, bool doBinM, vector<RooUnfoldResponse*> &BinM_eff, int bin_true = 0, int type = 0) {
  int Bin0 = ToyEff[0]->FindBin(fabs(E0_eta), E0_pt);
  int Bin1 = ToyEff[0]->FindBin(fabs(E1_eta), E1_pt);
  double Norm0 = ToyEff[0]->GetBinContent(Bin0);
  double Norm1 = ToyEff[0]->GetBinContent(Bin1);
  for (int t = 0; t < Ntoys2; t++) {
    double New0 = ToyEff[t + 1]->GetBinContent(Bin0);
    double New1 = ToyEff[t + 1]->GetBinContent(Bin1);
    double new_weight = weight * New0 * New1 / (Norm0 * Norm1);
    if (type == 1) new_weight = weight * New1 / Norm1;
    else if (type == 2) new_weight = weight * New0 / Norm0;
    // cout<<type<<" "<<weight<<" "<<new_weight<<endl;
    if (fabs((weight / new_weight) - 1) > 0.5) cout << Norm0 << " " << Norm1 << " " << New0 << " " << New1 << " " << weight / new_weight << endl;
    if (!doBinM) {
      h_eff[idx] ->Fill(bin, new_weight);
    } else {
      h_eff[idx] ->Fill(bin_true, new_weight);
      BinM_eff[idx]->Fill(bin, bin_true, new_weight);
    }
    idx++;
  }
}

void FillTrigEffFluc(int &idx, vector<TH2D*> ToyEffMC, vector<TH2D*> ToyEffData, vector<TH1D*> &h_eff, double weight, int bin, double E0_pt, double E0_eta, double E1_pt, double E1_eta, bool doBinM, vector<RooUnfoldResponse*> &BinM_eff, int bin_true = 0) {
  int Bin0 = ToyEffMC[0]->FindBin(E0_eta, E0_pt);
  int Bin1 = ToyEffMC[0]->FindBin(E1_eta, E1_pt);
  double MCNorm0 = ToyEffMC[0]->GetBinContent(Bin0);
  double MCNorm1 = ToyEffMC[0]->GetBinContent(Bin1);
  double DataNorm0 = ToyEffData[0]->GetBinContent(Bin0);
  double DataNorm1 = ToyEffData[0]->GetBinContent(Bin1);
  for (int t = 0; t < Ntoys2; t++) {
    double MCNew0 = ToyEffMC[t + 1]->GetBinContent(Bin0);
    double MCNew1 = ToyEffMC[t + 1]->GetBinContent(Bin1);
    double new_weight = weight;
    if (fabs(E0_eta) > 2.1 || E0_pt < 30) {
      new_weight = new_weight * (MCNorm1 / MCNew1);
      // cout<<"2:"<< E1_pt<<" "<<E1_eta<<" "<<MCNorm1<<" "<<BinX1<<" "<<BinY1<<endl;
    } else if (fabs(E1_eta) > 2.1 || E1_pt < 30) {
      new_weight = new_weight * (MCNorm0 / MCNew0);
      // cout<<"2:"<< E0_pt<<" "<<E0_eta<<" "<<MCNorm0<<" "<<BinX0<<" "<<BinY0<<endl;
    } else {
      double old_w = (1. - (1. - DataNorm0)*(1. - DataNorm1)) / (1. - (1. - MCNorm0)*(1. - MCNorm1));
      double new_w = (1 - (1 - DataNorm0)*(1 - DataNorm1)) / (1 - (1 - MCNew0)*(1 - MCNew1));
      new_weight = new_weight * new_w / old_w;
      // cout<<"2:"<< E0_pt<<" "<<E0_eta<<" "<<MCNorm0<<" "<<BinX1<<" "<<BinY1<<endl;
      // cout<<"2:"<< E1_pt<<" "<<E1_eta<<" "<<MCNorm1<<" "<<BinX1<<" "<<BinY1<<endl;
    }
    if (fabs((weight / new_weight) - 1) > 0.5) cout << weight / new_weight << endl;
    if (!doBinM) {
      h_eff[idx]->Fill(bin, new_weight);
    } else {
      h_eff[idx] ->Fill(bin_true, new_weight);
      BinM_eff[idx]->Fill(bin, bin_true, new_weight);
    }
    idx++;
  }
  for (int t = 0; t < Ntoys2; t++) {
    double DataNew0 = ToyEffData[t + 1]->GetBinContent(Bin0);
    double DataNew1 = ToyEffData[t + 1]->GetBinContent(Bin1);
    double new_weight = weight;
    if (fabs(E0_eta) > 2.1 || E0_pt < 30) {
      // cout<<"3:"<< E1_pt<<" "<<E1_eta<<" "<<DataNorm1<<endl;
      new_weight = new_weight * (DataNew1 / DataNorm1);
    } else if (fabs(E1_eta) > 2.1 || E1_pt < 30) {
      // cout<<"3:"<< E0_pt<<" "<<E0_eta<<" "<<DataNorm1<<endl;
      new_weight = new_weight * (DataNew0 / DataNorm0);
    } else {
      double old_w = (1. - (1. - DataNorm0)*(1. - DataNorm1)) / (1. - (1. - MCNorm0)*(1. - MCNorm1));
      double new_w = (1 - (1 - DataNew0)*(1 - DataNew1)) / (1 - (1 - MCNorm0)*(1 - MCNorm1));
      new_weight = new_weight * new_w / old_w;
      // cout<<"3:"<< E0_pt<<" "<<E0_eta<<" "<<DataNorm1<<endl;
      // cout<<"3:"<< E1_pt<<" "<<E1_eta<<" "<<DataNorm1<<endl;
    }
    if (fabs((weight / new_weight) - 1) > 0.5) cout << weight / new_weight << endl;
    if (!doBinM) {
      h_eff[idx]->Fill(bin, new_weight);
    } else {
      h_eff[idx] ->Fill(bin_true, new_weight);
      BinM_eff[idx]->Fill(bin, bin_true, new_weight);
    }
    idx++;
  }
}

void GetBG(std::string FileName, double sampleweight, vector<TH1D*> &h_eff, vector<TH1D*> &h_fsr_pileup, vector<TH2D*> h_ToyEffSF, vector<TH2D*> h_ToyEffMSF, vector<TH2D*> h_ToyEffTSF, vector<TH2D*> h_ToyEffTMC, vector<TH2D*> h_ToyEffTData) {
  gErrorIgnoreLevel = kError;
  cout << "reading data for " << FileName << endl;
  TChain* t = new TChain(reco_name.c_str(), reco_name.c_str());
  t->Add(FileName.c_str());
  t->SetBranchStatus("event_info", 0); //to disable all branches
  t->SetBranchStatus("truth", 0); //to disable all branches
  TBranch *b_reco = t->GetBranch("reco");
  TLeaf *l_phistar = b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_e0_pt = b_reco->GetLeaf("e_pt0");
  TLeaf *l_e0_eta = b_reco->GetLeaf("e_eta0");
  TLeaf *l_e1_pt = b_reco->GetLeaf("e_pt1");
  TLeaf *l_e1_eta = b_reco->GetLeaf("e_eta1");
  TLeaf *l_YZ = b_reco->GetLeaf("z_y");
  int nweights;
  t->SetBranchAddress("weight_size", &nweights);
  t->GetEntry(0);
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights", &weights);
  t->SetBranchAddress("weight_ids", &weightid);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr", &weight_fsr);
  h_eff = EmptyHVec(Ntoys2 * 5 + 1);
  h_fsr_pileup = EmptyHVec(3);
  cout << "Entries: " << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
    // if (!(i % 10000))cout << " Entry number " << i << endl;
    t->GetEntry(i);
    double E0_pt = l_e0_pt ->GetValue();
    double E0_eta = l_e0_eta ->GetValue();
    double E1_pt = l_e1_pt ->GetValue();
    double E1_eta = l_e1_eta ->GetValue();
    double phistar = l_phistar->GetValue();
    double Z_Y = l_YZ->GetValue();
    int bin=GetBin(phistar,Z_Y);
    double weight = sampleweight;
    double weightpu_0 = 0;
    double weightpu_p = 0;
    double weightpu_m = 0;
    for (int w = 0; w < nweights; w++) {
      if (weightid[w] == 1 || weightid[w] == 2 || weightid[w] == 12 || weightid[w] == 13 || weightid[w] == 20 || weightid[w] == 30) {
        weight = weight * weights[w];
      }
      if (weightid[w] == 2) weightpu_0 = weights[w];
      if (weightid[w] == 3) weightpu_p = weights[w];
      if (weightid[w] == 4) weightpu_m = weights[w];
    }
    // if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;
    h_eff[0]->Fill(bin, weight);
    h_fsr_pileup[0]->Fill(bin, weight * weight_fsr);
    if (weightpu_0 != 0) {
      h_fsr_pileup[1]->Fill(bin, weight * weightpu_p / weightpu_0);
      h_fsr_pileup[2]->Fill(bin, weight * weightpu_m / weightpu_0);
    }
    int idx = 1;
    vector<RooUnfoldResponse*> dummy;
    FillRecoEffFluc(idx, h_ToyEffSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy);
    FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy, 0, 1);
    FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy, 0, 2);
    if (idx == 1) cout << "idx did not increase when it should have" << endl;
    FillTrigEffFluc(idx, h_ToyEffTMC, h_ToyEffTData, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy);
  }
  cout << "done reading data for " << FileName << endl;
  // delete l_phistar;
  // delete l_e0_pt;
  // delete l_e0_eta;
  // delete l_e1_pt;
  // delete l_e1_eta;
  // delete l_YZ;
  // delete b_reco;
  delete t;
}

TH1D * GetData(double sampleweight, int up = 0) {
  cout << "reading data" << ny<<" "<<nphistar<<" "<<nbins<< endl;
  TH1D *h_phistar = new TH1D("phistar", "phistar", nbins, 0, nbins);
  h_phistar->Sumw2();
  std::string name = reco_name;
  if (up < 0) name = reco_name_en_l;
  if (up > 0) name = reco_name_en_h;
  TChain* t = new TChain(name.c_str(), name.c_str());
  if (up == 0) t->Add(File_Data.c_str());
  else if (up > 0) t->Add(File_Data_Pt_H.c_str());
  else t->Add(File_Data_Pt_L.c_str());
  TBranch *b_reco = t->GetBranch("reco");
  TLeaf *l_phistar = b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_ZY = b_reco->GetLeaf("z_y");
  cout << "Entries: " << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
    // if (!(i % 10000))cout << " entry2 number " << i << endl;
    //    if (i%100000==0)cout<<i<<endl;
    t->GetEntry(i);
    double Z_Y = l_ZY->GetValue();
    double phistar = l_phistar->GetValue();
    int bin = GetBin(phistar,Z_Y);
    h_phistar->Fill(bin, sampleweight);
  }
  // delete l_phistar;
  // delete l_ZY;
  // delete b_reco;
  delete t;
  cout << "filled data phistar histogram" << endl;
  return h_phistar;
}

TH1D * GetDataPt(double sampleweight, int up = 0) {
  cout << "reading data" << endl;
  TH1D *h_phistar = new TH1D("phistar", "phistar", nbins, 0, nbins);
  h_phistar->Sumw2();
  std::string name = reco_name;
  if (up < 0) name = reco_name_en_l;
  // if (up>0) name=reco_name_en_h;
  TChain* t = new TChain(name.c_str(), name.c_str());
  if (up >= 0) t->Add(File_Data.c_str());
    // else if (up > 0) t->Add(File_Data_Pt_H.c_str());
  else t->Add(File_Data_Pt_L.c_str());
  TBranch *b_reco = t->GetBranch("reco");
  TLeaf *l_phistar = b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_e0pt = b_reco->GetLeaf("e_pt0");
  TLeaf *l_e0eta = b_reco->GetLeaf("e_eta0");
  TLeaf *l_e1pt = b_reco->GetLeaf("e_pt1");
  TLeaf *l_e1eta = b_reco->GetLeaf("e_eta1");
  TLeaf *l_e0_tight = b_reco->GetLeaf("t0tight");
  TLeaf *l_e1_tight = b_reco->GetLeaf("t1tight");
  TLeaf *l_ZY = b_reco->GetLeaf("z_y");
  cout << "Entries1: " << t->GetEntries() << endl;
  double pth = 30;
  double ptl = 20;
  double etah = 2.1;
  if (up > 0) {
    pth = pth * 1.003;
    ptl = ptl * 1.003;
  } else if (up < 0) {
    pth = pth * 0.997;
    ptl = ptl * 0.997;
  }
  int idx = 0;
  for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
    // if (!(i % 10000))cout << " entry number " << i << endl;
    t->GetEntry(i);
    bool E0Tight = l_e0_tight ->GetValue();
    bool E1Tight = l_e1_tight ->GetValue();
    double pt0 = l_e0pt->GetValue();
    double pt1 = l_e1pt->GetValue();
    double eta0 = l_e0eta->GetValue();
    double eta1 = l_e1eta->GetValue();
    double phistar = l_phistar->GetValue();
    double Z_Y = l_ZY->GetValue();
    int bin= GetBin(phistar, Z_Y);
    if (pt0 < ptl || pt1 < ptl) continue;
    if ((pt0 < pth || fabs(eta0) > etah || !E0Tight) && (pt1 < pth || fabs(eta1) > etah || !E1Tight)) continue;
    h_phistar->Fill(bin, sampleweight);
    idx++;
  }
  //  cout << "filled data phistar histogram, selected:" << idx << endl;
  // delete l_phistar;
  // delete l_e0pt;
  // delete l_e0eta;
  // delete l_e1pt;
  // delete l_e1eta;
  // delete l_e0_tight;
  // delete l_e1_tight;
  // delete l_ZY;
  // cout << "filled data phistar histogram, selected:" << idx << endl;
  // delete b_reco;
  // cout << "filled data phistar histogram, selected:" << idx << endl;
  delete t;
  cout << "filled data phistar histogram, selected:" << idx << endl;
  return h_phistar;
}

void GetGen(double sampleweight, TH1D* &h_phistar, vector<TH1D*> &h_cteq, vector<TH1D*> &h_fsr_pileup) {
  gErrorIgnoreLevel = kError;
  cout << "reading signal gen" << endl;
  h_phistar = new TH1D("phistar", "phistar", nbins, 0, nbins);
  h_phistar->Sumw2();
  TChain* t = new TChain(gen_name.c_str(), gen_name.c_str());
  if (doMG) {
    if (elec == 0) t->Add(File_Signal_gen.c_str());
    else if (elec == 1) t->Add(File_Signal_gen_born.c_str());
    else t->Add(File_Signal_gen_bare.c_str());
  } else {
    if (elec == 0) t->Add(File_Powheg_gen.c_str());
    else if (elec == 1) t->Add(File_Powheg_gen_born.c_str());
    else t->Add(File_Powheg_gen_bare.c_str());
  }
  ///t->SetBranchStatus("event_info", 0); //to disable all branches
  //  t->SetBranchStatus("reco", 0); //to disable all branches
  TBranch *b_truth = t->GetBranch("truth");
  TBranch *b_reco = t->GetBranch("reco");
  TLeaf * l_ZY = b_truth->GetLeaf("z_y");
  if (elec == 1) l_ZY = b_reco->GetLeaf("z_yBorn"); //TO BE CHANGED
  if (elec == 2) l_ZY = b_reco->GetLeaf("z_yNaked");
  TLeaf *l_phistar = b_truth->GetLeaf("z_phistar_dressed");
  if (elec == 1) l_phistar = b_truth->GetLeaf("z_phistar_born");
  if (elec == 2) l_phistar = b_truth->GetLeaf("z_phistar_naked");
  int nweights;
  int nwcteq;
  t->SetBranchAddress("weight_size", &nweights);
  t->SetBranchAddress("weight_cteq_size", &nwcteq);
  t->GetEntry(0);
  double weights[nweights];
  int weightid[nweights];
  double weights_cteq[nwcteq];
  t->SetBranchAddress("weights", &weights);
  t->SetBranchAddress("weight_ids", &weightid);
  t->SetBranchAddress("weights_cteq", &weights_cteq);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr", &weight_fsr);
  cout << "Entries: " << t->GetEntries() << endl;
  h_cteq = EmptyHVec(nwcteq);
  h_fsr_pileup = EmptyHVec(3);

  TH1D* h_PDFWeight = new TH1D("h_PDFWeight", "h_PDFWeight", 10000, 0, 100);
  h_PDFWeight->Sumw2();
  // cout << "Nweightcteq=" << nwcteq << endl;

  for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
    // for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double weight = sampleweight;
    double pdfnorm = weights_cteq[0];
    double weightpu_0 = 0;
    double weightpu_p = 0;
    double weightpu_m = 0;
    double phistar = l_phistar->GetValue();
    double Z_Y = l_ZY->GetValue();
    int bin=GetBin(phistar,Z_Y);
    //cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;
    for (int w = 0; w < nweights; w++) {
      if (weightid[w] == 1 || weightid[w] == 2) {
        weight = weight * weights[w];
      }
      //if (weightid[w]==1) {weight=weight*weights[w];}
      if (weightid[w] == 2) weightpu_0 = weights[w];
      if (weightid[w] == 3) weightpu_p = weights[w];
      if (weightid[w] == 4) weightpu_m = weights[w];
    }
    // if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;
    h_phistar->Fill(bin, weight);
    if (pdfnorm != 0) {
      for (int w = 0; w < nwcteq; w++) {
        h_cteq[w] ->Fill(bin, weight * weights_cteq[w] / pdfnorm);
        if (w != 0) {
          h_PDFWeight->Fill(weights_cteq[w] / pdfnorm, 1);
          // cout<<weights_cteq[w] / pdfnorm<<endl;
        }
      }
    }
    if (doMG) weight_fsr = 1;
    h_fsr_pileup[0]->Fill(bin, weight * weight_fsr);
    if (weightpu_0 != 0) {
      h_fsr_pileup[1]->Fill(bin, weight * weightpu_p / weightpu_0);
      h_fsr_pileup[2]->Fill(bin, weight * weightpu_m / weightpu_0);
    }
  }
  cout << "done reading signal gen" << endl;
  // TCanvas* PDFW = new TCanvas("PDFW", "PDFW", 800, 900);
  // PDFW->cd();
  // h_PDFWeight->SetStats(0);
  // h_PDFWeight->SetBit(TH1::kNoTitle, true);
  // h_PDFWeight->SetLineColor(1);
  // h_PDFWeight->Draw();
  // delete l_ZY;
  // delete l_phistar;
  // delete b_truth;
  delete t;
  return;
}

void GetBinM(double sampleweight, vector<RooUnfoldResponse*> &BinM_eff, vector<RooUnfoldResponse*> &BinM_mcstat, vector<RooUnfoldResponse*> &BinM_cteq, vector<RooUnfoldResponse*> &BinM_fsr_pileup, vector<TH1D*> &h_eff, vector<TH1D*> &h_cteq, vector<TH1D*> &h_fsr_pileup, vector<TH2D*> h_ToyEffSF, vector<TH2D*> h_ToyEffMSF, vector<TH2D*> h_ToyEffTSF, vector<TH2D*> h_ToyEffTMC, vector<TH2D*> h_ToyEffTData) {
  gErrorIgnoreLevel = kError;
  cout << "reading data for " << File_Signal_reco << " " << reco_name << endl;
  TChain* t;
  t = new TChain(reco_name.c_str(), reco_name.c_str());
  // else t= new TChain(reco_name2.c_str(), reco_name2.c_str());
  if (doMG){
    if (elec==0)t->Add(File_Signal_reco.c_str());
    if (elec==1)t->Add(File_Signal_reco_born.c_str());
    if (elec==2)t->Add(File_Signal_reco_bare.c_str());
  }  
  else {
    if (elec==0)t->Add(File_Powheg_reco.c_str());
    if (elec==1)t->Add(File_Powheg_reco_born.c_str());
    if (elec==2)t->Add(File_Powheg_reco_bare.c_str());
  }
  t->SetBranchStatus("event_info", 0); //to disable all branches
  TBranch *b_reco = t->GetBranch("reco");
  TBranch *b_truth = t->GetBranch("truth");
  TLeaf *l_phistar = b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_e0_pt = b_reco->GetLeaf("e_pt0");
  TLeaf *l_e0_eta = b_reco->GetLeaf("e_eta0");
  TLeaf *l_e1_pt = b_reco->GetLeaf("e_pt1");
  TLeaf *l_e1_eta = b_reco->GetLeaf("e_eta1");
  TLeaf *l_e0_tight = b_reco->GetLeaf("t0tight");
  TLeaf *l_ZY = b_reco->GetLeaf("z_y");
  TLeaf *l_YZTruth = b_truth->GetLeaf("z_y");
  if (elec == 1) l_YZTruth = b_reco->GetLeaf("z_yBorn"); //TO BE CHANGED
  if (elec == 2) l_YZTruth = b_reco->GetLeaf("z_yNaked");
  TLeaf *l_phistar_true = b_truth->GetLeaf("z_phistar_dressed");
  if (elec == 1) l_phistar_true = b_truth->GetLeaf("z_phistar_born");
  if (elec == 2) l_phistar_true = b_truth->GetLeaf("z_phistar_naked");
  int nweights;
  int nwcteq;
  t->SetBranchAddress("weight_size", &nweights);
  t->SetBranchAddress("weight_cteq_size", &nwcteq);
  t->GetEntry(0);
  double weights[nweights];
  int weightid[nweights];
  double weights_cteq[nwcteq];
  t->SetBranchAddress("weights", &weights);
  t->SetBranchAddress("weight_ids", &weightid);
  t->SetBranchAddress("weights_cteq", &weights_cteq);
  double weight_fsr;
  t->SetBranchAddress("weight_fsr", &weight_fsr);
  cout << "Entries: " << t->GetEntries() << endl;
  h_eff = EmptyHVec(Ntoys2 * 5 + 1);
  h_cteq = EmptyHVec(nwcteq);
  h_fsr_pileup = EmptyHVec(3);
  BinM_eff = EmptyBinMVec(Ntoys2 * 5 + 1);
  BinM_cteq = EmptyBinMVec(nwcteq);
  BinM_fsr_pileup = EmptyBinMVec(3);
  // TH1D *h_rec = new TH1D("phistar", "phistar", nbins, 0, nbins);
  // h_rec->Sumw2();
  TH2D* BinMigration = new TH2D("BinMigration", "BinMigration", nbins, 0, nbins, nbins, 0, nbins);
  BinMigration->Sumw2();

  // TH1D* h_PDFWeight = new TH1D("h_PDFWeight", "h_PDFWeight", 10000, 0, 100);
  // h_PDFWeight->Sumw2();

  for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
    // if (i%10000==0) cout<<i<<endl;
    // if (i>4740000) cout<<i<<endl;
    // Percent = 100 * Percent / (double(t->GetEntries()));
    t->GetEntry(i);
    bool E0Tight = 1;
    E0Tight = l_e0_tight ->GetValue();
    double E0_pt = l_e0_pt ->GetValue();
    double E0_eta = l_e0_eta ->GetValue();
    double E1_pt = l_e1_pt ->GetValue();
    double E1_eta = l_e1_eta ->GetValue();
    double phistar = l_phistar->GetValue();
    double phistar_true = l_phistar_true->GetValue();
    double Z_Y = l_ZY->GetValue();
    double Z_YTruth = l_YZTruth->GetValue();
    int bin=GetBin(phistar,Z_Y);
    int bin_true=GetBin(phistar_true,Z_YTruth);
    double pdfnorm = 0;
    if (!doMG)pdfnorm = weights_cteq[0];
    double weight = sampleweight;
    double weightpu_0 = 0;
    double weightpu_p = 0;
    double weightpu_m = 0;
    //cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;
    for (int w = 0; w < nweights; w++) {
      if (weightid[w] == 1 || weightid[w] == 2 || weightid[w] == 12 || weightid[w] == 13 || weightid[w] == 20 || weightid[w] == 30) {
        weight = weight * weights[w];
      }
      if (weightid[w] == 2) weightpu_0 = weights[w];
      if (weightid[w] == 3) weightpu_p = weights[w];
      if (weightid[w] == 4) weightpu_m = weights[w];
    }
    // if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;
    // if (weightpu_0 == 0) cout << "pile-up nominal weight 0" << endl;
    // int ii = 0; //just using it as a counter for a print out statement
    h_eff[0] ->Fill(bin_true, weight);
    BinM_eff[0] ->Fill(bin, bin_true, weight);
    //h_rec ->Fill(bin, weight);
    BinMigration ->Fill(bin, bin_true, weight);
    if (pdfnorm != 0) {
      for (int w = 0; w < nwcteq; w++) {
        //cout << w << endl;
        h_cteq[w] ->Fill(bin_true, weight * weights_cteq[w] / pdfnorm);
        BinM_cteq[w]->Fill(bin, bin_true, weight * weights_cteq[w] / pdfnorm);
        // if (w != 0) {
        //   h_PDFWeight->Fill(weights_cteq[w] / pdfnorm, 1);
        //   // cout<<weights_cteq[w] / pdfnorm<<endl;
        // }
      }
    }
    if (doMG) weight_fsr = 1;
    h_fsr_pileup[0] ->Fill(bin_true, weight * weight_fsr);
    BinM_fsr_pileup[0]->Fill(bin, bin_true, weight * weight_fsr);
    if (weightpu_0 != 0) {
      h_fsr_pileup[1] ->Fill(bin_true, weight * weightpu_p / weightpu_0);
      BinM_fsr_pileup[1]->Fill(bin, bin_true, weight * weightpu_p / weightpu_0);
      h_fsr_pileup[2] ->Fill(bin_true, weight * weightpu_m / weightpu_0);
      BinM_fsr_pileup[2]->Fill(bin, bin_true, weight * weightpu_m / weightpu_0);
    }
    int idx = 1;
    //cout<<"Fillreco start"<<endl;
    FillRecoEffFluc(idx, h_ToyEffSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, bin_true);
    if (E0Tight) {
      FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, bin_true, 1);
      FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, bin_true, 2);
    } else {
      FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, bin, E1_pt, E1_eta, E0_pt, E0_eta, 1, BinM_eff, bin_true, 1);
      FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, bin, E1_pt, E1_eta, E0_pt, E0_eta, 1, BinM_eff, bin_true, 2);
    }
    FillTrigEffFluc(idx, h_ToyEffTMC, h_ToyEffTData, h_eff, weight, bin, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, bin_true);
    //cout<<"Fillreco start"<<endl;
  }
  cout << "done reading data for " << File_Signal_reco << " " << reco_name << endl;
  GetToyResponse(BinM_mcstat, BinMigration);
  cout << "done reading data for " << File_Signal_reco << " " << reco_name << endl;
  delete BinMigration;
  // delete l_phistar;
  // delete l_e0_pt;
  // delete l_e0_eta;
  // delete l_e1_pt;
  // delete l_e1_eta;
  // delete l_e0_tight;
  // delete l_ZY;
  // delete l_YZTruth;
  // delete l_phistar_true;
  // delete b_reco;
  // delete b_truth;
  delete t;

  return;
}

void NTupleZShape() {
  double Lumi = 19712.;
  double ttbar_weight = 23.64 / 4246440.;
  double tautau_weight = 1966.7 / 47271600.;
  double tbarw_weight = 11.1 / 493460.;
  double tw_weight = 11.1 / 497658.;
  double ww_weight = 54.84 / 10000430.;
  double wz_weight = 33.21 / 10000280.;
  double zz_weight = 17.0 / 9799908.;
  double signal_weight = 3531.89 / 30459500.; //
  if (!doMG) signal_weight = 1966.7 / 42705454.; //
 
  TH1D* Data = GetData(1. / Lumi);
  TH1D* Data_down = GetDataPt(1. / Lumi, -1);
  TH1D* Data_up = GetDataPt(1. / Lumi, 1);

  vector<TH2D*> h_ToyEffSF = GetEffSFToys();
  vector<TH2D*> h_ToyEffMSF = GetEffSFToys(1);
  vector<TH2D*> h_ToyEffTSF = GetEffSFToys(2);
  vector<TH2D*> h_ToyEffTMC = GetEffTMCToys(1);
  vector<TH2D*> h_ToyEffTData = GetEffTMCToys(0);
  vector<TH1D*> h_tt_eff, h_tt_fsr_pileup;
  vector<TH1D*> h_tautau_eff, h_tautau_fsr_pileup;
  vector<TH1D*> h_tbarw_eff, h_tbarw_fsr_pileup;
  vector<TH1D*> h_tw_eff, h_tw_fsr_pileup;
  vector<TH1D*> h_ww_eff, h_ww_fsr_pileup;
  vector<TH1D*> h_wz_eff, h_wz_fsr_pileup;
  vector<TH1D*> h_zz_eff, h_zz_fsr_pileup;
  GetBG(File_tt, ttbar_weight, h_tt_eff, h_tt_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_tautau, tautau_weight, h_tautau_eff, h_tautau_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_tbarw, tbarw_weight, h_tbarw_eff, h_tbarw_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_tw, tw_weight, h_tw_eff, h_tw_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_ww, ww_weight, h_ww_eff, h_ww_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_wz, wz_weight, h_wz_eff, h_wz_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
  GetBG(File_zz, zz_weight, h_zz_eff, h_zz_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);

  vector<TH1D*> bg_sf = GetToyBg();
  TH1D* bg_ss = GetBgSS();
  PrintBG(Data, h_tt_eff[0], h_tautau_eff[0], h_tbarw_eff[0], h_tw_eff[0], h_ww_eff[0], h_wz_eff[0], h_zz_eff[0], bg_ss);

  vector<TH1D*> h_data_eff = RemoveBG(Data, bg_sf[0], bg_ss, h_tt_eff, h_tautau_eff, h_tbarw_eff, h_tw_eff, h_ww_eff, h_wz_eff, h_zz_eff);
  cout<<"check eff"<<endl;
  NoNeg(h_data_eff);
  cout<<"check again"<<endl;
  NoNeg(h_data_eff);
  vector<TH1D*> h_data_fsr_pileup = RemoveBG(Data, bg_sf[0], bg_ss, h_tt_fsr_pileup, h_tautau_fsr_pileup, h_tbarw_fsr_pileup, h_tw_fsr_pileup, h_ww_fsr_pileup, h_wz_fsr_pileup, h_zz_fsr_pileup);
  cout<<"check fsr-pileup"<<endl;
  NoNeg(h_data_fsr_pileup);
  vector<TH1D*> h_data_bg;
  for (size_t idx = 0; idx < Ntoys + 7; idx++) {
    //TH1D* bgtemp = (TH1D*) h_tt_eff[0]->Clone();
    //bgtemp->Add(h_tautau_eff[0], 1.0);
      TH1D* bgtemp = (TH1D*) h_tautau_eff[0]->Clone();
    bgtemp->Add(h_tbarw_eff[0], 1.0);
    bgtemp->Add(h_tw_eff[0], 1.0);
    bgtemp->Add(h_ww_eff[0], 1.0);
    if (idx == Ntoys + 1 || idx == Ntoys + 2 || idx == Ntoys + 3 || idx == Ntoys + 4|| idx == Ntoys + 5|| idx == Ntoys + 6) {
      bgtemp->Multiply(bg_sf[0]);
    }
    else {
      bgtemp->Multiply(bg_sf[idx]);
    }
    TH1D* datatemp = (TH1D*) Data->Clone();
    datatemp->Add(bgtemp, -1.0);
    double scale = -1;
    double scaless = -2;
    double scalett=-1;
    if (idx == Ntoys + 1) {
      scale = -0.8;
    }
    if (idx == Ntoys + 2) {
      scale = -1.2;
    }
    if (idx == Ntoys + 3) {
      scaless = 0.0;
    }
    if (idx == Ntoys + 4) {
      scaless = -2.0;
    }
    if (idx == Ntoys + 5) {
      scalett = -.9;
    }
    if (idx == Ntoys + 6) {
      scalett = -1.1;
    }
    datatemp->Add(h_wz_eff[0], scale);
    datatemp->Add(h_zz_eff[0], scale);
    datatemp->Add(bg_ss, scaless);
    datatemp->Add(h_tt_eff[0], scalett);
    h_data_bg.push_back(datatemp);
  }
  cout<<"check bg"<<endl;
  NoNeg(h_data_bg);
  vector<TH1D*> h_data_pt;
  TH1D* bgtemp = (TH1D*) h_tt_eff[0]->Clone();
  bgtemp->Add(h_tautau_eff[0], 1.0);
  bgtemp->Add(h_tbarw_eff[0], 1.0);
  bgtemp->Add(h_tw_eff[0], 1.0);
  bgtemp->Add(h_ww_eff[0], 1.0);
  bgtemp->Multiply(bg_sf[0]);
  bgtemp->Add(h_wz_eff[0], 1.0);
  bgtemp->Add(h_zz_eff[0], 1.0);
  bgtemp->Add(bg_ss, 2.0);
  Data_up->Add(bgtemp, -1.0);
  Data_down->Add(bgtemp, -1.0);
  h_data_pt.push_back(Data_up);
  h_data_pt.push_back(Data_down);
  cout<<"check pt"<<endl;
  NoNeg(h_data_pt);

  delete Data;

  DeleteVec(h_tt_eff);
  DeleteVec(h_tt_fsr_pileup);
  DeleteVec(h_tautau_eff);
  DeleteVec(h_tautau_fsr_pileup);
  DeleteVec(h_tbarw_eff);
  DeleteVec(h_tbarw_fsr_pileup);
  DeleteVec(h_tw_eff);
  DeleteVec(h_tw_fsr_pileup);
  DeleteVec(h_ww_eff);
  DeleteVec(h_ww_fsr_pileup);
  DeleteVec(h_wz_eff);
  DeleteVec(h_wz_fsr_pileup);
  DeleteVec(h_zz_eff);
  DeleteVec(h_zz_fsr_pileup);

  vector<RooUnfoldResponse*> BinM_eff, BinM_mcstat, BinM_cteq, BinM_fsr_pileup;
  vector<TH1D*> mc_truereco_eff, mc_truereco_cteq, mc_truereco_fsr_pileup;
  GetBinM(signal_weight, BinM_eff, BinM_mcstat, BinM_cteq, BinM_fsr_pileup, mc_truereco_eff, mc_truereco_cteq, mc_truereco_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);

  TH1D* mc_truegen;
  vector<TH1D*> mc_truegen_cteq, mc_truegen_fsr_pileup;
  GetGen(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_fsr_pileup);

  vector<TH1D*> h_eff_eff = GetEff(mc_truereco_eff, mc_truegen);
  vector<TH1D*> h_eff_cteq = GetEff(mc_truereco_cteq, mc_truegen_cteq);

  vector<TH1D*> h_eff_fsr_pileup = GetEff(mc_truereco_fsr_pileup, mc_truegen_fsr_pileup);
  delete mc_truegen;
  DeleteVec(mc_truegen_cteq);
  DeleteVec(mc_truegen_fsr_pileup);
  DeleteVec(mc_truereco_eff);
  DeleteVec(mc_truereco_cteq);
  DeleteVec(mc_truereco_fsr_pileup);
  TH1D* eff_0 = (TH1D*) h_eff_eff[0]->Clone();
  for (uint idx = 0; idx < nbins; idx++) {
    eff_0->SetBinError(idx + 1, 0);
  }

  TMatrixD CovM_stat(nbins,nbins);
  vector<TGraphAsymmErrors *> g_data_phistar_unf = GetUnfoldedData(BinM_eff[0], h_data_eff[0], eff_0, CovM_stat);
  vector<TGraphAsymmErrors *> g_data_phistar_eff = GetUnfoldedData(BinM_eff, h_data_eff, h_eff_eff);
  vector<TGraphAsymmErrors *> g_data_phistar_bg = GetUnfoldedData(BinM_eff, h_data_bg, h_eff_eff, 1);
  vector<TGraphAsymmErrors *> g_data_phistar_pt = GetUnfoldedData(BinM_eff, h_data_pt, h_eff_eff, 1);
  vector<TGraphAsymmErrors *> g_data_phistar_cteq = GetUnfoldedData(BinM_cteq, h_data_eff[0], h_eff_cteq);
  vector<TGraphAsymmErrors *> g_data_phistar_fsr_pileup = GetUnfoldedData(BinM_fsr_pileup, h_data_fsr_pileup, h_eff_fsr_pileup);
  vector<TGraphAsymmErrors *> g_data_phistar_mcstat = GetUnfoldedData(BinM_mcstat, h_data_eff, h_eff_eff, 2);

  DeleteVec(BinM_eff);
  DeleteVec(BinM_cteq);
  DeleteVec(BinM_fsr_pileup);
  DeleteVec(BinM_mcstat);

  //copy graphs to make seperate absolute and normalised distributions
  vector<TGraphAsymmErrors *> g_data_norm_unf = CreateCopy(g_data_phistar_unf);
  vector<TGraphAsymmErrors *> g_data_norm_eff = CreateCopy(g_data_phistar_eff);
  vector<TGraphAsymmErrors *> g_data_norm_bg = CreateCopy(g_data_phistar_bg);
  vector<TGraphAsymmErrors *> g_data_norm_pt = CreateCopy(g_data_phistar_pt);
  vector<TGraphAsymmErrors *> g_data_norm_cteq = CreateCopy(g_data_phistar_cteq);
  vector<TGraphAsymmErrors *> g_data_norm_fsr_pileup = CreateCopy(g_data_phistar_fsr_pileup);
  vector<TGraphAsymmErrors *> g_data_norm_mcstat = CreateCopy(g_data_phistar_mcstat);
  //first get absolute distributions
  cout << "done unfolding, going to normalise" << endl;
  TMatrixD CovM_Abs_stat=*((TMatrixD*) (&CovM_stat)->Clone());
  NormalizeGraph(g_data_phistar_unf,0,CovM_Abs_stat); //empty
  NormalizeGraph(g_data_phistar_eff);
  NormalizeGraph(g_data_phistar_bg);
  NormalizeGraph(g_data_phistar_pt);
  NormalizeGraph(g_data_phistar_cteq);
  NormalizeGraph(g_data_phistar_fsr_pileup);
  NormalizeGraph(g_data_phistar_mcstat);

  TMatrixD CovM_Abs_unf=CalcCovM_unfsys(g_data_phistar_unf[0],0);//,1 for normalised
  TMatrixD CovM_Abs_mcstat_toy=CalcCovM_toy(g_data_phistar_mcstat);
  TMatrixD CovM_Abs_mcstat_event=CalcCovM_uncor(g_data_phistar_unf[0], g_data_phistar_eff[0]);
  TMatrixD CovM_Abs_mcstat=CovM_Abs_mcstat_toy+CovM_Abs_mcstat_event;
  TMatrixD CovM_Abs_reff=CalcCovM_toy(g_data_phistar_eff,1,1);
  TMatrixD CovM_Abs_meff=CalcCovM_toy(g_data_phistar_eff,2,1);
  TMatrixD CovM_Abs_teff=CalcCovM_toy(g_data_phistar_eff,3,1);
  TMatrixD CovM_Abs_teff_m=CalcCovM_toy(g_data_phistar_eff,4,1);
  TMatrixD CovM_Abs_teff_d=CalcCovM_toy(g_data_phistar_eff,5,1);
  TMatrixD CovM_Abs_eff=CovM_Abs_reff+CovM_Abs_meff+CovM_Abs_teff+CovM_Abs_teff_m+CovM_Abs_teff_d;
  TMatrixD CovM_Abs_bg_toy=CalcCovM_toy(g_data_phistar_bg);
  TMatrixD CovM_Abs_bg_ZZWZ=CalcCovM_downup(g_data_phistar_bg[Ntoys+1],g_data_phistar_bg[Ntoys+2],g_data_phistar_bg[0]);
  TMatrixD CovM_Abs_bg_QCD=CalcCovM_downup(g_data_phistar_bg[Ntoys+3],g_data_phistar_bg[Ntoys+4],g_data_phistar_bg[0]);
  TMatrixD CovM_Abs_bg_TTBar=CalcCovM_downup(g_data_phistar_bg[Ntoys+5],g_data_phistar_bg[Ntoys+6],g_data_phistar_bg[0]);
  TMatrixD CovM_Abs_bg=CovM_Abs_bg_toy+CovM_Abs_bg_ZZWZ+CovM_Abs_bg_QCD+CovM_Abs_bg_TTBar;
  TMatrixD CovM_Abs_pileup=CalcCovM_downup(g_data_phistar_fsr_pileup[1], g_data_phistar_fsr_pileup[2], g_data_phistar_eff[0]);
  TMatrixD CovM_Abs_enscale=CalcCovM_downup(g_data_phistar_pt[0], g_data_phistar_pt[1], g_data_phistar_eff[0]);
  TMatrixD CovM_Abs_fsr=CalcCovM_down(g_data_phistar_fsr_pileup[0], g_data_phistar_eff[0]);
  TMatrixD CovM_Abs_cteq=CalcCovM_cteq(g_data_phistar_cteq);
  TMatrixD CovM_Abs_lumi=CalcCovM_lumi(g_data_phistar_unf[0]);
  TMatrixD CovM_Abs_syst=CovM_Abs_unf+CovM_Abs_mcstat+CovM_Abs_eff+CovM_Abs_bg+CovM_Abs_pileup+CovM_Abs_enscale+CovM_Abs_fsr+CovM_Abs_lumi;
  if (!doMG) CovM_Abs_syst+=CovM_Abs_cteq;
  TMatrixD CovM_Abs_tot=CovM_Abs_syst+CovM_Abs_stat;

  TGraphAsymmErrors* g_syst_phistar_mctoy = CalcTotalSysU_toymc(g_data_phistar_mcstat, g_data_phistar_eff[0]);
  TGraphAsymmErrors* g_syst_phistar_bg = CalcTotalSysU_toymc(g_data_phistar_bg, g_data_phistar_bg[0], 1);
  TGraphAsymmErrors* g_syst_phistar_reff = CalcTotalSysU_toymc(g_data_phistar_eff, g_data_phistar_eff[0], 0, 1, 1);
  TGraphAsymmErrors* g_syst_phistar_meff = CalcTotalSysU_toymc(g_data_phistar_eff, g_data_phistar_eff[0], 0, 2, 1);
  TGraphAsymmErrors* g_syst_phistar_teff = CalcTotalSysU_toymc(g_data_phistar_eff, g_data_phistar_eff[0], 0, 3, 1);
  TGraphAsymmErrors* g_syst_phistar_teff_m = CalcTotalSysU_toymc(g_data_phistar_eff, g_data_phistar_eff[0], 0, 4, 1);
  TGraphAsymmErrors* g_syst_phistar_teff_d = CalcTotalSysU_toymc(g_data_phistar_eff, g_data_phistar_eff[0], 0, 5, 1);
  TGraphAsymmErrors* g_syst_phistar_eff = CalcTotalSysU_comb5(g_syst_phistar_reff, g_syst_phistar_meff, g_syst_phistar_teff, g_syst_phistar_teff_m, g_syst_phistar_teff_d);
  TGraphAsymmErrors* g_syst_phistar_mcstat = CalcTotalSysU_comb3(g_data_phistar_eff[0], g_data_phistar_unf[0], g_syst_phistar_mctoy, 0);
  TGraphAsymmErrors* g_syst_phistar_cteq = CalcTotalSysU_updown(g_data_phistar_cteq, g_data_phistar_cteq[0], 1, 1);
  TGraphAsymmErrors* g_syst_phistar_fsr = CalcTotalSysU_fsr(g_data_phistar_fsr_pileup[0], g_data_phistar_eff[0]); //commented out since it was not used any more
  TGraphAsymmErrors* g_syst_phistar_pu_Temp = CalcTotalSysU_pileup(g_data_phistar_fsr_pileup[1], g_data_phistar_fsr_pileup[2], g_data_phistar_eff[0]);
  TGraphAsymmErrors* g_syst_phistar_pt_Temp = CalcTotalSysU_pileup(g_data_phistar_pt[0], g_data_phistar_pt[1], g_data_phistar_eff[0]);
  //so take the Pileup and and Pt
  ///ADD IT HERE
  TGraphAsymmErrors* g_syst_phistar_pu = ErrorSmoother(g_syst_phistar_pu_Temp);
  TGraphAsymmErrors* g_syst_phistar_pt = ErrorSmoother( g_syst_phistar_pt_Temp);
  TGraphAsymmErrors* g_syst_phistar_unsyst = GetUnfoldingSysUn(0, g_data_phistar_unf[0]);
 
  // DeleteVec(g_data_phistar_eff);
  // DeleteVec(g_data_phistar_bg);
  // DeleteVec(g_data_phistar_pt);
  // DeleteVec(g_data_phistar_cteq);
  // DeleteVec(g_data_phistar_mcstat);

  vector<TGraphAsymmErrors *> g_data_syst_muon;
  g_data_syst_muon.push_back(g_data_phistar_unf[0]);
  g_data_syst_muon.push_back(g_syst_phistar_eff);
  g_data_syst_muon.push_back(g_syst_phistar_mcstat);
  g_data_syst_muon.push_back(g_syst_phistar_pt);
  vector<std::string> syst_list_muon;
  syst_list_muon.push_back("unfolding");
  syst_list_muon.push_back("eff");
  syst_list_muon.push_back("mcstat");
  syst_list_muon.push_back("pt");

  TGraphAsymmErrors* g_data_final_muon = GetDataFinal(g_data_syst_muon, syst_list_muon, 0, 0);
  // DeleteVec(g_data_syst_muon);
  TH1D* h_data_elec = ConvertToHist(g_data_final_muon, "h_data_elec");
  TH1D* h_data_bgnd = ConvertToHist(g_syst_phistar_bg, "h_data_bgnd");
  TH1D* h_data_pdf;
  if (!doMG) h_data_pdf = ConvertToHist(g_syst_phistar_cteq, "h_data_pdf");
  TH1D* h_data_fsr = ConvertToHist(g_data_phistar_fsr_pileup[0], "h_data_fsr");
  TH1D* h_data_pup = ConvertToHist(g_data_phistar_fsr_pileup[1], "h_data_pup");
  TH1D* h_data_pum = ConvertToHist(g_data_phistar_fsr_pileup[2], "h_data_pum");
  TH1D* h_data_uns = ConvertToHist(g_syst_phistar_unsyst, "h_data_uns");

  vector<TGraphAsymmErrors *> g_data_syst;
  g_data_syst.push_back(g_data_phistar_unf[0]);
  g_data_syst.push_back(g_syst_phistar_eff);
  g_data_syst.push_back(g_syst_phistar_mcstat);
  g_data_syst.push_back(g_syst_phistar_bg);
  g_data_syst.push_back(g_syst_phistar_fsr);
  g_data_syst.push_back(g_syst_phistar_pu);
  g_data_syst.push_back(g_syst_phistar_pt);
  g_data_syst.push_back(g_syst_phistar_unsyst);
  if (!doMG) g_data_syst.push_back(g_syst_phistar_cteq);
  vector<TMatrixD> m_data_syst;
  m_data_syst.push_back(CovM_Abs_stat);
  m_data_syst.push_back(CovM_Abs_eff);
  m_data_syst.push_back(CovM_Abs_mcstat);
  m_data_syst.push_back(CovM_Abs_bg);
  m_data_syst.push_back(CovM_Abs_fsr);
  m_data_syst.push_back(CovM_Abs_pileup);
  m_data_syst.push_back(CovM_Abs_enscale);
  m_data_syst.push_back(CovM_Abs_unf);
  if (!doMG) m_data_syst.push_back(CovM_Abs_cteq);
  vector<std::string> syst_list;
  syst_list.push_back("unfolding");
  syst_list.push_back("eff");
  syst_list.push_back("mcstat");
  syst_list.push_back("bg");
  syst_list.push_back("fsr");
  syst_list.push_back("pileup");
  syst_list.push_back("pt");
  syst_list.push_back("unsyst");
  if (!doMG) syst_list.push_back("cteq");

  //So currently I am just sticking in the part that makes the covariant matrix here, right after everything is normalized.
  // cout << " test 1" << endl;
  TGraphAsymmErrors* g_data_final = GetDataFinal(g_data_syst, syst_list);
  TGraphAsymmErrors* g_data_final_m = GetDataFinal(m_data_syst, syst_list, g_syst_phistar_eff);
  // DeleteVec(g_data_syst);
  std::string textn = "Output/Final_Hist_";
  textn += Tag;
  textn += "Abs_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr(textn.c_str(), "RECREATE");
  h_data_elec->Write();
  h_data_bgnd->Write();
  if (!doMG)h_data_pdf->Write();
  if (!doMG)h_data_fsr->Write();
  h_data_pup->Write();
  h_data_pum->Write();
  h_data_uns->Write();

  textn = "Output/Data_Graph_";
  textn += Tag;
  textn += "Abs_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr2(textn.c_str(), "RECREATE");
  g_data_final->Write();

  textn = "Output/Data_Graph_MCStat_";
  textn += Tag;
  textn += "Abs_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr2b(textn.c_str(), "RECREATE");
  g_syst_phistar_mcstat->Write();

  textn = "Output/CovM_";
  textn += Tag;
  textn += "Abs_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr3(textn.c_str(), "RECREATE");
  CovM_Abs_stat.Write("CovM_Abs_stat");
  CovM_Abs_meff.Write("CovM_Abs_meff");
  CovM_Abs_reff.Write("CovM_Abs_reff");
  CovM_Abs_teff.Write("CovM_Abs_teff");
  CovM_Abs_teff_m.Write("CovM_Abs_teff_m");
  CovM_Abs_teff_d.Write("CovM_Abs_teff_d");
  CovM_Abs_eff.Write("CovM_Abs_eff");
  CovM_Abs_mcstat_toy.Write("CovM_Abs_mcstat_toy");
  CovM_Abs_mcstat_event.Write("CovM_Abs_mcstat_event");
  CovM_Abs_mcstat.Write("CovM_Abs_mcstat");
  CovM_Abs_bg_QCD.Write("CovM_Abs_bg_QCD");
  CovM_Abs_bg_TTBar.Write("CovM_Abs_bg_TTBar");
  CovM_Abs_bg_ZZWZ.Write("CovM_Abs_bg_ZZWZ");
  CovM_Abs_bg_toy.Write("CovM_Abs_bg_toy");
  CovM_Abs_bg.Write("CovM_Abs_bg");
  CovM_Abs_fsr.Write("CovM_Abs_fsr");
  CovM_Abs_pileup.Write("CovM_Abs_pu");
  CovM_Abs_enscale.Write("CovM_Abs_pt");
  CovM_Abs_unf.Write("CovM_Abs_unf");
  CovM_Abs_lumi.Write("CovM_Abs_lumi");
  if (!doMG) CovM_Abs_cteq.Write("CovM_Abs_pdf");
  CovM_Abs_syst.Write("CovM_Abs_syst");
  CovM_Abs_tot.Write("CovM_Abs_tot");
 
  TMatrixD CovM_Norm_stat=*((TMatrixD*) (&CovM_stat)->Clone());
  NormalizeGraph(g_data_norm_unf, 1, CovM_Norm_stat);
  NormalizeGraph(g_data_norm_eff, 1);
  NormalizeGraph(g_data_norm_bg, 1);
  NormalizeGraph(g_data_norm_pt, 1);
  if (!doMG) NormalizeGraph(g_data_norm_cteq, 1);
  NormalizeGraph(g_data_norm_fsr_pileup, 1);
  NormalizeGraph(g_data_norm_mcstat, 1);

  TMatrixD CovM_Norm_unf=CalcCovM_unfsys(g_data_norm_unf[0],1);//,1 for normalised
  TMatrixD CovM_Norm_mcstat_toy=CalcCovM_toy(g_data_norm_mcstat);
  TMatrixD CovM_Norm_mcstat_event=CalcCovM_uncor(g_data_norm_unf[0], g_data_norm_eff[0]);
  TMatrixD CovM_Norm_mcstat=CovM_Norm_mcstat_toy+CovM_Norm_mcstat_event;
  TMatrixD CovM_Norm_reff=CalcCovM_toy(g_data_norm_eff,1,1);
  TMatrixD CovM_Norm_meff=CalcCovM_toy(g_data_norm_eff,2,1);
  TMatrixD CovM_Norm_teff=CalcCovM_toy(g_data_norm_eff,3,1);
  TMatrixD CovM_Norm_teff_m=CalcCovM_toy(g_data_norm_eff,4,1);
  TMatrixD CovM_Norm_teff_d=CalcCovM_toy(g_data_norm_eff,5,1);
  TMatrixD CovM_Norm_eff=CovM_Norm_reff+CovM_Norm_meff+CovM_Norm_teff+CovM_Norm_teff_m+CovM_Norm_teff_d;
  TMatrixD CovM_Norm_bg_toy=CalcCovM_toy(g_data_norm_bg);
  TMatrixD CovM_Norm_bg_ZZWZ=CalcCovM_downup(g_data_norm_bg[Ntoys+1],g_data_norm_bg[Ntoys+2],g_data_norm_bg[0]);
  TMatrixD CovM_Norm_bg_QCD=CalcCovM_downup(g_data_norm_bg[Ntoys+3],g_data_norm_bg[Ntoys+4],g_data_norm_bg[0]);
  TMatrixD CovM_Norm_bg_TTBar=CalcCovM_downup(g_data_norm_bg[Ntoys+5],g_data_norm_bg[Ntoys+6],g_data_norm_bg[0]);
  TMatrixD CovM_Norm_bg=CovM_Norm_bg_toy+CovM_Norm_bg_ZZWZ+CovM_Norm_bg_QCD+CovM_Norm_bg_TTBar;
  TMatrixD CovM_Norm_pileup=CalcCovM_downup(g_data_norm_fsr_pileup[1], g_data_norm_fsr_pileup[2], g_data_norm_eff[0]);
  TMatrixD CovM_Norm_enscale=CalcCovM_downup(g_data_norm_pt[0], g_data_norm_pt[1], g_data_norm_eff[0]);
  TMatrixD CovM_Norm_fsr=CalcCovM_down(g_data_norm_fsr_pileup[0], g_data_norm_eff[0]);
  TMatrixD CovM_Norm_cteq=CalcCovM_cteq(g_data_norm_cteq);
  TMatrixD CovM_Norm_syst=CovM_Norm_unf+CovM_Norm_mcstat+CovM_Norm_eff+CovM_Norm_bg+CovM_Norm_pileup+CovM_Norm_enscale+CovM_Norm_fsr;
  if (!doMG) CovM_Norm_syst+=CovM_Norm_cteq;
  TMatrixD CovM_Norm_tot=CovM_Norm_syst+CovM_Norm_stat;

  TGraphAsymmErrors* g_syst_norm_mctoy = CalcTotalSysU_toymc(g_data_norm_mcstat, g_data_norm_eff[0]);
  TGraphAsymmErrors* g_syst_norm_bg = CalcTotalSysU_toymc(g_data_norm_bg, g_data_norm_bg[0], 1);
  TGraphAsymmErrors* g_syst_norm_reff = CalcTotalSysU_toymc(g_data_norm_eff, g_data_norm_eff[0], 0, 1, 1);
  TGraphAsymmErrors* g_syst_norm_meff = CalcTotalSysU_toymc(g_data_norm_eff, g_data_norm_eff[0], 0, 2, 1);
  TGraphAsymmErrors* g_syst_norm_teff = CalcTotalSysU_toymc(g_data_norm_eff, g_data_norm_eff[0], 0, 3, 1);
  TGraphAsymmErrors* g_syst_norm_teff_m = CalcTotalSysU_toymc(g_data_norm_eff, g_data_norm_eff[0], 0, 4, 1);
  TGraphAsymmErrors* g_syst_norm_teff_d = CalcTotalSysU_toymc(g_data_norm_eff, g_data_norm_eff[0], 0, 5, 1);
  TGraphAsymmErrors* g_syst_norm_eff = CalcTotalSysU_comb5(g_syst_norm_reff, g_syst_norm_meff, g_syst_norm_teff, g_syst_norm_teff_m, g_syst_norm_teff_d);
  TGraphAsymmErrors* g_syst_norm_mcstat = CalcTotalSysU_comb3(g_data_norm_eff[0], g_data_norm_unf[0], g_syst_norm_mctoy, 0);
  TGraphAsymmErrors* g_syst_norm_cteq = CalcTotalSysU_updown(g_data_norm_cteq, g_data_norm_cteq[0], 1, 1);
  TGraphAsymmErrors* g_syst_norm_fsr = CalcTotalSysU_fsr(g_data_norm_fsr_pileup[0], g_data_norm_eff[0]); //not used so I commented it out
  TGraphAsymmErrors* g_syst_norm_pu_temp = CalcTotalSysU_pileup(g_data_norm_fsr_pileup[1], g_data_norm_fsr_pileup[2], g_data_norm_eff[0]);
  TGraphAsymmErrors* g_syst_norm_pt_temp = CalcTotalSysU_pileup(g_data_norm_pt[0], g_data_norm_pt[1], g_data_norm_eff[0]);
  TGraphAsymmErrors* g_syst_norm_unsyst = GetUnfoldingSysUn(1, g_data_norm_unf[0]);
  // and here is where I need to do everything AGAIN
  TGraphAsymmErrors* g_syst_norm_pu = ErrorSmoother(g_syst_norm_pu_temp);
  TGraphAsymmErrors* g_syst_norm_pt = ErrorSmoother(g_syst_norm_pt_temp);

  vector<TGraphAsymmErrors *> g_data_syst_norm_muon;
  g_data_syst_norm_muon.push_back(g_data_norm_unf[0]);
  g_data_syst_norm_muon.push_back(g_syst_norm_eff);
  g_data_syst_norm_muon.push_back(g_syst_norm_mcstat);
  g_data_syst_norm_muon.push_back(g_syst_norm_pt);
  vector<std::string> syst_list_norm_muon;
  syst_list_norm_muon.push_back("unfolding");
  syst_list_norm_muon.push_back("eff");
  syst_list_norm_muon.push_back("mcstat");
  syst_list_norm_muon.push_back("pt");
  TGraphAsymmErrors* g_data_final_norm_muon = GetDataFinal(g_data_syst_norm_muon, syst_list_norm_muon, 1, 0);

  TH1D* h_data_norm_elec = ConvertToHist(g_data_final_norm_muon, "h_data_elec");
  TH1D* h_data_norm_bgnd = ConvertToHist(g_syst_norm_bg, "h_data_bgnd");
  TH1D* h_data_norm_pdf;
  if (!doMG) h_data_norm_pdf = ConvertToHist(g_syst_norm_cteq, "h_data_pdf");
  TH1D * h_data_norm_fsr = ConvertToHist(g_data_norm_fsr_pileup[0], "h_data_fsr");
  TH1D * h_data_norm_pup = ConvertToHist(g_data_norm_fsr_pileup[1], "h_data_pup");
  TH1D * h_data_norm_pum = ConvertToHist(g_data_norm_fsr_pileup[2], "h_data_pum");
  TH1D * h_data_norm_uns = ConvertToHist(g_syst_norm_unsyst, "h_data_uns");

  vector<TGraphAsymmErrors *> g_data_syst_norm;
  g_data_syst_norm.push_back(g_data_norm_unf[0]);
  g_data_syst_norm.push_back(g_syst_norm_eff);
  g_data_syst_norm.push_back(g_syst_norm_mcstat);
  g_data_syst_norm.push_back(g_syst_norm_bg);
  if (!doMG)g_data_syst_norm.push_back(g_syst_norm_fsr);
  g_data_syst_norm.push_back(g_syst_norm_pu);
  g_data_syst_norm.push_back(g_syst_norm_pt);
  g_data_syst_norm.push_back(g_syst_norm_unsyst);
  if (!doMG) g_data_syst_norm.push_back(g_syst_norm_cteq);
  vector<TMatrixD> m_data_syst_norm;
  m_data_syst_norm.push_back(CovM_Norm_stat);
  m_data_syst_norm.push_back(CovM_Norm_eff);
  m_data_syst_norm.push_back(CovM_Norm_mcstat);
  m_data_syst_norm.push_back(CovM_Norm_bg);
  m_data_syst_norm.push_back(CovM_Norm_fsr);
  m_data_syst_norm.push_back(CovM_Norm_pileup);
  m_data_syst_norm.push_back(CovM_Norm_enscale);
  m_data_syst_norm.push_back(CovM_Norm_unf);
  if (!doMG) m_data_syst_norm.push_back(CovM_Norm_cteq);
  vector<std::string> syst_list_norm;
  syst_list_norm.push_back("unfolding");
  syst_list_norm.push_back("eff");
  syst_list_norm.push_back("mcstat");
  syst_list_norm.push_back("bg");
  if (!doMG)syst_list_norm.push_back("fsr");
  syst_list_norm.push_back("pileup");
  syst_list_norm.push_back("pt");
  syst_list_norm.push_back("unsyst");
  if (!doMG) syst_list_norm.push_back("cteq");

  TGraphAsymmErrors * g_data_final_norm = GetDataFinal(g_data_syst_norm, syst_list_norm, 1);
  TGraphAsymmErrors * g_data_final_norm_m = GetDataFinal(m_data_syst_norm, syst_list_norm, g_syst_norm_eff, 1);

  textn = "Output/Final_Hist_";
  textn += Tag;
  textn += "Norm_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr4(textn.c_str(), "RECREATE");
  h_data_norm_elec->Write();
  h_data_norm_bgnd->Write();
  if (!doMG)h_data_norm_pdf->Write();
  if (!doMG)h_data_norm_fsr->Write();
  h_data_norm_pup->Write();
  h_data_norm_pum->Write();
  h_data_norm_uns->Write();

  textn = "Output/Data_Graph_";
  textn += Tag;
  textn += "Norm_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr5(textn.c_str(), "RECREATE");
  cout << "file name " << textn << endl;
  g_data_final_norm->Write();

  textn = "Output/Data_Graph_MCStat_";
  textn += Tag;
  textn += "Norm_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr5b(textn.c_str(), "RECREATE");
  cout << "file name " << textn << endl;
  g_syst_norm_mcstat->Write();

  textn = "Output/CovM_";
  textn += Tag;
  textn += "Norm_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "BornUnChanged.root";
  if (elec == 2)textn += "Naked.root";
  TFile tr6(textn.c_str(), "RECREATE");
  CovM_Norm_stat.Write("CovM_Norm_stat");
  CovM_Norm_meff.Write("CovM_Norm_meff");
  CovM_Norm_reff.Write("CovM_Norm_reff");
  CovM_Norm_teff.Write("CovM_Norm_teff");
  CovM_Norm_teff_m.Write("CovM_Norm_teff_m");
  CovM_Norm_teff_d.Write("CovM_Norm_teff_d");
  CovM_Norm_eff.Write("CovM_Norm_eff");
  CovM_Norm_mcstat_toy.Write("CovM_Norm_mcstat_toy");
  CovM_Norm_mcstat_event.Write("CovM_Norm_mcstat_event");
  CovM_Norm_mcstat.Write("CovM_Norm_mcstat");
  CovM_Norm_bg_QCD.Write("CovM_Norm_bg_QCD");
  CovM_Norm_bg_ZZWZ.Write("CovM_Norm_bg_ZZWZ");
  CovM_Norm_bg_toy.Write("CovM_Norm_bg_toy");
  CovM_Norm_bg_TTBar.Write("CovM_Norm_bg_TTBar");
  CovM_Norm_bg.Write("CovM_Norm_bg");
  CovM_Norm_fsr.Write("CovM_Norm_fsr");
  CovM_Norm_pileup.Write("CovM_Norm_pu");
  CovM_Norm_enscale.Write("CovM_Norm_pt");
  CovM_Norm_unf.Write("CovM_Norm_unf");
  if (!doMG) CovM_Norm_cteq.Write("CovM_Norm_pdf"); 
  CovM_Norm_syst.Write("CovM_Norm_syst");
  CovM_Norm_tot.Write("CovM_Norm_tot");
}
