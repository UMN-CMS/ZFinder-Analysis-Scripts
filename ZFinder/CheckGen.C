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

const double TEfficiencyEtaBins[11] = {-2.1, -2.0, -1.556, -1.442, -0.8, 0.0, 0.8, 1.442, 1.556, 2.0, 2.1};
const int TEfficiencyETBins[5] = {30, 40, 50, 70, 250};
const double TEfficiencyData[10][4][2] = {{{0.741,  0.003}, {0.773,  0.003}, {0.780,  0.005}, {0.790,  0.010}},
					  {{0.734,  0.001}, {0.772,  0.001}, {0.786,  0.002}, {0.792,  0.005}},
					  {{0.725,  0.003}, {0.821,  0.002}, {0.809,  0.004}, {0.848,  0.010}},
					  {{0.893,  0.0005},{0.9396, 0.0003},{0.9509, 0.0006},{0.966,  0.001}},
					  {{0.9213, 0.0004},{0.9528, 0.0002},{0.9601, 0.0004},{0.969 , 0.001}},
					  {{0.9174, 0.0004},{0.9473, 0.0004},{0.9561, 0.0004},{0.963,  0.001}},
					  {{0.8964, 0.0005},{0.9424, 0.0003},{0.9533, 0.0006},{0.966,  0.001}},
					  {{0.714,  0.003}, {0.823,  0.002}, {0.827,  0.004}, {0.861,  0.010}},
					  {{0.758,  0.001}, {0.800,  0.001}, {0.811,  0.002}, {0.823,  0.005}},
					  {{0.764,  0.003}, {0.792,  0.002}, {0.797,  0.005}, {0.820,  0.010}}};
const double TEfficiencyMC[10][4][2]   = {{{0.734,  0.004}, {0.769,  0.003}, {0.771,  0.004}, {0.760,  0.020}},
					  {{0.736,  0.002}, {0.768,  0.004}, {0.779,  0.003}, {0.789,  0.008}},
					  {{0.791,  0.004}, {0.847,  0.003}, {0.850,  0.006}, {0.870,  0.020}},
					  {{0.9395, 0.0006},{0.9612, 0.0004},{0.9690, 0.0007},{0.980,  0.002}},
					  {{0.9469, 0.0005},{0.9670, 0.0003},{0.9745, 0.0005},{0.982,  0.001}},
					  {{0.9466, 0.0005},{0.9665, 0.0003},{0.9739, 0.0005},{0.982,  0.001}},
					  {{0.9364, 0.0007},{0.9597, 0.0004},{0.9668, 0.0008},{0.979,  0.002}},
					  {{0.779,  0.004}, {0.841,  0.003}, {0.842,  0.006}, {0.860,  0.020}},
					  {{0.749,  0.002}, {0.786,  0.002}, {0.798,  0.003}, {0.810,  0.008}},
					  {{0.737,  0.004}, {0.769,  0.004}, {0.779,  0.007}, {0.820,  0.020}}};

const double EfficiencyEtaBins[6] = {0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};
const int EfficiencyETBins[7] = {10, 15, 20, 30, 40, 50, 1000000};
const double EfficiencySF[5][6][3] = {{{0.977,0.024,0.054}, {0.997,0.009,0.031}, {0.982,0.003,0.012}, {0.988,0.001,0.008}, {0.990,0.001,0.004}, {0.990,0.001,0.004}},
                                      {{0.977,0.024,0.054}, {0.997,0.009,0.031}, {0.993,0.002,0.012}, {0.993,0.001,0.008}, {0.993,0.001,0.004}, {0.991,0.001,0.004}},
                                      {{1.076,0.152,0.095}, {0.952,0.025,0.070}, {1.016,0.012,0.020}, {0.985,0.004,0.009}, {0.987,0.004,0.004}, {0.974,0.009,0.006}},
                                      {{1.096,0.036,0.060}, {1.008,0.010,0.031}, {0.988,0.003,0.012}, {0.993,0.002,0.008}, {0.992,0.001,0.004}, {0.990,0.003,0.004}},
                                      {{1.096,0.036,0.060}, {1.008,0.010,0.031}, {1.002,0.004,0.012}, {1.004,0.002,0.008}, {1.005,0.002,0.004}, {0.998,0.004,0.004}}};

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
RooBinning phistarbin(nphistar,phistarBins);

vector<std::string> Get_File_Data(){
  vector<std::string> vec;
  vec.push_back("/afs/cern.ch/work/r/ruckstuh/Data.root");
  return vec;
}
// std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/MadGraph_reco.root";
// std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/MadGraph_gen.root";
// std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_cteq_reco.root";
// std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_cteq_gen.root";
std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco.root";
std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Powheg_reco= "/afs/cern.ch/work/r/ruckstuh/public/Powheg_reco.root";
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


void GetGenPhiStar(){
  gErrorIgnoreLevel = kError;
  cout<<"reading signal gen"<<endl;

  TChain* t = new TChain(gen_name.c_str(),gen_name.c_str());
  int nfiles;
  if (doMG) nfiles=t->Add(File_Signal_gen.c_str());
  else nfiles=t->Add(File_Powheg_gen.c_str());
  t->SetBranchStatus("reco",0); //to disable all branches
  TBranch *b_truth=t->GetBranch("truth");
  TBranch *b_event=t->GetBranch("event_info");
  // TLeaf *l_phistar=b_truth->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar=b_truth->GetLeaf("z_phistar_naked");
  TLeaf *l_m=b_truth->GetLeaf("z_m");
  TLeaf *l_e0_q =b_truth->GetLeaf("e_charge0");
  TLeaf *l_e1_q =b_truth->GetLeaf("e_charge1");
  TLeaf *l_e0_pt =b_truth->GetLeaf("e_pt0");
  TLeaf *l_e1_pt =b_truth->GetLeaf("e_pt1");
  TLeaf *l_e0_eta =b_truth->GetLeaf("e_eta0");
  TLeaf *l_e1_eta =b_truth->GetLeaf("e_eta1");
  TLeaf *l_number =b_event->GetLeaf("event_number");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
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
  t->SetBranchAddress("weight_cteq6",&weight_cteq6);

  int ntot=0;
  double npu=0;
  int ntot2=0;
  double npu2=0;

  double cte=0;
  double nnp=0;
  double mst=0;
  double ct6=0;

  // vector<int > eventnvec;
  // for (int i=0; i<t->GetEntries();i++){
  for (int i=0; i<10000;i++){
    if (i%100000==0) cout<<i<<endl;
    t->GetEntry(i);
    double weight=1;
    double phistar=l_phistar->GetValue();
    double mz=l_m->GetValue();
    int q0=l_e0_q->GetValue();
    int q1=l_e1_q->GetValue();
    double eta0=l_e0_eta->GetValue();
    double eta1=l_e1_eta->GetValue();
    double pt0=l_e0_pt->GetValue();
    double pt1=l_e1_pt->GetValue();
    int eventN=l_number->GetValue();

    // for (int j=0; j<eventnvec.size(); j++){
    //   if (eventN==eventnvec[j]) cout<<"Event number dublicated:"<<eventN<<endl;
    // }
    // eventnvec.push_back(l_number->GetValue());

    if (mz<60 || mz>120) cout<<"Mass:"<<mz<<" Electron0, pt:"<<pt0<<", eta:"<<eta0<<", q:"<<q0<<" Electron1, pt:"<<pt1<<", eta:"<<eta1<<", q:"<<q1<<endl;
    else if (!((q0==-1 && q1==1)||(q0==1 && q1==-1))) cout<<"Mass:"<<mz<<" Electron0, pt:"<<pt0<<", eta:"<<eta0<<", q:"<<q0<<" Electron1, pt:"<<pt1<<", eta:"<<eta1<<", q:"<<q1<<endl;
    else if (pt0 <20 ||pt1<20) cout<<"Mass:"<<mz<<" Electron0, pt:"<<pt0<<", eta:"<<eta0<<", q:"<<q0<<" Electron1, pt:"<<pt1<<", eta:"<<eta1<<", q:"<<q1<<endl;
    else if (fabs(eta0)>2.4 ||fabs(eta1)>2.4) cout<<"Mass:"<<mz<<" Electron0, pt:"<<pt0<<", eta:"<<eta0<<", q:"<<q0<<" Electron1, pt:"<<pt1<<", eta:"<<eta1<<", q:"<<q1<<endl;
    else if (!((pt0>30 && fabs(eta0)<2.1) || (pt1>30 && fabs(eta1)<2.1))) cout<<"Mass:"<<mz<<" Electron0, pt:"<<pt0<<", eta:"<<eta0<<", q:"<<q0<<" Electron1, pt:"<<pt1<<", eta:"<<eta1<<", q:"<<q1<<endl;

    // cout<<weights_cteq[0]<<endl;
    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2) {weight=weight*weights[w];}
    }
    ntot++;
    npu+=weight;
    mst+=weights_mstw[0];
    cte+=weights_cteq[0];
    nnp+=weights_nnpdf[0];
    ct6+=weight_cteq6;
    if (pt0>pt1){
      if (fabs(eta0)<2.1) {
	ntot2++; 
	npu2+=weight;
      }
    }
    else {
      if (fabs(eta1)<2.1) {
	ntot2++; 
	npu2+=weight;
      }
    }
  }
  cout<<ntot<<"  "<<npu<<" "<<npu/double(ntot)<<"  "<<double(ntot2)/double(ntot)<<"  "<<npu2/npu<<endl;
  cout<<ct6<<"  "<<cte<<" "<<mst<<"  "<<nnp<<"  "<<double(cte)/double(ct6)<<"  "<<double(mst)/double(ct6)<<"  "<<double(nnp)/double(ct6)<<endl;
  return;
}

