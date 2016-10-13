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

const int elec=1;
const int doMG=1;
const std::string Tag="";

const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 5, 10};
const double yBins[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
size_t ny=(sizeof(yBins)/sizeof(yBins[0]))-1;
size_t nbins=nphistar*ny;

std::string File_Signal_gen = "/afs/cern.ch/work/r/ruckstuh/public/MG_dressed_gen.root";
std::string File_Signal_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/MG_born_gen.root";
std::string File_Signal_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/MG_bare_gen.root";
std::string File_Powheg_gen = "/afs/cern.ch/work/r/ruckstuh/public/PH_dressed_gen.root";
std::string File_Powheg_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/PH_born_gen.root";
std::string File_Powheg_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/PH_bare_gen.root";

std::string gen_name ="Combined Gen Cuts Reco";

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

TGraphAsymmErrors * GetMCFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist, bool doNorm=0, bool doLumi = 1) {
  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nbins);
  uint nun = graph.size();
  for (size_t ibin = 0; ibin < nbins; ibin++) {
    double x, y, errorl, errorh;
    graph[0]->GetPoint(ibin, x, y);
    double x_nom = x;
    double y_nom = y;
    double error_sys_max = 0;
    double error_sys_min = 0;
    if (doMG && !doNorm && doLumi) {
      error_sys_max=0.033*0.033;
      error_sys_min=0.033*0.033;
    }
    for (uint i = 0; i < nun; i++) {
      graph[i]->GetPoint(ibin, x, y);
      errorl = graph[i]->GetErrorYlow(ibin);
      errorh = graph[i]->GetErrorYhigh(ibin);
      double error_max = errorh / y;
      double error_min = errorl / y;
      error_sys_max += error_max*error_max;
      error_sys_min += error_min*error_min;
    }
    g_un->SetPoint(ibin, x_nom, y_nom);
    g_un->SetPointError(ibin, 0, 0, y_nom * sqrt(error_sys_min), y_nom * sqrt(error_sys_max));
  }
  return g_un;
}

TGraphAsymmErrors* ConvertToTGraph(TH1D* h){
  TGraphAsymmErrors* g=new TGraphAsymmErrors(nbins);
  for (size_t ibin=0; ibin<nbins; ibin++){
    g->SetPoint(ibin, ibin, h->GetBinContent(ibin+1));
    g->SetPointError(ibin, 0, 0, h->GetBinError(ibin+1), h->GetBinError(ibin+1));
  }
  return g;
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

vector<TH1D*> EmptyHVec(int n) {
  vector<TH1D*> h;
  for (int i = 0; i < n; i++) {
    TH1D *phistartemp = new TH1D("phistar", "phistar", nbins, 0, nbins);
    phistartemp->Sumw2();
    h.push_back(phistartemp);
  }
  return h;
}

vector<TGraphAsymmErrors *> CreateCopy(vector<TGraphAsymmErrors *> graphvec) {
  vector<TGraphAsymmErrors *> newvec;
  for (int i = 0; i < ((int) graphvec.size()); i++) {//just getting rid of the warnings
    TGraphAsymmErrors *temp = (TGraphAsymmErrors *) graphvec[i]->Clone();
    newvec.push_back(temp);
  }
  return newvec;
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

  for (int i = 0; i < t->GetEntries(); i++) {
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
  delete t;
  return;
}


void ZShape2dMC(){
  double signal_weight = 3531.89/30459500.;//3504
  if (!doMG) signal_weight = 1966.7 / 42705454.; //

  TH1D* mc_truegen;
  vector<TH1D*> mc_truegen_cteq, mc_truegen_fsr_pileup;
  GetGen(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_fsr_pileup);
  
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

  vector<TGraphAsymmErrors *> g_MC_norm = CreateCopy(g_MC_phistar);
  vector<TGraphAsymmErrors *> g_MC_norm_cteq = CreateCopy(g_MC_phistar_cteq);
  vector<TGraphAsymmErrors *> g_MC_norm_fsr_pileup = CreateCopy(g_MC_phistar_fsr_pileup);
  NormalizeGraph(g_MC_phistar, 0);
  NormalizeGraph(g_MC_phistar_cteq, 0);
  NormalizeGraph(g_MC_phistar_fsr_pileup, 0);
  NormalizeGraph(g_MC_norm, 1);
  NormalizeGraph(g_MC_norm_cteq, 1);
  NormalizeGraph(g_MC_norm_fsr_pileup, 1);

  TGraphAsymmErrors* g_syst_mc_cteq = CalcTotalSysU_updown(g_MC_phistar_cteq,g_MC_phistar_cteq[0],1,1);
  TGraphAsymmErrors* g_syst_mc_fsr  = CalcTotalSysU_fsr(g_MC_phistar_fsr_pileup[0],g_MC_phistar[0]);
  TGraphAsymmErrors* g_syst_norm_mc_cteq = CalcTotalSysU_updown(g_MC_norm_cteq,g_MC_norm_cteq[0],1,1);
  TGraphAsymmErrors* g_syst_norm_mc_fsr  = CalcTotalSysU_fsr(g_MC_norm_fsr_pileup[0],g_MC_norm[0]);
  
  vector<TGraphAsymmErrors *> g_mc_syst;
  g_mc_syst.push_back(g_MC_phistar[0]);
  g_mc_syst.push_back(g_syst_mc_fsr);
  if (!doMG) g_mc_syst.push_back(g_syst_mc_cteq);
  vector<std::string> syst_mc_list;
  syst_mc_list.push_back("unfolding");
  syst_mc_list.push_back("fsr");
  if (!doMG) syst_mc_list.push_back("cteq");
  TGraphAsymmErrors* g_mc_final=GetMCFinal(g_mc_syst, syst_mc_list, 0);
  
  vector<TGraphAsymmErrors *> g_mc_syst_norm;
  g_mc_syst_norm.push_back(g_MC_norm[0]);
  g_mc_syst_norm.push_back(g_syst_norm_mc_fsr);
  if (!doMG) g_mc_syst_norm.push_back(g_syst_norm_mc_cteq);
  vector<std::string> syst_mc_list_norm;
  syst_mc_list_norm.push_back("unfolding");
  syst_mc_list_norm.push_back("fsr");
  if (!doMG) syst_mc_list_norm.push_back("cteq");
  TGraphAsymmErrors* g_mc_final_norm=GetMCFinal(g_mc_syst_norm, syst_mc_list_norm, 1);
  
  string textn="Output/TESTMC_Graph_Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  TFile tr2(textn.c_str(),"RECREATE");
  g_mc_final->Write();

  textn="Output/TESTMC_Graph_Norm_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  if (elec==0)textn+="Dressed.root";
  if (elec==1)textn+="Born.root";
  if (elec==2)textn+="Naked.root";
  TFile tr3(textn.c_str(),"RECREATE");
  g_mc_final_norm->Write();
}
