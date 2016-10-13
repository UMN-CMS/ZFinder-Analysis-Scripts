#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TMatrixD.h"
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

const int nbins=6*36;
const std::string Tag="";

void PlotCovM(bool norm=0, int elec=1, bool doMG=0){
  string list[]={"stat","meff","reff","teff","teff_m","teff_d","eff","mcstat_toy","mcstat_event","mcstat","bg_QCD","bg_ZZWZ","bg_toy","bg","pu","pt","fsr","unf","lumi","pdf","syst","tot"};
  std::string textn = "Output/CovM_";
  textn += Tag;
  if (norm) textn += "Norm_";
  else textn += "Abs_";
  if (doMG) textn += "MG_";
  else textn += "PH_";
  if (elec == 0)textn += "Dressed.root";
  if (elec == 1)textn += "Born.root";
  if (elec == 2)textn += "Naked.root";

  TFile f(textn.c_str());

  for (uint i=0; i<22;i++){
    if (norm && list[i]=="lumi")continue;
    if (doMG && list[i]=="pdf")continue;
    std::string matrixn = "CovM_";
    if (norm) matrixn += "Norm_"; 
    else  matrixn += "Abs_";
    matrixn +=list[i];
    TMatrixD CovM=*((TMatrixD*)f.Get(matrixn.c_str()));
    TH2D* Plot=new TH2D(matrixn.c_str(),matrixn.c_str(),nbins,0,nbins,nbins,0,nbins);
    for (uint j=0; j<nbins; j++){
      for (uint k=0; k<nbins; k++){
	Plot->SetBinContent(j+1,k+1,CovM(j,k));
      }
    }
    std::string plotname="CovM_2D_";
    plotname += Tag;
    if (norm) plotname += "Norm_";
    else plotname += "Abs_";
    if (doMG) plotname += "MG_";
    else plotname += "PH_";
    if (elec == 0)plotname += "Dressed_";
    if (elec == 1)plotname += "Born_";
    if (elec == 2)plotname += "Naked_";
    plotname+=list[i];
    TCanvas* C = new TCanvas(plotname.c_str(),plotname.c_str(),800,900);
    C->cd();
    gPad->SetPad("p1","p1",0,0, 1,1,kWhite,0,0);
    gPad->SetBottomMargin(0.10);
    gPad->SetTopMargin(0.10);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.15);
    Plot->SetStats(0);
    Plot->GetYaxis()->SetTitle("(#phi*,y)-bin");
    Plot->GetYaxis()->SetTitleOffset(1.5); 
    Plot->GetXaxis()->SetTitle("(#phi*,y)-bin");
    if (list[i]=="stat")Plot->SetTitle("Statistical uncertainty");
    if (list[i]=="unf")Plot->SetTitle("Unfolding uncertainty");
    if (list[i]=="eff")Plot->SetTitle("Efficiency uncertainty");
    if (list[i]=="meff")Plot->SetTitle("Medium electron ID eff. uncertainty");
    if (list[i]=="reff")Plot->SetTitle("GSF electron reconstrction eff. uncertainty");
    if (list[i]=="teff")Plot->SetTitle("Tight electron ID eff.  uncertainty");
    if (list[i]=="teff_m")Plot->SetTitle("Trigger eff. in data uncertainty");
    if (list[i]=="teff_d")Plot->SetTitle("Trigger eff. in MC uncertainty");
    if (list[i]=="mcstat")Plot->SetTitle("MC stat uncertainty");
    if (list[i]=="mcstat_toy")Plot->SetTitle("MC stat unfolding uncertainty");
    if (list[i]=="mcstat_event")Plot->SetTitle("MC stat event eff uncertainty");
    if (list[i]=="bg")Plot->SetTitle("Background uncertainty");
    if (list[i]=="bg_QCD")Plot->SetTitle("QCD background uncertainty");
    if (list[i]=="bg_ZZWZ")Plot->SetTitle("ZZ+WZ background uncertainty");
    if (list[i]=="bg_toy")Plot->SetTitle("Flavour independent background uncertainty");
    if (list[i]=="pt")Plot->SetTitle("Electron energy scale uncertainty");
    if (list[i]=="pu")Plot->SetTitle("PileUp uncertainty");
    if (list[i]=="fsr")Plot->SetTitle("FSR uncertainty");
    if (list[i]=="lumi")Plot->SetTitle("Luminosity uncertainty");
    if (list[i]=="pdf")Plot->SetTitle("PDF uncertainty");
    if (list[i]=="tot")Plot->SetTitle("Total uncertainty");
    if (list[i]=="syst")Plot->SetTitle("Total systematic uncertainty");
    Plot->Draw("COLZ");
    C->SaveAs(("Plots/"+plotname+".png").c_str());
    delete C;
  }
}
