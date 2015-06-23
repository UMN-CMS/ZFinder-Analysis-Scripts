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
const int Ntoys = 500 ;
const int Ntoys2 = 100 ;
const int elec = 0;
const int doMG = 1;
const Int_t RooUnfoldExtraInfo = 0; //set it to one to print the chi^2 stuff I think. at zero it supposedly makes the covariance matrix though the Roo code has a comment that it might be broke
const std::string Tag = "";

string MatrixName = "";
const double TEfficiencyEtaBins[11] = {-2.1, -2.0, -1.556, -1.442, -0.8, 0.0, 0.8, 1.442, 1.556, 2.0, 2.1};
const double TEfficiencyETBins[5] = {30, 40, 50, 70, 5000};
const double TEfficiencyData[10][4][2] = {
    {
        {0.741, 0.003},
        {0.773, 0.003},
        {0.780, 0.005},
        {0.790, 0.010}
    },
    {
        {0.734, 0.001},
        {0.772, 0.001},
        {0.786, 0.002},
        {0.792, 0.005}
    },
    {
        {0.725, 0.003},
        {0.821, 0.002},
        {0.809, 0.004},
        {0.848, 0.010}
    },
    {
        {0.893, 0.0005},
        {0.9396, 0.0003},
        {0.9509, 0.0006},
        {0.966, 0.001}
    },
    {
        {0.9213, 0.0004},
        {0.9528, 0.0002},
        {0.9601, 0.0004},
        {0.969, 0.001}
    },
    {
        {0.9174, 0.0004},
        {0.9473, 0.0004},
        {0.9561, 0.0004},
        {0.963, 0.001}
    },
    {
        {0.8964, 0.0005},
        {0.9424, 0.0003},
        {0.9533, 0.0006},
        {0.966, 0.001}
    },
    {
        {0.714, 0.003},
        {0.823, 0.002},
        {0.827, 0.004},
        {0.861, 0.010}
    },
    {
        {0.758, 0.001},
        {0.800, 0.001},
        {0.811, 0.002},
        {0.823, 0.005}
    },
    {
        {0.764, 0.003},
        {0.792, 0.002},
        {0.797, 0.005},
        {0.820, 0.010}
    }
};
const double TEfficiencyMC[10][4][2] = {
    {
        {0.734, 0.004},
        {0.769, 0.003},
        {0.771, 0.004},
        {0.760, 0.020}
    },
    {
        {0.736, 0.002},
        {0.768, 0.004},
        {0.779, 0.003},
        {0.789, 0.008}
    },
    {
        {0.791, 0.004},
        {0.847, 0.003},
        {0.850, 0.006},
        {0.870, 0.020}
    },
    {
        {0.9395, 0.0006},
        {0.9612, 0.0004},
        {0.9690, 0.0007},
        {0.980, 0.002}
    },
    {
        {0.9469, 0.0005},
        {0.9670, 0.0003},
        {0.9745, 0.0005},
        {0.982, 0.001}
    },
    {
        {0.9466, 0.0005},
        {0.9665, 0.0003},
        {0.9739, 0.0005},
        {0.982, 0.001}
    },
    {
        {0.9364, 0.0007},
        {0.9597, 0.0004},
        {0.9668, 0.0008},
        {0.979, 0.002}
    },
    {
        {0.779, 0.004},
        {0.841, 0.003},
        {0.842, 0.006},
        {0.860, 0.020}
    },
    {
        {0.749, 0.002},
        {0.786, 0.002},
        {0.798, 0.003},
        {0.810, 0.008}
    },
    {
        {0.737, 0.004},
        {0.769, 0.004},
        {0.779, 0.007},
        {0.820, 0.020}
    }
};

const double EfficiencyEtaBins[6] = {0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};
const double EfficiencyETBins[5] = {18, 30, 40, 50, 1000000};
const double EfficiencySF[5][4][3] = {
    {
        {0.982, 0.003, 0.012},
        {0.988, 0.001, 0.008},
        {0.990, 0.001, 0.004},
        {0.990, 0.001, 0.004}
    },
    {
        {0.993, 0.002, 0.012},
        {0.993, 0.001, 0.008},
        {0.993, 0.001, 0.004},
        {0.991, 0.001, 0.004}
    },
    {
        {1.016, 0.012, 0.020},
        {0.985, 0.004, 0.009},
        {0.987, 0.004, 0.004},
        {0.974, 0.009, 0.006}
    },
    {
        {0.988, 0.003, 0.012},
        {0.993, 0.002, 0.008},
        {0.992, 0.001, 0.004},
        {0.990, 0.003, 0.004}
    },
    {
        {1.002, 0.004, 0.012},
        {1.004, 0.002, 0.008},
        {1.005, 0.002, 0.004},
        {0.998, 0.004, 0.004}
    }
};
const double EfficiencyMediumSF[5][4][3] = {
    {
        {0.986, 0.002, 0.001},
        {1.002, 0.001, 0.001},
        {1.005, 0.001, 0.001},
        {1.004, 0.001, 0.001}
    },
    {
        {0.959, 0.003, 0.003},
        {0.980, 0.001, 0.001},
        {0.988, 0.001, 0.001},
        {0.988, 0.002, 0.002}
    },
    {
        {0.967, 0.007, 0.013},
        {0.950, 0.006, 0.007},
        {0.958, 0.005, 0.005},
        {0.966, 0.009, 0.009}
    },
    {
        {0.941, 0.005, 0.005},
        {0.967, 0.003, 0.003},
        {0.992, 0.002, 0.002},
        {1.000, 0.003, 0.003}
    },
    {
        {1.020, 0.003, 0.003},
        {1.021, 0.003, 0.003},
        {1.019, 0.002, 0.002},
        {1.022, 0.004, 0.004}
    }
};
const double EfficiencyTightSF[5][4][3] = {
    {
        {0.960, 0.003, 0.003},
        {0.978, 0.001, 0.001},
        {0.981, 0.001, 0.001},
        {0.982, 0.002, 0.002}
    },
    {
        {0.936, 0.004, 0.004},
        {0.958, 0.002, 0.002},
        {0.969, 0.001, 0.001},
        {0.969, 0.002, 0.002}
    },
    {
        {0.933, 0.015, 0.017},
        {0.907, 0.008, 0.008},
        {0.904, 0.004, 0.004},
        {0.929, 0.011, 0.011}
    },
    {
        {0.879, 0.007, 0.007},
        {0.909, 0.003, 0.003},
        {0.942, 0.002, 0.002},
        {0.957, 0.004, 0.004}
    },
    {
        {0.974, 0.004, 0.004},
        {0.987, 0.004, 0.004},
        {0.991, 0.003, 0.003},
        {0.999, 0.005, 0.005}
    }
};

const bool DoYSeperation = true;
const double YSeperation[] = {.1, .4}; //elements represent seperation between seperation would be 0<a<b<c<inf with abc elents of YSeperation
const double phistarBinsOne[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 10};
//To handle 2D measurements will be making a second one that handles more
const double* phistarBins;
const size_t nYSeper = (sizeof (YSeperation) / sizeof (phistarBinsOne[0]));

const size_t nphistar = ((sizeof (phistarBinsOne) / sizeof (phistarBinsOne[0])) - 1)*(nYSeper*DoYSeperation + 1); //this allows the 1D graphs to have the data for multiple Y values, with the first 35 being the 0<Y<a from Y seperation
//this makes it so if we don't want to look at YSeperation

std::string File_Data = "/afs/cern.ch/work/r/ruckstuh/public/Data_R9.root";
std::string File_Data_Pt_L = "/afs/cern.ch/work/r/ruckstuh/public/Data_Low_R9.root";
std::string File_Data_Pt_H = "/afs/cern.ch/work/r/ruckstuh/public/Data_High_R9.root";
std::string File_Signal_reco = "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco_R9.root";
std::string File_Signal_gen = "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Signal_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_born_gen.root";
std::string File_Signal_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_naked_gen.root";
std::string File_Powheg_reco = "/afs/cern.ch/work/r/ruckstuh/public/Powheg_reco.root";
std::string File_Powheg_gen = "/afs/cern.ch/work/r/ruckstuh/public/Powheg_gen.root";
std::string File_Powheg_gen_born = "/afs/cern.ch/work/r/ruckstuh/public/Powheg_born_gen.root";
std::string File_Powheg_gen_bare = "/afs/cern.ch/work/r/ruckstuh/public/Powheg_naked_gen.root";
std::string File_tt = "/afs/cern.ch/work/r/ruckstuh/public/TTbar.root";
std::string File_tautau = "/afs/cern.ch/work/r/ruckstuh/public/TauTau.root";
std::string File_tbarw = "/afs/cern.ch/work/r/ruckstuh/public/TbarW.root";
std::string File_tw = "/afs/cern.ch/work/r/ruckstuh/public/TW.root";
std::string File_ww = "/afs/cern.ch/work/r/ruckstuh/public/WW.root";
std::string File_wz = "/afs/cern.ch/work/r/ruckstuh/public/WZ.root";
std::string File_zz = "/afs/cern.ch/work/r/ruckstuh/public/ZZ.root";

std::string reco_name = "Combined Single Reco";
std::string reco_name_en_l = "Combined Single Lowered Threshold Reco";
std::string reco_name_en_h = "Combined Single Higher Threshold Reco";
//std::string reco_name_en="Combined Single Reco";
std::string gen_name = "Combined Gen Cuts Reco";

double OffSetter(double phistar, double Z_Y) {
    Z_Y = fabs(Z_Y);
    //cout<<"Z_Y :"<<Z_Y;
    for (size_t TwoDCheck = 0; TwoDCheck < nYSeper; TwoDCheck++) {
        if ((TwoDCheck + 1) < nYSeper) {

            if (YSeperation[TwoDCheck] < Z_Y && YSeperation[TwoDCheck + 1] >= Z_Y) phistar += (TwoDCheck + 1) *(1 + phistarBinsOne[(nphistar / (nYSeper + 1))]); //so for doing the checks for the phistar bin useful 
        } else if (YSeperation[TwoDCheck] < Z_Y) {
            phistar += (TwoDCheck + 1) *(1 + phistarBinsOne[(nphistar / (nYSeper + 1))]);
            //cout<<"end thingy :"<<phistarBinsOne[(nphistar / (nYSeper + 1))]<<" and before that :"<<phistarBinsOne[(nphistar / (nYSeper + 1))-1]<<endl;}
        }
    }
    //cout<<" and our phistar is :"<<phistar<<endl;
    return phistar;
}

void MakeCovTMatrix(TH2D* CovHisto, TMatrixD & CovMatrix)//this is for the Toy MC sample since I already made a histogram but it was asked I make a TMatrix 
{// since I decided histograms were useful in there own way and I wanted to save them this seemed the best way to deal with them

    for (size_t iphistarx = 0; iphistarx < nphistar; iphistarx++) {
        for (size_t iphistary = 0; iphistary < nphistar; iphistary++) {
            CovMatrix(iphistarx, iphistary) = CovHisto->GetBinContent(iphistarx + 1, iphistary + 1);
        }
    }

}

void MakeCovTMatrix(TGraphAsymmErrors* OrginalGraph, TMatrixD & CovMatrix)// so if we don't have any toys just grab the error from the TGraph
{
    double IntValue = 0; //this will be used to normalize the graph. 
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y;
        OrginalGraph->GetPoint(iphistar, x, y);
        IntValue += y;
    }

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {

        double yerror;
        yerror = OrginalGraph->GetErrorYhigh(iphistar);
        yerror = yerror / IntValue; //normalizing the error
        CovMatrix(iphistar, iphistar) = (yerror * yerror);
    }

}

void MakeCovTMatrix(TH1D* OrginalGraph, TMatrixD & CovMatrix)// so if we don't have any toys just grab the error from the TGraph
{

    for (size_t iphistar = 1; iphistar < nphistar; iphistar++) {

        double yerror;
        yerror = OrginalGraph->GetBinError(iphistar + 1); //plus one since the zero bun is overflow

        //yerror = yerror / (OrginalGraph->GetBinContent(iphistar + 1)); //Nicole said to do it this way, 
        yerror = yerror / (OrginalGraph->Integral()); //seemed to make more sense to me to change it by the normalized graph.
        CovMatrix(iphistar, iphistar) = (yerror * yerror);
    }

}

void MakeCovTMatrix(TGraphAsymmErrors* OrginalGraph, TMatrixD& CovMatrix, double PercentError)// Since there are a couple of things that all change by a constant percentage this just makes a matrix down the diagonal with that percentage 
{

    double IntValue = 0; //this will be used to normalize the graph. 
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y;
        OrginalGraph->GetPoint(iphistar, x, y);
        IntValue += y;
    }

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y;
        double yerror;
        OrginalGraph->GetPoint(iphistar, x, y);
        yerror = y * PercentError / IntValue;

        CovMatrix(iphistar, iphistar) = (yerror * yerror);
    }

}

double Return_RMS(double mean_sq, double mean) {
    return sqrt(mean_sq - mean * mean);
}

TH1D * ConvertToHist(TGraphAsymmErrors* g, std::string name) {
    TH1D* h_temp;

    h_temp = new TH1D(name.c_str(), name.c_str(), nphistar, phistarBins);
    h_temp->Sumw2();

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y;
        g->GetPoint(iphistar, x, y);

        double error = g->GetErrorYhigh(iphistar);

        h_temp->SetBinContent(iphistar + 1, y);
        h_temp->SetBinError(iphistar + 1, error);
    }
    return h_temp;
}

void Cov_Histo_Creator(vector<TH1D *> HistoVector, TH2D* Covar_hist, TH2D * Correl_hist) {//makes the covariance matrix if given TH1D.
    int Niter = HistoVector.size(); //number of toys
    const int nphistar_bins = nphistar;
    //const int nphistar_bins = 34;
    //double phistar_var[nphistar_bins] = {0.0, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.051, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
    double Mean_PROD_SUM[nphistar_bins ][nphistar_bins ] = {
        {0.}
    };

    double Mean_SUM[nphistar_bins ] = {0.}; //holds the sum of all the toys for the bin
    double Mean_SQ_SUM[nphistar_bins ] = {0.};

    double Mean_j = 0.;
    double Mean_k = 0;
    double RMS_j = 0.;
    double RMS_k = 0;
    double Correlation = 0.;
    double Covariance = 0.;
    //TH2D *Covar_hist =  new TH2D("Covar_hist","Covar_hist",nphistar_bins-1,phistar_var,nphistar_bins-1,phistar_var);
    //TH2D *Correl_hist = new TH2D("Correl_hist", "Correl_hist", nphistar_bins - 1, phistar_var, nphistar_bins - 1, phistar_var);
    for (int i = 0; i < Niter; i++) {
        //sprintf(name,"hReco_%i",i);
        TH1F* h = (TH1F*) HistoVector[i];
        h->Scale(1 / h->Integral());
        for (int j = 1; j < 1 + ((int) nphistar); j++) {
            if (i == 0) {
                Mean_SUM[j - 1] = 0;
                Mean_SQ_SUM[j - 1] = 0;
            }
            Mean_SUM[j - 1] += h->GetBinContent(j);
            //Mean_SQ_SUM[j - 1] += pow(h->GetBinContent(j), 2);
            Mean_SQ_SUM[j - 1] += h->GetBinContent(j) * h->GetBinContent(j);
            for (int k = 1; k < 1 + ((int) nphistar); k++) {
                if (i == 0) {
                    Mean_PROD_SUM[j - 1][k - 1] = 0;
                }
                Mean_PROD_SUM[j - 1][k - 1] += h->GetBinContent(j) * h->GetBinContent(k);
            }
        }
    }
    TH1F* testhisto = (TH1F*) HistoVector[0];
    testhisto->Scale(1 / testhisto->Integral());

    if (MatrixName.length() > 2)cout << "Matrix named " << MatrixName << endl;
    for (int j = 1; j < 1 + ((int) nphistar); j++) {
        Mean_j = Mean_SUM[j - 1] / Niter;
        RMS_j = Return_RMS(Mean_SQ_SUM[j - 1] / Niter, Mean_j);
        for (int k = 1; k < 1 + ((int) nphistar); k++) {
            Mean_k = Mean_SUM[k - 1] / Niter;
            RMS_k = Return_RMS(Mean_SQ_SUM[k - 1] / Niter, Mean_k);
            Covariance = ((Mean_PROD_SUM[j - 1][k - 1] / Niter) - Mean_j * Mean_k);
            Correlation = ((Mean_PROD_SUM[j - 1][k - 1] / Niter) - Mean_j * Mean_k) / (RMS_j * RMS_k);
            if (j == k && MatrixName.length() > 2 && debug) {
                //cout << " okay the mean product is :" << Mean_PROD_SUM[j - 1][k - 1] / Niter << "  and the other half is " << Mean_j * Mean_k;
                cout << " and the covariance  :" << Covariance << "for bin number :" << k;
                cout << " and the orginal is :" << testhisto->GetBinContent(j) << endl;
                double sqrtCov = sqrt(Covariance);
                cout << "this gives us " << sqrtCov << endl << endl;

            }
            Covar_hist->SetBinContent(j, k, Covariance);
            Correl_hist->SetBinContent(j, k, Correlation);
        }
        //return Covar_hist;
    }
}

void Cov_Histo_Creator(vector<TGraphAsymmErrors*> test, TH2D* Covar_hist, TH2D * Correl_hist) {
    int Niter = test.size(); //number of toys
    vector<TH1D*> HistoHolder;

    cout << "test Cov 1" << endl;
    for (int w = 0; w < Niter; w++) {
        
        //TH1F* h = (TH1F*) HistoVector[w];
        ostringstream namesnumber;
        namesnumber << w;
        string name = "Histo" + namesnumber.str();
        cout << "test Cov 2" << endl;
        TH1D* h = ConvertToHist(test[w], name);
        cout << "test Cov 2.1" << endl;
        //for (int i = 1; (i < ((int) nphistar)) && w == 0 && MatrixName.length() > 2 && debug; i++) {
       //     cout << "and this should be the same but isn't? :" << h->GetBinContent(i) << "  bin number " << i << endl;
        //}
        HistoHolder.push_back(h);
    }
    cout << "test Cov 3" << endl;

    Cov_Histo_Creator(HistoHolder, Covar_hist, Correl_hist);
        cout << "test Cov 4" << endl;
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
    return EffToys;
}

vector<TH1D*> GetToyBg(TH1F * bg_sf_full) {
    vector<TH1D*> bg_sf;
    for (int t = 0; t < Ntoys + 1; t++) {
        TH1D* bg_sf_temp = new TH1D("bg", "bg", nphistar, phistarBins);
        bg_sf_temp->Sumw2();
        for (uint i = 0; i < nphistar; i++) {
            if (i == nphistar - 1) bg_sf_temp->SetBinContent(i + 1, 1.);
            else {
                double val = bg_sf_full->GetBinContent(i + 1);
                if (t == 0)bg_sf_temp->SetBinContent(i + 1, val);
                else {
                    double error = bg_sf_full->GetBinError(i + 1);
                    double x = gRandom->Gaus(val, error);
                    bg_sf_temp->SetBinContent(i + 1, x);
                }
            }
            bg_sf_temp->SetBinError(i + 1, 0.);
        }
        bg_sf.push_back(bg_sf_temp);
    }
    return bg_sf;
}

TH1D * GetBgSS(TH1D * ss_full) {
    TH1D* bg_ss = new TH1D("bg", "bg", nphistar, phistarBins);
    bg_ss->Sumw2();
    for (uint i = 0; i < nphistar; i++) {
        if (i == nphistar - 1) bg_ss->SetBinContent(i + 1, 0.);
        else {
            double val = ss_full->GetBinContent(i + 1);
            bg_ss->SetBinContent(i + 1, val);
        }
        bg_ss->SetBinError(i + 1, 0.);
    }
    return bg_ss;
}

void NormalizeGraph(vector<TGraphAsymmErrors*> &graph, bool doNorm = 0) {
    for (size_t i = 0; i < graph.size(); i++) {
        double xstot = 0;
        double x, y, errorl, errorh;
        for (size_t iphistar = 0; iphistar < nphistar - 1; iphistar++) {
            graph[i]->GetPoint(iphistar, x, y);
            xstot += y;
        }
        //    cout<<i<<":Normalisation:"<<xstot<<endl;
        for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
            double dphistar = phistarBins[iphistar + 1] - phistarBins[iphistar];
            double norm = dphistar*xstot;
            double errorbin2_h = 0;
            double errorbin2_l = 0;
            if (doNorm) {
                for (int j = 0; j < ((int) nphistar); j++) {
                    double errorbin_j_on_i_h = 0;
                    double errorbin_j_on_i_l = 0;
                    if (((int) iphistar) == j) {
                        errorbin_j_on_i_h = graph[i]->GetErrorYhigh(j) * ((1 / xstot)-(1 / (xstot * xstot)));
                        errorbin_j_on_i_l = graph[i]->GetErrorYlow(j) * ((1 / xstot)-(1 / (xstot * xstot)));
                    } else {
                        errorbin_j_on_i_h = graph[i]->GetErrorYhigh(j)*(1 / (xstot * xstot));
                        errorbin_j_on_i_l = graph[i]->GetErrorYlow(j)*(1 / (xstot * xstot));
                    }
                    errorbin2_h += errorbin_j_on_i_h*errorbin_j_on_i_h;
                    errorbin2_l += errorbin_j_on_i_l*errorbin_j_on_i_l;
                }
            }
            if (!doNorm) norm = dphistar;
            graph[i]->GetPoint(iphistar, x, y);
            errorl = graph[i]->GetErrorYlow(iphistar);
            errorh = graph[i]->GetErrorYhigh(iphistar);
            graph[i]->SetPoint(iphistar, x, y / norm);
            graph[i]->SetPointError(iphistar, 0, 0, errorl / norm, errorh / norm);
            if (doNorm) graph[i]->SetPointError(iphistar, 0, 0, sqrt(errorbin2_l) / dphistar, sqrt(errorbin2_h) / dphistar);
        }
    }
}

TGraphAsymmErrors * CalcTotalSysU_updown(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal, bool same = 0, bool cteq = 0) { // adds everything in quadrature (up down)
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, errorl, errorh;
        graph_nominal->GetPoint(iphistar, x, y);
        g_un->SetPoint(iphistar, x, y);
        errorl = 0;
        errorh = 0;
        int last = 0;
        int now = 0;
        double diflast = 0;
        for (size_t i = same; i < graph.size(); i++) {
            now = 0;
            double xtemp, ytemp;
            graph[i]->GetPoint(iphistar, xtemp, ytemp);
            if (xtemp != x) cout << "This is really weird and wrong" << endl;
            double dif = ytemp - y;
            if (dif / y > 0.1)cout << "WOW very large error: " << dif / y << " " << y << "  " << ytemp << endl;
            if (dif > 0) {
                errorh = sqrt(errorh * errorh + dif * dif);
                now = 1.0;
            }
            if (dif < 0) {
                errorl = sqrt(errorl * errorl + dif * dif);
                now = -1.0;
            }
            if (i > 0 && i % 2 != same) {
                if (last != 0 && now != 0 && last == now && fabs(diflast / y) > 0.0001 && fabs(dif / y) > 0.0001) cout << "both change go in the same dirrection: " << i << "  " << iphistar << "  " << dif / y << " " << diflast / y << endl;
            }
            last = now;
            diflast = dif;
        }
        if (cteq) {
            errorl = errorl / 1.645;
            errorh = errorh / 1.645;
        }
        g_un->SetPointError(iphistar, 0, 0, errorl, errorh);
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_comb3(TGraphAsymmErrors* graph, TGraphAsymmErrors* graph_nominal, TGraphAsymmErrors* graph_mctoy, bool add = 1) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, errorl, errorh;
        double x2, y2, errorl2, errorh2;
        double errorl3, errorh3;
        graph_nominal->GetPoint(iphistar, x, y);
        graph->GetPoint(iphistar, x2, y2);
        if (x != x2 || y != y2) cout << "eff stat: nominal values don't agree: " << x << " " << x2 << " : " << y << " " << y2 << endl;
        g_un->SetPoint(iphistar, x, y);
        errorl = graph_nominal->GetErrorYlow(iphistar);
        errorh = graph_nominal->GetErrorYhigh(iphistar);
        errorl2 = graph->GetErrorYlow(iphistar);
        errorh2 = graph->GetErrorYhigh(iphistar);
        errorl3 = graph_mctoy->GetErrorYlow(iphistar);
        errorh3 = graph_mctoy->GetErrorYhigh(iphistar);
        if (!add) {
            if (errorl2 < errorl || errorh2 < errorh) cout << "eff stat: uncertainty size don't agree: " << errorl << " " << errorl2 << " : " << errorh << " " << errorh2 << endl;
            g_un->SetPointError(iphistar, 0, 0, sqrt((errorl2 * errorl2)-(errorl * errorl)+(errorl3 * errorl3)), sqrt((errorh2 * errorh2)-(errorh * errorh)+(errorh3 * errorh3)));
        } else {
            g_un->SetPointError(iphistar, 0, 0, sqrt((errorl2 * errorl2)+(errorl * errorl)+(errorl3 * errorl3)), sqrt((errorh2 * errorh2)+(errorh * errorh)+(errorh3 * errorh3)));
        }
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_comb5(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2, TGraphAsymmErrors* graph3, TGraphAsymmErrors* graph4, TGraphAsymmErrors * graph5) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, errorl, errorh;
        double errorl2, errorh2;
        double errorl3, errorh3;
        double errorl4, errorh4;
        double errorl5, errorh5;
        graph1->GetPoint(iphistar, x, y);
        g_un->SetPoint(iphistar, x, y);
        errorl = graph1->GetErrorYlow(iphistar);
        errorh = graph1->GetErrorYhigh(iphistar);
        errorl2 = graph2->GetErrorYlow(iphistar);
        errorh2 = graph2->GetErrorYhigh(iphistar);
        errorl3 = graph3->GetErrorYlow(iphistar);
        errorh3 = graph3->GetErrorYhigh(iphistar);
        errorl4 = graph4->GetErrorYlow(iphistar);
        errorh4 = graph4->GetErrorYhigh(iphistar);
        errorl5 = graph5->GetErrorYlow(iphistar);
        errorh5 = graph5->GetErrorYhigh(iphistar);
        double errorltot = sqrt((errorl * errorl)+(errorl2 * errorl2)+(errorl3 * errorl3)+(errorl4 * errorl4)+(errorl5 * errorl5));
        double errorhtot = sqrt((errorh * errorh)+(errorh2 * errorh2)+(errorh3 * errorh3)+(errorh4 * errorh4)+(errorh5 * errorh5));
        g_un->SetPointError(iphistar, 0, 0, errorltot, errorhtot);
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_toyvariation(vector<TGraphAsymmErrors*> graph, bool useMedian = 0, int part = 1, bool eff = 0) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y;
        vector<double> phis;
        graph[0]->GetPoint(iphistar, x, y);
        phis.push_back(y);
        int nt = Ntoys;
        if (eff) nt = Ntoys2;
        uint gs = 1 + (nt * (part - 1));
        for (size_t i = gs; i < gs + nt; i++) {
            graph[i]->GetPoint(iphistar, x, y);
            phis.push_back(y);
        }
        g_un->SetPoint(iphistar, x, TMath::Mean(phis.begin(), phis.end()));
        double rms = TMath::RMS(phis.begin(), phis.end());
        g_un->SetPointError(iphistar, 0, 0, rms, rms);
        if (useMedian && !eff) {
            std::sort(phis.begin(), phis.end());
            double med = phis[(nt - 1) / 2];
            int idx68 = 0.1587 * (nt) - 0.5;
            double min68 = phis[idx68];
            double max68 = phis[(nt) - idx68];
            g_un->SetPoint(iphistar, x, med);
            g_un->SetPointError(iphistar, 0, 0, med - min68, max68 - med);
        }
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_toymc(vector<TGraphAsymmErrors*> graph, TGraphAsymmErrors* graph_nominal, bool isBG = 0, int part = 1, bool eff = 0) {
    TGraphAsymmErrors* g_toy = CalcTotalSysU_toyvariation(graph, 1, part, eff);
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, xmc, ymc;
        graph_nominal->GetPoint(iphistar, x, y);
        g_un->SetPoint(iphistar, x, y);
        g_toy->GetPoint(iphistar, xmc, ymc);
        double eh = ymc + g_toy->GetErrorYhigh(iphistar);
        double el = ymc - g_toy->GetErrorYlow(iphistar);
        eh = eh - y;
        el = el - y;
        // cout << iphistar << " " << ymc << " " << eh << " " << el << " " << y << endl;
        if (eh < 0) eh = 0;
        if (el > 0) el = 0;
        //cout << iphistar << " " << ymc << " " << eh << " " << el << " " << y << endl;
        if (isBG) {
            double ymc1, ymc2;
            double ymc3, ymc4;
            graph[graph.size() - 1]->GetPoint(iphistar, x, ymc1);
            graph[graph.size() - 2]->GetPoint(iphistar, x, ymc2);
            graph[graph.size() - 3]->GetPoint(iphistar, x, ymc3);
            graph[graph.size() - 4]->GetPoint(iphistar, x, ymc4);
            if (ymc1 > y) eh = sqrt(pow(eh, 2) + pow((ymc1 - y), 2));
            if (ymc2 > y) eh = sqrt(pow(eh, 2) + pow((ymc2 - y), 2));
            if (ymc3 > y) eh = sqrt(pow(eh, 2) + pow((ymc3 - y), 2));
            if (ymc4 > y) eh = sqrt(pow(eh, 2) + pow((ymc4 - y), 2));
            if (ymc1 < y) el = sqrt(pow(el, 2) + pow((ymc1 - y), 2));
            if (ymc2 < y) el = sqrt(pow(el, 2) + pow((ymc2 - y), 2));
            if (ymc3 < y) el = sqrt(pow(el, 2) + pow((ymc3 - y), 2));
            if (ymc4 < y) el = sqrt(pow(el, 2) + pow((ymc4 - y), 2));
        }
        g_un->SetPointError(iphistar, 0, 0, el, eh);
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_fsr(TGraphAsymmErrors* graph, TGraphAsymmErrors * graph_nominal) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, xtemp, ytemp, error;
        graph_nominal->GetPoint(iphistar, x, y);
        // g_un->SetPoint(iphistar,x,y);
        g_un->SetPoint(iphistar, x, y);
        graph->GetPoint(iphistar, xtemp, ytemp);
        if (xtemp != x) cout << "This is really weird and wrong" << endl;
        error = fabs(ytemp - y);
        g_un->SetPointError(iphistar, 0, 0, error, error);
    }
    return g_un;
}

TGraphAsymmErrors * CalcTotalSysU_pileup(TGraphAsymmErrors* graph_1, TGraphAsymmErrors* graph_2, TGraphAsymmErrors * graph_nominal) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, xtemp1, ytemp1, xtemp2, ytemp2;
        double errorh = 0;
        double errorl = 0;
        graph_nominal->GetPoint(iphistar, x, y);
        // g_un->SetPoint(iphistar,x,y);
        g_un->SetPoint(iphistar, x, y);
        graph_1->GetPoint(iphistar, xtemp1, ytemp1);
        graph_2->GetPoint(iphistar, xtemp2, ytemp2);
        if (xtemp1 != x) cout << "This is really weird and wrong" << endl;
        if (xtemp2 != x) cout << "This is really weird and wrong" << endl;
        if (y - ytemp1 > 0) errorl = fabs(ytemp1 - y);
        else errorh = fabs(ytemp1 - y);
        if (y - ytemp2 > 0) {
            if (errorl == 0) errorl = fabs(ytemp2 - y);
            else {
                if (errorl < fabs(ytemp2 - y)) errorl = fabs(ytemp2 - y);
                cout << "pile-up: both errors on same side:" << y - ytemp1 << " " << y - ytemp2 << endl;
            }
        } else {
            if (errorh == 0) errorh = fabs(ytemp2 - y);
            else {
                if (errorh < fabs(ytemp2 - y)) errorh = fabs(ytemp2 - y);
                cout << "pile-up: both errors on same side:" << y - ytemp1 << " " << y - ytemp2 << endl;
            }
        }
        g_un->SetPointError(iphistar, 0, 0, errorl, errorh);
    }
    return g_un;
}

TGraphAsymmErrors * GetDataFinal(vector<TGraphAsymmErrors *> graph, vector<std::string> slist, bool doNorm = 0, bool doLumi = 1) {
    TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);
    uint nun = graph.size();
    ofstream outputfile;
    std::string textname = "Table_Un_";
    textname += Tag;
    if (doNorm) textname += "Norm_";
    else textname += "Abs_";
    if (doMG) textname += "MG_";
    else textname += "PH_";
    if (elec == 0)textname += "Dressed.txt";
    if (elec == 1)textname += "Born.txt";
    if (elec == 2)textname += "Naked.txt";
    outputfile.open(textname.c_str());
    outputfile << "$\\phi^*$ range & MC stat. & Pile-up & Background & Energy scale & Efficiencies ";
    if (!doMG) outputfile << "& PDF ";
    outputfile << "& Total syst. & Stat. & Total \\\\ \\hline" << "\n";
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double h_un = 0;
        double h_pdf = 0;
        double h_pt = 0;
        double l_pt = 0;
        double h_eff = 0;
        double h_bg = 0;
        double h_fsr = 0;
        double h_pu = 0;
        double h_mc = 0;
        double l_eff = 0;
        double l_bg = 0;
        double l_fsr = 0;
        double l_pu = 0;
        double l_mc = 0;
        double l_pdf = 0; //Set but not used, delete? 
        double x, y, errorl, errorh;
        graph[0]->GetPoint(iphistar, x, y);
        double x_nom = x;
        double y_nom = y;
        ;
        double error_sys_max = 0;
        double error_sys_min = 0;
        double error_sys2_max = 0;
        double error_sys2_min = 0;
        for (uint i = 0; i < nun; i++) {
            graph[i]->GetPoint(iphistar, x, y);
            errorl = graph[i]->GetErrorYlow(iphistar);
            errorh = graph[i]->GetErrorYhigh(iphistar);
            double error_max = errorh / y;
            double error_min = errorl / y;
            error_sys_max += error_max*error_max;
            error_sys_min += error_min*error_min;
            if (slist[i] != "unfolding") {
                error_sys2_max += error_max*error_max;
                error_sys2_min += error_min*error_min;
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
            if (slist[i] == "fsr") {
                h_fsr = error_max;
                l_fsr = error_min;
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
            //cout<<error_sys_min<<endl;
        }
        if (!doNorm && doLumi) error_sys_max += 0.026 * 0.026;
        if (!doNorm && doLumi) error_sys_min += 0.026 * 0.026;
        if (!doNorm && doLumi) error_sys2_max += 0.026 * 0.026;
        if (!doNorm && doLumi) error_sys2_min += 0.026 * 0.026;
        //cout<<error_sys_min<<endl;
        g_un->SetPoint(iphistar, x_nom, y_nom);
        g_un->SetPointError(iphistar, 0, 0, y_nom * sqrt(error_sys_min), y_nom * sqrt(error_sys_max));
        if (doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f : %9.2f : %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n", int(iphistar), y_nom, sqrt(error_sys_max) * y_nom, sqrt(error_sys_min) * y_nom, sqrt(error_sys_max)*100., h_un * 100., h_mc * 100., h_pt * 100., h_eff * 100., h_bg * 100., h_fsr * 100., h_pu * 100.);
        if (doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f : %9.2f : %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n", int(iphistar), y_nom, sqrt(error_sys_max) * y_nom, sqrt(error_sys_min) * y_nom, sqrt(error_sys_min)*100., h_un * 100., l_mc * 100., l_pt * 100., l_eff * 100., l_bg * 100., l_fsr * 100., l_pu * 100.);
        if (!doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f: %9.2f: %9.2f: %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n", int(iphistar), y_nom, sqrt(error_sys_max) * y_nom, sqrt(error_sys_min) * y_nom, sqrt(error_sys_max)*100., h_un * 100., 2.6, h_mc * 100., h_pt * 100., h_eff * 100., h_bg * 100., h_fsr * 100., h_pu * 100.);
        if (!doNorm) printf("%10d : %10.3f + %5.3f - %5.3f : %9.2f: %9.2f: %9.2f: %9.2f : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% : %9.2f%% \n", int(iphistar), y_nom, sqrt(error_sys_max) * y_nom, sqrt(error_sys_min) * y_nom, sqrt(error_sys_min)*100., h_un * 100., 2.6, l_mc * 100., l_pt * 100., l_eff * 100., l_bg * 100., l_fsr * 100., l_pu * 100.);

        outputfile << std::fixed;
        outputfile << std::setprecision(3) << phistarBins[iphistar] << "-" << phistarBins[iphistar + 1] << " & " << std::setprecision(2) << h_mc * 100. << " & " << h_pu * 100. << " & " << h_bg * 100. << " & " << h_pt * 100. << " & " << h_eff * 100.;
        if (!doMG) outputfile << " & " << h_pdf * 100.;
        outputfile << " & " << sqrt(error_sys2_max)*100 << " & " << h_un * 100. << " & " << sqrt(error_sys_max)*100 << "  \\\\ \\hline" << "\n";
    }
    outputfile.close();
    return g_un;
}

void PlotEff(TH1D * h_eff) {
    TCanvas* Efficiency = new TCanvas("Efficiency", "Efficiency", 800, 900);
    Efficiency->cd();
    Efficiency->SetLogx();
    h_eff->GetXaxis()->SetRangeUser(0.001, 3.2);
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
    h_eff->SetBit(TH1::kNoTitle, true);
    h_eff->SetLineColor(1);
    h_eff->Draw();
    TLegend* leg_eff = new TLegend(0.45, 0.77, 0.85, 0.91);
    leg_eff->SetFillStyle(0);
    leg_eff->SetBorderSize(0);
    leg_eff->SetLineWidth(1);
    leg_eff->SetNColumns(1);
    leg_eff->SetTextFont(42);
    leg_eff->SetTextSize(0.04);
    leg_eff->AddEntry(h_eff, "MadGraph", "PL");
    leg_eff->Draw();
    return;
}

void PrintBG(TH1D* Data, TH1D* h_tt, TH1D* h_tautau, TH1D* h_tbarw, TH1D* h_tw, TH1D* h_ww, TH1D* h_wz, TH1D* h_zz, TH1D * h_QCD) {
    double data_sel = Data->GetSumOfWeights();
    double tt_sel = h_tt->GetSumOfWeights();
    double tautau_sel = h_tautau->GetSumOfWeights();
    double tbarw_sel = h_tbarw->GetSumOfWeights();
    double tw_sel = h_tw->GetSumOfWeights();
    double ww_sel = h_ww->GetSumOfWeights();
    double wz_sel = h_wz->GetSumOfWeights();
    double zz_sel = h_zz->GetSumOfWeights();
    double qcd_sel = h_QCD->GetSumOfWeights()*2;

    double t_bg = tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel;

    cout << "Weighted number of events:" << endl;
    cout << "Data: " << data_sel << "  tt: " << tt_sel << "  tautau: " << tautau_sel << "  tbarw: " << tbarw_sel << "  tw: " << tw_sel << " singletop: " << tbarw_sel + tw_sel << "  ww: " << ww_sel << "  wz: " << wz_sel << "  zz: " << zz_sel << "  qcd: " << qcd_sel << endl;
    cout << "ratio:" << "  tt: " << tt_sel * 100. / data_sel << "  tautau: " << tautau_sel * 100. / data_sel << "  tbarw: " << tbarw_sel * 100. / data_sel << "  tw: " << tw_sel * 100. / data_sel << " singletop: " << (tbarw_sel + tw_sel)*100. / data_sel << "  ww: " << ww_sel * 100. / data_sel << "  wz: " << wz_sel * 100. / data_sel << "  zz: " << zz_sel * 100. / data_sel << "  qcd: " << qcd_sel * 100. / data_sel << "data: " << (data_sel - t_bg)*100. / data_sel << endl;
    cout << "ratio of bg:" << "  tt: " << tt_sel * 100. / t_bg << "  tautau: " << tautau_sel * 100. / t_bg << "  tbarw: " << tbarw_sel * 100. / t_bg << "  tw: " << tw_sel * 100. / t_bg << " singletop: " << (tbarw_sel + tw_sel)*100. / t_bg << "  ww: " << ww_sel * 100. / t_bg << "  wz: " << wz_sel * 100. / t_bg << "  zz: " << zz_sel * 100. / t_bg << "  qcd: " << qcd_sel * 100. / t_bg << endl;
    cout << "total: " << t_bg << " " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel) << "  persentage: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel + qcd_sel)*100. / data_sel << endl;
    cout << "wz+zz: " << (wz_sel + zz_sel) << "  persentage: " << (wz_sel + zz_sel)*100. / data_sel << " of background: " << (wz_sel + zz_sel)*100. / (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel) << endl;
    cout << "e+mu: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel) << "  persentage: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel)*100. / data_sel << " of background: " << (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel)*100. / (tt_sel + tautau_sel + tbarw_sel + tw_sel + ww_sel + wz_sel + zz_sel) << endl;
}

void GetToyResponse(vector<RooUnfoldResponse*> &BinM, TH2D * BinMigration) {
    // RooUnfoldResponse* BinM1  =new RooUnfoldResponse (h_reco,h_gen,BinMigration);
    RooUnfoldResponse* BinM1 = new RooUnfoldResponse(0, 0, BinMigration);
    BinM.push_back(BinM1);
    double x;
    for (int i = 0; i < Ntoys; i++) {
        TH2D* BinMigrationtemp = new TH2D("BinMigration", "BinMigration", nphistar, phistarBins, nphistar, phistarBins);
        BinMigrationtemp->Sumw2();
        for (uint j = 0; j < nphistar; j++) {
            for (uint k = 0; k < nphistar; k++) {
                double mean = BinMigration->GetBinContent(j + 1, k + 1);
                if (mean == 0) continue;
                double error = BinMigration->GetBinError(j + 1, k + 1);
                if (mean / error < 5) {
                    x = gRandom->Poisson(mean);
                } else {
                    x = gRandom->Gaus(mean, error);
                }
                BinMigrationtemp->SetBinContent(j + 1, k + 1, x);
                BinMigrationtemp->SetBinError(j + 1, k + 1, error);
            }
        }
        // TH1D* recotemp=BinMigrationtemp->ProjectionX();
        // TH1D* gentemp=BinMigrationtemp->ProjectionY();
        //    RooUnfoldResponse* BinM1temp  =new RooUnfoldResponse (recotemp,gentemp,BinMigrationtemp);
        RooUnfoldResponse* BinM1temp = new RooUnfoldResponse(0, 0, BinMigrationtemp);
        BinM.push_back(BinM1temp);
    }
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
        datatemp->Add(bg_ss, -2.0);
        h_data.push_back(datatemp);
    }
    return h_data;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, vector<TH1D*> mc_truegen) {
    vector<TH1D*> h_eff;
    for (uint i = 0; i < mc_truereco.size(); i++) {
        TH1D* h_efftemp = new TH1D("h_eff", "h_eff", nphistar, phistarBins);
        TH1D* gentemp = (TH1D*) mc_truegen[i]->Clone();
        h_efftemp->Divide(mc_truereco[i], gentemp, 1., 1., "B");
        h_eff.push_back(h_efftemp);
    }
    return h_eff;
}

vector<TH1D*> GetEff(vector<TH1D*> mc_truereco, TH1D * mc_truegen) {
    vector<TH1D*> h_eff;
    for (uint i = 0; i < mc_truereco.size(); i++) {
        TH1D* h_efftemp = new TH1D("h_eff", "h_eff", nphistar, phistarBins);
        TH1D* gentemp = (TH1D*) mc_truegen->Clone();
        h_efftemp->Divide(mc_truereco[i], gentemp, 1., 1., "B");
        h_eff.push_back(h_efftemp);
    }
    return h_eff;
}

TGraphAsymmErrors * ConvertToTGraph(TH1D * h) {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h->GetBinContent(iphistar + 1));
        g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar + 1), h->GetBinError(iphistar + 1));
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
        unfoldBay_data.SetVerbose(RooUnfoldExtraInfo);
        TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
        TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
        h_Unfolded_temp->Divide(h_BinMBay_data, h_eff[e]);

        TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
        g_data.push_back(g_data_temp);

    }
    return g_data;

}

vector<TGraphAsymmErrors *> GetUnfoldedData(vector<RooUnfoldResponse*> BinM, TH1D* h_data, vector<TH1D*> h_eff) {
    vector<TGraphAsymmErrors *> g_data;
    size_t n = BinM.size();
    for (size_t u = 0; u < n; u++) {
        RooUnfoldBayes unfoldBay_data(BinM[u], h_data, 4);
        unfoldBay_data.SetVerbose(RooUnfoldExtraInfo);
        TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
        TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
        h_Unfolded_temp->Divide(h_BinMBay_data, h_eff[u]);

        TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
        g_data.push_back(g_data_temp);
    }
    return g_data;
}

vector<TGraphAsymmErrors *> GetUnfoldedData(RooUnfoldResponse* BinM, TH1D* h_data, TH1D * h_eff) {
    vector<TGraphAsymmErrors *> g_data;
    RooUnfoldBayes unfoldBay_data(BinM, h_data, 4);
    unfoldBay_data.SetVerbose(RooUnfoldExtraInfo);
    TH1D* h_BinMBay_data = (TH1D*) unfoldBay_data.Hreco();
    TH1D* h_Unfolded_temp = (TH1D*) unfoldBay_data.Hreco();
    h_Unfolded_temp->Divide(h_BinMBay_data, h_eff);

    TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_Unfolded_temp);
    g_data.push_back(g_data_temp);

    return g_data;
}

vector<TH1D*> EmptyHVec(int n) {
    vector<TH1D*> h;
    for (int i = 0; i < n; i++) {
        TH1D *phistartemp = new TH1D("phistar", "phistar", nphistar, phistarBins);
        phistartemp->Sumw2();
        h.push_back(phistartemp);
    }
    return h;
}

vector<RooUnfoldResponse*> EmptyBinMVec(int n) {
    vector<RooUnfoldResponse*> BinM;
    for (int i = 0; i < n; i++) {
        TH1D *phistartemp = new TH1D("phistar", "phistar", nphistar, phistarBins);
        phistartemp->Sumw2();
        RooUnfoldResponse* responsetemp = new RooUnfoldResponse(phistartemp, phistartemp);
        BinM.push_back(responsetemp);
    }
    return BinM;
}

/*ask about the blue code to Jeremey*/
void FillRecoEffFluc(int &idx, vector<TH2D*> ToyEff, vector<TH1D*> &h_eff, double weight, double phistar, double E0_pt, double E0_eta, double E1_pt, double E1_eta, bool doBinM, vector<RooUnfoldResponse*> &BinM_eff, double phistar_true = 0, int type = 0) {
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
            h_eff[idx] ->Fill(phistar, new_weight);
        } else {
            h_eff[idx] ->Fill(phistar_true, new_weight);
            BinM_eff[idx]->Fill(phistar, phistar_true, new_weight);
        }
        idx++;
    }
}

void FillTrigEffFluc(int &idx, vector<TH2D*> ToyEffMC, vector<TH2D*> ToyEffData, vector<TH1D*> &h_eff, double weight, double phistar, double E0_pt, double E0_eta, double E1_pt, double E1_eta, bool doBinM, vector<RooUnfoldResponse*> &BinM_eff, double phistar_true = 0) {
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
            h_eff[idx]->Fill(phistar, new_weight);
        } else {
            h_eff[idx] ->Fill(phistar_true, new_weight);
            BinM_eff[idx]->Fill(phistar, phistar_true, new_weight);
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
            //   cout<<"3:"<< E0_pt<<" "<<E0_eta<<" "<<DataNorm1<<endl;
            // cout<<"3:"<< E1_pt<<" "<<E1_eta<<" "<<DataNorm1<<endl;
        }
        if (fabs((weight / new_weight) - 1) > 0.5) cout << weight / new_weight << endl;
        if (!doBinM) {
            h_eff[idx]->Fill(phistar, new_weight);
        } else {
            h_eff[idx] ->Fill(phistar_true, new_weight);
            BinM_eff[idx]->Fill(phistar, phistar_true, new_weight);
        }
        idx++;
    }
}

void GetBGPhiStar(std::string FileName, double sampleweight, vector<TH1D*> &h_eff, vector<TH1D*> &h_fsr_pileup, vector<TH2D*> h_ToyEffSF, vector<TH2D*> h_ToyEffMSF, vector<TH2D*> h_ToyEffTSF, vector<TH2D*> h_ToyEffTMC, vector<TH2D*> h_ToyEffTData) {
    gErrorIgnoreLevel = kError;
    cout << "reading data for " << FileName << endl;
    TChain* t = new TChain(reco_name.c_str(), reco_name.c_str());
    int nfiles; //nfiles set but not used
    nfiles = t->Add(FileName.c_str());

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
        if (!(i % 10000))cout << " Entry number " << i << endl;
        t->GetEntry(i);
        double E0_pt = l_e0_pt ->GetValue();
        double E0_eta = l_e0_eta ->GetValue();
        double E1_pt = l_e1_pt ->GetValue();
        double E1_eta = l_e1_eta ->GetValue();
        double phistar = l_phistar->GetValue();
        double Z_Y = l_YZ->GetValue();
        double weight = sampleweight;
        double weightpu_0 = 0;
        double weightpu_p = 0;
        double weightpu_m = 0;
        if (YSeperation)phistar = OffSetter(phistar, Z_Y);
        for (int w = 0; w < nweights; w++) {

            if (weightid[w] == 1 || weightid[w] == 2 || weightid[w] == 12 || weightid[w] == 13 || weightid[w] == 20 || weightid[w] == 30) {
                weight = weight * weights[w];
            }
            if (weightid[w] == 2) weightpu_0 = weights[w];
            if (weightid[w] == 3) weightpu_p = weights[w];
            if (weightid[w] == 4) weightpu_m = weights[w];
        }
        if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;

        h_eff[0]->Fill(phistar, weight);
        h_fsr_pileup[0]->Fill(phistar, weight * weight_fsr);
        if (weightpu_0 != 0) {
            h_fsr_pileup[1]->Fill(phistar, weight * weightpu_p / weightpu_0);
            h_fsr_pileup[2]->Fill(phistar, weight * weightpu_m / weightpu_0);
        }
        int idx = 1;
        vector<RooUnfoldResponse*> dummy;
        FillRecoEffFluc(idx, h_ToyEffSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy);
        FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy, 0, 1);
        FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy, 0, 2);
        if (idx == 1) cout << "idx did not increase when it should have" << endl;
        FillTrigEffFluc(idx, h_ToyEffTMC, h_ToyEffTData, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 0, dummy);

    }
    cout << "done reading data for " << FileName << endl;
}

TH1D * GetDataPhiStar(double sampleweight, int up = 0) {
    cout << "reading data" << endl;
    TH1D *h_phistar = new TH1D("phistar", "phistar", nphistar, phistarBins);
    h_phistar->Sumw2();
    std::string name = reco_name;
    if (up < 0) name = reco_name_en_l;
    if (up > 0) name = reco_name_en_h;
    TChain* t = new TChain(name.c_str(), name.c_str());
    int nfiles;
    if (up == 0) nfiles = t->Add(File_Data.c_str());
    else if (up > 0) nfiles = t->Add(File_Data_Pt_H.c_str());
    else nfiles = t->Add(File_Data_Pt_L.c_str());
    TBranch *b_reco = t->GetBranch("reco");
    TLeaf *l_phistar = b_reco->GetLeaf("z_phistar_dressed");
    TLeaf *l_ZY = b_reco->GetLeaf("z_y");
    cout << "Entries: " << t->GetEntries() << endl;



    for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
        if (!(i % 10000))cout << " entry2 number " << i << endl;
        t->GetEntry(i);
        double Z_Y = l_ZY->GetValue();
        double phistar = l_phistar->GetValue();
        if (YSeperation)phistar = OffSetter(phistar, Z_Y);

        h_phistar->Fill(phistar, sampleweight);
    }
    cout << "filled data phistar histogram" << endl;
    return h_phistar;
}

TH1D * GetDataPhiStarPt(double sampleweight, int up = 0) {
    cout << "reading data" << endl;
    TH1D *h_phistar = new TH1D("phistar", "phistar", nphistar, phistarBins);
    h_phistar->Sumw2();
    std::string name = reco_name;
    if (up < 0) name = reco_name_en_l;
    //  if (up>0) name=reco_name_en_h;
    TChain* t = new TChain(name.c_str(), name.c_str());
    int nfiles;
    if (up >= 0) nfiles = t->Add(File_Data.c_str());
        //  else if (up > 0) nfiles=t->Add(File_Data_Pt_H.c_str());
    else nfiles = t->Add(File_Data_Pt_L.c_str());
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
    for (int i = 0; i < t->GetEntries()&&(!debug || i < 50000); i++) {
        if (!(i % 10000))cout << " entry number " << i << endl;
        t->GetEntry(i);
        bool E0Tight = l_e0_tight ->GetValue();
        bool E1Tight = l_e1_tight ->GetValue();
        double pt0 = l_e0pt->GetValue();
        double pt1 = l_e1pt->GetValue();
        double eta0 = l_e0eta->GetValue();
        double eta1 = l_e1eta->GetValue();
        double phistar = l_phistar->GetValue();
        double Z_Y = l_ZY->GetValue();
        if (pt0 < ptl || pt1 < ptl) continue;
        if ((pt0 < pth || fabs(eta0) > etah || !E0Tight) && (pt1 < pth || fabs(eta1) > etah || !E1Tight)) continue;

        if (YSeperation)phistar = OffSetter(phistar, Z_Y);
        h_phistar->Fill(phistar, sampleweight);
        idx++;
    }
    cout << "filled data phistar histogram, selected:" << idx << endl;
    return h_phistar;
}

void GetGenPhiStar(double sampleweight, TH1D* &h_phistar, vector<TH1D*> &h_cteq, vector<TH1D*> &h_fsr_pileup) {
    gErrorIgnoreLevel = kError;
    cout << "reading signal gen" << endl;
    h_phistar = new TH1D("phistar", "phistar", nphistar, phistarBins);
    h_phistar->Sumw2();

    TChain* t = new TChain(gen_name.c_str(), gen_name.c_str());
    int nfiles;
    if (doMG) {
        if (elec == 0) nfiles = t->Add(File_Signal_gen.c_str());
        else if (elec == 1) nfiles = t->Add(File_Signal_gen_born.c_str());
        else nfiles = t->Add(File_Signal_gen_bare.c_str());
    } else {
        if (elec == 0) nfiles = t->Add(File_Powheg_gen.c_str());
        else if (elec == 1) nfiles = t->Add(File_Powheg_gen_born.c_str());
        else nfiles = t->Add(File_Powheg_gen_bare.c_str());
    }
    t->SetBranchStatus("event_info", 0); //to disable all branches
    t->SetBranchStatus("reco", 0); //to disable all branches
    TBranch *b_truth = t->GetBranch("truth");
    TLeaf * l_ZY = b_truth->GetLeaf("z_y");
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

    for (int i = 0; i < t->GetEntries()&&(!debug || i < 100000); i++) {
        //  for (int i=0; i<50000;i++){
        t->GetEntry(i);
        double weight = sampleweight;
        double pdfnorm = weights_cteq[0];
        double weightpu_0 = 0;
        double weightpu_p = 0;
        double weightpu_m = 0;
        double phistar = l_phistar->GetValue();
        double Z_Y = l_ZY->GetValue();
        //cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;

        if (YSeperation)phistar = OffSetter(phistar, Z_Y);

        for (int w = 0; w < nweights; w++) {
            if (weightid[w] == 1 || weightid[w] == 2) {
                weight = weight * weights[w];
            }
            //if (weightid[w]==1) {weight=weight*weights[w];}
            if (weightid[w] == 2) weightpu_0 = weights[w];
            if (weightid[w] == 3) weightpu_p = weights[w];
            if (weightid[w] == 4) weightpu_m = weights[w];
        }
        if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;
        h_phistar->Fill(phistar, weight);
        if (pdfnorm != 0) {
            for (int w = 0; w < nwcteq; w++) {
                h_cteq[w] ->Fill(phistar, weight * weights_cteq[w] / pdfnorm);
            }
        }
        h_fsr_pileup[0]->Fill(phistar, weight * weight_fsr);
        if (weightpu_0 != 0) {
            h_fsr_pileup[1]->Fill(phistar, weight * weightpu_p / weightpu_0);
            h_fsr_pileup[2]->Fill(phistar, weight * weightpu_m / weightpu_0);
        }
    }
    cout << "done reading signal gen" << endl;
    return;
}

void GetBinM(double sampleweight, vector<RooUnfoldResponse*> &BinM_eff, vector<RooUnfoldResponse*> &BinM_mcstat, vector<RooUnfoldResponse*> &BinM_cteq, vector<RooUnfoldResponse*> &BinM_fsr_pileup, vector<TH1D*> &h_eff, vector<TH1D*> &h_cteq, vector<TH1D*> &h_fsr_pileup, vector<TH2D*> h_ToyEffSF, vector<TH2D*> h_ToyEffMSF, vector<TH2D*> h_ToyEffTSF, vector<TH2D*> h_ToyEffTMC, vector<TH2D*> h_ToyEffTData) {
    gErrorIgnoreLevel = kError;
    cout << "reading data for " << File_Signal_reco << "  " << reco_name << endl;
    TChain* t = new TChain(reco_name.c_str(), reco_name.c_str());
    int nfiles;
    if (doMG) nfiles = t->Add(File_Signal_reco.c_str());
    else nfiles = t->Add(File_Powheg_reco.c_str());

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

    TH1D *h_rec = new TH1D("phistar", "phistar", nphistar, phistarBins);
    h_rec->Sumw2();
    TH2D* BinMigration = new TH2D("BinMigration", "BinMigration", nphistar, phistarBins, nphistar, phistarBins);
    BinMigration->Sumw2();

    for (int i = 0; i < t->GetEntries()&&(!debug || i < 10000); i++) {
        double Percent = i;
        Percent = 100 * Percent / (double(t->GetEntries()));
        if (!(i % 10000))cout << "percent done? " << Percent << endl;

        //for (int i=0; i<50000;i++){
        t->GetEntry(i);

        bool E0Tight = 1;

        if (doMG) E0Tight = l_e0_tight ->GetValue();
        double E0_pt = l_e0_pt ->GetValue();
        double E0_eta = l_e0_eta ->GetValue();
        double E1_pt = l_e1_pt ->GetValue();
        double E1_eta = l_e1_eta ->GetValue();
        double phistar = l_phistar->GetValue();

        double phistar_true = l_phistar_true->GetValue();
        double Z_Y = l_ZY->GetValue();
        double Z_YTruth = l_YZTruth->GetValue();
        double pdfnorm = 0;

        if (!doMG)pdfnorm = weights_cteq[0];
        double weight = sampleweight;
        double weightpu_0 = 0;
        double weightpu_p = 0;
        double weightpu_m = 0;
        //    cout<<pdfnorm<<" "<<weights_cteq[0]<<endl;
        if (YSeperation)phistar = OffSetter(phistar, Z_Y);
        if (YSeperation)phistar_true = OffSetter(phistar_true, Z_YTruth);


        for (int w = 0; w < nweights; w++) {
            if (weightid[w] == 1 || weightid[w] == 2 || weightid[w] == 12 || weightid[w] == 13 || weightid[w] == 20 || weightid[w] == 30) {
                weight = weight * weights[w];
            }
            if (weightid[w] == 2) weightpu_0 = weights[w];
            if (weightid[w] == 3) weightpu_p = weights[w];
            if (weightid[w] == 4) weightpu_m = weights[w];
        }

        if (weightpu_0 == 0 || weightpu_p == 0 || weightpu_m == 0) cout << "pile-up weights not there" << endl;
        if (weightpu_0 == 0) cout << "pile-up nominal weight 0" << endl;
        //       int ii = 0; //just using it as a counter for a print out statement
        h_eff[0] ->Fill(phistar_true, weight);
        BinM_eff[0] ->Fill(phistar, phistar_true, weight);
        h_rec ->Fill(phistar, weight);
        BinMigration ->Fill(phistar, phistar_true, weight);
        if (pdfnorm != 0) {
            for (int w = 0; w < nwcteq; w++) {
                //cout << w << endl;
                h_cteq[w] ->Fill(phistar_true, weight * weights_cteq[w] / pdfnorm);
                BinM_cteq[w]->Fill(phistar, phistar_true, weight * weights_cteq[w] / pdfnorm);
            }
        }
        h_fsr_pileup[0] ->Fill(phistar_true, weight * weight_fsr);
        BinM_fsr_pileup[0]->Fill(phistar, phistar_true, weight * weight_fsr);
        if (weightpu_0 != 0) {
            h_fsr_pileup[1] ->Fill(phistar_true, weight * weightpu_p / weightpu_0);
            BinM_fsr_pileup[1]->Fill(phistar, phistar_true, weight * weightpu_p / weightpu_0);
            h_fsr_pileup[2] ->Fill(phistar_true, weight * weightpu_m / weightpu_0);
            BinM_fsr_pileup[2]->Fill(phistar, phistar_true, weight * weightpu_m / weightpu_0);
        }
        int idx = 1;
        FillRecoEffFluc(idx, h_ToyEffSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, phistar_true);
        if (E0Tight) {
            FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, phistar_true, 1);
            FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, phistar_true, 2);
        } else {
            FillRecoEffFluc(idx, h_ToyEffMSF, h_eff, weight, phistar, E1_pt, E1_eta, E0_pt, E0_eta, 1, BinM_eff, phistar_true, 1);
            FillRecoEffFluc(idx, h_ToyEffTSF, h_eff, weight, phistar, E1_pt, E1_eta, E0_pt, E0_eta, 1, BinM_eff, phistar_true, 2);
        }
        FillTrigEffFluc(idx, h_ToyEffTMC, h_ToyEffTData, h_eff, weight, phistar, E0_pt, E0_eta, E1_pt, E1_eta, 1, BinM_eff, phistar_true);
    }
    GetToyResponse(BinM_mcstat, BinMigration);
    cout << "done reading data for " << File_Signal_reco << "  " << reco_name << endl;
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
    double signal_weight = 3531.89 / 30459500.; //3504
    double phistarBinsTest[360]; //test purpose, make it more elegent later. 


    if (DoYSeperation) {
        size_t nPhistarsingle = ((sizeof (phistarBinsOne) / sizeof (phistarBinsOne[0])) - 1);

        for (size_t i = 0; i < 10; i++) {
            int index;
            for (size_t j = 0; j <= nPhistarsingle; j++) {
                index = j + nPhistarsingle*i;
                phistarBinsTest[index] = phistarBinsOne[j] + i * (phistarBinsOne[nPhistarsingle] + 1);
                if (j != 0) {
                    if (phistarBinsTest[index] < phistarBinsTest[index - 1])cout << " index number " << index << " number -1 :" << phistarBinsTest[index - 1] << "  index :" << phistarBinsTest[index] << endl;
                }
            }

        }
        phistarBins = phistarBinsTest;
    } else {
        phistarBins = phistarBinsOne;
    }

    if (!doMG) signal_weight = 1966.7 / 3297045.;
    TH1D* Data = GetDataPhiStar(1. / Lumi);
    TH1D* Data_down = GetDataPhiStarPt(1. / Lumi, -1);
    TH1D* Data_up = GetDataPhiStarPt(1. / Lumi, 1);

    // TH1D* Data_down=(TH1D*)Data->Clone();
    // TH1D* Data_up  =(TH1D*)Data->Clone();
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

    GetBGPhiStar(File_tt, ttbar_weight, h_tt_eff, h_tt_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    GetBGPhiStar(File_tautau, tautau_weight, h_tautau_eff, h_tautau_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    GetBGPhiStar(File_tbarw, tbarw_weight, h_tbarw_eff, h_tbarw_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    GetBGPhiStar(File_tw, tw_weight, h_tw_eff, h_tw_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);

    GetBGPhiStar(File_ww, ww_weight, h_ww_eff, h_ww_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    GetBGPhiStar(File_wz, wz_weight, h_wz_eff, h_wz_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    GetBGPhiStar(File_zz, zz_weight, h_zz_eff, h_zz_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);

    TFile f_bg("ratio_data_mc_emu.root");
    TCanvas *c_bg = (TCanvas*) f_bg.Get("Canvas_1");
    TH1F *bg_sf_full = (TH1F*) c_bg->FindObject("hPull");
    bg_sf_full->Sumw2();
    vector<TH1D*> bg_sf = GetToyBg(bg_sf_full);
    TFile f_ss("SS.root");
    TH1D *ss_full = (TH1D*) f_ss.Get("qcd_phistar");
    ss_full->Sumw2();
    TH1D* bg_ss = GetBgSS(ss_full);

    PrintBG(Data, h_tt_eff[0], h_tautau_eff[0], h_tbarw_eff[0], h_tw_eff[0], h_ww_eff[0], h_wz_eff[0], h_zz_eff[0], bg_ss);

    vector<TH1D*> h_data_eff = RemoveBG(Data, bg_sf[0], bg_ss, h_tt_eff, h_tautau_eff, h_tbarw_eff, h_tw_eff, h_ww_eff, h_wz_eff, h_zz_eff);
    vector<TH1D*> h_data_fsr_pileup = RemoveBG(Data, bg_sf[0], bg_ss, h_tt_fsr_pileup, h_tautau_fsr_pileup, h_tbarw_fsr_pileup, h_tw_fsr_pileup, h_ww_fsr_pileup, h_wz_fsr_pileup, h_zz_fsr_pileup);
    vector<TH1D*> h_data_bg;
    for (size_t idx = 0; idx < Ntoys + 5; idx++) {
        TH1D* bgtemp = (TH1D*) h_tt_eff[0]->Clone();
        bgtemp->Add(h_tautau_eff[0], 1.0);
        bgtemp->Add(h_tbarw_eff[0], 1.0);
        bgtemp->Add(h_tw_eff[0], 1.0);
        bgtemp->Add(h_ww_eff[0], 1.0);
        if (idx == Ntoys + 1 || idx == Ntoys + 2 || idx == Ntoys + 3 || idx == Ntoys + 4) {
            bgtemp->Multiply(bg_sf[0]);
        } else {
            bgtemp->Multiply(bg_sf[idx]);
        }
        TH1D* datatemp = (TH1D*) Data->Clone();
        datatemp->Add(bgtemp, -1.0);
        double scale = -1;
        double scaless = -2;
        if (idx == Ntoys + 1) {
            scale = -0.8;
        }
        if (idx == Ntoys + 2) {
            scale = -1.2;
        }
        if (idx == Ntoys + 3) {
            scaless = 0;
        }
        if (idx == Ntoys + 4) {
            scaless = -4.0;
        }
        datatemp->Add(h_wz_eff[0], scale);
        datatemp->Add(h_zz_eff[0], scale);
        datatemp->Add(bg_ss, scaless);
        h_data_bg.push_back(datatemp);
    }
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

    vector<RooUnfoldResponse*> BinM_eff, BinM_mcstat, BinM_cteq, BinM_fsr_pileup;
    vector<TH1D*> mc_truereco_eff, mc_truereco_cteq, mc_truereco_fsr_pileup;
    GetBinM(signal_weight, BinM_eff, BinM_mcstat, BinM_cteq, BinM_fsr_pileup, mc_truereco_eff, mc_truereco_cteq, mc_truereco_fsr_pileup, h_ToyEffSF, h_ToyEffMSF, h_ToyEffTSF, h_ToyEffTMC, h_ToyEffTData);
    TH1D* mc_truegen;
    vector<TH1D*> mc_truegen_cteq, mc_truegen_fsr_pileup;
    GetGenPhiStar(signal_weight, mc_truegen, mc_truegen_cteq, mc_truegen_fsr_pileup);

    vector<TH1D*> h_eff_eff = GetEff(mc_truereco_eff, mc_truegen);
    vector<TH1D*> h_eff_cteq = GetEff(mc_truereco_cteq, mc_truegen_cteq);
    vector<TH1D*> h_eff_fsr_pileup = GetEff(mc_truereco_fsr_pileup, mc_truegen_fsr_pileup);

    TH1D* eff_0 = (TH1D*) h_eff_eff[0]->Clone();
    for (uint idx = 0; idx < nphistar; idx++) {
        eff_0->SetBinError(idx + 1, 0);
    }
    vector<TGraphAsymmErrors *> g_data_phistar_unf = GetUnfoldedData(BinM_eff[0], h_data_eff[0], eff_0);
    vector<TGraphAsymmErrors *> g_data_phistar_eff = GetUnfoldedData(BinM_eff, h_data_eff, h_eff_eff);
    vector<TGraphAsymmErrors *> g_data_phistar_bg = GetUnfoldedData(BinM_eff, h_data_bg, h_eff_eff, 1);
    vector<TGraphAsymmErrors *> g_data_phistar_pt = GetUnfoldedData(BinM_eff, h_data_pt, h_eff_eff, 1);
    vector<TGraphAsymmErrors *> g_data_phistar_cteq = GetUnfoldedData(BinM_cteq, h_data_eff[0], h_eff_cteq);
    vector<TGraphAsymmErrors *> g_data_phistar_fsr_pileup = GetUnfoldedData(BinM_fsr_pileup, h_data_fsr_pileup, h_eff_fsr_pileup);
    vector<TGraphAsymmErrors *> g_data_phistar_mcstat = GetUnfoldedData(BinM_mcstat, h_data_eff, h_eff_eff, 2);

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
    NormalizeGraph(g_data_phistar_unf); //empty
    NormalizeGraph(g_data_phistar_eff);
    NormalizeGraph(g_data_phistar_bg);
    NormalizeGraph(g_data_phistar_pt);
    NormalizeGraph(g_data_phistar_cteq);
    NormalizeGraph(g_data_phistar_fsr_pileup);
    NormalizeGraph(g_data_phistar_mcstat);


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
    //  TGraphAsymmErrors* g_syst_phistar_fsr = CalcTotalSysU_fsr(g_data_phistar_fsr_pileup[0], g_data_phistar_eff[0]); commented out since it was not used any more
    TGraphAsymmErrors* g_syst_phistar_pu = CalcTotalSysU_pileup(g_data_phistar_fsr_pileup[1], g_data_phistar_fsr_pileup[2], g_data_phistar_eff[0]);
    TGraphAsymmErrors* g_syst_phistar_pt = CalcTotalSysU_pileup(g_data_phistar_pt[0], g_data_phistar_pt[1], g_data_phistar_eff[0]);

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
    TH1D* h_data_elec = ConvertToHist(g_data_final_muon, "h_data_elec");
    TH1D* h_data_bgnd = ConvertToHist(g_syst_phistar_bg, "h_data_bgnd");
    TH1D* h_data_fsr = ConvertToHist(g_data_phistar_fsr_pileup[0], "h_data_fsr");
    TH1D* h_data_pup = ConvertToHist(g_data_phistar_fsr_pileup[1], "h_data_pup");
    TH1D* h_data_pum = ConvertToHist(g_data_phistar_fsr_pileup[2], "h_data_pum");



    vector<TGraphAsymmErrors *> g_data_syst;
    g_data_syst.push_back(g_data_phistar_unf[0]);
    g_data_syst.push_back(g_syst_phistar_eff);
    g_data_syst.push_back(g_syst_phistar_mcstat);
    g_data_syst.push_back(g_syst_phistar_bg);
    //g_data_syst.push_back(g_syst_phistar_fsr);
    g_data_syst.push_back(g_syst_phistar_pu);
    g_data_syst.push_back(g_syst_phistar_pt);
    if (!doMG) g_data_syst.push_back(g_syst_phistar_cteq);
    vector<std::string> syst_list;

    syst_list.push_back("unfolding");
    syst_list.push_back("eff");
    syst_list.push_back("mcstat");
    syst_list.push_back("bg");
    //syst_list.push_back("fsr");
    syst_list.push_back("pileup");
    syst_list.push_back("pt");
    if (!doMG) syst_list.push_back("cteq");

    //if (doNorm) printf("%10s : %26s : %10s : %10s : %10s : %10s : %10s : %10s : %10s : %10s \n", "bin", "phistar        ", "total", "unfolding", "MC stat", "pt", "efficiency", "background", "fsr", "pile-up");

    printf("%10s : %26s : %10s : %10s : %10s : %10s : %10s : %10s : %10s : %10s : %10s \n", "bin", "phistar        ", "total", "unfolding", "luminosity", "MC stat", "pt", "efficiency", "background", "fsr", "pile-up");
    //So currently I am just sticking in the part that makes the covariant matrix here, right after everything is normalized.

    cout << " test 1" << endl;
    std::string YSeperationname = "YSeperated";
    for (size_t i = 0; i < nYSeper; i++) {
        char buffer [10];
        sprintf(buffer, "-%.2f", YSeperation[i]);
        YSeperationname += buffer;
    }


    std::string textn = "CovarianceMatrix";

    textn += Tag;
    // if (doNorm) textn+="Norm_";
    textn += "Abs_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (DoYSeperation) textn += YSeperationname;
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    TFile tr7(textn.c_str(), "RECREATE");
    TMatrixD unf_Covariant_Matrix(nphistar, nphistar);
    if (g_data_phistar_unf.size() > 10) {
        TH2D *unf_Covariant_Matrix_Histo = new TH2D("unf_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* unf_Correl_hist = new TH2D("unf_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        Cov_Histo_Creator(g_data_phistar_unf, unf_Covariant_Matrix_Histo, unf_Correl_hist); //empty
        MakeCovTMatrix(unf_Covariant_Matrix_Histo, unf_Covariant_Matrix);
        unf_Covariant_Matrix_Histo->Write();
    } else {
        MakeCovTMatrix(g_data_phistar_unf[0], unf_Covariant_Matrix);
    }
    cout << " test 2" << endl;



    vector<TGraphAsymmErrors *> g_syst_phistar_reff_Vector;
    vector<TGraphAsymmErrors *> g_syst_phistar_meff_Vector;
    vector<TGraphAsymmErrors *> g_syst_phistar_teff_Vector;
    vector<TGraphAsymmErrors *> g_syst_phistar_teff_m_Vector;
    vector<TGraphAsymmErrors *> g_syst_phistar_teff_d_Vector;
    g_syst_phistar_reff_Vector.push_back(g_data_phistar_eff[0]);
    g_syst_phistar_meff_Vector.push_back(g_data_phistar_eff[0]);
    g_syst_phistar_teff_Vector.push_back(g_data_phistar_eff[0]);
    g_syst_phistar_teff_m_Vector.push_back(g_data_phistar_eff[0]);
    g_syst_phistar_teff_d_Vector.push_back(g_data_phistar_eff[0]); //probably should make a function relook later
    cout << " test 3" << endl;

    int nt = Ntoys2;
    double part = 1;
    uint gs = 1 + (nt * (part - 1));
    for (size_t i = gs; i < 101; i++) {
        g_syst_phistar_reff_Vector.push_back(g_data_phistar_eff[gs]);
        g_syst_phistar_meff_Vector.push_back(g_data_phistar_eff[gs + 100]);
        g_syst_phistar_teff_Vector.push_back(g_data_phistar_eff[gs + 200]);
        g_syst_phistar_teff_m_Vector.push_back(g_data_phistar_eff[gs + 300]);
        g_syst_phistar_teff_d_Vector.push_back(g_data_phistar_eff[gs + 400]);
    }


    TMatrixD reff_Covariant_Matrix(nphistar, nphistar);
    cout << " test 4" << endl;
    if (g_syst_phistar_reff_Vector.size() > 10) {
        TH2D *reff_Covariant_Matrix_Histo = new TH2D("reff_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* reff_Correl_hist = new TH2D("reff_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);

        Cov_Histo_Creator(g_syst_phistar_reff_Vector, reff_Covariant_Matrix_Histo, reff_Correl_hist);
        MakeCovTMatrix(reff_Covariant_Matrix_Histo, reff_Covariant_Matrix);
        reff_Covariant_Matrix_Histo->Write();
        reff_Correl_hist->Write();
    } else {
        MakeCovTMatrix(g_syst_phistar_eff, reff_Covariant_Matrix);
    }
    cout << " test 4.1" << endl;

    TMatrixD eff_Covariant_Matrix(reff_Covariant_Matrix);
    TMatrixD meff_Covariant_Matrix(nphistar, nphistar);
    cout << " test 5" << endl;
    if (g_syst_phistar_meff_Vector.size() > 10) {
        TH2D *meff_Covariant_Matrix_Histo = new TH2D("meff_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* meff_Correl_hist = new TH2D("meff_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        cout << " test 6" << endl;
        Cov_Histo_Creator(g_syst_phistar_meff_Vector, meff_Covariant_Matrix_Histo, meff_Correl_hist);
        cout << " test 7" << endl;
        MakeCovTMatrix(meff_Covariant_Matrix_Histo, meff_Covariant_Matrix);

        meff_Covariant_Matrix_Histo->Write();
        meff_Correl_hist->Write();
    } else {
        MakeCovTMatrix(g_syst_phistar_eff, meff_Covariant_Matrix);
    }
    eff_Covariant_Matrix = eff_Covariant_Matrix + meff_Covariant_Matrix;
    TMatrixD teff_m_Covariant_Matrix(nphistar, nphistar);

    if (g_syst_phistar_teff_m_Vector.size() > 10) {
        TH2D *teff_m_Covariant_Matrix_Histo = new TH2D("teff_m_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* teff_m_Correl_hist = new TH2D("teff_m_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        Cov_Histo_Creator(g_syst_phistar_teff_m_Vector, teff_m_Covariant_Matrix_Histo, teff_m_Correl_hist);
        MakeCovTMatrix(teff_m_Covariant_Matrix_Histo, teff_m_Covariant_Matrix);
        teff_m_Covariant_Matrix_Histo->Write();
        teff_m_Correl_hist->Write();
    } else {
        MakeCovTMatrix(g_syst_phistar_eff, teff_m_Covariant_Matrix);
    }
    eff_Covariant_Matrix = eff_Covariant_Matrix + teff_m_Covariant_Matrix;

    TMatrixD teff_d_Covariant_Matrix(nphistar, nphistar);
    if (g_syst_phistar_teff_d_Vector.size() > 10) {
        TH2D *teff_d_Covariant_Matrix_Histo = new TH2D("teff_d_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* teff_d_Correl_hist = new TH2D("teff_d_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        Cov_Histo_Creator(g_syst_phistar_teff_d_Vector, teff_d_Covariant_Matrix_Histo, teff_d_Correl_hist);
        MakeCovTMatrix(teff_d_Covariant_Matrix_Histo, teff_d_Covariant_Matrix);
        teff_d_Covariant_Matrix_Histo->Write();
        teff_d_Correl_hist->Write();
    } else {
        MakeCovTMatrix(g_syst_phistar_eff, teff_d_Covariant_Matrix);
    }

    eff_Covariant_Matrix = eff_Covariant_Matrix + teff_d_Covariant_Matrix;

    TMatrixD mcstat_Covariant_Matrix(nphistar, nphistar);

    if (g_data_phistar_mcstat.size() > 10) {
        TH2D *mcstat_Covariant_Matrix_Histo = new TH2D("mcstat_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* mcstat_Correl_hist = new TH2D("mcstat_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        Cov_Histo_Creator(g_data_phistar_mcstat, mcstat_Covariant_Matrix_Histo, mcstat_Correl_hist);
        MakeCovTMatrix(mcstat_Covariant_Matrix_Histo, mcstat_Covariant_Matrix);
        mcstat_Covariant_Matrix_Histo->Write();
        mcstat_Correl_hist->Write();
    } else {
        MakeCovTMatrix(g_syst_phistar_mcstat, mcstat_Covariant_Matrix);
    }

    TMatrixD bg_Covariant_Matrix(nphistar, nphistar);
    if (g_data_phistar_bg.size() > 10) {
        TH2D *bg_Covariant_Matrix_Histo = new TH2D("bg_Covariant_Matrix_Histo", "Covar_hist", nphistar, phistarBins, nphistar, phistarBins);
        TH2D* bg_Correl_hist = new TH2D("bg_Correl_hist", "Correlation", nphistar, phistarBins, nphistar, phistarBins);
        MatrixName = "BGMatrix";
        double x, y;
        g_data_phistar_bg[0]->GetPoint(0, x, y);
        cout << " and our value is :" << y << endl;
        Cov_Histo_Creator(g_data_phistar_bg, bg_Covariant_Matrix_Histo, bg_Correl_hist);
        MakeCovTMatrix(bg_Covariant_Matrix_Histo, bg_Covariant_Matrix);
        bg_Covariant_Matrix_Histo->Write();
        MatrixName = "";
    } else {
        MakeCovTMatrix(g_syst_phistar_bg, bg_Covariant_Matrix);
    }

    TMatrixD Pilup_Matrix(nphistar, nphistar);
    MakeCovTMatrix(g_syst_phistar_pu, Pilup_Matrix);

    TMatrixD PtMatrix(nphistar, nphistar);
    MakeCovTMatrix(g_syst_phistar_pt, PtMatrix);

    TMatrixD LumMatrix(nphistar, nphistar);
    double LumiError = .026; // error in luminosity
    MakeCovTMatrix(g_data_phistar_eff[0], LumMatrix, LumiError); //if I understand everything the first element in most of the systems should be a the orginal so I just chose this one

    TMatrixD h_eff_effMatrix(nphistar, nphistar);
    MakeCovTMatrix(h_eff_eff[0], h_eff_effMatrix);

    unf_Covariant_Matrix.Write("Unfolded");
    eff_Covariant_Matrix.Write("Efficency");
    h_eff_effMatrix.Write("MCStatBinbyBin");
    mcstat_Covariant_Matrix.Write("MCstat");
    bg_Covariant_Matrix.Write("BG");
    Pilup_Matrix.Write("PileUp");
    LumMatrix.Write("Lumi");
    PtMatrix.Write("Pt");
    //Cov_Histo_Creator(g_data_phistar_cteq, cteq_Covariant_Matrix_Histo); //empty
    //Cov_Histo_Creator(g_data_phistar_mcstat, mcstat_Covariant_Matrix_Histo);



    TGraphAsymmErrors* g_data_final = GetDataFinal(g_data_syst, syst_list);

    textn = "Final_Hist_";
    textn += Tag;
    // if (doNorm) textn+="Norm_";
    textn += "Abs_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (DoYSeperation) textn += YSeperationname;
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    TFile tr(textn.c_str(), "RECREATE");
    h_data_elec->Write();
    h_data_bgnd->Write();
    h_data_fsr->Write();
    h_data_pup->Write();
    h_data_pum->Write();

    textn = "Data_Graph_";
    textn += Tag;
    // if (doNorm) textn+="Norm_";
    textn += "Abs_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (DoYSeperation) textn += YSeperationname;
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    TFile tr2(textn.c_str(), "RECREATE");
    g_data_final->Write();

    for (int i = 0; i < ((int) nphistar) && debug; i++) {
        double x, y;
        g_data_final->GetPoint(i, x, y);
        cout << " and last we have :" << y << endl;
    }

    //Now for the normalised distribution:

    NormalizeGraph(g_data_norm_unf, 1);
    NormalizeGraph(g_data_norm_eff, 1);
    NormalizeGraph(g_data_norm_bg, 1);
    NormalizeGraph(g_data_norm_pt, 1);
    NormalizeGraph(g_data_norm_cteq, 1);
    NormalizeGraph(g_data_norm_fsr_pileup, 1);
    NormalizeGraph(g_data_norm_mcstat, 1);

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
    // TGraphAsymmErrors* g_syst_norm_fsr = CalcTotalSysU_fsr(g_data_norm_fsr_pileup[0], g_data_norm_eff[0]); not used so I commented it out
    TGraphAsymmErrors* g_syst_norm_pu = CalcTotalSysU_pileup(g_data_norm_fsr_pileup[1], g_data_norm_fsr_pileup[2], g_data_norm_eff[0]);
    TGraphAsymmErrors* g_syst_norm_pt = CalcTotalSysU_pileup(g_data_norm_pt[0], g_data_norm_pt[1], g_data_norm_eff[0]);

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
    TH1D* h_data_norm_fsr = ConvertToHist(g_data_norm_fsr_pileup[0], "h_data_fsr");
    TH1D* h_data_norm_pup = ConvertToHist(g_data_norm_fsr_pileup[1], "h_data_pup");
    TH1D* h_data_norm_pum = ConvertToHist(g_data_norm_fsr_pileup[2], "h_data_pum");

    vector<TGraphAsymmErrors *> g_data_syst_norm;
    g_data_syst_norm.push_back(g_data_norm_unf[0]);
    g_data_syst_norm.push_back(g_syst_norm_eff);
    g_data_syst_norm.push_back(g_syst_norm_mcstat);
    g_data_syst_norm.push_back(g_syst_norm_bg);
    //g_data_syst_norm.push_back(g_syst_norm_fsr);
    g_data_syst_norm.push_back(g_syst_norm_pu);
    g_data_syst_norm.push_back(g_syst_norm_pt);
    if (!doMG) g_data_syst_norm.push_back(g_syst_norm_cteq);
    vector<std::string> syst_list_norm;
    syst_list_norm.push_back("unfolding");
    syst_list_norm.push_back("eff");
    syst_list_norm.push_back("mcstat");
    syst_list_norm.push_back("bg");
    //syst_list_norm.push_back("fsr");
    syst_list_norm.push_back("pileup");
    syst_list_norm.push_back("pt");
    if (!doMG) syst_list_norm.push_back("cteq");
    printf("%10s : %26s : %10s : %10s : %10s : %10s : %10s : %10s : %10s : %10s \n", "bin", "phistar        ", "total", "unfolding", "MC stat", "pt", "efficiency", "background", "fsr", "pile-up");

    TGraphAsymmErrors* g_data_final_norm = GetDataFinal(g_data_syst_norm, syst_list_norm, 1);

    textn = "Final_Hist_";
    textn += Tag;
    textn += "Norm_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (DoYSeperation) textn += YSeperationname;
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    TFile tr3(textn.c_str(), "RECREATE");
    h_data_norm_elec->Write();
    h_data_norm_bgnd->Write();
    h_data_norm_fsr->Write();
    h_data_norm_pup->Write();
    h_data_norm_pum->Write();


    textn = "Data_Graph_";
    textn += Tag;
    textn += "Norm_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (DoYSeperation) textn += YSeperationname;
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    TFile tr4(textn.c_str(), "RECREATE");
    cout << "file name " << textn << endl;
    g_data_final_norm->Write();

}
