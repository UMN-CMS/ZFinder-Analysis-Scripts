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

const bool doPHcomp = 0;//Done for andrews to compare differnt tunes?
const bool doNorm = 0;
const int elec = 1;
const int doMG = 0;
const std::string Tag = "";
const std::string Type = "combined"; //elec, muon or combined

const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;

double ElectAbsTTBar[] = {0.001574, 0.001257, 0.001685, 0.001581, 0.001027, 0.001283, 0.001130, 0.001953, 0.001491, 0.002356, 0.002693, 0.002553, 0.002389, 0.002398, 0.004139, 0.003285, 0.002939, 0.005734, 0.005421, 0.005708, 0.008338, 0.009034, 0.010407, 0.015298, 0.019190, 0.030277, 0.042358, 0.079723, 0.128774, 0.169056, 0.201204, 0.190319, 0.186362, 0.190599, 0.181289, 0.181289};
double ElectNormTTBar[] = {0.012417, 0.012734, 0.012307, 0.012411, 0.012965, 0.012709, 0.012861, 0.012038, 0.012501, 0.011636, 0.011299, 0.011439, 0.011602, 0.011593, 0.009853, 0.010707, 0.011053, 0.008257, 0.008570, 0.008283, 0.005652, 0.004957, 0.003583, 0.001308, 0.005201, 0.016289, 0.028372, 0.065743, 0.114800, 0.155088, 0.187240, 0.176354, 0.172396, 0.176633, 0.167323, 0.167323};
double MuonAbsTTBar[] = {0.00151815, 0.00150637, 0.0015723, 0.00159348, 0.00170658, 0.00162156, 0.00183434, 0.00192313, 0.00208938, 0.00222283, 0.00240952, 0.00267002, 0.00289458, 0.00313999, 0.00345079, 0.00397061, 0.00456066, 0.00533871, 0.00616368, 0.0072238, 0.00836111, 0.010165, 0.0127692, 0.0164479, 0.0225485, 0.0335744, 0.055569, 0.098426, 0.155315, 0.209841, 0.236886, 0.23359, 0.217536, 0.211622, 0.226457};
double MuonNormTTBar[] = {0.000273672, 0.000259274, 0.000273625, 0.00027254, 0.000298313, 0.000267256, 0.000314293, 0.00033616, 0.000356221, 0.000380016, 0.000428403, 0.000457978, 0.000502123, 0.000552694, 0.000597325, 0.000694674, 0.000796269, 0.000933197, 0.00106172, 0.00126091, 0.00145946, 0.00177974, 0.00222539, 0.00287601, 0.00393807, 0.00585608, 0.00968358, 0.0170782, 0.0268222, 0.0360863, 0.0406359, 0.0400708, 0.0373819, 0.0363764, 0.0388815};

void NewError(TGraphAsymmErrors* g_data_final, double* array) {

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, yerror, xerror;
        g_data_final->GetPoint(iphistar, x, y);
        double NewTTError = y * array[iphistar] / 100;
        double OldError = g_data_final->GetErrorYhigh(iphistar);
        double newError = sqrt(OldError * OldError + NewTTError * NewTTError);
        g_data_final->SetPointError(iphistar, 0, 0, newError, newError);
    }

}

void PrintValue(ofstream &outputfile, int n, int p, double value, double error) {
    if (p > 1 || p<-1) {
        outputfile << std::setprecision(n + p) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
    } else {
        outputfile << std::setprecision(n) << value << "&" << error;
    }
}

void PrintFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, bool dataonly = 0, TGraphAsymmErrors* g_re_final = 0) {
    ofstream outputfile;
    std::string textname = "Table_Values_";
    textname += Tag;
    if (!dataonly) {
        textname += Type;
        textname += "_";
    }
    if (dataonly && Type == "combined")textname += "MuEl";
    if (dataonly && Type == "elec")textname += "PHMG";
    if (Type == "elec" && !doMG) textname += "PH_";
    if (doNorm) textname += "Norm_";
    else textname += "Abs_";
    if (elec == 0)textname += "Dressed.txt";
    if (elec == 1)textname += "Born.txt";
    if (elec == 2)textname += "Naked.txt";
    outputfile.open(textname.c_str());
    std::string tableheader = "";
    if (!dataonly) {
        if (elec == 1) {
            if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data (pb)} & \\multicolumn{2}{c|}{MadGraph (pb)} & \\multicolumn{2}{c|}{Powheg (pb)} & \\multicolumn{2}{c|}{Resbos (pb)}\\\\ \\hline";
            else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data} & \\multicolumn{2}{c|}{MadGraph} & \\multicolumn{2}{c|}{Powheg} & \\multicolumn{2}{c|}{Resbos}\\\\ \\hline";
        } else {
            if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data (pb)} & \\multicolumn{2}{c|}{MadGraph (pb)} & \\multicolumn{2}{c|}{Powheg (pb)}\\\\ \\hline";
            else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data} & \\multicolumn{2}{c|}{MadGraph} & \\multicolumn{2}{c|}{Powheg}\\\\ \\hline";
        }
    } else {
        if (Type == "combined") {
            if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data: Z \\rightarrow ee (pb)} & \\multicolumn{2}{c|}{Data: Z \\rightarrow \\mu\\mu (pb)} & \\multicolumn{2}{c|}{Data: Z \\rightarrow ll (pb)} \\\\ \\hline";
            else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data: Z \\rightarrow ee} & \\multicolumn{2}{c|}{Data: Z \\rightarrow \\mu\\mu} & \\multicolumn{2}{c|}{Data: Z \\rightarrow ll} \\\\ \\hline";
        }
        if (Type == "elec") {
            if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data unfolded with Powheg (pb)} & \\multicolumn{2}{c|}{Data unfolded with MadGraph (pb)} \\\\ \\hline";
            else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data unfolded with Powheg} & \\multicolumn{2}{c|}{Data unfolded with MadGraph}\\\\ \\hline";
        }
    }
    outputfile << tableheader << "\n";
    std::cout << tableheader << endl;
    for (size_t i = 0; i < nphistar; i++) {
        double x, y, xmg, ymg, xph, yph, xre, yre, temp_r, temp2_r;
        g_data_final->GetPoint(i, x, y);
        g_mg_final->GetPoint(i, xmg, ymg);
        g_ph_final->GetPoint(i, xph, yph);
        double temp_d = g_data_final->GetErrorYhigh(i);
        double temp_m = g_mg_final->GetErrorYhigh(i);
        double temp_p = g_ph_final->GetErrorYhigh(i);
        double temp2_d = y;
        double temp2_m = ymg;
        double temp2_p = yph;
        int n_d = 0;
        int n_m = 0;
        int n_p = 0;
        int n_r = 0;
        int p_d = 0;
        int p_m = 0;
        int p_p = 0;
        int p_r = 0;
        while (temp_d < 1) {
            temp_d = temp_d * 10.;
            n_d++;
        }
        while (temp_m < 1) {
            temp_m = temp_m * 10.;
            n_m++;
        }
        while (temp_p < 1) {
            temp_p = temp_p * 10.;
            n_p++;
        }
        while (temp_d > 10) {
            temp_d = temp_d / 10.;
            n_d = n_d - 1;
        }
        while (temp_m > 10) {
            temp_m = temp_m / 10.;
            n_m = n_m - 1;
        }
        while (temp_p > 10) {
            temp_p = temp_p / 10.;
            n_p = n_p - 1;
        }
        if (!dataonly && elec == 1) {
            g_re_final->GetPoint(i, xre, yre);
            temp_r = g_re_final->GetErrorYhigh(i);
            while (temp_r < 1) {
                temp_r = temp_r * 10.;
                n_r++;
            }
            while (temp_r > 10) {
                temp_r = temp_r / 10.;
                n_r = n_r - 1;
            }
        }
        while (temp2_d > 10) {
            temp2_d = temp2_d / 10.;
            p_d++;
        }
        while (temp2_m > 10) {
            temp2_m = temp2_m / 10.;
            p_m++;
        }
        while (temp2_p > 10) {
            temp2_p = temp2_p / 10.;
            p_p++;
        }
        while (temp2_d < 1) {
            temp2_d = temp2_d * 10.;
            p_d = p_d - 1;
        }
        while (temp2_m < 1) {
            temp2_m = temp2_m * 10.;
            p_m = p_m - 1;
        }
        while (temp2_p < 1) {
            temp2_p = temp2_p * 10.;
            p_p = p_p - 1;
        }
        if (!dataonly && elec == 1) {
            g_re_final->GetPoint(i, xre, yre);
            temp_r = g_re_final->GetErrorYhigh(i);
            temp2_r = yre;
            while (temp_r > 10) {
                temp_r = temp_r / 10.;
                p_r++;
            }
            while (temp2_r < 1) {
                temp2_r = temp2_r * 10.;
                p_r = p_r - 1;
            }
        }

        outputfile << std::fixed << std::setprecision(3) << phistarBins[i] << "-" << phistarBins[i + 1] << " & ";
        if (!(dataonly && Type == "combined")) {
            if (!(Type == "elec") || !(dataonly)) {
                PrintValue(outputfile, n_d + 1, p_d, y, g_data_final->GetErrorYhigh(i));
                outputfile << " & ";
            }
            PrintValue(outputfile, n_m + 1, p_m, ymg, g_mg_final->GetErrorYhigh(i));
            outputfile << " & ";
            PrintValue(outputfile, n_p + 1, p_p, yph, g_ph_final->GetErrorYhigh(i));
            if (!dataonly && elec == 1) {
                outputfile << " & ";
                PrintValue(outputfile, n_r + 1, p_r, yre, g_re_final->GetErrorYhigh(i));
            }
        } else {
            PrintValue(outputfile, n_m + 1, p_m, ymg, g_mg_final->GetErrorYhigh(i));
            outputfile << " & ";
            PrintValue(outputfile, n_p + 1, p_p, yph, g_ph_final->GetErrorYhigh(i));
            outputfile << " & ";
            PrintValue(outputfile, n_d + 1, p_d, y, g_data_final->GetErrorYhigh(i));
        }
        outputfile << " \\\\ \\hline" << "\n";
    }
    outputfile.close();
}

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_dummy_phistar, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2 = 0, TGraphAsymmErrors* g_re_final = 0, TGraphAsymmErrors* g_ratio_re_phistar = 0) {

    cout<<"test 1"<<endl;
    TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot", 800, 900);
    FinalPhiTot->Divide(1, 2);
    FinalPhiTot->cd(1);
    gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    g_dummy_phistar->GetXaxis()->SetRangeUser(0.001, 3.0);
    if (doNorm) g_dummy_phistar->GetYaxis()->SetRangeUser(0.0012, 200.0);
    else g_dummy_phistar->GetYaxis()->SetRangeUser(0.8, 70000.0);
    g_dummy_phistar->GetXaxis()->CenterTitle();
    g_dummy_phistar->GetYaxis()->CenterTitle();
    g_dummy_phistar->Draw("A2");
    cout<<"test 2"<<endl;
    g_mg_final->SetMarkerColor(kBlue - 7);
    g_mg_final->SetLineColor(kBlue - 7);
    g_mg_final->SetMarkerSize(1);
    g_mg_final->SetLineWidth(2);
    g_mg_final->SetMarkerStyle(21);
    g_mg_final->Draw("PEsame");
    cout<<"test 3"<<endl;
    g_ph_final->SetMarkerColor(kRed);
    g_ph_final->SetLineColor(kRed);
    g_ph_final->SetMarkerSize(1);
    g_ph_final->SetLineWidth(2);
    g_ph_final->SetMarkerStyle(22);
    g_ph_final->Draw("PEsame");
    cout<<"test 4"<<endl;
    if (!isPlot2 && elec == 1) {
        cout<<"test 5"<<endl;
        g_re_final->SetMarkerColor(kGreen + 1);
        g_re_final->SetLineColor(kGreen + 1);
        g_re_final->SetMarkerSize(1);
        g_re_final->SetLineWidth(2);
        g_re_final->SetMarkerStyle(23);
        
        g_re_final->Draw("PEsame");
        
    }
    cout<<"test 6"<<endl;
    g_data_final->SetFillColor(kYellow);
    cout<<"test 7"<<endl;
    g_data_final->SetMarkerSize(1);
    g_data_final->SetLineWidth(2);
    g_data_final->SetMarkerStyle(20);
    g_data_final->Draw("PEsame");
    g_data_final->SetFillColor(kYellow);

    TLegend* leg;
    if (isPlot2) leg = new TLegend(0.15, 0.06, 0.80, 0.27);
    else leg = new TLegend(0.15, 0.06, 0.80, 0.31); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);
cout<<"test 8"<<endl;
    if (!isPlot2) {
        leg->AddEntry(g_data_final, "2012 data", "PEF");
        if (Type == "elec" && doPHcomp) {
            leg->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ee (Powheg+Pythia8 4C*)", "P");
            leg->AddEntry(g_ph_final, "#gamma*/Z #rightarrow ee (Powheg+Pythia6 Z2*)", "P");
            if (elec == 1) leg->AddEntry(g_re_final, "#gamma*/Z #rightarrow ee (Powheg+Pythia8 Andrew)", "P");
        }
        if (Type == "elec" && !doPHcomp) {
            leg->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ee (MadGraph+Pythia6 Z2*)", "P");
            leg->AddEntry(g_ph_final, "#gamma*/Z #rightarrow ee (Powheg+Pythia6 Z2*)", "P");
            if (elec == 1) leg->AddEntry(g_re_final, "#gamma*/Z #rightarrow ee (Resbos)", "P");
        }
        if (Type == "muon") {
            leg->AddEntry(g_mg_final, "#gamma*/Z #rightarrow #mu#mu (MadGraph+Pythia6 Z2*)", "P");
            leg->AddEntry(g_ph_final, "#gamma*/Z #rightarrow #mu#mu (Powheg+Pythia6 Z2*)", "P");
            if (elec == 1) leg->AddEntry(g_re_final, "#gamma*/Z #rightarrow #mu#mu (Resbos)", "P");
        }
        cout<<"test 4"<<endl;
        if (Type == "combined") {
            leg->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "P");
            leg->AddEntry(g_ph_final, "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "P");
            if (elec == 1) leg->AddEntry(g_re_final, "#gamma*/Z #rightarrow ll (Resbos)", "P");
        }
    } else {
        if (Type == "combined") {
            leg->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ee (2012 data)", "P");
            leg->AddEntry(g_ph_final, "#gamma*/Z #rightarrow #mu#mu (2012 data)", "P");
            if (doMG)leg->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
            else leg->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
        }
        if (Type == "elec") {
            leg->AddEntry(g_mg_final, "data (unfolded with Powheg)", "P");
            leg->AddEntry(g_ph_final, "data (unfolded with MadGraph)", "P");
            if (doMG)leg->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
            else leg->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
        }
    }
cout<<"test 9"<<endl;
    leg->Draw();

    TLatex mark3;
    mark3.SetTextSize(0.043);
    mark3.SetTextFont(42);
    mark3.SetNDC(true);
    mark3.DrawLatex(0.71, 0.955, "19.7 fb^{-1} (8 TeV)");
    TLatex mark1;
    mark1.SetTextSize(0.05);
    mark1.SetTextFont(61);
    mark1.SetNDC(true);
    mark1.DrawLatex(0.15, 0.955, "CMS");
    TLatex mark2;
    mark2.SetTextSize(0.038);
    mark2.SetTextFont(52);
    mark2.SetNDC(true);
    mark2.DrawLatex(0.23, 0.955, " Preliminary");
    TLatex mark;
    mark.SetTextSize(0.043);
    mark.SetTextFont(42);
    //mark.SetTextSize(0.035);
    mark.SetNDC(true);
    // mark.DrawLatex(0.70,0.95,"19.7 fb^{-1} (8 TeV)");
    // mark.DrawLatex(0.19,0.88,"CMS Preliminary");
    if (Type == "elec" && !isPlot2) {
        mark.DrawLatex(0.53, 0.88, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
        mark.DrawLatex(0.53, 0.81, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
        mark.DrawLatex(0.53, 0.74, "60 GeV < M_{ee} < 120 GeV");
    }
    if (Type == "muon" && !isPlot2) {
        mark.DrawLatex(0.19, 0.20, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.13, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.06, "60 GeV < M_{#mu#mu} < 120 GeV");
    }
    if (Type == "combined" || isPlot2) {
        mark.DrawLatex(0.53, 0.88, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
        mark.DrawLatex(0.53, 0.81, "p_{T}^{l_{0}} > 30 GeV,  p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.53, 0.74, "60 GeV < M_{ll} < 120 GeV");
    }
    FinalPhiTot->cd(2);
    gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
    gPad->SetBottomMargin(0.37);
    gPad->SetTopMargin(0.01);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);

    if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC  ");
    else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data   ");
    g_ratio_phistar->GetXaxis()->SetRangeUser(0.001, 3.0);
    g_ratio_phistar->GetYaxis()->SetRangeUser(0.76, 1.24);
    if (isPlot2) g_ratio_phistar->GetYaxis()->SetRangeUser(0.75, 1.25);
    if (isPlot2 && (!doNorm)) g_ratio_phistar->GetYaxis()->SetRangeUser(0.7, 1.3);
    if (!isPlot2 && Type == "elec" && doPHcomp && !doNorm) g_ratio_phistar->GetYaxis()->SetRangeUser(0.3, 1.7);
    g_ratio_phistar->GetYaxis()->SetTitleOffset(0.45);
    g_ratio_phistar->SetFillColor(kYellow);
    g_ratio_phistar->GetXaxis()->SetTitleSize(0.15);
    g_ratio_phistar->GetXaxis()->CenterTitle();
    g_ratio_phistar->GetYaxis()->CenterTitle();
    g_ratio_phistar->Draw("AE2");
    g_ratio_mg_phistar->SetMarkerSize(1);
    g_ratio_mg_phistar->SetLineWidth(2);
    g_ratio_mg_phistar->SetMarkerStyle(21);
    g_ratio_mg_phistar->SetMarkerColor(kBlue - 7);
    g_ratio_mg_phistar->SetLineColor(kBlue - 7);
    g_ratio_mg_phistar->Draw("PEsame");
    g_ratio_ph_phistar->SetMarkerSize(1);
    g_ratio_ph_phistar->SetLineWidth(2);
    g_ratio_ph_phistar->SetMarkerStyle(22);
    g_ratio_ph_phistar->SetMarkerColor(kRed);
    g_ratio_ph_phistar->SetLineColor(kRed);
    g_ratio_ph_phistar->Draw("PEsame");
    if (!isPlot2 && elec == 1) {
        g_ratio_re_phistar->SetMarkerSize(1);
        g_ratio_re_phistar->SetLineWidth(2);
        g_ratio_re_phistar->SetMarkerStyle(23);
        g_ratio_re_phistar->SetMarkerColor(kGreen + 1);
        g_ratio_re_phistar->SetLineColor(kGreen + 1);
        g_ratio_re_phistar->Draw("PEsame");
    }

    std::string plotname = "ZShape_";
    plotname += Tag;
    if (doPHcomp)plotname += "Tunes_";
    if (isPlot2 && Type == "combined")plotname += "MuEl";
    if (isPlot2 && Type == "elec")plotname += "PHMG";
    else plotname += Type;
    plotname += "_";
    if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    FinalPhiTot->SaveAs((plotname + "pdf").c_str());
    FinalPhiTot->SaveAs((plotname + "png").c_str());
    FinalPhiTot->SaveAs((plotname + "C").c_str());

    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->cd();
    FinalPhiRatio->SetLogx();
    if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC");
    else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data");
    g_ratio_phistar->SetFillColor(kYellow);
    g_ratio_phistar->GetYaxis()->SetRangeUser(0.7, 1.3);
    g_ratio_phistar->GetYaxis()->SetTitleOffset(1.2);
    g_ratio_phistar->GetYaxis()->SetTitleSize(0.04);
    g_ratio_phistar->GetYaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetTitleSize(0.04);
    g_ratio_phistar->GetXaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetLabelOffset(-0.01);
    g_ratio_phistar->Draw("AE2");
    g_ratio_mg_phistar->Draw("PEsame");
    g_ratio_ph_phistar->Draw("PEsame");
    if (!isPlot2 && elec == 1) {
        g_ratio_re_phistar->Draw("PEsame");
    }
    // TLatex mark3;
    mark3.SetTextSize(0.043);
    mark3.SetTextFont(42);
    mark3.SetNDC(true);
    mark3.DrawLatex(0.71, 0.955, "19.7 fb^{-1} (8 TeV)");
    // TLatex mark1;
    mark1.SetTextSize(0.05);
    mark1.SetTextFont(61);
    mark1.SetNDC(true);
    mark1.DrawLatex(0.15, 0.955, "CMS");
    // TLatex mark2;
    mark2.SetTextSize(0.038);
    mark2.SetTextFont(52);
    mark2.SetNDC(true);
    mark2.DrawLatex(0.23, 0.955, " Preliminary");

    mark.SetTextSize(0.03);
    // mark.DrawLatex(0.7,0.907,"19.7 fb^{-1} (8 TeV)");
    // mark.DrawLatex(0.19,0.87,"CMS Preliminary");
    if (Type == "elec" && !isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ee} < 120 GeV");
    }
    if (Type == "muon" && !isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{#mu#mu} < 120 GeV");
    }
    if (Type == "combined" || isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
    }
    TLegend* leg2;
    if (Type == "elec" && isPlot2) leg2 = new TLegend(0.50, 0.72, 0.91, 0.91);
    else leg2 = new TLegend(0.50, 0.74, 0.84, 0.88);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetLineWidth(1);
    leg2->SetNColumns(1);
    leg2->SetTextFont(42);

    if (!isPlot2) {
        leg2->AddEntry(g_ratio_phistar, "2012 data", "F");
        if (Type == "elec") {
            leg2->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ee (MadGraph)", "P");
            leg2->AddEntry(g_ph_final, "#gamma*/Z #rightarrow ee (Powheg)", "P");
            if (elec == 1) leg2->AddEntry(g_re_final, "#gamma*/Z #rightarrow ee (Resbos)", "P");
        }
        if (Type == "muon") {
            leg2->AddEntry(g_mg_final, "#gamma*/Z #rightarrow #mu#mu (MadGraph)", "P");
            leg2->AddEntry(g_ph_final, "#gamma*/Z #rightarrow #mu#mu (Powheg)", "P");
            if (elec == 1) leg2->AddEntry(g_re_final, "#gamma*/Z #rightarrow #mu#mu (Resbos)", "P");
        }
        if (Type == "combined") {
            leg2->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ll (MadGraph)", "P");
            leg2->AddEntry(g_ph_final, "#gamma*/Z #rightarrow ll (Powheg)", "P");
            if (elec == 1) leg2->AddEntry(g_re_final, "#gamma*/Z #rightarrow ll (Resbos)", "P");
        }
    } else {
        if (Type == "combined") {
            leg2->AddEntry(g_mg_final, "#gamma*/Z #rightarrow ee (2012 data)", "P");
            leg2->AddEntry(g_ph_final, "#gamma*/Z #rightarrow #mu#mu (2012 data)", "P");
            if (doMG) leg2->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (MadGraph)", "PEF");
            else leg2->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (Powheg)", "PEF");
        }
        if (Type == "elec") {
            leg2->AddEntry(g_mg_final, "data (unfolded with Powheg)", "P");
            leg2->AddEntry(g_ph_final, "data (unfolded with MadGraph)", "P");
            if (doMG) leg2->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (MadGraph)", "PEF");
            else leg2->AddEntry(g_data_final, "#gamma*/Z #rightarrow ll (Powheg)", "PEF");
        }
    }
    leg2->Draw();
    plotname = "ZShape_Ratio";
    plotname += Tag;
    if (isPlot2 && Type == "combined")plotname += "MuEl";
    if (isPlot2 && Type == "elec")plotname += "PHMG";
    else plotname += Type;
    plotname += "_";
    if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    FinalPhiRatio->SaveAs((plotname + "pdf").c_str());
    FinalPhiRatio->SaveAs((plotname + "png").c_str());
    FinalPhiRatio->SaveAs((plotname + "C").c_str());
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph) {
    vector<TGraphAsymmErrors*> g_dummyvec;
    // double x,y,errorl,errorh;
    double x, y;
    TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        // errorl=graph->GetErrorYlow(iphistar);
        // errorh=graph->GetErrorYhigh(iphistar);
        g_dummy->SetPoint(iphistar, x, y);
        g_dummy->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], 0, 0);
    }
    g_dummy->GetXaxis()->SetRangeUser(0.001, 3.2);
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

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        errorl = graph->GetErrorYlow(iphistar);
        errorh = graph->GetErrorYhigh(iphistar);
        if (!isData) {
            graphmc->GetPoint(iphistar, xmc, ymc);
            errorlmc = graphmc->GetErrorYlow(iphistar);
            errorhmc = graphmc->GetErrorYhigh(iphistar);
            g_ratio->SetPoint(iphistar, x, ymc / y);
            g_ratio->SetPointError(iphistar, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            g_ratio->SetPoint(iphistar, x, 1);
            g_ratio->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], errorl / y, errorh / y);
        }
    }
    if (isData) {
        g_ratio->SetLineWidth(2);
        g_ratio->GetXaxis()->SetRangeUser(0.001, 3.2);
        g_ratio->GetXaxis()->SetTitle("#phi*");
        g_ratio->GetXaxis()->SetTitleOffset(1.05);
        g_ratio->GetXaxis()->SetTitleSize(0.12);
        g_ratio->GetXaxis()->SetLabelSize(0.12);
        g_ratio->GetYaxis()->SetTitle("MC/Data  ");
        g_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
        g_ratio->GetYaxis()->SetTitleOffset(0.32);
        g_ratio->GetYaxis()->SetTitleSize(0.12);
        g_ratio->GetYaxis()->SetLabelSize(0.12);
        g_ratio->GetYaxis()->SetNdivisions(3, 5, 0);
        g_ratio->SetTitle(0);
    }
    return g_ratio;
}

TGraphAsymmErrors* ConvertToTGraph(TH1D* h) {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h->GetBinContent(iphistar + 1));
        g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar + 1), h->GetBinError(iphistar + 1));
    }
    return g;
}

TGraphAsymmErrors* CreateCopy(TGraphAsymmErrors* graph, double scale = 1.0) {
    double x, y, errorl, errorh;
    TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        errorl = graph->GetErrorYlow(iphistar);
        errorh = graph->GetErrorYhigh(iphistar);
        g_copy->SetPoint(iphistar, x, y * scale);
        g_copy->SetPointError(iphistar, 0, 0, errorl*scale, errorh * scale);
    }
    return g_copy;
}

void MakeFinalPlots() {

    TGraphAsymmErrors* g_data_final;
    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_ph_final;
    TGraphAsymmErrors* g_re_final = 0;
    TGraphAsymmErrors* g_ratio_re_phistar = 0;

    if (Type == "elec" && !doPHcomp) {
        std::string textn = "Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        if (doMG) textn += "MG_";
        else textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dg(textn.c_str());
        TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
        g_data_final = CreateCopy(g_data_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile mg(textn.c_str());
        TGraphAsymmErrors* g_mg_temp = (TGraphAsymmErrors*) mg.Get("Graph");
        g_mg_final = CreateCopy(g_mg_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile fph(textn.c_str());
        TGraphAsymmErrors* g_ph_temp = (TGraphAsymmErrors*) fph.Get("Graph");
        g_ph_final = CreateCopy(g_ph_temp);
    }
    if (Type == "elec" && doPHcomp) {
        std::string textn = "Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        if (doMG) textn += "MG_";
        else textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dg(textn.c_str());
        TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
        g_data_final = CreateCopy(g_data_temp);
        textn = "PhiStarAbs_NormPythia8_4C_5_3_13.root";
        cout << textn.c_str() << endl;
        TFile mg(textn.c_str());
        std::string histname = "phistar_abs";
        if (doNorm) histname = "phistar_norm";
        cout << histname.c_str() << endl;
        TH1D* h_mg_temp = (TH1D*) mg.Get(histname.c_str());
        TGraphAsymmErrors* g_mg_temp = ConvertToTGraph(h_mg_temp);
        g_mg_final = CreateCopy(g_mg_temp);
        textn = "PhiStarAbs_NormPythia6_Z2Star_5_3_13v2.root";
        cout << textn.c_str() << endl;
        TFile ph(textn.c_str());
        TH1D* h_ph_temp = (TH1D*) ph.Get(histname.c_str());
        TGraphAsymmErrors* g_ph_temp = ConvertToTGraph(h_ph_temp);
        g_ph_final = CreateCopy(g_ph_temp);
        textn = "PhiStarAbs_NormPythia8_4CHard1.6ktSoft0.35_5_3_13.root";
        cout << textn.c_str() << endl;
        TFile re(textn.c_str());
        TH1D* h_re_temp = (TH1D*) re.Get(histname.c_str());
        TGraphAsymmErrors* g_re_temp = ConvertToTGraph(h_re_temp);
        g_re_final = CreateCopy(g_re_temp);
        g_ratio_re_phistar = CreateRatio(g_data_final, g_re_final, 0);
    }
    if (Type == "combined") {
        std::string textn = "Comb_Hist_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        textn += "Born.root";
        TFile tr(textn.c_str());
        TH1D* h_data_temp = (TH1D*) tr.Get("h_Comb");
        TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_data_temp);
        g_data_final = CreateCopy(g_data_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile mg(textn.c_str());
        TGraphAsymmErrors* g_mg_temp = (TGraphAsymmErrors*) mg.Get("Graph");

        g_mg_final = CreateCopy(g_mg_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile fph(textn.c_str());
        TGraphAsymmErrors* g_ph_temp = (TGraphAsymmErrors*) fph.Get("Graph");
        g_ph_final = CreateCopy(g_ph_temp);
    }
    if (elec == 1 && !doPHcomp) {
        std::string textn = "Resbos_Born_";
        if (doNorm) textn += "Normal";
        else textn += "Absolute";
        textn += "_Central_Full_Errors.root";
        TFile tr(textn.c_str());
        TH1D* h_re_temp = (TH1D*) tr.Get("hRecoClone1");
        TGraphAsymmErrors* g_re_temp = ConvertToTGraph(h_re_temp);
        if (doNorm) g_re_final = CreateCopy(g_re_temp);
        else g_re_final = CreateCopy(g_re_temp, 0.001);
        g_ratio_re_phistar = CreateRatio(g_data_final, g_re_final, 0);
    }

    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_dummy_phistar = CreateDummy(g_data_final);
    TGraphAsymmErrors* g_ratio_phistar = CreateRatio(g_data_final, g_data_final, 1);
    TGraphAsymmErrors* g_ratio_mg_phistar = CreateRatio(g_data_final, g_mg_final, 0);
    TGraphAsymmErrors* g_ratio_ph_phistar = CreateRatio(g_data_final, g_ph_final, 0);

    if (!doPHcomp) PrintFinal(g_data_final, g_mg_final, g_ph_final, 0, g_re_final);
    PlotFinal(g_data_final, g_mg_final, g_ph_final, g_dummy_phistar, g_ratio_phistar, g_ratio_mg_phistar, g_ratio_ph_phistar, 0, g_re_final, g_ratio_re_phistar);
}

TH1D* GetHistMuon(std::string name) {
    TH1D* h_Muon = new TH1D("h_Comb", "h_Comb", nphistar, phistarBins);
    h_Muon->Sumw2();
    TFile fm(name.c_str());
    TH1F *h_temp = (TH1F*) fm.Get("hRecoClone1");
    for (uint i = 0; i < nphistar; i++) {
        h_Muon->SetBinContent(i + 1, h_temp->GetBinContent(i + 1));
        h_Muon->SetBinError(i + 1, h_temp->GetBinError(i + 1));
    }
    return h_Muon;
}

void MakeFinalPlots2() {

    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_data_muon;
    TGraphAsymmErrors* g_data_elec;

    std::string textn = "Data_Graph_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "PH_";
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    std::cout << textn.c_str() << std::endl;
    TFile dg(textn.c_str());
    TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
    g_data_elec = CreateCopy(g_data_temp);
    if (doNorm)NewError(g_data_elec, ElectNormTTBar);
    else NewError(g_data_elec, ElectAbsTTBar);
    textn = "Data_Graph_MC_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    if (doMG) textn += "MG_";
    else textn += "PH_";
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    cout << textn.c_str() << endl;
    TFile mg(textn.c_str());
    TGraphAsymmErrors* g_mg_temp = (TGraphAsymmErrors*) mg.Get("Graph");
    g_mg_final = CreateCopy(g_mg_temp);
    // textn="Data_Graph_MC_";
    // textn+=Tag;
    // if (doNorm) textn+="Norm_";
    // else        textn+="Abs_";
    // textn="Comb_Hist_MC_";
    // textn+=Tag;
    // if (doNorm) textn+="Norm_";
    // else        textn+="Abs_";
    // if (doMG)   textn+="MG_";
    // else textn+="PH_";
    // textn+="Dressed.root";
    // TFile tr2(textn.c_str());
    // TH1D* h_mg_temp= (TH1D*)tr2.Get("h_Comb");
    // TGraphAsymmErrors* g_mg_temp=ConvertToTGraph(h_mg_temp);
    // g_mg_final=CreateCopy(g_mg_temp);  

    //   textn+="MG_";
    // if (elec==0)textn+="Dressed.root";
    // if (elec==1)textn+="Born.root";
    // if (elec==2)textn+="Naked.root";
    // cout<<textn.c_str()<<endl;
    // TFile mg(textn.c_str());  
    // TGraphAsymmErrors* g_mg_temp= (TGraphAsymmErrors*)mg.Get("Graph");
    // g_mg_final=CreateCopy(g_mg_temp);
    if (Type == "combined") {
        std::string muon_name = "Muon_";
        if (doNorm) muon_name += "Normalized_";
        else muon_name += "Absolute_";
        muon_name += "Files_For_Combination_Unfold_Born_with_Powheg/Powheg_Born_";
        if (doNorm) muon_name += "Normal_";
        else muon_name += "Absolute_";
        muon_name += "Central_Full_Errors.root";
        TH1D* h_Muon_temp = GetHistMuon(muon_name);
        if (!doNorm) {
            h_Muon_temp->Scale(1. / 1000.);
        }
        g_data_muon = ConvertToTGraph(h_Muon_temp);
        if (doNorm)NewError(g_data_muon, MuonNormTTBar);
        else NewError(g_data_muon, MuonAbsTTBar);
    }
    if (Type == "elec") {
        textn = "Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dm(textn.c_str());
        TGraphAsymmErrors* g_muon_temp = (TGraphAsymmErrors*) dm.Get("Graph");
        g_data_muon = CreateCopy(g_muon_temp);
        if (doNorm)NewError(g_data_muon, ElectNormTTBar);
        else NewError(g_data_muon, ElectAbsTTBar);

    }
    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_dummy_phistar = CreateDummy(g_mg_final);
    TGraphAsymmErrors* g_ratio_elec = CreateRatio(g_mg_final, g_data_elec, 0);
    TGraphAsymmErrors* g_ratio_mg = CreateRatio(g_mg_final, g_mg_final, 1);
    TGraphAsymmErrors* g_ratio_muon = CreateRatio(g_mg_final, g_data_muon, 0);
    //
    //  PrintFinal(g_data_final,g_mg_final,g_ph_final);
    PlotFinal(g_mg_final, g_data_elec, g_data_muon, g_dummy_phistar, g_ratio_mg, g_ratio_elec, g_ratio_muon, 1);

    textn = "Comb_Hist_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "PH_";
    textn += "BornTT.root";
    TFile tr_comb(textn.c_str());
    TH1D* h_comb_temp = (TH1D*) tr_comb.Get("h_Comb");
    TGraphAsymmErrors* g_comb_temp = ConvertToTGraph(h_comb_temp);
    TGraphAsymmErrors* g_comb_final = CreateCopy(g_comb_temp);


    PrintFinal(g_comb_final, g_data_elec, g_data_muon, 1);
}

