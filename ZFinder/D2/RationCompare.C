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
#include <TStyle.h>

const bool doNorm = true;
const int elec = 1;
const int doMG = 0;
const std::string Tag = "";
const std::string Type = "elec"; //elec, muon or combined
const std::string OutType = "png";
const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 5, 10};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double yBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
size_t ny = (sizeof (yBins) / sizeof (yBins[0])) - 1;
size_t nbins = nphistar*ny;

using namespace std;

void PrintValue(ofstream &outputfile, int n, int p, double value, double error) {
    if (p > 1 || p<-1) {
        if (n + p > 1)outputfile << std::setprecision(n + p) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
        else outputfile << std::setprecision(2) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
    } else {
        if (n > 1) outputfile << std::setprecision(n) << value << "&" << error;
        else outputfile << std::setprecision(2) << value << "&" << error;
    }
}

void PrintFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, bool dataonly = 0, TGraphAsymmErrors* g_re_final = 0) {
    ofstream outputfile;
    std::string textname = "Output/Table_2D_Values_";
    textname += Tag;
    if (!dataonly) {
        textname += Type;
        textname += "_";
    }
    if (dataonly && Type == "combined")textname += "MuEl";
    if (dataonly && Type == "elec")textname += "PHMG";
    if (Type == "elec" && !doMG) textname += "PH_";
    if (Type == "elec" && doMG) textname += "MG_";
    if (doNorm) textname += "Norm_";
    else textname += "Abs_";
    if (elec == 0)textname += "Dressed.txt";
    if (elec == 1)textname += "Born.txt";
    if (elec == 2)textname += "Naked.txt";
    outputfile.open(textname.c_str());
    std::string tableheader = "";
    if (!dataonly) {
        // if (elec==1){
        //   if (!doNorm) tableheader="$\\phi^*$ range & \\multicolumn{2}{c|}{Data (pb)} & \\multicolumn{2}{c|}{MadGraph (pb)} & \\multicolumn{2}{c|}{Powheg (pb)} & \\multicolumn{2}{c|}{Resbos (pb)}\\\\ \\hline";
        //   else         tableheader="$\\phi^*$ range & \\multicolumn{2}{c|}{Data} & \\multicolumn{2}{c|}{MadGraph} & \\multicolumn{2}{c|}{Powheg} & \\multicolumn{2}{c|}{Resbos}\\\\ \\hline";
        // }
        // else{
        if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data (pb)} & \\multicolumn{2}{c|}{MadGraph (pb)} & \\multicolumn{2}{c|}{Powheg (pb)}\\\\ \\hline";
        else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c|}{Data} & \\multicolumn{2}{c|}{MadGraph} & \\multicolumn{2}{c|}{Powheg}\\\\ \\hline";
        // }     
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
    // for (size_t i=0; i<nbins;i++){
    for (size_t i = 0; i < (nbins); i++) {
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
        cout << temp_d << " " << temp_m << " " << temp_p << " " << temp2_d << " " << temp2_m << " " << temp2_p << endl;
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
        if (!dataonly && elec == 1 && g_re_final) {
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
        if (!dataonly && elec == 1 && g_re_final) {
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
        cout << n_d << " " << n_m << " " << n_p << " " << n_d << " " << n_m << " " << n_p << endl;

        if (i % nphistar == 0) {
            outputfile << std::fixed;
            outputfile << std::setprecision(1);
            outputfile << "\n";
            outputfile << yBins[i / nphistar] << "|Z(y)|" << yBins[(i / nphistar) + 1] << "\n";
            outputfile << "\n";
            outputfile << tableheader << "\n";
            std::cout << tableheader << endl;
        }
        outputfile << std::fixed << std::setprecision(3) << phistarBins[i % nphistar] << "-" << phistarBins[i % nphistar + 1] << " & ";
        if (!(dataonly && Type == "combined")) {
            if (!(Type == "elec") || !(dataonly)) {
                cout << i << " " << y << " " << g_data_final->GetErrorYhigh(i) << endl;
                PrintValue(outputfile, n_d + 1, p_d, y, g_data_final->GetErrorYhigh(i));
                cout << i << " " << y << " " << g_data_final->GetErrorYhigh(i) << endl;
                outputfile << " & ";
            }
            PrintValue(outputfile, n_m + 1, p_m, ymg, g_mg_final->GetErrorYhigh(i));
            outputfile << " & ";
            PrintValue(outputfile, n_p + 1, p_p, yph, g_ph_final->GetErrorYhigh(i));
            if (!dataonly && elec == 1 && g_re_final) {
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

vector<TGraphAsymmErrors*> SplitGraph(TGraphAsymmErrors* graph, bool doXerrors = 0) {
    vector<TGraphAsymmErrors*> v;
    for (uint i = 0; i < ny; i++) {
        TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
        for (uint j = 0; j < nphistar; j++) {
            int bin = i * nphistar + j;
            double x, y;
            graph->GetPoint(bin, x, y);
            g->SetPoint(j, (phistarBins[j] + phistarBins[j + 1]) / 2., y);
            g->SetPointError(j, 0, 0, graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            if (doXerrors) g->SetPointError(j, (phistarBins[j + 1] - phistarBins[j]) / 2., (phistarBins[j + 1] - phistarBins[j]) / 2., graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            cout << i << " " << j << " " << bin << " " << (phistarBins[j] + phistarBins[j + 1]) / 2. << " " << (phistarBins[j + 1] - phistarBins[j]) / 2. << " " << y << " " << graph->GetErrorYlow(bin) << " " << graph->GetErrorYhigh(bin) << endl;
        }
        v.push_back(g);
    }
    return v;
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph) {
    vector<TGraphAsymmErrors*> g_dummyvec;
    // double x,y,errorl,errorh;
    double x, y;
    TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        g_dummy->SetPoint(iphistar, x, y);
        g_dummy->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], 0, 0);
    }
    g_dummy->GetXaxis()->SetRangeUser(0.001, 10);
    g_dummy->GetXaxis()->SetTitleOffset(1.05);
    g_dummy->GetXaxis()->SetTitle(0);
    g_dummy->GetXaxis()->SetTitleSize(0.05);
    g_dummy->GetXaxis()->SetLabelSize(0.05);
    g_dummy->GetYaxis()->SetTitleOffset(1.05);
    //  if (!doNorm)g_dummy->GetYaxis()->SetTitleOffset(1.05); 
    g_dummy->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi*dy ");
    if (!doNorm) g_dummy->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi*dy (pb)");
    g_dummy->GetYaxis()->SetTitleSize(0.05);
    g_dummy->GetYaxis()->SetLabelSize(0.05);
    g_dummy->SetFillColor(kWhite);
    g_dummy->SetTitle(0);
    return g_dummy;
}

vector<TGraphAsymmErrors*> CreateDummy(vector<TGraphAsymmErrors*> graphs) {
    vector<TGraphAsymmErrors*> v;
    for (uint i = 0; i < graphs.size(); i++) {
        TGraphAsymmErrors* d = CreateDummy(graphs[i]);
        v.push_back(d);
    }
    return v;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbins);
    for (size_t ibin = 0; ibin < nbins; ibin++) {
        graph->GetPoint(ibin, x, y);
        if (!isData) {
            graphmc->GetPoint(ibin, xmc, ymc);
            errorlmc = graph->GetErrorYlow(ibin);
            errorhmc = graph->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, ymc / y);
            g_ratio->SetPointError(ibin, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            errorl = graph->GetErrorYlow(ibin);
            errorh = graph->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, 1);
            g_ratio->SetPointError(ibin, 0, 0, errorl / y, errorh / y);
            cout << errorl / y << " " << errorh / y << endl;
        }
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
    TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nbins);
    for (size_t ibin = 0; ibin < nbins; ibin++) {
        graph->GetPoint(ibin, x, y);
        errorl = graph->GetErrorYlow(ibin);
        errorh = graph->GetErrorYhigh(ibin);
        g_copy->SetPoint(ibin, x, y * scale);
        g_copy->SetPointError(ibin, 0, 0, errorl*scale, errorh * scale);
    }
    return g_copy;
}

TGraphAsymmErrors* GetIntValues(vector<TGraphAsymmErrors*> All) {//get the values intigrated 
    //Unfolded Data  
    TGraphAsymmErrors* Total = new TGraphAsymmErrors(nphistar);

    for (uint j = 0; j < nphistar; j++) {

        double X = 0;
        double YTotal = 0;
        double NewError = 0;
        //double MCX = 0;
        for (uint i = 0; i < ny; i++) {
            double testx, testy;
            double HoldYError = 0;
            All[i]->GetPoint(j, testx, testy);
            double holdy = testy;
            X = testx;
            YTotal += holdy;
            HoldYError = All[i]->GetErrorY(j);
            NewError = sqrt(NewError * NewError + HoldYError * HoldYError);
        }
        Total->SetPoint(j, X, YTotal);
        Total->SetPointEYhigh(j, NewError);
        Total->SetPointEYlow(j, NewError);
    }


    return Total;
}

void PlotFinalZach(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2 = 0, TGraphAsymmErrors* g_re_final = 0, TGraphAsymmErrors* g_ratio_re_phistar = 0) {

    vector<TGraphAsymmErrors*> g_data = SplitGraph(g_data_final);
    vector<TGraphAsymmErrors*> g_mg = SplitGraph(g_mg_final);
    vector<TGraphAsymmErrors*> g_ph = SplitGraph(g_ph_final);
    vector<TGraphAsymmErrors*> r_data = SplitGraph(g_ratio_phistar, 1);
    vector<TGraphAsymmErrors*> r_mg = SplitGraph(g_ratio_mg_phistar);
    vector<TGraphAsymmErrors*> r_ph = SplitGraph(g_ratio_ph_phistar);
    vector<TGraphAsymmErrors*> g_re;
    vector<TGraphAsymmErrors*> r_re;
    if (g_re_final) {
        g_re = SplitGraph(g_re_final);
        r_re = SplitGraph(g_ratio_re_phistar);
    }



    for (uint i = 0; i < ny; i++) {
        //  for (uint i=0; i<5; i++){ 
        //Okay need to find new errors
        for (uint j = 0; j < nphistar; j++) {

            double MGSampleX, MGSampleY, PHSampleX, PHSampleY, DataX, DataY;
            double MGError, PHError, DataError;
            g_data[i]->GetPoint(j, DataX, DataY);
            g_mg[i]->GetPoint(j, MGSampleX, MGSampleY);
            g_ph[i]->GetPoint(j, PHSampleX, PHSampleY);
            DataError = g_data[i]->GetErrorY(j);
            PHError = g_ph[i]->GetErrorY(j);
            double NewErrorMGRatio = (sqrt(MGError * MGError / (DataY * DataY) + DataError * DataError * MGSampleY * MGSampleY / (DataY * DataY * DataY * DataY)));
            double NewErrorPHRatio = (sqrt(PHError * PHError / (DataY * DataY) + DataError * DataError * PHSampleY * PHSampleY / (DataY * DataY * DataY * DataY)));
            r_ph[i]->SetPointEYlow(j, NewErrorPHRatio);
            r_ph[i]->SetPointEYhigh(j, NewErrorPHRatio);
            r_mg[i]->SetPointEYlow(j, NewErrorMGRatio);
            r_mg[i]->SetPointEYhigh(j, NewErrorMGRatio);

            //OffSetter
            double x, y;
            r_mg[i]->GetPoint(j, x, y);
            r_mg[i]->SetPoint(j, x, y + (i + 1) * .1);
            r_ph[i]->GetPoint(j, x, y);
            r_ph[i]->SetPoint(j, x, y + (i + 1)* .1);


        }
    }

    TCanvas* C1 = new TCanvas("c1", "c1", 800, 900);
    C1->SetLogx();
    C1->SetRightMargin(.02);
    gStyle->SetOptStat(0);
    vector<TGraphAsymmErrors*> g_dummy = CreateDummy(g_data);

    TLine* line0 = new TLine(0, 1, 10, 1);
    TLine* line1 = new TLine(0, 1.1, 10, 1.1);
    TLine* line2 = new TLine(0, 1.2, 10, 1.2);
    TLine* line3 = new TLine(0, 1.3, 10, 1.3);
    TLine* line4 = new TLine(0, 1.4, 10, 1.4);
    TLine* line5 = new TLine(0, 1.5, 10, 1.5);
    TLine* line6 = new TLine(0, 1.6, 10, 1.6);
    //
    TGraphAsymmErrors* OrginPH = GetIntValues(g_ph); //contains the intigrated ratios
    TGraphAsymmErrors* OrginMG = GetIntValues(g_mg);
    ;
    TGraphAsymmErrors* OrginData = GetIntValues(g_data);

    TGraphAsymmErrors* OrginRatioPH = CreateRatio(OrginData, OrginPH, 0);
    TGraphAsymmErrors* OrginRatioMG = CreateRatio(OrginData, OrginMG, 0);
    for (uint j = 0; j < nphistar; j++) {

        double MGSampleX, MGSampleY, PHSampleX, PHSampleY, DataX, DataY;
        double MGError, PHError, DataError;
        OrginData->GetPoint(j, DataX, DataY);
        OrginMG->GetPoint(j, MGSampleX, MGSampleY);
        OrginPH->GetPoint(j, PHSampleX, PHSampleY);
        DataError = OrginData->GetErrorY(j);
        PHError = OrginPH->GetErrorY(j);
        double NewErrorMGRatio = (sqrt(MGError * MGError / (DataY * DataY) + DataError * DataError * MGSampleY * MGSampleY / (DataY * DataY * DataY * DataY)));
        double NewErrorPHRatio = (sqrt(PHError * PHError / (DataY * DataY) + DataError * DataError * PHSampleY * PHSampleY / (DataY * DataY * DataY * DataY)));
        OrginRatioPH->SetPointEYlow(j, NewErrorPHRatio);
        OrginRatioPH->SetPointEYhigh(j, NewErrorPHRatio);
        OrginMG->SetPointEYlow(j, NewErrorMGRatio);
        OrginMG->SetPointEYhigh(j, NewErrorMGRatio);

       
    }



    OrginRatioPH->SetMarkerStyle(kStar);
    OrginRatioMG->SetMarkerStyle(kStar);
    //
    r_ph[0]->SetFillStyle(0);
    r_ph[1]->SetFillStyle(0);
    r_ph[2]->SetFillStyle(0);
    r_ph[3]->SetFillStyle(0);
    r_ph[4]->SetFillStyle(0);
    r_ph[5]->SetFillStyle(0);




    r_ph[0]->SetMarkerStyle(kFullTriangleUp);
    //r_ph[0]->SetMarkerColor(kBlue);
    //r_ph[0]->SetLineColor(kBlue);
    r_ph[1]->SetMarkerStyle(kOpenTriangleUp);
    //r_ph[1]->SetMarkerColor(kGreen + 2);
    //r_ph[1]->SetLineColor(kGreen + 2);
    r_ph[2]->SetMarkerStyle(kFullTriangleDown);
    //r_ph[2]->SetMarkerColor(kMagenta + 2);
    //r_ph[2]->SetLineColor(kMagenta + 2);
    r_ph[3]->SetMarkerStyle(kOpenTriangleDown);
    r_ph[3]->SetMarkerColor(kBlack);
    r_ph[3]->SetLineColor(kBlack);
    r_ph[4]->SetMarkerStyle(kFullCircle);
    //r_ph[4]->SetMarkerColor(kRed);
    //r_ph[4]->SetLineColor(kRed);
    r_ph[5]->SetMarkerStyle(kOpenCircle);
    //r_ph[5]->SetMarkerSize(2);
    //r_ph[5]->SetMarkerColor(kCyan + 1);
    //r_ph[5]->SetLineColor(kCyan + 1);

    //r_ph[0]->GetYaxis()->SetRangeUser(0.5, 2);
    //r_ph[0]->GetYaxis()->SetTitleOffset(1.15);
    ////r_ph[0]->GetYaxis()->SetLabelColor(kWhite);
    //r_ph[0]->SetTitle("");
    //r_ph[0]->GetXaxis()->SetTitle("#phi*");
    //r_ph[0]->GetXaxis()->SetTitleOffset(.7);
    //r_ph[0]->GetXaxis()->SetTitleSize(0.05);
    ////r_ph[0]->GetXaxis()->SetLabelSize(0.12);
    //r_ph[0]->GetYaxis()->SetTitle("(POWHEG+Pythia6 Z2*)/Data  ");

    TH1D* Dummy = new (TH1D) ("test", "", 1, .001, 10);
    Dummy->GetYaxis()->SetRangeUser(.8, 2.4);
    Dummy->SetTitle("");
    Dummy->GetXaxis()->SetTitle("#phi*   ");
    Dummy->GetYaxis()->SetTitle("(POWHEG+Pythia6 Z2*)/Data  ");
    Dummy->GetYaxis()->SetTitleOffset(1.15);
    Dummy->Draw();

    r_ph[0]->Draw("PEsame");
    r_ph[1]->Draw("PEsame");
    r_ph[2]->Draw("PEsame");
    r_ph[3]->Draw("PEsame");
    r_ph[4]->Draw("PEsame");
    r_ph[5]->Draw("PEsame");
    OrginRatioPH->Draw("PEsame");
    TLegend* legPow1;
    TLegend* legPow2;
    TLegend* legPow3;
    //if (isPlot2) leg = new TLegend(0.15, 0.06, 0.80, 0.27);
    legPow1 = new TLegend(0.1, 0.73, 0.4, 0.92); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legPow1->SetLineStyle(3);
    legPow1->SetFillStyle(4000);
    legPow1->SetFillColor(0);
    legPow1->SetLineColor(0);
    legPow1->AddEntry(r_ph[0], "0.0 < |y_{ee}| < 0.4+0.1", "P");
    legPow1->AddEntry(r_ph[1], "(0.4 < |y_{ee}| < 0.8) +0.2", "P");
    legPow1->AddEntry(r_ph[2], "(0.8 < |y_{ee}| < 1.2) +0.3", "P");

    legPow2 = new TLegend(0.4, 0.73, 0.7, 0.92); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legPow2->SetLineStyle(3);
    legPow2->SetFillStyle(4000);
    legPow2->SetFillColor(0);
    legPow2->SetLineColor(0);
    legPow2->AddEntry(r_ph[3], "(1.2 < |y_{ee}| < 1.6) +0.4", "P");
    legPow2->AddEntry(r_ph[4], "(1.6 < |y_{ee}| < 2.0) +0.5", "P");
    legPow2->AddEntry(r_ph[5], "(2.0 < |y_{ee}| < 2.4) +0.6", "P");
    legPow3 = new TLegend(0.25, 0.7, 0.55, 0.73); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legPow3->SetLineStyle(3);
    legPow3->SetFillStyle(4000);
    legPow3->SetFillColor(0);
    legPow3->SetLineColor(0);
    //legPow3->SetTextSize(0.043);
    legPow3->AddEntry(OrginRatioPH, "(0.0<|y_{ee}|<2.4 )", "P");
    TLatex mark;
    mark.SetTextSize(0.043);
    mark.SetTextFont(42);
    mark.SetTextSize(0.035);
    mark.SetNDC(true);
    mark.DrawLatex(0.65, 0.905, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.905, "CMS Preliminary");

    line0->Draw("same");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");
    line5->Draw("same");
    line6->Draw("same");

    legPow1->Draw();
    legPow2->Draw();
    legPow3->Draw();
    string RatioName = "CompPH";
    if (doNorm)RatioName += "Norm.png";
    else RatioName += "Abs.png";
    C1->Print(RatioName.c_str());
    delete C1;

    TCanvas* C2 = new TCanvas("c2", "c1", 800, 900);
    C2->SetLogx();



    //r_mg[0]->GetXaxis()->SetRangeUser(.001,10);

    gStyle->SetOptStat(0);
    r_mg[0]->SetFillStyle(0);
    r_mg[1]->SetFillStyle(0);
    r_mg[2]->SetFillStyle(0);
    r_mg[3]->SetFillStyle(0);
    r_mg[4]->SetFillStyle(0);
    r_mg[5]->SetFillStyle(0);

    r_mg[0]->SetMarkerStyle(kFullTriangleUp);
    //r_mg[0]->SetMarkerColor(kBlue);
    //r_mg[0]->SetLineColor(kBlue);
    r_mg[1]->SetMarkerStyle(kOpenTriangleUp);
    //r_mg[1]->SetMarkerColor(kGreen + 2);
    //r_mg[1]->SetLineColor(kGreen + 2);
    r_mg[2]->SetMarkerStyle(kFullTriangleDown);
    //r_mg[2]->SetMarkerColor(kMagenta + 2);
    //r_mg[2]->SetLineColor(kMagenta + 2);
    r_mg[3]->SetMarkerStyle(kOpenTriangleDown);
    r_mg[3]->SetMarkerColor(kBlack);
    r_mg[3]->SetLineColor(kBlack);
    r_mg[4]->SetMarkerStyle(kFullCross);
    //r_mg[4]->SetMarkerColor(kRed);
    //r_mg[4]->SetLineColor(kRed);
    r_mg[5]->SetMarkerStyle(kOpenCircle);
    //r_mg[5]->SetMarkerSize(2);
    //r_mg[5]->SetMarkerColor(kCyan + 1);
    //r_mg[5]->SetLineColor(kCyan + 1);

    //r_mg[0]->GetYaxis()->SetRangeUser(0.5, 2);
    //r_mg[0]->GetYaxis()->SetLabelColor(kWhite);
    Dummy->GetYaxis()->SetTitle("(MADGRAPH+Pythia6 Z2*)/Data  ");
    Dummy->Draw();
    r_mg[0]->Draw("PEsame");
    r_mg[1]->Draw("PEsame");
    r_mg[2]->Draw("PEsame");
    r_mg[3]->Draw("PEsame");
    r_mg[4]->Draw("PEsame");
    r_mg[5]->Draw("PEsame");
    OrginRatioMG->Draw("PEsame");
    line0->Draw("same");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");
    line5->Draw("same");
    line6->Draw("same");
    TLegend* legMad1;
    TLegend* legMad2;
    TLegend* legMad3;
    //if (isPlot2) leg2 = new TLegend(0.15, 0.06, 0.80, 0.27);
    legMad1 = new TLegend(0.1, 0.73, 0.4, 0.92); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legMad1->SetLineStyle(3);
    legMad1->SetFillStyle(4000);
    legMad1->SetFillColor(0);
    legMad1->SetLineColor(0);
    legMad1->AddEntry(r_mg[0], "(0.0 < |y_{ee}| < 0.4", "P");
    legMad1->AddEntry(r_mg[1], "(0.4 < |y_{ee}| < 0.8) +0.1", "P");
    legMad1->AddEntry(r_mg[2], "(0.8 < |y_{ee}| < 1.2) +.2", "P");

    legMad2 = new TLegend(0.4, 0.73, 0.7, 0.92); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legMad2->SetLineStyle(3);
    legMad2->SetFillStyle(4000);
    legMad2->SetFillColor(0);
    legMad2->SetLineColor(0);
    legMad2->AddEntry(r_mg[3], "(1.2 < |y_{ee}| < 1.6) +0.3", "P");
    legMad2->AddEntry(r_mg[4], "(1.6 < |y_{ee}| < 2.0) +0.4", "P");
    legMad2->AddEntry(r_mg[5], "(2.0 < |y_{ee}| < 2.4) +0.5", "P");
    
    legMad3 = new TLegend(0.25, 0.67, 0.55, 0.73); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
    legMad3->SetLineStyle(3);
    legMad3->SetFillStyle(4000);
    legMad3->SetFillColor(0);
    legMad3->SetLineColor(0);
    legMad3->AddEntry(OrginRatioMG, "(0.0< |y_{ee}|<2.4 )", "P");

    mark.DrawLatex(0.65, 0.905, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.905, "CMS Preliminary");

    legMad1->Draw();
    legMad2->Draw();
    legMad3->Draw();

    RatioName = "CompMG";
    if (doNorm)RatioName += "Norm.png";
    else RatioName += "Abs.png";
    C2->Print(RatioName.c_str());
    delete C2;
    delete Dummy;
    //Okay this should contain everything I hope 

}

void CompareAllPlots() {

    TGraphAsymmErrors* g_data_final;
    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_ph_final;
    TGraphAsymmErrors* g_re_final = 0;
    TGraphAsymmErrors* g_ratio_re_phistar = 0;

    if (Type == "elec") {
        std::string textn = "Output/Data_Graph_";
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

        textn = "Output/MC_Graph_";
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

        textn = "Output/MC_Graph_";
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

    // if (Type=="combined"){
    // }
    // if (elec==1){
    // }

    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_ratio_phistar = CreateRatio(g_data_final, g_data_final, 1);
    TGraphAsymmErrors* g_ratio_mg_phistar = CreateRatio(g_data_final, g_mg_final, 0);
    TGraphAsymmErrors* g_ratio_ph_phistar = CreateRatio(g_data_final, g_ph_final, 0);

    //PrintFinal(g_data_final,g_mg_final,g_ph_final,0,g_re_final);
    PlotFinalZach(g_data_final, g_mg_final, g_ph_final, g_ratio_phistar, g_ratio_mg_phistar, g_ratio_ph_phistar, 0, g_re_final, g_ratio_re_phistar);
}


