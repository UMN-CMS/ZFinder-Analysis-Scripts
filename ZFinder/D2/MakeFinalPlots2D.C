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

void PrintValue(ofstream &outputfile, int n, int p, double value, double error) {
    if (p > 1 || p<-1) {
        if (n + p > 1)outputfile << std::setprecision(n + p) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
        else outputfile << std::setprecision(2) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
    } else {
        if (n > 1) outputfile << std::setprecision(n) << value << "&" << error;
        else outputfile << std::setprecision(2) << value << "&" << error;
    }
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
            //cout << i << " " << j << " " << bin << " " << (phistarBins[j] + phistarBins[j + 1]) / 2. << " " << (phistarBins[j + 1] - phistarBins[j]) / 2. << " " << y << " " << graph->GetErrorYlow(bin) << " " << graph->GetErrorYhigh(bin) << endl;
        }
        v.push_back(g);
    }
    return v;
}



TGraphAsymmErrors* ResbosFromRaj(int FType = 0) {
    TGraphAsymmErrors* ResHolder = new TGraphAsymmErrors(nphistar);
    string FName = "";
    if (doNorm) {
        if (FType == 0)FName = "Resbos_0_1000_8TeV_Normalized.root";
        else if (FType == 1)FName = "PowhegPythia8_NLO_AODSIM_8TeV_Normalized.root";
        else if (FType == 2)FName = "AMCAT_NLO_AODSIM_8TeV_Normalized.root";
    } else {
        if (FType == 0)FName = "Resbos_0_1000_8TeV_Absolute.root";
        else if (FType == 1)FName = "PhiStar_PowhegPythia8_phistar_Status3_Absolute.root";
        else if (FType == 2)
            FName = "PhiStar_PowhegPythia8_phistar_Status3_Absolute.root";
    }
    TFile ResFile(FName.c_str());
    TH1D* Bin0 = (TH1D*) ResFile.Get("PhiStar_YBin_0");
    TH1D* Bin1 = (TH1D*) ResFile.Get("PhiStar_YBin_1");
    TH1D* Bin2 = (TH1D*) ResFile.Get("PhiStar_YBin_2");
    TH1D* Bin3 = (TH1D*) ResFile.Get("PhiStar_YBin_3");
    TH1D* Bin4 = (TH1D*) ResFile.Get("PhiStar_YBin_4");
    TH1D* Bin5 = (TH1D*) ResFile.Get("PhiStar_YBin_5");
    //TH1D* Allstuff = (TH1D*) ResFile.Get("PhiStar")
    if (Bin0 == 0)cout << "missing Bin0" << endl;
    if (Bin1 == 0)cout << "missing Bin1" << endl;
    if (Bin2 == 0)cout << "missing Bin2" << endl;
    if (Bin3 == 0)cout << "missing Bin3" << endl;
    if (Bin4 == 0)cout << "missing Bin4" << endl;
    if (Bin5 == 0)cout << "missing Bin5" << endl;





    for (uint i = 1; i <= nphistar; i++) {
        uint iphistar = i - 1;
        ResHolder->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2, Bin0->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 1 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 1 * phistarBins[nphistar], Bin1->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 2 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 2 * phistarBins[nphistar], Bin2->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 3 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 3 * phistarBins[nphistar], Bin3->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 4 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 4 * phistarBins[nphistar], Bin4->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 5 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 5 * phistarBins[nphistar], Bin5->GetBinContent(i));

        //cout<<"our bin content is "<<Bin1->GetBinContent(i)<<endl;

        ResHolder->SetPointEYhigh(iphistar, Bin0->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 1 * nphistar, Bin1->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 2 * nphistar, Bin2->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 3 * nphistar, Bin3->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 4 * nphistar, Bin4->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 5 * nphistar, Bin5->GetBinError(i));

        ResHolder->SetPointEYlow(iphistar, Bin0->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 1 * nphistar, Bin1->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 2 * nphistar, Bin2->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 3 * nphistar, Bin3->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 4 * nphistar, Bin4->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 5 * nphistar, Bin5->GetBinError(i));

        ResHolder->SetPointEYhigh(iphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 1 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 2 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 3 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 4 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 5 * nphistar, .000001);

        ResHolder->SetPointEYlow(iphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 1 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 2 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 3 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 4 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 5 * nphistar, .000001);

    }
    return ResHolder;
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
    cout << "test 2" << endl;
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
    cout << "test 3" << endl;
    for (size_t i = 0; i < (nbins); i++) {
        if (i % 10 == 0)cout << "still going at " << i / nbins << "%" << endl;
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
        //cout << temp_d << " " << temp_m << " " << temp_p << " " << temp2_d << " " << temp2_m << " " << temp2_p << endl;
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
        cout << "test 4" << endl;
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
            cout << "test 4.1" << endl;
            cout << " and our point is " << i << endl;
            while (temp_r < 1) {
                temp_r = temp_r * 10.;
                n_r++;
            }
            cout << "test 4.2" << endl;
            while (temp_r > 10) {
                temp_r = temp_r / 10.;
                n_r = n_r - 1;
            }
        }
        cout << "test 5" << endl;
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
        cout << "test 6" << endl;
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
        cout << "test 7" << endl;
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
                //cout << i << " " << y << " " << g_data_final->GetErrorYhigh(i) << endl;
                PrintValue(outputfile, n_d + 1, p_d, y, g_data_final->GetErrorYhigh(i));
                //cout << i << " " << y << " " << g_data_final->GetErrorYhigh(i) << endl;
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

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2 = 0, TGraphAsymmErrors* g_re_final = 0, TGraphAsymmErrors* g_ratio_re_phistar = 0) {
    cout << "test 2.5" << endl;
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
    cout << "test 3" << endl;
    vector<TGraphAsymmErrors*> g_dummy = CreateDummy(g_data);
    for (uint i = 0; i < ny; i++) {
        //  for (uint i=0; i<5; i++){ 
        std::ostringstream strs;
        strs << i;
        std::string Canvasname = "EventEff_Bin" + strs.str();
        TCanvas* FinalPhiTot = new TCanvas(Canvasname.c_str(), Canvasname.c_str(), 800, 900);
        FinalPhiTot->Divide(1, 2);
        FinalPhiTot->cd(1);
        gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
        gPad->SetBottomMargin(0.01);
        gPad->SetTopMargin(0.06);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        g_dummy[i]->GetXaxis()->SetRangeUser(0.001, 10.0);
        if (doNorm) g_dummy[i]->GetYaxis()->SetRangeUser(0.00000005, 1000.0);
        else g_dummy[i]->GetYaxis()->SetRangeUser(0.00005, 1000000.0);
        g_dummy[i]->GetXaxis()->CenterTitle();
        g_dummy[i]->GetYaxis()->CenterTitle();
        g_dummy[i]->Draw("A2");
        g_mg[i]->SetMarkerColor(kBlue - 7);
        g_mg[i]->SetLineColor(kBlue - 7);
        g_mg[i]->SetMarkerSize(1);
        g_mg[i]->SetLineWidth(2);
        g_mg[i]->SetMarkerStyle(21);
        g_mg[i]->Draw("PEsame");
        g_ph[i]->SetMarkerColor(kRed);
        g_ph[i]->SetLineColor(kRed);
        g_ph[i]->SetMarkerSize(1);
        g_ph[i]->SetLineWidth(2);
        g_ph[i]->SetMarkerStyle(22);
        g_ph[i]->Draw("PEsame");
        if (!isPlot2 && elec == 1 && g_re_final) {
            g_re[i]->SetMarkerColor(kGreen + 1);
            g_re[i]->SetLineColor(kGreen + 1);
            g_re[i]->SetMarkerSize(1);
            g_re[i]->SetLineWidth(2);
            g_re[i]->SetMarkerStyle(23);
            g_re[i]->Draw("PEsame");
        }
        g_data[i]->SetFillColor(kYellow);
        g_data[i]->SetMarkerSize(1);
        g_data[i]->SetLineWidth(2);
        g_data[i]->SetMarkerStyle(20);
        g_data[i]->Draw("PEsame");
        g_data[i]->SetFillColor(kYellow);

        TLegend* leg;
        if (isPlot2) leg = new TLegend(0.15, 0.06, 0.80, 0.27);
        else leg = new TLegend(0.15, 0.06, 0.80, 0.31); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        if (!isPlot2) {
            leg->AddEntry(g_data[i], "2012 data", "PEF");
            if (Type == "elec") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ee (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow ee (Powheg+Pythia6 Z2*)", "P");
                //ToDo AMCAt decisions
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ee (Resbos)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ee (POWHEG+Pythia8)", "P");
            }
            if (Type == "muon") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow #mu#mu (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow #mu#mu (Powheg+Pythia6 Z2*)", "P");
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow #mu#mu (Resbos)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow #mu#mu (POWHEG+Pythia8 CUETP8M1)", "P");
            }
            if (Type == "combined") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "P");
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ll (Resbos)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ll (POWHEG+Pythia8 CUETP8M1)", "P");
            }
        } else {
            if (Type == "combined") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ee (2012 data)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow #mu#mu (2012 data)", "P");
                if (doMG)leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
                else leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
            }
            if (Type == "elec") {
                leg->AddEntry(g_mg[i], "data (unfolded with Powheg)", "P");
                leg->AddEntry(g_ph[i], "data (unfolded with MadGraph)", "P");
                if (doMG)leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
                else leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
            }
        }
        leg->Draw();

        TLatex mark3;
        mark3.SetTextSize(0.043);
        mark3.SetTextFont(42);
        mark3.SetNDC(true);
        mark3.DrawLatex(0.71, 0.955, "19.7 fb^{-1} (8 TeV)");
        // TLatex mark1;
        // mark1.SetTextSize(0.05);
        // mark1.SetTextFont(61);
        // mark1.SetNDC(true);
        // mark1.DrawLatex(0.15,0.955,"CMS");
        // TLatex mark2;
        // mark2.SetTextSize(0.038);
        // mark2.SetTextFont(52);
        // mark2.SetNDC(true);
        // mark2.DrawLatex(0.23,0.955," Preliminary");
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
        if (Type == "combined" || isPlot2) {
            mark.DrawLatex(0.53, 0.88, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
            mark.DrawLatex(0.53, 0.81, "p_{T}^{l_{0}} > 30 GeV,  p_{T}^{l_{1}} > 20 GeV");
            mark.DrawLatex(0.53, 0.74, "60 GeV < M_{ll} < 120 GeV");
        }
        if (i == 0) mark.DrawLatex(0.2, 0.88, "0.0 < |y_{ee}| < 0.4");
        if (i == 1) mark.DrawLatex(0.2, 0.88, "0.4 < |y_{ee}| < 0.8");
        if (i == 2) mark.DrawLatex(0.2, 0.88, "0.8 < |y_{ee}| < 1.2");
        if (i == 3) mark.DrawLatex(0.2, 0.88, "1.2 < |y_{ee}| < 1.6");
        if (i == 4) mark.DrawLatex(0.2, 0.88, "1.6 < |y_{ee}| < 2.0");
        if (i == 5) mark.DrawLatex(0.2, 0.88, "2.0 < |y_{ee}| < 2.4");
        FinalPhiTot->cd(2);
        gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
        gPad->SetBottomMargin(0.37);
        gPad->SetTopMargin(0.01);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetLogx(1);

        r_data[i]->SetLineWidth(2);
        r_data[i]->GetXaxis()->SetRangeUser(0.001, 10);
        r_data[i]->GetXaxis()->SetTitle("#phi*");
        r_data[i]->GetXaxis()->SetTitleOffset(1.05);
        r_data[i]->GetXaxis()->SetTitleSize(0.12);
        r_data[i]->GetXaxis()->SetLabelSize(0.12);
        r_data[i]->GetYaxis()->SetTitle("MC/Data  ");
        //   r_data[i]->GetYaxis()->SetRangeUser(0.8,1.2);
        r_data[i]->GetYaxis()->SetTitleOffset(0.32);
        r_data[i]->GetYaxis()->SetTitleSize(0.12);
        r_data[i]->GetYaxis()->SetLabelSize(0.12);
        r_data[i]->SetTitle(0);
        if (isPlot2) r_data[i]->GetYaxis()->SetTitle("Data/MC  ");
        else r_data[i]->GetYaxis()->SetTitle("MC/Data   ");
        r_data[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
        r_data[i]->GetYaxis()->SetNdivisions(2, 5, 0);
        if (isPlot2) r_data[i]->GetYaxis()->SetRangeUser(0.3, 1.35);
        if (isPlot2 && (!doNorm)) r_data[i]->GetYaxis()->SetRangeUser(0.7, 1.3);
        r_data[i]->GetYaxis()->SetNdivisions(2, 5, 0);
        r_data[i]->GetYaxis()->SetTitleOffset(0.45);
        r_data[i]->SetFillColor(kYellow);
        r_data[i]->GetXaxis()->SetTitleSize(0.15);
        r_data[i]->GetXaxis()->CenterTitle();
        r_data[i]->GetYaxis()->CenterTitle();
        r_data[i]->Draw("APE2");
        r_mg[i]->SetMarkerSize(1);
        r_mg[i]->SetLineWidth(2);
        r_mg[i]->SetMarkerStyle(21);
        r_mg[i]->SetMarkerColor(kBlue - 7);
        r_mg[i]->SetLineColor(kBlue - 7);
        r_mg[i]->Draw("PEsame");
        r_ph[i]->SetMarkerSize(1);
        r_ph[i]->SetLineWidth(2);
        r_ph[i]->SetMarkerStyle(22);
        r_ph[i]->SetMarkerColor(kRed);
        r_ph[i]->SetLineColor(kRed);
        r_ph[i]->Draw("PEsame");
        if (!isPlot2 && elec == 1 && g_re_final) {
            r_re[i]->SetMarkerSize(1);
            r_re[i]->SetLineWidth(2);
            r_re[i]->SetMarkerStyle(23);
            r_re[i]->SetMarkerColor(kGreen + 1);
            r_re[i]->SetLineColor(kGreen + 1);
            r_re[i]->Draw("PEsame");
        }

        std::string plotname = "Plots/ZShape_2D_";
        plotname += "Bin" + strs.str() + "_";
        plotname += Tag;
        if (isPlot2 && Type == "combined")plotname += "MuEl";
        else if (isPlot2 && Type == "elec")plotname += "PHMG";
        else plotname += Type;
        plotname += "_";
        if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
        if ((Type == "elec" || isPlot2) && doMG) plotname += "MG_";
        if (doNorm) plotname += "Norm_";
        else plotname += "Abs_";
        if (elec == 0)plotname += "Dressed.";
        if (elec == 1)plotname += "Born.";
        if (elec == 2)plotname += "Naked.";
        FinalPhiTot->SaveAs((plotname + OutType).c_str());
    }
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbins);
    for (size_t ibin = 0; ibin < nbins; ibin++) {
        graph->GetPoint(ibin, x, y);
        if (!isData) {
            graphmc->GetPoint(ibin, xmc, ymc);
            errorlmc = graphmc->GetErrorYlow(ibin);
            errorhmc = graphmc->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, ymc / y);
            g_ratio->SetPointError(ibin, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            errorl = graph->GetErrorYlow(ibin);
            errorh = graph->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, 1);
            g_ratio->SetPointError(ibin, 0, 0, errorl / y, errorh / y);
            //cout << errorl / y << " " << errorh / y << endl;
        }
    }
    return g_ratio;
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

void MakeFinalPlots() {

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
        if (elec == 1)textn += "BornTT.root";
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

    g_re_final = ResbosFromRaj();
    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_ratio_phistar = CreateRatio(g_data_final, g_data_final, 1);
    TGraphAsymmErrors* g_ratio_mg_phistar = CreateRatio(g_data_final, g_mg_final, 0);
    TGraphAsymmErrors* g_ratio_ph_phistar = CreateRatio(g_data_final, g_ph_final, 0);
    if (g_re_final != 0)g_ratio_re_phistar = CreateRatio(g_data_final, g_re_final, 0);


    cout << "test 1" << endl;
    //PrintFinal(g_data_final, g_mg_final, g_ph_final, 0, g_re_final);
    //cout << "test 2" << endl;
    PlotFinal(g_data_final, g_mg_final, g_ph_final, g_ratio_phistar, g_ratio_mg_phistar, g_ratio_ph_phistar, 0, g_re_final, g_ratio_re_phistar);
}

// TH1D* GetHistMuon(std::string name){
//   TH1D* h_Muon=new TH1D("h_Comb","h_Comb",nphistar,phistarBins);
//   h_Muon->Sumw2();
//   TFile fm(name.c_str());
//   TH1F *h_temp = (TH1F*)fm.Get("hRecoClone1");
//   for (uint i=0; i<nphistar; i++){
//     h_Muon->SetBinContent(i+1,h_temp->GetBinContent(i+1));
//     h_Muon->SetBinError(i+1,h_temp->GetBinError(i+1));
//   }
//   return h_Muon;
// } 

void MakeFinalPlots2() {

    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_data_muon;
    TGraphAsymmErrors* g_data_elec;

    std::string textn = "Output/Data_Graph_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "PH_";
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "BornTT.root";
    if (elec == 2)textn += "Naked.root";
    std::cout << textn.c_str() << std::endl;
    TFile dg(textn.c_str());
    TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
    g_data_elec = CreateCopy(g_data_temp);
    textn = "Output/MC_Graph_";
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
    if (Type == "combined") {
    }
    if (Type == "elec") {
        textn = "Output/Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dm(textn.c_str());
        TGraphAsymmErrors* g_muon_temp = (TGraphAsymmErrors*) dm.Get("Graph");
        g_data_muon = CreateCopy(g_muon_temp);
    }
    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_ratio_elec = CreateRatio(g_mg_final, g_data_elec, 0);
    TGraphAsymmErrors* g_ratio_mg = CreateRatio(g_mg_final, g_mg_final, 1);
    TGraphAsymmErrors* g_ratio_muon = CreateRatio(g_mg_final, g_data_muon, 0);
    //
    //  PrintFinal(g_data_final,g_mg_final,g_ph_final);
    PlotFinal(g_mg_final, g_data_elec, g_data_muon, g_ratio_mg, g_ratio_elec, g_ratio_muon, 1);

    // textn="Comb_Hist_";
    // textn+=Tag;
    // if (doNorm) textn+="Norm_";
    // else        textn+="Abs_";
    // textn+="PH_";
    // textn+="Born.root";
    // TFile tr_comb(textn.c_str());
    // TH1D* h_comb_temp= (TH1D*)tr_comb.Get("h_Comb");
    // TGraphAsymmErrors* g_comb_temp=ConvertToTGraph(h_comb_temp);
    // TGraphAsymmErrors* g_comb_final=CreateCopy(g_comb_temp);
    // PrintFinal(g_comb_final,g_data_elec,g_data_muon,1);
}

