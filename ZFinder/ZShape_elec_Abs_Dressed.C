{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Tue May  5 15:41:25 2015) by ROOT version5.34/28
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",1922,204,800,872);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-3.924951,-0.1474366,0.7623949,4.905219);
   p1->SetFillColor(0);
   p1->SetBorderMode(0);
   p1->SetBorderSize(0);
   p1->SetLogx();
   p1->SetLogy();
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.06);
   p1->SetTopMargin(0.06);
   p1->SetBottomMargin(0.01);
   p1->SetFrameBorderMode(0);
   p1->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph0");
   grae->SetTitle("");
   grae->SetPoint(0,0.002,4297.204);
   grae->SetPointError(0,0.002,0.002,0,0);
   grae->SetPoint(1,0.006,4302.073);
   grae->SetPointError(1,0.002,0.002,0,0);
   grae->SetPoint(2,0.01,4189.518);
   grae->SetPointError(2,0.002,0.002,0,0);
   grae->SetPoint(3,0.014,4137.143);
   grae->SetPointError(3,0.002,0.002,0,0);
   grae->SetPoint(4,0.018,3950.681);
   grae->SetPointError(4,0.002,0.002,0,0);
   grae->SetPoint(5,0.022,3816.286);
   grae->SetPointError(5,0.002,0.002,0,0);
   grae->SetPoint(6,0.0265,3631.916);
   grae->SetPointError(6,0.0025,0.0025,0,0);
   grae->SetPoint(7,0.0315,3404.196);
   grae->SetPointError(7,0.0025,0.0025,0,0);
   grae->SetPoint(8,0.0365,3181.777);
   grae->SetPointError(8,0.0025,0.0025,0,0);
   grae->SetPoint(9,0.042,2971.462);
   grae->SetPointError(9,0.003,0.003,0,0);
   grae->SetPoint(10,0.0485,2724.068);
   grae->SetPointError(10,0.0035,0.0035,0,0);
   grae->SetPoint(11,0.0545,2521.256);
   grae->SetPointError(11,0.0025,0.0025,0,0);
   grae->SetPoint(12,0.0605,2335.465);
   grae->SetPointError(12,0.0035,0.0035,0,0);
   grae->SetPoint(13,0.068,2098.982);
   grae->SetPointError(13,0.004,0.004,0,0);
   grae->SetPoint(14,0.0765,1903.899);
   grae->SetPointError(14,0.0045,0.0045,0,0);
   grae->SetPoint(15,0.086,1694.217);
   grae->SetPointError(15,0.005,0.005,0,0);
   grae->SetPoint(16,0.0965,1498.413);
   grae->SetPointError(16,0.0055,0.0055,0,0);
   grae->SetPoint(17,0.108,1319.693);
   grae->SetPointError(17,0.006,0.006,0,0);
   grae->SetPoint(18,0.121,1152.376);
   grae->SetPointError(18,0.007,0.007,0,0);
   grae->SetPoint(19,0.1365,984.1828);
   grae->SetPointError(19,0.0085,0.0085,0,0);
   grae->SetPoint(20,0.155,826.6441);
   grae->SetPointError(20,0.01,0.01,0,0);
   grae->SetPoint(21,0.177,677.7919);
   grae->SetPointError(21,0.012,0.012,0,0);
   grae->SetPoint(22,0.204,536.8369);
   grae->SetPointError(22,0.015,0.015,0,0);
   grae->SetPoint(23,0.2385,415.3886);
   grae->SetPointError(23,0.0195,0.0195,0,0);
   grae->SetPoint(24,0.285,304.1331);
   grae->SetPointError(24,0.027,0.027,0,0);
   grae->SetPoint(25,0.3515,203.7893);
   grae->SetPointError(25,0.0395,0.0395,0,0);
   grae->SetPoint(26,0.4575,119.952);
   grae->SetPointError(26,0.0665,0.0665,0,0);
   grae->SetPoint(27,0.6095,62.36899);
   grae->SetPointError(27,0.0855,0.0855,0,0);
   grae->SetPoint(28,0.8065,31.56173);
   grae->SetPointError(28,0.1115,0.1115,0,0);
   grae->SetPoint(29,1.0355,16.11513);
   grae->SetPointError(29,0.1175,0.1175,0,0);
   grae->SetPoint(30,1.3245,8.491314);
   grae->SetPointError(30,0.1715,0.1715,0,0);
   grae->SetPoint(31,1.7215,3.980523);
   grae->SetPointError(31,0.2255,0.2255,0,0);
   grae->SetPoint(32,2.2345,1.96241);
   grae->SetPointError(32,0.2875,0.2875,0,0);
   grae->SetPoint(33,2.8995,1.001341);
   grae->SetPointError(33,0.3775,0.3775,0,0);
   
   TH1F *Graph_Graph21 = new TH1F("Graph_Graph21","",100,0.0006,3.604633);
   Graph_Graph21->SetMinimum(0.8);
   Graph_Graph21->SetMaximum(40000);
   Graph_Graph21->SetDirectory(0);
   Graph_Graph21->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph21->SetLineColor(ci);
   Graph_Graph21->GetXaxis()->SetRange(1,84);
   Graph_Graph21->GetXaxis()->SetLabelFont(42);
   Graph_Graph21->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph21->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph21->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph21->GetXaxis()->SetTitleFont(42);
   Graph_Graph21->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi* (pb)");
   Graph_Graph21->GetYaxis()->SetLabelFont(42);
   Graph_Graph21->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph21->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph21->GetYaxis()->SetTitleOffset(1.05);
   Graph_Graph21->GetYaxis()->SetTitleFont(42);
   Graph_Graph21->GetZaxis()->SetLabelFont(42);
   Graph_Graph21->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph21->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph21->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph21);
   
   grae->Draw("a2");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.002,4154.677);
   grae->SetPointError(0,0,0,137.9876,137.9876);
   grae->SetPoint(1,0.006,4083.134);
   grae->SetPointError(1,0,0,135.6147,135.6147);
   grae->SetPoint(2,0.01,4037.129);
   grae->SetPointError(2,0,0,134.1118,134.1118);
   grae->SetPoint(3,0.014,3968.948);
   grae->SetPointError(3,0,0,131.8459,131.8459);
   grae->SetPoint(4,0.018,3879.055);
   grae->SetPointError(4,0,0,128.8688,128.8688);
   grae->SetPoint(5,0.022,3733.707);
   grae->SetPointError(5,0,0,124.072,124.072);
   grae->SetPoint(6,0.0265,3585.588);
   grae->SetPointError(6,0,0,119.0733,119.0733);
   grae->SetPoint(7,0.0315,3396.305);
   grae->SetPointError(7,0,0,112.8175,112.8175);
   grae->SetPoint(8,0.0365,3191.902);
   grae->SetPointError(8,0,0,106.0526,106.0526);
   grae->SetPoint(9,0.042,2986.599);
   grae->SetPointError(9,0,0,99.20193,99.20193);
   grae->SetPoint(10,0.0485,2749.509);
   grae->SetPointError(10,0,0,91.324,91.324);
   grae->SetPoint(11,0.0545,2513.574);
   grae->SetPointError(11,0,0,83.61512,83.61512);
   grae->SetPoint(12,0.0605,2335.24);
   grae->SetPointError(12,0,0,77.60991,77.60991);
   grae->SetPoint(13,0.068,2104.109);
   grae->SetPointError(13,0,0,69.92933,69.92933);
   grae->SetPoint(14,0.0765,1883.177);
   grae->SetPointError(14,0,0,62.59515,62.59515);
   grae->SetPoint(15,0.086,1677.629);
   grae->SetPointError(15,0,0,55.77028,55.77028);
   grae->SetPoint(16,0.0965,1471.188);
   grae->SetPointError(16,0,0,48.91121,48.91121);
   grae->SetPoint(17,0.108,1288.342);
   grae->SetPointError(17,0,0,42.84817,42.84817);
   grae->SetPoint(18,0.121,1118.149);
   grae->SetPointError(18,0,0,37.18851,37.18851);
   grae->SetPoint(19,0.1365,958.1795);
   grae->SetPointError(19,0,0,31.86742,31.86742);
   grae->SetPoint(20,0.155,797.3958);
   grae->SetPointError(20,0,0,26.52335,26.52335);
   grae->SetPoint(21,0.177,648.1036);
   grae->SetPointError(21,0,0,21.55964,21.55964);
   grae->SetPoint(22,0.204,516.6624);
   grae->SetPointError(22,0,0,17.18805,17.18805);
   grae->SetPoint(23,0.2385,394.2691);
   grae->SetPointError(23,0,0,13.11855,13.11855);
   grae->SetPoint(24,0.285,286.5257);
   grae->SetPointError(24,0,0,9.533377,9.533377);
   grae->SetPoint(25,0.3515,192.0471);
   grae->SetPointError(25,0,0,6.390185,6.390185);
   grae->SetPoint(26,0.4575,111.5901);
   grae->SetPointError(26,0,0,3.71396,3.71396);
   grae->SetPoint(27,0.6095,58.56676);
   grae->SetPointError(27,0,0,1.95222,1.95222);
   grae->SetPoint(28,0.8065,29.28843);
   grae->SetPointError(28,0,0,0.9789815,0.9789815);
   grae->SetPoint(29,1.0355,15.23137);
   grae->SetPointError(29,0,0,0.5125173,0.5125173);
   grae->SetPoint(30,1.3245,8.027508);
   grae->SetPointError(30,0,0,0.271293,0.271293);
   grae->SetPoint(31,1.7215,3.986646);
   grae->SetPointError(31,0,0,0.1360319,0.1360319);
   grae->SetPoint(32,2.2345,1.934107);
   grae->SetPointError(32,0,0,0.06713188,0.06713188);
   grae->SetPoint(33,2.8995,0.9723212);
   grae->SetPointError(33,0,0,0.0344906,0.0344906);
   
   TH1F *Graph_Graph22 = new TH1F("Graph_Graph22","Graph",100,0.0018,3.18925);
   Graph_Graph22->SetMinimum(0.8440475);
   Graph_Graph22->SetMaximum(4721.837);
   Graph_Graph22->SetDirectory(0);
   Graph_Graph22->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph22->SetLineColor(ci);
   Graph_Graph22->GetXaxis()->SetLabelFont(42);
   Graph_Graph22->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph22->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph22->GetXaxis()->SetTitleFont(42);
   Graph_Graph22->GetYaxis()->SetLabelFont(42);
   Graph_Graph22->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph22->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph22->GetYaxis()->SetTitleFont(42);
   Graph_Graph22->GetZaxis()->SetLabelFont(42);
   Graph_Graph22->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph22->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph22->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph22);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetPoint(0,0.002,3871.012);
   grae->SetPointError(0,0,0,104.1892,105.1251);
   grae->SetPoint(1,0.006,3880.854);
   grae->SetPointError(1,0,0,104.5812,105.3478);
   grae->SetPoint(2,0.01,3805.844);
   grae->SetPointError(2,0,0,101.7995,101.92);
   grae->SetPoint(3,0.014,3746.193);
   grae->SetPointError(3,0,0,101.6182,102.9882);
   grae->SetPoint(4,0.018,3660.446);
   grae->SetPointError(4,0,0,98.75869,99.38061);
   grae->SetPoint(5,0.022,3641.954);
   grae->SetPointError(5,0,0,98.63588,99.92567);
   grae->SetPoint(6,0.0265,3517.822);
   grae->SetPointError(6,0,0,94.86181,95.87792);
   grae->SetPoint(7,0.0315,3393.033);
   grae->SetPointError(7,0,0,90.99152,91.69483);
   grae->SetPoint(8,0.0365,3228.637);
   grae->SetPointError(8,0,0,86.78341,87.81123);
   grae->SetPoint(9,0.042,3071.199);
   grae->SetPointError(9,0,0,82.07134,82.59198);
   grae->SetPoint(10,0.0485,2877.875);
   grae->SetPointError(10,0,0,76.90605,77.67196);
   grae->SetPoint(11,0.0545,2661.719);
   grae->SetPointError(11,0,0,72.14259,73.26733);
   grae->SetPoint(12,0.0605,2479.24);
   grae->SetPointError(12,0,0,66.67986,67.81527);
   grae->SetPoint(13,0.068,2257.32);
   grae->SetPointError(13,0,0,60.78242,62.05906);
   grae->SetPoint(14,0.0765,2056.613);
   grae->SetPointError(14,0,0,54.94883,55.66435);
   grae->SetPoint(15,0.086,1830.73);
   grae->SetPointError(15,0,0,48.83455,49.4012);
   grae->SetPoint(16,0.0965,1588.761);
   grae->SetPointError(16,0,0,42.70708,43.73673);
   grae->SetPoint(17,0.108,1391.693);
   grae->SetPointError(17,0,0,37.20393,37.80053);
   grae->SetPoint(18,0.121,1198.248);
   grae->SetPointError(18,0,0,31.99189,32.62879);
   grae->SetPoint(19,0.1365,1008.127);
   grae->SetPointError(19,0,0,26.92092,27.35844);
   grae->SetPoint(20,0.155,824.3385);
   grae->SetPointError(20,0,0,21.9491,22.33665);
   grae->SetPoint(21,0.177,675.7394);
   grae->SetPointError(21,0,0,18.07732,18.47589);
   grae->SetPoint(22,0.204,535.9388);
   grae->SetPointError(22,0,0,14.26476,14.63632);
   grae->SetPoint(23,0.2385,412.1994);
   grae->SetPointError(23,0,0,10.8855,11.02243);
   grae->SetPoint(24,0.285,298.3513);
   grae->SetPointError(24,0,0,7.8578,7.948862);
   grae->SetPoint(25,0.3515,197.186);
   grae->SetPointError(25,0,0,5.139431,5.154013);
   grae->SetPoint(26,0.4575,114.161);
   grae->SetPointError(26,0,0,2.953508,2.937482);
   grae->SetPoint(27,0.6095,57.50506);
   grae->SetPointError(27,0,0,1.516287,1.521268);
   grae->SetPoint(28,0.8065,28.27642);
   grae->SetPointError(28,0,0,0.752416,0.7486204);
   grae->SetPoint(29,1.0355,14.12064);
   grae->SetPointError(29,0,0,0.3942329,0.3933234);
   grae->SetPoint(30,1.3245,7.091961);
   grae->SetPointError(30,0,0,0.2057675,0.2060806);
   grae->SetPoint(31,1.7215,3.295981);
   grae->SetPointError(31,0,0,0.1044981,0.1055117);
   grae->SetPoint(32,2.2345,1.597593);
   grae->SetPointError(32,0,0,0.05677508,0.05782978);
   grae->SetPoint(33,2.8995,0.8524869);
   grae->SetPointError(33,0,0,0.03356333,0.03435786);
   
   TH1F *Graph_Graph23 = new TH1F("Graph_Graph23","Graph",100,0.0018,3.18925);
   Graph_Graph23->SetMinimum(0.7370313);
   Graph_Graph23->SetMaximum(4384.74);
   Graph_Graph23->SetDirectory(0);
   Graph_Graph23->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph23->SetLineColor(ci);
   Graph_Graph23->GetXaxis()->SetLabelFont(42);
   Graph_Graph23->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetXaxis()->SetTitleFont(42);
   Graph_Graph23->GetYaxis()->SetLabelFont(42);
   Graph_Graph23->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetYaxis()->SetTitleFont(42);
   Graph_Graph23->GetZaxis()->SetLabelFont(42);
   Graph_Graph23->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph23);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph3");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.002,4297.204);
   grae->SetPointError(0,0,0,116.6324,116.7365);
   grae->SetPoint(1,0.006,4302.073);
   grae->SetPointError(1,0,0,116.7317,116.9427);
   grae->SetPoint(2,0.01,4189.518);
   grae->SetPointError(2,0,0,113.7089,113.4655);
   grae->SetPoint(3,0.014,4137.143);
   grae->SetPointError(3,0,0,112.2763,112.5655);
   grae->SetPoint(4,0.018,3950.681);
   grae->SetPointError(4,0,0,108.3704,108.5136);
   grae->SetPoint(5,0.022,3816.286);
   grae->SetPointError(5,0,0,103.1852,103.3199);
   grae->SetPoint(6,0.0265,3631.916);
   grae->SetPointError(6,0,0,98.74501,98.69488);
   grae->SetPoint(7,0.0315,3404.196);
   grae->SetPointError(7,0,0,92.40891,92.62502);
   grae->SetPoint(8,0.0365,3181.777);
   grae->SetPointError(8,0,0,86.60188,86.69261);
   grae->SetPoint(9,0.042,2971.462);
   grae->SetPointError(9,0,0,80.63014,80.58662);
   grae->SetPoint(10,0.0485,2724.068);
   grae->SetPointError(10,0,0,73.71232,73.81946);
   grae->SetPoint(11,0.0545,2521.256);
   grae->SetPointError(11,0,0,68.82099,68.73831);
   grae->SetPoint(12,0.0605,2335.465);
   grae->SetPointError(12,0,0,63.24248,63.29069);
   grae->SetPoint(13,0.068,2098.982);
   grae->SetPointError(13,0,0,57.3449,57.36593);
   grae->SetPoint(14,0.0765,1903.899);
   grae->SetPointError(14,0,0,51.52481,51.6476);
   grae->SetPoint(15,0.086,1694.217);
   grae->SetPointError(15,0,0,45.93894,45.95309);
   grae->SetPoint(16,0.0965,1498.413);
   grae->SetPointError(16,0,0,40.74544,40.78842);
   grae->SetPoint(17,0.108,1319.693);
   grae->SetPointError(17,0,0,35.83975,35.89686);
   grae->SetPoint(18,0.121,1152.376);
   grae->SetPointError(18,0,0,31.25338,31.30127);
   grae->SetPoint(19,0.1365,984.1828);
   grae->SetPointError(19,0,0,26.49915,26.55286);
   grae->SetPoint(20,0.155,826.6441);
   grae->SetPointError(20,0,0,22.3267,22.35187);
   grae->SetPoint(21,0.177,677.7919);
   grae->SetPointError(21,0,0,18.27912,18.30116);
   grae->SetPoint(22,0.204,536.8369);
   grae->SetPointError(22,0,0,14.49433,14.51594);
   grae->SetPoint(23,0.2385,415.3886);
   grae->SetPointError(23,0,0,11.26514,11.28101);
   grae->SetPoint(24,0.285,304.1331);
   grae->SetPointError(24,0,0,8.215319,8.229753);
   grae->SetPoint(25,0.3515,203.7893);
   grae->SetPointError(25,0,0,5.506072,5.509655);
   grae->SetPoint(26,0.4575,119.952);
   grae->SetPointError(26,0,0,3.253207,3.246873);
   grae->SetPoint(27,0.6095,62.36899);
   grae->SetPointError(27,0,0,1.708501,1.705341);
   grae->SetPoint(28,0.8065,31.56173);
   grae->SetPointError(28,0,0,0.8649183,0.8624879);
   grae->SetPoint(29,1.0355,16.11513);
   grae->SetPointError(29,0,0,0.4531494,0.453447);
   grae->SetPoint(30,1.3245,8.491314);
   grae->SetPointError(30,0,0,0.2443236,0.2421763);
   grae->SetPoint(31,1.7215,3.980523);
   grae->SetPointError(31,0,0,0.1178618,0.1172956);
   grae->SetPoint(32,2.2345,1.96241);
   grae->SetPointError(32,0,0,0.06091066,0.0606875);
   grae->SetPoint(33,2.8995,1.001341);
   grae->SetPointError(33,0,0,0.03279666,0.03274627);
   
   TH1F *Graph_Graph24 = new TH1F("Graph_Graph24","Graph",100,0.0018,3.18925);
   Graph_Graph24->SetMinimum(0.8716902);
   Graph_Graph24->SetMaximum(4860.82);
   Graph_Graph24->SetDirectory(0);
   Graph_Graph24->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph24->SetLineColor(ci);
   Graph_Graph24->GetXaxis()->SetLabelFont(42);
   Graph_Graph24->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph24->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph24->GetXaxis()->SetTitleFont(42);
   Graph_Graph24->GetYaxis()->SetLabelFont(42);
   Graph_Graph24->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph24->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph24->GetYaxis()->SetTitleFont(42);
   Graph_Graph24->GetZaxis()->SetLabelFont(42);
   Graph_Graph24->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph24->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph24->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph24);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.53,0.72,0.85,0.91,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph3","2012 data","PEF");

   ci = TColor::GetColor("#999999");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Z #rightarrow ee MadGraph","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Z #rightarrow ee Powheg","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.745,0.95,"19.7 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.895,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.2,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.13,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.06,"60 GeV < M_{ee} < 120 GeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
   p1->Modified();
   FinalPhiTot->cd();
  
// ------------>Primitives in pad: p2
   p2 = new TPad("p2", "p2",0,0,1,0.2777778);
   p2->Draw();
   p2->cd();
   p2->Range(-3.924951,0.4735484,0.7623949,1.247742);
   p2->SetFillColor(0);
   p2->SetBorderMode(0);
   p2->SetBorderSize(0);
   p2->SetLogx();
   p2->SetLeftMargin(0.15);
   p2->SetRightMargin(0.06);
   p2->SetTopMargin(0.01);
   p2->SetBottomMargin(0.37);
   p2->SetFrameBorderMode(0);
   p2->SetFrameBorderMode(0);
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph0");
   grae->SetTitle("");

   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetPoint(0,0.002,1);
   grae->SetPointError(0,0.002,0.002,0.02714147,0.02716569);
   grae->SetPoint(1,0.006,1);
   grae->SetPointError(1,0.002,0.002,0.02713383,0.02718287);
   grae->SetPoint(2,0.01,1);
   grae->SetPointError(2,0.002,0.002,0.02714129,0.02708319);
   grae->SetPoint(3,0.014,1);
   grae->SetPointError(3,0.002,0.002,0.02713861,0.02720851);
   grae->SetPoint(4,0.018,1);
   grae->SetPointError(4,0.002,0.002,0.0274308,0.02746705);
   grae->SetPoint(5,0.022,1);
   grae->SetPointError(5,0.002,0.002,0.02703811,0.02707342);
   grae->SetPoint(6,0.0265,1);
   grae->SetPointError(6,0.0025,0.0025,0.02718813,0.02717433);
   grae->SetPoint(7,0.0315,1);
   grae->SetPointError(7,0.0025,0.0025,0.02714559,0.02720908);
   grae->SetPoint(8,0.0365,1);
   grae->SetPointError(8,0.0025,0.0025,0.02721808,0.0272466);
   grae->SetPoint(9,0.042,1);
   grae->SetPointError(9,0.003,0.003,0.02713484,0.02712019);
   grae->SetPoint(10,0.0485,1);
   grae->SetPointError(10,0.0035,0.0035,0.02705965,0.02709898);
   grae->SetPoint(11,0.0545,1);
   grae->SetPointError(11,0.0025,0.0025,0.02729631,0.02726352);
   grae->SetPoint(12,0.0605,1);
   grae->SetPointError(12,0.0035,0.0035,0.02707919,0.02709983);
   grae->SetPoint(13,0.068,1);
   grae->SetPointError(13,0.004,0.004,0.02732034,0.02733035);
   grae->SetPoint(14,0.0765,1);
   grae->SetPointError(14,0.0045,0.0045,0.02706278,0.02712728);
   grae->SetPoint(15,0.086,1);
   grae->SetPointError(15,0.005,0.005,0.02711514,0.02712349);
   grae->SetPoint(16,0.0965,1);
   grae->SetPointError(16,0.0055,0.0055,0.0271924,0.02722108);
   grae->SetPoint(17,0.108,1);
   grae->SetPointError(17,0.006,0.006,0.02715763,0.02720091);
   grae->SetPoint(18,0.121,1);
   grae->SetPointError(18,0.007,0.007,0.02712081,0.02716237);
   grae->SetPoint(19,0.1365,1);
   grae->SetPointError(19,0.0085,0.0085,0.02692502,0.0269796);
   grae->SetPoint(20,0.155,1);
   grae->SetPointError(20,0.01,0.01,0.02700884,0.02703929);
   grae->SetPoint(21,0.177,1);
   grae->SetPointError(21,0.012,0.012,0.02696863,0.02700115);
   grae->SetPoint(22,0.204,1);
   grae->SetPointError(22,0.015,0.015,0.0269995,0.02703976);
   grae->SetPoint(23,0.2385,1);
   grae->SetPointError(23,0.0195,0.0195,0.02711952,0.02715774);
   grae->SetPoint(24,0.285,1);
   grae->SetPointError(24,0.027,0.027,0.02701225,0.02705971);
   grae->SetPoint(25,0.3515,1);
   grae->SetPointError(25,0.0395,0.0395,0.02701845,0.02703603);
   grae->SetPoint(26,0.4575,1);
   grae->SetPointError(26,0.0665,0.0665,0.02712092,0.02706811);
   grae->SetPoint(27,0.6095,1);
   grae->SetPointError(27,0.0855,0.0855,0.02739343,0.02734277);
   grae->SetPoint(28,0.8065,1);
   grae->SetPointError(28,0.1115,0.1115,0.02740402,0.02732701);
   grae->SetPoint(29,1.0355,1);
   grae->SetPointError(29,0.1175,0.1175,0.0281195,0.02813796);
   grae->SetPoint(30,1.3245,1);
   grae->SetPointError(30,0.1715,0.1715,0.02877336,0.02852047);
   grae->SetPoint(31,1.7215,1);
   grae->SetPointError(31,0.2255,0.2255,0.02960962,0.02946737);
   grae->SetPoint(32,2.2345,1);
   grae->SetPointError(32,0.2875,0.2875,0.03103871,0.03092499);
   grae->SetPoint(33,2.8995,1);
   grae->SetPointError(33,0.3775,0.3775,0.03275273,0.03270241);
   
   TH1F *Graph_Graph25 = new TH1F("Graph_Graph25","",100,0.0006,3.604633);
   Graph_Graph25->SetMinimum(0.76);
   Graph_Graph25->SetMaximum(1.24);
   Graph_Graph25->SetDirectory(0);
   Graph_Graph25->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph25->SetLineColor(ci);
   Graph_Graph25->GetXaxis()->SetTitle("#phi*");
   Graph_Graph25->GetXaxis()->SetRange(1,84);
   Graph_Graph25->GetXaxis()->SetLabelFont(42);
   Graph_Graph25->GetXaxis()->SetLabelSize(0.12);
   Graph_Graph25->GetXaxis()->SetTitleSize(0.12);
   Graph_Graph25->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph25->GetXaxis()->SetTitleFont(42);
   Graph_Graph25->GetYaxis()->SetTitle("MC/Data   ");
   Graph_Graph25->GetYaxis()->SetNdivisions(503);
   Graph_Graph25->GetYaxis()->SetLabelFont(42);
   Graph_Graph25->GetYaxis()->SetLabelSize(0.12);
   Graph_Graph25->GetYaxis()->SetTitleSize(0.12);
   Graph_Graph25->GetYaxis()->SetTitleOffset(0.45);
   Graph_Graph25->GetYaxis()->SetTitleFont(42);
   Graph_Graph25->GetZaxis()->SetLabelFont(42);
   Graph_Graph25->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph25->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph25->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph25);
   
   grae->Draw("ae2");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.002,0.9668325);
   grae->SetPointError(0,0,0,0.02714147,0.02716569);
   grae->SetPoint(1,0.006,0.9491084);
   grae->SetPointError(1,0,0,0.02713383,0.02718287);
   grae->SetPoint(2,0.01,0.963626);
   grae->SetPointError(2,0,0,0.02714129,0.02708319);
   grae->SetPoint(3,0.014,0.9593452);
   grae->SetPointError(3,0,0,0.02713861,0.02720851);
   grae->SetPoint(4,0.018,0.98187);
   grae->SetPointError(4,0,0,0.0274308,0.02746705);
   grae->SetPoint(5,0.022,0.9783616);
   grae->SetPointError(5,0,0,0.02703811,0.02707342);
   grae->SetPoint(6,0.0265,0.9872442);
   grae->SetPointError(6,0,0,0.02718813,0.02717433);
   grae->SetPoint(7,0.0315,0.9976821);
   grae->SetPointError(7,0,0,0.02714559,0.02720908);
   grae->SetPoint(8,0.0365,1.003182);
   grae->SetPointError(8,0,0,0.02721808,0.0272466);
   grae->SetPoint(9,0.042,1.005094);
   grae->SetPointError(9,0,0,0.02713484,0.02712019);
   grae->SetPoint(10,0.0485,1.00934);
   grae->SetPointError(10,0,0,0.02705965,0.02709898);
   grae->SetPoint(11,0.0545,0.9969529);
   grae->SetPointError(11,0,0,0.02729631,0.02726352);
   grae->SetPoint(12,0.0605,0.9999038);
   grae->SetPointError(12,0,0,0.02707919,0.02709983);
   grae->SetPoint(13,0.068,1.002443);
   grae->SetPointError(13,0,0,0.02732034,0.02733035);
   grae->SetPoint(14,0.0765,0.9891161);
   grae->SetPointError(14,0,0,0.02706278,0.02712728);
   grae->SetPoint(15,0.086,0.9902089);
   grae->SetPointError(15,0,0,0.02711514,0.02712349);
   grae->SetPoint(16,0.0965,0.981831);
   grae->SetPointError(16,0,0,0.0271924,0.02722108);
   grae->SetPoint(17,0.108,0.9762431);
   grae->SetPointError(17,0,0,0.02715763,0.02720091);
   grae->SetPoint(18,0.121,0.9702988);
   grae->SetPointError(18,0,0,0.02712081,0.02716237);
   grae->SetPoint(19,0.1365,0.9735788);
   grae->SetPointError(19,0,0,0.02692502,0.0269796);
   grae->SetPoint(20,0.155,0.964618);
   grae->SetPointError(20,0,0,0.02700884,0.02703929);
   grae->SetPoint(21,0.177,0.9561985);
   grae->SetPointError(21,0,0,0.02696863,0.02700115);
   grae->SetPoint(22,0.204,0.9624196);
   grae->SetPointError(22,0,0,0.0269995,0.02703976);
   grae->SetPoint(23,0.2385,0.9491574);
   grae->SetPointError(23,0,0,0.02711952,0.02715774);
   grae->SetPoint(24,0.285,0.9421065);
   grae->SetPointError(24,0,0,0.02701225,0.02705971);
   grae->SetPoint(25,0.3515,0.9423807);
   grae->SetPointError(25,0,0,0.02701845,0.02703603);
   grae->SetPoint(26,0.4575,0.9302901);
   grae->SetPointError(26,0,0,0.02712092,0.02706811);
   grae->SetPoint(27,0.6095,0.9390366);
   grae->SetPointError(27,0,0,0.02739343,0.02734277);
   grae->SetPoint(28,0.8065,0.927973);
   grae->SetPointError(28,0,0,0.02740402,0.02732701);
   grae->SetPoint(29,1.0355,0.9451597);
   grae->SetPointError(29,0,0,0.0281195,0.02813796);
   grae->SetPoint(30,1.3245,0.9453787);
   grae->SetPointError(30,0,0,0.02877336,0.02852047);
   grae->SetPoint(31,1.7215,1.001538);
   grae->SetPointError(31,0,0,0.02960962,0.02946737);
   grae->SetPoint(32,2.2345,0.9855775);
   grae->SetPointError(32,0,0,0.03103871,0.03092499);
   grae->SetPoint(33,2.8995,0.9710188);
   grae->SetPointError(33,0,0,0.03275273,0.03270241);
   
   TH1F *Graph_Graph26 = new TH1F("Graph_Graph26","Graph",100,0.0018,3.18925);
   Graph_Graph26->SetMinimum(0.886982);
   Graph_Graph26->SetMaximum(1.050026);
   Graph_Graph26->SetDirectory(0);
   Graph_Graph26->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph26->SetLineColor(ci);
   Graph_Graph26->GetXaxis()->SetLabelFont(42);
   Graph_Graph26->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph26->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph26->GetXaxis()->SetTitleFont(42);
   Graph_Graph26->GetYaxis()->SetLabelFont(42);
   Graph_Graph26->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph26->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph26->GetYaxis()->SetTitleFont(42);
   Graph_Graph26->GetZaxis()->SetLabelFont(42);
   Graph_Graph26->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph26->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph26->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph26);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetPoint(0,0.002,0.9008209);
   grae->SetPointError(0,0,0,0.02714147,0.02716569);
   grae->SetPoint(1,0.006,0.9020893);
   grae->SetPointError(1,0,0,0.02713383,0.02718287);
   grae->SetPoint(2,0.01,0.9084203);
   grae->SetPointError(2,0,0,0.02714129,0.02708319);
   grae->SetPoint(3,0.014,0.9055025);
   grae->SetPointError(3,0,0,0.02713861,0.02720851);
   grae->SetPoint(4,0.018,0.9265354);
   grae->SetPointError(4,0,0,0.0274308,0.02746705);
   grae->SetPoint(5,0.022,0.9543191);
   grae->SetPointError(5,0,0,0.02703811,0.02707342);
   grae->SetPoint(6,0.0265,0.9685856);
   grae->SetPointError(6,0,0,0.02718813,0.02717433);
   grae->SetPoint(7,0.0315,0.9967209);
   grae->SetPointError(7,0,0,0.02714559,0.02720908);
   grae->SetPoint(8,0.0365,1.014727);
   grae->SetPointError(8,0,0,0.02721808,0.0272466);
   grae->SetPoint(9,0.042,1.033565);
   grae->SetPointError(9,0,0,0.02713484,0.02712019);
   grae->SetPoint(10,0.0485,1.056463);
   grae->SetPointError(10,0,0,0.02705965,0.02709898);
   grae->SetPoint(11,0.0545,1.055712);
   grae->SetPointError(11,0,0,0.02729631,0.02726352);
   grae->SetPoint(12,0.0605,1.061562);
   grae->SetPointError(12,0,0,0.02707919,0.02709983);
   grae->SetPoint(13,0.068,1.075435);
   grae->SetPointError(13,0,0,0.02732034,0.02733035);
   grae->SetPoint(14,0.0765,1.080211);
   grae->SetPointError(14,0,0,0.02706278,0.02712728);
   grae->SetPoint(15,0.086,1.080576);
   grae->SetPointError(15,0,0,0.02711514,0.02712349);
   grae->SetPoint(16,0.0965,1.060296);
   grae->SetPointError(16,0,0,0.0271924,0.02722108);
   grae->SetPoint(17,0.108,1.054558);
   grae->SetPointError(17,0,0,0.02715763,0.02720091);
   grae->SetPoint(18,0.121,1.039807);
   grae->SetPointError(18,0,0,0.02712081,0.02716237);
   grae->SetPoint(19,0.1365,1.024329);
   grae->SetPointError(19,0,0,0.02692502,0.0269796);
   grae->SetPoint(20,0.155,0.9972109);
   grae->SetPointError(20,0,0,0.02700884,0.02703929);
   grae->SetPoint(21,0.177,0.9969719);
   grae->SetPointError(21,0,0,0.02696863,0.02700115);
   grae->SetPoint(22,0.204,0.9983271);
   grae->SetPointError(22,0,0,0.0269995,0.02703976);
   grae->SetPoint(23,0.2385,0.9923225);
   grae->SetPointError(23,0,0,0.02711952,0.02715774);
   grae->SetPoint(24,0.285,0.9809895);
   grae->SetPointError(24,0,0,0.02701225,0.02705971);
   grae->SetPoint(25,0.3515,0.9675976);
   grae->SetPointError(25,0,0,0.02701845,0.02703603);
   grae->SetPoint(26,0.4575,0.9517224);
   grae->SetPointError(26,0,0,0.02712092,0.02706811);
   grae->SetPoint(27,0.6095,0.9220137);
   grae->SetPointError(27,0,0,0.02739343,0.02734277);
   grae->SetPoint(28,0.8065,0.8959085);
   grae->SetPointError(28,0,0,0.02740402,0.02732701);
   grae->SetPoint(29,1.0355,0.8762349);
   grae->SetPointError(29,0,0,0.0281195,0.02813796);
   grae->SetPoint(30,1.3245,0.8352018);
   grae->SetPointError(30,0,0,0.02877336,0.02852047);
   grae->SetPoint(31,1.7215,0.8280271);
   grae->SetPointError(31,0,0,0.02960962,0.02946737);
   grae->SetPoint(32,2.2345,0.8140976);
   grae->SetPointError(32,0,0,0.03103871,0.03092499);
   grae->SetPoint(33,2.8995,0.8513451);
   grae->SetPointError(33,0,0,0.03275273,0.03270241);
   
   TH1F *Graph_Graph27 = new TH1F("Graph_Graph27","Graph",100,0.0018,3.18925);
   Graph_Graph27->SetMinimum(0.7505948);
   Graph_Graph27->SetMaximum(1.140163);
   Graph_Graph27->SetDirectory(0);
   Graph_Graph27->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph27->SetLineColor(ci);
   Graph_Graph27->GetXaxis()->SetLabelFont(42);
   Graph_Graph27->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph27->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph27->GetXaxis()->SetTitleFont(42);
   Graph_Graph27->GetYaxis()->SetLabelFont(42);
   Graph_Graph27->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph27->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph27->GetYaxis()->SetTitleFont(42);
   Graph_Graph27->GetZaxis()->SetLabelFont(42);
   Graph_Graph27->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph27->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph27->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph27);
   
   grae->Draw("pe");
   p2->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->SetSelected(FinalPhiTot);
}
