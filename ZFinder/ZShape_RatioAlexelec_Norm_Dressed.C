{
//=========Macro generated from canvas: FinalPhiRatio/FinalPhiRatio
//=========  (Thu May 21 03:08:07 2015) by ROOT version5.34/18
   TCanvas *FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio",51,66,800,900);
   FinalPhiRatio->Range(-3.684724,0.625,0.9440295,1.375);
   FinalPhiRatio->SetFillColor(0);
   FinalPhiRatio->SetBorderMode(0);
   FinalPhiRatio->SetBorderSize(2);
   FinalPhiRatio->SetLogx();
   FinalPhiRatio->SetFrameBorderMode(0);
   FinalPhiRatio->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph0");
   grae->SetTitle("");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetPoint(0,0.002,1);
   grae->SetPointError(0,0.002,0.002,0.004420763,0.004298288);
   grae->SetPoint(1,0.006,1);
   grae->SetPointError(1,0.002,0.002,0.004801765,0.00496134);
   grae->SetPoint(2,0.01,1);
   grae->SetPointError(2,0.002,0.002,0.004808034,0.004665582);
   grae->SetPoint(3,0.014,1);
   grae->SetPointError(3,0.002,0.002,0.004803764,0.004992469);
   grae->SetPoint(4,0.018,1);
   grae->SetPointError(4,0.002,0.002,0.004993568,0.005127593);
   grae->SetPoint(5,0.022,1);
   grae->SetPointError(5,0.002,0.002,0.004937079,0.005079358);
   grae->SetPoint(6,0.0265,1);
   grae->SetPointError(6,0.0025,0.0025,0.004505659,0.00437283);
   grae->SetPoint(7,0.0315,1);
   grae->SetPointError(7,0.0025,0.0025,0.004493937,0.004597712);
   grae->SetPoint(8,0.0365,1);
   grae->SetPointError(8,0.0025,0.0025,0.004694736,0.004759421);
   grae->SetPoint(9,0.042,1);
   grae->SetPointError(9,0.003,0.003,0.00438431,0.004330876);
   grae->SetPoint(10,0.0485,1);
   grae->SetPointError(10,0.0035,0.0035,0.003982059,0.004227181);
   grae->SetPoint(11,0.0545,1);
   grae->SetPointError(11,0.0025,0.0025,0.005326253,0.005163359);
   grae->SetPoint(12,0.0605,1);
   grae->SetPointError(12,0.0035,0.0035,0.004393695,0.004475009);
   grae->SetPoint(13,0.068,1);
   grae->SetPointError(13,0.004,0.004,0.004455648,0.004482539);
   grae->SetPoint(14,0.0765,1);
   grae->SetPointError(14,0.0045,0.0045,0.004190541,0.004267437);
   grae->SetPoint(15,0.086,1);
   grae->SetPointError(15,0.005,0.005,0.004260442,0.004205803);
   grae->SetPoint(16,0.0965,1);
   grae->SetPointError(16,0.0055,0.0055,0.004285788,0.004312357);
   grae->SetPoint(17,0.108,1);
   grae->SetPointError(17,0.006,0.006,0.004305399,0.004401398);
   grae->SetPoint(18,0.121,1);
   grae->SetPointError(18,0.007,0.007,0.00429067,0.004263173);
   grae->SetPoint(19,0.1365,1);
   grae->SetPointError(19,0.0085,0.0085,0.004226132,0.004271211);
   grae->SetPoint(20,0.155,1);
   grae->SetPointError(20,0.01,0.01,0.004189642,0.004189911);
   grae->SetPoint(21,0.177,1);
   grae->SetPointError(21,0.012,0.012,0.004202421,0.004227179);
   grae->SetPoint(22,0.204,1);
   grae->SetPointError(22,0.015,0.015,0.004223922,0.004259339);
   grae->SetPoint(23,0.2385,1);
   grae->SetPointError(23,0.0195,0.0195,0.004200511,0.004229068);
   grae->SetPoint(24,0.285,1);
   grae->SetPointError(24,0.027,0.027,0.004215942,0.004230344);
   grae->SetPoint(25,0.3515,1);
   grae->SetPointError(25,0.0395,0.0395,0.004337619,0.004322881);
   grae->SetPoint(26,0.4575,1);
   grae->SetPointError(26,0.0665,0.0665,0.004584233,0.004529651);
   grae->SetPoint(27,0.6095,1);
   grae->SetPointError(27,0.0855,0.0855,0.005746868,0.005619484);
   grae->SetPoint(28,0.8065,1);
   grae->SetPointError(28,0.1115,0.1115,0.007527976,0.00734401);
   grae->SetPoint(29,1.0355,1);
   grae->SetPointError(29,0.1175,0.1175,0.009430274,0.00953222);
   grae->SetPoint(30,1.3245,1);
   grae->SetPointError(30,0.1715,0.1715,0.01119774,0.01065549);
   grae->SetPoint(31,1.7215,1);
   grae->SetPointError(31,0.2255,0.2255,0.01291288,0.01273643);
   grae->SetPoint(32,2.2345,1);
   grae->SetPointError(32,0.2875,0.2875,0.01742002,0.01734582);
   grae->SetPoint(33,2.8995,1);
   grae->SetPointError(33,0.3775,0.3775,0.02023863,0.02020248);
   
   TH1F *Graph_Graph_Graph2328 = new TH1F("Graph_Graph_Graph2328","",100,0.0006,3.604633);
   Graph_Graph_Graph2328->SetMinimum(0.7);
   Graph_Graph_Graph2328->SetMaximum(1.3);
   Graph_Graph_Graph2328->SetDirectory(0);
   Graph_Graph_Graph2328->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph2328->SetLineColor(ci);
   Graph_Graph_Graph2328->GetXaxis()->SetTitle("#phi*");
   Graph_Graph_Graph2328->GetXaxis()->SetRange(1,84);
   Graph_Graph_Graph2328->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph2328->GetXaxis()->SetLabelOffset(-0.01);
   Graph_Graph_Graph2328->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph_Graph2328->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph2328->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph_Graph2328->GetYaxis()->SetNdivisions(503);
   Graph_Graph_Graph2328->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph2328->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph2328->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph2328->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph2328->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2328->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2328->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph2328);
   
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
   grae->SetMarkerStyle(4);
   grae->SetPoint(0,0.002,0.9930977);
   grae->SetPointError(0,0,0,0.004420763,0.004298288);
   grae->SetPoint(1,0.006,0.974892);
   grae->SetPointError(1,0,0,0.004801765,0.00496134);
   grae->SetPoint(2,0.01,0.9898041);
   grae->SetPointError(2,0,0,0.004808034,0.004665582);
   grae->SetPoint(3,0.014,0.9854069);
   grae->SetPointError(3,0,0,0.004803764,0.004992469);
   grae->SetPoint(4,0.018,1.008544);
   grae->SetPointError(4,0,0,0.004993568,0.005127593);
   grae->SetPoint(5,0.022,1.00494);
   grae->SetPointError(5,0,0,0.004937079,0.005079358);
   grae->SetPoint(6,0.0265,1.014064);
   grae->SetPointError(6,0,0,0.004505659,0.00437283);
   grae->SetPoint(7,0.0315,1.024785);
   grae->SetPointError(7,0,0,0.004493937,0.004597712);
   grae->SetPoint(8,0.0365,1.030435);
   grae->SetPointError(8,0,0,0.004694736,0.004759421);
   grae->SetPoint(9,0.042,1.032399);
   grae->SetPointError(9,0,0,0.00438431,0.004330876);
   grae->SetPoint(10,0.0485,1.036759);
   grae->SetPointError(10,0,0,0.003982059,0.004227181);
   grae->SetPoint(11,0.0545,1.024036);
   grae->SetPointError(11,0,0,0.005326253,0.005163359);
   grae->SetPoint(12,0.0605,1.027067);
   grae->SetPointError(12,0,0,0.004393695,0.004475009);
   grae->SetPoint(13,0.068,1.029675);
   grae->SetPointError(13,0,0,0.004455648,0.004482539);
   grae->SetPoint(14,0.0765,1.015987);
   grae->SetPointError(14,0,0,0.004190541,0.004267437);
   grae->SetPoint(15,0.086,1.017109);
   grae->SetPointError(15,0,0,0.004260442,0.004205803);
   grae->SetPoint(16,0.0965,1.008504);
   grae->SetPointError(16,0,0,0.004285788,0.004312357);
   grae->SetPoint(17,0.108,1.002764);
   grae->SetPointError(17,0,0,0.004305399,0.004401398);
   grae->SetPoint(18,0.121,0.9966581);
   grae->SetPointError(18,0,0,0.00429067,0.004263173);
   grae->SetPoint(19,0.1365,1.000027);
   grae->SetPointError(19,0,0,0.004226132,0.004271211);
   grae->SetPoint(20,0.155,0.990823);
   grae->SetPointError(20,0,0,0.004189642,0.004189911);
   grae->SetPoint(21,0.177,0.9821747);
   grae->SetPointError(21,0,0,0.004202421,0.004227179);
   grae->SetPoint(22,0.204,0.9885649);
   grae->SetPointError(22,0,0,0.004223922,0.004259339);
   grae->SetPoint(23,0.2385,0.9749423);
   grae->SetPointError(23,0,0,0.004200511,0.004229068);
   grae->SetPoint(24,0.285,0.9677);
   grae->SetPointError(24,0,0,0.004215942,0.004230344);
   grae->SetPoint(25,0.3515,0.9679816);
   grae->SetPointError(25,0,0,0.004337619,0.004322881);
   grae->SetPoint(26,0.4575,0.9555626);
   grae->SetPointError(26,0,0,0.004584233,0.004529651);
   grae->SetPoint(27,0.6095,0.9645466);
   grae->SetPointError(27,0,0,0.005746868,0.005619484);
   grae->SetPoint(28,0.8065,0.9531825);
   grae->SetPointError(28,0,0,0.007527976,0.00734401);
   grae->SetPoint(29,1.0355,0.970836);
   grae->SetPointError(29,0,0,0.009430274,0.00953222);
   grae->SetPoint(30,1.3245,0.971061);
   grae->SetPointError(30,0,0,0.01119774,0.01065549);
   grae->SetPoint(31,1.7215,1.028746);
   grae->SetPointError(31,0,0,0.01291288,0.01273643);
   grae->SetPoint(32,2.2345,1.012352);
   grae->SetPointError(32,0,0,0.01742002,0.01734582);
   grae->SetPoint(33,2.8995,0.9973977);
   grae->SetPointError(33,0,0,0.02023863,0.02020248);
   
   TH1F *Graph_Graph_Graph2429 = new TH1F("Graph_Graph_Graph2429","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph2429->SetMinimum(0.9360717);
   Graph_Graph_Graph2429->SetMaximum(1.051065);
   Graph_Graph_Graph2429->SetDirectory(0);
   Graph_Graph_Graph2429->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph2429->SetLineColor(ci);
   Graph_Graph_Graph2429->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph2429->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2429->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2429->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph2429->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph2429->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2429->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2429->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph2429->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph2429->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2429->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2429->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph2429);
   
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
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.002,0.9074065);
   grae->SetPointError(0,0,0,0.004420763,0.004298288);
   grae->SetPoint(1,0.006,0.9086841);
   grae->SetPointError(1,0,0,0.004801765,0.00496134);
   grae->SetPoint(2,0.01,0.9150614);
   grae->SetPointError(2,0,0,0.004808034,0.004665582);
   grae->SetPoint(3,0.014,0.9121222);
   grae->SetPointError(3,0,0,0.004803764,0.004992469);
   grae->SetPoint(4,0.018,0.9333089);
   grae->SetPointError(4,0,0,0.004993568,0.005127593);
   grae->SetPoint(5,0.022,0.9612957);
   grae->SetPointError(5,0,0,0.004937079,0.005079358);
   grae->SetPoint(6,0.0265,0.9756665);
   grae->SetPointError(6,0,0,0.004505659,0.00437283);
   grae->SetPoint(7,0.0315,1.004007);
   grae->SetPointError(7,0,0,0.004493937,0.004597712);
   grae->SetPoint(8,0.0365,1.022146);
   grae->SetPointError(8,0,0,0.004694736,0.004759421);
   grae->SetPoint(9,0.042,1.041121);
   grae->SetPointError(9,0,0,0.00438431,0.004330876);
   grae->SetPoint(10,0.0485,1.064186);
   grae->SetPointError(10,0,0,0.003982059,0.004227181);
   grae->SetPoint(11,0.0545,1.06343);
   grae->SetPointError(11,0,0,0.005326253,0.005163359);
   grae->SetPoint(12,0.0605,1.069323);
   grae->SetPointError(12,0,0,0.004393695,0.004475009);
   grae->SetPoint(13,0.068,1.083297);
   grae->SetPointError(13,0,0,0.004455648,0.004482539);
   grae->SetPoint(14,0.0765,1.088108);
   grae->SetPointError(14,0,0,0.004190541,0.004267437);
   grae->SetPoint(15,0.086,1.088475);
   grae->SetPointError(15,0,0,0.004260442,0.004205803);
   grae->SetPoint(16,0.0965,1.068047);
   grae->SetPointError(16,0,0,0.004285788,0.004312357);
   grae->SetPoint(17,0.108,1.062267);
   grae->SetPointError(17,0,0,0.004305399,0.004401398);
   grae->SetPoint(18,0.121,1.047408);
   grae->SetPointError(18,0,0,0.00429067,0.004263173);
   grae->SetPoint(19,0.1365,1.031817);
   grae->SetPointError(19,0,0,0.004226132,0.004271211);
   grae->SetPoint(20,0.155,1.004501);
   grae->SetPointError(20,0,0,0.004189642,0.004189911);
   grae->SetPoint(21,0.177,1.00426);
   grae->SetPointError(21,0,0,0.004202421,0.004227179);
   grae->SetPoint(22,0.204,1.005625);
   grae->SetPointError(22,0,0,0.004223922,0.004259339);
   grae->SetPoint(23,0.2385,0.999577);
   grae->SetPointError(23,0,0,0.004200511,0.004229068);
   grae->SetPoint(24,0.285,0.9881611);
   grae->SetPointError(24,0,0,0.004215942,0.004230344);
   grae->SetPoint(25,0.3515,0.9746713);
   grae->SetPointError(25,0,0,0.004337619,0.004322881);
   grae->SetPoint(26,0.4575,0.9586801);
   grae->SetPointError(26,0,0,0.004584233,0.004529651);
   grae->SetPoint(27,0.6095,0.9287542);
   grae->SetPointError(27,0,0,0.005746868,0.005619484);
   grae->SetPoint(28,0.8065,0.9024581);
   grae->SetPointError(28,0,0,0.007527976,0.00734401);
   grae->SetPoint(29,1.0355,0.8826407);
   grae->SetPointError(29,0,0,0.009430274,0.00953222);
   grae->SetPoint(30,1.3245,0.8413076);
   grae->SetPointError(30,0,0,0.01119774,0.01065549);
   grae->SetPoint(31,1.7215,0.8340805);
   grae->SetPointError(31,0,0,0.01291288,0.01273643);
   grae->SetPoint(32,2.2345,0.8200491);
   grae->SetPointError(32,0,0,0.01742002,0.01734582);
   grae->SetPoint(33,2.8995,0.8575689);
   grae->SetPointError(33,0,0,0.02023863,0.02020248);
   
   TH1F *Graph_Graph_Graph2530 = new TH1F("Graph_Graph_Graph2530","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph2530->SetMinimum(0.7736239);
   Graph_Graph_Graph2530->SetMaximum(1.121686);
   Graph_Graph_Graph2530->SetDirectory(0);
   Graph_Graph_Graph2530->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph2530->SetLineColor(ci);
   Graph_Graph_Graph2530->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph2530->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2530->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2530->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph2530->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph2530->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2530->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2530->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph2530->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph2530->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2530->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2530->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph2530);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph3");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#009900");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.002,1.061826);
   grae->SetPointError(0,0,0,0.004420763,0.004298288);
   grae->SetPoint(1,0.006,1.060938);
   grae->SetPointError(1,0,0,0.004801765,0.00496134);
   grae->SetPoint(2,0.01,1.075712);
   grae->SetPointError(2,0,0,0.004808034,0.004665582);
   grae->SetPoint(3,0.014,1.035353);
   grae->SetPointError(3,0,0,0.004803764,0.004992469);
   grae->SetPoint(4,0.018,1.04555);
   grae->SetPointError(4,0,0,0.004993568,0.005127593);
   grae->SetPoint(5,0.022,1.069066);
   grae->SetPointError(5,0,0,0.004937079,0.005079358);
   grae->SetPoint(6,0.0265,1.059003);
   grae->SetPointError(6,0,0,0.004505659,0.00437283);
   grae->SetPoint(7,0.0315,1.07124);
   grae->SetPointError(7,0,0,0.004493937,0.004597712);
   grae->SetPoint(8,0.0365,1.044673);
   grae->SetPointError(8,0,0,0.004694736,0.004759421);
   grae->SetPoint(9,0.042,1.042866);
   grae->SetPointError(9,0,0,0.00438431,0.004330876);
   grae->SetPoint(10,0.0485,1.026906);
   grae->SetPointError(10,0,0,0.003982059,0.004227181);
   grae->SetPoint(11,0.0545,1.023805);
   grae->SetPointError(11,0,0,0.005326253,0.005163359);
   grae->SetPoint(12,0.0605,1.0111);
   grae->SetPointError(12,0,0,0.004393695,0.004475009);
   grae->SetPoint(13,0.068,1.013074);
   grae->SetPointError(13,0,0,0.004455648,0.004482539);
   grae->SetPoint(14,0.0765,0.9955082);
   grae->SetPointError(14,0,0,0.004190541,0.004267437);
   grae->SetPoint(15,0.086,0.9994183);
   grae->SetPointError(15,0,0,0.004260442,0.004205803);
   grae->SetPoint(16,0.0965,0.9927216);
   grae->SetPointError(16,0,0,0.004285788,0.004312357);
   grae->SetPoint(17,0.108,0.972626);
   grae->SetPointError(17,0,0,0.004305399,0.004401398);
   grae->SetPoint(18,0.121,0.9700977);
   grae->SetPointError(18,0,0,0.00429067,0.004263173);
   grae->SetPoint(19,0.1365,0.9582074);
   grae->SetPointError(19,0,0,0.004226132,0.004271211);
   grae->SetPoint(20,0.155,0.9674132);
   grae->SetPointError(20,0,0,0.004189642,0.004189911);
   grae->SetPoint(21,0.177,0.9646633);
   grae->SetPointError(21,0,0,0.004202421,0.004227179);
   grae->SetPoint(22,0.204,0.9804985);
   grae->SetPointError(22,0,0,0.004223922,0.004259339);
   grae->SetPoint(23,0.2385,0.9678674);
   grae->SetPointError(23,0,0,0.004200511,0.004229068);
   grae->SetPoint(24,0.285,0.9314136);
   grae->SetPointError(24,0,0,0.004215942,0.004230344);
   grae->SetPoint(25,0.3515,0.9392079);
   grae->SetPointError(25,0,0,0.004337619,0.004322881);
   grae->SetPoint(26,0.4575,0.9070061);
   grae->SetPointError(26,0,0,0.004584233,0.004529651);
   grae->SetPoint(27,0.6095,0.8983865);
   grae->SetPointError(27,0,0,0.005746868,0.005619484);
   grae->SetPoint(28,0.8065,0.8730775);
   grae->SetPointError(28,0,0,0.007527976,0.00734401);
   grae->SetPoint(29,1.0355,0.8376494);
   grae->SetPointError(29,0,0,0.009430274,0.00953222);
   grae->SetPoint(30,1.3245,0.8645467);
   grae->SetPointError(30,0,0,0.01119774,0.01065549);
   grae->SetPoint(31,1.7215,0.8930329);
   grae->SetPointError(31,0,0,0.01291288,0.01273643);
   grae->SetPoint(32,2.2345,0.8983404);
   grae->SetPointError(32,0,0,0.01742002,0.01734582);
   grae->SetPoint(33,2.8995,0.9200284);
   grae->SetPointError(33,0,0,0.02023863,0.02020248);
   
   TH1F *Graph_Graph_Graph2631 = new TH1F("Graph_Graph_Graph2631","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph2631->SetMinimum(0.8030033);
   Graph_Graph_Graph2631->SetMaximum(1.105594);
   Graph_Graph_Graph2631->SetDirectory(0);
   Graph_Graph_Graph2631->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph2631->SetLineColor(ci);
   Graph_Graph_Graph2631->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph2631->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2631->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2631->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph2631->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph2631->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2631->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2631->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph2631->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph2631->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2631->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2631->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph2631);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph4");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff6600");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff6600");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetPoint(0,0.002,1.056243);
   grae->SetPointError(0,0,0,0.004420763,0.004298288);
   grae->SetPoint(1,0.006,1.055603);
   grae->SetPointError(1,0,0,0.004801765,0.00496134);
   grae->SetPoint(2,0.01,1.048985);
   grae->SetPointError(2,0,0,0.004808034,0.004665582);
   grae->SetPoint(3,0.014,1.045766);
   grae->SetPointError(3,0,0,0.004803764,0.004992469);
   grae->SetPoint(4,0.018,1.056944);
   grae->SetPointError(4,0,0,0.004993568,0.005127593);
   grae->SetPoint(5,0.022,1.064272);
   grae->SetPointError(5,0,0,0.004937079,0.005079358);
   grae->SetPoint(6,0.0265,1.043833);
   grae->SetPointError(6,0,0,0.004505659,0.00437283);
   grae->SetPoint(7,0.0315,1.074286);
   grae->SetPointError(7,0,0,0.004493937,0.004597712);
   grae->SetPoint(8,0.0365,1.075934);
   grae->SetPointError(8,0,0,0.004694736,0.004759421);
   grae->SetPoint(9,0.042,1.046602);
   grae->SetPointError(9,0,0,0.00438431,0.004330876);
   grae->SetPoint(10,0.0485,1.03413);
   grae->SetPointError(10,0,0,0.003982059,0.004227181);
   grae->SetPoint(11,0.0545,1.020969);
   grae->SetPointError(11,0,0,0.005326253,0.005163359);
   grae->SetPoint(12,0.0605,1.023193);
   grae->SetPointError(12,0,0,0.004393695,0.004475009);
   grae->SetPoint(13,0.068,1.027772);
   grae->SetPointError(13,0,0,0.004455648,0.004482539);
   grae->SetPoint(14,0.0765,0.9972674);
   grae->SetPointError(14,0,0,0.004190541,0.004267437);
   grae->SetPoint(15,0.086,0.9856335);
   grae->SetPointError(15,0,0,0.004260442,0.004205803);
   grae->SetPoint(16,0.0965,0.9919204);
   grae->SetPointError(16,0,0,0.004285788,0.004312357);
   grae->SetPoint(17,0.108,0.9710293);
   grae->SetPointError(17,0,0,0.004305399,0.004401398);
   grae->SetPoint(18,0.121,0.9558788);
   grae->SetPointError(18,0,0,0.00429067,0.004263173);
   grae->SetPoint(19,0.1365,0.9605227);
   grae->SetPointError(19,0,0,0.004226132,0.004271211);
   grae->SetPoint(20,0.155,0.9589894);
   grae->SetPointError(20,0,0,0.004189642,0.004189911);
   grae->SetPoint(21,0.177,0.9564709);
   grae->SetPointError(21,0,0,0.004202421,0.004227179);
   grae->SetPoint(22,0.204,0.9630035);
   grae->SetPointError(22,0,0,0.004223922,0.004259339);
   grae->SetPoint(23,0.2385,0.9582963);
   grae->SetPointError(23,0,0,0.004200511,0.004229068);
   grae->SetPoint(24,0.285,0.9513278);
   grae->SetPointError(24,0,0,0.004215942,0.004230344);
   grae->SetPoint(25,0.3515,0.9457955);
   grae->SetPointError(25,0,0,0.004337619,0.004322881);
   grae->SetPoint(26,0.4575,0.8997428);
   grae->SetPointError(26,0,0,0.004584233,0.004529651);
   grae->SetPoint(27,0.6095,0.896371);
   grae->SetPointError(27,0,0,0.005746868,0.005619484);
   grae->SetPoint(28,0.8065,0.8656377);
   grae->SetPointError(28,0,0,0.007527976,0.00734401);
   grae->SetPoint(29,1.0355,0.8767441);
   grae->SetPointError(29,0,0,0.009430274,0.00953222);
   grae->SetPoint(30,1.3245,0.8450043);
   grae->SetPointError(30,0,0,0.01119774,0.01065549);
   grae->SetPoint(31,1.7215,0.9486875);
   grae->SetPointError(31,0,0,0.01291288,0.01273643);
   grae->SetPoint(32,2.2345,0.9949353);
   grae->SetPointError(32,0,0,0.01742002,0.01734582);
   grae->SetPoint(33,2.8995,1.010197);
   grae->SetPointError(33,0,0,0.02023863,0.02020248);
   
   TH1F *Graph_Graph_Graph2732 = new TH1F("Graph_Graph_Graph2732","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph2732->SetMinimum(0.8091179);
   Graph_Graph_Graph2732->SetMaximum(1.105382);
   Graph_Graph_Graph2732->SetDirectory(0);
   Graph_Graph_Graph2732->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph2732->SetLineColor(ci);
   Graph_Graph_Graph2732->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph2732->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2732->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2732->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph2732->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph2732->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2732->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2732->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph2732->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph2732->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph2732->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph2732->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph2732);
   
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.7,0.907,"19.7 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.87,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.25,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.2,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.15,"60 GeV < M_{ee} < 120 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.13,0.62,0.95,0.84,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","2012 data","F");

   ci = TColor::GetColor("#999999");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Z #rightarrow ee MadGraph+Pythia6 (Z2star)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(4);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Z #rightarrow ee POWHEG+Pythia6 (Z2star)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","Z #rightarrow ee POWHEG+Pythia8  (Tunepp 5)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","Z #rightarrow ee POWHEG+Pythia8 (Tunepp 14)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff6600");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   FinalPhiRatio->Modified();
   FinalPhiRatio->cd();
   FinalPhiRatio->SetSelected(FinalPhiRatio);
}
