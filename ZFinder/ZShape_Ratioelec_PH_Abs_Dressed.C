{
//=========Macro generated from canvas: FinalPhiRatio/FinalPhiRatio
//=========  (Tue May  5 15:36:53 2015) by ROOT version5.34/28
   TCanvas *FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio",6,97,800,900);
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

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetPoint(0,0.002,1);
   grae->SetPointError(0,0.002,0.002,0.02809126,0.02790952);
   grae->SetPoint(1,0.006,1);
   grae->SetPointError(1,0.002,0.002,0.02851094,0.02845184);
   grae->SetPoint(2,0.01,1);
   grae->SetPointError(2,0.002,0.002,0.02868182,0.02882308);
   grae->SetPoint(3,0.014,1);
   grae->SetPointError(3,0.002,0.002,0.02854511,0.0288731);
   grae->SetPoint(4,0.018,1);
   grae->SetPointError(4,0.002,0.002,0.02851053,0.02816792);
   grae->SetPoint(5,0.022,1);
   grae->SetPointError(5,0.002,0.002,0.02870409,0.02870386);
   grae->SetPoint(6,0.0265,1);
   grae->SetPointError(6,0.0025,0.0025,0.02791722,0.02792448);
   grae->SetPoint(7,0.0315,1);
   grae->SetPointError(7,0.0025,0.0025,0.02792113,0.02787675);
   grae->SetPoint(8,0.0365,1);
   grae->SetPointError(8,0.0025,0.0025,0.02836554,0.02828455);
   grae->SetPoint(9,0.042,1);
   grae->SetPointError(9,0.003,0.003,0.02791216,0.02785795);
   grae->SetPoint(10,0.0485,1);
   grae->SetPointError(10,0.0035,0.0035,0.02823525,0.02813367);
   grae->SetPoint(11,0.0545,1);
   grae->SetPointError(11,0.0025,0.0025,0.02827969,0.02848347);
   grae->SetPoint(12,0.0605,1);
   grae->SetPointError(12,0.0035,0.0035,0.0279947,0.02810981);
   grae->SetPoint(13,0.068,1);
   grae->SetPointError(13,0.004,0.004,0.0283045,0.02809505);
   grae->SetPoint(14,0.0765,1);
   grae->SetPointError(14,0.0045,0.0045,0.02776272,0.02787052);
   grae->SetPoint(15,0.086,1);
   grae->SetPointError(15,0.005,0.005,0.02797781,0.02791598);
   grae->SetPoint(16,0.0965,1);
   grae->SetPointError(16,0.0055,0.0055,0.02786719,0.02779345);
   grae->SetPoint(17,0.108,1);
   grae->SetPointError(17,0.006,0.006,0.02804569,0.02802867);
   grae->SetPoint(18,0.121,1);
   grae->SetPointError(18,0.007,0.007,0.02787522,0.02787028);
   grae->SetPoint(19,0.1365,1);
   grae->SetPointError(19,0.0085,0.0085,0.02796356,0.02799082);
   grae->SetPoint(20,0.155,1);
   grae->SetPointError(20,0.01,0.01,0.02761253,0.02769821);
   grae->SetPoint(21,0.177,1);
   grae->SetPointError(21,0.012,0.012,0.02773791,0.0277793);
   grae->SetPoint(22,0.204,1);
   grae->SetPointError(22,0.015,0.015,0.02788096,0.02784412);
   grae->SetPoint(23,0.2385,1);
   grae->SetPointError(23,0.0195,0.0195,0.0277587,0.02774809);
   grae->SetPoint(24,0.285,1);
   grae->SetPointError(24,0.027,0.027,0.02788039,0.02783703);
   grae->SetPoint(25,0.3515,1);
   grae->SetPointError(25,0.0395,0.0395,0.02783791,0.02777281);
   grae->SetPoint(26,0.4575,1);
   grae->SetPointError(26,0.0665,0.0665,0.02784636,0.02797055);
   grae->SetPoint(27,0.6095,1);
   grae->SetPointError(27,0.0855,0.0855,0.02857616,0.02809262);
   grae->SetPoint(28,0.8065,1);
   grae->SetPointError(28,0.1115,0.1115,0.02904505,0.02894753);
   grae->SetPoint(29,1.0355,1);
   grae->SetPointError(29,0.1175,0.1175,0.03112653,0.0310173);
   grae->SetPoint(30,1.3245,1);
   grae->SetPointError(30,0.1715,0.1715,0.03314886,0.03275989);
   grae->SetPoint(31,1.7215,1);
   grae->SetPointError(31,0.2255,0.2255,0.03542199,0.03578039);
   grae->SetPoint(32,2.2345,1);
   grae->SetPointError(32,0.2875,0.2875,0.04021734,0.0400026);
   grae->SetPoint(33,2.8995,1);
   grae->SetPointError(33,0.3775,0.3775,0.04467076,0.04399739);
   
   TH1F *Graph_Graph_Graph1518 = new TH1F("Graph_Graph_Graph1518","",100,0.0006,3.604633);
   Graph_Graph_Graph1518->SetMinimum(0.7);
   Graph_Graph_Graph1518->SetMaximum(1.3);
   Graph_Graph_Graph1518->SetDirectory(0);
   Graph_Graph_Graph1518->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph1518->SetLineColor(ci);
   Graph_Graph_Graph1518->GetXaxis()->SetTitle("#phi*");
   Graph_Graph_Graph1518->GetXaxis()->SetRange(1,84);
   Graph_Graph_Graph1518->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph1518->GetXaxis()->SetLabelOffset(-0.01);
   Graph_Graph_Graph1518->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph_Graph1518->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph1518->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph_Graph1518->GetYaxis()->SetNdivisions(503);
   Graph_Graph_Graph1518->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph1518->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph1518->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph1518->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph1518->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1518->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1518->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph1518);
   
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
   grae->SetPoint(0,0.002,0.9484236);
   grae->SetPointError(0,0,0,0.02809126,0.02790952);
   grae->SetPoint(1,0.006,0.94912);
   grae->SetPointError(1,0,0,0.02851094,0.02845184);
   grae->SetPoint(2,0.01,0.9632997);
   grae->SetPointError(2,0,0,0.02868182,0.02882308);
   grae->SetPoint(3,0.014,0.9690533);
   grae->SetPointError(3,0,0,0.02854511,0.0288731);
   grae->SetPoint(4,0.018,0.9659571);
   grae->SetPointError(4,0,0,0.02851053,0.02816792);
   grae->SetPoint(5,0.022,0.9953682);
   grae->SetPointError(5,0,0,0.02870409,0.02870386);
   grae->SetPoint(6,0.0265,0.9812144);
   grae->SetPointError(6,0,0,0.02791722,0.02792448);
   grae->SetPoint(7,0.0315,0.9953846);
   grae->SetPointError(7,0,0,0.02792113,0.02787675);
   grae->SetPoint(8,0.0365,0.9933999);
   grae->SetPointError(8,0,0,0.02836554,0.02828455);
   grae->SetPoint(9,0.042,1.000372);
   grae->SetPointError(9,0,0,0.02791216,0.02785795);
   grae->SetPoint(10,0.0485,1.005812);
   grae->SetPointError(10,0,0,0.02823525,0.02813367);
   grae->SetPoint(11,0.0545,1.001951);
   grae->SetPointError(11,0,0,0.02827969,0.02848347);
   grae->SetPoint(12,0.0605,1.017649);
   grae->SetPointError(12,0,0,0.0279947,0.02810981);
   grae->SetPoint(13,0.068,0.9988122);
   grae->SetPointError(13,0,0,0.0283045,0.02809505);
   grae->SetPoint(14,0.0765,0.9959211);
   grae->SetPointError(14,0,0,0.02776272,0.02787052);
   grae->SetPoint(15,0.086,1.001263);
   grae->SetPointError(15,0,0,0.02797781,0.02791598);
   grae->SetPoint(16,0.0965,0.9888672);
   grae->SetPointError(16,0,0,0.02786719,0.02779345);
   grae->SetPoint(17,0.108,0.9711499);
   grae->SetPointError(17,0,0,0.02804569,0.02802867);
   grae->SetPoint(18,0.121,0.9851277);
   grae->SetPointError(18,0,0,0.02787522,0.02787028);
   grae->SetPoint(19,0.1365,0.9726936);
   grae->SetPointError(19,0,0,0.02796356,0.02799082);
   grae->SetPoint(20,0.155,0.9673759);
   grae->SetPointError(20,0,0,0.02761253,0.02769821);
   grae->SetPoint(21,0.177,0.9654481);
   grae->SetPointError(21,0,0,0.02773791,0.0277793);
   grae->SetPoint(22,0.204,0.9604686);
   grae->SetPointError(22,0,0,0.02788096,0.02784412);
   grae->SetPoint(23,0.2385,0.9563618);
   grae->SetPointError(23,0,0,0.0277587,0.02774809);
   grae->SetPoint(24,0.285,0.961544);
   grae->SetPointError(24,0,0,0.02788039,0.02783703);
   grae->SetPoint(25,0.3515,0.9516155);
   grae->SetPointError(25,0,0,0.02783791,0.02777281);
   grae->SetPoint(26,0.4575,0.9592932);
   grae->SetPointError(26,0,0,0.02784636,0.02797055);
   grae->SetPoint(27,0.6095,0.9524333);
   grae->SetPointError(27,0,0,0.02857616,0.02809262);
   grae->SetPoint(28,0.8065,0.9366437);
   grae->SetPointError(28,0,0,0.02904505,0.02894753);
   grae->SetPoint(29,1.0355,0.9314586);
   grae->SetPointError(29,0,0,0.03112653,0.0310173);
   grae->SetPoint(30,1.3245,0.9590798);
   grae->SetPointError(30,0,0,0.03314886,0.03275989);
   grae->SetPoint(31,1.7215,1.033778);
   grae->SetPointError(31,0,0,0.03542199,0.03578039);
   grae->SetPoint(32,2.2345,0.9923798);
   grae->SetPointError(32,0,0,0.04021734,0.0400026);
   grae->SetPoint(33,2.8995,1.005431);
   grae->SetPointError(33,0,0,0.04467076,0.04399739);
   
   TH1F *Graph_Graph_Graph1619 = new TH1F("Graph_Graph_Graph1619","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph1619->SetMinimum(0.8834094);
   Graph_Graph_Graph1619->SetMaximum(1.086481);
   Graph_Graph_Graph1619->SetDirectory(0);
   Graph_Graph_Graph1619->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph1619->SetLineColor(ci);
   Graph_Graph_Graph1619->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph1619->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1619->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1619->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph1619->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph1619->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1619->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1619->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph1619->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph1619->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1619->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1619->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph1619);
   
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
   grae->SetPoint(0,0.002,0.8836689);
   grae->SetPointError(0,0,0,0.02809126,0.02790952);
   grae->SetPoint(1,0.006,0.9021003);
   grae->SetPointError(1,0,0,0.02851094,0.02845184);
   grae->SetPoint(2,0.01,0.9081127);
   grae->SetPointError(2,0,0,0.02868182,0.02882308);
   grae->SetPoint(3,0.014,0.9146657);
   grae->SetPointError(3,0,0,0.02854511,0.0288731);
   grae->SetPoint(4,0.018,0.9115193);
   grae->SetPointError(4,0,0,0.02851053,0.02816792);
   grae->SetPoint(5,0.022,0.9709077);
   grae->SetPointError(5,0,0,0.02870409,0.02870386);
   grae->SetPoint(6,0.0265,0.9626698);
   grae->SetPointError(6,0,0,0.02791722,0.02792448);
   grae->SetPoint(7,0.0315,0.9944255);
   grae->SetPointError(7,0,0,0.02792113,0.02787675);
   grae->SetPoint(8,0.0365,1.004833);
   grae->SetPointError(8,0,0,0.02836554,0.02828455);
   grae->SetPoint(9,0.042,1.028709);
   grae->SetPointError(9,0,0,0.02791216,0.02785795);
   grae->SetPoint(10,0.0485,1.052771);
   grae->SetPointError(10,0,0,0.02823525,0.02813367);
   grae->SetPoint(11,0.0545,1.061005);
   grae->SetPointError(11,0,0,0.02827969,0.02848347);
   grae->SetPoint(12,0.0605,1.080402);
   grae->SetPointError(12,0,0,0.0279947,0.02810981);
   grae->SetPoint(13,0.068,1.07154);
   grae->SetPointError(13,0,0,0.0283045,0.02809505);
   grae->SetPoint(14,0.0765,1.087643);
   grae->SetPointError(14,0,0,0.02776272,0.02787052);
   grae->SetPoint(15,0.086,1.092638);
   grae->SetPointError(15,0,0,0.02797781,0.02791598);
   grae->SetPoint(16,0.0965,1.067894);
   grae->SetPointError(16,0,0,0.02786719,0.02779345);
   grae->SetPoint(17,0.108,1.049056);
   grae->SetPointError(17,0,0,0.02804569,0.02802867);
   grae->SetPoint(18,0.121,1.055698);
   grae->SetPointError(18,0,0,0.02787522,0.02787028);
   grae->SetPoint(19,0.1365,1.023397);
   grae->SetPointError(19,0,0,0.02796356,0.02799082);
   grae->SetPoint(20,0.155,1.000062);
   grae->SetPointError(20,0,0,0.02761253,0.02769821);
   grae->SetPoint(21,0.177,1.006616);
   grae->SetPointError(21,0,0,0.02773791,0.0277793);
   grae->SetPoint(22,0.204,0.9963032);
   grae->SetPointError(22,0,0,0.02788096,0.02784412);
   grae->SetPoint(23,0.2385,0.9998545);
   grae->SetPointError(23,0,0,0.0277587,0.02774809);
   grae->SetPoint(24,0.285,1.001229);
   grae->SetPointError(24,0,0,0.02788039,0.02783703);
   grae->SetPoint(25,0.3515,0.9770795);
   grae->SetPointError(25,0,0,0.02783791,0.02777281);
   grae->SetPoint(26,0.4575,0.9813936);
   grae->SetPointError(26,0,0,0.02784636,0.02797055);
   grae->SetPoint(27,0.6095,0.9351676);
   grae->SetPointError(27,0,0,0.02857616,0.02809262);
   grae->SetPoint(28,0.8065,0.9042795);
   grae->SetPointError(28,0,0,0.02904505,0.02894753);
   grae->SetPoint(29,1.0355,0.8635329);
   grae->SetPointError(29,0,0,0.03112653,0.0310173);
   grae->SetPoint(30,1.3245,0.8473062);
   grae->SetPointError(30,0,0,0.03314886,0.03275989);
   grae->SetPoint(31,1.7215,0.8546817);
   grae->SetPointError(31,0,0,0.03542199,0.03578039);
   grae->SetPoint(32,2.2345,0.8197163);
   grae->SetPointError(32,0,0,0.04021734,0.0400026);
   grae->SetPoint(33,2.8995,0.8815159);
   grae->SetPointError(33,0,0,0.04467076,0.04399739);
   
   TH1F *Graph_Graph_Graph1720 = new TH1F("Graph_Graph_Graph1720","Graph",100,0.0018,3.18925);
   Graph_Graph_Graph1720->SetMinimum(0.7453934);
   Graph_Graph_Graph1720->SetMaximum(1.15466);
   Graph_Graph_Graph1720->SetDirectory(0);
   Graph_Graph_Graph1720->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph1720->SetLineColor(ci);
   Graph_Graph_Graph1720->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph1720->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1720->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1720->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph1720->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph1720->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1720->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1720->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph1720->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph1720->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph1720->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph1720->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph1720);
   
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
   
   TLegend *leg = new TLegend(0.53,0.72,0.85,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
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
   FinalPhiRatio->Modified();
   FinalPhiRatio->cd();
   FinalPhiRatio->SetSelected(FinalPhiRatio);
}
