{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Tue Apr 14 13:38:56 2015) by ROOT version5.34/18
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",4,79,800,900);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-3.924951,-0.1077066,0.7623949,10.66295);
   p1->SetFillColor(0);
   p1->SetBorderMode(0);
   p1->SetBorderSize(0);
   p1->SetLogx();
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.06);
   p1->SetTopMargin(0.06);
   p1->SetBottomMargin(0.01);
   p1->SetFrameBorderMode(0);
   p1->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph0");
   grae->SetTitle("");
   grae->SetPoint(0,0.002,9.106297);
   grae->SetPointError(0,0.002,0.002,0,0);
   grae->SetPoint(1,0.006,9.078291);
   grae->SetPointError(1,0.002,0.002,0,0);
   grae->SetPoint(2,0.01,8.860786);
   grae->SetPointError(2,0.002,0.002,0,0);
   grae->SetPoint(3,0.014,8.676274);
   grae->SetPointError(3,0.002,0.002,0,0);
   grae->SetPoint(4,0.018,8.337609);
   grae->SetPointError(4,0.002,0.002,0,0);
   grae->SetPoint(5,0.022,8.044422);
   grae->SetPointError(5,0.002,0.002,0,0);
   grae->SetPoint(6,0.0265,7.659918);
   grae->SetPointError(6,0.0025,0.0025,0,0);
   grae->SetPoint(7,0.0315,7.17912);
   grae->SetPointError(7,0.0025,0.0025,0,0);
   grae->SetPoint(8,0.0365,6.728084);
   grae->SetPointError(8,0.0025,0.0025,0,0);
   grae->SetPoint(9,0.042,6.282668);
   grae->SetPointError(9,0.003,0.003,0,0);
   grae->SetPoint(10,0.0485,5.775932);
   grae->SetPointError(10,0.0035,0.0035,0,0);
   grae->SetPoint(11,0.0545,5.332089);
   grae->SetPointError(11,0.0025,0.0025,0,0);
   grae->SetPoint(12,0.0605,4.923704);
   grae->SetPointError(12,0.0035,0.0035,0,0);
   grae->SetPoint(13,0.068,4.449297);
   grae->SetPointError(13,0.004,0.004,0,0);
   grae->SetPoint(14,0.0765,4.011493);
   grae->SetPointError(14,0.0045,0.0045,0,0);
   grae->SetPoint(15,0.086,3.564145);
   grae->SetPointError(15,0.005,0.005,0,0);
   grae->SetPoint(16,0.0965,3.159969);
   grae->SetPointError(16,0.0055,0.0055,0,0);
   grae->SetPoint(17,0.108,2.77742);
   grae->SetPointError(17,0.006,0.006,0,0);
   grae->SetPoint(18,0.121,2.420708);
   grae->SetPointError(18,0.007,0.007,0,0);
   grae->SetPoint(19,0.1365,2.071115);
   grae->SetPointError(19,0.0085,0.0085,0,0);
   grae->SetPoint(20,0.155,1.734207);
   grae->SetPointError(20,0.01,0.01,0,0);
   grae->SetPoint(21,0.177,1.424907);
   grae->SetPointError(21,0.012,0.012,0,0);
   grae->SetPoint(22,0.204,1.133518);
   grae->SetPointError(22,0.015,0.015,0,0);
   grae->SetPoint(23,0.2385,0.876128);
   grae->SetPointError(23,0.0195,0.0195,0,0);
   grae->SetPoint(24,0.285,0.6394784);
   grae->SetPointError(24,0.027,0.027,0,0);
   grae->SetPoint(25,0.3515,0.4301066);
   grae->SetPointError(25,0.0395,0.0395,0,0);
   grae->SetPoint(26,0.4575,0.2519978);
   grae->SetPointError(26,0.0665,0.0665,0,0);
   grae->SetPoint(27,0.6095,0.1319807);
   grae->SetPointError(27,0.0855,0.0855,0,0);
   grae->SetPoint(28,0.8065,0.06668263);
   grae->SetPointError(28,0.1115,0.1115,0,0);
   grae->SetPoint(29,1.0355,0.0345134);
   grae->SetPointError(29,0.1175,0.1175,0,0);
   grae->SetPoint(30,1.3245,0.01800679);
   grae->SetPointError(30,0.1715,0.1715,0,0);
   grae->SetPoint(31,1.7215,0.00856106);
   grae->SetPointError(31,0.2255,0.2255,0,0);
   grae->SetPoint(32,2.2345,0.00417938);
   grae->SetPointError(32,0.2875,0.2875,0,0);
   grae->SetPoint(33,2.8995,0.002152797);
   grae->SetPointError(33,0.3775,0.3775,0,0);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","",100,0.0006,3.604633);
   Graph_Graph6->SetMinimum(0);
   Graph_Graph6->SetMaximum(10.01671);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph6->SetLineColor(ci);
   Graph_Graph6->GetXaxis()->SetRange(1,84);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph6->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi*_{#eta} ");
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph6->GetYaxis()->SetTitleOffset(1.05);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph6);
   
   grae->Draw("a2");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#6666ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6666ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.002,9.014114);
   grae->SetPointError(0,0,0,0.0239678,0.0239678);
   grae->SetPoint(1,0.006,8.858891);
   grae->SetPointError(1,0,0,0.02377138,0.02377138);
   grae->SetPoint(2,0.01,8.759078);
   grae->SetPointError(2,0,0,0.0235779,0.0235779);
   grae->SetPoint(3,0.014,8.611151);
   grae->SetPointError(3,0,0,0.02339807,0.02339807);
   grae->SetPoint(4,0.018,8.416117);
   grae->SetPointError(4,0,0,0.0231308,0.0231308);
   grae->SetPoint(5,0.022,8.100765);
   grae->SetPointError(5,0,0,0.02265084,0.02265084);
   grae->SetPoint(6,0.0265,7.779402);
   grae->SetPointError(6,0,0,0.01988439,0.01988439);
   grae->SetPoint(7,0.0315,7.368728);
   grae->SetPointError(7,0,0,0.01931979,0.01931979);
   grae->SetPoint(8,0.0365,6.925249);
   grae->SetPointError(8,0,0,0.01871791,0.01871791);
   grae->SetPoint(9,0.042,6.479816);
   grae->SetPointError(9,0,0,0.01652182,0.01652182);
   grae->SetPoint(10,0.0485,5.96542);
   grae->SetPointError(10,0,0,0.01465632,0.01465632);
   grae->SetPoint(11,0.0545,5.453526);
   grae->SetPointError(11,0,0,0.01657029,0.01657029);
   grae->SetPoint(12,0.0605,5.066608);
   grae->SetPointError(12,0,0,0.01349848,0.01349848);
   grae->SetPoint(13,0.068,4.56514);
   grae->SetPointError(13,0,0,0.01198234,0.01198234);
   grae->SetPoint(14,0.0765,4.085799);
   grae->SetPointError(14,0,0,0.01068761,0.01068761);
   grae->SetPoint(15,0.086,3.639835);
   grae->SetPointError(15,0,0,0.009572266,0.009572266);
   grae->SetPoint(16,0.0965,3.191935);
   grae->SetPointError(16,0,0,0.008545742,0.008545742);
   grae->SetPoint(17,0.108,2.795225);
   grae->SetPointError(17,0,0,0.007663977,0.007663977);
   grae->SetPoint(18,0.121,2.425971);
   grae->SetPointError(18,0,0,0.00661298,0.00661298);
   grae->SetPoint(19,0.1365,2.078896);
   grae->SetPointError(19,0,0,0.005560853,0.005560853);
   grae->SetPoint(20,0.155,1.730054);
   grae->SetPointError(20,0,0,0.004679248,0.004679248);
   grae->SetPoint(21,0.177,1.406145);
   grae->SetPointError(21,0,0,0.003851259,0.003851259);
   grae->SetPoint(22,0.204,1.120966);
   grae->SetPointError(22,0,0,0.003076839,0.003076839);
   grae->SetPoint(23,0.2385,0.8554184);
   grae->SetPointError(23,0,0,0.002361113,0.002361113);
   grae->SetPoint(24,0.285,0.621655);
   grae->SetPointError(24,0,0,0.001710575,0.001710575);
   grae->SetPoint(25,0.3515,0.4166713);
   grae->SetPointError(25,0,0,0.001157509,0.001157509);
   grae->SetPoint(26,0.4575,0.2421093);
   grae->SetPointError(26,0,0,0.0006814084,0.0006814084);
   grae->SetPoint(27,0.6095,0.1270682);
   grae->SetPointError(27,0,0,0.0004337708,0.0004337708);
   grae->SetPoint(28,0.8065,0.06354508);
   grae->SetPointError(28,0,0,0.0002683031,0.0002683031);
   grae->SetPoint(29,1.0355,0.03304646);
   grae->SetPointError(29,0,0,0.0001883912,0.0001883912);
   grae->SetPoint(30,1.3245,0.01741673);
   grae->SetPointError(30,0,0,0.0001131872,0.0001131872);
   grae->SetPoint(31,1.7215,0.008649548);
   grae->SetPointError(31,0,0,6.951573e-05,6.951573e-05);
   grae->SetPoint(32,2.2345,0.004196297);
   grae->SetPointError(32,0,0,4.288426e-05,4.288426e-05);
   grae->SetPoint(33,2.8995,0.002109578);
   grae->SetPointError(33,0,0,2.652624e-05,2.652624e-05);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,0.0018,3.18925);
   Graph_Graph7->SetMinimum(0);
   Graph_Graph7->SetMaximum(9.941682);
   Graph_Graph7->SetDirectory(0);
   Graph_Graph7->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph7->SetLineColor(ci);
   Graph_Graph7->GetXaxis()->SetLabelFont(42);
   Graph_Graph7->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleFont(42);
   Graph_Graph7->GetYaxis()->SetLabelFont(42);
   Graph_Graph7->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleFont(42);
   Graph_Graph7->GetZaxis()->SetLabelFont(42);
   Graph_Graph7->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph7);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.002,9.106297);
   grae->SetPointError(0,0,0,0.040844,0.040844);
   grae->SetPoint(1,0.006,9.078291);
   grae->SetPointError(1,0,0,0.04207356,0.04207356);
   grae->SetPoint(2,0.01,8.860786);
   grae->SetPointError(2,0,0,0.04106303,0.04106303);
   grae->SetPoint(3,0.014,8.676274);
   grae->SetPointError(3,0,0,0.04062779,0.04062779);
   grae->SetPoint(4,0.018,8.337609);
   grae->SetPointError(4,0,0,0.04030721,0.04030721);
   grae->SetPoint(5,0.022,8.044422);
   grae->SetPointError(5,0,0,0.0396631,0.0396631);
   grae->SetPoint(6,0.0265,7.659918);
   grae->SetPointError(6,0,0,0.03422848,0.03422848);
   grae->SetPoint(7,0.0315,7.17912);
   grae->SetPointError(7,0,0,0.03234194,0.03234194);
   grae->SetPoint(8,0.0365,6.728084);
   grae->SetPointError(8,0,0,0.03201601,0.03201601);
   grae->SetPoint(9,0.042,6.282668);
   grae->SetPointError(9,0,0,0.02743313,0.02743313);
   grae->SetPoint(10,0.0485,5.775932);
   grae->SetPointError(10,0,0,0.02472142,0.02472142);
   grae->SetPoint(11,0.0545,5.332089);
   grae->SetPointError(11,0,0,0.0278476,0.0278476);
   grae->SetPoint(12,0.0605,4.923704);
   grae->SetPointError(12,0,0,0.02277378,0.02277378);
   grae->SetPoint(13,0.068,4.449297);
   grae->SetPointError(13,0,0,0.0205291,0.0205291);
   grae->SetPoint(14,0.0765,4.011493);
   grae->SetPointError(14,0,0,0.01784722,0.01784722);
   grae->SetPoint(15,0.086,3.564145);
   grae->SetPointError(15,0,0,0.01598776,0.01598776);
   grae->SetPoint(16,0.0965,3.159969);
   grae->SetPointError(16,0,0,0.01433522,0.01433522);
   grae->SetPoint(17,0.108,2.77742);
   grae->SetPointError(17,0,0,0.01291623,0.01291623);
   grae->SetPoint(18,0.121,2.420708);
   grae->SetPointError(18,0,0,0.01126212,0.01126212);
   grae->SetPoint(19,0.1365,2.071115);
   grae->SetPointError(19,0,0,0.00939609,0.00939609);
   grae->SetPoint(20,0.155,1.734207);
   grae->SetPointError(20,0,0,0.007871789,0.007871789);
   grae->SetPoint(21,0.177,1.424907);
   grae->SetPointError(21,0,0,0.006580181,0.006580181);
   grae->SetPoint(22,0.204,1.133518);
   grae->SetPointError(22,0,0,0.005234261,0.005234261);
   grae->SetPoint(23,0.2385,0.876128);
   grae->SetPointError(23,0,0,0.004004305,0.004004305);
   grae->SetPoint(24,0.285,0.6394784);
   grae->SetPointError(24,0,0,0.002925678,0.002925678);
   grae->SetPoint(25,0.3515,0.4301066);
   grae->SetPointError(25,0,0,0.001966129,0.001966129);
   grae->SetPoint(26,0.4575,0.2519978);
   grae->SetPointError(26,0,0,0.001147313,0.001147313);
   grae->SetPoint(27,0.6095,0.1319807);
   grae->SetPointError(27,0,0,0.000738591,0.000738591);
   grae->SetPoint(28,0.8065,0.06668263);
   grae->SetPointError(28,0,0,0.000472941,0.000472941);
   grae->SetPoint(29,1.0355,0.0345134);
   grae->SetPointError(29,0,0,0.0003363745,0.0003363745);
   grae->SetPoint(30,1.3245,0.01800679);
   grae->SetPointError(30,0,0,0.0001993801,0.0001993801);
   grae->SetPoint(31,1.7215,0.00856106);
   grae->SetPointError(31,0,0,0.0001151668,0.0001151668);
   grae->SetPoint(32,2.2345,0.00417938);
   grae->SetPointError(32,0,0,7.548669e-05,7.548669e-05);
   grae->SetPoint(33,2.8995,0.002152797);
   grae->SetPointError(33,0,0,5.124535e-05,5.124535e-05);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Graph",100,0.0018,3.18925);
   Graph_Graph8->SetMinimum(0);
   Graph_Graph8->SetMaximum(10.06165);
   Graph_Graph8->SetDirectory(0);
   Graph_Graph8->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph8->SetLineColor(ci);
   Graph_Graph8->GetXaxis()->SetLabelFont(42);
   Graph_Graph8->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleFont(42);
   Graph_Graph8->GetYaxis()->SetLabelFont(42);
   Graph_Graph8->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetYaxis()->SetTitleFont(42);
   Graph_Graph8->GetZaxis()->SetLabelFont(42);
   Graph_Graph8->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph8);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.53,0.77,0.9,0.91,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph2","2012 data","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Z#rightarrow ll MadGraph","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6666ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.745,0.95,"19.8 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.895,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.3,"#sqrt{s} = 8 TeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.23,"|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.16,"p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.09,"60 GeV < M_{ll} < 120 GeV");
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
   p2->Range(-3.924951,0.7367742,0.7623949,1.123871);
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

   ci = TColor::GetColor("#6666ff");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#6666ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6666ff");
   grae->SetMarkerColor(ci);
   grae->SetPoint(0,0.002,1);
   grae->SetPointError(0,0.002,0.002,0.002658919,0.002658919);
   grae->SetPoint(1,0.006,1);
   grae->SetPointError(1,0.002,0.002,0.002683336,0.002683336);
   grae->SetPoint(2,0.01,1);
   grae->SetPointError(2,0.002,0.002,0.002691824,0.002691824);
   grae->SetPoint(3,0.014,1);
   grae->SetPointError(3,0.002,0.002,0.002717182,0.002717182);
   grae->SetPoint(4,0.018,1);
   grae->SetPointError(4,0.002,0.002,0.002748393,0.002748393);
   grae->SetPoint(5,0.022,1);
   grae->SetPointError(5,0.002,0.002,0.002796135,0.002796135);
   grae->SetPoint(6,0.0265,1);
   grae->SetPointError(6,0.0025,0.0025,0.00255603,0.00255603);
   grae->SetPoint(7,0.0315,1);
   grae->SetPointError(7,0.0025,0.0025,0.002621862,0.002621862);
   grae->SetPoint(8,0.0365,1);
   grae->SetPointError(8,0.0025,0.0025,0.00270285,0.00270285);
   grae->SetPoint(9,0.042,1);
   grae->SetPointError(9,0.003,0.003,0.002549737,0.002549737);
   grae->SetPoint(10,0.0485,1);
   grae->SetPointError(10,0.0035,0.0035,0.00245688,0.00245688);
   grae->SetPoint(11,0.0545,1);
   grae->SetPointError(11,0.0025,0.0025,0.003038455,0.003038455);
   grae->SetPoint(12,0.0605,1);
   grae->SetPointError(12,0.0035,0.0035,0.002664204,0.002664204);
   grae->SetPoint(13,0.068,1);
   grae->SetPointError(13,0.004,0.004,0.002624748,0.002624748);
   grae->SetPoint(14,0.0765,1);
   grae->SetPointError(14,0.0045,0.0045,0.002615794,0.002615794);
   grae->SetPoint(15,0.086,1);
   grae->SetPointError(15,0.005,0.005,0.002629863,0.002629863);
   grae->SetPoint(16,0.0965,1);
   grae->SetPointError(16,0.0055,0.0055,0.002677292,0.002677292);
   grae->SetPoint(17,0.108,1);
   grae->SetPointError(17,0.006,0.006,0.00274181,0.00274181);
   grae->SetPoint(18,0.121,1);
   grae->SetPointError(18,0.007,0.007,0.002725911,0.002725911);
   grae->SetPoint(19,0.1365,1);
   grae->SetPointError(19,0.0085,0.0085,0.002674907,0.002674907);
   grae->SetPoint(20,0.155,1);
   grae->SetPointError(20,0.01,0.01,0.002704683,0.002704683);
   grae->SetPoint(21,0.177,1);
   grae->SetPointError(21,0.012,0.012,0.002738877,0.002738877);
   grae->SetPoint(22,0.204,1);
   grae->SetPointError(22,0.015,0.015,0.002744809,0.002744809);
   grae->SetPoint(23,0.2385,1);
   grae->SetPointError(23,0.0195,0.0195,0.002760186,0.002760186);
   grae->SetPoint(24,0.285,1);
   grae->SetPointError(24,0.027,0.027,0.002751646,0.002751646);
   grae->SetPoint(25,0.3515,1);
   grae->SetPointError(25,0.0395,0.0395,0.002777991,0.002777991);
   grae->SetPoint(26,0.4575,1);
   grae->SetPointError(26,0.0665,0.0665,0.002814465,0.002814465);
   grae->SetPoint(27,0.6095,1);
   grae->SetPointError(27,0.0855,0.0855,0.003413684,0.003413684);
   grae->SetPoint(28,0.8065,1);
   grae->SetPointError(28,0.1115,0.1115,0.004222249,0.004222249);
   grae->SetPoint(29,1.0355,1);
   grae->SetPointError(29,0.1175,0.1175,0.0057008,0.0057008);
   grae->SetPoint(30,1.3245,1);
   grae->SetPointError(30,0.1715,0.1715,0.006498765,0.006498765);
   grae->SetPoint(31,1.7215,1);
   grae->SetPointError(31,0.2255,0.2255,0.00803692,0.00803692);
   grae->SetPoint(32,2.2345,1);
   grae->SetPointError(32,0.2875,0.2875,0.01021955,0.01021955);
   grae->SetPoint(33,2.8995,1);
   grae->SetPointError(33,0.3775,0.3775,0.01257419,0.01257419);
   
   TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","",100,0.0006,3.604633);
   Graph_Graph9->SetMinimum(0.88);
   Graph_Graph9->SetMaximum(1.12);
   Graph_Graph9->SetDirectory(0);
   Graph_Graph9->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph9->SetLineColor(ci);
   Graph_Graph9->GetXaxis()->SetTitle("#phi*");
   Graph_Graph9->GetXaxis()->SetRange(1,84);
   Graph_Graph9->GetXaxis()->SetLabelFont(42);
   Graph_Graph9->GetXaxis()->SetLabelSize(0.12);
   Graph_Graph9->GetXaxis()->SetTitleSize(0.12);
   Graph_Graph9->GetXaxis()->SetTitleOffset(1.05);
   Graph_Graph9->GetXaxis()->SetTitleFont(42);
   Graph_Graph9->GetYaxis()->SetTitle("Data/MC");
   Graph_Graph9->GetYaxis()->SetNdivisions(503);
   Graph_Graph9->GetYaxis()->SetLabelFont(42);
   Graph_Graph9->GetYaxis()->SetLabelSize(0.12);
   Graph_Graph9->GetYaxis()->SetTitleSize(0.12);
   Graph_Graph9->GetYaxis()->SetTitleOffset(0.4);
   Graph_Graph9->GetYaxis()->SetTitleFont(42);
   Graph_Graph9->GetZaxis()->SetLabelFont(42);
   Graph_Graph9->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph9);
   
   grae->Draw("ae2");
   
   grae = new TGraphAsymmErrors(34);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.002,1.010227);
   grae->SetPointError(0,0,0,0.004531116,0.004531116);
   grae->SetPoint(1,0.006,1.024766);
   grae->SetPointError(1,0,0,0.004749303,0.004749303);
   grae->SetPoint(2,0.01,1.011612);
   grae->SetPointError(2,0,0,0.004688053,0.004688053);
   grae->SetPoint(3,0.014,1.007563);
   grae->SetPointError(3,0,0,0.004718044,0.004718044);
   grae->SetPoint(4,0.018,0.9906717);
   grae->SetPointError(4,0,0,0.004789289,0.004789289);
   grae->SetPoint(5,0.022,0.9930446);
   grae->SetPointError(5,0,0,0.004896217,0.004896217);
   grae->SetPoint(6,0.0265,0.984641);
   grae->SetPointError(6,0,0,0.004399887,0.004399887);
   grae->SetPoint(7,0.0315,0.9742685);
   grae->SetPointError(7,0,0,0.00438908,0.00438908);
   grae->SetPoint(8,0.0365,0.9715296);
   grae->SetPointError(8,0,0,0.004623085,0.004623085);
   grae->SetPoint(9,0.042,0.969575);
   grae->SetPointError(9,0,0,0.004233628,0.004233628);
   grae->SetPoint(10,0.0485,0.9682357);
   grae->SetPointError(10,0,0,0.00414412,0.00414412);
   grae->SetPoint(11,0.0545,0.9777324);
   grae->SetPointError(11,0,0,0.005106347,0.005106347);
   grae->SetPoint(12,0.0605,0.971795);
   grae->SetPointError(12,0,0,0.004494877,0.004494877);
   grae->SetPoint(13,0.068,0.9746244);
   grae->SetPointError(13,0,0,0.004496926,0.004496926);
   grae->SetPoint(14,0.0765,0.9818136);
   grae->SetPointError(14,0,0,0.00436811,0.00436811);
   grae->SetPoint(15,0.086,0.979205);
   grae->SetPointError(15,0,0,0.004392441,0.004392441);
   grae->SetPoint(16,0.0965,0.9899855);
   grae->SetPointError(16,0,0,0.004491074,0.004491074);
   grae->SetPoint(17,0.108,0.99363);
   grae->SetPointError(17,0,0,0.004620821,0.004620821);
   grae->SetPoint(18,0.121,0.9978308);
   grae->SetPointError(18,0,0,0.004642313,0.004642313);
   grae->SetPoint(19,0.1365,0.9962573);
   grae->SetPointError(19,0,0,0.004519751,0.004519751);
   grae->SetPoint(20,0.155,1.0024);
   grae->SetPointError(20,0,0,0.004550024,0.004550024);
   grae->SetPoint(21,0.177,1.013342);
   grae->SetPointError(21,0,0,0.004679588,0.004679588);
   grae->SetPoint(22,0.204,1.011197);
   grae->SetPointError(22,0,0,0.004669418,0.004669418);
   grae->SetPoint(23,0.2385,1.02421);
   grae->SetPointError(23,0,0,0.004681108,0.004681108);
   grae->SetPoint(24,0.285,1.028671);
   grae->SetPointError(24,0,0,0.004706273,0.004706273);
   grae->SetPoint(25,0.3515,1.032244);
   grae->SetPointError(25,0,0,0.004718657,0.004718657);
   grae->SetPoint(26,0.4575,1.040843);
   grae->SetPointError(26,0,0,0.004738822,0.004738822);
   grae->SetPoint(27,0.6095,1.03866);
   grae->SetPointError(27,0,0,0.005812554,0.005812554);
   grae->SetPoint(28,0.8065,1.049375);
   grae->SetPointError(28,0,0,0.007442606,0.007442606);
   grae->SetPoint(29,1.0355,1.04439);
   grae->SetPointError(29,0,0,0.01017884,0.01017884);
   grae->SetPoint(30,1.3245,1.033879);
   grae->SetPointError(30,0,0,0.01144762,0.01144762);
   grae->SetPoint(31,1.7215,0.9897696);
   grae->SetPointError(31,0,0,0.01331478,0.01331478);
   grae->SetPoint(32,2.2345,0.9959687);
   grae->SetPointError(32,0,0,0.01798888,0.01798888);
   grae->SetPoint(33,2.8995,1.020487);
   grae->SetPointError(33,0,0,0.02429175,0.02429175);
   
   TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","Graph",100,0.0018,3.18925);
   Graph_Graph10->SetMinimum(0.954819);
   Graph_Graph10->SetMaximum(1.06609);
   Graph_Graph10->SetDirectory(0);
   Graph_Graph10->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph10->SetLineColor(ci);
   Graph_Graph10->GetXaxis()->SetLabelFont(42);
   Graph_Graph10->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetXaxis()->SetTitleFont(42);
   Graph_Graph10->GetYaxis()->SetLabelFont(42);
   Graph_Graph10->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetYaxis()->SetTitleFont(42);
   Graph_Graph10->GetZaxis()->SetLabelFont(42);
   Graph_Graph10->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph10);
   
   grae->Draw("pe");
   p2->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->SetSelected(FinalPhiTot);
}
