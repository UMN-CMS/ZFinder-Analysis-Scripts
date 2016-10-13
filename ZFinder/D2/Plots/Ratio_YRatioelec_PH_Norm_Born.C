{
//=========Macro generated from canvas: FinalPhiRatio/FinalPhiRatio
//=========  (Wed Oct  5 11:18:58 2016) by ROOT version5.34/18
   TCanvas *FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio",52,80,450,896);
   FinalPhiRatio->Range(-0.3003,0.40625,2.7027,1.34375);
   FinalPhiRatio->SetFillColor(0);
   FinalPhiRatio->SetBorderMode(0);
   FinalPhiRatio->SetBorderSize(2);
   FinalPhiRatio->SetFrameBorderMode(0);
   FinalPhiRatio->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetPoint(0,0.2,1);
   grae->SetPointError(0,0.2,0.2,0.001471333,0.001471333);
   grae->SetPoint(1,0.6,1);
   grae->SetPointError(1,0.2,0.2,0.00145628,0.00145628);
   grae->SetPoint(2,1,1);
   grae->SetPointError(2,0.2,0.2,0.001587423,0.001587423);
   grae->SetPoint(3,1.4,1);
   grae->SetPointError(3,0.2,0.2,0.001854211,0.001854211);
   grae->SetPoint(4,1.8,1);
   grae->SetPointError(4,0.2,0.2,0.002667337,0.002667337);
   grae->SetPoint(5,2.2,1);
   grae->SetPointError(5,0.2,0.2,0.005826074,0.005826074);
   
   TH1F *Graph_Graph_Graph58 = new TH1F("Graph_Graph_Graph58","",100,0,2.64);
   Graph_Graph_Graph58->SetMinimum(0.5);
   Graph_Graph_Graph58->SetMaximum(1.25);
   Graph_Graph_Graph58->SetDirectory(0);
   Graph_Graph_Graph58->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph58->SetLineColor(ci);
   Graph_Graph_Graph58->GetXaxis()->SetTitle("y_{ll}");
   Graph_Graph_Graph58->GetXaxis()->SetRange(1,91);
   Graph_Graph_Graph58->GetXaxis()->CenterTitle(true);
   Graph_Graph_Graph58->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetXaxis()->SetLabelOffset(0.012);
   Graph_Graph_Graph58->GetXaxis()->SetTitleSize(0.043);
   Graph_Graph_Graph58->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph58->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph_Graph58->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph58->GetYaxis()->SetNdivisions(505);
   Graph_Graph_Graph58->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph58->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph58->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph58->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph58->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph58);
   
   grae->Draw("ae2");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.2,0.9972603);
   grae->SetPointError(0,0,0,0.002529146,0.002529146);
   grae->SetPoint(1,0.6,1.001735);
   grae->SetPointError(1,0,0,0.001771153,0.001771153);
   grae->SetPoint(2,1,0.9974508);
   grae->SetPointError(2,0,0,0.0008238954,0.0008238954);
   grae->SetPoint(3,1.4,0.9954442);
   grae->SetPointError(3,0,0,0.001322895,0.001322895);
   grae->SetPoint(4,1.8,1.01187);
   grae->SetPointError(4,0,0,0.002476914,0.002476914);
   grae->SetPoint(5,2.2,1.006533);
   grae->SetPointError(5,0,0,0.003735747,0.003735747);
   
   TH1F *Graph_Graph_Graph79 = new TH1F("Graph_Graph_Graph79","Graph",100,0,2.4);
   Graph_Graph_Graph79->SetMinimum(0.9920988);
   Graph_Graph_Graph79->SetMaximum(1.016369);
   Graph_Graph_Graph79->SetDirectory(0);
   Graph_Graph_Graph79->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph79->SetLineColor(ci);
   Graph_Graph_Graph79->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph79->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph79->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph79);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#6666ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6666ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetPoint(0,0.2,0.9727288);
   grae->SetPointError(0,0,0,0.001100537,0.001100537);
   grae->SetPoint(1,0.6,0.9841751);
   grae->SetPointError(1,0,0,0.001115194,0.001115194);
   grae->SetPoint(2,1,0.9921097);
   grae->SetPointError(2,0,0,0.00113392,0.00113392);
   grae->SetPoint(3,1.4,1.013662);
   grae->SetPointError(3,0,0,0.001241168,0.001241168);
   grae->SetPoint(4,1.8,1.054109);
   grae->SetPointError(4,0,0,0.001566015,0.001566015);
   grae->SetPoint(5,2.2,1.079226);
   grae->SetPointError(5,0,0,0.003238944,0.003238944);
   
   TH1F *Graph_Graph_Graph610 = new TH1F("Graph_Graph_Graph610","Graph",100,0,2.4);
   Graph_Graph_Graph610->SetMinimum(0.9605445);
   Graph_Graph_Graph610->SetMaximum(1.093549);
   Graph_Graph_Graph610->SetDirectory(0);
   Graph_Graph_Graph610->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph610->SetLineColor(ci);
   Graph_Graph_Graph610->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph610->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph610->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph610);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#00cc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00cc00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(3);
   grae->SetPoint(0,0.2,1.053594);
   grae->SetPointError(0,0,0,9.601556e-06,9.601556e-06);
   grae->SetPoint(1,0.6,1.056796);
   grae->SetPointError(1,0,0,9.745491e-06,9.745491e-06);
   grae->SetPoint(2,1,1.051799);
   grae->SetPointError(2,0,0,1.000135e-05,1.000135e-05);
   grae->SetPoint(3,1.4,0.9937728);
   grae->SetPointError(3,0,0,1.171517e-05,1.171517e-05);
   grae->SetPoint(4,1.8,0.8282453);
   grae->SetPointError(4,0,0,1.78902e-05,1.78902e-05);
   grae->SetPoint(5,2.2,0.5347542);
   grae->SetPointError(5,0,0,7.45738e-05,7.45738e-05);
   
   TH1F *Graph_Graph11 = new TH1F("Graph_Graph11","Graph",100,0,2.4);
   Graph_Graph11->SetMinimum(0.482467);
   Graph_Graph11->SetMaximum(1.109018);
   Graph_Graph11->SetDirectory(0);
   Graph_Graph11->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph11->SetLineColor(ci);
   Graph_Graph11->GetXaxis()->SetLabelFont(42);
   Graph_Graph11->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleFont(42);
   Graph_Graph11->GetYaxis()->SetLabelFont(42);
   Graph_Graph11->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleFont(42);
   Graph_Graph11->GetZaxis()->SetLabelFont(42);
   Graph_Graph11->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph11);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#009999");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#009999");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(24);
   grae->SetPoint(0,0.2,0.9885504);
   grae->SetPointError(0,0,0,9.601556e-06,9.601556e-06);
   grae->SetPoint(1,0.6,0.9952626);
   grae->SetPointError(1,0,0,9.745491e-06,9.745491e-06);
   grae->SetPoint(2,1,0.9962789);
   grae->SetPointError(2,0,0,1.000135e-05,1.000135e-05);
   grae->SetPoint(3,1.4,0.9979535);
   grae->SetPointError(3,0,0,1.171517e-05,1.171517e-05);
   grae->SetPoint(4,1.8,1.028673);
   grae->SetPointError(4,0,0,1.78902e-05,1.78902e-05);
   grae->SetPoint(5,2.2,1.036644);
   grae->SetPointError(5,0,0,7.45738e-05,7.45738e-05);
   
   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,0,2.4);
   Graph_Graph12->SetMinimum(0.983723);
   Graph_Graph12->SetMaximum(1.041536);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph12);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(25);
   grae->SetPoint(0,0.2,0.9964778);
   grae->SetPointError(0,0,0,9.601556e-06,9.601556e-06);
   grae->SetPoint(1,0.6,1.001509);
   grae->SetPointError(1,0,0,9.745491e-06,9.745491e-06);
   grae->SetPoint(2,1,0.9973628);
   grae->SetPointError(2,0,0,1.000135e-05,1.000135e-05);
   grae->SetPoint(3,1.4,0.997847);
   grae->SetPointError(3,0,0,1.171517e-05,1.171517e-05);
   grae->SetPoint(4,1.8,1.013896);
   grae->SetPointError(4,0,0,1.78902e-05,1.78902e-05);
   grae->SetPoint(5,2.2,1.007439);
   grae->SetPointError(5,0,0,7.45738e-05,7.45738e-05);
   
   TH1F *Graph_Graph13 = new TH1F("Graph_Graph13","Graph",100,0,2.4);
   Graph_Graph13->SetMinimum(0.9947237);
   Graph_Graph13->SetMaximum(1.015658);
   Graph_Graph13->SetDirectory(0);
   Graph_Graph13->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph13->SetLineColor(ci);
   Graph_Graph13->GetXaxis()->SetLabelFont(42);
   Graph_Graph13->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleFont(42);
   Graph_Graph13->GetYaxis()->SetLabelFont(42);
   Graph_Graph13->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleFont(42);
   Graph_Graph13->GetZaxis()->SetLabelFont(42);
   Graph_Graph13->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph13);
   
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.7,0.907,"19.7 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.907,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.25,"|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.2,"p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.15,"60 GeV < M_{ll} < 120 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.13,0.72,0.7,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(22);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","2012 data","F");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph1"," MadGraph+Pythia6 (Z2*)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6666ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph2","POWHEG+Pythia6 (Z2*)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph0","Resbos","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00cc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(3);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph1","AMC@nlo+Pythia8(CUETP8M1)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009999");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph2","POWHEG+Pythia8 (CT10)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   leg->Draw();
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","",100,0,2.64);
   Graph_Graph5->SetMinimum(0.5);
   Graph_Graph5->SetMaximum(1.25);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph5->SetLineColor(ci);
   Graph_Graph5->GetXaxis()->SetTitle("y_{ll}");
   Graph_Graph5->GetXaxis()->SetRange(1,91);
   Graph_Graph5->GetXaxis()->CenterTitle(true);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelOffset(0.012);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.043);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph5->GetYaxis()->CenterTitle(true);
   Graph_Graph5->GetYaxis()->SetNdivisions(505);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   Graph_Graph5->Draw("sameaxis");
   FinalPhiRatio->Modified();
   FinalPhiRatio->cd();
   FinalPhiRatio->SetSelected(FinalPhiRatio);
}
