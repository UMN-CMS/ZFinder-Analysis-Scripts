{
//=========Macro generated from canvas: BinMigration_MP/BinMigration_MP
//=========  (Wed Apr  1 12:33:03 2015) by ROOT version5.34/18
   TCanvas *BinMigration_MP = new TCanvas("BinMigration_MP", "BinMigration_MP",56,79,800,900);
   BinMigration_MP->Range(-3.151198,0.75,0.9135106,1.25);
   BinMigration_MP->SetFillColor(0);
   BinMigration_MP->SetBorderMode(0);
   BinMigration_MP->SetBorderSize(2);
   BinMigration_MP->SetLogx();
   BinMigration_MP->SetFrameBorderMode(0);
   BinMigration_MP->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(35);
   grae->SetName("Graph0");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.002,1.00325);
   grae->SetPointError(0,0,0,0.02839551,0.02839551);
   grae->SetPoint(1,0.006,0.9952978);
   grae->SetPointError(1,0,0,0.02993856,0.02993856);
   grae->SetPoint(2,0.01,1.018079);
   grae->SetPointError(2,0,0,0.02988458,0.02988458);
   grae->SetPoint(3,0.014,0.9956378);
   grae->SetPointError(3,0,0,0.0300566,0.0300566);
   grae->SetPoint(4,0.018,0.982187);
   grae->SetPointError(4,0,0,0.03082209,0.03082209);
   grae->SetPoint(5,0.022,0.9902705);
   grae->SetPointError(5,0,0,0.03018221,0.03018221);
   grae->SetPoint(6,0.0265,1.019899);
   grae->SetPointError(6,0,0,0.02768277,0.02768277);
   grae->SetPoint(7,0.0315,0.9849935);
   grae->SetPointError(7,0,0,0.02703851,0.02703851);
   grae->SetPoint(8,0.0365,1.001824);
   grae->SetPointError(8,0,0,0.02906803,0.02906803);
   grae->SetPoint(9,0.042,0.99336);
   grae->SetPointError(9,0,0,0.02624488,0.02624488);
   grae->SetPoint(10,0.0485,1.00638);
   grae->SetPointError(10,0,0,0.02459282,0.02459282);
   grae->SetPoint(11,0.0545,1.00328);
   grae->SetPointError(11,0,0,0.03128612,0.03128612);
   grae->SetPoint(12,0.0605,1.006335);
   grae->SetPointError(12,0,0,0.02661772,0.02661772);
   grae->SetPoint(13,0.068,0.9896944);
   grae->SetPointError(13,0,0,0.02567866,0.02567866);
   grae->SetPoint(14,0.0765,1.012878);
   grae->SetPointError(14,0,0,0.02578888,0.02578888);
   grae->SetPoint(15,0.086,0.9945516);
   grae->SetPointError(15,0,0,0.02532306,0.02532306);
   grae->SetPoint(16,0.0965,1.007028);
   grae->SetPointError(16,0,0,0.02665868,0.02665868);
   grae->SetPoint(17,0.108,0.9875582);
   grae->SetPointError(17,0,0,0.02656513,0.02656513);
   grae->SetPoint(18,0.121,1.000975);
   grae->SetPointError(18,0,0,0.02648634,0.02648634);
   grae->SetPoint(19,0.1365,1.004351);
   grae->SetPointError(19,0,0,0.02633074,0.02633074);
   grae->SetPoint(20,0.155,1.005233);
   grae->SetPointError(20,0,0,0.0261118,0.0261118);
   grae->SetPoint(21,0.177,0.9932978);
   grae->SetPointError(21,0,0,0.02607367,0.02607367);
   grae->SetPoint(22,0.204,1.003765);
   grae->SetPointError(22,0,0,0.02622969,0.02622969);
   grae->SetPoint(23,0.2385,0.9933049);
   grae->SetPointError(23,0,0,0.02621693,0.02621693);
   grae->SetPoint(24,0.285,1.002918);
   grae->SetPointError(24,0,0,0.02608952,0.02608952);
   grae->SetPoint(25,0.3515,1.001113);
   grae->SetPointError(25,0,0,0.02621762,0.02621762);
   grae->SetPoint(26,0.4575,1.001119);
   grae->SetPointError(26,0,0,0.02659177,0.02659177);
   grae->SetPoint(27,0.6095,0.9988394);
   grae->SetPointError(27,0,0,0.03274974,0.03274974);
   grae->SetPoint(28,0.8065,0.9958707);
   grae->SetPointError(28,0,0,0.04067625,0.04067625);
   grae->SetPoint(29,1.0355,1.002392);
   grae->SetPointError(29,0,0,0.05886637,0.05886637);
   grae->SetPoint(30,1.3245,1.001956);
   grae->SetPointError(30,0,0,0.0635891,0.0635891);
   grae->SetPoint(31,1.7215,1.006134);
   grae->SetPointError(31,0,0,0.08337461,0.08337461);
   grae->SetPoint(32,2.2345,0.9914842);
   grae->SetPointError(32,0,0,0.1054447,0.1054447);
   grae->SetPoint(33,2.8995,0.9986668);
   grae->SetPointError(33,0,0,0.1367861,0.1367861);
   grae->SetPoint(34,6.6385,1.00008);
   grae->SetPointError(34,0,0,0.08064454,0.08064454);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0.0018,7.30215);
   Graph_Graph1->SetMinimum(0.8);
   Graph_Graph1->SetMaximum(1.2);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetTitle("#phi*");
   Graph_Graph1->GetXaxis()->SetRange(0,44);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(-0.01);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.8);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("Unfolded/Generated");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.03);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1);
   
   grae->Draw("pa");
   
   grae = new TGraphAsymmErrors(35);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineColor(2);
   grae->SetMarkerColor(2);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,0.002,1.005501);
   grae->SetPointError(0,0,0,0.02939741,0.03003413);
   grae->SetPoint(1,0.006,0.9946167);
   grae->SetPointError(1,0,0,0.02753389,0.02831055);
   grae->SetPoint(2,0.01,1.018079);
   grae->SetPointError(2,0,0,0.02947236,0.02887356);
   grae->SetPoint(3,0.014,0.995831);
   grae->SetPointError(3,0,0,0.02559223,0.0273592);
   grae->SetPoint(4,0.018,0.9818425);
   grae->SetPointError(4,0,0,0.031652,0.02889034);
   grae->SetPoint(5,0.022,0.9899925);
   grae->SetPointError(5,0,0,0.02947727,0.02790075);
   grae->SetPoint(6,0.0265,1.02277);
   grae->SetPointError(6,0,0,0.03042974,0.02558654);
   grae->SetPoint(7,0.0315,0.9855716);
   grae->SetPointError(7,0,0,0.02801516,0.02864083);
   grae->SetPoint(8,0.0365,1.002003);
   grae->SetPointError(8,0,0,0.02832047,0.02883147);
   grae->SetPoint(9,0.042,0.9939141);
   grae->SetPointError(9,0,0,0.02522137,0.02639494);
   grae->SetPoint(10,0.0485,1.009867);
   grae->SetPointError(10,0,0,0.02345891,0.02264355);
   grae->SetPoint(11,0.0545,1.000315);
   grae->SetPointError(11,0,0,0.03114619,0.03037806);
   grae->SetPoint(12,0.0605,1.006324);
   grae->SetPointError(12,0,0,0.02519001,0.02889961);
   grae->SetPoint(13,0.068,0.9903491);
   grae->SetPointError(13,0,0,0.0259888,0.02459971);
   grae->SetPoint(14,0.0765,1.012911);
   grae->SetPointError(14,0,0,0.02433978,0.02779281);
   grae->SetPoint(15,0.086,0.9921953);
   grae->SetPointError(15,0,0,0.0248095,0.02471185);
   grae->SetPoint(16,0.0965,1.007497);
   grae->SetPointError(16,0,0,0.02349752,0.02865718);
   grae->SetPoint(17,0.108,0.9895114);
   grae->SetPointError(17,0,0,0.02547316,0.02594816);
   grae->SetPoint(18,0.121,1.001336);
   grae->SetPointError(18,0,0,0.02618233,0.02937608);
   grae->SetPoint(19,0.1365,1.003574);
   grae->SetPointError(19,0,0,0.02745523,0.02555906);
   grae->SetPoint(20,0.155,1.00444);
   grae->SetPointError(20,0,0,0.02437915,0.02404638);
   grae->SetPoint(21,0.177,0.9934361);
   grae->SetPointError(21,0,0,0.02671766,0.02679438);
   grae->SetPoint(22,0.204,1.004926);
   grae->SetPointError(22,0,0,0.02692655,0.02333326);
   grae->SetPoint(23,0.2385,0.9923631);
   grae->SetPointError(23,0,0,0.02432763,0.02671954);
   grae->SetPoint(24,0.285,1.003234);
   grae->SetPointError(24,0,0,0.03021711,0.02776293);
   grae->SetPoint(25,0.3515,1.001353);
   grae->SetPointError(25,0,0,0.0256878,0.02758108);
   grae->SetPoint(26,0.4575,1.001261);
   grae->SetPointError(26,0,0,0.0259895,0.02808484);
   grae->SetPoint(27,0.6095,0.9978241);
   grae->SetPointError(27,0,0,0.03223191,0.03293131);
   grae->SetPoint(28,0.8065,0.9960565);
   grae->SetPointError(28,0,0,0.04282531,0.03905719);
   grae->SetPoint(29,1.0355,1.000828);
   grae->SetPointError(29,0,0,0.06159048,0.05995788);
   grae->SetPoint(30,1.3245,0.99688);
   grae->SetPointError(30,0,0,0.06136981,0.06500153);
   grae->SetPoint(31,1.7215,1.004182);
   grae->SetPointError(31,0,0,0.08695831,0.08473954);
   grae->SetPoint(32,2.2345,0.9942198);
   grae->SetPointError(32,0,0,0.1054604,0.1033279);
   grae->SetPoint(33,2.8995,1.015166);
   grae->SetPointError(33,0,0,0.1546953,0.1313852);
   grae->SetPoint(34,6.6385,1.000994);
   grae->SetPointError(34,0,0,0.07688533,0.07975831);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0.0018,7.30215);
   Graph_Graph2->SetMinimum(0.8318623);
   Graph_Graph2->SetMaximum(1.175159);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph2);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.45,0.77,0.85,0.91,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","RooUnfold","PL");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Toy MC","PL");
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.19,0.21,"50000 Powheg events");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.17,"unfolded using MadGraph");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   BinMigration_MP->Modified();
   BinMigration_MP->cd();
   BinMigration_MP->SetSelected(BinMigration_MP);
}
