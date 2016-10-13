{
//=========Macro generated from canvas: c1/c1
//=========  (Mon May  2 15:28:26 2016) by ROOT version5.34/18
   TCanvas *c1 = new TCanvas("c1", "c1",12,230,700,500);
   gStyle->SetOptStat(0);
   c1->Range(-1.127451,-0.133687,10.2451,1.180902);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetRightMargin(0.02155172);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *h = new TH1F("h","h",1,0,10);
   h->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   h->SetLineColor(ci);
   h->GetXaxis()->SetLabelFont(42);
   h->GetXaxis()->SetLabelSize(0.035);
   h->GetXaxis()->SetTitleSize(0.035);
   h->GetXaxis()->SetTitleFont(42);
   h->GetYaxis()->SetLabelFont(42);
   h->GetYaxis()->SetLabelSize(0.035);
   h->GetYaxis()->SetTitleSize(0.035);
   h->GetYaxis()->SetTitleFont(42);
   h->GetZaxis()->SetLabelFont(42);
   h->GetZaxis()->SetLabelSize(0.035);
   h->GetZaxis()->SetTitleSize(0.035);
   h->GetZaxis()->SetTitleFont(42);
   h->Draw("");
   
   TPaveText *pt = new TPaveText(0.4799425,0.94,0.5200575,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("h");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
