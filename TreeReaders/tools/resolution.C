{
  //  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile*  hFile = new TFile ("HiggsAnalysis120GeV.root");

 TCanvas *c1 = new TCanvas("c1", "c1",430, 10, 600, 600);
 c1->SetBorderMode(0);
 h_mass_2gammaSelEB->SetLineWidth(2);
 h_mass_2gammaSelEB->Draw();
 h_mass_2gammaGoldenSelEB->SetLineColor(2);
 h_mass_2gammaGoldenSelEB->SetLineWidth(2);
 h_mass_2gammaGoldenSelEB->Draw("same");
 h_mass_2gamma1goodconvSelEB->SetLineColor(3);
 h_mass_2gamma1goodconvSelEB->SetLineWidth(2);
 h_mass_2gamma1goodconvSelEB->Draw("same");
 h_mass_2gamma2convSelEB->SetLineColor(4);
 h_mass_2gamma2convSelEB->SetLineWidth(2);
 h_mass_2gamma2convSelEB->Draw("same");
 h_mass_2gamma1poorconvSelEB->SetLineColor(6);
 h_mass_2gamma1poorconvSelEB->SetLineWidth(2);
 h_mass_2gamma1poorconvSelEB->Draw("same");
 h_mass_2gammaleftoverSelEB->SetLineColor(7);
 h_mass_2gammaleftoverSelEB->SetLineWidth(2);
 h_mass_2gammaleftoverSelEB->Draw("same");



 TCanvas *c2 = new TCanvas("c2", "c2",430, 10, 600, 600);
 c2->SetBorderMode(0);
 h_mass_2gammaSelEE->SetLineWidth(2);
 h_mass_2gammaSelEE->Draw();
 h_mass_2gammaGoldenSelEE->SetLineColor(2);
 h_mass_2gammaGoldenSelEE->SetLineWidth(2);
 h_mass_2gammaGoldenSelEE->Draw("same");
 h_mass_2gamma1goodconvSelEE->SetLineColor(3);
 h_mass_2gamma1goodconvSelEE->SetLineWidth(2);
 h_mass_2gamma1goodconvSelEE->Draw("same");
 h_mass_2gamma2convSelEE->SetLineColor(4);
 h_mass_2gamma2convSelEE->SetLineWidth(2);
 h_mass_2gamma2convSelEE->Draw("same");
 h_mass_2gamma1poorconvSelEE->SetLineColor(6);
 h_mass_2gamma1poorconvSelEE->SetLineWidth(2);
 h_mass_2gamma1poorconvSelEE->Draw("same");
 h_mass_2gammaleftoverSelEE->SetLineColor(7);
 h_mass_2gammaleftoverSelEE->SetLineWidth(2);
 h_mass_2gammaleftoverSelEE->Draw("same");






 TCanvas *c3 = new TCanvas("c3", "c3",430, 10, 600, 600);
 c3->SetBorderMode(0);
 c3->Divide(3,2);
 c3->cd(1);
 h_mass_2gammaSelEB->SetLineWidth(2);
 h_mass_2gammaSelEB->Draw();
 c3->cd(2);
 h_mass_2gammaGoldenSelEB->SetLineColor(2);
 h_mass_2gammaGoldenSelEB->SetLineWidth(2);
 h_mass_2gammaGoldenSelEB->Draw();
 c3->cd(3);
 h_mass_2gamma1poorconvSelEB->SetLineColor(6);
 h_mass_2gamma1poorconvSelEB->SetLineWidth(2);
 h_mass_2gamma1poorconvSelEB->Draw();
 c3->cd(4);
 h_mass_2gamma1goodconvSelEB->SetLineColor(3);
 h_mass_2gamma1goodconvSelEB->SetLineWidth(2);
 h_mass_2gamma1goodconvSelEB->Draw();
 c3->cd(5);
 h_mass_2gamma2convSelEB->SetLineColor(4);
 h_mass_2gamma2convSelEB->SetLineWidth(2);
 h_mass_2gamma2convSelEB->Draw();
 c3->cd(6);
 h_mass_2gammaleftoverSelEB->SetLineColor(7);
 h_mass_2gammaleftoverSelEB->SetLineWidth(2);
 h_mass_2gammaleftoverSelEB->Draw();



 TCanvas *c4 = new TCanvas("c4", "c4",430, 10, 600, 600);
 c4->SetBorderMode(0);
 c4->Divide(3,2);
 c4->cd(1);
 h_mass_2gammaSelEE->SetLineWidth(2);
 h_mass_2gammaSelEE->Draw();
 c4->cd(2);
 h_mass_2gammaGoldenSelEE->SetLineColor(2);
 h_mass_2gammaGoldenSelEE->SetLineWidth(2);
 h_mass_2gammaGoldenSelEE->Draw();
 c4->cd(3);
 h_mass_2gamma1poorconvSelEE->SetLineColor(6);
 h_mass_2gamma1poorconvSelEE->SetLineWidth(2);
 h_mass_2gamma1poorconvSelEE->Draw();
 c4->cd(4);
 h_mass_2gamma1goodconvSelEE->SetLineColor(3);
 h_mass_2gamma1goodconvSelEE->SetLineWidth(2);
 h_mass_2gamma1goodconvSelEE->Draw();
 c4->cd(5);
 h_mass_2gamma2convSelEE->SetLineColor(4);
 h_mass_2gamma2convSelEE->SetLineWidth(2);
 h_mass_2gamma2convSelEE->Draw();
 c4->cd(6);
 h_mass_2gammaleftoverSelEE->SetLineColor(7);
 h_mass_2gammaleftoverSelEE->SetLineWidth(2);
 h_mass_2gammaleftoverSelEE->Draw();







}
