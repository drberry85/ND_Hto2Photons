{
  using namespace RooFit;
  // Parameter used for mass in both CrystalBall or BreitWigner
  RooRealVar mass("mass","Mass_{2#gamma}", 80, 160,"GeV/c^{2}");

  //  Parameters for Crystal Ball Lineshape 
  RooRealVar m0("#Delta m_{0}", "Bias", 0.0, -100.0, 100.0);//,"GeV/c^{2}"); 
  RooRealVar sigma("#sigma_{CB}","Width", 2,0,10.0);//,"GeV/c^{2}"); 
  RooRealVar cut("#alpha","Cut", 1,0,10.0); 
  RooRealVar power("#gamma","Power", 1, -10, 10.0); 
  RooCBShape CrystalBall("CrystalBall", "A  Crystal Ball Lineshape", mass, m0,sigma, cut, power);
  
  //  Parameters for Breit-Wigner Distribution
  RooRealVar mRes("M_{Higgs}", "Higgs Resonance  Mass", 110, 70, 200);//,"GeV/c^{2}"); 
  RooRealVar Gamma("#Gamma", "#Gamma", 1, 0,10.0);//,"GeV/c^{2}"); 
  RooBreitWigner BreitWigner("BreitWigner","A Breit-Wigner Distribution",mass,mRes,Gamma);

  // Convoluve the BreitWigner and Crystal Ball
  RooFFTConvPdf Convolution("Convolution","Convolution", mass, BreitWigner, CrystalBall);

  //  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile * hFile = new TFile ("../HiggsAnalysis120GeV.root");
  TCanvas *c1 = new TCanvas("c1", "c1",430, 10, 600, 600);
  TH1F * h_mass_2gammaSelEB = hFile->Get("h_mass_2gammaSelEB");
  c1->SetBorderMode(0);
  h_mass_2gammaSelEB->SetLineWidth(2);
  h_mass_2gammaSelEB->Draw();

  RooDataHist data("h_mass_2gammaSelEBfit","h_mass_2gammaSelEBfit",mass,h_mass_2gammaSelEB,1.0);
  Convolution.fitTo(data);
  RooPlot * plot = mass.frame(80,160,160);
  data->plotOn(plot);
  Convolution.plotOn(plot);
  //Convolution.paramOn(plot);
  Convolution.paramOn(plot,Format("NEU",AutoPrecision(2)), Parameters(RooArgSet(m0, sigma, cut, power, Gamma, mRes)), Layout(.67, 0.97, 0.97), ShowConstants(kFALSE));
  plot->Draw();

}
