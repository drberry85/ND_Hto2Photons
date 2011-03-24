void fitter (char * filename = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/UnweightedHiggsAnalysis130GeV.root", float startmass = 130) {
  using namespace RooFit;

  //  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleOffset(1.8,"y");
  TCanvas * c1 = new TCanvas("c1", "c1",430, 10, 800, 800);
  c1->SetBorderMode(0);
  
  TFile * hFile = new TFile(filename);
  TList * HistList = hFile->GetListOfKeys();
  
  for (Int_t i=1; i<HistList->GetSize(); ++i) {
    TString HistName(HistList->At(i)->GetName());
    //cout << "Looking at Hist: " << HistName << endl; 
    if (HistName.Contains("mass_2gammaGoldenSel")==1) {
      // Parameter used for mass in both CrystalBall or BreitWigner
      RooRealVar mass("mass","Mass_{2#gamma}", 80, 160,"GeV/c^{2}");

      //  Parameters for Crystal Ball Lineshape 
      RooRealVar m0("M_{Higgs}", "Bias", startmass, startmass-5, startmass+5);//,"GeV/c^{2}"); 
      RooRealVar sigma("#sigma_{CB}","Width", 2,0,10.0);//,"GeV/c^{2}"); 
      RooRealVar cut("#alpha","Cut", 2,0.0,20.0); 
      RooRealVar power("#gamma","Power", 2, 0, 100000.0); 
      RooCBShape CrystalBall("CrystalBall", "A  Crystal Ball Lineshape", mass, m0,sigma, cut, power);
            
      //  Parameters for Breit-Wigner Distribution
      //RooRealVar mRes("M_{Higgs}", "Higgs Resonance  Mass", startmass, 70, 200);//,"GeV/c^{2}"); 
      //RooRealVar Gamma("#Gamma", "#Gamma", 1, 0,10.0);//,"GeV/c^{2}"); 
      //RooBreitWigner BreitWigner("BreitWigner","A Breit-Wigner Distribution",mass,mRes,Gamma);

      // Convoluve the BreitWigner and Crystal Ball
      //RooFFTConvPdf Convolution("Convolution","Convolution", mass, BreitWigner, CrystalBall);

      //  Introduce a resolution model
      // Gaussian Sigma
      RooRealVar GaussianSigma("#sigma_{G}","Core Width", 1.0,0.99,1.01,"GeV/c^{2}"); 
      
      //  RooBifurGauss  res("res", "A Bifurcated Gaussian Distribution", deltam, mu,sigL,sigR);
      RooGaussian GaussianResolution("GaussianResolution", "A  Gaussian Lineshape", mass, m0,GaussianSigma);
      RooRealVar fracG("f_{G}", "Gaussian Fraction", 0.0,0.0,1.0);
      RooAddPdf ResolutionModel("ResolutionModel", "Resolution Model", GaussianResolution, CrystalBall, fracG); 
      //fracG.setConstant();
      //RooFFTConvPdf ResolutionModel("ResolutionModel","Resolution Model", mass, CrystalBall, GaussianResolution);
      
      cout << HistName << " is being fit." << endl;
      TH1F * DataHist = (TH1F*) hFile->Get(HistName);
      DataHist->UseCurrentStyle();
      TString DataHistName = HistName;
      DataHistName += "fit";
      RooDataHist data(DataHistName,DataHistName,mass,DataHist,1.0);
      //CrystalBall.fitTo(data,Range(startmass-20, startmass+10));
      ResolutionModel.fitTo(data,Range(startmass-20, startmass+10));
      RooPlot * plot = mass.frame(80,160,80);
      data.plotOn(plot,DataError(RooAbsData::None));
      plot->SetMaximum(1.1*plot->GetMaximum());
      //CrystalBall.plotOn(plot);
      //CrystalBall.paramOn(plot,Format("NEL",FixedPrecision(4)), Parameters(RooArgSet(m0, sigma, cut, power)), Layout(.55, 0.99, 0.99), ShowConstants(kFALSE));
      ResolutionModel.plotOn(plot);
      ResolutionModel.paramOn(plot,Format("NEL",FixedPrecision(4)), Parameters(RooArgSet(m0, sigma, cut, power, GaussianSigma)), Layout(.15, 0.58, 0.9), ShowConstants(kFALSE));
      TString PlotTitle(DataHist->GetTitle());
      PlotTitle += ";Mass_{2#gamma} (GeV/c^{2});Number of Weighted Events";
      plot->SetTitle(PlotTitle);
      plot->GetXaxis()->SetRangeUser(startmass-20,startmass+10);
      gStyle->SetTitleW(.54);
      plot->Draw();
      TString OutPutName = HistName + "fit.gif";
      c1->SaveAs(OutPutName);
    }
  }
  
}
