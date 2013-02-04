#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include <ctime>
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <TProof.h>
#include <TROOT.h>
//#include "ND_Hto2Photons/TreeReaders/interface/sdaReader.h"
#include "ND_Hto2Photons/TreeReaders/interface/sdaReaderFast_5X_new.h"
#include "ND_Hto2Photons/TreeReaders/interface/Selections.h"
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

bool PhotonPreSelectionPFbased(unsigned int PhotonIndex);
string GetPhotonCat(unsigned int index);
map<TString,double> GetkFactor();
map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight);
template <class type> string makestring(type value);
void BookHistograms(HistoContainer *histoContainer);
void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  gROOT->ProcessLine(".L $CMSSW_BASE/src/ND_Hto2Photons/TreeReaders/interface/link_def.h+");
  //gROOT->ProcessLine("TProof")

  TString InputArgs(input[1]);
  bool bar = false;
  bool debug = false;
  bool mc = true;
  bool unweighted = false;
  float globalweight = 5100.8;

  int FirstFileNum = 0;

  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Debug")) debug=true;

  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;

  map<TString,double> kFactor=GetkFactor();
  map<TString,double> WeightsMap=GetWeightsMap(kFactor, globalweight);
  if (debug) cout << "argc is: " << argc << endl;
  if (argc==2) MakeFilesAndWeights(InputArgs, filesAndWeights, filelist, kFactor, WeightsMap);
  else if (argc>2) MakeFilesAndWeights(string(input[2]), InputArgs, filesAndWeights, filelist, kFactor, WeightsMap);
  
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, PromptReco, 90GeV, 95GeV, 100GeV, 105GeV, 110GeV, 115GeV, 120GeV, 130GeV, 140GeV, PJet, QCD, DY, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }

  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = "";
    if (argc>3) {
      outfilename = TString(input[3]);
      outfilename += "/";
    }
    outfilename += "Vertex_";
    outfilename += itFilePair->first;
    if (argc==2) {
      if (unweighted) outfilename.ReplaceAll(".root","_Unweighted.root");
    }

    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << outfilename << " created." << endl;

    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float fileweight = filesAndWeights[itFile].second * globalweight;
      if (itFilePair->first=="PhotonPlusJetData.root" || itFilePair->first=="Photon_Data.root" || unweighted) fileweight=1;
      
      TChain* filechain = new TChain("event");
      filechain->Add(file.c_str());

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "Reading the tree in file " << file << endl;

      map <int, double> PileUpMap;
      MakePileUpWeights(InputArgs, PileUpMap);

      Long64_t nentries = filechain->GetEntries();
      cout << "TreeRead Entries " << nentries << endl;     

      outfile->cd();

      int percent = 0;
      time_t start,now;
      time(&start);

      for ( Long64_t i = 0; i < nentries; i++ ) {

        if (bar && nentries>1000 && i % (nentries/1000) == 0) {
          time (&now);
          double elapsed = difftime(now,start);
          float fracdone = float(i)/float(nentries);
          double estimate = elapsed / fracdone;
          estimate -= elapsed;
          if (percent == 1000) percent = 999;
          ProgressBar(percent,estimate);
        }

        int ientry = filechain->LoadTree(i);
        if (ientry==0) connect_variables(filechain->GetTree());
        m_entry_number = ientry;

      }

      delete filechain;

    }

    histoContainer->Save();
    delete histoContainer;
    outfile->Write();
    outfile->Close();
    delete outfile;
    FirstFileNum+=itFilePair->second;
  }

}

bool PhotonPreSelectionPFbased(unsigned int PhotonIndex) {

  TLorentzVector PhotonP4 = *((TLorentzVector*) pho_p4()->At(PhotonIndex));
  double EtCorrEcalIso = pho_ecalsumetconedr03()[PhotonIndex] - 0.012*PhotonP4.Pt();
  double EtCorrHcalIso = pho_hcalsumetconedr03()[PhotonIndex] - 0.005*PhotonP4.Pt();
  double EtCorrTrkIso = pho_trksumpthollowconedr03()[PhotonIndex] - 0.002*PhotonP4.Pt();

  if (pho_r9()[PhotonIndex]<=0.9) {
    if (pho_isEB()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.014)) return false;
    if (pho_isEE()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.034)) return false;
    if (EtCorrEcalIso>4.0) return false;
    if (EtCorrHcalIso>4.0) return false;
    if (EtCorrTrkIso>4.0) return false;
    if ( pho_pfiso_mycharged02()->at(PhotonIndex).at(0) > 4.0 )  return false;
    if (pho_isconv()[PhotonIndex]!=1) return false;
    return true;
  } else {
    if (pho_isEB()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.082 || pho_sieie()[PhotonIndex]>0.014)) return false;
    if (pho_isEE()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.034)) return false;
    if (EtCorrEcalIso>50.0) return false;
    if (EtCorrHcalIso>50.0) return false;
    if (EtCorrTrkIso>50.0) return false;
    if ( pho_pfiso_mycharged02()->at(PhotonIndex).at(0) > 4.0 )  return false;
    if (pho_isconv()[PhotonIndex]!=1) return false;
    return true;
  }

}

string GetPhotonCat(unsigned int index) {

  string ReturnValue = "PlaceHolder";
  if (pho_r9()[index]>=0.94 && pho_isEB()[index]) ReturnValue="_cat0";
  if (pho_r9()[index]< 0.94 && pho_isEB()[index]) ReturnValue="_cat1";
  if (pho_r9()[index]>=0.94 && pho_isEE()[index]) ReturnValue="_cat2";
  if (pho_r9()[index]< 0.94 && pho_isEE()[index]) ReturnValue="_cat3";
  return ReturnValue;
  
}

map<TString,double> GetkFactor() {

  map<TString,double> kFactor;

  kFactor["DiPhotonJets"]=1.0;
  kFactor["PJet"]=1.0;
  kFactor["QCD"]=1.0;
  kFactor["WJets"]=1.0;
  kFactor["ZJets"]=1.0;
  return kFactor;
  
}

map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight) {

  map<TString,double> WeightsMap;
  WeightsMap["None"]=1/globalweight;
  WeightsMap["DiPhotonJets"]=kFactor["DiPhotonJets"]*1/1154970*75.39;
  WeightsMap["PJet15to30"]=kFactor["PJet"]*1/1970745.0*200061.7;
  WeightsMap["PJet30to50"]=kFactor["PJet"]*1/1993325.0*19931.62;
  WeightsMap["PJet50to80"]=kFactor["PJet"]*1/1995062.0*3322.309;
  WeightsMap["PJet80to120"]=kFactor["PJet"]*1/1992627.0*558.2865;
  WeightsMap["PJet120to170"]=kFactor["PJet"]*1/1989603.0*108.0068;
  WeightsMap["PJet170to300"]=kFactor["PJet"]*1/2000069.0*30.12207;
  WeightsMap["PJet300to470"]=kFactor["PJet"]*1/2000130.0*2.138632;
  WeightsMap["PJet470to800"]=kFactor["PJet"]*1/1975231.0*0.2119244;
  WeightsMap["QCD20to30"]=kFactor["QCD"]*1/35040695.0*288600000.0*0.0101;
  WeightsMap["QCD30to80"]=kFactor["QCD"]*1/33088888.0*74330000.0*0.0621;
  WeightsMap["QCD80to170"]=kFactor["QCD"]*1/34542763.0*1191000.0*0.1539;
  WeightsMap["QCD170to250"]=kFactor["QCD"]*1/31697066.0*30990.0*0.148;
  WeightsMap["QCD250to350"]=kFactor["QCD"]*1/34611322.0*4250.0*0.131;
  WeightsMap["QCD350"]=kFactor["QCD"]*1/34083442.0*810.0*0.11;
  WeightsMap["WJets"]=kFactor["WJets"]*1/18393090*30400.0;
  WeightsMap["ZJets"]=kFactor["ZJets"]*1/30421222*2950.0;
  return WeightsMap;
  
}

template <class type> string makestring(type value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

void BookHistograms(HistoContainer *histoContainer) {

  histoContainer->Add("Numvtx","Number of Primary Verticies",100,0,100);
  histoContainer->Add("NumPhoton","Number of Photons",10,0,10);
  histoContainer->Add("NumConversion","Number of Conversions",100,0,100);
  
}

void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  if (inputstring.Contains("Test") && !inputstring.Contains("_")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/PhotonPlusJetMC.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Test_Data")) {
    inputfilelist.push_back(pair<string,int> ("TestData.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/PhotonPlusJetData.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Test_NoPull")) {
    inputfilelist.push_back(pair<string,int> ("Test_NoPull.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/PhotonPlusJetMC_nopull.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("PJet50to80")) {
    inputfilelist.push_back(pair<string,int> ("PJet50to80.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-50to80_TuneZ2star_8TeV.root",kFactor["PJet"]*WeightsMap["PJet50to80"]));
  } 
  if (inputstring.Contains("Data") && !inputstring.Contains("_")) {
    inputfilelist.push_back(pair<string,int> ("DataTest.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/Photon-Run2012A-29Jun2012_0.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Higgs") && !inputstring.Contains("_")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_GluGluToHToGG_M-125_8TeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-124_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1_0.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Higgs_NoPull")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_GluGluToHToGG_M-125_8TeV_NoPull.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-125_8TeV_nopull.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Higgs_newBS")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_GluGluToHToGG_M-125_8TeV_newBS.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-125_8TeV_newBS.root",WeightsMap["None"]));
  }
    
}

void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  string outfile;
  if (infile.rfind("/")!=string::npos) outfile=infile.substr(infile.rfind("/")+1);
  else outfile=infile;
  
  inputfilelist.push_back(pair<string,int> (outfile,1));
  if (inputstring.Contains("Run2012A")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("Run2012B")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("Run2012C")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("DiPhotonJets")) inputvector.push_back(pair<string,float> (infile,kFactor["DiPhotonJets"]*WeightsMap["DiPhotonJets"]));
  if (inputstring.Contains("PJet15to30")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet15to30"]));
  if (inputstring.Contains("PJet30to50")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet30to50"]));
  if (inputstring.Contains("PJet50to80")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet50to80"]));
  if (inputstring.Contains("PJet80to120")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet80to120"]));
  if (inputstring.Contains("PJet120to170")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet120to170"]));
  if (inputstring.Contains("PJet170to300")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet170to300"]));
  if (inputstring.Contains("PJet300to470")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet300to470"]));
  if (inputstring.Contains("PJet470to800")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet470to800"]));
  if (inputstring.Contains("QCD20to30")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD20to30"]));
  if (inputstring.Contains("QCD30to80")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD30to80"]));
  if (inputstring.Contains("QCD80to170")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD80to170"]));
  if (inputstring.Contains("QCD170to250")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD170to250"]));
  if (inputstring.Contains("QCD250to350")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD250to350"]));
  if (inputstring.Contains("QCD350")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD350"]));
  if (inputstring.Contains("WJets")) inputvector.push_back(pair<string,float> (infile,kFactor["WJets"]*WeightsMap["WJets"]));
  if (inputstring.Contains("ZJets")) inputvector.push_back(pair<string,float> (infile,kFactor["ZJets"]*WeightsMap["ZJets"]));
  if (inputstring.Contains("Test") && !inputstring.Contains("TestMC")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
    
}

void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap) {

   if (inputstring.Contains("PJet")) {
     #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/PhotonPlusJet_53X.h"
   } else if (inputstring.Contains("QCD")) {
     #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/QCDEMEnriched_53X.h"
   } else {
     #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Dummy.h"
   }
   //#include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Dummy.h"
  
}

void ProgressBar(int &percent, double estimate) {
  
  cout << "\033[100m";
  cout << "\r[";
  cout << "\033[42m";
  for (int ctr = 0; ctr <= percent / 20; ++ctr) cout << "-";
  cout << "\b>";
  cout << "\033[100m";
  for (int ctr = percent / 20; ctr < 49; ++ctr) cout << " ";
  cout << "]\033[0m";
  cout << " " << (float) ++percent/10;
  if (percent%10==0) cout << ".0";
  cout << "%";
  if (int(estimate)>0) {
    cout << " Remaining " << int(estimate) / 60 << ":";
    if (int(estimate) % 60 < 10) cout << "0";
    cout << int(estimate) % 60 << " ";
  }
  cout << flush;
}
