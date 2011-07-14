#include "TCanvas.h"
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

#include "ND_Hto2Photons/TreeReaders/interface/sdaReader.h"
//#include "ND_Hto2Photons/TreeReaders/interface/sdaReaderfast.h"
#include "ND_Hto2Photons/TreeReaders/interface/Selections.h"
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

bool GoodConverion(sdaReader *currentTree, int index);
bool sortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, double &subleadpt, unsigned int &leadindex, unsigned int &subleadindex);
double etaTransformation(double EtaParticle, double Zvertex);
double DeltaPhi(double Phi1, double Phi2);
string DetectorPosition(sdaReader *currentTree, unsigned int index);
TString MakeFileName(string filename, bool unweighted, bool dataweight, double RCut);
unsigned int DoGenMatching(sdaReader *currentTree, TVector3 Photon);
int gettrackerconvindex(sdaReader *currentTree, TVector3 Photonxyz);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool dataweight = false;
  bool debug = false;
  bool unweighted = false;
  double RCut = 999999;
  float globalWeight = 1000;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;
    
  TString InputArgs(input[1]);

  //PrintWeights();
  MakeFilesAndWeights(InputArgs, filesAndWeights, filelist, data);

  if (InputArgs.Contains("Dataweight")) {
    globalWeight = 94.4;
    dataweight="Dataweight";
  }
  if (InputArgs.Contains("PromptRecoweight")) {
    globalWeight = 228.6;
    dataweight="PromptRecoweight";
  }
  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Debug")) debug=true;
  if (InputArgs.Contains("DiPhoton")) debug=true;
  if (InputArgs.Contains("RCut")) RCut=40;
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, PromptReco, 90GeV, 95GeV, 100GeV, 105GeV, 110GeV, 115GeV, 120GeV, 130GeV, 140GeV, PJet, QCD, DY, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = "Conversion_";
    outfilename += MakeFileName(itFilePair->first, unweighted, dataweight, RCut);
    
    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << "\n" << outfilename << " created." << endl;

    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float weight = filesAndWeights[itFile].second * globalWeight;

      if (unweighted) weight=1;
      if (itFilePair->first=="Data.root" || itFilePair->first=="YousiData.root") weight=1;
      
      TFile * currentFile = new TFile(file.c_str());
      currentFile->cd();

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      sdaReader currentTree(currentFile);
      //sdaReaderfast currentFilefast(currentFilefast);
      
      Long64_t nentries = currentTree.fChain->GetEntries();
      cout << "TreeRead Entries " << nentries << endl;     

      outfile->cd();

      int percent = 0;
      time_t start,now;
      time(&start);
      
      for ( Long64_t i = 0; i < nentries; i++ ) {

        if (i % (nentries/100) == 0 && bar) {
          time (&now);
          double elapsed = difftime(now,start);
          float fracdone = float(i)/float(nentries);
          double estimate = elapsed / fracdone;
          estimate -= elapsed;
          if (percent == 100) percent = 99;
          ProgressBar(percent,estimate);
        }

        currentTree.GetEntry(i);

        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionRefittedPairMomentum;

        vector<TVector3> Photonxyz;
        vector<TLorentzVector> Photonp4;
        vector<TLorentzVector> SuperClusterp4;
        map <double, unsigned int> ptindex;

        if (currentTree.pho_n<2) continue;

        
        for (unsigned int j=0; j<(unsigned int) currentTree.sc_p4->GetSize(); j++) SuperClusterp4.push_back(*((TLorentzVector*) currentTree.sc_p4->At(j)));
        for (unsigned int j=0; j<(unsigned int) currentTree.pho_n; j++) {
          Photonxyz.push_back(*((TVector3*) currentTree.pho_calopos->At(j)));
          Photonp4.push_back(*((TLorentzVector*) currentTree.pho_p4->At(j)));
          ptindex[Photonp4[j].Pt()]=j;
        }
        for (unsigned int j=0; j<(unsigned int) currentTree.conv_n; j++) {
          ConversionVertex.push_back(*((TVector3*) currentTree.conv_vtx->At(j)));
          ConversionRefittedPairMomentum.push_back(*((TVector3*) currentTree.conv_refitted_momentum->At(j)));
        }

        double leadpt = -1;
        double subleadpt = -1;
        unsigned int leadindex = 0;
        unsigned int subleadindex = 0;

        if (debug) cout << "Doing pt sorting." << endl;
        bool sorted = sortpt(ptindex, currentTree.pho_n, leadpt, subleadpt, leadindex, subleadindex);
        if (debug && !sorted) cout << "Warning: Photons not pt sorted." << endl;
        
        if (Photonp4[leadindex].Pt()<38.33) continue; // leading photon
        if (Photonp4[subleadindex].Pt()<28.75) continue; // subleading photon
        if (fabs(Photonxyz[leadindex].Eta())>2.5 || fabs(Photonxyz[subleadindex].Eta())>2.5) continue;
        if (currentTree.pho_cic4cutlevel_lead->at(leadindex)<4 && currentTree.pho_cic4cutlevel_sublead->at(subleadindex)<4) continue;
        
        for ( unsigned int j = 0; j < (unsigned int) currentTree.pho_n; j++) {

          if (j!=leadindex && j!=subleadindex) continue;
          string PhotonDetector = DetectorPosition(&currentTree, j);
          TLorentzVector MatchedSuperCluster = SuperClusterp4[currentTree.pho_scind[j]];
          int nconv = 0;
          
          for ( unsigned int k = 0; k < (unsigned int) currentTree.conv_n; k++) {

            double deltaeta = Photonxyz[j].Eta() - etaTransformation(ConversionRefittedPairMomentum[k].Eta(),currentTree.conv_zofprimvtxfromtrks[k]);
            double deltaphi = DeltaPhi(Photonxyz[j].Phi(),ConversionRefittedPairMomentum[k].Phi());
            if (abs(deltaphi)<0.1 && abs(deltaeta)<0.1) nconv++;
            double EoP = MatchedSuperCluster.E()/ConversionRefittedPairMomentum[k].Mag();

            histoContainer->Fill("convr",PhotonDetector,ConversionVertex[k].Perp(),weight);
            histoContainer->Fill("convdeta",PhotonDetector,deltaeta,weight);
            histoContainer->Fill("convdetauncorrected",PhotonDetector,Photonxyz[j].Eta()-ConversionRefittedPairMomentum[k].Eta(),weight);
            histoContainer->Fill("convdphi",PhotonDetector,deltaphi,weight);
            histoContainer->Fill("convdr",PhotonDetector,sqrt(deltaeta*deltaeta+deltaphi*deltaphi),weight);
            histoContainer->Fill("conveop",PhotonDetector,EoP,weight);

            histoContainer->Fill("conveta",ConversionRefittedPairMomentum[k].Eta(),weight);
            histoContainer->Fill("convetacorrected",etaTransformation(ConversionRefittedPairMomentum[k].Eta(),currentTree.conv_zofprimvtxfromtrks[k]),weight);
            histoContainer->Fill("convphi",PhotonDetector,ConversionRefittedPairMomentum[k].Phi(),weight);
            histoContainer->Fill("convpt",PhotonDetector,ConversionRefittedPairMomentum[k].Pt(),weight);
            
          }
          histoContainer->Fill("nconv",PhotonDetector,nconv,weight);
          
          //if (Photonp4[j].Pt()<30) continue;
          //if (abs(Photonxyz[j].Eta())>2.5) continue;
          //if (currentTree.pho_isEBEEGap[j]) continue;
          
          histoContainer->Fill("photoneta",Photonp4[j].Eta(),weight);
          histoContainer->Fill("photonr9",PhotonDetector,currentTree.pho_r9[j],weight);

          if (!data) {
            unsigned int GenMatchedIndex = 999999;
            GenMatchedIndex = DoGenMatching(&currentTree, Photonxyz[j]);
            if (GenMatchedIndex == 999999) continue;
            histoContainer->Fill("matchedphotoneta",Photonp4[j].Eta(),weight);
            histoContainer->Fill("matchedphotonr9",PhotonDetector,currentTree.pho_r9[j],weight);
          }
          
          int convindex = gettrackerconvindex(&currentTree,Photonxyz[j]);
          //int convindex = currentTree.pho_matchingConv->at(j);
          if (convindex==-1 || currentTree.conv_validvtx[convindex]==0) continue;
          histoContainer->Fill("convphotoneta",Photonp4[j].Eta(),weight);
          histoContainer->Fill("convphotonr9",PhotonDetector,currentTree.pho_r9[j],weight);

          bool goodconversion = GoodConverion(&currentTree, convindex);
          if (!goodconversion) continue;
          histoContainer->Fill("goodconvphotoneta",Photonp4[j].Eta(),weight);
          histoContainer->Fill("goodconvphotonr9",PhotonDetector,currentTree.pho_r9[j],weight);

          double deltaeta = Photonxyz[j].Eta() - etaTransformation(ConversionRefittedPairMomentum[convindex].Eta(),currentTree.conv_zofprimvtxfromtrks[convindex]);
          double deltaphi = DeltaPhi(Photonxyz[j].Phi(),ConversionRefittedPairMomentum[convindex].Phi());
          double EoP = MatchedSuperCluster.E()/ConversionRefittedPairMomentum[convindex].Mag();
          
          histoContainer->Fill("goodconvr",PhotonDetector,ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("goodconvchi2",PhotonDetector,currentTree.conv_chi2_probability[convindex],weight);
          histoContainer->Fill("goodconvdistmin",PhotonDetector,currentTree.conv_distofminapproach[convindex],weight);
          histoContainer->Fill("goodconveop",PhotonDetector,EoP,weight);
          histoContainer->Fill("goodconvdeta",PhotonDetector,deltaeta,weight);
          histoContainer->Fill("goodconvdphi",PhotonDetector,deltaphi,weight);
          histoContainer->Fill("goodconveta",Photonxyz[j].Eta(),weight);
          histoContainer->Fill("goodconvphi",Photonxyz[j].Phi(),weight);
        }
        
      }
      currentFile->Close();
      delete currentFile;

    }
    
    histoContainer->Save();
    delete histoContainer;
    outfile->Write();
    outfile->Close();
    delete outfile;
    FirstFileNum+=itFilePair->second;
  }
}

bool GoodConverion(sdaReader *currentTree, int index) {

  bool ReturnBool = false;
  if (currentTree->conv_ntracks[index]==2 && currentTree->conv_validvtx[index]==1 && currentTree->conv_chi2_probability[index]>0.0005) ReturnBool = true;
  return ReturnBool;
  
}

bool sortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, double &subleadpt, unsigned int &leadindex, unsigned int &subleadindex) {

  int count=0;
  for (map<double,unsigned int>::iterator it_ptindex=ptindex.begin(); it_ptindex!=ptindex.end(); ++it_ptindex) {
    //cout << "Pt: " << it_ptindex->first << " Index: " << it_ptindex->second << endl;
    if (count==NumPhotons-1) leadindex=it_ptindex->second;
    if (count==NumPhotons-1) leadpt=it_ptindex->first;
    if (count==NumPhotons-2) subleadindex=it_ptindex->second;
    if (count==NumPhotons-2) subleadpt=it_ptindex->first;
    count++;
  }

  if (leadpt<=subleadpt) {
    return false;
  } else {
    return true;
  }
  
}

double etaTransformation(double EtaParticle, double Zvertex)  {

  //---Definitions
  const float PI    = 3.1415927;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction
  float Theta = 0.0; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+PI;
  double ETA = -log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+PI ;
      ETA = - log(tan(0.5*Theta));                    
    } 
  //---Return the result
  return ETA;
  //---end
}

double DeltaPhi(double Phi1, double Phi2) {
  double Pi = 3.14159265;
  double deltaphi = Phi1 - Phi2;
  if (deltaphi>Pi) deltaphi = 2*Pi-deltaphi;
  if (deltaphi<-Pi) deltaphi = 2*Pi+deltaphi;
  return deltaphi;
}

string DetectorPosition(sdaReader *currentTree, unsigned int index) {

  string ReturnValue="";
  if (currentTree->pho_isEB[index]) ReturnValue="Barrel";
  if (currentTree->pho_isEE[index]) ReturnValue="Endcap";
  return ReturnValue;
 
}

TString MakeFileName(string filename, bool unweighted, bool dataweight, double RCut) {

  TString outfilename = "";
  
  if (unweighted) {
    outfilename = "Unweighted";
    outfilename += filename;
  } else if (dataweight) {
    outfilename = "Dataweight";
    outfilename += filename;
  } else {
    outfilename = filename;
  }
  if (RCut!=999999) {
    
    outfilename.ReplaceAll(".root","RCut.root");
  }
  return outfilename;
  
}

unsigned int DoGenMatching(sdaReader *currentTree, TVector3 Photon) {
  unsigned int ReturnValue = 999999;

  for (unsigned int i=0; i< (unsigned int) currentTree->gp_n; i++) {
    if (currentTree->gp_pdgid[i]!=22) continue;
    TLorentzVector GenParticlep4 = (*((TLorentzVector*) currentTree->gp_p4->At(i)));
    if (GenParticlep4.Pt()<20) continue;
    double deltaeta = abs(Photon.Eta() - GenParticlep4.Eta());
    double deltaphi = DeltaPhi(Photon.Phi(),GenParticlep4.Phi());
    double DeltaR = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);

    if (DeltaR<.1) {
      ReturnValue=i;
      break;
    }
  }

  return ReturnValue;

}

int gettrackerconvindex(sdaReader *currentTree, TVector3 Photonxyz) {

  int ReturnIndex = -1;
  double Mindeltaeta = 999999;
  double Mindeltaphi = 999999;
  
  for (int i=0; i<currentTree->conv_n; i++) {
    TVector3 ConversionRefittedPairMomentum = *((TVector3*) currentTree->conv_refitted_momentum->At(i));
    if (ConversionRefittedPairMomentum.Pt()<1) continue;
    double deltaphi = DeltaPhi(Photonxyz.Phi(),ConversionRefittedPairMomentum.Phi());

    double ConvEta = etaTransformation(ConversionRefittedPairMomentum.Eta(),currentTree->conv_zofprimvtxfromtrks[i]);
    double deltaeta = abs(Photonxyz.Eta()-ConvEta);

    if (abs(deltaeta)<abs(Mindeltaeta) && abs(deltaphi)<abs(Mindeltaphi)) {
      Mindeltaphi=abs(deltaphi);
      Mindeltaeta=abs(deltaeta);
      ReturnIndex = i;
    }
  }

  if (abs(Mindeltaeta)<0.1 && abs(Mindeltaphi)<0.1) {
    return ReturnIndex;
  } else {
    return -1;
  }

}

void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {
  TString Regions[2] = {"Barrel","Endcap"};
  for (int i=0; i<2; i++) {
    TString histnametemp = histname;
    TString histtitletemp = histtitle;
    histnametemp += Regions[i];
    histtitletemp.ReplaceAll("region",Regions[i]);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
  }
}

void BookHistograms(HistoContainer *histoContainer) {

  BookBarrelAndEndcap(histoContainer,"convr","R of conversion; R (cm): region; Counts",100,0,100);
  BookBarrelAndEndcap(histoContainer,"nconv","Number of Matched Conversions; Number of Conversions; Counts",20,0,20);
  BookBarrelAndEndcap(histoContainer,"conveop","E over P of conversion: region; E over P; Counts",100,0,3);
  BookBarrelAndEndcap(histoContainer,"convdphi","#Delta#phi between the refitted tracks and supercluster; #Delta#phi; Counts",100,-0.5,0.5);
  BookBarrelAndEndcap(histoContainer,"convdeta","#Delta#eta between the refitted tracks and supercluster; #Delta#eta; Counts",100,-0.5,0.5);
  BookBarrelAndEndcap(histoContainer,"convdetauncorrected","#Delta#eta uncorrected between the refitted tracks and supercluster; #Delta#eta; Counts",100,-0.5,0.5);
  BookBarrelAndEndcap(histoContainer,"convdr","#DeltaR between the refitted tracks and supercluster; #DeltaR; Counts",100,0,1);
  BookBarrelAndEndcap(histoContainer,"convphi","#phi of refitted tracks; #phi; Counts",100,-3.2,3.2);
  BookBarrelAndEndcap(histoContainer,"convpt","Pt of refitted tracks; Pt (GeV); Counts",300,0,300);

  histoContainer->Add("conveta","#eta of refitted tracks; #eta; Counts",100,-3,3);
  histoContainer->Add("convetacorrected","Corrected #eta of refitted tracks; #eta; Counts",100,-3,3);
  histoContainer->Add("photoneta","#eta of photons; #eta; Counts",100,-3,3);
  histoContainer->Add("matchedphotoneta","#eta of gen matched photons; #eta; Counts",100,-3,3);
  histoContainer->Add("convphotoneta","#eta of photons with valid conversion vertexes; #eta; Counts",100,-3,3);
  histoContainer->Add("goodconvphotoneta","#eta of photons with good conversions; #eta; Counts",100,-3,3);

  BookBarrelAndEndcap(histoContainer,"photonr9","R9 of photons; R9; Counts",100,0,1);
  BookBarrelAndEndcap(histoContainer,"matchedphotonr9","R9 of matched photons; R9; Counts",100,0,1);
  BookBarrelAndEndcap(histoContainer,"convphotonr9","R9 of photons with valid conversion vertexes; R9; Counts",100,0,1);
  BookBarrelAndEndcap(histoContainer,"goodconvphotonr9","r9 of photons with good conversions; #eta; Counts",100,0,1);

  BookBarrelAndEndcap(histoContainer,"goodconvr","R of conversion; R (cm): region;",100,0,100);
  BookBarrelAndEndcap(histoContainer,"goodconvchi2","#Chi^2 probability of vertex; Pt (GeV);",100,0,1);
  BookBarrelAndEndcap(histoContainer,"goodconvdistmin","Distantance of minimum approach of conversion tracks; Pt (GeV);",100,0,1);
  BookBarrelAndEndcap(histoContainer,"goodconveop","E over P of conversion: region; E over P;",100,0,3);
  BookBarrelAndEndcap(histoContainer,"goodconvdeta","#Delta#eta between the refitted tracks and supercluster; #Delta#eta;",100,-0.2,0.2);
  BookBarrelAndEndcap(histoContainer,"goodconvdphi","#Delta#phi between the refitted tracks and supercluster; #Delta#phi;",100,-0.2,0.2);
  histoContainer->Add("goodconveta","#eta of matched good conversion; #eta;",100,-3,3);
  histoContainer->Add("goodconvphi","#phi of matched good conversion; #phi;",100,-3.2,3.2);

}

void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData) {

  float BranchingFraction = 0;

  if (inputstring.Contains("Data") && !inputstring.Contains("Dataweight")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/paoloReRecoWithLaserV5ESFixed_v1.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/paoloReRecoWithLaserV5_v2.root",1));
    isData=true;
  }
  if (inputstring.Contains("PromptReco") && !inputstring.Contains("PromptRecoTrial") && !inputstring.Contains("PromptRecoweight")) {
    inputfilelist.push_back(pair<string,int> ("PromptReco.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_160404_161474.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_161475_163754.root",1));
    isData=true;
  }
  if (inputstring.Contains("PromptRecoTrial")) {
    inputfilelist.push_back(pair<string,int> ("PromptReco.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_160404_161474.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_161475_163754.root",1));
    isData=true;
  }
  if (inputstring.Contains("90GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00123;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis90GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M90_2011PU.root",29.47*BranchingFraction/109995));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M90_2011PU.root",1.710*BranchingFraction/108813));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M90_2011PU.root",(1.640+0.8597+0.2162)*BranchingFraction/110000));
  }
  if (inputstring.Contains("95GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00140;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis95GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M95_2011PU.root",26.58*BranchingFraction/109988));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M95_2011PU.root",1.628*BranchingFraction/109579));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M95_2011PU.root",(1.392+0.7348+0.1880)*BranchingFraction/110000));
  }
  if (inputstring.Contains("100GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00159;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis100GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M100_2011PU.root",24.02*BranchingFraction/109989));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M100_2011PU.root",1.546*BranchingFraction/109826));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M100_2011PU.root",(1.186+0.6313+0.1638)*BranchingFraction/110000));
  }
  if (inputstring.Contains("105GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00178;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis105GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M105_2011PU.root",21.78*BranchingFraction/109470));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M105_2011PU.root",1.472*BranchingFraction/109835));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M105_2011PU.root",(1.018+0.5449+0.1433 )*BranchingFraction/110000));
  }
  if (inputstring.Contains("110GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00197;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis110GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M110_2011PU.root",19.84*BranchingFraction/109994));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M110_2011PU.root",1.398*BranchingFraction/108563));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M110_2011PU.root",(0.8754+0.4721+0.1257)*BranchingFraction/110000));
  }
  if (inputstring.Contains("115GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00213;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis115GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M115_2011PU.root",18.13*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M115_2011PU.root",1.332*BranchingFraction/109834));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M115_2011PU.root",(0.7546+0.4107+0.1106)*BranchingFraction/110000));
  }
  if (inputstring.Contains("120GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00225;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M120_2011PU.root",16.63*BranchingFraction/109992));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M120_2011PU.root",1.269*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M120_2011PU.root",(0.6561+0.3598+0.09756)*BranchingFraction/110000));
  }
  if (inputstring.Contains("130GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00226;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M130_2011PU.root",14.12*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M130_2011PU.root",1.154*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M130_2011PU.root",(0.5008+0.2778+0.07658)*BranchingFraction/110000));
  }
  if (inputstring.Contains("140GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.00194;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis140GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M140_2011PU.root",12.13*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/VBF_M140_2011PU.root",1.052*BranchingFraction/109842));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/WHZHTTH_M140_2011PU.root",(0.3857+0.2172+0.06072 )*BranchingFraction/110000));
  }
  if (inputstring.Contains("PJet") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("PhotonJetEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GJet_Pt20_2011PU.root",77100/(1182075/0.0064)));
  }
  if (inputstring.Contains("QCD") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDDoubleEMEnriched.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/QCD_Pt30to40_2011PU.root",41800000/(3550408/0.00023)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/QCD_Pt40_2011PU.root",18700000/(21276029/0.00216)));
  }
  if (inputstring.Contains("Born") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Born.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt10to25_2011PU.root",236.4/522865));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt25to250_2011PU.root",22.37/537445));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt250toinf_2011PU.root",0.008072/546355));
  }
  if (inputstring.Contains("Box") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt10to25_2011PU.root",358.2/797975));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt25to250_2011PU.root",12.37/777725));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt250toinf_2011PU.root",0.000208/789470));
  }
  if (inputstring.Contains("DY") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("DrellYan.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DYToEE_M10To20_2011PU.root",2659.0/1933000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DYToEE_powheg_M20_2011PU.root",1614.0/1998990));
  }
  if (inputstring.Contains("Test")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M120_2011PUTest.root",1));
  }

}

void ProgressBar(int &percent, double estimate) {
  cout << "\033[100m";
  cout << "\r[";
  cout << "\033[42m";
  for (int ctr = 0; ctr <= percent / 2; ++ctr) cout << "-";
  cout << "\b>";
  cout << "\033[100m";
  for (int ctr = percent / 2; ctr < 49; ++ctr) cout << " ";
  cout << "]\033[0m";
  cout << " " << ++percent << "%";
  if (int(estimate)>0) {
    cout << " Remaining " << int(estimate) / 60 << ":";
    if (int(estimate) % 60 < 10) cout << "0";
    cout << int(estimate) % 60 << " ";
  }
  cout << flush;
}
