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
#include "ND_Hto2Photons/TreeReaders/interface/sdaReaderFast.h"
#include "ND_Hto2Photons/TreeReaders/interface/Selections.h"
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

bool GoodConverion(int index);
bool leadsortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, unsigned int &leadindex);
double etaTransformation(double EtaParticle, double Zvertex);
double DeltaPhi(double Phi1, double Phi2);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
string DetectorPosition(unsigned int index);
TString MakeFileName(string filename, bool unweighted, bool onevertex);
int gettrackerconvindex(TVector3 Photonxyz, TVector3 BeamSpot);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void BookTrackerdZ(HistoContainer *histoContainer, TString histname, TString options="None");
void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ, double R, float weight, bool limit=false);
void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double R, float weight);
void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ, double Z, float weight, bool limit=false);
void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double Z, float weight);
void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  gROOT->ProcessLine(".L /data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/interface/link_def.h+");
  //gROOT->ProcessLine("TProof")
  
  bool bar = false;
  bool debug = false;
  bool mc = true;
  bool nojet = false;
  bool onevertex = false;
  bool unweighted = false;
  float globalWeight = 1206.786849;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;

  TString InputArgs(input[1]);

  //PrintWeights();
  MakeFilesAndWeights(InputArgs, filesAndWeights, filelist);

  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Debug")) debug=true;
  if (InputArgs.Contains("NoJet")) nojet=true;
  if (InputArgs.Contains("OneVertex")) onevertex=true;
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, PromptReco, 90GeV, 95GeV, 100GeV, 105GeV, 110GeV, 115GeV, 120GeV, 130GeV, 140GeV, PJet, QCD, DY, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = "Vertex_";
    outfilename += MakeFileName(itFilePair->first, unweighted, onevertex);
    
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
      if (itFilePair->first=="PhotonPlusJetData.root") weight=1;
      
      //TFile * currentFile = new TFile(file.c_str());
      //currentFile->cd();

      TChain* filechain = new TChain("event");
      filechain->Add(file.c_str());

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      Long64_t nentries = filechain->GetEntries();
      cout << "TreeRead Entries " << nentries << endl;     

      outfile->cd();

      int percent = 0;
      time_t start,now;
      time(&start);

      for ( Long64_t i = 0; i < nentries; i++ ) {

        if (bar && nentries>100 && i % (nentries/100) == 0) {
          time (&now);
          double elapsed = difftime(now,start);
          float fracdone = float(i)/float(nentries);
          double estimate = elapsed / fracdone;
          estimate -= elapsed;
          if (percent == 100) percent = 99;
          ProgressBar(percent,estimate);
        }

        int ientry = filechain->LoadTree(i);
        if (ientry==0) connect_variables(filechain->GetTree());
        m_entry_number = ientry;
        if (ientry==0 && gp_n()==0) mc=false;

        TVector3 SimVertex(0,0,0);
        vector<TVector3> BeamSpot;
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TVector3> SuperClusterxyz;
        vector<TLorentzVector> Photonp4;
        map <double, unsigned int> ptindex;
        if (pho_n()<1) continue;
        if (mc) SimVertex = *((TVector3*) gv_pos()->At(0));
        for (unsigned int j=0; j!=(unsigned int) vtx_std_xyz()->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) vtx_std_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) sc_xyz()->GetSize(); j++) SuperClusterxyz.push_back(*((TVector3*) sc_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) bs_xyz()->GetSize(); j++) BeamSpot.push_back(*((TVector3*) bs_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) pho_n(); j++) {
          Photonp4.push_back(*((TLorentzVector*) pho_p4()->At(j)));
          Photonxyz.push_back(*((TVector3*) pho_calopos()->At(j)));
          ptindex[Photonp4[j].Pt()]=j;
        }
        for (unsigned int j=0; j<(unsigned int) conv_n(); j++) {
          ConversionVertex.push_back(*((TVector3*) conv_vtx()->At(j)));
          ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum()->At(j)));
        }

        double leadpt = -1;
        unsigned int leadindex = 0;

        if (debug) cout << "Doing pt sorting." << endl;
        bool sorted = leadsortpt(ptindex, pho_n(), leadpt, leadindex);
        if (debug && !sorted) cout << "Warning: Photons not pt sorted." << endl;
        int scleadindex = pho_scind()[leadindex];

        double MaxJetPt = 0;
        if (!onevertex) {
          for (unsigned int j=0; j<(unsigned int) jet_algoPF1_n(); j++) {
            TLorentzVector JetP4 = *((TLorentzVector*) jet_algoPF1_p4()->At(j));
            double dR = (JetP4.Eta()-Photonp4[leadindex].Eta())*(JetP4.Eta()-Photonp4[leadindex].Eta())+DeltaPhi(JetP4.Phi(),Photonp4[leadindex].Phi())*DeltaPhi(JetP4.Phi(),Photonp4[leadindex].Phi());
            if (JetP4.Eta()>2.4 || dR<0.4) continue;
            if (JetP4.Pt()>MaxJetPt) MaxJetPt=JetP4.Pt();
          }
        }
        histoContainer->Fill("MaxJetPt",MaxJetPt,weight);
        if (Photonp4[leadindex].Pt()<30) continue;
        if (fabs(SuperClusterxyz[scleadindex].Eta())>2.5) continue;
        if (pho_cic4cutlevel_lead()->at(leadindex).at(0)<4) continue;
        if (onevertex && PrimaryVertex.size()!=1) continue;
        if (!onevertex && !nojet && MaxJetPt<40) continue;
        int convindex = gettrackerconvindex(SuperClusterxyz[scleadindex],BeamSpot[0]);
        if (convindex==-1) continue;

        if (debug) cout << "Lead pt is: " << Photonp4[leadindex].Pt() << " and CICatagory is: " << pho_cic4cutlevel_lead()->at(leadindex).at(0) << endl;
        
        double superclusterdz = FindNewZConvLinear(ConversionVertex[convindex],SuperClusterxyz[scleadindex],BeamSpot[0]);
        double convvertexdz = FindNewdZ(ConversionVertex[convindex], ConversionRefittedPairMomentum[convindex], BeamSpot[0]);
        string PhotonDetector = DetectorPosition(leadindex);

        histoContainer->Fill("ZPV",PrimaryVertex[0].Z(),weight);
        if (mc) histoContainer->Fill("deltaZcheck",SimVertex.Z()-PrimaryVertex[0].Z(),weight);
        if (mc) histoContainer->Fill("deltaZcheckzoom",(SimVertex.Z()-PrimaryVertex[0].Z())*10*1000,weight);
        histoContainer->Fill("deltaEta",SuperClusterxyz[scleadindex].Eta() - etaTransformation(ConversionRefittedPairMomentum[convindex].Eta(),superclusterdz),weight);

        histoContainer->Fill("SuperZPV",superclusterdz,weight);
        histoContainer->Fill("ConvZPV",convvertexdz,weight);
        histoContainer->Fill("SuperdeltaZ",PrimaryVertex[0].Z()-superclusterdz,weight);
        histoContainer->Fill("ConvdeltaZ",PrimaryVertex[0].Z()-convvertexdz,weight);
        if (mc) histoContainer->Fill("SuperdeltaZsim",SimVertex.Z()-superclusterdz,weight);
        if (mc) histoContainer->Fill("ConvdeltaZsim",SimVertex.Z()-convvertexdz,weight);

        histoContainer->Fill("SuperZPV",PhotonDetector,superclusterdz,weight);
        histoContainer->Fill("ConvZPV",PhotonDetector,convvertexdz,weight);
        histoContainer->Fill("SuperdeltaZ",PhotonDetector,PrimaryVertex[0].Z()-superclusterdz,weight);
        histoContainer->Fill("ConvdeltaZ",PhotonDetector,PrimaryVertex[0].Z()-convvertexdz,weight);
        if (mc) histoContainer->Fill("SuperdeltaZsim",PhotonDetector,SimVertex.Z()-superclusterdz,weight);
        if (mc) histoContainer->Fill("ConvdeltaZsim",PhotonDetector,SimVertex.Z()-convvertexdz,weight);

        if (PhotonDetector=="Barrel") {
          FilldZTrackerBarrel(histoContainer, "SuperdZ", PrimaryVertex[0].Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight, true);
          FilldZTrackerBarrel(histoContainer, "SuperdZEff", PrimaryVertex[0].Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight);
          FilldZTrackerBarrel(histoContainer, "ConvdZ", PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight, true);
          FilldZTrackerBarrel(histoContainer, "ConvdZEff", PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
          FilldZTrackerBarrel(histoContainer, "ConvdZCompairison", PrimaryVertex[0].Z()-superclusterdz, PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
          if (mc) FilldZTrackerBarrel(histoContainer, "SuperdZsim", SimVertex.Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight);
          if (mc) FilldZTrackerBarrel(histoContainer, "ConvdZsim", SimVertex.Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
        }
        if (PhotonDetector=="Endcap") {
          FilldZTrackerEndcap(histoContainer, "SuperdZ", PrimaryVertex[0].Z()-superclusterdz, abs(ConversionVertex[convindex].Z()), weight, true);
          FilldZTrackerEndcap(histoContainer, "SuperdZEff", PrimaryVertex[0].Z()-superclusterdz, abs(ConversionVertex[convindex].Z()), weight);
          FilldZTrackerEndcap(histoContainer, "ConvdZ", PrimaryVertex[0].Z()-convvertexdz, abs(ConversionVertex[convindex].Z()), weight, true);
          FilldZTrackerEndcap(histoContainer, "ConvdZEff", PrimaryVertex[0].Z()-convvertexdz, abs(ConversionVertex[convindex].Z()), weight);
          FilldZTrackerEndcap(histoContainer, "ConvdZCompairison", PrimaryVertex[0].Z()-superclusterdz, PrimaryVertex[0].Z()-convvertexdz, abs(ConversionVertex[convindex].Z()), weight);
          if (mc) FilldZTrackerEndcap(histoContainer, "SuperdZsim", SimVertex.Z()-superclusterdz, abs(ConversionVertex[convindex].Z()), weight);
          if (mc) FilldZTrackerEndcap(histoContainer, "ConvdZsim", SimVertex.Z()-convvertexdz, abs(ConversionVertex[convindex].Z()), weight);
        }

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

bool GoodConverion(int index) {

  bool ReturnBool = false;
  if (conv_ntracks()[index]==2 && conv_validvtx()[index]==1 && conv_chi2_probability()[index]>0.0005) ReturnBool = true;
  return ReturnBool;
  
}

bool leadsortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, unsigned int &leadindex) {

  int count=0;
  for (map<double,unsigned int>::iterator it_ptindex=ptindex.begin(); it_ptindex!=ptindex.end(); ++it_ptindex) {
    //cout << "Pt: " << it_ptindex->first << " Index: " << it_ptindex->second << endl;
    if (count==NumPhotons-1) leadindex=it_ptindex->second;
    if (count==NumPhotons-1) leadpt=it_ptindex->first;
    count++;
  }

  return true;
  
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

double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot) {

  double dz = (vtx.z()-myBeamSpot.z()) - ((vtx.x()-myBeamSpot.x())*mom.x()+(vtx.y()-myBeamSpot.y())*mom.y())/mom.Perp() * mom.z()/mom.Perp();
  return dz + myBeamSpot.z();
  
}

double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 beamSpot) {
  
  double deltaX1 = superclustervtx.X()-convvtx.X();
  double deltaY1 = superclustervtx.Y()-convvtx.Y();
  double deltaZ1 = superclustervtx.Z()-convvtx.Z();
  double R1 = sqrt(deltaX1*deltaX1+deltaY1*deltaY1);
  double tantheta = R1/deltaZ1;
  
  double deltaX2 = convvtx.X()-beamSpot.X();
  double deltaY2 = convvtx.Y()-beamSpot.Y();
  double R2 = sqrt(deltaX2*deltaX2+deltaY2*deltaY2);
  double deltaZ2 = R2/tantheta;
  double primaryvertexZ = superclustervtx.Z()-deltaZ1-deltaZ2;
  return primaryvertexZ;

}


string DetectorPosition(unsigned int index) {

  string ReturnValue="";
  if (pho_isEB()[index]) ReturnValue="Barrel";
  if (pho_isEE()[index]) ReturnValue="Endcap";
  return ReturnValue;
 
}

TString MakeFileName(string filename, bool unweighted, bool onevertex) {

  TString outfilename = "";
  
  if (unweighted) {
    outfilename = "Unweighted";
    outfilename += filename;
  } else if (onevertex) {
    outfilename = "OneVertex_";
    outfilename += filename;
  } else {
    outfilename = filename;
  }
  return outfilename;
  
}

int gettrackerconvindex(TVector3 Photonxyz, TVector3 BeamSpot) {

  int ReturnIndex = -1;
  double Mindeltaeta = 999999;
  double Mindeltaphi = 999999;
  
  for (int i=0; i<conv_n(); i++) {
    TVector3 ConversionRefittedPairMomentum = *((TVector3*) conv_refitted_momentum()->At(i));
    TVector3 ConversionVertex = *((TVector3*) conv_vtx()->At(i));
    if (ConversionRefittedPairMomentum.Pt()<1) continue;
    if (conv_chi2_probability()[i]<0.000001 || !conv_validvtx()[i] || conv_ntracks()[i]!=2) continue;

    double deltaphi = DeltaPhi(Photonxyz.Phi(),ConversionVertex.Phi());
    double zfromconv = FindNewZConvLinear(ConversionVertex,Photonxyz,BeamSpot);
    double deltaeta = abs(Photonxyz.Eta() - etaTransformation(ConversionRefittedPairMomentum.Eta(),zfromconv));

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

  histoContainer->Add("MaxJetPt","Highest Pt Jet in the Events;Max Pt Jet;Counts",100,0,300);
  
  histoContainer->Add("deltaEta","Delta Eta Between Supercluster and Conversion Vertex;#Delta#eta;Counts",101,-0.1,0.1);
  histoContainer->Add("ZPV","Z of Primary Vertex;Z (cm);Counts",101,-20,20);
  histoContainer->Add("deltaZcheck","#DeltaZ of Primary Vertex and the Sim Vertex;#DeltaZ (cm);Counts",101,-2,2);
  histoContainer->Add("deltaZcheckzoom","#DeltaZ of Primary Vertex and the Sim Vertex;#DeltaZ (#mum);Counts",201,-100,100);

  histoContainer->Add("SuperZPV","Z of Primary Vertex from SuperCluster Conversion;Z (cm);Counts",101,-20,20);
  histoContainer->Add("SuperdeltaZ","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion;#DeltaZ (cm);Counts",101,-2,2);
  histoContainer->Add("SuperdeltaZsim","#DeltaZ of Sim Vertex and Primary Vertex from SuperCluster Conversion;#DeltaZ (cm);Counts",101,-2,2);
  BookBarrelAndEndcap(histoContainer,"SuperZPV","Z of Primary Vertex from SuperCluster Conversion: region;Z (cm);Counts",101,-20,20);
  BookBarrelAndEndcap(histoContainer,"SuperdeltaZ","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion: region;#DeltaZ (cm);Counts",101,-2,2);
  BookBarrelAndEndcap(histoContainer,"SuperdeltaZsim","#DeltaZ of Sim Vertex and Primary Vertex from SuperCluster Conversion: region;#DeltaZ (cm);Counts",101,-2,2);
  BookTrackerdZ(histoContainer,"SuperdZ");
  BookTrackerdZ(histoContainer,"SuperdZEff","efficiency");
  BookTrackerdZ(histoContainer,"SuperdZsim");

  histoContainer->Add("ConvZPV","Z of Primary Vertex from Conversion Vertex;Z (cm);Counts",101,-20,20);
  histoContainer->Add("ConvdeltaZ","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex;#DeltaZ (cm);Counts",101,-2,2);
  histoContainer->Add("ConvdeltaZsim","#DeltaZ of Sim Vertex and Primary Vertex from Conversion Vertex;#DeltaZ (cm);Counts",101,-2,2);
  BookBarrelAndEndcap(histoContainer,"ConvZPV","Z of Primary Vertex from Conversion Vertex: region;Z (cm);Counts",101,-20,20);
  BookBarrelAndEndcap(histoContainer,"ConvdeltaZ","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex: region;#DeltaZ (cm);Counts",101,-2,2);
  BookBarrelAndEndcap(histoContainer,"ConvdeltaZsim","#DeltaZ of Sim Vertex and Primary Vertex from Conversion Vertex: region;#DeltaZ (cm);Counts",101,-2,2);
  BookTrackerdZ(histoContainer,"ConvdZ");
  BookTrackerdZ(histoContainer,"ConvdZEff","efficiency");
  BookTrackerdZ(histoContainer,"ConvdZsim");
  BookTrackerdZ(histoContainer,"ConvdZCompairison","twoD");

}

void BookTrackerdZ(HistoContainer *histoContainer, TString histname, TString options) {

  int bins = 100;
  float lowerlimit = -5;
  float upperlimit = 5;

  TString histtitle = "#DeltaZ between the Z of the Primary Vertex using method and the Primary Vertex: region;#DeltaZ (cm);";
  if (options.Contains("twoD")) histtitle =  "#DeltaZ between the Z of the Primary Vertex using method and the Primary Vertex: region;SuperCluster #DeltaZ (cm);Conversion #DeltaZ (cm)";
  histtitle.ReplaceAll("method",histname+" method");
    
  vector<pair<TString, TString> > histarray;
  histarray.push_back(pair<TString,TString> ("PixelBarrel","Pixel Barrel"));
  histarray.push_back(pair<TString,TString> ("TIB","Tracker Inner Barrel"));
  histarray.push_back(pair<TString,TString> ("TOB","Tracker Outer Barrel"));
  histarray.push_back(pair<TString,TString> ("PixelFwd","Pixel Foward"));
  histarray.push_back(pair<TString,TString> ("TID","Tracker Inner Disk"));
  histarray.push_back(pair<TString,TString> ("TEC","Tracker EndCaps"));

  for (unsigned int i=0; i<histarray.size(); i++) {
    TString histnametemp = histarray[i].first+histname;
    TString histtitletemp = histtitle;
    histtitletemp.ReplaceAll("region",histarray[i].second); 
    if (!options.Contains("efficiency")) {
      bins=50;
      if ( i==0 ) {
        lowerlimit = -0.2;
        upperlimit = 0.2;
      }
      if ( i==1 ) {
        lowerlimit = -3;
        upperlimit = 3;
      }
      if ( i==2 ) {
        lowerlimit = -5;
        upperlimit = 5;
      }
      if ( i==3 ) {
        lowerlimit = -0.5;
        upperlimit = 0.5;
      }
      if ( i==4 ) {
        lowerlimit = -3;
        upperlimit = 3;
      }
      if ( i==5 ) {
        lowerlimit = -5;
        upperlimit = 5;
      }
    }
    if (!options.Contains("twoD")) histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
    if (options.Contains("twoD")) histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit,bins,lowerlimit,upperlimit);
  }
  
}

void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ, double R, float weight, bool limit) {

  TString histnametemp = "";
  if (R<=15.0) {
    histnametemp = "PixelBarrel"+histname;
    if( fabs(deltaZ) > 0.2 && limit ) deltaZ=(0.2-0.00009)*deltaZ/fabs(deltaZ);
  }
  else if  (R>15 && R<=60.0) {
    histnametemp = "TIB"+histname;
    if( fabs(deltaZ) >3 && limit ) deltaZ=(3.-0.009)*deltaZ/fabs(deltaZ);
  }
  else if (R>60.0) {
    histnametemp = "TOB"+histname;
    if( fabs(deltaZ) >5 && limit ) deltaZ=(5.-0.009)*deltaZ/fabs(deltaZ);
  }
  histoContainer->Fill(histnametemp.Data(),deltaZ,weight);
  
}

void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ, double Z, float weight, bool limit) {

  TString histnametemp = "";
  if (Z<=50.0) {
    histnametemp = "PixelFwd"+histname;
    if( fabs(deltaZ) >0.5 && limit ) deltaZ=(0.5-0.009)*deltaZ/fabs(deltaZ);
  }
  else if (Z>50 && Z<=100.0) {
    histnametemp = "TID"+histname;
    if( fabs(deltaZ) >3 && limit ) deltaZ=(3.-0.009)*deltaZ/fabs(deltaZ);
  }
  else if (Z>100.0) {
    histnametemp = "TEC"+histname;
    if( fabs(deltaZ) >5 && limit ) deltaZ=(5.-0.009)*deltaZ/fabs(deltaZ);
  }
  histoContainer->Fill(histnametemp.Data(),deltaZ,weight);

}

void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double R, float weight) {

  TString histnametemp = "";
  if (R<=15.0) {
    histnametemp = "PixelBarrel"+histname;
    if( fabs(deltaZ1) > 0.2 ) deltaZ1=(0.2-0.00009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) > 0.2 ) deltaZ2=(0.2-0.00009)*deltaZ2/fabs(deltaZ2);
  }
  else if  (R>15 && R<=60.0) {
    histnametemp = "TIB"+histname;
    if( fabs(deltaZ1) >3 ) deltaZ1=(3.-0.009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) >3 ) deltaZ2=(3.-0.009)*deltaZ2/fabs(deltaZ2);
  }
  else if (R>60.0) {
    histnametemp = "TOB"+histname;
    if( fabs(deltaZ1) >5  ) deltaZ1=(5.-0.009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) >5  ) deltaZ2=(5.-0.009)*deltaZ2/fabs(deltaZ2);
  }
  histoContainer->Fill(histnametemp.Data(),deltaZ1,deltaZ2,weight);
  
}

void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double Z, float weight) {

  TString histnametemp = "";
  if (Z<=50.0) {
    histnametemp = "PixelFwd"+histname;
    if( fabs(deltaZ1) >0.5 ) deltaZ1=(0.5-0.009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) >0.5 ) deltaZ2=(0.5-0.009)*deltaZ2/fabs(deltaZ2);
  }
  else if (Z>50 && Z<=100.0) {
    histnametemp = "TID"+histname;
    if( fabs(deltaZ1) >3  ) deltaZ1=(3.-0.009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) >3  ) deltaZ2=(3.-0.009)*deltaZ2/fabs(deltaZ2);
  }
  else if (Z>100.0) {
    histnametemp = "TEC"+histname;
    if( fabs(deltaZ1) >5 ) deltaZ1=(5.-0.009)*deltaZ1/fabs(deltaZ1);
    if( fabs(deltaZ2) >5 ) deltaZ2=(5.-0.009)*deltaZ2/fabs(deltaZ2);
  }
  histoContainer->Fill(histnametemp.Data(),deltaZ1,deltaZ2,weight);

}

void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist) {

  if (inputstring.Contains("Data")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetData.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/PhotonPlusJetData.root",1));
  }
  if (inputstring.Contains("PJet") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetMC.root",12));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-0to15_TuneZ2.root",1/2080768.0*84200000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-15to30_TuneZ2.root",1/2046119.0*172000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-30to50_TuneZ2.root",1/2187260.0*16700));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-50to80_TuneZ2.root",1/2036704.0*2720));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-80to120_TuneZ2.root",1/2046637.0*447));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-120to170_TuneZ2.root",1/2088216.0*84.2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-170to300_TuneZ2.root",1/2069161.0*22.6));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-300to470_TuneZ2.root",1/2076880.0*1.49));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-470to800_TuneZ2.root",1/2087212.0*0.132));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-800to1400_TuneZ2.root",1/2131800.0*0.00348));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1400to1800_TuneZ2.root",1/2198160.0*0.0000127));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1800_TuneZ2.root",1/2188301.0*0.000000294));
  }
  if (inputstring.Contains("QCD") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDEMEnriched.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/QCD_Pt-20to30_EMEnriched_TuneZ2.root",1/35729669.0*236100000/0.0106));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/QCD_Pt-30to80_EMEnriched_TuneZ2.root",1/70392060.0*59440000/0.061));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/QCD_Pt-80to170_EMEnriched_TuneZ2.root",1/8150672.0*898200/0.159));
  }
  if (inputstring.Contains("Signal")) {
    inputfilelist.push_back(pair<string,int> ("120GeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-120.root",1/1206.786849));
  }
  if (inputstring.Contains("Test") && !inputstring.Contains("TestMC")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/Test.root",1/1206.786849));
  }
  if (inputstring.Contains("TestMC")) {
    inputfilelist.push_back(pair<string,int> ("TestMC.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/TestMC.root",1/1206.786849));
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
