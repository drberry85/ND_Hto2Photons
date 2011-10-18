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

bool GenMatch(TVector3 Photon);
double etaTransformation(double EtaParticle, double Zvertex);
double DeltaPhi(double Phi1, double Phi2);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
string DetectorPosition(unsigned int index);
template <class type> string makestring(type value);
string GetConversionRegion(HistoContainer *histoContainer, double Z, double R, bool isEB);
TString MakeFileName(string filename, bool unweighted, bool onevertex);
int gettrackerconvindex(TVector3 Photonxyz, TVector3 BeamSpot);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void BookEtadZ(HistoContainer *histoContainer, TString histname);
void BookPtdZ(HistoContainer *histoContainer, TString histname);
void BookTrackerdZ(HistoContainer *histoContainer, TString histname, TString options="None");
void FilldZEta(HistoContainer *histoContainer, TString histname, double PhotonEta, double deltaZ, float weight);
void FilldZPt(HistoContainer *histoContainer, TString histname, double PhotonPt, double deltaZ, float weight);
void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ, double R, float weight, bool limit=false);
void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double R, float weight);
void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ, double Z, float weight, bool limit=false);
void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ1, double deltaZ2, double Z, float weight);
void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist);
void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap);
void MakeEtWeights(TString inputstring, map<int,double> &EtMap);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  gROOT->ProcessLine(".L /data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/interface/link_def.h+");
  //gROOT->ProcessLine("TProof")

  bool background = false;
  bool bar = false;
  bool debug = false;
  bool fake = false;
  bool mc = true;
  bool nojet = false;
  bool onevertex = false;
  bool onephoton = false;
  bool prompt = false;
  bool singleleg = false;
  bool doubleleg = false;
  bool unweighted = false;
  bool usesimvertex = false;
  float globalWeight = 1078.306616393;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;

  TString InputArgs(input[1]);

  //PrintWeights();
  MakeFilesAndWeights(InputArgs, filesAndWeights, filelist);

  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Background")) background=true;
  if (InputArgs.Contains("Debug")) debug=true;
  if (InputArgs.Contains("Fake")) fake=true;
  if (InputArgs.Contains("NoJet")) nojet=true;
  if (InputArgs.Contains("OneVertex")) onevertex=true;
  if (InputArgs.Contains("OnePhoton")) onephoton=true;
  if (InputArgs.Contains("Prompt")) prompt=true;
  if (InputArgs.Contains("SingleLeg")) singleleg=true;
  if (InputArgs.Contains("DoubleLeg")) doubleleg=true; 
  if (InputArgs.Contains("SimVertex")) usesimvertex=true;
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, PromptReco, 90GeV, 95GeV, 100GeV, 105GeV, 110GeV, 115GeV, 120GeV, 130GeV, 140GeV, PJet, QCD, DY, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = "Vertex_";
    outfilename += MakeFileName(itFilePair->first, unweighted, onevertex);
    if (usesimvertex) outfilename.ReplaceAll(".root","_SimVertex.root");
    if (background) outfilename.ReplaceAll(".root","_Background.root");
    if (onephoton) outfilename.ReplaceAll(".root","_OnePhoton.root");
    if (fake) outfilename.ReplaceAll(".root","_Fake.root");
    if (prompt) outfilename.ReplaceAll(".root","_Prompt.root");
    if (singleleg) outfilename.ReplaceAll(".root","_SingleLeg.root");
    if (doubleleg) outfilename.ReplaceAll(".root","_DoubleLeg.root");
    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << "\n" << outfilename << " created." << endl;

    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float fileweight = filesAndWeights[itFile].second * globalWeight;
      if (itFilePair->first=="PhotonPlusJetData.root" || itFilePair->first=="Photon_Data.root" || unweighted) fileweight=1;
      
      TChain* filechain = new TChain("event");
      filechain->Add(file.c_str());

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      map <int, double> PileUpMap;
      MakePileUpWeights(InputArgs, PileUpMap);
      map <int, double> EtMap;
      MakeEtWeights(InputArgs, EtMap);

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
        if (ientry==0 && gp_n()==0) mc=false;
        if (debug) cout << "Looking at event: " << i << endl;
        TVector3 SimVertex = mc ? *((TVector3*) gv_pos()->At(0)) : TVector3(0,0,0);
        vector<TVector3> BeamSpot;
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TVector3> SuperClusterxyz;
        vector<TLorentzVector> SuperClusterp4;
        vector<TLorentzVector> Photonp4;
        if (pho_n()<1) continue;
        if (jet_algo2_n()<1 && !nojet) continue;
        unsigned int leadphotonindex = 0;
        for (unsigned int j=0; j!=(unsigned int) vtx_std_xyz()->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) vtx_std_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) sc_xyz()->GetSize(); j++) {
          SuperClusterxyz.push_back(*((TVector3*) sc_xyz()->At(j)));
          SuperClusterp4.push_back(*((TLorentzVector*) sc_p4()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) bs_xyz()->GetSize(); j++) BeamSpot.push_back(*((TVector3*) bs_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) pho_n(); j++) {
          Photonp4.push_back(*((TLorentzVector*) pho_p4()->At(j)));
          if (Photonp4[j].Pt()>Photonp4[leadphotonindex].Pt()) leadphotonindex=j;
          Photonxyz.push_back(*((TVector3*) pho_calopos()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) conv_n(); j++) {
          ConversionVertex.push_back(*((TVector3*) conv_vtx()->At(j)));
          if (conv_singleleg_momentum()!=NULL && conv_ntracks()[j]==1) ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_singleleg_momentum()->At(j))); else ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum()->At(j)));
        }

        float PUWeight = PileUpMap[PrimaryVertex.size()];

        for (unsigned int photonindex=0; photonindex<(unsigned int) pho_n(); photonindex++) {
          if (onephoton && photonindex!=leadphotonindex) continue;
          bool GenMatched = false;
          if (mc) GenMatched = GenMatch(Photonxyz[photonindex]);
          if (debug && GenMatched) cout << "Photon is Gen Matched" << endl;
          if (prompt && !GenMatched) continue;
          if (fake && GenMatched) continue;

          int scphotonindex = pho_scind()[photonindex];

          if (Photonp4[photonindex].Pt()<30) continue;
          if (fabs(SuperClusterxyz[scphotonindex].Eta())>2.5) continue;
          if (pho_cic4cutlevel_lead()->at(photonindex).at(0)<4 && !background) continue;
          if (pho_cic4cutlevel_lead()->at(photonindex).at(0)>=4 && background) continue;
          if (onevertex && PrimaryVertex.size()!=1) continue;
          
          float EtWeight = EtMap[(int) floor(Photonp4[photonindex].Pt()/5)];
          float weight = fileweight*PUWeight*EtWeight;
          histoContainer->Fill("PhotonPt",Photonp4[photonindex].Pt(),weight);
          histoContainer->Fill("PhotonEta",SuperClusterxyz[scphotonindex].Eta(),weight);
          if (debug) cout << "FileWeight: " << fileweight << " PileUpWeight: " << PUWeight << " EtWeight: " << EtWeight << " TotalWeight: " << weight << endl;
          
          if (!mc || unweighted) weight=1;
          histoContainer->Fill("Numvtx",PrimaryVertex.size(),weight);
          
          double MaxJetPt = 0;
          double MaxJetTrackPt = 0;
          if (!onevertex) {
            for (unsigned int j=0; j<(unsigned int) jet_algo2_n(); j++) {
              TLorentzVector JetP4 = *((TLorentzVector*) jet_algo2_p4()->At(j));
              TLorentzVector TrackSumP4(0,0,0,0);
              double dR = (JetP4.Eta()-SuperClusterxyz[scphotonindex].Eta())*(JetP4.Eta()-SuperClusterxyz[scphotonindex].Eta())+DeltaPhi(JetP4.Phi(),SuperClusterxyz[scphotonindex].Phi())*DeltaPhi(JetP4.Phi(),SuperClusterxyz[scphotonindex].Phi());
              if (JetP4.Eta()>2.4 || dR<0.4 || jet_algo2_ntk()[j]<3) continue;
              if (debug) cout << "Good jet found with: " << jet_algo2_ntk()[j] << " tracks" << endl;
              for (unsigned int k=0; k<(unsigned int) jet_algo2_ntk()[j]; k++) {
                unsigned short trackindex = jet_algo2_tkind()->at(j).at(k);
                //if (debug) cout << "Track Index is: " << trackindex << endl;
                if (trackindex>=1500) continue;
                TLorentzVector TrackP4 = *((TLorentzVector*) tk_p4()->At(trackindex));
                if (TrackP4.Pt()>1.0) TrackSumP4+=TrackP4;
              }
              if (JetP4.Pt()>MaxJetPt && TrackSumP4.Pt()>30) {
                MaxJetPt=JetP4.Pt();
                if (TrackSumP4.Pt()>MaxJetTrackPt) MaxJetTrackPt=TrackSumP4.Pt();
              }
            }
          }

          histoContainer->Fill("MaxJetPt",MaxJetPt,weight);
          histoContainer->Fill("MaxJetTrackPt",MaxJetTrackPt,weight);
          if (vtx_std_sumpt2()->size()>0) histoContainer->Fill("sumpt2",vtx_std_sumpt2()->at(0).at(0),weight);
        
          if (!onevertex && !nojet && MaxJetPt<30) continue;
          if (usesimvertex) PrimaryVertex[0]=SimVertex;
          int convindex = gettrackerconvindex(SuperClusterxyz[scphotonindex],BeamSpot[0]);
          if (convindex==-1) continue;

          double superclusterdz = FindNewZConvLinear(ConversionVertex[convindex],SuperClusterxyz[scphotonindex],BeamSpot[0]);
          double convvertexdz = FindNewdZ(ConversionVertex[convindex],ConversionRefittedPairMomentum[convindex],BeamSpot[0]);
          string PhotonDetector = DetectorPosition(photonindex);
          if (singleleg && conv_ntracks()[convindex]!=1) continue;
          if (doubleleg && conv_ntracks()[convindex]!=2) continue;
          histoContainer->Fill("PhotonPtWeight",Photonp4[photonindex].Pt(),weight);
          histoContainer->Fill("PhotonEoP",SuperClusterp4[scphotonindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          if (GenMatched) histoContainer->Fill("PhotonEoPGen",SuperClusterp4[scphotonindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          histoContainer->Fill("ZPV",PrimaryVertex[0].Z(),weight);
          if (mc) histoContainer->Fill("dZcheck",SimVertex.Z()-PrimaryVertex[0].Z(),weight);
          if (mc) histoContainer->Fill("dZcheckzoom",(SimVertex.Z()-PrimaryVertex[0].Z())*10*1000,weight);
          histoContainer->Fill("deltaEta",SuperClusterxyz[scphotonindex].Eta() - etaTransformation(ConversionRefittedPairMomentum[convindex].Eta(),superclusterdz),weight);

          histoContainer->Fill("SuperZPV",superclusterdz,weight);
          histoContainer->Fill("ConvZPV",convvertexdz,weight);
          histoContainer->Fill("SuperdZ",PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("ConvdZ",PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("SuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("ConvdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("SuperdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("ConvdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-convvertexdz,weight);
          FilldZPt(histoContainer,"SuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
          FilldZPt(histoContainer,"ConvdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
          FilldZEta(histoContainer,"SuperdZEta",abs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-superclusterdz,weight);
          FilldZEta(histoContainer,"ConvdZEta",abs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-convvertexdz,weight);

          if (mc) histoContainer->Fill("SuperdZsim",SimVertex.Z()-superclusterdz,weight);
          if (mc) histoContainer->Fill("ConvdZsim",SimVertex.Z()-convvertexdz,weight);

          histoContainer->Fill("SuperZPV",PhotonDetector,superclusterdz,weight);
          histoContainer->Fill("ConvZPV",PhotonDetector,convvertexdz,weight);
          histoContainer->Fill("SuperdZ",PhotonDetector,PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("ConvdZ",PhotonDetector,PrimaryVertex[0].Z()-convvertexdz,weight);
          if (mc) histoContainer->Fill("SuperdZsim",PhotonDetector,SimVertex.Z()-superclusterdz,weight);
          if (mc) histoContainer->Fill("ConvdZsim",PhotonDetector,SimVertex.Z()-convvertexdz,weight);

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

          string ConversionRegion = GetConversionRegion(histoContainer, ConversionVertex[convindex].Z(), ConversionVertex[convindex].Perp(), pho_isEB()[photonindex]);
          if (ConversionRegion=="PixelBarrel" || ConversionRegion=="PixelFwd" || ConversionRegion=="TID") {
            histoContainer->Fill("MixdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
            histoContainer->Fill("MixdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-convvertexdz,weight);
          }
          if (ConversionRegion=="TIB" || ConversionRegion=="TOB" || ConversionRegion=="TEC") {
            histoContainer->Fill("MixdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
            histoContainer->Fill("MixdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-superclusterdz,weight);
          }
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

bool GenMatch(TVector3 Photon) {

  for (unsigned int i=0; i<(unsigned int) gp_n(); i++) {
    if (gp_pdgid()[i]!=22 || gp_status()[i]!=1 || gp_pdgid()[gp_mother()[i]]>22) continue;
    TLorentzVector GenParticlep4 = *((TLorentzVector*) gp_p4()->At(i));
    if (GenParticlep4.Pt()<20) continue;
    double deltaeta = abs(Photon.Eta() - GenParticlep4.Eta());
    double deltaphi = DeltaPhi(Photon.Phi(),GenParticlep4.Phi());
    double DeltaR = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);

    if (DeltaR<.1) {
      return true;
    }
  }


  return false;
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

template <class type> string makestring(type value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

string GetConversionRegion(HistoContainer *histoContainer, double Z, double R, bool isEB) {

  string ReturnString = "";
  if (isEB) {
    if (R<=15.0) {
      ReturnString = "PixelBarrel";
    } else if  (R>15 && R<=60.0) {
      ReturnString = "TIB";
    } else if (R>60.0) {
      ReturnString = "TOB";
    }
  } else {
    if (Z<=50.0) {
      ReturnString = "PixelFwd";
    }  else if (Z>50 && Z<=100.0) {
      ReturnString = "TID";
    }  else if (Z>100.0) {
      ReturnString = "TEC";
    }
  }

  return ReturnString;
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
    TVector3 ConversionRefittedPairMomentum = conv_singleleg_momentum()!=NULL && conv_ntracks()[i]==1 ? *((TVector3*) conv_singleleg_momentum()->At(i)) : *((TVector3*) conv_refitted_momentum()->At(i));
    TVector3 ConversionVertex = *((TVector3*) conv_vtx()->At(i));
    if (conv_ntracks()[i]==1) 
    if (conv_ntracks()[i]!=2 && conv_ntracks()[i]!=1) continue;
    if (conv_ntracks()[i]==2 && (conv_chi2_probability()[i]<0.000001 || ConversionRefittedPairMomentum.Pt()<1)) continue;
    if (conv_ntracks()[i]==1 && ConversionRefittedPairMomentum.Pt()<1) continue;
    
    double deltaphi = DeltaPhi(Photonxyz.Phi(),ConversionVertex.Phi());
    double zfromconv = FindNewZConvLinear(ConversionVertex,Photonxyz,BeamSpot);
    //cout << "NTracks: " << conv_ntracks()[i] << " Conerion Refitted Pair Momentum Eta: " << ConversionRefittedPairMomentum.Eta() << endl;
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

  histoContainer->Add("Numvtx","Number of Primary Verticies",100,0,100);
  
  histoContainer->Add("PhotonPt","Pt of Photon;Pt (GeV);Counts",100,0,200);
  histoContainer->Add("PhotonEta","Eta of Photon;#eta;Counts",50,-2.5,2.5);
  histoContainer->Add("PhotonPtWeight","Pt of Photon;Pt (GeV);Counts",100,0,500);
  histoContainer->Add("PhotonEoP","EoP of Photon;EoP;Counts",100,0,3);
  histoContainer->Add("PhotonEoPGen","EoP of Gen Matched Photon;EoP;Counts",100,0,3);
  histoContainer->Add("MaxJetPt","Highest Pt Jet in the Events;Max Pt Jet;Counts",100,0,500);
  histoContainer->Add("MaxJetTrackPt","Highest Track Pt Jet in the Events;Max Pt Jet;Counts",100,0,500);

  histoContainer->Add("sumpt2","Sum Pt^2 of Primary Vertex;Sum Pt^2;Counts",100,0,1000);

  histoContainer->Add("deltaEta","Delta Eta Between Supercluster and Conversion Vertex;#Delta#eta;Counts",100,-0.1,0.1);
  histoContainer->Add("ZPV","Z of Primary Vertex;Z (cm);Counts",100,-20,20);
  histoContainer->Add("dZcheck","#DeltaZ of Primary Vertex and the Sim Vertex;#DeltaZ (cm);Counts",100,-5,5);
  histoContainer->Add("dZcheckzoom","#DeltaZ of Primary Vertex and the Sim Vertex;#DeltaZ (#mum);Counts",200,-100,100);

  histoContainer->Add("SuperZPV","Z of Primary Vertex from SuperCluster Conversion;Z (cm);Counts",100,-20,20);
  histoContainer->Add("SuperdZ","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion;#DeltaZ (cm);Counts",100,-2,2);
  histoContainer->Add("SuperdZsim","#DeltaZ of Sim Vertex and Primary Vertex from SuperCluster Conversion;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"SuperZPV","Z of Primary Vertex from SuperCluster Conversion: region;Z (cm);Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"SuperdZ","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion: region;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"SuperdZsim","#DeltaZ of Sim Vertex and Primary Vertex from SuperCluster Conversion: region;#DeltaZ (cm);Counts",100,-2,2);
  BookTrackerdZ(histoContainer,"SuperdZ");
  BookTrackerdZ(histoContainer,"SuperdZEff","efficiency");
  BookTrackerdZ(histoContainer,"SuperdZsim");
  BookEtadZ(histoContainer,"SuperdZEta");
  BookPtdZ(histoContainer,"SuperdZPt");
  histoContainer->Add("SuperdZPt","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("SuperdZEta","#DeltaZ of Primary Vertex and Primary Vertex from SuperCluster Conversion vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);

  histoContainer->Add("ConvZPV","Z of Primary Vertex from Conversion Vertex;Z (cm);Counts",100,-20,20);
  histoContainer->Add("ConvdZ","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex;#DeltaZ (cm);Counts",100,-2,2);
  histoContainer->Add("ConvdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Conversion Vertex;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"ConvZPV","Z of Primary Vertex from Conversion Vertex: region;Z (cm);Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"ConvdZ","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"ConvdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Conversion Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookTrackerdZ(histoContainer,"ConvdZ");
  BookTrackerdZ(histoContainer,"ConvdZEff","efficiency");
  BookTrackerdZ(histoContainer,"ConvdZsim");
  BookTrackerdZ(histoContainer,"ConvdZCompairison","twoD");
  BookEtadZ(histoContainer,"ConvdZEta");
  BookPtdZ(histoContainer,"ConvdZPt");
  histoContainer->Add("ConvdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("ConvdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);

  histoContainer->Add("MixdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("MixdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);

}

void BookEtadZ(HistoContainer *histoContainer, TString histname) {

  int bins=100;
  int lowerlimit=-5.0;
  int upperlimit=5.0;
  float etabins[] = {0,0.9,1.479,2.1,2.5};
  for (int i=0; i<4; i++) {
    TString histtitle = "#DeltaZ between the Z of the Primary Vertex and the Primary Vertex of the Conversion for Photons with an eta of eta1-eta2;#DeltaZ (cm);counts";
    TString histnametemp=histname;
    histnametemp+=makestring(etabins[i]);
    histtitle.ReplaceAll("eta1",makestring(etabins[i]));
    histtitle.ReplaceAll("eta2",makestring(etabins[i+1]));
    histoContainer->Add(histnametemp.Data(),histtitle.Data(),bins,lowerlimit,upperlimit);
  }
 
}

void BookPtdZ(HistoContainer *histoContainer, TString histname) {

  int bins=100;
  int lowerlimit=-5;
  int upperlimit=5;
  int ptbins[] = {30,40,60,100};
  for (int i=0; i<3; i++) {
    TString histtitle = "#DeltaZ between the Z of the Primary Vertex and the Primary Vertex of the Conversion for Photons in the replaceGeV Bin;#DeltaZ (cm);counts";
    TString histnametemp=histname;
    histnametemp+=makestring(ptbins[i]);
    histtitle.ReplaceAll("replace",makestring(ptbins[i]));
    histoContainer->Add(histnametemp.Data(),histtitle.Data(),bins,lowerlimit,upperlimit);
  }
  
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

void FilldZEta(HistoContainer *histoContainer, TString histname, double PhotonEta, double deltaZ, float weight) {

  float etabins[] = {0,0.9,1.479,2.1,2.5};
  for (int i=0; i<4; i++) if (PhotonEta>etabins[i] && PhotonEta<=etabins[i+1]) histname+=makestring(etabins[i]);
  if (histname=="SuperdZEta" || histname=="ConvdZEta") cout << "PhotonEta: " << PhotonEta << " HistName: " << histname << endl;
  histoContainer->Fill(histname.Data(),deltaZ,weight);

}

void FilldZPt(HistoContainer *histoContainer, TString histname, double PhotonPt, double deltaZ, float weight) {

  int ptbins[] = {30,40,60,999999};
  if (PhotonPt<ptbins[0]) return;
  for (int i=0; i<3; i++) if (PhotonPt>ptbins[i] && PhotonPt<ptbins[i+1]) histname+=makestring(ptbins[i]);
  //cout << "Photon Pt: " << PhotonPt << " HistName: " << histname << endl;
  histoContainer->Fill(histname.Data(),deltaZ,weight);
  
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

void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist) {

  float kFactor = 1.0;
  
  if (inputstring.Contains("Data") && !inputstring.Contains("Data44X")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetData.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/PhotonPlusJetData.root",1));
  }
  if (inputstring.Contains("Data44X")) {
    inputfilelist.push_back(pair<string,int> ("Photon_Data.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/Photon_Data.root",1));
  }
  if (inputstring.Contains("PJet")) {
    kFactor=1.3;
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetMC.root",12));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-0to15_TuneZ2.root",kFactor*1/2080768.0*84200000.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-15to30_TuneZ2.root",kFactor*1/2046119.0*172000.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-30to50_TuneZ2.root",kFactor*1/2187260.0*16700.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-50to80_TuneZ2.root",kFactor*1/2036704.0*2720.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-80to120_TuneZ2.root",kFactor*1/2046637.0*447.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-120to170_TuneZ2.root",kFactor*1/2088216.0*84.2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-170to300_TuneZ2.root",kFactor*1/2069161.0*22.6));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-300to470_TuneZ2.root",kFactor*1/2076880.0*1.49));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-470to800_TuneZ2.root",kFactor*1/2087212.0*0.132));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-800to1400_TuneZ2.root",kFactor*1/2131800.0*0.00348));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1400to1800_TuneZ2.root",kFactor*1/2198160.0*0.0000127));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1800_TuneZ2.root",kFactor*1/2188301.0*0.000000294));
  }
  if (inputstring.Contains("QCD")) {
    inputfilelist.push_back(pair<string,int> ("QCDEMEnriched.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/QCD_Pt-20to30_EMEnriched_TuneZ2.root",1/35729669.0*236100000.0*0.0106));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/QCD_Pt-30to80_EMEnriched_TuneZ2.root",1/70392060.0*59440000.0*0.061));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/QCD_Pt-80to170_EMEnriched_TuneZ2.root",1/8150672.0*898200.0*0.159));
  }
  if (inputstring.Contains("Diphoton")) {
    kFactor=2.0;
    inputfilelist.push_back(pair<string,int> ("Diphoton.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/DiPhotonsJets.root",kFactor*1/1150800.0*154.7));
  }
  if (inputstring.Contains("Box")) {
    kFactor=2.0;
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Box10to25.root",kFactor*1/528400.0*358.2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Box25to250.root",kFactor*1/515028.0*12.37));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Box250.root",kFactor*1/518288.0*0.000208));
  }
  if (inputstring.Contains("WJets")) {
    inputfilelist.push_back(pair<string,int> ("WJets.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/WJets.root",1/81352581.0*27770.0));
  }
  if (inputstring.Contains("ZJets")) {
    kFactor=1.15;
    inputfilelist.push_back(pair<string,int> ("ZJets.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/ZJets.root",kFactor*1/36277961.0*2475.0));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/ZJets_1.root",kFactor*1/36277961.0*2475.0));
  }
  if (inputstring.Contains("Signal")) {
    inputfilelist.push_back(pair<string,int> ("120GeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-120.root",1/1078.306616393));
  }
  if (inputstring.Contains("130GeV")) {
    inputfilelist.push_back(pair<string,int> ("130GeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/130GeVHiggs.root",1/1078.306616393));
  }
  if (inputstring.Contains("Test") && !inputstring.Contains("TestMC")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/Test.root",1/1078.306616393));
  }
  if (inputstring.Contains("TestMC")) {
    inputfilelist.push_back(pair<string,int> ("TestMC.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-50to80_TuneZ2.root",1/1078.306616393));
  }

}

void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap) {

  if (inputstring.Contains("PJet")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/PhotonPlusJet.h"
  } else if (inputstring.Contains("QCD")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/QCDEMEnriched.h"
  } else if (inputstring.Contains("Diphoton")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Diphoton.h"
  } else if (inputstring.Contains("WJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/WJets.h"
  } else if (inputstring.Contains("ZJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZJets.h"
  } else if (inputstring.Contains("Signal")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/120GeV.h"
  } else {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Dummy.h"
  }

}

void MakeEtWeights(TString inputstring, map<int,double> &EtMap) {

  if (inputstring.Contains("PJet")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/PhotonPlusJet.h"
  } else if (inputstring.Contains("QCD")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/QCDEMEnriched.h"
  } else if (inputstring.Contains("Diphoton")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/Diphoton.h"
  } else if (inputstring.Contains("WJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/WJets.h"
  } else if (inputstring.Contains("ZJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/ZJets.h"
  } else {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/Dummy.h"
  }

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
