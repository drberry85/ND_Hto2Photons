#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"
#include "TString.h"

#include <ctime>
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <sstream>
#include <string>

#include "ND_Hto2Photons/TreeReaders/interface/sdaReader.h"
#include "ND_Hto2Photons/TreeReaders/interface/Selections.h"
#include "ND_Hto2Photons/TreeReaders/interface/lumi.h"
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

unsigned int DoGenMatching(sdaReader *currentTree, TVector3 Photon);
unsigned int getconvindex(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
int gettrackerconvindex(sdaReader *currentTree, TVector3 Photonxyz);
bool DoVertexMatching(TVector3 TestVertex, TVector3 SimVertex);
bool sortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, double &subleadpt, unsigned int &leadindex, unsigned int &subleadindex);
double CosThetaStar(TLorentzVector VLead, TLorentzVector VSum);
double DeltaPhi(double Phi1, double Phi2);
double etaTransformation(double EtaParticle, double Zvertex);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
double FindRefittedZ(TVector3 ConversionVertex, TVector3 ConversionRefittedPairMomentum);
string HiggsDetectorPosition(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
string DetectorPosition(sdaReader *currentTree, unsigned int index);
TString MakeFileName(string filename, bool unweighted, string dataweight, double RCut);
TLorentzVector CalcDiffVertex(TLorentzVector p4, TVector3 xyz, float newz);
void BookConversionPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void BookFourHists(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookMassPlots(HistoContainer *histoContainer, TString histname);
void BookPhysicsPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookRCutsMassPlots(HistoContainer *histoContainer, TString histname);
void BookRCutsdZPlots(HistoContainer *histoContainer, TString histname);
void BookTrackerdZ(HistoContainer *histoContainer, TString histname);
void FillCatHists(sdaReader *currentTree, HistoContainer *histoContainer, float weight, unsigned int leadindex, unsigned int subleadindex, int leadPhoCategory, int subleadPhoCategory);
void FillPhotonHists(sdaReader *currentTree, HistoContainer *histoContainer, string photon, unsigned int index, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, TVector3 PrimaryVertex, TVector3 SimVertex, vector<TLorentzVector> Photonp4, bool data);
void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, unsigned int photonconvindex, int convindex, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop);
void FillConvHistsInDiPho(sdaReader *currentTree, HistoContainer *histoContainer, unsigned int photonconvindex, int convindex, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop);
void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double SimInvMass, double cos_thetastar);
void FillMassNewVertexHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, string vertextype, TLorentzVector VSum, double InvMass, double cos_thetastar);
void FillMassRcut(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, double R, double InvMass);
void FilldZRcut(HistoContainer *histoContainer, string selection, float weight, string iConvDetector, double R, double deltaz);
void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ, double R, float weight);
void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ, double Z, float weight);
void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData);
void ProgressBar(int &percent, double estimate);
void PrintWeights();

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool debug = false;
  bool pjet = false;
  bool unweighted = false;
  double RCut = 999999;
  float globalWeight = 1000;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;
    
  string dataweight = "";
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
  if (InputArgs.Contains("RCut")) RCut=40;
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, PromptReco, 90GeV, 95GeV, 100GeV, 105GeV, 110GeV, 115GeV, 120GeV, 130GeV, 140GeV, PJet, QCD, DY, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {
    TString outfilename = MakeFileName(itFilePair->first, unweighted, dataweight, RCut);
    
    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << "\n" << outfilename << " created." << endl;

    if (outfilename.Contains("PhotonJetEMEnriched.root")) pjet=true; else pjet=false;
    
    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float weight = filesAndWeights[itFile].second * globalWeight;

      if (unweighted) weight=1;
      if (itFilePair->first=="Data.root" || itFilePair->first=="PromptReco.root") weight=1;
      
      TFile * currentFile = new TFile(file.c_str());
      currentFile->cd();

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      sdaReader currentTree(currentFile);

      Long64_t nentries = currentTree.fChain->GetEntries();
      cout << "TreeRead Entries " << nentries << endl;     

      if (debug) cout << "Going to Outfile" << endl;
      outfile->cd();

      int percent = 0;
      time_t start,now;
      time(&start);
      if (debug) cout << "Timer Started" << endl;

      for ( Long64_t i = 0; i < nentries; i++ ) {

        if (debug) cout << "Calling Bar Function" << endl;
        if (bar && nentries>=100 && i%(nentries/100)==0) {
          time (&now);
          double elapsed = difftime(now,start);
          float fracdone = float(i)/float(nentries);
          double estimate = elapsed / fracdone;
          estimate -= elapsed;
          if (percent == 100) percent = 99;
          ProgressBar(percent, estimate);
        }

        if (debug) cout << "Getting Tree Entry: " << i << endl;
        currentTree.GetEntry(i);
        if (pjet && currentTree.process_id==18) continue;
        if (data && !lumiRunSelection(currentTree.run, currentTree.lumis)) continue;
        
        TVector3 SimVertex(0,0,0);
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionPairMomentum;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TLorentzVector> Photonp4;
        vector<TLorentzVector> SuperClusterp4;
        if (debug) cout << "Vectors Defined" << endl;
        if (debug && !data) cout << "Number Gen Particles: " << currentTree.gp_n << endl;
        if (debug && !data) cout << "Number Gen Vertexes: " << currentTree.gp_vtx->GetSize() << endl;
        if (debug && !data) cout << "Number Gen P4: " << currentTree.gp_p4->GetSize() << endl;
        //Finding Vertex of Higgs

        int HiggsIndex = -1;
        TLorentzVector Higgs_p4(0,0,0,0);
        if (!data) {
          for (int j=0; j<currentTree.gp_n; j++) {
            if (currentTree.gp_pdgid[j]==25 && currentTree.gp_status[j]==2) {
              HiggsIndex = j;
              SimVertex = (*((TVector3*) currentTree.gp_vtx->At(j)));
              Higgs_p4 = (*((TLorentzVector*) currentTree.gp_p4->At(j)));
            }
          }
        }
        if (debug) cout << "Higgs Vertexing Done" << endl;
        if (HiggsIndex != -1) {
          histoContainer->Fill("GenHiggsPt",Higgs_p4.Pt(),weight);
          histoContainer->Fill("GenHiggsEta",Higgs_p4.Eta(),weight);
        }

        for (unsigned int j=0; j!=(unsigned int) currentTree.vtx_std_xyz->GetSize(); j++) {
          PrimaryVertex.push_back(*((TVector3*) currentTree.vtx_std_xyz->At(j)));
          histoContainer->Fill("ZVertices",PrimaryVertex[j].Z(),weight);
        }        
        histoContainer->Fill("ZVertex",PrimaryVertex[0].Z(),weight);

        if (!data) {
          histoContainer->Fill("ZSimVertex",SimVertex.Z(),weight);
          histoContainer->Fill("ZdZ",PrimaryVertex[0].Z()-SimVertex.Z(),weight);
          histoContainer->Fill("ZdZZoom",PrimaryVertex[0].Z()-SimVertex.Z(),weight);
          for (unsigned int j=0; j<PrimaryVertex.size(); j++) histoContainer->Fill("ZdZAll",PrimaryVertex[j].Z()-SimVertex.Z(),weight);
        }

        if (currentTree.pho_n<1) continue;
        if (debug) cout << "One Photon Found" << endl;
        if (debug) cout << "There are " << currentTree.pho_n << " photons" << endl;
        map <double, unsigned int> ptindex;
        for (unsigned int j=0; j<(unsigned int) currentTree.pho_n; j++) {
          Photonxyz.push_back(*((TVector3*) currentTree.pho_calopos->At(j)));
          Photonp4.push_back(*((TLorentzVector*) currentTree.pho_p4->At(j)));
          ptindex[Photonp4[j].Pt()]=j;
        }
        for (unsigned int j=0; j<(unsigned int) currentTree.sc_p4->GetSize(); j++) SuperClusterp4.push_back(*((TLorentzVector*) currentTree.sc_p4->At(j)));
        for (unsigned int j=0; j<(unsigned int) currentTree.conv_n; j++) {
          ConversionVertex.push_back(*((TVector3*) currentTree.conv_vtx->At(j)));
          ConversionPairMomentum.push_back(*((TVector3*) currentTree.conv_pair_momentum->At(j)));
          ConversionRefittedPairMomentum.push_back(*((TVector3*) currentTree.conv_refitted_momentum->At(j)));
        }
        if (debug) cout << "Conversions Vectors Filled" << endl;
        
        double leadpt = -1;
        double subleadpt = -1;
        unsigned int leadindex = 0;
        unsigned int subleadindex = 0;

        if (debug) cout << "Doing pt sorting." << endl;
        bool sorted = sortpt(ptindex, currentTree.pho_n, leadpt, subleadpt, leadindex, subleadindex);
        
        if (!sorted) cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        if (debug) cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        histoContainer->Fill("NumberVertices",float(currentTree.vtx_std_xyz->GetSize()),weight);

        if (!data) {
          double genleadpt = -1;
          double gensubleadpt = -1;
          unsigned int genleadindex = 0;
          unsigned int gensubleadindex = 0;

          map <double, unsigned int> genptindex;
          //cout << "Filling Gen Tree with " << currentTree.gp_n << " entries." << endl;
          for (unsigned int j=0; j<(unsigned int) currentTree.gp_n; j++) {
            //cout << "J is: " << j << " and pdgid is: " << currentTree.gp_pdgid[j] << endl;
            if (currentTree.gp_pdgid[j]==22 && currentTree.gp_status[j]==1) {
              //cout << "Filling pt: " << ((TLorentzVector*) currentTree.gp_p4->At(j))->Pt() << " with index " << counter << endl;
              genptindex[((TLorentzVector*) currentTree.gp_p4->At(j))->Pt()]=j;
            }
          }

          if (debug) cout << "Doing Gen Sorting." << endl;
          bool gensorted = sortpt(genptindex, genptindex.size(), genleadpt, gensubleadpt, genleadindex, gensubleadindex);
          if (!gensorted) cout << "GenParticles not sorted!" << endl;

          TLorentzVector LeadGenPhotonp4 = (*((TLorentzVector*) currentTree.gp_p4->At(genleadindex)));
          TLorentzVector SubleadGenPhotonp4 = (*((TLorentzVector*) currentTree.gp_p4->At(gensubleadindex)));

          histoContainer->Fill("GenPhotonEta",LeadGenPhotonp4.Eta(),weight);
          histoContainer->Fill("GenPhotonEta",SubleadGenPhotonp4.Eta(),weight);

          if (abs(LeadGenPhotonp4.Eta())<2.5 && abs(SubleadGenPhotonp4.Eta())<2.5) {
            histoContainer->Fill("GenPhotonEtaBool",1,weight);
          } else {
            histoContainer->Fill("GenPhotonEtaBool",0,weight);
          }
        }

        //////////////// basic selection
        if (!preselection(Photonp4[leadindex].Pt(), Photonp4[subleadindex].Pt(), Photonxyz[leadindex].Eta(), Photonxyz[subleadindex].Eta(), currentTree.pho_isEBEEGap[leadindex], currentTree.pho_isEBEEGap[subleadindex])) continue;

        ////////////////////////////////////
        //unsigned int convindex = getconvindex(&currentTree,leadindex,subleadindex);
        if (debug) cout << "Matching Conversions" << endl;
        unsigned int photonconvindex = 0;
        bool validconversion = false;
        int convindex = gettrackerconvindex(&currentTree,Photonxyz[leadindex]);
        if (convindex!=-1) photonconvindex=leadindex;
        if (convindex==-1) {
          convindex = gettrackerconvindex(&currentTree,Photonxyz[subleadindex]);
          if (convindex!=-1) photonconvindex=subleadindex;
        }
        if (convindex!=-1) validconversion=true;
        unsigned int convscindex = (unsigned int) currentTree.pho_scind[photonconvindex];
        unsigned int nearvertexindex = 0;
        double scEnergy = Photonp4[photonconvindex].E();
        scEnergy = SuperClusterp4[convscindex].E();
        double EoP = 0;
        if (validconversion) EoP=scEnergy/ConversionRefittedPairMomentum[convindex].Mag();

        if (debug) cout << "Getting Detector Positions" << endl;
        string iLeadDetector = DetectorPosition(&currentTree, leadindex);
        string iSubleadDetector = DetectorPosition(&currentTree, subleadindex);
        string iConvDetector = DetectorPosition(&currentTree, photonconvindex);
        string selection = "All";

        if (debug) cout << "Getting Photon Category" << endl;
        int leadPhoCategory = photonCategory(currentTree.pho_r9[leadindex], currentTree.pho_isEB[leadindex], currentTree.pho_isEE[leadindex]);
        int subleadPhoCategory = photonCategory(currentTree.pho_r9[subleadindex], currentTree.pho_isEB[subleadindex], currentTree.pho_isEE[subleadindex]);

        bool convsel1 = false;
        bool convsel2 = false;

        if (debug) cout << "Checking for Good Conversions" << endl;
        if (photonconvindex==leadindex && validconversion) convsel1 = convSel(currentTree.conv_ntracks[convindex],
                                                           currentTree.conv_validvtx[convindex] ,  
                                                           currentTree.conv_chi2_probability[convindex], 
                                                           currentTree.conv_dphitrksatvtx[convindex], 
                                                           currentTree.conv_paircotthetasep[convindex], 
                                                           EoP,
                                                           ConversionVertex[convindex].Perp());
	
        if (photonconvindex==subleadindex && validconversion) convsel2 = convSel(currentTree.conv_ntracks[convindex],
                                                              currentTree.conv_validvtx[convindex] ,  
                                                              currentTree.conv_chi2_probability[convindex], 
                                                              currentTree.conv_dphitrksatvtx[convindex], 
                                                              currentTree.conv_paircotthetasep[convindex], 
                                                              EoP,
                                                              ConversionVertex[convindex].Perp());

        int diPhoCategory = diPhotonCategory( leadPhoCategory, subleadPhoCategory );

        double NewZconvlinear = 0;
        ////////////////////////////////////
        if (debug) cout << "Filling Conversion Plots" << endl;
        if (convsel1 || convsel2) {
          if (!data) FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
          histoContainer->Fill("convdZvsEtaAll",Photonxyz[photonconvindex].Eta(),currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          histoContainer->Fill("convEoPAll",iConvDetector,SuperClusterp4[convscindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          histoContainer->Fill("convr",iConvDetector,ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("Zconv",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex],weight);
          if (!data) {
            if (debug) cout << "Filling DeltaZ resolution plots." << endl;
            histoContainer->Fill("convdZPrimary",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),PrimaryVertex[0].Z()-SimVertex.Z(),weight);
            histoContainer->Fill("ZconvdZ",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            histoContainer->Fill("ZconvdZVsNumOfVertices",iConvDetector,float(currentTree.vtx_std_xyz->GetSize()),currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
	
            if ( fabs(currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < 1 ) histoContainer->Fill("ZconvdZLessOnecm",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            if ( fabs(currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < .5 ) histoContainer->Fill("ZconvdZLessHalfcm",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            if ( fabs(currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < .3 ) histoContainer->Fill("ZconvdZLess3mm",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
	    
          }
          if (debug) cout << "Filling Refitted Tracks" << endl;
          histoContainer->Fill("RefittedPt",iConvDetector, ConversionRefittedPairMomentum[convindex].Pt(), weight);
          histoContainer->Fill("RefittedEta",iConvDetector, ConversionRefittedPairMomentum[convindex].Eta(), weight);
          histoContainer->Fill("RefittedPhi",iConvDetector, ConversionRefittedPairMomentum[convindex].Phi(), weight);

          if (debug) cout << "Finding new vertex Z." << endl;
          TVector3 myBeamSpot(0,0,0);
          double NewZ = FindNewdZ(Photonxyz[photonconvindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          double RefittedZ = FindRefittedZ(ConversionVertex[convindex], ConversionRefittedPairMomentum[convindex]);
          NewZconvlinear = FindNewZConvLinear(ConversionVertex[convindex],Photonxyz[photonconvindex],PrimaryVertex[0]);
          myBeamSpot.SetXYZ(PrimaryVertex[0].X(),PrimaryVertex[0].Y(),0);
          double NewZPV = FindNewdZ(Photonxyz[photonconvindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          
          float deltaz = 100000;
          for (unsigned int j=0; j<PrimaryVertex.size(); j++) {
            if (deltaz > fabs(NewZconvlinear-PrimaryVertex[j].Z())) {
              deltaz = fabs(NewZconvlinear-PrimaryVertex[j].Z());
              nearvertexindex = j;
            }
          }

          if (debug) cout << "Filling Z Plots" << endl;
          histoContainer->Fill("RefittedZconv",iConvDetector, RefittedZ, weight);
          histoContainer->Fill("NewZconv",iConvDetector, NewZ, weight);
          histoContainer->Fill("NewZPVconv",iConvDetector, NewZPV, weight);
          histoContainer->Fill("NewZconvlinear",iConvDetector, NewZconvlinear, weight);

          if (!data) {
            if (debug) cout << "Filling DeltaZ Plots" << endl;
            histoContainer->Fill("ZconvdZNearest",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-PrimaryVertex[nearvertexindex].Z(),weight);
            histoContainer->Fill("ZdZNear",iConvDetector,PrimaryVertex[nearvertexindex].Z()-SimVertex.Z(),weight);
            histoContainer->Fill("RefittedZconvdZ",iConvDetector, RefittedZ-SimVertex.Z(), weight);
            histoContainer->Fill("NewZconvdZ",iConvDetector, NewZ-SimVertex.Z(), weight);
            histoContainer->Fill("NewZPVconvdZ",iConvDetector, NewZPV-SimVertex.Z(), weight);
            histoContainer->Fill("NewZconvlineardZ",iConvDetector, NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<1) histoContainer->Fill("ZconvdZLinearLessOnecm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<0.5) histoContainer->Fill("ZconvdZLinearLessHalfcm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<0.3) histoContainer->Fill("ZconvdZLinearLess3mm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
            if (iConvDetector=="Barrel") {
              FilldZTrackerBarrel(histoContainer, "NewPVdZ", NewZPV-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZPV-SimVertex.Z())<1.0) FilldZTrackerBarrel(histoContainer, "Less10mmNewPVdZ", NewZPV-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZPV-SimVertex.Z())<0.5) FilldZTrackerBarrel(histoContainer, "Less5mmNewPVdZ", NewZPV-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZPV-SimVertex.Z())<0.3) FilldZTrackerBarrel(histoContainer, "Less3mmNewPVdZ", NewZPV-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              FilldZTrackerBarrel(histoContainer, "LineardZ", NewZconvlinear-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<1.0) FilldZTrackerBarrel(histoContainer, "Less10mmLineardZ", NewZconvlinear-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<0.5) FilldZTrackerBarrel(histoContainer, "Less5mmLineardZ", NewZconvlinear-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<0.3) FilldZTrackerBarrel(histoContainer, "Less3mmLineardZ", NewZconvlinear-SimVertex.Z(), ConversionVertex[convindex].Perp(), weight);
            }
            if (iConvDetector=="Endcap") {
              FilldZTrackerEndcap(histoContainer, "NewPVdZ", NewZPV-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZPV-SimVertex.Z())<1.0) FilldZTrackerEndcap(histoContainer, "Less10mmNewPVdZ", NewZPV-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZPV-SimVertex.Z())<0.5) FilldZTrackerEndcap(histoContainer, "Less5mmNewPVdZ", NewZPV-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZPV-SimVertex.Z())<0.3) FilldZTrackerEndcap(histoContainer, "Less3mmNewPVdZ", NewZPV-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              FilldZTrackerEndcap(histoContainer, "LineardZ", NewZconvlinear-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<1.0) FilldZTrackerEndcap(histoContainer, "Less10mmLineardZ", NewZconvlinear-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<0.5) FilldZTrackerEndcap(histoContainer, "Less5mmLineardZ", NewZconvlinear-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
              if (abs(NewZconvlinear-SimVertex.Z())<0.3) FilldZTrackerEndcap(histoContainer, "Less3mmLineardZ", NewZconvlinear-SimVertex.Z(), abs(ConversionVertex[convindex].Z()), weight);
            }
          }
          
          if (validconversion) {
            if (debug) cout << "Filling Combined Conversion Plots" << endl;
            histoContainer->Fill("trackerconvr",iConvDetector,ConversionVertex[convindex].Perp(),weight);
            histoContainer->Fill("trackerconvEoPAll",iConvDetector,SuperClusterp4[convscindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
            histoContainer->Fill("Ztrackerconv",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex],weight);
            histoContainer->Fill("ZtrackerconvdZ",iConvDetector,currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          }
          
        }

        if (debug) cout << "Doing Vertex Matching" << endl;
        TVector3 NearVertex = PrimaryVertex[nearvertexindex];
        bool PrimaryVertexMatched = DoVertexMatching(PrimaryVertex[0],SimVertex);
        bool NearVertexMatched = DoVertexMatching(NearVertex,SimVertex);

        if (debug) cout << "Filling Vertex Matched Plots" << endl;
        histoContainer->Fill("PrimaryVertexMatchedAllEcal", PrimaryVertexMatched, weight);
        histoContainer->Fill("PrimaryVertexMatched", iConvDetector, PrimaryVertexMatched, weight);
        if (validconversion) histoContainer->Fill("NearVertexMatchedAllEcal", NearVertexMatched, weight);
        if (validconversion) histoContainer->Fill("NearVertexMatched", iConvDetector, NearVertexMatched, weight);

        /*if (Photonp4[subleadindex].Pt()>Photonp4[leadindex].Pt()) {
          cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
          cout << "LeadPt is " << Photonp4[leadindex].Pt() << endl;
          cout << "SubLeadPt is " << Photonp4[subleadindex].Pt() << endl;
          }*/

        if (debug) cout << "Filling Conversion Hists" << endl;
        histoContainer->Fill("NPhotonsAll",currentTree.pho_n,weight);
        //// Fill histograms hists only for conversions satisfying the conversion selection
        if (convsel1 || convsel2) FillConvHists(&currentTree, histoContainer, photonconvindex, convindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);

        if (debug) cout << "Filling Photon Hists" << endl;
        FillPhotonHists(&currentTree, histoContainer, string("lead"), leadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        FillPhotonHists(&currentTree, histoContainer, string("sublead"), subleadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);

        if (validconversion && currentTree.conv_ntracks[convindex]==2 && currentTree.conv_chi2_probability[convindex]>0.0005 && (bool) currentTree.pho_isEB[photonconvindex])
          histoContainer->Fill("convVtxRvsZBarrelAll",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
        
        if (debug) cout << "Making Vertex Vectors" << endl;
        TLorentzVector VLead(Photonp4[leadindex]);
        TLorentzVector VSubLead(Photonp4[subleadindex]);
        TLorentzVector VSum=VLead+VSubLead;
        double InvMass=fabs(VSum.M());
        double cos_thetastar = CosThetaStar(VLead,VSum);
        TLorentzVector SimVLead = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], SimVertex.Z());
        TLorentzVector SimVSubLead = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], SimVertex.Z());
        TLorentzVector SimVSum=SimVLead+SimVSubLead;
        double SimInvMass=fabs(SimVSum.M());
        
        //With New Vertex
        TLorentzVector VSum_newvertex(0,0,0,0);
        TLorentzVector VLead_newvertex(0,0,0,0);
        TLorentzVector VSubLead_newvertex(0,0,0,0);
        double InvMass_newvertex = 0;
        double cos_thetastar_newvertex = 0;

        //With New Linear Vertex Fit
        TLorentzVector VSum_linearvertex(0,0,0,0);
        TLorentzVector VLead_linearvertex(0,0,0,0);
        TLorentzVector VSubLead_linearvertex(0,0,0,0);
        double InvMass_linearvertex = 0;
        double cos_thetastar_linearvertex = 0;
        
        //With Sim Vertex
        TLorentzVector VSum_simvertex(0,0,0,0);
        TLorentzVector VLead_simvertex(0,0,0,0);
        TLorentzVector VSubLead_simvertex(0,0,0,0);
        double InvMass_simvertex = 0;
        double cos_thetastar_simvertex = 0;

        //With Nearest Vertex
        TLorentzVector VSum_nearvertex(0,0,0,0);
        TLorentzVector VLead_nearvertex(0,0,0,0);
        TLorentzVector VSubLead_nearvertex(0,0,0,0);
        double InvMass_nearvertex = 0;
        double cos_thetastar_nearvertex = 0;

        if (convsel1 || convsel2) {
          if (debug) cout << "Filling Conversion Vertex Hists" << endl;
          VLead_newvertex = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], currentTree.conv_zofprimvtxfromtrks[convindex]);
          VSubLead_newvertex = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], currentTree.conv_zofprimvtxfromtrks[convindex]);
          VSum_newvertex=VLead_newvertex+VSubLead_newvertex;
          InvMass_newvertex=fabs(VSum_newvertex.M());
          cos_thetastar_newvertex= CosThetaStar(VLead_newvertex,VSum_newvertex);

          VLead_linearvertex = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], NewZconvlinear);
          VSubLead_linearvertex = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], NewZconvlinear);
          VSum_linearvertex=VLead_linearvertex+VSubLead_linearvertex;
          InvMass_linearvertex=fabs(VSum_linearvertex.M());
          cos_thetastar_linearvertex= CosThetaStar(VLead_linearvertex,VSum_linearvertex);
          
          if (!data) {
            VLead_simvertex = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], SimVertex.Z());
            VSubLead_simvertex = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], SimVertex.Z());
            VSum_simvertex=VLead_simvertex+VSubLead_simvertex;
            InvMass_simvertex=fabs(VSum_simvertex.M());
            cos_thetastar_simvertex= CosThetaStar(VLead_simvertex,VSum_simvertex);
          }
          
          VLead_nearvertex = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], NearVertex.Z());
          VSubLead_nearvertex = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], NearVertex.Z());
          VSum_nearvertex=VLead_nearvertex+VSubLead_nearvertex;
          InvMass_nearvertex=fabs(VSum_nearvertex.M());
          cos_thetastar_nearvertex= CosThetaStar(VLead_nearvertex,VSum_nearvertex);
          
        }
        
        if (Photonp4[leadindex].Pt()>40 && Photonp4[subleadindex].Pt()>30 && InvMass>90 && InvMass<250
            && MarcosCut(Photonp4[leadindex].Pt(), currentTree.pho_ecalsumetconedr04[leadindex], currentTree.pho_hcalsumetconedr04[leadindex], currentTree.pho_trksumpthollowconedr04[leadindex], currentTree.pho_haspixseed[leadindex], currentTree.pho_isEB[leadindex], currentTree.pho_isEE[leadindex], currentTree.pho_sieie[leadindex], currentTree.pho_hoe[leadindex])
            && MarcosCut(Photonp4[subleadindex].Pt(), currentTree.pho_ecalsumetconedr04[subleadindex], currentTree.pho_hcalsumetconedr04[subleadindex], currentTree.pho_trksumpthollowconedr04[subleadindex], currentTree.pho_haspixseed[subleadindex], currentTree.pho_isEB[subleadindex], currentTree.pho_isEE[subleadindex], currentTree.pho_sieie[subleadindex], currentTree.pho_hoe[subleadindex])
            ) {
          histoContainer->Fill("leadEtMarco_allEcal",Photonp4[leadindex].Pt(),weight);
          histoContainer->Fill("subleadEtMarco_allEcal",Photonp4[subleadindex].Pt(),weight);
          histoContainer->Fill("mass_Marco_allEcal",InvMass,weight);

          int MarcosCategory = MarcosCutCategory(currentTree.pho_r9[leadindex], currentTree.pho_isEB[leadindex], currentTree.pho_isEE[leadindex], currentTree.pho_r9[subleadindex], currentTree.pho_isEB[subleadindex], currentTree.pho_isEE[subleadindex]);
          TString MarcosHistName = "leadEtMarco_cat";
          MarcosHistName += MarcosCategory;
          histoContainer->Fill(MarcosHistName.Data(),Photonp4[leadindex].Pt(),weight);
          MarcosHistName = "subleadEtMarco_cat";
          MarcosHistName += MarcosCategory;
          histoContainer->Fill(MarcosHistName.Data(),Photonp4[subleadindex].Pt(),weight);
          MarcosHistName = "mass_Marco_cat";
          MarcosHistName += MarcosCategory;
          histoContainer->Fill(MarcosHistName.Data(),InvMass,weight);

        }

        /// di-photon system before event selection
        string HiggsInWhichDetector = HiggsDetectorPosition(&currentTree, leadindex, subleadindex);

        if (convsel1 || convsel2) {
          FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass_nearvertex, SimInvMass, cos_thetastar);
          FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        } else {
          FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, SimInvMass, cos_thetastar);
        }
          
        ///////////////////////////////  Event selection ///////////////////////////////////////
        ///////////// build the di-photon system
        if (debug) cout << "Doing Selection" << endl;
        if (Photonp4[leadindex].Pt()<40.0) continue; // leading photon
        if (Photonp4[subleadindex].Pt()<30.0) continue; // subleading photon
        if (fabs(Photonxyz[leadindex].Eta())>2.5 || fabs(Photonxyz[subleadindex].Eta())>2.5) continue;

        //isolation 
//         if (!(tightId(Photonp4[leadindex].Pt(),
//                       currentTree.pho_ecalsumetconedr04[leadindex],
//                       currentTree.pho_hcalsumetconedr04[leadindex],
//                       currentTree.pho_trksumpthollowconedr04[leadindex],
//                       (bool) currentTree.pho_isEB[leadindex],
//                       (bool) currentTree.pho_isEE[leadindex],
//                       currentTree.pho_sieie[leadindex],
//                       currentTree.pho_hoe[leadindex]))) continue;

//         if (!(tightId(Photonp4[subleadindex].Pt(),
//                       currentTree.pho_ecalsumetconedr04[subleadindex],
//                       currentTree.pho_hcalsumetconedr04[subleadindex],
//                       currentTree.pho_trksumpthollowconedr04[subleadindex],
//                       (bool) currentTree.pho_isEB[subleadindex],
//                       (bool) currentTree.pho_isEE[subleadindex],
//                       currentTree.pho_sieie[subleadindex],
//                       currentTree.pho_hoe[subleadindex]))) continue;

         //New CiC Cut 4 is SuperTight
        if (debug) for (unsigned int j=0; j<currentTree.pho_cic4cutlevel_lead->size(); j++) cout << "Index is: " << j << " and lead cut level is " << currentTree.pho_cic4cutlevel_lead->at(j) << endl;
        if (debug) for (unsigned int j=0; j<currentTree.pho_cic4cutlevel_sublead->size(); j++) cout << "Index is: " << j << " and sublead cut level is " << currentTree.pho_cic4cutlevel_sublead->at(j) << endl;
        if (currentTree.pho_cic4cutlevel_lead->at(leadindex)<4 && currentTree.pho_cic4cutlevel_sublead->at(subleadindex)<4) continue;
        
        histoContainer->Fill("NPhotonsSel",currentTree.pho_n,weight);
        if (currentTree.pho_haspixseed[leadindex]==1 && convsel1==false) continue;
        if (currentTree.pho_haspixseed[subleadindex]==1 && convsel2==false) continue;
        selection = "Sel";

        if (debug) cout << "Filling Photon Hists" << endl;
        FillPhotonHists(&currentTree, histoContainer, string("lead"), leadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        FillPhotonHists(&currentTree, histoContainer, string("sublead"), subleadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        if (convsel1 || convsel2) FillConvHists(&currentTree, histoContainer, photonconvindex, convindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);

        if (convsel1 || convsel2) {
          string diPhoSel = "DiPho";	  
          if (convsel1 || convsel2) FillConvHistsInDiPho(&currentTree, histoContainer, photonconvindex, convindex, diPhoSel, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
	  
          histoContainer->Fill("convdZvsEtaSel",Photonxyz[photonconvindex].Eta(),currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          if (!data) {
            histoContainer->Fill("convdZvsR",iConvDetector,ConversionVertex[convindex].Perp(),currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            histoContainer->Fill("convdZvsZ",iConvDetector,ConversionVertex[convindex].z(),currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            histoContainer->Fill("convdZvsRlinear",iConvDetector,ConversionVertex[convindex].Perp(),NewZconvlinear-SimVertex.Z(),weight);
            histoContainer->Fill("convdZvsZlinear",iConvDetector,ConversionVertex[convindex].z(),NewZconvlinear-SimVertex.Z(),weight);
            FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
          }
        }

        if (validconversion && currentTree.conv_ntracks[convindex]==2 && currentTree.conv_chi2_probability[convindex]>0.0005 && (bool) currentTree.pho_isEB[photonconvindex])
          histoContainer->Fill("convVtxRvsZBarrelSel",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);

        if (convsel1) histoContainer->Fill("selconvVtxRvsZBarrelSel",ConversionVertex[leadindex].Z(),sqrt(ConversionVertex[leadindex].X()*ConversionVertex[leadindex].X()+ConversionVertex[leadindex].Y()*ConversionVertex[leadindex].Y()));
        if (convsel2) histoContainer->Fill("selconvVtxRvsZBarrelSel",ConversionVertex[subleadindex].Z(),sqrt(ConversionVertex[subleadindex].X()*ConversionVertex[subleadindex].X()+ConversionVertex[subleadindex].Y()*ConversionVertex[subleadindex].Y()));

        FillCatHists(&currentTree, histoContainer, weight, leadindex, subleadindex, leadPhoCategory, subleadPhoCategory);

        if (convsel1 || convsel2) {
          FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass_nearvertex, SimInvMass, cos_thetastar);
          FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        } else {
          FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, SimInvMass, cos_thetastar);
        }

        if ((convsel1 || convsel2) && InvMass>100) {
          TLorentzVector MatchedSuperCluster = SuperClusterp4[currentTree.pho_scind[photonconvindex]];
          double deltaeta = Photonxyz[photonconvindex].Eta() - etaTransformation(ConversionRefittedPairMomentum[convindex].Eta(),currentTree.conv_zofprimvtxfromtrks[convindex]);
          double deltaphi = DeltaPhi(Photonxyz[photonconvindex].Phi(),ConversionRefittedPairMomentum[convindex].Phi());
          double EoP = MatchedSuperCluster.E()/ConversionRefittedPairMomentum[convindex].Mag();
         
          histoContainer->Fill("goodconveta",ConversionVertex[convindex].Eta(),weight);
          histoContainer->Fill("goodconvr",iConvDetector,ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("goodconvchi2",iConvDetector,currentTree.conv_chi2_probability[convindex],weight);
          histoContainer->Fill("goodconvz",iConvDetector,ConversionVertex[convindex].Z(),weight);
          histoContainer->Fill("goodconvdistmin",iConvDetector,currentTree.conv_distofminapproach[convindex],weight);
          histoContainer->Fill("goodconveop",iConvDetector,EoP,weight);
          histoContainer->Fill("goodconvdeta",iConvDetector,deltaeta,weight);
          histoContainer->Fill("goodconvdphi",iConvDetector,deltaphi,weight);
          if (!data) histoContainer->Fill("higgsPtWhenAtLeastOneGoodConv",Higgs_p4.Pt(),weight);
          float sigma_z=0;
          if (iConvDetector=="Barrel") {

            if ( ConversionVertex[convindex].Perp() <=15 ) {
              sigma_z = 0.06;
            } else if (ConversionVertex[convindex].Perp()>15 && ConversionVertex[convindex].Perp()<=60) {
              sigma_z= 0.67;  
            } else if (ConversionVertex[convindex].Perp()>60 ) {
              sigma_z=2;
            }

          } else if ( iConvDetector=="Endcap") {

            if ( fabs(ConversionVertex[convindex].Z()) <=50 ) {
              sigma_z = 0.18;
            } else if ( fabs(ConversionVertex[convindex].Z())>50 && fabs(ConversionVertex[convindex].Z())<=100) {
              sigma_z = 0.61;
            } else if ( fabs(ConversionVertex[convindex].Z())>100 ) {
              sigma_z = 0.99;
            }

          }

          if ( !data && fabs(NewZconvlinear-SimVertex.Z()) < 1 )	  
            histoContainer->Fill("higgsPtWhenAtLeastOneGoodConvZInOneSigma",Higgs_p4.Pt(),weight);

        }
        
        //GenMatching
        if (data) continue;
        if (debug) cout << "Doing Gen Matching" << endl;
        unsigned int LeadGenMatchedIndex = 999999;
        unsigned int SubleadGenMatchedIndex = 999999;

        LeadGenMatchedIndex = DoGenMatching(&currentTree, Photonxyz[leadindex]);
        SubleadGenMatchedIndex = DoGenMatching(&currentTree, Photonxyz[subleadindex]);

        if (LeadGenMatchedIndex!=999999) {
          TVector3 GenPhotonxyz = (*((TVector3*) currentTree.gp_vtx->At(LeadGenMatchedIndex)));
          TLorentzVector GenPhotonp4 = (*((TLorentzVector*) currentTree.gp_p4->At(LeadGenMatchedIndex)));
          TLorentzVector Correctedp4 = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], GenPhotonp4.Z());
          histoContainer->Fill("MomentumResolution",iLeadDetector,(Correctedp4.E()-GenPhotonp4.E())/GenPhotonp4.E(),weight);
        }
        if (SubleadGenMatchedIndex!=999999) {
          TVector3 GenPhotonxyz = (*((TVector3*) currentTree.gp_vtx->At(SubleadGenMatchedIndex)));
          TLorentzVector GenPhotonp4 = (*((TLorentzVector*) currentTree.gp_p4->At(SubleadGenMatchedIndex)));
          TLorentzVector Correctedp4 = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], GenPhotonp4.Z());
          histoContainer->Fill("MomentumResolution",iSubleadDetector,(Correctedp4.E()-GenPhotonp4.E())/GenPhotonp4.E(),weight);
        }

        if (LeadGenMatchedIndex!=999999 && SubleadGenMatchedIndex!=999999) {
          selection="Matched";
          if (convsel1 || convsel2) {
            FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass_nearvertex, SimInvMass, cos_thetastar);
            FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
            FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
            FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
            if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
            FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
          } else {
            FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, SimInvMass, cos_thetastar);
          }
          
          TLorentzVector LeadGenPhoton = (*((TLorentzVector*) currentTree.gp_p4->At(LeadGenMatchedIndex)));
          TLorentzVector SubleadGenPhoton = (*((TLorentzVector*) currentTree.gp_p4->At(SubleadGenMatchedIndex)));
          TLorentzVector GenSum = LeadGenPhoton + SubleadGenPhoton;
          double GenInvMass = GenSum.M();
          
          histoContainer->Fill("MassResolution",HiggsInWhichDetector,InvMass-GenInvMass,weight);
          histoContainer->Fill("MassResolutionSim",HiggsInWhichDetector,SimInvMass-GenInvMass,weight);
          histoContainer->Fill("MassResolutionGen",HiggsInWhichDetector,GenInvMass,weight);
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

bool DoVertexMatching(TVector3 TestVertex, TVector3 SimVertex) {

  bool ReturnValue = false;
  if (abs(TestVertex.Z()-SimVertex.Z()) < 0.2) ReturnValue=true;
  return ReturnValue;
  
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

unsigned int getconvindex(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex) {
  if (currentTree->pho_conv_validvtx[leadindex]) return leadindex;
  if (currentTree->pho_conv_validvtx[subleadindex]) return subleadindex;
  return 0;
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

double CosThetaStar(TLorentzVector VLead, TLorentzVector VSum) {

  double beta_b = VSum.Beta();
  double gamma_b = VSum.Gamma();
  TVector3 directionV = VSum.Vect().Unit();

  TVector3 CrossVLead = VLead.Vect().Cross(directionV);
  double DotVLeadValue = VLead.Vect().Dot(directionV);
  double CrossVLeadValue = sqrt(CrossVLead.x()*CrossVLead.x()+CrossVLead.y()*CrossVLead.y()+CrossVLead.z()*CrossVLead.z());
  double sin_theta =  CrossVLeadValue/VLead.E();
  double cos_theta =  DotVLeadValue/VLead.E(); 
  double tg_thetas = sin_theta/(gamma_b*(cos_theta-beta_b));
  double cos_thetastar = 1.0/sqrt(1.0+tg_thetas*tg_thetas);
  return cos_thetastar;

}

double DeltaPhi(double Phi1, double Phi2) {
  double Pi = 3.14159265;
  double deltaphi = Phi1 - Phi2;
  if (deltaphi>Pi) deltaphi = 2*Pi-deltaphi;
  if (deltaphi<-Pi) deltaphi = 2*Pi+deltaphi;
  return deltaphi;
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


double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot) {

  double dz = (vtx.z()-myBeamSpot.z()) - ((vtx.x()-myBeamSpot.x())*mom.x()+(vtx.y()-myBeamSpot.y())*mom.y())/mom.Perp() * mom.z()/mom.Perp();
  return dz + myBeamSpot.z();
  
}

double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex) {
  
  double deltaX1 = superclustervtx.X()-convvtx.X();
  double deltaY1 = superclustervtx.Y()-convvtx.Y();
  double deltaZ1 = superclustervtx.Z()-convvtx.Z();
  double R1 = sqrt(deltaX1*deltaX1+deltaY1*deltaY1);
  double tantheta = R1/deltaZ1;
  
  double deltaX2 = convvtx.X()-primaryvertex.X();
  double deltaY2 = convvtx.Y()-primaryvertex.Y();
  double R2 = sqrt(deltaX2*deltaX2+deltaY2*deltaY2);
  double deltaZ2 = R2/tantheta;
  double primaryvertexZ = superclustervtx.Z()-deltaZ1-deltaZ2;
  return primaryvertexZ;

}

double FindRefittedZ(TVector3 ConversionVertex, TVector3 ConversionRefittedPairMomentum) {

  double theZOfPrimaryVertexFromTracks=-999999.;
  double theta=ConversionRefittedPairMomentum.Theta();

  theZOfPrimaryVertexFromTracks = ConversionVertex.z() - sqrt(ConversionVertex.Perp2())*(1./tan(theta));

  return  theZOfPrimaryVertexFromTracks;
  
}

string HiggsDetectorPosition(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex) {

  string ReturnValue = "";
  if (currentTree->pho_isEB[leadindex] && currentTree->pho_isEB[subleadindex]) ReturnValue="Barrel";
  if ((currentTree->pho_isEB[leadindex] && currentTree->pho_isEE[subleadindex]) ||
      (currentTree->pho_isEE[leadindex] && currentTree->pho_isEB[subleadindex]) ||
      (currentTree->pho_isEE[leadindex] && currentTree->pho_isEE[subleadindex]) ) ReturnValue="Endcap";
  return ReturnValue;
  
}

string DetectorPosition(sdaReader *currentTree, unsigned int index) {

  string ReturnValue="";
  if (currentTree->pho_isEB[index]) ReturnValue="Barrel";
  if (currentTree->pho_isEE[index]) ReturnValue="Endcap";
  return ReturnValue;
 
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

TString MakeFileName(string filename, bool unweighted, string dataweight, double RCut) {

  TString outfilename = "";
  
  if (unweighted) {
    outfilename = "Unweighted";
    outfilename += filename;
  } else {
    outfilename = dataweight;
    outfilename += filename;
  }
  if (RCut!=999999) {
    
    outfilename.ReplaceAll(".root","RCut.root");
  }
  return outfilename;
  
}

TLorentzVector CalcDiffVertex(TLorentzVector p4, TVector3 xyz, float newz) {

  TLorentzVector Return4Vector(0,0,0,0);
  double theta = atan2(xyz.Perp(),xyz.Z()-newz);
  double eta = -log(tan(theta/2));
  double pt = fabs(p4.E()*sin(theta));

  Return4Vector.SetPtEtaPhiE(pt,eta,xyz.Phi(),p4.E());
  return Return4Vector;

}

void BookConversionPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {
  vector <TString> region;
  region.push_back("allEcal");
  region.push_back("Barrel");
  region.push_back("Endcap");
  for (unsigned int i=0; i<2; i++) {
    for (unsigned int j=0; j<region.size(); j++) {
      TString histnametemp(histname);
      TString histtitletemp(histtitle);
      histnametemp.ReplaceAll("region",region[j]);
      histtitletemp.ReplaceAll("region",region[j]);
      if (i==1) {
        histnametemp.ReplaceAll("All","Sel");
        histtitletemp.ReplaceAll("Photon","Selected Photon");
      }
      histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
    }
  }
}

void BookConversionPlotsDiPho(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {
  vector <TString> region;
  region.push_back("allEcal");
  region.push_back("Barrel");
  region.push_back("Endcap");

  for (unsigned int j=0; j<region.size(); j++) {
    TString histnametemp(histname);
    TString histtitletemp(histtitle);
    histnametemp.ReplaceAll("region",region[j]);
    histtitletemp.ReplaceAll("region",region[j]);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
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

void BookFourHists(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {
  
  for (int i=0; i<4; i++) {
    char Catagory[255];
    TString histnametemp(histname);
    TString histtitletemp(histtitle);
    sprintf(Catagory,"%d",i);
    histnametemp.ReplaceAll("dummy",Catagory);
    histtitletemp.ReplaceAll("dummy",Catagory);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
  }
  
}

void BookHistograms(HistoContainer *histoContainer) {

  histoContainer->Add("NumberVertices","Number of Reconstructed Vertices;Number of Vertices; Counts",20,-0.5,19.5);
  histoContainer->Add("NumberSimVertices","Number of Simulated Vertices;Number of Sim Vertices; Counts",20,-0.5,19.5);

  histoContainer->Add("NumberSimVerticesVsRecoVertices","Number of Simulated Vertices;RecoVertices;Number of Sim Vertices;",20,-0.5,19.5,10,0,10);

  histoContainer->Add("ZVertex","Z of Primary Vertex;Z (cm); Counts",100,-20,20);
  histoContainer->Add("ZVertices","Z of all Vertecies;Z (cm); Counts",100,-20,20);
  histoContainer->Add("ZSimVertex","Z of Simulated Vertex;Z (com); Counts",100,-20,20);
    
  histoContainer->Add("ZdZ","#deltaZ between the Z of the Primary Vertex and the Sim Vertex;Z (cm); Counts",100,-1,1);
  histoContainer->Add("ZdZZoom","#deltaZ between the Z of the Primary Vertex and the Sim Vertex: Zoom;Z (cm); Counts",100,-0.02,0.02);
  histoContainer->Add("ZdZAll","#deltaZ between the Z of the all vertecies and the Sim Vertex;Z (cm); Counts",100,-1,1);

  BookBarrelAndEndcap(histoContainer,"RefittedPt","Pt distrobution of refitted momentum;Pt (GeV);Counts",100,0,300);
  BookBarrelAndEndcap(histoContainer,"RefittedEta","#eta distrobution of refitted momentum;#eta (GeV);Counts",60,-3,3);
  BookBarrelAndEndcap(histoContainer,"RefittedPhi","#phi distrobution of refitted momentum;#phi (GeV);Counts",64,-3.2,3.2);
  
  BookRCutsdZPlots(histoContainer,"dz_conv_");
  BookTrackerdZ(histoContainer,"NewPVdZ");
  BookTrackerdZ(histoContainer,"Less10mmNewPVdZ");
  BookTrackerdZ(histoContainer,"Less5mmNewPVdZ");
  BookTrackerdZ(histoContainer,"Less3mmNewPVdZ");
  BookTrackerdZ(histoContainer,"LineardZ");
  BookTrackerdZ(histoContainer,"Less10mmLineardZ");
  BookTrackerdZ(histoContainer,"Less5mmLineardZ");
  BookTrackerdZ(histoContainer,"Less3mmLineardZ");
  
  histoContainer->Add("PrimaryVertexMatchedAllEcal","#deltaZ between the Primary Vertex and SimVertex less than 2mm: AllEcal",2,0,2);
  histoContainer->Add("NearVertexMatchedAllEcal","#deltaZ between the Near Vertex and SimVertex less than 2mm: AllEcal",2,0,2);
  BookBarrelAndEndcap(histoContainer,"PrimaryVertexMatched","#deltaZ between the Primary Vertex and SimVertex less than 2mm: AllEcal",2,0,2);
  BookBarrelAndEndcap(histoContainer,"NearVertexMatched","#deltaZ between the Near Vertex and SimVertex less than 2mm: region",2,0,2);
  
  BookBarrelAndEndcap(histoContainer,"convr","R of conversion; R (cm): region; Counts",100,0,100);
  BookBarrelAndEndcap(histoContainer,"convEoPAll","E over P of Conversion; E over P; Counts",100,0,3);
  BookBarrelAndEndcap(histoContainer,"Zconv","Z of Primary Vertex from Conversion: region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"ZconvdZ","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-5,5);

  histoContainer->Add("ZconvdZVsNumOfVerticesBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: Barrel;# of PU vertices;dZ (cm)",20,0,20,100,-10,10);
  histoContainer->Add("ZconvdZVsNumOfVerticesEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: Endcap;# of PU vertices;dZ (cm)",20,0,20,100,-10,10);
   
  BookBarrelAndEndcap(histoContainer,"ZconvdZLessOnecm","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-10,10);
  BookBarrelAndEndcap(histoContainer,"ZconvdZLessHalfcm","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"ZconvdZLess3mm","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-5,5);

  BookBarrelAndEndcap(histoContainer,"ZconvdZLinearLessOnecm","#deltaZ between the Z of the Primary Vertex from Conversion Linear Method and Sim Vertex: region;Z (cm); Counts",100,-10,10);
  BookBarrelAndEndcap(histoContainer,"ZconvdZLinearLessHalfcm","#deltaZ between the Z of the Primary Vertex from Conversion Linear Method and Sim Vertex: region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"ZconvdZLinearLess3mm","#deltaZ between the Z of the Primary Vertex from Conversion Linear Method and Sim Vertex: region;Z (cm); Counts",100,-5,5);

  BookBarrelAndEndcap(histoContainer,"ZconvdZNearest","#deltaZ between the Z of the Conversion and the nearest vertex: region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"ZdZNear","#deltaZ between the Z of the vertex nearest the conversion Z position and the Sim Vertex: region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"NewZconv","Josh's Z of Primary Vertex from Conversion (0,0,0): region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"NewZconvdZ","#deltaZ between the Josh's Z of the Primary Vertex from Conversion and Sim Vertex (0,0,0): region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"NewZPVconv","Josh's Z of Primary Vertex from Conversion (PV): region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"NewZPVconvdZ","#deltaZ between the Josh's Z of the Primary Vertex from Conversion and Sim Vertex (PV): region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"RefittedZconv","Z of Primary Vertex from Refitted Conversion: region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"RefittedZconvdZ","#deltaZ between the Refitted Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-5,5);
  BookBarrelAndEndcap(histoContainer,"NewZconvlinear","Z of Primary Vertex using linear method: region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"NewZconvlineardZ","#deltaZ between the Z of the Primary Vertex using linear method and the SimVertex: region;Z (cm); Counts",100,-5,5);

  BookBarrelAndEndcap(histoContainer,"goodconvr","R of conversion; R (cm): region; Counts",100,0,100);
  BookBarrelAndEndcap(histoContainer,"goodconvz","Z of conversion; Z (cm): region; Counts",100,0,200);
  BookBarrelAndEndcap(histoContainer,"goodconvchi2","#Chi^2 probability of vertex; #Chi^2 Probability; Counts",100,0,1);
  BookBarrelAndEndcap(histoContainer,"goodconvdistmin","Distantance of minimum approach of conversion tracks; cm; Counts",100,-1,1);
  BookBarrelAndEndcap(histoContainer,"goodconveop","E over P of conversion: region; E over P; Counts",100,0,3);
  BookBarrelAndEndcap(histoContainer,"goodconvdeta","#Delta#eta between the refitted tracks and supercluster; #Delta#eta; Counts",100,-0.1,0.1);
  BookBarrelAndEndcap(histoContainer,"goodconvdphi","#Delta#phi between the refitted tracks and supercluster; #Delta#phi; Counts",100,-0.1,0.1);
  histoContainer->Add("goodconveta","Conversion #eta; #eta; Counts",50,-2.5,2.5); 
  
  histoContainer->Add("convdZvsEtaAll","#deltaZ of the Primary Vertex from the Conversion and the Sim Vertex vs #eta;#eta of Conversion;#deltaZ of the Primary Vertex from the Conversion from the SimVertex(cm)",60, -3.0, 3.0, 100, -5, 5);
  histoContainer->Add("convdZvsEtaSel","#deltaZ of the Primary Vertex from the Conversion and the Sim Vertex vs #eta;#eta of Conversion;#deltaZ of the Primary Vertex from the Conversion from the SimVertex(cm)",60, -3.0, 3.0, 100, -5, 5);
 
  histoContainer->Add("convdZvsRBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Barrel;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsREndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Endcap;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsZBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus Z of Conversion: Barrel;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);
  histoContainer->Add("convdZvsZEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus Z of Conversion: Endcap;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);

  histoContainer->Add("convdZvsRlinearBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex Linear Fit versus R of Conversion: Barrel;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsRlinearEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex Linear Fit versus R of Conversion: Endcap;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsZlinearBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex Linear Fit versus Z of Conversion: Barrel;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);
  histoContainer->Add("convdZvsZlinearEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex Linear Fit versus Z of Conversion: Endcap;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);

  histoContainer->Add("convdZPrimaryBarrel","#deltaZ of the Primary Vertex from Conversion and Sim Vertex versus #deltaZ of the Primary Vertex and the Sim Vertex: Barrel;#deltaZ of Primary Vertex from Conversion and SimVertex (cm);#deltaZ of Primary Vertex and SimVertex (cm)",100, -5, 5, 100, -5, 5);
  histoContainer->Add("convdZPrimaryEndcap","#deltaZ of the Primary Vertex from Conversion and Sim Vertex versus #deltaZ of the Primary Vertex and the Sim Vertex: Endcap;#deltaZ of Primary Vertex from Conversion and SimVertex (cm);#deltaZ of Primary Vertex and SimVertex (cm)",100, -5, 5, 100, -5, 5);

  BookBarrelAndEndcap(histoContainer,"trackerconvr","R of tracker conversion; R (cm): region; Counts",100,0,100);
  BookBarrelAndEndcap(histoContainer,"trackerconvEoPAll","E over P of Tracker Conversion; E over P; Counts",100,0,3);
  BookBarrelAndEndcap(histoContainer,"Ztrackerconv","Z of Primary Vertex from Tracker Conversion: region;Z (cm); Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"ZtrackerconvdZ","#deltaZ between the Z of the Primary Vertex from Tracker Conversion and Sim Vertex: region;Z (cm); Counts",100,-5,5);
  
  histoContainer->Add("leadEtMarco_allEcal","Leading Photon Et with Marco's Cuts, Et (GeV); Counts",30,0.,150.);
  histoContainer->Add("subleadEtMarco_allEcal","Subleading Photon Et with Marco's Cuts, Et (GeV); Counts",30,0.,150.);
  histoContainer->Add("mass_Marco_allEcal","Invariant Mass of Photons with Marco's Cuts; Mass (GeV); Counts",20,100.,200.);

  BookFourHists(histoContainer,"leadEtMarco_catdummy","Leading Photon Et in Catagory dummy, Et (GeV); Counts",30,0,150);
  BookFourHists(histoContainer,"subleadEtMarco_catdummy","SubLeading Photon Et in Catagory dummy, Et (GeV); Counts",30,0,150);
  BookFourHists(histoContainer,"mass_Marco_catdummy","Invariant Mass of Photons in Catagory dummy; Mass (GeV); Counts",20,100,200);
    
  histoContainer->Add("NPhotonsAll","Num of photons in the event: all candidates",20,-0.5,19.5);
  histoContainer->Add("NPhotonsSel","Num of photons in the event: selected candidates",20,-0.5,19.5);
    
  BookPhysicsPlots(histoContainer,"leadPhoEtAll_region","leading photon Et, all candidates: region",100,0.,300.);
  BookPhysicsPlots(histoContainer,"leadPhoEtaAll_region","leading photon Eta, all candidates: region",100,-3.,3.);
  BookPhysicsPlots(histoContainer,"leadPhoPhiAll_region","leading photon Phi, all candidates: region",32,-3.2,3.2);
  BookPhysicsPlots(histoContainer,"leadPhoR9All_region","leading photon R9, all candidates: region",100,0.,1.1);
  BookPhysicsPlots(histoContainer,"leadPhoHoEAll_region","leading photon HoE, all candidates: region",100,0.,0.1);
  BookPhysicsPlots(histoContainer,"leadPhoTrkPtSumHollow04All_region","leading photon trk pt sum dr=04, all candidates: region",100,0.,2);
  BookPhysicsPlots(histoContainer,"leadPhoEcalPtSumHollow04All_region","leading photon ecal pt sum dr=04, all candidates: region",100,0.,5);
  BookPhysicsPlots(histoContainer,"leadPhoHcalPtSumHollow04All_region","leading photon hcal pt sum dr=04, all candidates: region",100,0.,5);
  BookPhysicsPlots(histoContainer,"leadPhoSigmaIetaIetaAll_region","leading photon SigmaIetaIeta, all candidates: region",100,0.,0.15);
  BookPhysicsPlots(histoContainer,"leadPhoZPVAll_region","leading photon Z(PV) (cm), all candidates: region",100,0.,20);
  BookPhysicsPlots(histoContainer,"leadPhoDZPVAll_region","leading photon #Deltaz_{Zpho - Ztrue} (cm), all candidates: region",100,-5.,5);
  BookPhysicsPlots(histoContainer,"leadPhoDZPVconvAll_region","leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: region",100,-5.,5);

  BookConversionPlots(histoContainer,"phi_conv_All_region","#phi of Photon Conversion; #phi of Conversion: region", 64, -3.2, 3.2);
  BookConversionPlots(histoContainer,"eta_conv_All_region","#eta of Photon Conversion; #eta of Conversion: region", 60, -3.0, 3.0);
  BookConversionPlots(histoContainer,"pt_conv_All_region","pt of Photon Conversion; pt of Conversion: region", 200, 0, 200);
  BookConversionPlots(histoContainer,"z_conv_All_region","z of Photon Conversion; z (cm): region", 200, -100, 100);
  BookConversionPlots(histoContainer,"r_conv_All_region","r of Photon Conversion; r (cm): region", 100, 0, 100);
  BookConversionPlots(histoContainer,"dz_conv_All_region","dZ of Photon Conversion; dZ (cm): region", 100, -5, 5);
  BookConversionPlots(histoContainer,"eop_conv_All_region","E/p of Photon Conversion; E/p: region", 100, 0, 3);

  BookConversionPlotsDiPho(histoContainer,"phi_conv_DiPho_region","#phi of Photon Conversion; #phi of Conversion: region", 64, -3.2, 3.2);
  BookConversionPlotsDiPho(histoContainer,"eta_conv_DiPho_region","#eta of Photon Conversion; #eta of Conversion: region", 60, -3.0, 3.0);
  BookConversionPlotsDiPho(histoContainer,"pt_conv_DiPho_region","pt of Photon Conversion; pt of Conversion: region", 200, 0, 200);
  BookConversionPlotsDiPho(histoContainer,"z_conv_DiPho_region","z of Photon Conversion; z (cm): region", 200, -100, 100);
  BookConversionPlotsDiPho(histoContainer,"r_conv_DiPho_region","r of Photon Conversion; r (cm): region", 100, 0, 100);
  BookConversionPlotsDiPho(histoContainer,"dz_conv_DiPho_region","dZ of Photon Conversion; dZ (cm): region", 100, -5, 5);
  BookConversionPlotsDiPho(histoContainer,"eop_conv_DiPho_region","E/p of Photon Conversion; E/p: region", 100, 0, 3);

  // diphoton system
  BookBarrelAndEndcap(histoContainer,"MomentumResolution","#deltaEnergy between the GenParticle and Measured Energy: region;#DeltaE;Counts",100,-0.2,0.2);
  BookBarrelAndEndcap(histoContainer,"MassResolution","#deltaInvMass between the GenParticle Invariant Mass and Measured Invariant Mass: region;#DeltaE;Counts",100,-10,10);
  BookBarrelAndEndcap(histoContainer,"MassResolutionSim","#deltaInvMass between the GenParticle Invariant Mass and Measured Invariant Mass using Sim Vertex: region;#DeltaE;Counts",100,-10,10);
  BookBarrelAndEndcap(histoContainer,"MassResolutionGen","#deltaInvMass between the GenParticle Invariant Mass and Measured Invariant Mass using Sim Vertex: region;#DeltaE;Counts",80,80,160);
  BookMassPlots(histoContainer,"2gamma");
  BookMassPlots(histoContainer,"2gammaGolden");
  BookMassPlots(histoContainer,"2gamma1goodconv");
  BookRCutsMassPlots(histoContainer,"2gamma1goodconv");
  BookMassPlots(histoContainer,"2gamma1goodconvnewvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvlinearvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvsimvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvnearvertex");
  //BookMassPlots(histoContainer,"2gamma1poorconv");
  //BookMassPlots(histoContainer,"2gamma2conv");
  BookMassPlots(histoContainer,"2gammaleftover");
    
  BookFourHists(histoContainer,"lead_r9_catdummy_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
  BookFourHists(histoContainer,"lead_r9_catdummy_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
  BookFourHists(histoContainer,"lead_r9_catdummy_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_allEcal","subleading photon R9, selected candidates: all ECAL",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_Barrel","subleading photon R9, selected candidates: Barrel",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_Endcap","subleading photon R9, selected candidates: Endcap",100,0.,1.1);

  histoContainer->Add("convVtxRvsZBarrelAll"," Photon conversion vtx position all candidates Barrel",400, -200., 200., 200, 0., 80.);
  histoContainer->Add("convVtxRvsZBarrelSel"," Photon conversion vtx position selected candidates Barrel",400, -200., 200., 200, 0., 80.);
  histoContainer->Add("selconvVtxRvsZBarrelSel"," Photon conversion vtx position selected conversion Barrel",400, -200., 200., 200, 0., 80.);
  
  histoContainer->Add("GenPhotonEta","Eta of Gen Daughter Photons;#eta;Counts;",120,-6,6);
  histoContainer->Add("GenPhotonEtaBool","Is the daughter of the Higgs in the Fiducial Region;Bool Value;Counts;",2,0,2);
  histoContainer->Add("GenHiggsEta","Eta of Gen Higgs Boson;",120,-6,6);
  histoContainer->Add("GenHiggsPt","Pt of Gen Daughter Photons;#eta;Counts;",200,0,200);

  histoContainer->Add("higgsPtWhenAtLeastOneGoodConv","Higgs pt when at least one good conversion;#Pt;Counts;",50,0,200);
  histoContainer->Add("higgsPtWhenAtLeastOneGoodConvZInOneSigma","Higgs pt when at least one good conversion;#Pt;Counts;",50,0,200);

}

void BookMassPlots(HistoContainer *histoContainer, TString histname) {

  vector<pair<TString, TString> > parameter;
  parameter.push_back(pair<TString,TString> ("mass_","M_{#gamma#gamma}"));
  parameter.push_back(pair<TString,TString> ("masssimvertex_","M_{#gamma#gamma} SimVertex"));
  parameter.push_back(pair<TString,TString> ("pt_","PT_{2#gamma}"));
  parameter.push_back(pair<TString,TString> ("pz_","Pz_{2#gamma}"));
  parameter.push_back(pair<TString,TString> ("eta_","#eta(2#gamma)"));
  parameter.push_back(pair<TString,TString> ("phi_","#hpi(2#gamma)"));
  parameter.push_back(pair<TString,TString> ("CosThetaStar_","cos#theta^{*}"));

  vector<pair<TString, TString> > region;
  region.push_back(pair<TString,TString> ("Barrel"," barrel "));
  region.push_back(pair<TString,TString> ("Endcap"," endcap "));
    
  vector<pair<TString, TString> > selection;
  selection.push_back(pair<TString,TString> ("All","all"));
  selection.push_back(pair<TString,TString> ("Sel","selected"));
  selection.push_back(pair<TString,TString> ("Matched","matched"));

  int bins[7] = {80, 80, 200, 100, 160, 64, 60};
  float lowerlimit[7] = {80.0, 80.0, 0.0, -1000.0, -8.0, -3.2, 0.0};
  float upperlimit[7] = {160.0, 160.0, 200.0, 1000.0, 8.0, 3.2, 1.0};
    
  for (unsigned int i=0; i<selection.size(); i++) {
    for (unsigned int j=0; j<region.size(); j++) {
      for (unsigned int k=0; k<parameter.size(); k++) {
        TString histnametemp(parameter[k].first);
        histnametemp += histname;
        histnametemp += selection[i].first;
        histnametemp += region[j].first;

        TString histtitletemp = "Di-photon ";
        histtitletemp += parameter[k].second;
        if (histname.Contains("newvertex")) histtitletemp += " with New Vertex";
        if (histname.Contains("linearvertex")) histtitletemp += " with New Linear Fit Vertex";
        if (histname.Contains("simvertex")) histtitletemp += " with Sim Vertex";
        if (histname.Contains("nearvertex") && !histname.Contains("linearvertex")) histtitletemp += " with Near Vertex";
        histtitletemp += ";";
        histtitletemp += parameter[k].second;
        if (k==1 || k==2) {
          histtitletemp += " (GeV) ";
        } else {
          histtitletemp += " ";
        }
        histtitletemp += selection[i].second;
        histtitletemp += region[j].second;
        histtitletemp += "candidates";
          
        if (!(histname.Contains("vertex") && parameter[k].first.Contains("simvertex"))) histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins[k],lowerlimit[k],upperlimit[k]);
      }
    }
  }
    
}

void BookPhysicsPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {

  vector <TString> region;
  region.push_back("allEcal");
  region.push_back("Barrel");
  region.push_back("Endcap");

  for (unsigned int i=0; i<2; i++){
    for (unsigned int j=0; j<2; j++) {
      for (unsigned int k=0; k<region.size(); k++) {
        TString histnametemp(histname);
        TString histtitletemp(histtitle);
        histnametemp.ReplaceAll("region",region[k]);
        histtitletemp.ReplaceAll("region",region[k]);
        if (i==1) {
          histnametemp.ReplaceAll("All_","Sel_");
          histnametemp.ReplaceAll(" all "," selected ");
        }
        if (j==1) {
          histnametemp.ReplaceAll("lead","sublead");
          histtitletemp.ReplaceAll("lead","sublead");
        }
        histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
      }
    }
  }

}

void BookRCutsMassPlots(HistoContainer *histoContainer, TString histname) {

  vector<pair<TString, TString> > region;
  region.push_back(pair<TString,TString> ("Barrel","barrel "));
  region.push_back(pair<TString,TString> ("Endcap","endcap "));
    
  vector<pair<TString, TString> > selection;
  selection.push_back(pair<TString,TString> ("All","all "));
  selection.push_back(pair<TString,TString> ("Sel","selected "));
  selection.push_back(pair<TString,TString> ("Matched","matched "));

  int RCut[8] = {20, 30, 40, 50, 60, 70, 80, 90};
  int bins = 80;
  float lowerlimit = 80.0;
  float upperlimit = 160.0;
    
  for (unsigned int i=0; i<selection.size(); i++) {
    for (unsigned int j=0; j<region.size(); j++) {
      for (unsigned int k=0; k<8; k++) {
        TString histnametemp("mass_");
        histnametemp += histname;
        histnametemp += selection[i].first;
        histnametemp += RCut[k];
        histnametemp += region[j].first;

        TString histtitletemp = "Di-photon invariant mass with an R cut of ";
        histtitletemp += RCut[k];
        histtitletemp += "cm;M_{#gamma#gamma} (GeV) ";
        histtitletemp += selection[i].second;
        histtitletemp += region[j].second;
        histtitletemp += "candidates with one conversion";
          
        histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
      }
    }
  }
    
}

void BookRCutsdZPlots(HistoContainer *histoContainer, TString histname) {

  vector<pair<TString, TString> > region;
  region.push_back(pair<TString,TString> ("Barrel","barrel "));
  region.push_back(pair<TString,TString> ("Endcap","endcap "));
    
  vector<pair<TString, TString> > selection;
  selection.push_back(pair<TString,TString> ("All","all "));
  selection.push_back(pair<TString,TString> ("Sel","selected "));

  int RCut[8] = {20, 30, 40, 50, 60, 70, 80, 90};
  int bins = 100;
  float lowerlimit = -5;
  float upperlimit = 5;
    
  for (unsigned int i=0; i<selection.size(); i++) {
    for (unsigned int j=0; j<region.size(); j++) {
      for (unsigned int k=0; k<8; k++) {
        TString histnametemp = histname;
        histnametemp += selection[i].first;
        histnametemp += RCut[k];
        histnametemp += region[j].first;

        TString histtitletemp = "#deltaZ with an R cut of ";
        histtitletemp += RCut[k];
        histtitletemp += "cm;#delta Z of ";
        histtitletemp += selection[i].second;
        histtitletemp += region[j].second;
        histtitletemp += "candidates with one conversion";
          
        histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
      }
    }
  }
    
}

void BookTrackerdZ(HistoContainer *histoContainer, TString histname) {

  int bins = 100;
  int lowerlimit = -5;
  int upperlimit = 5;

  TString histtitle = "#DeltaZ between the Z of the Primary Vertex using method and the SimVertex: region;#DeltaZ (cm);";
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
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
  }

}

void FillCatHists(sdaReader *currentTree, HistoContainer *histoContainer, float weight, unsigned int leadindex, unsigned int subleadindex, int leadPhoCategory, int subleadPhoCategory) {

  for (int i=0; i<4; i++) {
    for (unsigned int j=0; j<2; j++) {
      TString histname = "lead_r9_cat";
      histname += i;
      string region = "";
      if (j==0) region = "_allEcal";
      if (j==1 && currentTree->pho_isEB[leadindex]) region = "_Barrel";
      if (j==1 && currentTree->pho_isEE[leadindex]) region = "_Endcap";
      histname += region;
      
      if (leadPhoCategory==i) histoContainer->Fill(histname.Data(),currentTree->pho_r9[leadindex],weight);
      histname.Prepend("sub");
      if (j==1 && currentTree->pho_isEB[subleadindex]) region = "_Barrel";
      if (j==1 && currentTree->pho_isEE[subleadindex]) region = "_Endcap";
      if (subleadPhoCategory==i) histoContainer->Fill(histname.Data(),currentTree->pho_r9[subleadindex],weight);
      
    }
  }
  
}

void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, unsigned int photonconvindex, int convindex, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop) {

  for (unsigned int j=0; j<2; j++) {
    TString histnames[] = {"phi_conv_", "eta_conv_", "pt_conv_", "z_conv_", "r_conv_", "dz_conv_","eop_conv_"};
    string region = "";
    if (j==0) region = "_allEcal";
    if (j==1 && currentTree->pho_isEB[photonconvindex]) region = "_Barrel";
    if (j==1 && currentTree->pho_isEE[photonconvindex]) region = "_Endcap";

    for (unsigned int k=0; k<7; k++) {
      histnames[k] += selection;
      histnames[k] += region;
    }

    histoContainer->Fill(histnames[0].Data(),Photonxyz[photonconvindex].Phi(), weight);
    histoContainer->Fill(histnames[1].Data(),Photonxyz[photonconvindex].Eta(), weight);
    histoContainer->Fill(histnames[2].Data(),Photonp4[photonconvindex].Pt(), weight);
    histoContainer->Fill(histnames[3].Data(),ConversionVertex[convindex].z(), weight);
    histoContainer->Fill(histnames[4].Data(),ConversionVertex[convindex].Perp(), weight);
    histoContainer->Fill(histnames[5].Data(),zPVFromConv - SimVertex.Z(), weight);
    histoContainer->Fill(histnames[6].Data(),eop, weight);

  }

}

void FillConvHistsInDiPho(sdaReader *currentTree, HistoContainer *histoContainer, unsigned int photonconvindex, int convindex, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop) {

  for (unsigned int j=0; j<2; j++) {
    TString histnames[] = {"phi_conv_", "eta_conv_", "pt_conv_", "z_conv_", "r_conv_", "dz_conv_","eop_conv_"};
    string region = "";
    if (j==0) region = "_allEcal";
    if (j==1 && currentTree->pho_isEB[photonconvindex]) region = "_Barrel";
    if (j==1 && currentTree->pho_isEE[photonconvindex]) region = "_Endcap";

    for (unsigned int k=0; k<7; k++) {
      histnames[k] += selection;
      histnames[k] += region;
    }

    histoContainer->Fill(histnames[0].Data(),Photonxyz[photonconvindex].Phi(), weight);
    histoContainer->Fill(histnames[1].Data(),Photonxyz[photonconvindex].Eta(), weight);
    histoContainer->Fill(histnames[2].Data(),Photonp4[photonconvindex].Pt(), weight);
    histoContainer->Fill(histnames[3].Data(),ConversionVertex[convindex].z(), weight);
    histoContainer->Fill(histnames[4].Data(),ConversionVertex[convindex].Perp(), weight);
    histoContainer->Fill(histnames[5].Data(),zPVFromConv - SimVertex.Z(), weight);
    histoContainer->Fill(histnames[6].Data(),eop, weight);

  }

}

void FilldZRcut(HistoContainer *histoContainer, string selection, float weight, string iConvDetector, double R, double deltaz) {

  float RCuts[] = {90, 80, 70, 60, 50, 40 ,30, 20};
  for (unsigned int i=0; i<8; i++) {
    TString histname = "dz_conv_";
    histname += selection;
    histname += (int) floor(RCuts[i]);
    histname += iConvDetector;
    if (R<RCuts[i]) histoContainer->Fill(histname.Data(),deltaz,weight);
  }

}

void FilldZTrackerBarrel(HistoContainer *histoContainer, TString histname, double deltaZ, double R, float weight) {

  TString histnametemp = "";
  if (R<=15.0) histnametemp = "PixelBarrel"+histname;
  else if (R>15 && R<=60.0) histnametemp = "TIB"+histname;
  else if (R>60.0) histnametemp = "TOB"+histname;
  histoContainer->Fill(histnametemp.Data(),deltaZ,weight);
  
}

void FilldZTrackerEndcap(HistoContainer *histoContainer, TString histname, double deltaZ, double Z, float weight) {

  TString histnametemp = "";
  if (Z<=50.0) histnametemp = "PixelFwd"+histname;
  else if (Z>50 && Z<=100.0) histnametemp = "TID"+histname;
  else if (Z>100.0) histnametemp = "TEC"+histname;
  histoContainer->Fill(histnametemp.Data(),deltaZ,weight);

}

void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double SimInvMass, double cos_thetastar) {

  TString histnames[] = {"mass_", "masssimvertex_", "pt_", "pz_", "eta_", "phi_", "CosThetaStar_"};
  unsigned int arraylength = (unsigned int) sizeof(histnames)/sizeof(histnames[0]);
  for (unsigned int i=0; i<arraylength; i++) {
    histnames[i] += "2gamma";
    histnames[i] += selection;
    histnames[i] += HiggsInWhichDetector;
  }

  if ( InvMass > 100 ) { // get rid of the Z
    histoContainer->Fill(histnames[0].Data(),InvMass,weight);
    histoContainer->Fill(histnames[1].Data(),SimInvMass,weight);
    histoContainer->Fill(histnames[2].Data(),VSum.Pt(),weight);
    histoContainer->Fill(histnames[3].Data(),VSum.Pz(),weight);
    histoContainer->Fill(histnames[4].Data(),VSum.Eta(),weight);
    histoContainer->Fill(histnames[5].Data(),VSum.Phi(),weight);
    histoContainer->Fill(histnames[6].Data(),cos_thetastar,weight);
    
    if (diPhoCategory==1) for (unsigned int i=0; i<arraylength; i++) histnames[i].ReplaceAll("2gamma","2gammaGolden");
    if (diPhoCategory==2) for (unsigned int i=0; i<arraylength; i++) histnames[i].ReplaceAll("2gamma","2gammaleftover");
    //if (diPhoCategory==3) for (unsigned int i=0; i<arraylength; i++) histnames[i].ReplaceAll("2gamma","2gamma1poorconv");
    //if (diPhoCategory==4) for (unsigned int i=0; i<arraylength; i++) histnames[i].ReplaceAll("2gamma","2gamma2conv");

    if (diPhoCategory==1 || diPhoCategory==2) {
      histoContainer->Fill(histnames[0].Data(),InvMass,weight);
      histoContainer->Fill(histnames[1].Data(),SimInvMass,weight);
      histoContainer->Fill(histnames[2].Data(),VSum.Pt(),weight);
      histoContainer->Fill(histnames[3].Data(),VSum.Pz(),weight);
      histoContainer->Fill(histnames[4].Data(),VSum.Eta(),weight);
      histoContainer->Fill(histnames[5].Data(),VSum.Phi(),weight);
      histoContainer->Fill(histnames[6].Data(),cos_thetastar,weight);
    }

  }
          
}

void FillMassNewVertexHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, string vertextype, TLorentzVector VSum, double InvMass, double cos_thetastar) {

  TString histnames[] = {"mass_", "pt_", "pz_", "eta_", "phi_", "CosThetaStar_"};
  for (unsigned int i=0; i<6; i++) {
    histnames[i] += "2gamma1goodconv";
    histnames[i] += vertextype;
    histnames[i] += selection;
    histnames[i] += HiggsInWhichDetector;
  }
  
  histoContainer->Fill(histnames[0].Data(),InvMass,weight);
  histoContainer->Fill(histnames[1].Data(),VSum.Pt(),weight);
  histoContainer->Fill(histnames[2].Data(),VSum.Pz(),weight);
  histoContainer->Fill(histnames[3].Data(),VSum.Eta(),weight);
  histoContainer->Fill(histnames[4].Data(),VSum.Phi(),weight);
  histoContainer->Fill(histnames[5].Data(),cos_thetastar,weight);

}

void FillMassRcut(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, double R, double InvMass) {

  float RCuts[] = {90, 80, 70, 60, 50, 40 ,30, 20};
  for (unsigned int i=0; i<8; i++) {
    TString histname = "mass_2gamma1goodconv";
    histname += selection;
    histname += (int) floor(RCuts[i]);
    histname += HiggsInWhichDetector;
    if (R<RCuts[i]) histoContainer->Fill(histname.Data(),InvMass,weight);
  }

}

void FillPhotonHists(sdaReader *currentTree, HistoContainer *histoContainer, string photon, unsigned int index, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, TVector3 PrimaryVertex, TVector3 SimVertex, vector<TLorentzVector> Photonp4, bool data) {

  for (unsigned int i=0; i<2; i++) {
    TString histnames[] = {"PhoEt", "PhoEta", "PhoPhi", "PhoR9", "PhoHoE", "PhoTrkPtSumHollow04", "PhoEcalPtSumHollow04", "PhoHcalPtSumHollow04", "PhoSigmaIetaIeta", "PhoZPV", "PhoDZPV", "PhoDZPVconv"};
    string region = "";
    if (i==0) region = "_allEcal";
    if (i==1 && currentTree->pho_isEB[index]) region = "_Barrel";
    if (i==1 && currentTree->pho_isEE[index]) region = "_Endcap";
      
    for (unsigned int k=0; k<12; k++) {
      histnames[k].Prepend(photon);
      histnames[k].Append(selection);
      histnames[k].Append(region);
    }
    
    histoContainer->Fill(histnames[0].Data(),Photonp4[index].Et(),weight);
    histoContainer->Fill(histnames[1].Data(),Photonxyz[index].Eta(),weight);
    histoContainer->Fill(histnames[2].Data(),Photonxyz[index].Phi(),weight);
    histoContainer->Fill(histnames[3].Data(),currentTree->pho_r9[index],weight);
    histoContainer->Fill(histnames[4].Data(),currentTree->pho_hoe[index],weight);
    histoContainer->Fill(histnames[5].Data(),currentTree->pho_trksumpthollowconedr04[index],weight);
    histoContainer->Fill(histnames[6].Data(),currentTree->pho_ecalsumetconedr04[index],weight);
    histoContainer->Fill(histnames[7].Data(),currentTree->pho_hcalsumetconedr04[index],weight);
    histoContainer->Fill(histnames[8].Data(),currentTree->pho_sieie[index],weight);
    histoContainer->Fill(histnames[9].Data(),PrimaryVertex.Z(),weight);
    if (!data) {
      histoContainer->Fill(histnames[10].Data(),PrimaryVertex.Z()-SimVertex.Z(),weight);
      histoContainer->Fill(histnames[11].Data(),currentTree->pho_conv_zofprimvtxfromtrks[index]-SimVertex.Z(),weight);
    }
  }

}

void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData) {

  float BranchingFraction = 0;
  float kFactor = 1.0;
  
  if (inputstring.Contains("Data") && !inputstring.Contains("Dataweight")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/paoloReRecoWithLaserV5ESFixed_v1.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/paoloReRecoWithLaserV5_v2.root",1));
    isData=true;
  }
  if (inputstring.Contains("PromptReco") && !inputstring.Contains("PromptRecoweight")) {
    inputfilelist.push_back(pair<string,int> ("PromptReco.root",7));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/Run2010B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_160404_161474.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_162000_163100.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_163100_163400.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_163400_163700.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/promptreco_163700_164000.root",1));
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
    kFactor=1.2;
    inputfilelist.push_back(pair<string,int> ("PhotonJetEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GJet_Pt20_2011PU.root",kFactor*77100/(1182075/0.0064)));
  }
  if (inputstring.Contains("QCD") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDDoubleEMEnriched.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/QCD_Pt30to40_2011PU.root",41800000/(3550408/0.00023)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/QCD_Pt40_2011PU.root",18700000/(21276029/0.00216)));
  }
  if (inputstring.Contains("Born") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    kFactor = 1.15;
    inputfilelist.push_back(pair<string,int> ("Born.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt10to25_2011PU.root",kFactor*236.4/522865));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt25to250_2011PU.root",kFactor*22.37/537445));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBorn_Pt250toinf_2011PU.root",kFactor*0.008072/546355));
  }
  if (inputstring.Contains("Box") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    kFactor = 1.15;
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt10to25_2011PU.root",kFactor*358.2/797975));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt25to250_2011PU.root",kFactor*12.37/777725));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DiPhotonBox_Pt250toinf_2011PU.root",kFactor*0.000208/789470));
  }
  if (inputstring.Contains("DY") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("DrellYan.root",2));
    kFactor = 1.15;
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DYToEE_M10To20_2011PU.root",kFactor*2659.0/1933000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/DYToEE_powheg_M20_2011PU.root",kFactor*1614.0/1998990));
  }
  if (inputstring.Contains("Test")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/CMSSW414/GluGlu_M120_2011PU.root",1));
  }

}

void PrintWeights() {

  //Print Weights
  const int nMassPoints=3;
  const int nProdProcess=3;
  double br[nMassPoints];
  br[0]= 0.001939;
  br[1]= 0.002219;
  br[2]= 0.001363;
  double SignalXSec[nMassPoints][nProdProcess];
  double SignalW[nMassPoints][nProdProcess];
  SignalXSec[0][0]=20.493;
  SignalXSec[0][1]=1.4405;
  SignalXSec[0][2]=1.4421;
  SignalXSec[1][0]=17.173;
  SignalXSec[1][1]=1.3062;
  SignalXSec[1][2]=1.0921;
  SignalXSec[2][0]=10.863;
  SignalXSec[2][1]=0.9868;
  SignalXSec[2][2]=0.5515;

  for ( int j=0; j<nProdProcess; j++ ) {
    for ( int i=0; i<nMassPoints; i++ ) {
      SignalW[i][j] = br[i]*SignalXSec[i][j];
      cout << " Process " << j << " Mass Point "<< i << " Weight " <<
        SignalW[i][j] << endl;
    }
  }

}
