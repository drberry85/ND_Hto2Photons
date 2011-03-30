#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"
#include "TString.h"

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
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

bool DoGenMatching(sdaReader *currentTree, TVector3 LeadPhoton, TVector3 SubleadPhoton);
bool DoVertexMatching(TVector3 TestVertex, TVector3 SimVertex);
bool sortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, double &subleadpt, unsigned int &leadindex, unsigned int &subleadindex);
unsigned int getconvindex(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
double CosThetaStar(TLorentzVector VLead, TLorentzVector VSum);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
double FindRefittedZ(TVector3 ConversionVertex, TVector3 ConversionRefittedPairMomentum);
string HiggsDetectorPosition(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
string DetectorPosition(sdaReader *currentTree, unsigned int index);
TString MakeFileName(string filename, bool unweighted, bool dataweight, double RCut);
TLorentzVector CalcDiffVertex(TLorentzVector p4, TVector3 xyz, float newz);
void BookConversionPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void BookFourHists(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookMassPlots(HistoContainer *histoContainer, TString histname);
void BookPhysicsPlots(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookRCutsMassPlots(HistoContainer *histoContainer, TString histname);
void BookRCutsdZPlots(HistoContainer *histoContainer, TString histname);
void FillCatHists(sdaReader *currentTree, HistoContainer *histoContainer, float weight, unsigned int leadindex, unsigned int subleadindex, int leadPhoCategory, int subleadPhoCategory);
void FillPhotonHists(sdaReader *currentTree, HistoContainer *histoContainer, string photon, unsigned int index, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, TVector3 PrimaryVertex, TVector3 SimVertex, vector<TLorentzVector> Photonp4, bool data);
void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, int iPho, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop);
void FillConvHistsInDiPho(sdaReader *currentTree, HistoContainer *histoContainer, int iPho, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop);
void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double cos_thetastar, double R);
void FillMassNewVertexHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, string vertextype, TLorentzVector VSum, double InvMass, double cos_thetastar);
void FillMassRcut(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, double R, double InvMass);
void FilldZRcut(HistoContainer *histoContainer, string selection, float weight, string iConvDetector, double R, double deltaz);

void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData);
void ProgressBar(int &percent);
void PrintWeights();

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool dataweight = false;
  bool unweighted = false;
  double RCut = 999999;
  float globalWeight = 5000;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;
    
  TString InputArgs(input[1]);

  //PrintWeights();
  MakeFilesAndWeights(InputArgs, filesAndWeights, filelist, data);
  
  if (InputArgs.Contains("Dataweight")) {
    globalWeight = 36.0;
    dataweight=true;
  }
  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("RCut")) RCut=40;
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: 90GeV, 110GeV, 120GeV, 150GeV, PhotonPlusJet, EMEnriched, DoubleEMEnriched, QCDBEtoE, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = MakeFileName(itFilePair->first, unweighted, dataweight, RCut);
    
    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << "\n" << outfilename << " created." << endl;

    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      bool HasSimVertex = false;
      bool HasGenParticles = false;
      string file = filesAndWeights[itFile].first;
      float weight = filesAndWeights[itFile].second * globalWeight;

      if (unweighted) weight=1;
      if (itFilePair->first=="Data.root" || itFilePair->first=="YousiData.root") weight=1;
      
      TFile * currentFile = new TFile(file.c_str());
      currentFile->cd();

      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      sdaReader currentTree(currentFile);

      Long64_t nentries = currentTree.fChain->GetEntries();
      cout << "TreeRead Entries " << nentries << endl;     

      outfile->cd();

      int percent = 0;
      
      for ( Long64_t i = 0; i < nentries; i++ ) {

        if (i % (nentries/100) == 0 && bar) {
          if (percent == 100) percent = 99;
          ProgressBar(percent);
        }

        currentTree.GetEntry(i);

        if (i==0) HasSimVertex = currentTree.FindLeaf("simvtx");
        if (i==0) HasGenParticles = currentTree.FindLeaf("gp_vtx");
        
        TVector3 SimVertex(0,0,0);
        if (!data && HasSimVertex) SimVertex = (*((TVector3*) currentTree.simvtx->At(0)));
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionPairMomentum;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TLorentzVector> Photonp4;
        vector<TLorentzVector> SuperClusterp4;

        //Finding Vertex of Higgs
        if (!data && HasGenParticles) {
          for (int j=0; j<currentTree.gp_n; j++) {
            if (currentTree.gp_pdgid[j]==25 && currentTree.gp_status[j]==2) {
              SimVertex = (*((TVector3*) currentTree.gp_vtx->At(j)));
            }
          }
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

        map <double, unsigned int> ptindex;
        for (unsigned int j=0; j<(unsigned int) currentTree.pho_n; j++) {
          ConversionVertex.push_back(*((TVector3*) currentTree.pho_conv_vtx->At(j)));
          ConversionPairMomentum.push_back(*((TVector3*) currentTree.pho_conv_pair_momentum->At(j)));
          Photonxyz.push_back(*((TVector3*) currentTree.pho_calopos->At(j)));
          ConversionRefittedPairMomentum.push_back(*((TVector3*) currentTree.pho_conv_refitted_momentum->At(j)));
          Photonp4.push_back(*((TLorentzVector*) currentTree.pho_p4->At(j)));
          ptindex[Photonp4[j].Pt()]=j;
        }
        for (unsigned int j=0; j<(unsigned int) currentTree.sc_p4->GetSize(); j++) SuperClusterp4.push_back(*((TLorentzVector*) currentTree.sc_p4->At(j)));

        double leadpt = -1;
        double subleadpt = -1;
        unsigned int leadindex = 0;
        unsigned int subleadindex = 0;

        bool sorted = sortpt(ptindex, currentTree.pho_n, leadpt, subleadpt, leadindex, subleadindex);
        
        if (!sorted) cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        //cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        histoContainer->Fill("NumberVertices",float(currentTree.vtx_std_xyz->GetSize()),weight);
        if (HasSimVertex) histoContainer->Fill("NumberSimVerticesVsRecoVertices",float(currentTree.vtx_std_xyz->GetSize()), float(currentTree.simvtx->GetSize()));
        if (!data && HasSimVertex) histoContainer->Fill("NumberSimVertices",float(currentTree.simvtx->GetSize()),weight);

        //////////////// basic selection
        if (!preselection(Photonp4[leadindex].Pt(), Photonp4[subleadindex].Pt(), Photonxyz[leadindex].Eta(), Photonxyz[subleadindex].Eta(), currentTree.pho_isEBEEGap[leadindex], currentTree.pho_isEBEEGap[subleadindex])) continue;

        ////////////////////////////////////
        unsigned int convindex = getconvindex(&currentTree,leadindex,subleadindex);
        unsigned int convscindex = (unsigned int) currentTree.pho_scind[convindex];
        unsigned int nearvertexindex = 0;
        float scEnergy = Photonp4[convindex].E();
        if (currentTree.sc_p4->GetSize()>0) scEnergy = SuperClusterp4[convscindex].E();
        float EoP = scEnergy/ConversionRefittedPairMomentum[convindex].Mag();
        
        string iLeadDetector = DetectorPosition(&currentTree, leadindex);
        string iSubleadDetector = DetectorPosition(&currentTree, subleadindex);
        string iConvDetector = DetectorPosition(&currentTree, convindex);
        string selection = "All";

        int leadPhoCategory = photonCategory ( currentTree.pho_haspixseed[leadindex], currentTree.pho_r9[leadindex],  currentTree.pho_conv_ntracks[leadindex], currentTree.pho_conv_chi2_probability[leadindex] , Photonp4[leadindex].Pt()/ConversionPairMomentum[leadindex].Perp(), ConversionVertex[leadindex].Perp());
        int subleadPhoCategory = photonCategory (currentTree.pho_haspixseed[subleadindex], currentTree.pho_r9[subleadindex],  currentTree.pho_conv_ntracks[subleadindex], currentTree.pho_conv_chi2_probability[subleadindex] ,  Photonp4[subleadindex].Pt()/ConversionPairMomentum[subleadindex].Perp(), ConversionVertex[subleadindex].Perp());


        bool convsel1 = convSel(currentTree.pho_conv_ntracks[leadindex],
                                currentTree.pho_conv_validvtx[leadindex] ,  
                                currentTree.pho_conv_chi2_probability[leadindex], 
                                currentTree.pho_conv_dphitrksatvtx[leadindex], 
                                currentTree.pho_conv_paircotthetasep[leadindex], 
                                EoP,
                                ConversionVertex[leadindex].Perp());
	
        bool convsel2 = convSel(currentTree.pho_conv_ntracks[subleadindex],
                                currentTree.pho_conv_validvtx[subleadindex] ,  
                                currentTree.pho_conv_chi2_probability[subleadindex], 
                                currentTree.pho_conv_dphitrksatvtx[subleadindex], 
                                currentTree.pho_conv_paircotthetasep[subleadindex], 
                                EoP,
                                ConversionVertex[subleadindex].Perp());

        int diPhoCategory = diPhotonCategory( leadPhoCategory, subleadPhoCategory );
        if (leadPhoCategory==2 && ConversionVertex[convindex].Perp()>RCut) leadPhoCategory=3;
        if (subleadPhoCategory==2 && ConversionVertex[convindex].Perp()>RCut) subleadPhoCategory=3;

        double NewZconvlinear = 0;
        ////////////////////////////////////

        if (leadPhoCategory==2 || subleadPhoCategory==2) {

          if (!data) FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
          histoContainer->Fill("convdZvsEtaAll",Photonxyz[convindex].Eta(),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          histoContainer->Fill("convEoPAll",iConvDetector,SuperClusterp4[convscindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          histoContainer->Fill("convr",iConvDetector,ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("Zconv",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex],weight);
          if (!data) {
            histoContainer->Fill("convdZPrimary",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),PrimaryVertex[0].Z()-SimVertex.Z(),weight);
            histoContainer->Fill("ZconvdZ",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            histoContainer->Fill("ZconvdZVsNumOfVertices",iConvDetector,float(currentTree.vtx_std_xyz->GetSize()),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
	
            if ( fabs(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < 1 ) histoContainer->Fill("ZconvdZLessOnecm",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            if ( fabs(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < .5 ) histoContainer->Fill("ZconvdZLessHalfcm",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            if ( fabs(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z() ) < .3 ) histoContainer->Fill("ZconvdZLess3mm",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
	    
          }

          TVector3 myBeamSpot(0,0,0);
          double NewZ = FindNewdZ(Photonxyz[convindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          double RefittedZ = FindRefittedZ(ConversionVertex[convindex], ConversionRefittedPairMomentum[convindex]);
          NewZconvlinear = FindNewZConvLinear(ConversionVertex[convindex],Photonxyz[convindex],PrimaryVertex[0]);

          float deltaz = 100000;
          for (unsigned int j=0; j<PrimaryVertex.size(); j++) {
            if (deltaz > fabs(NewZconvlinear-PrimaryVertex[j].Z())) {
              deltaz = fabs(NewZconvlinear-PrimaryVertex[j].Z());
              nearvertexindex = j;
            }
          }

          histoContainer->Fill("RefittedZconv",iConvDetector, RefittedZ, weight);
          histoContainer->Fill("NewZconv",iConvDetector, NewZ, weight);
          myBeamSpot.SetXYZ(PrimaryVertex[0].X(),PrimaryVertex[0].Y(),0);
          NewZ = FindNewdZ(Photonxyz[convindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          histoContainer->Fill("NewZPVconv",iConvDetector, NewZ, weight);
          histoContainer->Fill("NewZconvlinear",iConvDetector, NewZconvlinear, weight);
          
          if (!data) {
            histoContainer->Fill("ZconvdZNearest",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-PrimaryVertex[nearvertexindex].Z(),weight);
            histoContainer->Fill("ZdZNear",iConvDetector,PrimaryVertex[nearvertexindex].Z()-SimVertex.Z(),weight);
            histoContainer->Fill("RefittedZconvdZ",iConvDetector, RefittedZ-SimVertex.Z(), weight);
            histoContainer->Fill("NewZconvdZ",iConvDetector, NewZ-SimVertex.Z(), weight);
            histoContainer->Fill("NewZPVconvdZ",iConvDetector, NewZ-SimVertex.Z(), weight);
            histoContainer->Fill("NewZconvlineardZ",iConvDetector, NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<1) histoContainer->Fill("ZconvdZLinearLessOnecm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<0.5) histoContainer->Fill("ZconvdZLinearLessHalfcm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
            if (abs(NewZconvlinear-SimVertex.Z())<0.3) histoContainer->Fill("ZconvdZLinearLess3mm",iConvDetector,NewZconvlinear-SimVertex.Z(), weight);
          }
        }

        TVector3 NearVertex = PrimaryVertex[nearvertexindex];
        bool PrimaryVertexMatched = DoVertexMatching(PrimaryVertex[0],SimVertex);
        bool NearVertexMatched = DoVertexMatching(NearVertex,SimVertex);

        histoContainer->Fill("PrimaryVertexMatchedAllEcal", PrimaryVertexMatched, weight);
        histoContainer->Fill("NearVertexMatchedAllEcal", NearVertexMatched, weight);
        histoContainer->Fill("PrimaryVertexMatched", iConvDetector, PrimaryVertexMatched, weight);
        histoContainer->Fill("NearVertexMatched", iConvDetector, NearVertexMatched, weight);

        /*if (Photonp4[subleadindex].Pt()>Photonp4[leadindex].Pt()) {
          cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
          cout << "LeadPt is " << Photonp4[leadindex].Pt() << endl;
          cout << "SubLeadPt is " << Photonp4[subleadindex].Pt() << endl;
          }*/

        histoContainer->Fill("NPhotonsAll",currentTree.pho_n,weight);
        //// Fill histograms hists only for conversions satisfying the conversion selection
        if ( convsel1) FillConvHists(&currentTree, histoContainer, leadindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
        if ( convsel2) FillConvHists(&currentTree, histoContainer, subleadindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
	
        FillPhotonHists(&currentTree, histoContainer, string("lead"), leadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        FillPhotonHists(&currentTree, histoContainer, string("sublead"), subleadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        
        if (currentTree.pho_conv_ntracks[leadindex]==2 && currentTree.pho_conv_chi2_probability[leadindex]>0.0005 && (bool) currentTree.pho_isEB[leadindex])
          histoContainer->Fill("convVtxRvsZBarrelAll",ConversionVertex[leadindex].Z(),ConversionVertex[leadindex].Perp(),weight);
        if (currentTree.pho_conv_ntracks[subleadindex]==2 && currentTree.pho_conv_chi2_probability[subleadindex]>0.0005 && (bool) currentTree.pho_isEB[subleadindex])
          histoContainer->Fill("convVtxRvsZBarrelAll",ConversionVertex[leadindex].Z(),ConversionVertex[leadindex].Perp(),weight);
        
        TLorentzVector VLead(Photonp4[leadindex]);
        TLorentzVector VSubLead(Photonp4[subleadindex]);
        TLorentzVector VSum=VLead+VSubLead;
        double InvMass=fabs(VSum.M());
        double cos_thetastar = CosThetaStar(VLead,VSum);
        
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
        
        if (diPhoCategory==2) {

          VLead_newvertex = CalcDiffVertex(Photonp4[leadindex], Photonxyz[leadindex], currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          VSubLead_newvertex = CalcDiffVertex(Photonp4[subleadindex], Photonxyz[subleadindex], currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
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
        FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, cos_thetastar, ConversionVertex[convindex].Perp());
        if (diPhoCategory==2) {
          FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        }
          
        ///////////////////////////////  Event selection ///////////////////////////////////////
        ///////////// build the di-photon system
        if (Photonp4[leadindex].Pt()<40) continue; // leading photon
        if (Photonp4[subleadindex].Pt()<30) continue; // subleading photon
        if (fabs(Photonxyz[leadindex].Eta())>2.5 || fabs(Photonxyz[subleadindex].Eta())>2.5) continue;

        //isolation 
        if (!(tightId(Photonp4[leadindex].Pt(),
                      currentTree.pho_ecalsumetconedr04[leadindex],
                      currentTree.pho_hcalsumetconedr04[leadindex],
                      currentTree.pho_trksumpthollowconedr04[leadindex],
                      (bool) currentTree.pho_isEB[leadindex],
                      (bool) currentTree.pho_isEE[leadindex],
                      currentTree.pho_sieie[leadindex],
                      currentTree.pho_hoe[leadindex]))) continue;

        if (!(tightId(Photonp4[subleadindex].Pt(),
                      currentTree.pho_ecalsumetconedr04[subleadindex],
                      currentTree.pho_hcalsumetconedr04[subleadindex],
                      currentTree.pho_trksumpthollowconedr04[subleadindex],
                      (bool) currentTree.pho_isEB[subleadindex],
                      (bool) currentTree.pho_isEE[subleadindex],
                      currentTree.pho_sieie[subleadindex],
                      currentTree.pho_hoe[subleadindex]))) continue;

        histoContainer->Fill("NPhotonsSel",currentTree.pho_n,weight);
        selection = "Sel";

        FillPhotonHists(&currentTree, histoContainer, string("lead"), leadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        FillPhotonHists(&currentTree, histoContainer, string("sublead"), subleadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        if ( convsel1 ) 
          FillConvHists(&currentTree, histoContainer, leadindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
        if ( convsel2 ) 
          FillConvHists(&currentTree, histoContainer, subleadindex, selection, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);

        if (diPhoCategory==2) {
          string diPhoSel = "DiPho";	  
          if ( convsel1 ) 
            FillConvHistsInDiPho(&currentTree, histoContainer, leadindex, diPhoSel, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
          if ( convsel2 ) 
            FillConvHistsInDiPho(&currentTree, histoContainer, subleadindex, diPhoSel, weight, Photonxyz, ConversionVertex, Photonp4,currentTree.pho_conv_zofprimvtxfromtrks[convindex],SimVertex, EoP);
	  
          histoContainer->Fill("convdZvsEtaSel",Photonxyz[convindex].Eta(),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          if (!data) {
            histoContainer->Fill("convdZvsR",iConvDetector,ConversionVertex[convindex].Perp(),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            histoContainer->Fill("convdZvsZ",iConvDetector,ConversionVertex[convindex].z(),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
            FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
          }
        }
        
        if (currentTree.pho_conv_ntracks[leadindex]==2 && currentTree.pho_conv_chi2_probability[leadindex]>0.0005 && (bool) currentTree.pho_isEB[leadindex])
          histoContainer->Fill("convVtxRvsZBarrelSel",ConversionVertex[leadindex].Z(),ConversionVertex[leadindex].Perp(),weight);
        if (currentTree.pho_conv_ntracks[subleadindex]==2 && currentTree.pho_conv_chi2_probability[subleadindex]>0.0005 && (bool) currentTree.pho_isEB[subleadindex])
          histoContainer->Fill("convVtxRvsZBarrelSel",ConversionVertex[subleadindex].Z(),ConversionVertex[subleadindex].Perp(),weight);

        if (convsel1) histoContainer->Fill("selconvVtxRvsZBarrelSel",ConversionVertex[leadindex].Z(),sqrt(ConversionVertex[leadindex].X()*ConversionVertex[leadindex].X()+ConversionVertex[leadindex].Y()*ConversionVertex[leadindex].Y()));
        if (convsel2) histoContainer->Fill("selconvVtxRvsZBarrelSel",ConversionVertex[subleadindex].Z(),sqrt(ConversionVertex[subleadindex].X()*ConversionVertex[subleadindex].X()+ConversionVertex[subleadindex].Y()*ConversionVertex[subleadindex].Y()));

        FillCatHists(&currentTree, histoContainer, weight, leadindex, subleadindex, leadPhoCategory, subleadPhoCategory);

        FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, cos_thetastar, ConversionVertex[convindex].Perp());
        if (diPhoCategory==2) {
          FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        }
        
        //GenMatching
//         bool GenMatched = false;
//         if (HasGenParticles) GenMatched = DoGenMatching(&currentTree, Photonxyz[leadindex], Photonxyz[subleadindex]);
//         if (!data && GenMatched) {
//           selection="Matched";
//           FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, cos_thetastar, ConversionVertex[convindex].Perp());
//           if (diPhoCategory==2) {
//             FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
//             FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
//             FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "linearvertex", VSum_linearvertex, InvMass_linearvertex, cos_thetastar_linearvertex);
//             if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
//             FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
//           }
//         }
        
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

bool DoGenMatching(sdaReader *currentTree, TVector3 LeadPhoton, TVector3 SubleadPhoton) {
  bool ReturnBool = false;
  bool LeadMatched = false;
  bool SubleadMatched = false;
  double Pi = 3.14159265;
  
  for (int i=0; i<currentTree->gp_n; i++) {
    if (currentTree->gp_pdgid[i]!=22) continue;
    TVector3 gp_vtx = (*((TVector3*) currentTree->gp_vtx->At(i)));
    double LeadDeltaEta = abs(LeadPhoton.Eta() - gp_vtx.Eta());
    double LeadDeltaPhi = abs(LeadPhoton.Phi() - gp_vtx.Phi());
    double SubleadDeltaEta = abs(SubleadPhoton.Eta() - gp_vtx.Eta());
    double SubleadDeltaPhi = abs(SubleadPhoton.Phi() - gp_vtx.Phi());

    if (LeadDeltaPhi>Pi) LeadDeltaPhi = 2*Pi-LeadDeltaPhi;
    if (SubleadDeltaPhi>Pi) SubleadDeltaPhi = 2*Pi-SubleadDeltaPhi;

    double LeadDeltaR = sqrt(LeadDeltaPhi*LeadDeltaPhi + LeadDeltaEta*LeadDeltaEta);
    double SubleadDeltaR = sqrt(SubleadDeltaPhi*SubleadDeltaPhi + SubleadDeltaEta*SubleadDeltaEta);

    if (LeadDeltaR<.1) LeadMatched=true;
    if (SubleadDeltaR<.1) SubleadMatched=true;
    if (LeadMatched && SubleadMatched) {
      ReturnBool=true;
      break;
    }
  }

  return ReturnBool;
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

  double theZOfPrimaryVertexFromTracks=-9999.;
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

void ProgressBar(int &percent) {
  cout << "\033[100m";
  cout << "\r[";
  cout << "\033[42m";
  for (int ctr = 0; ctr <= percent / 2; ++ctr) cout << "-";
  cout << "\b>";
  cout << "\033[100m";
  for (int ctr = percent / 2; ctr < 49; ++ctr) cout << " ";
  cout << "]\033[0m";
  cout << " " << ++percent << "%";
  cout << flush;
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

  BookRCutsdZPlots(histoContainer,"dz_conv_");

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
    
  histoContainer->Add("convdZvsEtaAll","#deltaZ of the Primary Vertex from the Conversion and the Sim Vertex vs #eta;#eta of Conversion;#deltaZ of the Primary Vertex from the Conversion from the SimVertex(cm)",60, -3.0, 3.0, 100, -5, 5);
  histoContainer->Add("convdZvsEtaSel","#deltaZ of the Primary Vertex from the Conversion and the Sim Vertex vs #eta;#eta of Conversion;#deltaZ of the Primary Vertex from the Conversion from the SimVertex(cm)",60, -3.0, 3.0, 100, -5, 5);
 
  histoContainer->Add("convdZvsRBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Barrel;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsREndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Endcap;R  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("convdZvsZBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus Z of Conversion: Barrel;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);
  histoContainer->Add("convdZvsZEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus Z of Conversion: Endcap;Z  (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-100,100,100, -5, 5);
    
  histoContainer->Add("convdZPrimaryBarrel","#deltaZ of the Primary Vertex from Conversion and Sim Vertex versus #deltaZ of the Primary Vertex and the Sim Vertex: Barrel;#deltaZ of Primary Vertex from Conversion and SimVertex (cm);#deltaZ of Primary Vertex and SimVertex (cm)",100, -5, 5, 100, -5, 5);
  histoContainer->Add("convdZPrimaryEndcap","#deltaZ of the Primary Vertex from Conversion and Sim Vertex versus #deltaZ of the Primary Vertex and the Sim Vertex: Endcap;#deltaZ of Primary Vertex from Conversion and SimVertex (cm);#deltaZ of Primary Vertex and SimVertex (cm)",100, -5, 5, 100, -5, 5);

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
  BookPhysicsPlots(histoContainer,"leadPhoTrkPtSumSolid03All_region","leading photon trk pt sum dr=03, all candidates: region",100,0.,2);
  BookPhysicsPlots(histoContainer,"leadPhoEcalPtSumSolid03All_region","leading photon ecal pt sum dr=03, all candidates: region",100,0.,5);
  BookPhysicsPlots(histoContainer,"leadPhoHcalPtSumSolid03All_region","leading photon hcal pt sum dr=03, all candidates: region",100,0.,5);
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
  BookMassPlots(histoContainer,"2gamma");
  BookMassPlots(histoContainer,"2gammaGolden");
  BookMassPlots(histoContainer,"2gamma1goodconv");
  BookRCutsMassPlots(histoContainer,"2gamma1goodconv");
  BookMassPlots(histoContainer,"2gamma1goodconvnewvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvlinearvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvsimvertex");
  BookMassPlots(histoContainer,"2gamma1goodconvnearvertex");
  BookMassPlots(histoContainer,"2gamma1poorconv");
  BookMassPlots(histoContainer,"2gamma2conv");
  BookMassPlots(histoContainer,"2gammaleftover");
    
  BookFourHists(histoContainer,"lead_r9_catdummy_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
  BookFourHists(histoContainer,"lead_r9_catdummy_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
  BookFourHists(histoContainer,"lead_r9_catdummy_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_allEcal","subleading photon R9, selected candidates: all ECAL",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_Barrel","subleading photon R9, selected candidates: Barrel",100,0.,1.1);
  BookFourHists(histoContainer,"sublead_r9_catdummy_Endcap","subleading photon R9, selected candidates: Endcap",100,0.,1.1);
    
  histoContainer->Add("convVtxRvsZBarrelAll"," Photon conversion vtx position all candidates Barrel",200, 0., 280., 200, 0., 80.);
  histoContainer->Add("convVtxRvsZBarrelSel"," Photon conversion vtx position selected candidates Barrel",200, 0., 280., 200, 0., 80.);
  histoContainer->Add("selconvVtxRvsZBarrelSel"," Photon conversion vtx position selected conversion Barrel",200, 0., 280., 200, 0., 80.);

}

void BookMassPlots(HistoContainer *histoContainer, TString histname) {

  vector<pair<TString, TString> > parameter;
  parameter.push_back(pair<TString,TString> ("mass_","M_{#gamma#gamma}"));
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

  int bins [6] = {80, 200, 100, 160, 64, 60};
  float lowerlimit [6] = {80.0, 0.0, -1000.0, -8.0, -3.2, 0.0};
  float upperlimit [6] = {160.0, 200.0, 1000.0, 8.0, 3.2, 1.0};
    
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
          
        histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins[k],lowerlimit[k],upperlimit[k]);
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

void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, int i, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop) {

  for (unsigned int j=0; j<2; j++) {
    TString histnames[] = {"phi_conv_", "eta_conv_", "pt_conv_", "z_conv_", "r_conv_", "dz_conv_","eop_conv_"};
    string region = "";
    if (j==0) region = "_allEcal";
    if (j==1 && currentTree->pho_isEB[i]) region = "_Barrel";
    if (j==1 && currentTree->pho_isEE[i]) region = "_Endcap";

    for (unsigned int k=0; k<7; k++) {
      histnames[k] += selection;
      histnames[k] += region;
    }

    histoContainer->Fill(histnames[0].Data(),Photonxyz[i].Phi(), weight);
    histoContainer->Fill(histnames[1].Data(),Photonxyz[i].Eta(), weight);
    histoContainer->Fill(histnames[2].Data(),Photonp4[i].Pt(), weight);
    histoContainer->Fill(histnames[3].Data(),ConversionVertex[i].z(), weight);
    histoContainer->Fill(histnames[4].Data(),ConversionVertex[i].Perp(), weight);
    histoContainer->Fill(histnames[5].Data(),zPVFromConv - SimVertex.Z(), weight);
    histoContainer->Fill(histnames[6].Data(),eop, weight);

  }

}

void FillConvHistsInDiPho(sdaReader *currentTree, HistoContainer *histoContainer, int i, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4, float zPVFromConv, TVector3 SimVertex, float eop) {

  for (unsigned int j=0; j<2; j++) {
    TString histnames[] = {"phi_conv_", "eta_conv_", "pt_conv_", "z_conv_", "r_conv_", "dz_conv_","eop_conv_"};
    string region = "";
    if (j==0) region = "_allEcal";
    if (j==1 && currentTree->pho_isEB[i]) region = "_Barrel";
    if (j==1 && currentTree->pho_isEE[i]) region = "_Endcap";

    for (unsigned int k=0; k<7; k++) {
      histnames[k] += selection;
      histnames[k] += region;
    }

    histoContainer->Fill(histnames[0].Data(),Photonxyz[i].Phi(), weight);
    histoContainer->Fill(histnames[1].Data(),Photonxyz[i].Eta(), weight);
    histoContainer->Fill(histnames[2].Data(),Photonp4[i].Pt(), weight);
    histoContainer->Fill(histnames[3].Data(),ConversionVertex[i].z(), weight);
    histoContainer->Fill(histnames[4].Data(),ConversionVertex[i].Perp(), weight);
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

void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double cos_thetastar, double R) {

  TString histnames[] = {"mass_", "pt_", "pz_", "eta_", "phi_", "CosThetaStar_"};
  for (unsigned int i=0; i<6; i++) {
    histnames[i] += "2gamma";
    histnames[i] += selection;
    histnames[i] += HiggsInWhichDetector;
  }

  if ( InvMass > 100 ) { // get rid of the Z
    histoContainer->Fill(histnames[0].Data(),InvMass,weight);
    histoContainer->Fill(histnames[1].Data(),VSum.Pt(),weight);
    histoContainer->Fill(histnames[2].Data(),VSum.Pz(),weight);
    histoContainer->Fill(histnames[3].Data(),VSum.Eta(),weight);
    histoContainer->Fill(histnames[4].Data(),VSum.Phi(),weight);
    histoContainer->Fill(histnames[5].Data(),cos_thetastar,weight);
    
    if (diPhoCategory==1) for (unsigned int i=0; i<6; i++) histnames[i].ReplaceAll("2gamma","2gammaGolden");
    if (diPhoCategory==2) for (unsigned int i=0; i<6; i++) histnames[i].ReplaceAll("2gamma","2gamma1goodconv");
    if (diPhoCategory==3) for (unsigned int i=0; i<6; i++) histnames[i].ReplaceAll("2gamma","2gamma1poorconv");
    if (diPhoCategory==4) for (unsigned int i=0; i<6; i++) histnames[i].ReplaceAll("2gamma","2gamma2conv");
    if (diPhoCategory==0 || diPhoCategory==5 || diPhoCategory==6) for (unsigned int i=0; i<6; i++) histnames[i].ReplaceAll("2gamma","2gammaleftover");
    
    histoContainer->Fill(histnames[0].Data(),InvMass,weight);
    histoContainer->Fill(histnames[1].Data(),VSum.Pt(),weight);
    histoContainer->Fill(histnames[2].Data(),VSum.Pz(),weight);
    histoContainer->Fill(histnames[3].Data(),VSum.Eta(),weight);
    histoContainer->Fill(histnames[4].Data(),VSum.Phi(),weight);
    histoContainer->Fill(histnames[5].Data(),cos_thetastar,weight);
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
    TString histnames[] = {"PhoEt", "PhoEta", "PhoPhi", "PhoR9", "PhoHoE", "PhoTrkPtSumSolid03", "PhoEcalPtSumSolid03", "PhoHcalPtSumSolid03", "PhoSigmaIetaIeta", "PhoZPV", "PhoDZPV", "PhoDZPVconv"};
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
    histoContainer->Fill(histnames[5].Data(),currentTree->pho_trksumptsolidconedr03[index],weight);
    histoContainer->Fill(histnames[6].Data(),currentTree->pho_ecalsumetconedr03[index],weight);
    histoContainer->Fill(histnames[7].Data(),currentTree->pho_hcalsumetconedr03[index],weight);
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

  if (inputstring.Contains("Data") && !inputstring.Contains("Dataweight")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DataRunA.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DataRunB.root",1));
    isData=true;
  }
  if ((inputstring.Contains("115GeV") && !inputstring.Contains("Pileup")) || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002101;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis115GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH115.root",18.735*BranchingFraction/109991));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF115.root",1.3712*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH115.root",1.2524*BranchingFraction/110000));
  }
  if ((inputstring.Contains("120GeV") && !inputstring.Contains("Pileup")) || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002219;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH120.root",17.173*BranchingFraction/109992));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF120.root",1.3062*BranchingFraction/109842));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH120.root",1.0921*BranchingFraction/110000));
  }
  if ((inputstring.Contains("130GeV") && !inputstring.Contains("Pileup")) || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002240;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH130.root",14.579*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF130.root",1.1866*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH130.root",0.8395*BranchingFraction/110000));
  }
  if (inputstring.Contains("120GeVPileup") || (inputstring.Contains("Pileup") && !inputstring.Contains("GeV")) || inputstring.Contains("All")) {
    BranchingFraction = 0.002219;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG120.root",17.173*BranchingFraction/109991));
  }
  if (inputstring.Contains("130GeVPileup") || (inputstring.Contains("Pileup") && !inputstring.Contains("GeV")) || inputstring.Contains("All")) {
    BranchingFraction = 0.002240;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG130.root",14.579*BranchingFraction/109991));
  }
  if (inputstring.Contains("140GeVPileup") || (inputstring.Contains("Pileup") && !inputstring.Contains("GeV")) || inputstring.Contains("All")) {
    BranchingFraction = 0.001929;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis140Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG140.root",12.525*BranchingFraction/109991));
  }
  if (inputstring.Contains("PJet") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("PhotonJetEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GJet_Pt-20.root",77100/(1182075/0.0064)));
  }
  if (inputstring.Contains("QCD") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDDoubleEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/QCD_Pt-40.root",18700000/(21301935/0.00216)));
  }
  if (inputstring.Contains("Born") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Born.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Born10to25.root",236.4/522865));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Born25to250.root",22.37/537445));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Born250toInf.root",0.008072/546355));
  }
  if (inputstring.Contains("Box") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Box10to25.root",358.2/797975));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Box25to250.root",12.37/777725));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Box250toInf.root",0.000208/789470));
  }
  if (inputstring.Contains("DY") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("DrellYan.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DYEEM1020.root",2659.0/1933000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DYEEM20.root",1300.0/2127607));
  }
  if (inputstring.Contains("Test")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/GGH130.root",1));
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
