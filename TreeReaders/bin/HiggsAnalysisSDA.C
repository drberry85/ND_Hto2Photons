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

bool sortpt(map <double, unsigned int> ptindex, int NumPhotons, double &leadpt, double &subleadpt, unsigned int &leadindex, unsigned int &subleadindex);
unsigned int getconvindex(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
double CosThetaStar(TLorentzVector VLead, TLorentzVector VSum);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindRefittedZ(TVector3 ConversionVertex, TVector3 ConversionRefittedPairMomentum);
string HiggsDetectorPosition(sdaReader *currentTree, unsigned int leadindex, unsigned int subleadindex);
string DetectorPosition(sdaReader *currentTree, unsigned int index);
string MakeFileName(string filename, bool unweighted, bool dataweight, double RCut);
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
void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4);
void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double cos_thetastar, double R);
void FillMassNewVertexHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, string vertextype, TLorentzVector VSum, double InvMass, double cos_thetastar);
void FillMassRcut(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, double R, double InvMass);
void FilldZRcut(HistoContainer *histoContainer, string selection, float weight, string iConvDetector, double R, double deltaz);
void FillPhotonHists(sdaReader *currentTree, HistoContainer *histoContainer, string photon, unsigned int index, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, TVector3 PrimaryVertex, TVector3 SimVertex, vector<TLorentzVector> Photonp4, bool data);
void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData);
void ProgressBar(int &percent);
void PrintWeights();

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool dataweight = false;
  bool unweighted = false;
  double RCut = 999999;
  //float globalWeight = 29.19;
  float globalWeight = 1000;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;
    
  TString InputArgs(input[1]);
  //cout << "Mass is: " << Mass << endl;

  //PrintWeights();
  MakeFilesAndWeights(InputArgs, filesAndWeights, filelist, data);
  
  if (InputArgs.Contains("dataweight")) {
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

    string outfilename = MakeFileName(itFilePair->first, unweighted, dataweight, RCut);
    
    TFile* outfile = new TFile(outfilename.c_str(),"RECREATE");
    outfile->cd();
    cout << outfilename << " created." << endl;

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
      cout << "TreeRead" << endl;
      Long64_t nentries = currentTree.fChain->GetEntries();

      outfile->cd();

      int percent = 0;
      
      for ( Long64_t i = 0; i < nentries; i++ ) {

        TVector3 conversionVertex;

        if (i % (nentries/100) == 0 && bar) {
          if (percent == 100) percent = 99;
          ProgressBar(percent);
        }

        currentTree.GetEntry(i);

        TVector3 SimVertex(0,0,0);
        if (!data) SimVertex.SetXYZ(((TVector3*) currentTree.simvtx->At(0))->x(),((TVector3*) currentTree.simvtx->At(0))->y(),((TVector3*) currentTree.simvtx->At(0))->z());
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionPairMomentum;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TLorentzVector> Photonp4;
        vector<TLorentzVector> SuperClusterp4;

        for (unsigned int j=0; j!=(unsigned int) currentTree.vtx_std_xyz->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) currentTree.vtx_std_xyz->At(j)));
        
        histoContainer->Fill("ZVertex",PrimaryVertex[0].Z(),weight);
        if (!data) {
          histoContainer->Fill("ZSimVertex",SimVertex.Z(),weight);
          histoContainer->Fill("ZdZ",SimVertex.Z()-PrimaryVertex[0].Z(),weight);
          histoContainer->Fill("ZdZZoom",SimVertex.Z()-PrimaryVertex[0].Z(),weight);
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
        histoContainer->Fill("NumberVertices",currentTree.vtx_std_xyz->GetSize(),weight);
        if (!data) histoContainer->Fill("NumberSimVertices",currentTree.simvtx->GetSize(),weight);

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
        int diPhoCategory = diPhotonCategory( leadPhoCategory, subleadPhoCategory );
        if (leadPhoCategory==2 && ConversionVertex[convindex].Perp()>RCut) leadPhoCategory=3;
        if (subleadPhoCategory==2 && ConversionVertex[convindex].Perp()>RCut) subleadPhoCategory=3;
        ////////////////////////////////////
        if (leadPhoCategory==2 || subleadPhoCategory==2) {
          if (!data) FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
          histoContainer->Fill("convetadZallEcal",Photonxyz[convindex].Eta(),currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);
          histoContainer->Fill("convEoPAll",iConvDetector,SuperClusterp4[convscindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          histoContainer->Fill("convr",iConvDetector,ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("Zconv",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex],weight);
          if (!data) {
            histoContainer->Fill("convrdZ",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),ConversionVertex[convindex].Perp(),weight);
            histoContainer->Fill("convdZPrimary",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),PrimaryVertex[0].Z()-SimVertex.Z(),weight);
            histoContainer->Fill("ZconvdZ",iConvDetector,SimVertex.Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex],weight);
          }
          float deltaz = 100000;
          for (unsigned int j=0; j<PrimaryVertex.size(); j++) {
            if (deltaz > abs(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-PrimaryVertex[j].Z())) {
              deltaz = abs(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-PrimaryVertex[j].Z());
              nearvertexindex = j;
            }
          }

          histoContainer->Fill("ZconvdZNearest",iConvDetector,PrimaryVertex[nearvertexindex].Z()-SimVertex.Z(),weight);
          if (!data) histoContainer->Fill("ZdZNear",iConvDetector,currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),weight);

          double RefittedZ = FindRefittedZ(ConversionVertex[convindex], ConversionRefittedPairMomentum[convindex]);
          histoContainer->Fill("RefittedZconv",iConvDetector, RefittedZ, weight);
          if (!data) histoContainer->Fill("RefittedZconvdZ",iConvDetector, RefittedZ-SimVertex.Z(), weight);
          
          TVector3 myBeamSpot(0,0,0);
          double NewZ = FindNewdZ(Photonxyz[convindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          histoContainer->Fill("NewZconv",iConvDetector, NewZ, weight);
          if (!data) histoContainer->Fill("NewZconvdZ",iConvDetector, NewZ-SimVertex.Z(), weight);
          myBeamSpot.SetXYZ(PrimaryVertex[0].x(),PrimaryVertex[0].y(),0);
          NewZ = FindNewdZ(Photonxyz[convindex], ConversionRefittedPairMomentum[convindex], myBeamSpot);
          histoContainer->Fill("NewZPVconv",iConvDetector, NewZ, weight);
          if (!data) histoContainer->Fill("NewZPVconvdZ",iConvDetector, NewZ-SimVertex.Z(), weight);
        }

        TVector3 NearVertex = PrimaryVertex[nearvertexindex];

        /*if (Photonp4[subleadindex].Pt()>Photonp4[leadindex].Pt()) {
          cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
          cout << "LeadPt is " << Photonp4[leadindex].Pt() << endl;
          cout << "SubLeadPt is " << Photonp4[subleadindex].Pt() << endl;
        }*/

        histoContainer->Fill("NPhotonsAll",currentTree.pho_n,weight);
        FillConvHists(&currentTree, histoContainer, selection, weight, Photonxyz, ConversionVertex, Photonp4);
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
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        }
          
        ///////////////////////////////  Event selection ///////////////////////////////////////
        if (Photonp4[leadindex].Pt()<40) continue; // leading photon
        if (Photonp4[subleadindex].Pt()<30) continue; // subleading photon
        //isolation 
        
        if (fabs(Photonxyz[leadindex].Eta())>2.5 || fabs(Photonxyz[subleadindex].Eta())>2.5) continue;
        if (!(looseId(Photonp4[leadindex].Pt(),
		      currentTree.pho_ecalsumetconedr04[leadindex],
		      currentTree.pho_hcalsumetconedr04[leadindex],
		      currentTree.pho_trksumpthollowconedr04[leadindex],
		      (bool) currentTree.pho_isEB[leadindex],
		      (bool) currentTree.pho_isEE[leadindex],
		      currentTree.pho_sieie[leadindex],
		      currentTree.pho_hoe[leadindex]))) continue;

        if (!(looseId(Photonp4[subleadindex].Pt(),
		      currentTree.pho_ecalsumetconedr04[subleadindex],
		      currentTree.pho_hcalsumetconedr04[subleadindex],
		      currentTree.pho_trksumpthollowconedr04[subleadindex],
		      (bool) currentTree.pho_isEB[subleadindex],
		      (bool) currentTree.pho_isEE[subleadindex],
		      currentTree.pho_sieie[subleadindex],
		      currentTree.pho_hoe[subleadindex]))) continue;

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

        histoContainer->Fill("NPhotonsSel",currentTree.pho_n,weight);
        selection = "Sel";
        FillConvHists(&currentTree, histoContainer, selection, weight, Photonxyz, ConversionVertex, Photonp4);
        FillPhotonHists(&currentTree, histoContainer, string("lead"), leadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);
        FillPhotonHists(&currentTree, histoContainer, string("sublead"), subleadindex, selection, weight, Photonxyz, ConversionVertex, PrimaryVertex[0], SimVertex, Photonp4, data);

        if (leadPhoCategory==2 || subleadPhoCategory==2) {
          histoContainer->Fill("convEoPSel",iConvDetector,SuperClusterp4[convscindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          if (!data) FilldZRcut(histoContainer, selection, weight, iConvDetector, ConversionVertex[convindex].Perp(), currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.z());
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
          if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
          FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
        }
        
        //GenMatching
        if (!data) {
          selection="Matched";
          FillMassHists(histoContainer, selection, weight, HiggsInWhichDetector, diPhoCategory, VSum, InvMass, cos_thetastar, ConversionVertex[convindex].Perp());
          if (diPhoCategory==2) {
            FillMassRcut(histoContainer, selection, weight, HiggsInWhichDetector, ConversionVertex[convindex].Perp(), InvMass_newvertex);
            FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "newvertex", VSum_newvertex, InvMass_newvertex, cos_thetastar_newvertex);
            if (!data) FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "simvertex", VSum_simvertex, InvMass_simvertex, cos_thetastar_simvertex);
            FillMassNewVertexHists(histoContainer, selection, weight, HiggsInWhichDetector, "nearvertex", VSum_nearvertex, InvMass_nearvertex, cos_thetastar_nearvertex);
          }
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

string MakeFileName(string filename, bool unweighted, bool dataweight, double RCut) {

  string outfilename = "";
  
  if (unweighted) {
      outfilename = "Unweighted";
      outfilename += filename;
    } else if (dataweight) {
      outfilename = "Dataweight";
      outfilename += filename;
    } else {
      outfilename = filename;
    }
  if (RCut!=999999) outfilename+="RCut";
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
    
    histoContainer->Add("NumberVertices","Number of Reconstructed Vertices;Number of Vertices; Counts",20,0,20);
    histoContainer->Add("NumberSimVertices","Number of Simulated Vertices;Number of Sim Vertices; Counts",20,0,20);

    histoContainer->Add("ZVertex","Z of Primary Vertex;Z (cm); Counts",100,-20,20);
    histoContainer->Add("ZSimVertex","Z of Simulated Vertex;Z (com); Counts",100,-20,20);
    
    histoContainer->Add("ZdZ","#deltaZ between the Z of the Primary Vertex and the Sim Vertex;Z (cm); Counts",100,-1,1);
    histoContainer->Add("ZdZZoom","#deltaZ between the Z of the Primary Vertex and the Sim Vertex: Zoom;Z (cm); Counts",100,-0.02,0.02);
    histoContainer->Add("ZdZAll","#deltaZ between the Z of the all vertecies and the Sim Vertex;Z (cm); Counts",100,-1,1);

    //histoContainer->Add("convr","R of conversion; R (cm); Counts",100,0,100);
    BookRCutsdZPlots(histoContainer,"ZconvdZ");

    BookBarrelAndEndcap(histoContainer,"convr","R of conversion; R (cm): region; Counts",100,0,100);
    BookBarrelAndEndcap(histoContainer,"convEoPAll","E over P of Conversion; E over P; Counts",100,0,3);
    BookBarrelAndEndcap(histoContainer,"convEoPSel","E over P of Conversion; E over P; Counts",100,0,3);
    BookBarrelAndEndcap(histoContainer,"Zconv","Z of Primary Vertex from Conversion: region;Z (cm); Counts",100,-20,20);
    BookBarrelAndEndcap(histoContainer,"ZconvdZ","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-1,1);
    BookBarrelAndEndcap(histoContainer,"ZconvdZNearest","#deltaZ between the Z of the Conversion and the nearest vertex: region;Z (cm); Counts",100,-1,1);
    BookBarrelAndEndcap(histoContainer,"ZdZNear","#deltaZ between the Z of the vertex nearest the conversion Z position and the Sim Vertex: region;Z (cm); Counts",100,-1,1);
    BookBarrelAndEndcap(histoContainer,"NewZconv","Josh's Z of Primary Vertex from Conversion (0,0,0): region;Z (cm); Counts",100,-20,20);
    BookBarrelAndEndcap(histoContainer,"NewZconvdZ","#deltaZ between the Josh's Z of the Primary Vertex from Conversion and Sim Vertex (0,0,0): region;Z (cm); Counts",100,-1,1);
    BookBarrelAndEndcap(histoContainer,"NewZPVconv","Josh's Z of Primary Vertex from Conversion (PV): region;Z (cm); Counts",100,-20,20);
    BookBarrelAndEndcap(histoContainer,"NewZPVconvdZ","#deltaZ between the Josh's Z of the Primary Vertex from Conversion and Sim Vertex (PV): region;Z (cm); Counts",100,-1,1);
    BookBarrelAndEndcap(histoContainer,"RefittedZconv","Z of Primary Vertex from Refitted Conversion: region;Z (cm); Counts",100,-20,20);
    BookBarrelAndEndcap(histoContainer,"RefittedZconvdZ","#deltaZ between the Refitted Z of the Primary Vertex from Conversion and Sim Vertex: region;Z (cm); Counts",100,-1,1);
    
    histoContainer->Add("convetadZallEcal","#eta of the conversion versus #deltaZ of the Primary Vertex from the Conversion and the Sim Vertex;#eta of Conversion;#deltaZ of the Primary Vertex from the Conversion from the SimVertex(cm)",60, -3.0, 3.0, 100, -5, 5);
    
    histoContainer->Add("convrdZBarrel","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Barrel;#deltaZ of Primary Vertex from Conversion (cm);R of Conversion (cm)",100, -5, 5, 0, 100);
    histoContainer->Add("convdZPrimaryBarrel","#deltaZ of the Primary Vertex from Conversion and Sim Vertex versus #deltaZ of the Primary Vertex and the Sim Vertex: Barrel;#deltaZ of Primary Vertex from Conversion and SimVertex (cm);#deltaZ of Primary Vertex and SimVertex (cm)",100, -5, 5, 100, -5, 5);
     
    histoContainer->Add("convrdZEndcap","#deltaZ between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion: Endcap;#deltaZ of Primary Vertex from Conversion (cm);R of Conversion (cm)",100, -5, 5, 0, 100);
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
    BookPhysicsPlots(histoContainer,"leadPhoDZPVAll_region","leading photon #Deltaz_{Zpho - Ztrue} (cm), all candidates: region",100,0.,1);
    BookPhysicsPlots(histoContainer,"leadPhoDZPVconvAll_region","leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: region",100,0.,1);

    BookConversionPlots(histoContainer,"phi_conv_All_region","#phi of Photon Conversion All ECAL; #phi of Conversion: region", 64, -3.2, 3.2);
    BookConversionPlots(histoContainer,"eta_conv_All_region","#eta of Photon Conversion All ECAL; #eta of Conversion: region", 60, -3.0, 3.0);
    BookConversionPlots(histoContainer,"pt_conv_All_region","pt of Photon Conversion All ECAL; pt of Conversion: region", 200, 0, 200);
    BookConversionPlots(histoContainer,"z_conv_All_region","z of Photon Conversion All ECAL; z of Conversion: region", 200, 0, 200);
    BookConversionPlots(histoContainer,"r_conv_All_region","r of Photon Conversion All ECAL; r of Conversion: region", 100, 0, 100);

    // diphoton system
    BookMassPlots(histoContainer,"2gamma");
    BookMassPlots(histoContainer,"2gammaGolden");
    BookMassPlots(histoContainer,"2gamma1goodconv");
    BookRCutsMassPlots(histoContainer,"2gamma1goodconv");
    BookMassPlots(histoContainer,"2gamma1goodconvnewvertex");
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
          if (histname.Contains("simvertex")) histtitletemp += " with Sim Vertex";
          if (histname.Contains("nearvertex")) histtitletemp += " with Near Vertex";
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
    float lowerlimit = -1;
    float upperlimit = 1;
    
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

void FillConvHists(sdaReader *currentTree, HistoContainer *histoContainer, string selection, float weight, vector<TVector3> Photonxyz, vector<TVector3> ConversionVertex, vector<TLorentzVector> Photonp4) {

  for (unsigned int i=0; i<(unsigned int) currentTree->pho_n; i++) {
    if (!currentTree->pho_conv_validvtx[i]) continue;

    for (unsigned int j=0; j<2; j++) {
      TString histnames[] = {"phi_conv_", "eta_conv_", "pt_conv_", "z_conv_", "r_conv_"};
      string region = "";
      if (j==0) region = "_allEcal";
      if (j==1 && currentTree->pho_isEB[i]) region = "_Barrel";
      if (j==1 && currentTree->pho_isEE[i]) region = "_Endcap";

      for (unsigned int k=0; k<5; k++) {
        histnames[k] += selection;
        histnames[k] += region;
      }

      histoContainer->Fill(histnames[0].Data(),Photonxyz[i].Phi(), weight);
      histoContainer->Fill(histnames[1].Data(),Photonxyz[i].Eta(), weight);
      histoContainer->Fill(histnames[2].Data(),Photonp4[i].Pt(), weight);
      histoContainer->Fill(histnames[3].Data(),ConversionVertex[i].z(), weight);
      histoContainer->Fill(histnames[4].Data(),ConversionVertex[i].Perp(), weight);
    }
  }

}

void FillMassHists(HistoContainer *histoContainer, string selection, float weight, string HiggsInWhichDetector, int diPhoCategory, TLorentzVector VSum, double InvMass, double cos_thetastar, double R) {

  TString histnames[] = {"mass_", "pt_", "pz_", "eta_", "phi_", "CosThetaStar_"};
  for (unsigned int i=0; i<6; i++) {
    histnames[i] += "2gamma";
    histnames[i] += selection;
    histnames[i] += HiggsInWhichDetector;
  }

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

void FilldZRcut(HistoContainer *histoContainer, string selection, float weight, string iConvDetector, double R, double deltaz) {

  float RCuts[] = {90, 80, 70, 60, 50, 40 ,30, 20};
  for (unsigned int i=0; i<8; i++) {
    TString histname = "ZconvdZ";
    histname += selection;
    histname += (int) floor(RCuts[i]);
    histname += iConvDetector;
    if (R<RCuts[i]) histoContainer->Fill(histname.Data(),deltaz,weight);
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

  if (inputstring.Contains("Data")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DataRunA.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/DataRunB.root",1));
    isData=true;
  }
  if (inputstring.Contains("115GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002101;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis115GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH115.root",18.735*BranchingFraction/109991));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF115.root",1.3712*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH115.root",1.2524*BranchingFraction/110000));
  }
  if (inputstring.Contains("120GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002219;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH120.root",17.173*BranchingFraction/109992));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF120.root",1.3062*BranchingFraction/109842));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH120.root",1.0921*BranchingFraction/110000));
  }
  if (inputstring.Contains("130GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002240;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GGH130.root",14.579*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/VBF130.root",1.1866*BranchingFraction/108813));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/WZTTH130.root",0.8395*BranchingFraction/110000));
  }
  if (inputstring.Contains("120Pileup") || inputstring.Contains("pileup") || inputstring.Contains("All")) {
    BranchingFraction = 0.002219;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG120.root",17.173*BranchingFraction/109991));
  }
  if (inputstring.Contains("130PileUp") || inputstring.Contains("pileup") || inputstring.Contains("All")) {
    BranchingFraction = 0.002240;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG130.root",14.579*BranchingFraction/109991));
  }
  if (inputstring.Contains("140PileUp") || inputstring.Contains("pileup") || inputstring.Contains("All")) {
    BranchingFraction = 0.001929;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis140Pileup.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/Pileup398/GluGluToHToGG140.root",12.525*BranchingFraction/109991));
  }
  if (inputstring.Contains("PJet") || inputstring.Contains("PJetEMEnriched") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("PhotonJetEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/GJet_Pt-20.root",18700000/(21301935/0.00216)));
  }
  if (inputstring.Contains("QCDEMenriched") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDDoubleEMEnriched.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/SDA/QCD_Pt-40_doubleEMEnriched.root",77100/(1182075/0.0064)));
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
