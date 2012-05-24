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
bool PhotonPreSelection(unsigned int photonindex);
double etaTransformation(double EtaParticle, double Zvertex);
double DeltaPhi(double Phi1, double Phi2);
double MomentumTrackdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double SuperclusterdZ(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
pair<float,float> CombineddZ(TVector3 BeamSpot, unsigned int photonindex, int ConversionIndex);
pair<float,float> MomentumTrackdZErr(unsigned int photonindex, int ConversionIndex);
pair<float,float> SuperClusterdZErr(TVector3 BeamSpot, unsigned int photonindex, int ConversionIndex);
double JacksonAngle(TLorentzVector p1, TLorentzVector p2);
string DetectorPosition(unsigned int index);
string GetPhotonCat(unsigned int index);
map<TString,double> GetkFactor();
map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight);
template <class type> string makestring(type value);
string GetConversionRegion(HistoContainer *histoContainer, double Z, double R, bool isEB);
TLorentzVector GetP4NewVertex(unsigned int photonindex, unsigned int scindex, TVector3 Vertex);
TLorentzVector dgieuler(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector dgloren(TLorentzVector p, double b, double g, double ikey);
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
void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap);
void MakeEtWeights(TString inputstring, map<int,double> &EtMap);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  gROOT->ProcessLine(".L $CMSSW_BASE/src/ND_Hto2Photons/TreeReaders/interface/link_def.h+");
  //gROOT->ProcessLine("TProof")

  TString InputArgs(input[1]);
  bool background = false;
  bool bar = false;
  bool debug = false;
  bool fake = false;
  bool highpt = false;
  bool higgs = false;
  bool mc = true;
  bool nocuts = false;
  bool nojet = false;
  bool onevertex = false;
  bool onephoton = false;
  bool prompt = false;
  bool singleleg = false;
  bool doubleleg = false;
  bool trigger = false;
  bool unweighted = false;
  bool usesimvertex = false;
  float globalweight = 4778.0;

  int FirstFileNum = 0;

  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Background")) background=true;
  if (InputArgs.Contains("Debug")) debug=true;
  if (InputArgs.Contains("Fake")) fake=true;
  if (InputArgs.Contains("HighPt")) highpt=true;
  if (InputArgs.Contains("NoCuts")) nocuts=true;
  if (InputArgs.Contains("NoJet")) nojet=true;
  if (InputArgs.Contains("OneVertex")) onevertex=true;
  if (InputArgs.Contains("OnePhoton")) onephoton=true;
  if (InputArgs.Contains("Prompt")) prompt=true;
  if (InputArgs.Contains("SingleLeg")) singleleg=true;
  if (InputArgs.Contains("DoubleLeg")) doubleleg=true; 
  if (InputArgs.Contains("SimVertex")) usesimvertex=true;
  if (InputArgs.Contains("Trigger")) trigger=true;

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
      if (onevertex) outfilename.ReplaceAll(".root","_OneVertex.root");
      if (usesimvertex) outfilename.ReplaceAll(".root","_SimVertex.root");
      if (background) outfilename.ReplaceAll(".root","_Background.root");
      if (onephoton) outfilename.ReplaceAll(".root","_OnePhoton.root");
      if (fake) outfilename.ReplaceAll(".root","_Fake.root");
      if (highpt) outfilename.ReplaceAll(".root","_HighPt.root");
      if (nocuts) outfilename.ReplaceAll(".root","_NoCuts.root");
      if (prompt) outfilename.ReplaceAll(".root","_Prompt.root");
      if (singleleg) outfilename.ReplaceAll(".root","_SingleLeg.root");
      if (doubleleg) outfilename.ReplaceAll(".root","_DoubleLeg.root");
    }

    if (outfilename.Contains("Higgs")) higgs=true;
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
        if (ientry==0 && gv_pos()==0) mc=false;

        if (ientry==0 && trigger) {
          cout << "Trigger Names:" << endl; 
          for (unsigned int j=0; j<(unsigned int) hlt_path_names_HLT1()->size(); j++) {
            TString TriggerName(hlt_path_names_HLT1()->at(j));
            //if (TriggerName.Contains("Photon")) cout << TriggerName << endl;
            cout << TriggerName << endl;
          }
        }

        if (debug) cout << "Looking at event: " << i << endl;
        TVector3 SimVertex = mc ? *((TVector3*) gv_pos()->At(0)) : TVector3(0,0,0);
        TVector3 DetectorOffset = mc ? TVector3(0,0,0) : TVector3(-0.147,-0.378,-0.485);
        vector<TVector3> BeamSpot;
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionRefittedPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TVector3> SuperClusterxyz;
        vector<TLorentzVector> SuperClusterp4;
        vector<TLorentzVector> Photonp4;
        map<double,int> PhotonPtMap;
        if (pho_n()<1) continue;
        if (pho_n()<2 && nojet) continue;
        if (jet_algo2_n()<1 && !nojet) continue;
        unsigned int leadphotonindex = 0;
        unsigned int subleadphotonindex = 1;
        for (unsigned int j=0; j!=(unsigned int) vtx_std_xyz()->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) vtx_std_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) sc_xyz()->GetSize(); j++) {
          SuperClusterxyz.push_back(*((TVector3*) sc_xyz()->At(j)));
          SuperClusterp4.push_back(*((TLorentzVector*) sc_p4()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) bs_xyz()->GetSize(); j++) BeamSpot.push_back(*((TVector3*) bs_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) pho_n(); j++) {
          Photonp4.push_back(*((TLorentzVector*) pho_p4()->At(j)));
          Photonxyz.push_back(*((TVector3*) pho_calopos()->At(j)));
          PhotonPtMap[Photonp4[j].Pt()]=j;
        }
        map<double,int>::reverse_iterator itLeadPhoton=PhotonPtMap.rbegin();
        map<double,int>::reverse_iterator itSubLeadPhoton=PhotonPtMap.rbegin();
        if (nojet) ++itSubLeadPhoton;
        leadphotonindex=itLeadPhoton->second;
        subleadphotonindex=itSubLeadPhoton->second;
        if (nojet && Photonp4[leadphotonindex].Pt()<Photonp4[subleadphotonindex].Pt()) cout << "Warning Pt is not sorted!!!!!!!!!!!!" << endl;  
        for (unsigned int j=0; j<(unsigned int) conv_n(); j++) {
          //ConversionVertex.push_back(*((TVector3*) conv_vtx()->At(j))+DetectorOffset);
          ConversionVertex.push_back(*((TVector3*) conv_vtx()->At(j)));
          if (conv_singleleg_momentum()!=NULL && conv_ntracks()[j]==1) ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_singleleg_momentum()->At(j))); else ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum()->At(j)));
        }
        float PUWeight = PileUpMap[PrimaryVertex.size()];
        histoContainer->Fill("NumPhoton",pho_n());
        histoContainer->Fill("NumConversion",conv_n());
        
        for (unsigned int photonindex=0; photonindex<(unsigned int) pho_n(); photonindex++) {

          if (onephoton && photonindex!=leadphotonindex) continue;
          if (higgs && photonindex!=leadphotonindex && photonindex!=subleadphotonindex) continue;
          bool GenMatched = false;
          if (mc) GenMatched = GenMatch(Photonxyz[photonindex]);
          if (debug && GenMatched) cout << "Photon is Gen Matched" << endl;
          if (prompt && !GenMatched) continue;
          if (fake && GenMatched) continue;

          int scphotonindex = pho_scind()[photonindex];

          unsigned int BestVertexIndex = vtx_std_sel();
          if (mc && usesimvertex) {
            double mindeltaZtoSimVertex = 999.0;
            for (unsigned int VertexIndex=0; VertexIndex<PrimaryVertex.size(); VertexIndex++) {
              double deltaZtoSimVertex = fabs(PrimaryVertex[VertexIndex].Z()-SimVertex.Z());
              if (deltaZtoSimVertex<mindeltaZtoSimVertex) {
                BestVertexIndex=VertexIndex;
                mindeltaZtoSimVertex=deltaZtoSimVertex;
              }
            }
          }
          
          if (Photonp4[photonindex].Pt()<30) continue;
          if (highpt && Photonp4[photonindex].Pt()<60) continue;
          if (fabs(SuperClusterxyz[scphotonindex].Eta())>2.5) continue;
          //if (!nocuts && pho_cic4cutlevel_lead()->at(photonindex).at(BestVertexIndex)<4 && !background) continue;
          //if (!nocuts && pho_cic4cutlevel_lead()->at(photonindex).at(BestVertexIndex)>=4 && background) continue;
          if (!nocuts && !PhotonPreSelection(photonindex)) continue;
          if (onevertex && PrimaryVertex.size()!=1) continue;
          if (usesimvertex && mc) PrimaryVertex[0]=SimVertex;

          float EtWeight = EtMap.find((int) floor(Photonp4[photonindex].Pt()/5))!=EtMap.end() ? EtMap[(int) floor(Photonp4[photonindex].Pt()/5)] : 1;
          float weight = fileweight*PUWeight*EtWeight;
          
          if (!mc || unweighted) weight=1;
          if (debug) cout << "FileWeight: " << fileweight << " PileUpWeight: " << PUWeight << " EtWeight: " << EtWeight << " TotalWeight: " << weight << endl;
          histoContainer->Fill("PhotonPt",Photonp4[photonindex].Pt(),weight);
          histoContainer->Fill("PhotonPtWeight",Photonp4[photonindex].Pt(),weight);
          histoContainer->Fill("PhotonEta",SuperClusterxyz[scphotonindex].Eta(),weight);
          histoContainer->Fill("PhotonPhi",SuperClusterxyz[scphotonindex].Phi(),weight);
          histoContainer->Fill("Numvtx",PrimaryVertex.size(),weight);

          if (PhotonPreSelection(photonindex)) {
            double tmva_id_mit_tiso1 = pho_tkiso_recvtx_030_002_0000_10_01()->at(photonindex).at(BestVertexIndex) + pho_ecalsumetconedr03()[photonindex] + pho_hcalsumetconedr04()[photonindex];
            double tmva_id_mit_tiso2 = pho_tkiso_badvtx_040_002_0000_10_01()[photonindex] + pho_ecalsumetconedr04()[photonindex] + pho_hcalsumetconedr04()[photonindex];
            histoContainer->Fill("tmva_id_mit_tiso1",rho_algo1(),(float) tmva_id_mit_tiso1);
            histoContainer->Fill("tmva_id_mit_tiso2",rho_algo1(),(float) tmva_id_mit_tiso2);
            histoContainer->Fill("tmva_id_mit_tiso1pro",rho_algo1(),(float) tmva_id_mit_tiso1);
            histoContainer->Fill("tmva_id_mit_tiso2pro",rho_algo1(),(float) tmva_id_mit_tiso2);

            float val_tkiso        = pho_tkiso_recvtx_030_002_0000_10_01()->at(photonindex)[BestVertexIndex];
            float val_ecaliso      = pho_ecalsumetconedr03()[photonindex];
            float val_hcaliso      = pho_hcalsumetconedr04()[photonindex];
            float val_ecalisobad   = pho_ecalsumetconedr04()[photonindex];
            float val_hcalisobad   = pho_hcalsumetconedr04()[photonindex];
            float val_tkisobad     = pho_tkiso_badvtx_040_002_0000_10_01()[photonindex];
            float val_sieie        = pho_sieie()[photonindex];
            float val_hoe          = pho_hoe()[photonindex];
            float val_r9           = pho_r9()[photonindex];
            float val_drtotk_25_99 = pho_drtotk_25_99()[photonindex];
            //float val_pixel        = (float)pho_haspixseed()[photonindex];

            float isosumconst = 5.;
            float isosumconstbad = 7.;
            string Category = GetPhotonCat(photonindex);
            if(Category=="_cat3") {
              isosumconst = 0.;
              isosumconstbad = 0.;
            }
            float rhofacbad=0.354293, rhofac=0.171688;

            TLorentzVector phop4_badvtx = GetP4NewVertex(photonindex,scphotonindex,PrimaryVertex[pho_tkiso_badvtx_id()[photonindex]]);
            float val_isosumoet    = (val_tkiso    + val_ecaliso    + val_hcaliso    + isosumconst    - rho_algo1() * rhofac )   * 50. / Photonp4[photonindex].Et();
            float val_isosumoetbad = (val_tkisobad + val_ecalisobad + val_hcalisobad + isosumconstbad - rho_algo1() * rhofacbad) * 50. / phop4_badvtx.Et();
            float val_trkisooet    = (val_tkiso) * 50. / Photonp4[photonindex].Et();

            histoContainer->Fill("isosumoet",val_isosumoet,weight);
            histoContainer->Fill("isosumoetbad",val_isosumoetbad,weight);
            histoContainer->Fill("trkisooet",val_trkisooet,weight);
            histoContainer->Fill("sieie",val_sieie,weight);
            histoContainer->Fill("hoe",val_hoe,weight);
            histoContainer->Fill("r9",val_r9,weight);
            histoContainer->Fill("drtotk_25_99",val_drtotk_25_99,weight);
            histoContainer->Fill("pho_cic4cutlevel_lead",pho_cic4cutlevel_lead()->at(photonindex).at(BestVertexIndex),weight);
            
            histoContainer->Fill("isosumoet",Category,val_isosumoet,weight);
            histoContainer->Fill("isosumoetbad",Category,val_isosumoetbad,weight);
            histoContainer->Fill("trkisooet",Category,val_trkisooet,weight);
            histoContainer->Fill("sieie",Category,val_sieie,weight);
            histoContainer->Fill("hoe",Category,val_hoe,weight);
            histoContainer->Fill("r9",Category,val_r9,weight);
            histoContainer->Fill("drtotk_25_99",Category,val_drtotk_25_99,weight);
            histoContainer->Fill("pho_cic4cutlevel_lead",Category,pho_cic4cutlevel_lead()->at(photonindex).at(BestVertexIndex),weight);

          }
          
          TLorentzVector MaxJetP4(0,0,0,0);
          TLorentzVector MaxJetTrackP4(0,0,0,0);
          if (!onevertex || !nojet) {
            for (unsigned int j=0; j<(unsigned int) jet_algo2_n(); j++) {
              TLorentzVector JetP4 = *((TLorentzVector*) jet_algo2_p4()->At(j));
              TLorentzVector TrackSumP4(0,0,0,0);
              double dR = (JetP4.Eta()-SuperClusterxyz[scphotonindex].Eta())*(JetP4.Eta()-SuperClusterxyz[scphotonindex].Eta())+DeltaPhi(JetP4.Phi(),SuperClusterxyz[scphotonindex].Phi())*DeltaPhi(JetP4.Phi(),SuperClusterxyz[scphotonindex].Phi());
              if (JetP4.Eta()>2.4 || dR<0.4 || jet_algo2_ntk()[j]<3) continue;
              if (debug) cout << "Good jet found with: " << jet_algo2_ntk()[j] << " tracks" << endl;
              for (unsigned int k=0; k<(unsigned int) jet_algo2_ntk()[j]; k++) {
                //if (debug) cout << "Track Index: " << k << endl;
                unsigned short trackindex = jet_algo2_tkind()->at(j).at(k);
                //if (debug) cout << "Track Index is: " << trackindex << endl;
                if (trackindex>=1500) continue;
                TLorentzVector TrackP4 = *((TLorentzVector*) tk_p4()->At(trackindex));
                if (TrackP4.Pt()>1.0) TrackSumP4+=TrackP4;
              }
              if (JetP4.Pt()>MaxJetP4.Pt() && TrackSumP4.Pt()>30) {
                MaxJetP4=JetP4;
                if (TrackSumP4.Pt()>MaxJetTrackP4.Pt()) MaxJetTrackP4=TrackSumP4;
              }
            }
          }

          if (nojet && onephoton) MaxJetP4=Photonp4[subleadphotonindex];
          if (nojet && photonindex==leadphotonindex) MaxJetP4=Photonp4[subleadphotonindex];
          if (nojet && photonindex==subleadphotonindex) MaxJetP4=Photonp4[leadphotonindex];
          TLorentzVector PJet = MaxJetP4+Photonp4[photonindex];

          histoContainer->Fill("MaxJetPt",MaxJetP4.Pt(),weight);
          histoContainer->Fill("MaxJetTrackPt",MaxJetTrackP4.Pt(),weight);

          histoContainer->Fill("MVAdZNVtxAll",PrimaryVertex.size(),PrimaryVertex[0].Z()-PrimaryVertex[vtx_std_sel()].Z(),weight);
          histoContainer->Fill("MVAdZVtxPtAll",PJet.Pt(),PrimaryVertex[0].Z()-PrimaryVertex[vtx_std_sel()].Z(),weight);
          
          int convindex = gettrackerconvindex(SuperClusterxyz[scphotonindex],BeamSpot[0]);
          if (convindex==-1) continue;
          if (singleleg && conv_ntracks()[convindex]!=1) continue;
          if (doubleleg && conv_ntracks()[convindex]!=2) continue;
          histoContainer->Fill("ConvPt",Photonp4[photonindex].Pt(),weight);
          histoContainer->Fill("ConvPairPt",ConversionRefittedPairMomentum[convindex].Pt(),weight);
          histoContainer->Fill("ConvEta",ConversionVertex[convindex].Eta(),weight);
          histoContainer->Fill("ConvPhi",ConversionVertex[convindex].Phi(),weight);
          histoContainer->Fill("ConvR",ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("ConvXvsZ",ConversionVertex[convindex].Z(),ConversionVertex[convindex].X(),weight);
          histoContainer->Fill("ConvYvsZ",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Y(),weight);
          histoContainer->Fill("ConvRvsEta",SuperClusterxyz[scphotonindex].Eta(),ConversionVertex[convindex].Perp(),weight);
          histoContainer->Fill("ConvRvsPhi",SuperClusterxyz[scphotonindex].Phi(),ConversionVertex[convindex].Perp(),weight);

          if (fabs(ConversionVertex[convindex].Z())<26.0) histoContainer->Fill("ConvRTracker",ConversionVertex[convindex].Perp(),weight);
          if (ConversionVertex[convindex].Perp()>3.5 && ConversionVertex[convindex].Perp()<19.0) histoContainer->Fill("ConvZPixel",ConversionVertex[convindex].Z(),weight);
          
          if (ConversionVertex[convindex].Phi()<0) histoContainer->Fill("ConvRvsZ",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (ConversionVertex[convindex].Phi()>=0) histoContainer->Fill("ConvRvsZ",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (-3.15<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-2.80) histoContainer->Fill("ConvRvsZslice1",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-2.80<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-2.45) histoContainer->Fill("ConvRvsZslice2",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-2.45<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-2.10) histoContainer->Fill("ConvRvsZslice3",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-2.10<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-1.75) histoContainer->Fill("ConvRvsZslice4",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-1.75<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-1.40) histoContainer->Fill("ConvRvsZslice5",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-1.40<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-1.05) histoContainer->Fill("ConvRvsZslice6",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-1.05<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-0.70) histoContainer->Fill("ConvRvsZslice7",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-0.70<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-0.35) histoContainer->Fill("ConvRvsZslice8",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (-0.35<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<-0.00) histoContainer->Fill("ConvRvsZslice9",ConversionVertex[convindex].Z(),-ConversionVertex[convindex].Perp(),weight);
          if (0.00<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<0.35) histoContainer->Fill("ConvRvsZslice1",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (0.35<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<0.70) histoContainer->Fill("ConvRvsZslice2",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (0.70<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<1.05) histoContainer->Fill("ConvRvsZslice3",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (1.05<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<1.40) histoContainer->Fill("ConvRvsZslice4",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (1.40<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<1.75) histoContainer->Fill("ConvRvsZslice5",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (1.75<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<2.10) histoContainer->Fill("ConvRvsZslice6",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (2.10<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<2.45) histoContainer->Fill("ConvRvsZslice7",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (2.45<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<2.80) histoContainer->Fill("ConvRvsZslice8",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);
          if (2.80<=ConversionVertex[convindex].Phi() && ConversionVertex[convindex].Phi()<3.15) histoContainer->Fill("ConvRvsZslice9",ConversionVertex[convindex].Z(),ConversionVertex[convindex].Perp(),weight);

          if (-1.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.1) histoContainer->Fill("SCConvR-Mod4",ConversionVertex[convindex].Perp(),weight);
          if (-1.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.3) histoContainer->Fill("SCConvR-1.5",ConversionVertex[convindex].Perp(),weight);
          if (-1.3<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.1) histoContainer->Fill("SCConvR-1.3",ConversionVertex[convindex].Perp(),weight);
          if (-1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-0.8) histoContainer->Fill("SCConvR-1.1",ConversionVertex[convindex].Perp(),weight);
          if (-0.8<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-0.5) histoContainer->Fill("SCConvR-0.8",ConversionVertex[convindex].Perp(),weight);
          if (-0.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.0 ) histoContainer->Fill("SCConvR-0.5",ConversionVertex[convindex].Perp(),weight);
          if ( 0.0<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.5 ) histoContainer->Fill("SCConvR0.0",ConversionVertex[convindex].Perp(),weight);
          if ( 0.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.8 ) histoContainer->Fill("SCConvR0.5",ConversionVertex[convindex].Perp(),weight);
          if ( 0.8<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.1 ) histoContainer->Fill("SCConvR0.8",ConversionVertex[convindex].Perp(),weight);
          if ( 1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.3 ) histoContainer->Fill("SCConvR1.1",ConversionVertex[convindex].Perp(),weight);
          if ( 1.3<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.5 ) histoContainer->Fill("SCConvR1.3",ConversionVertex[convindex].Perp(),weight);
          if ( 1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.5 ) histoContainer->Fill("SCConvRMod4",ConversionVertex[convindex].Perp(),weight);

          if (ConversionVertex[convindex].Perp()<25.0) {
            if (-1.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.1) histoContainer->Fill("SCConvRFine-Mod4",ConversionVertex[convindex].Perp(),weight);
            if (-1.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.3) histoContainer->Fill("SCConvRFine-1.5",ConversionVertex[convindex].Perp(),weight);
            if (-1.3<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-1.1) histoContainer->Fill("SCConvRFine-1.3",ConversionVertex[convindex].Perp(),weight);
            if (-1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-0.8) histoContainer->Fill("SCConvRFine-1.1",ConversionVertex[convindex].Perp(),weight);
            if (-0.8<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<-0.5) histoContainer->Fill("SCConvRFine-0.8",ConversionVertex[convindex].Perp(),weight);
            if (-0.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.0 ) histoContainer->Fill("SCConvRFine-0.5",ConversionVertex[convindex].Perp(),weight);
            if ( 0.0<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.5 ) histoContainer->Fill("SCConvRFine0.0",ConversionVertex[convindex].Perp(),weight);
            if ( 0.5<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<0.8 ) histoContainer->Fill("SCConvRFine0.5",ConversionVertex[convindex].Perp(),weight);
            if ( 0.8<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.1 ) histoContainer->Fill("SCConvRFine0.8",ConversionVertex[convindex].Perp(),weight);
            if ( 1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.3 ) histoContainer->Fill("SCConvRFine1.1",ConversionVertex[convindex].Perp(),weight);
            if ( 1.3<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.5 ) histoContainer->Fill("SCConvRFine1.3",ConversionVertex[convindex].Perp(),weight);
            if ( 1.1<=SuperClusterxyz[scphotonindex].Eta() && SuperClusterxyz[scphotonindex].Eta()<1.5 ) histoContainer->Fill("SCConvRFineMod4",ConversionVertex[convindex].Perp(),weight);
          }

          if (-1.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.1) histoContainer->Fill("ConvR-Mod4",ConversionVertex[convindex].Perp(),weight);
          if (-1.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.3) histoContainer->Fill("ConvR-1.5",ConversionVertex[convindex].Perp(),weight);
          if (-1.3<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.1) histoContainer->Fill("ConvR-1.3",ConversionVertex[convindex].Perp(),weight);
          if (-1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-0.8) histoContainer->Fill("ConvR-1.1",ConversionVertex[convindex].Perp(),weight);
          if (-0.8<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-0.5) histoContainer->Fill("ConvR-0.8",ConversionVertex[convindex].Perp(),weight);
          if (-0.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.0 ) histoContainer->Fill("ConvR-0.5",ConversionVertex[convindex].Perp(),weight);
          if ( 0.0<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.5 ) histoContainer->Fill("ConvR0.0",ConversionVertex[convindex].Perp(),weight);
          if ( 0.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.8 ) histoContainer->Fill("ConvR0.5",ConversionVertex[convindex].Perp(),weight);
          if ( 0.8<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.1 ) histoContainer->Fill("ConvR0.8",ConversionVertex[convindex].Perp(),weight);
          if ( 1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.3 ) histoContainer->Fill("ConvR1.1",ConversionVertex[convindex].Perp(),weight);
          if ( 1.3<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.5 ) histoContainer->Fill("ConvR1.3",ConversionVertex[convindex].Perp(),weight);
          if ( 1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.5 ) histoContainer->Fill("ConvRMod4",ConversionVertex[convindex].Perp(),weight);

          if (ConversionVertex[convindex].Perp()<25.0) {
            if (-1.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.1) histoContainer->Fill("ConvRFine-Mod4",ConversionVertex[convindex].Perp(),weight);
            if (-1.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.3) histoContainer->Fill("ConvRFine-1.5",ConversionVertex[convindex].Perp(),weight);
            if (-1.3<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-1.1) histoContainer->Fill("ConvRFine-1.3",ConversionVertex[convindex].Perp(),weight);
            if (-1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-0.8) histoContainer->Fill("ConvRFine-1.1",ConversionVertex[convindex].Perp(),weight);
            if (-0.8<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<-0.5) histoContainer->Fill("ConvRFine-0.8",ConversionVertex[convindex].Perp(),weight);
            if (-0.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.0 ) histoContainer->Fill("ConvRFine-0.5",ConversionVertex[convindex].Perp(),weight);
            if ( 0.0<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.5 ) histoContainer->Fill("ConvRFine0.0",ConversionVertex[convindex].Perp(),weight);
            if ( 0.5<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<0.8 ) histoContainer->Fill("ConvRFine0.5",ConversionVertex[convindex].Perp(),weight);
            if ( 0.8<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.1 ) histoContainer->Fill("ConvRFine0.8",ConversionVertex[convindex].Perp(),weight);
            if ( 1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.3 ) histoContainer->Fill("ConvRFine1.1",ConversionVertex[convindex].Perp(),weight);
            if ( 1.3<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.5 ) histoContainer->Fill("ConvRFine1.3",ConversionVertex[convindex].Perp(),weight);
            if ( 1.1<=ConversionVertex[convindex].Eta() && ConversionVertex[convindex].Eta()<1.5 ) histoContainer->Fill("ConvRFineMod4",ConversionVertex[convindex].Perp(),weight);
          }

          if (!onevertex && !nojet && MaxJetP4.Pt()<30) continue;

          histoContainer->Fill("PairMass",PJet.M(),weight);
          histoContainer->Fill("PairPt",PJet.Pt(),weight);
          histoContainer->Fill("PhotonJetJacksonAngle",JacksonAngle(Photonp4[photonindex],MaxJetP4),weight);

          unsigned int mvavertexindex = photonindex;
          if (nojet) for (unsigned int j=0; j<(unsigned int) dipho_n(); j++) if ((unsigned int) dipho_leadind()[j]==leadphotonindex && (unsigned int) dipho_subleadind()[j]==subleadphotonindex) mvavertexindex=j;
          if (debug) cout << "MVAVertex Index: " << mvavertexindex << " MVA Length is: " << vtx_std_mva()->size() << endl;
          if (vtx_std_mva()->size()==0 || vtx_std_mva()->size()<=mvavertexindex) continue;

          //unsigned int MVAVertexIndex = vtx_std_ranked_list()->at(mvavertexindex).at(0);
          unsigned int MVAVertexIndex = vtx_std_sel();

          histoContainer->Fill("MVAValueFirst",vtx_std_mva()->at(mvavertexindex).at(MVAVertexIndex),weight);
          histoContainer->Fill("MVAdZFirst",PrimaryVertex[MVAVertexIndex].Z()-PrimaryVertex[0].Z(),weight);
          if (vtx_std_ranked_list()->at(mvavertexindex).size()>1) {
            histoContainer->Fill("MVAValueSecond",vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(1)),weight);
            histoContainer->Fill("MVAdZSecond",PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(1)].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          }
          if (vtx_std_ranked_list()->at(mvavertexindex).size()>2) {
            histoContainer->Fill("MVAValueThird",vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(2)),weight);
            histoContainer->Fill("MVAdZThird",PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(2)].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          }

          for (unsigned int j=0; j<vtx_std_mva()->at(mvavertexindex).size(); j++) {
            histoContainer->Fill("logsumpt2All",log(vtx_std_sumpt2()->at(mvavertexindex).at(j)),weight);
            if (j==(unsigned int) MVAVertexIndex && fabs(PrimaryVertex[0].Z()-PrimaryVertex[j].Z())<1.0) {
              histoContainer->Fill("sumpt2",vtx_std_sumpt2()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("logsumpt2",log(vtx_std_sumpt2()->at(mvavertexindex).at(j)),weight);
              histoContainer->Fill("ptasymm",vtx_std_ptasym()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("ptbal",vtx_std_ptbal()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("pulltoconv",vtx_std_pulltoconv()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("limpulltoconv",vtx_std_limpulltoconv()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("MVAValue",vtx_std_mva()->at(mvavertexindex).at(j),weight);
            } else if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[j].Z())>1.0) {
              histoContainer->Fill("sumpt2bad",vtx_std_sumpt2()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("logsumpt2bad",log(vtx_std_sumpt2()->at(mvavertexindex).at(j)),weight);
              histoContainer->Fill("ptasymmbad",vtx_std_ptasym()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("ptbalbad",vtx_std_ptbal()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("pulltoconvbad",vtx_std_pulltoconv()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("limpulltoconvbad",vtx_std_limpulltoconv()->at(mvavertexindex).at(j),weight);
              histoContainer->Fill("MVAValuebad",vtx_std_mva()->at(mvavertexindex).at(j),weight);
            }
          }

          double superclusterdz = SuperclusterdZ(ConversionVertex[convindex],SuperClusterxyz[scphotonindex],BeamSpot[0]);
          double convvertexdz = MomentumTrackdZ(ConversionVertex[convindex],ConversionRefittedPairMomentum[convindex],BeamSpot[0]);
          pair<double,double> rishisuperclusterdzpair(0,0);
          pair<double,double> rishitrackprojectiondzpair(0,0);
          pair<double,double> rishicombineddzpair(0,0);
          if (conv_tk1_thetaerr()!=0) {
            pair<double,double> rishisuperclusterdzpair = SuperClusterdZErr(BeamSpot[0],photonindex,convindex);
            pair<double,double> rishitrackprojectiondzpair = MomentumTrackdZErr(photonindex,convindex);
            pair<double,double> rishicombineddzpair = CombineddZ(BeamSpot[0],photonindex,convindex);
          }
          string PhotonDetector = DetectorPosition(photonindex);
          histoContainer->Fill("PhotonEoP",SuperClusterp4[scphotonindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          if (GenMatched) histoContainer->Fill("PhotonEoPGen",SuperClusterp4[scphotonindex].E()/ConversionRefittedPairMomentum[convindex].Mag(),weight);
          histoContainer->Fill("ZPV",PrimaryVertex[0].Z(),weight);
          if (mc) histoContainer->Fill("dZcheck",SimVertex.Z()-PrimaryVertex[0].Z(),weight);
          if (mc) histoContainer->Fill("dZcheckzoom",(SimVertex.Z()-PrimaryVertex[0].Z())*10*1000,weight);
          histoContainer->Fill("deltaEta",SuperClusterxyz[scphotonindex].Eta() - etaTransformation(ConversionRefittedPairMomentum[convindex].Eta(),superclusterdz),weight);

          histoContainer->Fill("SuperZPV",superclusterdz,weight);
          histoContainer->Fill("SuperdZ",PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("SuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("SuperdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-superclusterdz,weight);
          FilldZPt(histoContainer,"SuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
          FilldZEta(histoContainer,"SuperdZEta",fabs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-superclusterdz,weight);
          
          histoContainer->Fill("ConvZPV",convvertexdz,weight);
          histoContainer->Fill("ConvdZ",PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("ConvdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("ConvdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-convvertexdz,weight);
          FilldZPt(histoContainer,"ConvdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
          FilldZEta(histoContainer,"ConvdZEta",fabs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-convvertexdz,weight);

          histoContainer->Fill("RishiSuperZPV",rishisuperclusterdzpair.first,weight);
          histoContainer->Fill("RishiSuperdZ",PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);
          histoContainer->Fill("RishiSuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);
          histoContainer->Fill("RishiSuperdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);
          FilldZPt(histoContainer,"RishiSuperdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);
          FilldZEta(histoContainer,"RishiSuperdZEta",fabs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);

          histoContainer->Fill("RishiTrackZPV",rishitrackprojectiondzpair.first,weight);
          histoContainer->Fill("RishiTrackdZ",PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);
          histoContainer->Fill("RishiTrackdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);
          histoContainer->Fill("RishiTrackdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);
          FilldZPt(histoContainer,"RishiTrackdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);
          FilldZEta(histoContainer,"RishiTrackdZEta",fabs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);

          histoContainer->Fill("RishiCombinedZPV",rishicombineddzpair.first,weight);
          histoContainer->Fill("RishiCombineddZ",PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          histoContainer->Fill("RishiCombineddZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          histoContainer->Fill("RishiCombineddZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          FilldZPt(histoContainer,"RishiCombineddZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          FilldZEta(histoContainer,"RishiCombineddZEta",fabs(SuperClusterxyz[scphotonindex].Eta()),PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          
          histoContainer->Fill("MVAdZ",PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          histoContainer->Fill("MVAdZNVtx",PrimaryVertex.size(),PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          histoContainer->Fill("MVAdZNVtx_short",PrimaryVertex.size(),PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          histoContainer->Fill("MVAdZVtxPt",PJet.Pt(),PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
          if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())<1.0) histoContainer->Fill("MVARes",vtx_std_evt_mva()->at(mvavertexindex),weight);
          else histoContainer->Fill("MVAResbad",vtx_std_evt_mva()->at(mvavertexindex),weight);
          
          if (mc) {
            histoContainer->Fill("MVAdZsim",PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z(),weight);
            if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0) histoContainer->Fill("MVAdZsimbad",PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z(),weight);
            if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0 && vtx_std_evt_mva()->at(mvavertexindex)<-0.9) histoContainer->Fill("MVAdZsimfunky",PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z(),weight);
            if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0 && fabs(PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z())<1.0) histoContainer->Fill("MVAResfunky",vtx_std_evt_mva()->at(mvavertexindex),weight);
            if (fabs(PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z())>1.0) histoContainer->Fill("MVAResbadsim",vtx_std_evt_mva()->at(mvavertexindex),weight);
            if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0 && fabs(PrimaryVertex[MVAVertexIndex].Z()-SimVertex.Z())<1.0) histoContainer->Fill("MVAResgoodbadsim",vtx_std_evt_mva()->at(mvavertexindex),weight);
          }
          
          if (debug && fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0 && vtx_std_evt_mva()->at(mvavertexindex)<-0.9 ) {
            cout << "Delta Z is: " << fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z()) << " and Per Event MVA is: " << vtx_std_evt_mva()->at(mvavertexindex) << " and MVA Vertex Value is: " << vtx_std_mva()->at(mvavertexindex).at(MVAVertexIndex) << " MVAVertex Index: " << mvavertexindex << " Vertex Index: " << MVAVertexIndex << endl;
            cout << "SumPt2:                   " << vtx_std_sumpt2()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "LogSumPt2:                " << log(vtx_std_sumpt2()->at(mvavertexindex).at(MVAVertexIndex)) << endl;
            cout << "Pt Asymmetry              " << vtx_std_ptasym()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "Pt Balance:               " << vtx_std_ptbal()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "Pull to Conversion:       " << vtx_std_pulltoconv()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "Limit Pull to Conversion: " << vtx_std_limpulltoconv()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "MVA Value:                " << vtx_std_mva()->at(mvavertexindex).at(MVAVertexIndex) << endl;
            cout << "Second MVA Value:         " << vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(1)) << endl;
            cout << "Third MVA Value:          " << vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(2)) << endl;
            cout << "Pair Pt:                  " << PJet.Pt() << endl;
            cout << "Second Delta Z:           " << abs(PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(1)].Z()-PrimaryVertex[MVAVertexIndex].Z()) << endl;
            cout << "Third Delta Z:            " << abs(PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(2)].Z()-PrimaryVertex[MVAVertexIndex].Z()) << endl;
          }

          if (fabs(PrimaryVertex[0].Z()-PrimaryVertex[MVAVertexIndex].Z())>1.0 && vtx_std_evt_mva()->at(mvavertexindex)<-0.95) {
            histoContainer->Fill("sumpt2busted",vtx_std_sumpt2()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("logsumpt2busted",log(vtx_std_sumpt2()->at(mvavertexindex).at(MVAVertexIndex)),weight);
            histoContainer->Fill("ptasymmbusted",vtx_std_ptasym()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("ptbalbusted",vtx_std_ptbal()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("pulltoconvbusted",vtx_std_pulltoconv()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("limpulltoconvbusted",vtx_std_limpulltoconv()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("MVAValuebusted",vtx_std_mva()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("PairPtbusted",PJet.Pt(),weight);
            histoContainer->Fill("MVAValueFirstbusted",vtx_std_mva()->at(mvavertexindex).at(MVAVertexIndex),weight);
            histoContainer->Fill("MVAdZFirstbusted",PrimaryVertex[MVAVertexIndex].Z()-PrimaryVertex[0].Z(),weight);
            if (vtx_std_ranked_list()->at(mvavertexindex).size()>1) {
              histoContainer->Fill("MVAValueSecondbusted",vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(1)),weight);
              histoContainer->Fill("MVAdZSecondbusted",PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(1)].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
            }
            if (vtx_std_ranked_list()->at(mvavertexindex).size()>2) {
              histoContainer->Fill("MVAValueThirdbusted",vtx_std_mva()->at(mvavertexindex).at(vtx_std_ranked_list()->at(mvavertexindex).at(2)),weight);
              histoContainer->Fill("MVAdZThirdbusted",PrimaryVertex[vtx_std_ranked_list()->at(mvavertexindex).at(2)].Z()-PrimaryVertex[MVAVertexIndex].Z(),weight);
            }
          }
          
          histoContainer->Fill("ConvdZvsR",ConversionVertex[convindex].Perp(),PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("ConvdZvsZ",ConversionVertex[convindex].Z(),PrimaryVertex[0].Z()-convvertexdz,weight);
          histoContainer->Fill("SuperdZvsR",ConversionVertex[convindex].Perp(),PrimaryVertex[0].Z()-superclusterdz,weight);
          histoContainer->Fill("SuperdZvsZ",ConversionVertex[convindex].Z(),PrimaryVertex[0].Z()-superclusterdz,weight);
          
          if (mc) histoContainer->Fill("SuperdZsim",SimVertex.Z()-superclusterdz,weight);
          if (mc) histoContainer->Fill("ConvdZsim",SimVertex.Z()-convvertexdz,weight);

          histoContainer->Fill("SuperZPV",PhotonDetector,superclusterdz,weight);
          histoContainer->Fill("SuperdZ",PhotonDetector,PrimaryVertex[0].Z()-superclusterdz,weight);
          if (mc) histoContainer->Fill("SuperdZsim",PhotonDetector,SimVertex.Z()-superclusterdz,weight);

          histoContainer->Fill("ConvZPV",PhotonDetector,convvertexdz,weight);
          histoContainer->Fill("ConvdZ",PhotonDetector,PrimaryVertex[0].Z()-convvertexdz,weight);
          if (mc) histoContainer->Fill("ConvdZsim",PhotonDetector,SimVertex.Z()-convvertexdz,weight);

          histoContainer->Fill("RishiSuperZPV",PhotonDetector,rishisuperclusterdzpair.first,weight);
          histoContainer->Fill("RishiSuperdZ",PhotonDetector,PrimaryVertex[0].Z()-rishisuperclusterdzpair.first,weight);
          if (mc) histoContainer->Fill("RishiSuperdZsim",PhotonDetector,SimVertex.Z()-rishisuperclusterdzpair.first,weight);

          histoContainer->Fill("RishiTrackZPV",PhotonDetector,rishitrackprojectiondzpair.first,weight);
          histoContainer->Fill("RishiTrackdZ",PhotonDetector,PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first,weight);
          if (mc) histoContainer->Fill("RishiTrackdZsim",PhotonDetector,SimVertex.Z()-rishitrackprojectiondzpair.first,weight);

          histoContainer->Fill("RishiCombinedZPV",PhotonDetector,rishicombineddzpair.first,weight);
          histoContainer->Fill("RishiCombineddZ",PhotonDetector,PrimaryVertex[0].Z()-rishicombineddzpair.first,weight);
          if (mc) histoContainer->Fill("RishiCombineddZsim",PhotonDetector,SimVertex.Z()-rishicombineddzpair.first,weight);
          
          if (PhotonDetector=="Barrel") {
            FilldZTrackerBarrel(histoContainer, "SuperdZ", PrimaryVertex[0].Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "SuperdZEff", PrimaryVertex[0].Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight, true);
            FilldZTrackerBarrel(histoContainer, "SuperdZRes", PrimaryVertex[0].Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight);
            if (mc) FilldZTrackerBarrel(histoContainer, "SuperdZsim", SimVertex.Z()-superclusterdz, ConversionVertex[convindex].Perp(), weight);

            FilldZTrackerBarrel(histoContainer, "RishiSuperdZ", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "RishiSuperdZEff", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, ConversionVertex[convindex].Perp(), weight, true);
            FilldZTrackerBarrel(histoContainer, "RishiSuperdZRes", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, ConversionVertex[convindex].Perp(), weight);
            if (mc) FilldZTrackerBarrel(histoContainer, "RishiSuperdZsim", SimVertex.Z()-rishisuperclusterdzpair.first, ConversionVertex[convindex].Perp(), weight);
            
            FilldZTrackerBarrel(histoContainer, "ConvdZ", PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "ConvdZEff", PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight, true);
            FilldZTrackerBarrel(histoContainer, "ConvdZRes", PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "ConvdZCompairison", PrimaryVertex[0].Z()-superclusterdz, PrimaryVertex[0].Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);
            if (mc) FilldZTrackerBarrel(histoContainer, "ConvdZsim", SimVertex.Z()-convvertexdz, ConversionVertex[convindex].Perp(), weight);

            FilldZTrackerBarrel(histoContainer, "RishiTrackdZ", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "RishiTrackdZEff", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, ConversionVertex[convindex].Perp(), weight, true);
            FilldZTrackerBarrel(histoContainer, "RishiTrackdZRes", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "RishiTrackdZCompairison", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, ConversionVertex[convindex].Perp(), weight);
            if (mc) FilldZTrackerBarrel(histoContainer, "RishiTrackdZsim", SimVertex.Z()-rishitrackprojectiondzpair.first, ConversionVertex[convindex].Perp(), weight);

            FilldZTrackerBarrel(histoContainer, "RishiCombineddZ", PrimaryVertex[0].Z()-rishicombineddzpair.first, ConversionVertex[convindex].Perp(), weight);
            FilldZTrackerBarrel(histoContainer, "RishiCombineddZEff", PrimaryVertex[0].Z()-rishicombineddzpair.first, ConversionVertex[convindex].Perp(), weight, true);
            FilldZTrackerBarrel(histoContainer, "RishiCombineddZRes", PrimaryVertex[0].Z()-rishicombineddzpair.first, ConversionVertex[convindex].Perp(), weight);
            if (mc) FilldZTrackerBarrel(histoContainer, "RishiCombineddZsim", SimVertex.Z()-rishicombineddzpair.first, ConversionVertex[convindex].Perp(), weight);
          }
          if (PhotonDetector=="Endcap") {
            FilldZTrackerEndcap(histoContainer, "SuperdZ", PrimaryVertex[0].Z()-superclusterdz, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "SuperdZEff", PrimaryVertex[0].Z()-superclusterdz, fabs(ConversionVertex[convindex].Z()), weight, true);
            FilldZTrackerEndcap(histoContainer, "SuperdZRes", PrimaryVertex[0].Z()-superclusterdz, fabs(ConversionVertex[convindex].Z()), weight);
            if (mc) FilldZTrackerEndcap(histoContainer, "SuperdZsim", SimVertex.Z()-superclusterdz, fabs(ConversionVertex[convindex].Z()), weight);

            FilldZTrackerEndcap(histoContainer, "RishiSuperdZ", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "RishiSuperdZEff", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, fabs(ConversionVertex[convindex].Z()), weight, true);
            FilldZTrackerEndcap(histoContainer, "RishiSuperdZRes", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            if (mc) FilldZTrackerEndcap(histoContainer, "RishiSuperdZsim", SimVertex.Z()-rishisuperclusterdzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            
            FilldZTrackerEndcap(histoContainer, "ConvdZ", PrimaryVertex[0].Z()-convvertexdz, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "ConvdZEff", PrimaryVertex[0].Z()-convvertexdz, fabs(ConversionVertex[convindex].Z()), weight, true);
            FilldZTrackerEndcap(histoContainer, "ConvdZRes", PrimaryVertex[0].Z()-convvertexdz, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "ConvdZCompairison", PrimaryVertex[0].Z()-superclusterdz, PrimaryVertex[0].Z()-convvertexdz, fabs(ConversionVertex[convindex].Z()), weight);
            if (mc) FilldZTrackerEndcap(histoContainer, "ConvdZsim", SimVertex.Z()-convvertexdz, fabs(ConversionVertex[convindex].Z()), weight);

            FilldZTrackerEndcap(histoContainer, "RishiTrackdZ", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "RishiTrackdZEff", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, fabs(ConversionVertex[convindex].Z()), weight, true);
            FilldZTrackerEndcap(histoContainer, "RishiTrackdZRes", PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "RishiTrackdZCompairison", PrimaryVertex[0].Z()-rishisuperclusterdzpair.first, PrimaryVertex[0].Z()-rishitrackprojectiondzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            if (mc) FilldZTrackerEndcap(histoContainer, "RishiTrackdZsim", SimVertex.Z()-rishitrackprojectiondzpair.first, fabs(ConversionVertex[convindex].Z()), weight);

            FilldZTrackerEndcap(histoContainer, "RishiCombineddZ", PrimaryVertex[0].Z()-rishicombineddzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            FilldZTrackerEndcap(histoContainer, "RishiCombineddZEff", PrimaryVertex[0].Z()-rishicombineddzpair.first, fabs(ConversionVertex[convindex].Z()), weight, true);
            FilldZTrackerEndcap(histoContainer, "RishiCombineddZRes", PrimaryVertex[0].Z()-rishicombineddzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
            if (mc) FilldZTrackerEndcap(histoContainer, "RishiCombineddZsim", SimVertex.Z()-rishicombineddzpair.first, fabs(ConversionVertex[convindex].Z()), weight);
          }

          string ConversionRegion = GetConversionRegion(histoContainer, ConversionVertex[convindex].Z(), ConversionVertex[convindex].Perp(), pho_isEB()[photonindex]);
          if ((ConversionRegion=="PixelBarrel" || ConversionRegion=="PixelFwd" || ConversionRegion=="TID") && conv_ntracks()[convindex]!=1) {
            histoContainer->Fill("MixdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
            histoContainer->Fill("MixdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-convvertexdz,weight);
            histoContainer->Fill("MixdZNVtx",PrimaryVertex.size(),PrimaryVertex[0].Z()-convvertexdz,weight);
            histoContainer->Fill("MixdZVtxPt",PJet.Pt(),PrimaryVertex[0].Z()-convvertexdz,weight);
          }
          if ((ConversionRegion=="TIB" || ConversionRegion=="TOB" || ConversionRegion=="TEC") || conv_ntracks()[convindex]==1) {
            histoContainer->Fill("MixdZPt",Photonp4[photonindex].Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
            histoContainer->Fill("MixdZEta",SuperClusterxyz[scphotonindex].Eta(),PrimaryVertex[0].Z()-superclusterdz,weight);
            histoContainer->Fill("MixdZNVtx",PrimaryVertex.size(),PrimaryVertex[0].Z()-superclusterdz,weight);
            histoContainer->Fill("MixdZVtxPt",PJet.Pt(),PrimaryVertex[0].Z()-superclusterdz,weight);
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
    double deltaeta = fabs(Photon.Eta() - GenParticlep4.Eta());
    double deltaphi = DeltaPhi(Photon.Phi(),GenParticlep4.Phi());
    double DeltaR = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);

    if (DeltaR<.1) {
      return true;
    }
  }


  return false;
}

bool PhotonPreSelection(unsigned int photonindex) {

  TLorentzVector PhotonP4 = *((TLorentzVector*) pho_p4()->At(photonindex));
  double EtCorrEcalIso = pho_ecalsumetconedr03()[photonindex] - 0.012*PhotonP4.Pt();
  double EtCorrHcalIso = pho_hcalsumetconedr03()[photonindex] - 0.005*PhotonP4.Pt();
  double EtCorrTrkIso = pho_trksumpthollowconedr03()[photonindex] - 0.002*PhotonP4.Pt();
  //double PuCorrHcalEcal = pho_ecalsumetconedr03()[photonindex] + pho_hcalsumetconedr03()[photonindex] - rho()*0.17;
  double AbsTrkIsoCIC = pho_tkiso_recvtx_030_002_0000_10_01()->at(photonindex).at(vtx_std_sel());
  
  if (pho_r9()[photonindex]<=0.9) {
    if (pho_isEB()[photonindex] && (pho_hoe()[photonindex]>0.075 || pho_sieie()[photonindex]>0.014)) return false;
    if (pho_isEE()[photonindex] && (pho_hoe()[photonindex]>0.075 || pho_sieie()[photonindex]>0.034)) return false;
    if (EtCorrEcalIso>4.0) return false;
    if (EtCorrHcalIso>4.0) return false;
    if (EtCorrTrkIso>4.0) return false;
    //if (PuCorrHcalEcal>3.0) return false;
    if (AbsTrkIsoCIC>2.8) return false;
    if (pho_trksumpthollowconedr03()[photonindex]>4.0) return false;
    if (pho_isconv()[photonindex]!=1) return false;
    return true;
  } else {
    if (pho_isEB()[photonindex] && (pho_hoe()[photonindex]>0.082 || pho_sieie()[photonindex]>0.014)) return false;
    if (pho_isEE()[photonindex] && (pho_hoe()[photonindex]>0.075 || pho_sieie()[photonindex]>0.034)) return false;
    if (EtCorrEcalIso>50.0) return false;
    if (EtCorrHcalIso>50.0) return false;
    if (EtCorrTrkIso>50.0) return false;
    //if (PuCorrHcalEcal>3.0) return false;
    if (AbsTrkIsoCIC>2.8) return false;
    if (pho_trksumpthollowconedr03()[photonindex]>4.0) return false;
    if (pho_isconv()[photonindex]!=1) return false;
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

double MomentumTrackdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot) {

  double dz = (vtx.z()-myBeamSpot.z()) - ((vtx.x()-myBeamSpot.x())*mom.x()+(vtx.y()-myBeamSpot.y())*mom.y())/mom.Perp() * mom.z()/mom.Perp();
  return dz + myBeamSpot.z();
  
}

double SuperclusterdZ(TVector3 convvtx, TVector3 superclustervtx, TVector3 beamSpot) {
  
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

pair<float,float> CombineddZ(TVector3 BeamSpot, unsigned int photonindex, int ConversionIndex){

  pair<float,float> combZ(0,0);

  pair<float,float> TrackProjection=MomentumTrackdZErr(photonindex, ConversionIndex);
  pair<float,float> SuperCluster=SuperClusterdZErr(BeamSpot, photonindex, ConversionIndex);

  //weighted avg of the two methods, where weights are based on the error
  float Z=((SuperCluster.first/(SuperCluster.second*SuperCluster.second))+(TrackProjection.first/(TrackProjection.second*TrackProjection.second)))/(1/(SuperCluster.second*SuperCluster.second)+1/(TrackProjection.second*TrackProjection.second));

  //total error sum of the two errors in quadrature
  float sigZ=sqrt((SuperCluster.second*SuperCluster.second)+(TrackProjection.second*TrackProjection.second));

  combZ.first=Z;
  combZ.second=sigZ;
  return combZ;

}

pair<float,float> MomentumTrackdZErr(unsigned int photonindex, int ConversionIndex){

  pair<float,float> ZProjection(0,0);
  TVector3 conv_momentum = *((TVector3*) conv_singleleg_momentum()->At(ConversionIndex));
  TVector3 conv_vtxpos = *((TVector3*) conv_vtx()->At(ConversionIndex));

  float theta=conv_momentum.Theta();
  float tkz=conv_vtxpos.z();
  float tkR=sqrt(conv_vtxpos.x()*conv_vtxpos.x()+conv_vtxpos.y()*conv_vtxpos.y());
  float thetaErr=conv_tk1_thetaerr()[ConversionIndex];
  float Z=tkz-tkR/tan(theta);//track projection
  float Zerr=((-1*(cos(theta)*cos(theta))/(sin(theta)* sin(theta))-1)*tkR*thetaErr); //projection error based on theta err of track theta derivative of the Zproj
  //for early tracks theta error is very small so just hard code an error 
  //that is measured by looking at the track Proj resolutin in MC
  if(tkR<39 && pho_isEB()[photonindex]) Zerr=0.234; 
  if(tkR<39 && !pho_isEB()[photonindex]) Zerr=0.341;
  ZProjection.first=Z;
  ZProjection.second=Zerr;
  return ZProjection;

}

pair<float,float> SuperClusterdZErr(TVector3 BeamSpot, unsigned int photonindex, int ConversionIndex){
  pair<float,float> Zint(0,0);
  TVector3 SuperCluster = *((TVector3*) sc_xyz()->At(pho_scind()[photonindex]));
  TVector3 SCPos(SuperCluster.x()-BeamSpot.x(),SuperCluster.y()-BeamSpot.y(),SuperCluster.z()-BeamSpot.z());

  TVector3 conv_vtxpos = *((TVector3*) conv_vtx()->At(ConversionIndex));
  TVector3 TkPos(conv_vtxpos.x()-BeamSpot.x(),conv_vtxpos.y()-BeamSpot.y(),conv_vtxpos.z()-BeamSpot.z());

  //Intersection fromt the two points:
  float R1=sqrt(SCPos.X()*SCPos.X() + SCPos.Y()*SCPos.Y()); 
  float R2=sqrt(TkPos.X()*TkPos.X() + TkPos.Y()*TkPos.Y());
  float Z1=SCPos.Z();
  float Z2=TkPos.Z();
  float slope=(Z1-Z2)/(R1-R2);
  Zint.first=Z2-R2*slope;

  //error of SL tracks of Conversion based on EB/EE and R>39 R<39 (4 cat)
  float EBLR=0.24;
  float EBHR=0.478;
  float EELR=0.416;
  float EEHR=0.888;

  if(pho_isEB()[photonindex]){
    if(conv_vtxpos.Perp()<39)Zint.second=EBLR;
    else Zint.second=EBHR;
  }
  else if (pho_isEE()[photonindex]){
    if(conv_vtxpos.Perp()<39)Zint.second=EELR;
    else Zint.second=EEHR;
  }

  return Zint;//return intersection at beamline and error in the pointing based on tracking region
}

double JacksonAngle(TLorentzVector p1, TLorentzVector p2) {

  TLorentzVector ppar = p1 + p2;

  // rotate so z axis is parent direction for p1
  TLorentzVector newp1 = dgieuler(ppar,p1);

  // boost to cm, along z axis
  TLorentzVector  pb = dgloren(newp1, ppar.Beta(), ppar.Gamma(), -1.0);

  // resulting angle is Jackson angle: return cosine theta
  return pb.Pz() / pb.P();
}

string DetectorPosition(unsigned int index) {

  string ReturnValue="";
  if (pho_isEB()[index]) ReturnValue="Barrel";
  if (pho_isEE()[index]) ReturnValue="Endcap";
  return ReturnValue;
 
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

  kFactor["PJet"]=1.3;
  kFactor["QCD"]=1.0;
  kFactor["Diphoton"]=2.0;
  kFactor["Box"]=2.0;
  kFactor["WJets"]=1.0;
  kFactor["ZJets"]=1.15;
  return kFactor;
  
}

map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight) {

  map<TString,double> WeightsMap;
  WeightsMap["None"]=1/globalweight;
  WeightsMap["PJet15to3000"]=kFactor["PJet"]*1/10566516.0*191975.0;
  WeightsMap["PJet0to15"]=kFactor["PJet"]*1/2080768.0*84200000.0;
  WeightsMap["PJet15to30"]=kFactor["PJet"]*1/2046119.0*172000.0;
  WeightsMap["PJet30to50"]=kFactor["PJet"]*1/2187260.0*16700.0;
  WeightsMap["PJet50to80"]=kFactor["PJet"]*1/2036704.0*2720.0;
  WeightsMap["PJet80to120"]=kFactor["PJet"]*1/2046637.0*447.0;
  WeightsMap["PJet120to170"]=kFactor["PJet"]*1/2088216.0*84.2;
  WeightsMap["PJet170to300"]=kFactor["PJet"]*1/2069161.0*22.6;
  WeightsMap["PJet300to470"]=kFactor["PJet"]*1/2076880.0*1.49;
  WeightsMap["PJet470to800"]=kFactor["PJet"]*1/2087212.0*0.132;
  WeightsMap["PJet800to1400"]=kFactor["PJet"]*1/2131800.0*0.00348;
  WeightsMap["PJet1400to1800"]=kFactor["PJet"]*1/2198160.0*0.0000127;
  WeightsMap["PJet1800"]=kFactor["PJet"]*1/2188301.0*0.000000294;
  WeightsMap["PJet15to3000_32PU"]=kFactor["PJet"]*1/624799.0*191975.0;
  WeightsMap["QCD20to30"]=kFactor["QCD"]*1/35721833.0*236100000.0*0.0106;
  WeightsMap["QCD30to80"]=kFactor["QCD"]*1/69968509.0*59440000.0*0.061;
  WeightsMap["QCD80to170"]=kFactor["QCD"]*1/8150672.0*898200.0*0.159;
  WeightsMap["Diphoton"]=kFactor["Diphoton"]*1/1150800.0*154.7;
  WeightsMap["Box10to25"]=kFactor["Box"]*1/528400.0*358.2;
  WeightsMap["Box25to250"]=kFactor["Box"]*1/518288.0*12.37;
  WeightsMap["Box250"]=kFactor["Box"]*1/518112.0*0.000208;
  WeightsMap["WJets"]=kFactor["WJets"]*1/78219917.0*27770.0;
  WeightsMap["ZJets"]=kFactor["ZJets"]*1/36195654.0*2475.0;
  return WeightsMap;
  
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

TLorentzVector GetP4NewVertex(unsigned int photonindex, unsigned int scindex, TVector3 Vertex) {

  TLorentzVector PhotonP4 = *((TLorentzVector*) pho_p4()->At(photonindex));
  TVector3 CaloPos = *((TVector3*) sc_xyz()->At(scindex));
  TVector3 direction = CaloPos - Vertex;
  TVector3 p = direction.Unit() * PhotonP4.Energy();
  TLorentzVector p4(p.x(),p.y(),p.z(),PhotonP4.Energy());
  return p4;
  
}

TLorentzVector dgieuler( TLorentzVector parent, TLorentzVector daughter) {

  // inverse euler angle rotations on daughter b from direction of parent a

  TVector3 a = parent.Vect(); 
  TVector3 b = daughter.Vect();

  // Parent's: costheta, sintheta, cosphi, sinphi
  double ct = a.CosTheta();
  double st = sqrt(1. - ct*ct);
  double cp  = cos(a.Phi());
  double sp  = sin(a.Phi());
 
  //bx, by, pz
  double bx = b.Px();
  double by = b.Py();
  double bz = b.Pz();

  // Euler's angle
  double rx =  ct*cp*bx + ct*sp*by - st*bz;
  double ry =  -sp*bx + cp*by ;
  double rz =  st*cp*bx + st*sp*by + ct*bz;
  double rt =  daughter.E();

  TLorentzVector r ( rx, ry, rz, rt);
  return r;
}

TLorentzVector  dgloren(TLorentzVector p, double b, double g, double ikey) {

  // does lorentz trans by beta, gamma (sense ikey)
  // on p vector (px,py,pz,e)

  double rx =  p.Px();
  double ry =  p.Py();
  double rz = g *( p.Pz() + ikey * b *p.E() );
  double rt = sqrt( p.M2() + rx*rx + ry*ry + rz*rz );

  TLorentzVector r ( rx, ry, rz, rt);
  return r;
}


int gettrackerconvindex(TVector3 Photonxyz, TVector3 BeamSpot) {

  int ReturnIndex = -1;
  double MindeltaR = 999999;
  for (int i=0; i<conv_n(); i++) {
    TVector3 ConversionRefittedPairMomentum = conv_singleleg_momentum()!=NULL && conv_ntracks()[i]==1 ? *((TVector3*) conv_singleleg_momentum()->At(i)) : *((TVector3*) conv_refitted_momentum()->At(i));
    TVector3 ConversionVertex = *((TVector3*) conv_vtx()->At(i));
    //cout << "Trying Conversion Index: " << i << " NTracks: " << conv_ntracks()[i] << " Chi2: " << conv_chi2_probability()[i] << " Conversion Refitted Pair Momentum Pt: " << ConversionRefittedPairMomentum.Pt() << endl;
    if (conv_ntracks()[i]!=2 && conv_ntracks()[i]!=1) continue;
    if (conv_ntracks()[i]==2 && (conv_chi2_probability()[i]<0.000001 || ConversionRefittedPairMomentum.Pt()<1)) continue;
    if (conv_ntracks()[i]==1 && ConversionRefittedPairMomentum.Pt()<1) continue;

    double deltaphi = DeltaPhi(Photonxyz.Phi(),ConversionVertex.Phi());
    double zfromconv = SuperclusterdZ(ConversionVertex,Photonxyz,BeamSpot);
    double deltaeta = fabs(Photonxyz.Eta() - etaTransformation(ConversionRefittedPairMomentum.Eta(),zfromconv));
    double deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    
    if (deltaR<MindeltaR) {
      MindeltaR=deltaR;
      ReturnIndex = i;
    }
    //cout << "Passing Conversion Index: " << i << " NTracks: " << conv_ntracks()[i] << " Conversion Refitted Pair Momentum Pt: " << ConversionRefittedPairMomentum.Pt() << " DeltaR: " << deltaR << endl;

  }

  if (MindeltaR<0.1) {
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

void BookCategories(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {

  TString Categories[5] = {"","_cat0","_cat1","_cat2","_cat3"};
  TString CatLabel[5] = {""," R9>0.94 Barrel"," R9<0.94 Barrel"," R9>0.94 Endcap"," R9<0.94 Endcap"};
  
  for (int cat=0; cat<5; cat++) {
    TString histnametemp = histname;
    histnametemp += Categories[cat];
    TString histtitletemp = histtitle;
    histtitletemp.ReplaceAll("Cat",CatLabel[cat]);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
  }


}

void BookHistograms(HistoContainer *histoContainer) {

  histoContainer->Add("Numvtx","Number of Primary Verticies",100,0,100);
  histoContainer->Add("NumPhoton","Number of Photons",10,0,10);
  histoContainer->Add("NumConversion","Number of Conversions",100,0,100);

  histoContainer->Add("tmva_id_mit_tiso1","tmva_id_mit_tiso1;rho_algo1;tmva_id_mit_tiso1",65,0,65,85,-10,75);
  histoContainer->Add("tmva_id_mit_tiso2","tmva_id_mit_tiso2;rho_algo1;tmva_id_mit_tiso2",65,0,65,85,-10,75);
  histoContainer->Add("tmva_id_mit_tiso1pro","tmva_id_mit_tiso1;rho_algo1;tmva_id_mit_tiso1",65,0,65,-10,75);
  histoContainer->Add("tmva_id_mit_tiso2pro","tmva_id_mit_tiso2;rho_algo1;tmva_id_mit_tiso2",65,0,65,-10,75);
  
  histoContainer->Add("PhotonPt","Pt of Photon;Pt (GeV);Counts",100,0,200);
  histoContainer->Add("PhotonEta","#eta of Photon;#eta;Counts",50,-2.5,2.5);
  histoContainer->Add("PhotonPhi","#phi of Photon;#phi;Counts",64,-3.2,3.2);
  histoContainer->Add("PhotonPtWeight","Pt of Photon;Pt (GeV);Counts",1000,0,5000);
  histoContainer->Add("PhotonEoP","EoP of Photon;EoP;Counts",100,0,3);
  histoContainer->Add("PhotonEoPGen","EoP of Gen Matched Photon;EoP;Counts",100,0,3);
  histoContainer->Add("MaxJetPt","Highest Pt Jet in the Events;Max Pt Jet;Counts",100,0,500);
  histoContainer->Add("MaxJetTrackPt","Highest Track Pt Jet in the Events;Max Pt Jet;Counts",100,0,500);
  histoContainer->Add("PairMass","Mass of Photon-Jet System;Pt (GeV);Counts",100,0,600);
  histoContainer->Add("PairPt","Pt of Photon-Jet System;Pt (GeV);Counts",60,0,120);
  histoContainer->Add("PairPtbusted","Pt of Photon-Jet System;Pt (GeV);Counts",60,0,120);
  histoContainer->Add("PhotonJetJacksonAngle","Jackson Angle between Jet and Photon;Jackson Angle;Counts",50,-1,1);

  histoContainer->Add("ConvPt","Pt of Conversion;Pt (GeV);Counts",100,0,200);
  histoContainer->Add("ConvPairPt","Pt of Conversion Refitted Tracks;Pt (GeV);Counts",100,0,200);
  histoContainer->Add("ConvEta","#eta of Conversion;#eta;Counts",50,-2.5,2.5);
  histoContainer->Add("ConvPhi","#phi of Conversion;#phi;Counts",64,-3.2,3.2);
  histoContainer->Add("ConvR","Radius of Conversion;Radius (cm);Counts",100,0,100);

  histoContainer->Add("ConvRTracker","Radius of Conversion;Radius (cm) |z|<26cm;Conversions/0.2cm",300,0,60);
  histoContainer->Add("ConvZPixel","Z Position of Conversion;Radius (cm) 3.5cm<R<19cm;Converisons/1.0cm",120,-60,60);

  histoContainer->Add("ConvR-0.5","Radius of Conversion -0.5<=#eta<0.0;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR-0.8","Radius of Conversion -0.8<=#eta<-0.5;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR-1.1","Radius of Conversion -1.1<=#eta<-0.8;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR-1.3","Radius of Conversion -1.3<=#eta<-1.1;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR-1.5","Radius of Conversion -1.5<=#eta<-1.3;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR0.0","Radius of Conversion 0.0<=#eta<0.5;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR0.5","Radius of Conversion 0.5<=#eta<0.8;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR0.8","Radius of Conversion 0.8<=#eta<1.1;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR1.1","Radius of Conversion 1.1<=#eta<1.3;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR1.3","Radius of Conversion 1.3<=#eta<1.5;#eta;Counts",100,0,100);

  histoContainer->Add("ConvRMod4","Radius of Conversion 1.1<=#eta<1.5;#eta;Counts",100,0,100);
  histoContainer->Add("ConvR-Mod4","Radius of Conversion -1.1<=#eta<-1.5;#eta;Counts",100,0,100);

  histoContainer->Add("ConvRFine-0.5","Radius of Conversion -0.5<=#eta<0.0;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine-0.8","Radius of Conversion -0.8<=#eta<-0.5;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine-1.1","Radius of Conversion -1.1<=#eta<-0.8;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine-1.3","Radius of Conversion -1.3<=#eta<-1.1;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine-1.5","Radius of Conversion -1.5<=#eta<-1.3;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine0.0","Radius of Conversion 0.0<=#eta<0.5;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine0.5","Radius of Conversion 0.5<=#eta<0.8;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine0.8","Radius of Conversion 0.8<=#eta<1.1;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine1.1","Radius of Conversion 1.1<=#eta<1.3;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine1.3","Radius of Conversion 1.3<=#eta<1.5;#eta;Counts",100,0,25);

  histoContainer->Add("ConvRFineMod4","Radius of Conversion 1.1<=#eta<1.5;#eta;Counts",100,0,25);
  histoContainer->Add("ConvRFine-Mod4","Radius of Conversion -1.1<=#eta<-1.5;#eta;Counts",100,0,25);

  histoContainer->Add("SCConvR-0.5","Radius of Conversion -0.5<=SC_{#eta}<0.0;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR-0.8","Radius of Conversion -0.8<=SC_{#eta}<-0.5;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR-1.1","Radius of Conversion -1.1<=SC_{#eta}<-0.8;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR-1.3","Radius of Conversion -1.3<=SC_{#eta}<-1.1;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR-1.5","Radius of Conversion -1.5<=SC_{#eta}<-1.3;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR0.0","Radius of Conversion 0.0<=SC_{#eta}<0.5;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR0.5","Radius of Conversion 0.5<=SC_{#eta}<0.8;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR0.8","Radius of Conversion 0.8<=SC_{#eta}<1.1;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR1.1","Radius of Conversion 1.1<=SC_{#eta}<1.3;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR1.3","Radius of Conversion 1.3<=SC_{#eta}<1.5;SC_{#eta};Counts",100,0,100);

  histoContainer->Add("SCConvRMod4","Radius of Conversion 1.1<=#eta<1.5;SC_{#eta};Counts",100,0,100);
  histoContainer->Add("SCConvR-Mod4","Radius of Conversion -1.1<=#eta<-1.5;SC_{#eta};Counts",100,0,100);

  histoContainer->Add("SCConvRFine-0.5","Radius of Conversion -0.5<=SC_{#eta}<0.0;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine-0.8","Radius of Conversion -0.8<=SC_{#eta}<-0.5;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine-1.1","Radius of Conversion -1.1<=SC_{#eta}<-0.8;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine-1.3","Radius of Conversion -1.3<=SC_{#eta}<-1.1;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine-1.5","Radius of Conversion -1.5<=SC_{#eta}<-1.3;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine0.0","Radius of Conversion 0.0<=SC_{#eta}<0.5;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine0.5","Radius of Conversion 0.5<=SC_{#eta}<0.8;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine0.8","Radius of Conversion 0.8<=SC_{#eta}<1.1;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine1.1","Radius of Conversion 1.1<=SC_{#eta}<1.3;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine1.3","Radius of Conversion 1.3<=SC_{#eta}<1.5;SC_{#eta};Counts",100,0,25);

  histoContainer->Add("SCConvRFineMod4","Radius of Conversion 1.1<=#eta<1.5;SC_{#eta};Counts",100,0,25);
  histoContainer->Add("SCConvRFine-Mod4","Radius of Conversion -1.1<=#eta<-1.5;SC_{#eta};Counts",100,0,25);

  histoContainer->Add("ConvXvsZ","X of Conversion vs Z;Z (cm);X (cm)",400,-200,200,200,-100,100);
  histoContainer->Add("ConvYvsZ","Y of Conversion vs Z;Z (cm);Y (cm)",400,-200,200,200,-100,100);

  histoContainer->Add("ConvRvsZ","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice1","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice2","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice3","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice4","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice5","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice6","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice7","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice8","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);
  histoContainer->Add("ConvRvsZslice9","R of Conversion vs Z;Z (cm);R (cm)",200,-200,200,100,-100,100);

  histoContainer->Add("ConvRvsEta","R of Conversion vs Eta;#eta;Radius (cm)",64,-3.2,3.2,100,0,100);
  histoContainer->Add("ConvRvsPhi","R of Conversion vs phi;#phi;Radius (cm)",60,-3.0,3.0,100,0,100);

  histoContainer->Add("sumpt2","Sum Pt^2 of Primary Vertex;Sum Pt^2;Counts",100,0,1000);
  histoContainer->Add("logsumpt2","Log Sum Pt^2 of Primary Vertex;log(Sum Pt^2);Counts",60,-5,15);
  histoContainer->Add("logsumpt2All","Log Sum Pt^2 of Primary Vertex;log(Sum Pt^2);Counts",60,-5,15);
  histoContainer->Add("ptasymm","Pt Asymmetry of Primary Vertex;ptasymm;Counts",50,-1,1);
  histoContainer->Add("ptbal","Pt Balance of Primary Vertex;ptbal;Counts",60,-50,100);
  histoContainer->Add("pulltoconv","Pull to Conversion;Pull;Counts",40,0,10);
  histoContainer->Add("limpulltoconv","Pull to Conversion;Pull;Counts",40,0,10);
  histoContainer->Add("MVAValue","Lead BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueFirst","First BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueSecond","Second BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueThird","Third BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVARes","Resolution MVA Value;Per Event ;Counts",50,-1,1);

  histoContainer->Add("sumpt2bad","Sum Pt^2 of other Verticies;Sum Pt^2;Counts",100,0,1000);
  histoContainer->Add("logsumpt2bad","Log Sum Pt^2 of other Verticies;log(Sum Pt^2);Counts",60,-5,15);
  histoContainer->Add("ptasymmbad","Pt Asymmetry of other Verticies;ptasymm;Counts",50,-1,1);
  histoContainer->Add("ptbalbad","Pt Balance of other Verticies;ptbal;Counts",60,-50,100);
  histoContainer->Add("pulltoconvbad","Pull to Converion of other Verticies;Pull;Counts",40,0,10);
  histoContainer->Add("limpulltoconvbad","Pull to Converison of other Verticies;Pull;Counts",40,0,10);
  histoContainer->Add("MVAValuebad","Output BDT value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAResbad","Resolution MVA Value;Per Event ;Counts",50,-1,1);
  histoContainer->Add("MVAResbadsim","Resolution MVA Value;Per Event ;Counts",50,-1,1);
  histoContainer->Add("MVAResgoodbadsim","Resolution MVA Value;Per Event ;Counts",50,-1,1);
  histoContainer->Add("MVAResfunky","Resolution MVA Value;Per Event ;Counts",50,-1,1);
  
  histoContainer->Add("sumpt2busted","Sum Pt^2 of Primary Vertex;Sum Pt^2;Counts",100,0,1000);
  histoContainer->Add("logsumpt2busted","Log Sum Pt^2 of Primary Vertex;log(Sum Pt^2);Counts",60,-5,15);
  histoContainer->Add("ptasymmbusted","Pt Asymmetry of Primary Vertex;ptasymm;Counts",50,-1,1);
  histoContainer->Add("ptbalbusted","Pt Balance of Primary Vertex;ptbal;Counts",60,-50,100);
  histoContainer->Add("pulltoconvbusted","Pull to Conversion;Pull;Counts",40,0,10);
  histoContainer->Add("limpulltoconvbusted","Pull to Conversion;Pull;Counts",40,0,10);
  histoContainer->Add("MVAValuebusted","Lead BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueFirstbusted","First BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueSecondbusted","Second BDT Output value;BDT Value;Counts",50,-1,1);
  histoContainer->Add("MVAValueThirdbusted","Third BDT Output value;BDT Value;Counts",50,-1,1);
  
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
  BookTrackerdZ(histoContainer,"SuperdZEff","visual");
  BookTrackerdZ(histoContainer,"SuperdZRes","resolution");
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
  BookTrackerdZ(histoContainer,"ConvdZEff","visual");
  BookTrackerdZ(histoContainer,"ConvdZRes","resolution");
  BookTrackerdZ(histoContainer,"ConvdZsim");
  BookTrackerdZ(histoContainer,"ConvdZCompairison","twoD");
  BookEtadZ(histoContainer,"ConvdZEta");
  BookPtdZ(histoContainer,"ConvdZPt");
  histoContainer->Add("ConvdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("ConvdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Conversion Vertex vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);

  histoContainer->Add("RishiSuperZPV","Z of Primary Vertex from Rishi's SuperCluster Projection Vertex;Z (cm);Counts",100,-20,20);
  histoContainer->Add("RishiSuperdZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  histoContainer->Add("RishiSuperdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiSuperZPV","Z of Primary Vertex from Rishi's SuperCluster Projection Vertex: region;Z (cm);Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"RishiSuperdZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiSuperdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookTrackerdZ(histoContainer,"RishiSuperdZ");
  BookTrackerdZ(histoContainer,"RishiSuperdZEff","visual");
  BookTrackerdZ(histoContainer,"RishiSuperdZRes","resolution");
  BookTrackerdZ(histoContainer,"RishiSuperdZsim");
  BookTrackerdZ(histoContainer,"RishiSuperdZCompairison","twoD");
  BookEtadZ(histoContainer,"RishiSuperdZEta");
  BookPtdZ(histoContainer,"RishiSuperdZPt");
  histoContainer->Add("RishiSuperdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("RishiSuperdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);
  
  histoContainer->Add("RishiTrackZPV","Z of Primary Vertex from Rishi's Track Projection Vertex;Z (cm);Counts",100,-20,20);
  histoContainer->Add("RishiTrackdZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's Track Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  histoContainer->Add("RishiTrackdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's Track Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiTrackZPV","Z of Primary Vertex from Rishi's Track Projection Vertex: region;Z (cm);Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"RishiTrackdZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's Track Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiTrackdZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's Track Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookTrackerdZ(histoContainer,"RishiTrackdZ");
  BookTrackerdZ(histoContainer,"RishiTrackdZEff","visual");
  BookTrackerdZ(histoContainer,"RishiTrackdZRes","resolution");
  BookTrackerdZ(histoContainer,"RishiTrackdZsim");
  BookTrackerdZ(histoContainer,"RishiTrackdZCompairison","twoD");
  BookEtadZ(histoContainer,"RishiTrackdZEta");
  BookPtdZ(histoContainer,"RishiTrackdZPt");
  histoContainer->Add("RishiTrackdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's Track Projection Vertex vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("RishiTrackdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's Track Projection Vertex vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);

  histoContainer->Add("RishiCombinedZPV","Z of Primary Vertex from Rishi's SuperCluster Projection Vertex;Z (cm);Counts",100,-20,20);
  histoContainer->Add("RishiCombineddZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  histoContainer->Add("RishiCombineddZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiCombinedZPV","Z of Primary Vertex from Rishi's SuperCluster Projection Vertex: region;Z (cm);Counts",100,-20,20);
  BookBarrelAndEndcap(histoContainer,"RishiCombineddZ","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookBarrelAndEndcap(histoContainer,"RishiCombineddZsim","#DeltaZ of Sim Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex: region;#DeltaZ (cm);Counts",100,-2,2);
  BookTrackerdZ(histoContainer,"RishiCombineddZ");
  BookTrackerdZ(histoContainer,"RishiCombineddZEff","visual");
  BookTrackerdZ(histoContainer,"RishiCombineddZRes","resolution");
  BookTrackerdZ(histoContainer,"RishiCombineddZsim");
  BookEtadZ(histoContainer,"RishiCombineddZEta");
  BookPtdZ(histoContainer,"RishiCombineddZPt");
  histoContainer->Add("RishiCombineddZPt","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("RishiCombineddZEta","#DeltaZ of Primary Vertex and Primary Vertex from Rishi's SuperCluster Projection Vertex vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);
  
  histoContainer->Add("MixdZPt","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Photon Pt;Pt (GeV);#DeltaZ (cm)",20,0,200,100,-5,5);
  histoContainer->Add("MixdZEta","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Photon #eta;#eta;#DeltaZ (cm)",25,0,2.5,100,-5,5);
  histoContainer->Add("MixdZNVtx","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Number of Verticies;Number of Verticies;#DeltaZ (cm)",50,0,50,100,-5,5);
  histoContainer->Add("MixdZVtxPt","#DeltaZ of Primary Vertex and Primary Vertex from Mixed Method vs Pt of Photon-Jet System;PJet Pt (GeV);#DeltaZ (cm)",40,0,200,100,-5,5);

  histoContainer->Add("MVAdZ","#DeltaZ of Primary Vertex and Primary Vertex MVA;#DeltaZ (cm);Counts",100,-5,5);
  histoContainer->Add("MVAdZFirst","#DeltaZ of Primary Vertex and First Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZSecond","#DeltaZ of Primary Vertex and Secondary Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZThird","#DeltaZ of Primary Vertex and Ternary Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZFirstbusted","#DeltaZ of Primary Vertex and First Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZSecondbusted","#DeltaZ of Primary Vertex and Secondary Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZThirdbusted","#DeltaZ of Primary Vertex and Ternary Vertex MVA;#DeltaZ (cm);Counts",40,-10,10);
  histoContainer->Add("MVAdZNVtx","#DeltaZ of Primary Vertex and Primary Vertex MVA vs Number of Verticies;Number of Verticies;#DeltaZ (cm)",50,0,50,100,-5,5);
  histoContainer->Add("MVAdZNVtx_short","#DeltaZ of Primary Vertex and Primary Vertex MVA vs Number of Verticies;Number of Verticies;#DeltaZ (cm)",30,0,30,100,-5,5);
  histoContainer->Add("MVAdZVtxPt","#DeltaZ of Primary Vertex and Primary Vertex MVA vs Pt of Photon-Jet System;PJet Pt (GeV);#DeltaZ (cm)",40,0,200,100,-5,5);

  histoContainer->Add("MVAdZsim","#DeltaZ of the SimVertex and Primary Vertex MVA;#DeltaZ (cm);Counts",100,-5,5);
  histoContainer->Add("MVAdZsimbad","#DeltaZ of the SimVertex and Primary Vertex MVA;#DeltaZ (cm);Counts",100,-5,5);
  histoContainer->Add("MVAdZsimfunky","#DeltaZ of the SimVertex and Primary Vertex MVA;#DeltaZ (cm);Counts",100,-5,5);

  histoContainer->Add("MVAdZNVtxAll","#DeltaZ of Primary Vertex and Primary Vertex MVA vs Number of Verticies;Number of Verticies;#DeltaZ (cm)",50,0,50,100,-5,5);
  histoContainer->Add("MVAdZVtxPtAll","#DeltaZ of Primary Vertex and Primary Vertex MVA vs Pt of Photon-Jet System;PJet Pt (GeV);#DeltaZ (cm)",40,0,200,100,-5,5);
  
  histoContainer->Add("ConvdZvsR","#deltaZ between the Z of the Primary Vertex from Conversion and the Primary Vertex versus R of Conversion;R (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("ConvdZvsZ","#deltaZ between the Z of the Primary Vertex from Conversion and the Primary Vertex versus Z of Conversion;Z (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-220,220,100, -5, 5);

  histoContainer->Add("SuperdZvsR","#deltaZ between the Z of the Primary Vertex from Conversion and the Primary Vertex versus R of Conversion;R (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,0,100,100, -5, 5);
  histoContainer->Add("SuperdZvsZ","#deltaZ between the Z of the Primary Vertex from Conversion and the Primary Vertex versus Z of Conversion;Z (cm);#deltaZ of Primary Vertex from Conversion (cm)",100,-220,220,100, -5, 5);

  BookCategories(histoContainer,"isosumoet","IsoSumEt:Cat;IsoSumEt;Counts",50,0,50);
  BookCategories(histoContainer,"isosumoetbad","IsoSumBadEt:Cat;IsoSumBadEt;Counts",50,0,50);
  BookCategories(histoContainer,"trkisooet","TrackIsoEt:Cat;TrackIsoEt;Counts",50,0,50);
  BookCategories(histoContainer,"sieie","#sigmai#etai#eta:Cat;#sigmai#etai#eta;Counts",100,0,0.04);
  BookCategories(histoContainer,"hoe","HoverE:Cat;HoverE;Counts",100,0,0.1);
  BookCategories(histoContainer,"r9","R9:Cat;R9;Counts",110,0,1.1);
  BookCategories(histoContainer,"drtotk_25_99","#DeltaR to Track:Cat;drtotk_25_99;Counts",50,0,0.5);
  BookCategories(histoContainer,"pho_cic4cutlevel_lead","CiC Cut Level Lead Photon:Cat;pho_cic4cutlevel_lead;Counts",9,-0.5,8.5);
  
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
    if (options.Contains("visual")) {
      bins=50;
      if ( i==0 ) {
        lowerlimit = -0.2;
        upperlimit = 0.2;
      }
      if ( i==1 || i==4 ) {
        lowerlimit = -3;
        upperlimit = 3;
      }
      if ( i==2 || i==5 ) {
        lowerlimit = -5;
        upperlimit = 5;
      }
      if ( i==3 ) {
        lowerlimit = -0.5;
        upperlimit = 0.5;
      }
    }
    if (options.Contains("resolution")) {
      if ( i==0 ) {
        bins=10000;
        lowerlimit = -20;
        upperlimit = 20;
      }
      if ( i==1 || i==4 ) {
        bins=1600;
        lowerlimit = -40;
        upperlimit = 40;
      }
      if ( i==2 || i==5 ) {
        bins=800;
        lowerlimit = -40;
        upperlimit = 40;
      }
      if ( i==3 ) {
        bins=4000;
        lowerlimit = -20;
        upperlimit = 20;
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

void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  if (inputstring.Contains("Run2011A") && !inputstring.Contains("Summer11")) {
    inputfilelist.push_back(pair<string,int> ("Run2011A.root",6));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_0.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_0_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_1_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_2.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011A_2_1.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2011B") && !inputstring.Contains("Summer11")) {
    inputfilelist.push_back(pair<string,int> ("Run2011B.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011B_0.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011B_0_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011B_0_2.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Run2011B_1.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2011A_Summer11")) {
    inputfilelist.push_back(pair<string,int> ("Run2011A_Summer11.root",5));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet/July5thReReco.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet/July5thReReco_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet/July5thReReco_2.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet/Aug5thReReco.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet/Oct3rdReReco.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2011B_Summer11")) {
    inputfilelist.push_back(pair<string,int> ("Run2011B_Summer11.root",4));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Run2011B.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Run2011B_1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Run2011B_2.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2011BSummer11Test")) {
    inputfilelist.push_back(pair<string,int> ("Run2011B_Summer11.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/b/drberry/PhotonPlusJet/Run2011B_2.root",WeightsMap["None"]));
  } 
  if (inputstring.Contains("PJet_Summer11")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetMC.root",10));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-0to15_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet0to15"]));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-15to30_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet15to30"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-30to50_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet30to50"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-50to80_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet50to80"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-80to120_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet80to120"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-120to170_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet120to170"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-170to300_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet170to300"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-300to470_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet300to470"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-470to800_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet470to800"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-800to1400_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet800to1400"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1400to1800_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet1400to1800"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-1800_TuneZ2.root",kFactor["PJet"]*WeightsMap["PJet1800"]));
  }
  if (inputstring.Contains("PJet15to3000")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetMC.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/G_Pt-15to3000_1.root",kFactor["PJet"]*WeightsMap["PJet15to3000"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/G_Pt-15to3000_2.root",kFactor["PJet"]*WeightsMap["PJet15to3000"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/G_Pt-15to3000_3.root",kFactor["PJet"]*WeightsMap["PJet15to3000"]));
  }
  if (inputstring.Contains("PJet15to3000_32PU")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJetMC_32PU.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/G_Pt-15to3000_32PU.root",kFactor["PJet"]*WeightsMap["PJet15to3000_32PU"]));
  }
  if (inputstring.Contains("QCD")) {
    inputfilelist.push_back(pair<string,int> ("QCDEMEnriched.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/QCD_Pt-20to30_EMEnriched.root",kFactor["QCD"]*WeightsMap["QCD20to30"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/QCD_Pt-30to80_EMEnriched.root",kFactor["QCD"]*WeightsMap["QCD30to80"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/QCD_Pt-80to170_EMEnriched.root",kFactor["QCD"]*WeightsMap["QCD80to170"]));
  }
  if (inputstring.Contains("Diphoton")) {
    inputfilelist.push_back(pair<string,int> ("Diphoton.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/DiPhotonsJets.root",kFactor["Diphoton"]*WeightsMap["Diphoton"]));
  }
  if (inputstring.Contains("Box")) {
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Box10to25.root",kFactor["Box"]*WeightsMap["Box10to25"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Box25to250.root",kFactor["Box"]*WeightsMap["Box25to250"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Box250.root",kFactor["Box"]*WeightsMap["Box250"]));
  }
  if (inputstring.Contains("WJets")) {
    inputfilelist.push_back(pair<string,int> ("WJets.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/WJets.root",kFactor["WJets"]*WeightsMap["WJets"]));
  }
  if (inputstring.Contains("ZJets")) {
    inputfilelist.push_back(pair<string,int> ("ZJets.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/ZJets.root",kFactor["ZJets"]*WeightsMap["ZJets"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/ZJets_1.root",kFactor["ZJets"]*WeightsMap["ZJets"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/ZJets_2.root",kFactor["ZJets"]*WeightsMap["ZJets"]));
  }
  if (inputstring.Contains("Signal")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_120GeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc3/b/drberry/PhotonPlusJet_S6/Higgs_120GeV.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("HighPU")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_145GeV_HighPU.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/HiggsGluGluToHToGG_M-145_7TeV_HighPU.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("145GeV")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_145GeV.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/HiggsGluGluToHToGG_M-145_7TeV.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("125GeV52XSL")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_125GeV52XSL.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-125_8TeV_52X_SingleLeg.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("125GeV52XDL")) {
    inputfilelist.push_back(pair<string,int> ("Higgs_125GeV52XDL.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/GluGluToHToGG_M-125_8TeV_52X_DoubleLeg.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Test") && !inputstring.Contains("TestMC") && !inputstring.Contains("Run2011")) {
    inputfilelist.push_back(pair<string,int> ("Test.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/DataTest.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("TestMC")) {
    inputfilelist.push_back(pair<string,int> ("TestMC.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/G_Pt-30to50_TuneZ2.root",WeightsMap["None"]));
  }

}

void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  string outfile;
  if (infile.rfind("/")!=string::npos) outfile=infile.substr(infile.rfind("/")+1);
  else outfile=infile;
  
  inputfilelist.push_back(pair<string,int> (outfile,1));
  if (inputstring.Contains("Run2011A")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("Run2011B")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("PJet15to3000")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet15to3000"]));
  if (inputstring.Contains("PJet0to15")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet0to15"]));
  if (inputstring.Contains("PJet15to30")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet15to30"]));
  if (inputstring.Contains("PJet30to50")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet30to50"]));
  if (inputstring.Contains("PJet50to80")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet50to80"]));
  if (inputstring.Contains("PJet80to120")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet80to120"]));
  if (inputstring.Contains("PJet120to170")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet120to170"]));
  if (inputstring.Contains("PJet170to300")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet170to300"]));
  if (inputstring.Contains("PJet300to470")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet300to470"]));
  if (inputstring.Contains("PJet470to800")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet470to800"]));
  if (inputstring.Contains("PJet800to1400")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet800to1400"]));
  if (inputstring.Contains("PJet1400to1800")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet1400to1800"]));
  if (inputstring.Contains("PJet1800")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet1800"]));
  if (inputstring.Contains("PJet15to3000_32PU")) inputvector.push_back(pair<string,float> (infile,kFactor["PJet"]*WeightsMap["PJet15to3000_32PU"]));
  if (inputstring.Contains("QCD20to30")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD20to30"]));
  if (inputstring.Contains("QCD30to80")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD30to80"]));
  if (inputstring.Contains("QCD80to170")) inputvector.push_back(pair<string,float> (infile,kFactor["QCD"]*WeightsMap["QCD80to170"]));
  if (inputstring.Contains("Diphoton")) inputvector.push_back(pair<string,float> (infile,kFactor["Diphoton"]*WeightsMap["Diphoton"]));
  if (inputstring.Contains("Box10to25")) inputvector.push_back(pair<string,float> (infile,kFactor["Box"]*WeightsMap["Box10to25"]));
  if (inputstring.Contains("Box25to250")) inputvector.push_back(pair<string,float> (infile,kFactor["Box"]*WeightsMap["Box25to250"]));
  if (inputstring.Contains("Box250")) inputvector.push_back(pair<string,float> (infile,kFactor["Box"]*WeightsMap["Box250"]));
  if (inputstring.Contains("WJets")) inputvector.push_back(pair<string,float> (infile,kFactor["WJets"]*WeightsMap["WJets"]));
  if (inputstring.Contains("ZJets")) inputvector.push_back(pair<string,float> (infile,kFactor["ZJets"]*WeightsMap["ZJets"]));
  if (inputstring.Contains("Signal")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("HighPU")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("145GeV")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("120GeV52X")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("Test") && !inputstring.Contains("TestMC")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("TestMC")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
    
}

void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap) {

  if (inputstring.Contains("PJet15to3000")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/PhotonPlusJet.h"
  } else if (inputstring.Contains("PJet") && (inputstring.Contains("to") || inputstring.Contains("1800"))) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/PhotonPlusJet_Summer11.h"
  } else if (inputstring.Contains("QCD")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/QCDEMEnriched.h"
  } else if (inputstring.Contains("Box")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Box.h"
  } else if (inputstring.Contains("Diphoton")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Diphoton.h"
  } else if (inputstring.Contains("WJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/WJets.h"
  } else if (inputstring.Contains("ZJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZJets.h"
//   } else if (inputstring.Contains("125GeV52X")) {
//     #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/125GeV_52X.h"
//   } else if (inputstring.Contains("125GeV42X")) {
//     #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/125GeV_42X.h"
  } else {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Dummy.h"
  }

}

void MakeEtWeights(TString inputstring, map<int,double> &EtMap) {

  if (inputstring.Contains("PJet15to3000")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/PhotonPlusJet.h"
  } else if (inputstring.Contains("PJet") && (inputstring.Contains("to") || inputstring.Contains("1800"))) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/PhotonPlusJet_Summer11.h"
  } else if (inputstring.Contains("PJet15to3000_32PU")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/PhotonPlusJet_32PU.h"
  } else if (inputstring.Contains("QCD")) {
    #include "ND_Hto2Photons/TreeReaders/interface/EtWeights/QCDEMEnriched.h"
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
