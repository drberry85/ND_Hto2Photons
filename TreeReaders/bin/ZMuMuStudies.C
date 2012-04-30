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
#include "ND_Hto2Photons/TreeReaders/interface/sdaReaderFast_ZMuMu.h"
#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.cc"

using namespace std;

bool GenMatch(TVector3 Photon);
bool PhotonPreSelection(unsigned int PhotonIndex);
double etaTransformation(double EtaParticle, double Zvertex);
double DeltaPhi(double Phi1, double Phi2);
double FindNewdZ(TVector3 vtx, TVector3 mom, TVector3 myBeamSpot);
double FindNewZConvLinear(TVector3 convvtx, TVector3 superclustervtx, TVector3 primaryvertex);
double JacksonAngle(TLorentzVector p1, TLorentzVector p2);
string DetectorPosition(unsigned int index);
string GetPhotonCat(unsigned int index);
map<TString,double> GetkFactor();
map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight);
template <class type> string makestring(type value);
TLorentzVector dgieuler(TLorentzVector parent, TLorentzVector daughter);
TLorentzVector dgloren(TLorentzVector p, double b, double g, double ikey);
int gettrackerconvindex(TVector3 Photonxyz, TVector3 BeamSpot);
void BookBarrelAndEndcap(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit);
void BookHistograms(HistoContainer *histoContainer);
void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap);
void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap);
void ProgressBar(int &percent, double estimate);

int main(int argc, char * input[]) {

  gROOT->ProcessLine(".L $CMSSW_BASE/src/ND_Hto2Photons/TreeReaders/interface/link_def.h+");

  TString InputArgs(input[1]);
  bool background = false;
  bool bar = false;
  bool debug = false;
  bool fake = false;
  bool higgs = false;
  bool highpt = false;
  bool presel = false;
  bool mc = true;
  bool nocuts = false;
  bool  nearMuonDRcut = false; 
  bool onevertex = false;
  bool onephoton = false;
  bool pixel = false;
  bool prompt = false;
  bool trigger = false;
  bool unweighted = false;
  bool usesimvertex = false;
  float globalweight = 5025.213;

  int FirstFileNum = 0;

  double RealZMass = 91.1876;
  
  if (InputArgs.Contains("Unweighted")) unweighted=true;
  if (InputArgs.Contains("Bar")) bar=true;
  if (InputArgs.Contains("Background")) background=true;
  if (InputArgs.Contains("Debug")) debug=true;
  if (InputArgs.Contains("Fake")) fake=true;
  if (InputArgs.Contains("Higgs")) higgs=true;
  if (InputArgs.Contains("HighPt")) highpt=true;
  if (InputArgs.Contains("PreSel")) presel=true;
  if (InputArgs.Contains("NearMuonDRCut")) nearMuonDRcut=true;
  if (InputArgs.Contains("NoCuts")) nocuts=true;
  if (InputArgs.Contains("OneVertex")) onevertex=true;
  if (InputArgs.Contains("OnePhoton")) onephoton=true;
  if (InputArgs.Contains("Pixel")) pixel=true;
  if (InputArgs.Contains("Prompt")) prompt=true;
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
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Data, ZMuMu, or TTJets." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }

  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    TString outfilename = "";
    if (argc>3) {
      outfilename = TString(input[3]);
      outfilename += "/";
    }
    outfilename += "ZMuMu_";
    outfilename += itFilePair->first;
    if (argc==2) {
      if (unweighted) outfilename.ReplaceAll(".root","_Unweighted.root");
      if (onevertex) outfilename.ReplaceAll(".root","_OneVertex.root");
      if (usesimvertex) outfilename.ReplaceAll(".root","_SimVertex.root");
      if (background) outfilename.ReplaceAll(".root","_Background.root");
      if (nocuts) outfilename.ReplaceAll(".root","_NoCuts.root");
      if (presel) outfilename.ReplaceAll(".root","_PhoPresel.root");
      if (nearMuonDRcut ) outfilename.ReplaceAll(".root","_NearMuonDRCut.root");
      if (onephoton) outfilename.ReplaceAll(".root","_OnePhoton.root");
      if (fake) outfilename.ReplaceAll(".root","_Fake.root");
      if (highpt) outfilename.ReplaceAll(".root","_HighPt.root");
      if (prompt) outfilename.ReplaceAll(".root","_Prompt.root");
      if (pixel) outfilename.ReplaceAll(".root","_Pixel.root");
    }
    
    TFile* outfile = new TFile(outfilename.Data(),"RECREATE");
    outfile->cd();
    cout << "\n" << outfilename << " created." << endl;

    HistoContainer* histoContainer;
    histoContainer = new HistoContainer();

    BookHistograms(histoContainer);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float fileweight = filesAndWeights[itFile].second * globalweight;
      if (unweighted) fileweight=1;
      
      TChain* filechain = new TChain("event");
      filechain->Add(file.c_str());

      if (debug) cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

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
        if (ientry==0 && gp_n()==0) mc=false;

        if (ientry==0 && trigger) {
          cout << "Trigger Names:" << endl; 
          for (unsigned int j=0; j<(unsigned int) hlt_path_names_HLT1()->size(); j++) {
            TString TriggerName(hlt_path_names_HLT1()->at(j));
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
        vector<TLorentzVector> MuonP4;
        map<double,unsigned int> MuonPtMap;
        map<double,unsigned int> PhotonPtMap;
        map<double,vector<unsigned int> > TripletMassIndex;
        if (pho_n()<1) continue;
        for (unsigned int j=0; j!=(unsigned int) vtx_std_xyz()->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) vtx_std_xyz()->At(j)));
        float PUWeight = PileUpMap[PrimaryVertex.size()];
        float weight = fileweight*PUWeight;
        if (debug) cout << "Number of Vertices: " << PrimaryVertex.size() << " PileUp Weight: " << PUWeight << endl;
        for (unsigned int j=0; j<(unsigned int) sc_xyz()->GetSize(); j++) {
          SuperClusterxyz.push_back(*((TVector3*) sc_xyz()->At(j)));
          SuperClusterp4.push_back(*((TLorentzVector*) sc_p4()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) bs_xyz()->GetSize(); j++) BeamSpot.push_back(*((TVector3*) bs_xyz()->At(j)));
        for (unsigned int j=0; j<(unsigned int) pho_n(); j++) {
          Photonp4.push_back(*((TLorentzVector*) pho_p4()->At(j)));
          PhotonPtMap[((TLorentzVector*) pho_p4()->At(j))->Pt()]=j;
          Photonxyz.push_back(*((TVector3*) pho_calopos()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) conv_n(); j++) {
          ConversionVertex.push_back(*((TVector3*) conv_vtx()->At(j))+DetectorOffset);
          ConversionRefittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum()->At(j)));
        }
        for (unsigned int j=0; j<(unsigned int) mu_glo_n(); j++) {
          TLorentzVector Muon_p4 = *((TLorentzVector*) mu_glo_p4()->At(j));
          MuonP4.push_back(Muon_p4);
          if (mu_glo_type()[j]<1100) continue;
          if (mu_glo_chi2()[j]/mu_glo_dof()[j]>10.0) continue;
          if (mu_glo_validhits()[j]<11) continue;
          if (mu_glo_pixelhits()[j]<1) continue;
          //if (fabs(mu_glo_d0()[j])>0.2) continue;
          if (mu_glo_tkiso03()[j]>3.0) continue;
          //if (debug) cout << "Has GSF Track: " << mu_glo_hasgsftrack()[j] << endl;
          //if (mu_glo_hasgsftrack()[j]) continue;
          if (Muon_p4.Pt()<10.0) continue;
          if (fabs(Muon_p4.Eta())>2.4) continue;
          MuonPtMap[Muon_p4.Pt()]=j;
        }
        map<double,unsigned int>::reverse_iterator PhotonIterator;
        unsigned int LeadMuonIndex,SubLeadMuonIndex,ZMuMuPhotonIndex,LeadPhotonIndex,SubLeadPhotonIndex = 0;
        unsigned int NearMuonIndex,FarMuonIndex =0 ;
        if (!higgs) {
          if (MuonPtMap.size()<2) continue;
          if (debug) cout << "MuonPt Map Filled: " << MuonPtMap.size() << endl;
          for (unsigned int j=0; j<(unsigned int) mu_glo_n(); j++) {
            TLorentzVector LeadMuon_p4 = *((TLorentzVector*) mu_glo_p4()->At(j));
            if (debug) cout << "Lead Trial Muon: " << j << " Pt: " << LeadMuon_p4.Pt() << endl;
            if (mu_glo_type()[j]<1100) continue;
            if (mu_glo_chi2()[j]/mu_glo_dof()[j]>10.0) continue;
            if (mu_glo_validhits()[j]<11) continue;
            if (mu_glo_pixelhits()[j]<1) continue;
            if (mu_glo_tkiso03()[j]>3.0) continue;
            //if (mu_glo_hasgsftrack()[j]) continue;
            if (LeadMuon_p4.Pt()<10.0) continue;
            if (fabs(LeadMuon_p4.Eta())>2.4) continue;
            if (debug) cout << "Lead Passing Muon: " << j << " Pt: " << LeadMuon_p4.Pt() << endl;
            for (unsigned int k=0; k<(unsigned int) mu_glo_n(); k++) {
              TLorentzVector SubLeadMuon_p4 = *((TLorentzVector*) mu_glo_p4()->At(k));
              if (debug) cout << "SubLead Trial Muon: " << k << " Pt: " << SubLeadMuon_p4.Pt() << endl;
              if (j<=k) continue;
              if (mu_glo_type()[k]<1100) continue;
              if (mu_glo_chi2()[k]/mu_glo_dof()[k]>10.0) continue;
              if (mu_glo_validhits()[k]<11) continue;
              if (mu_glo_pixelhits()[k]<1) continue;
              if (mu_glo_tkiso03()[k]>3.0) continue;
              //if (mu_glo_hasgsftrack()[k]) continue;
              if (SubLeadMuon_p4.Pt()<10.0) continue;
              if (fabs(SubLeadMuon_p4.Eta())>2.4) continue;
              if (mu_glo_charge()[j]+mu_glo_charge()[k]!=0) continue;
              TLorentzVector MuMuSystem = MuonP4[j]+MuonP4[k];
              if (MuMuSystem.M()<40 || MuMuSystem.M()>80) continue;
              if (debug) cout << "SubLead Passing Muon: " << k << " Pt: " << SubLeadMuon_p4.Pt() << endl;
              for (unsigned int l=0; l<(unsigned int) pho_n(); l++) {
                TLorentzVector TestCandidate = LeadMuon_p4+SubLeadMuon_p4+Photonp4[l];
                if (debug) cout << "Trial Photon: " << l << " Pt: " << Photonp4[l].Pt() << " Mass: " << TestCandidate.M() << endl;
                if (LeadMuon_p4.Pt() > SubLeadMuon_p4.Pt()) {
                  vector<unsigned int> indexvector;
                  indexvector.push_back(j);
                  indexvector.push_back(k);
                  indexvector.push_back(l);
                  TripletMassIndex[TestCandidate.M()] = indexvector;
                }
                else if (LeadMuon_p4.Pt() < SubLeadMuon_p4.Pt()) {
                  vector<unsigned int> indexvector;
                  indexvector.push_back(k);
                  indexvector.push_back(j);
                  indexvector.push_back(l);
                  TripletMassIndex[TestCandidate.M()] = indexvector;
                }
              }
            }
          }
          if (debug) cout << "Triplete Map Filled: " << TripletMassIndex.size() << endl;
          if (TripletMassIndex.size()==0) continue;
          double ZMassKey = 0;
          double ZMassDelta = 9999;
          for (map<double,vector<unsigned int> >::iterator MassIterator=TripletMassIndex.begin(); MassIterator!=TripletMassIndex.end(); ++MassIterator) {
            double TrialZMassDelta = fabs(MassIterator->first-RealZMass);
            if (TrialZMassDelta<ZMassDelta) {
              ZMassDelta=TrialZMassDelta;
              ZMassKey=MassIterator->first;
            }
          }
          if (debug) cout << "Triplete Selected: " << ZMassKey << endl;
          LeadMuonIndex = TripletMassIndex[ZMassKey][0];
          SubLeadMuonIndex = TripletMassIndex[ZMassKey][1];
          ZMuMuPhotonIndex = TripletMassIndex[ZMassKey][2];
          //unsigned int SCIndex = pho_scind()[PhotonIndex];
          if (debug) cout << "  LeadMuonIndex " <<  LeadMuonIndex << "  SubLeadMuonIndex " << SubLeadMuonIndex << endl;
          if (debug) cout << " Lead  pt " << MuonP4[LeadMuonIndex].Pt() << " SubLead Pt: " << MuonP4[SubLeadMuonIndex].Pt() << endl;
          if (debug) cout << "Filling MuonCharge" << endl;
          int NetCharge = mu_glo_charge()[LeadMuonIndex]+mu_glo_charge()[SubLeadMuonIndex];
          histoContainer->Fill("NetCharge",NetCharge);
          TLorentzVector MuMuSystem = MuonP4[LeadMuonIndex]+MuonP4[SubLeadMuonIndex];
          histoContainer->Fill("MuMuMass",MuMuSystem.M());
        }else { // higgs selection
          if (PhotonPtMap.size()<2) continue;
          if (debug) for (map<double,unsigned int>::reverse_iterator itPhotonPtMap=PhotonPtMap.rbegin(); itPhotonPtMap!=PhotonPtMap.rend(); itPhotonPtMap++) cout << "Photon Pt:" << itPhotonPtMap->first <<  "Photon Index: " << itPhotonPtMap->second << endl;
          PhotonIterator=PhotonPtMap.rbegin();
          LeadPhotonIndex=PhotonIterator->second;
          ++PhotonIterator;
          SubLeadPhotonIndex=PhotonIterator->second;
        }

        if (debug)  cout << " Lead " << LeadPhotonIndex << " Sub " << SubLeadPhotonIndex << "Higgs " << higgs<< endl;        
        for (unsigned int PhotonIndex=0; PhotonIndex<(unsigned int) pho_n(); PhotonIndex++) {

          if (!higgs && PhotonIndex!=ZMuMuPhotonIndex) continue;
          if (debug)  cout << " pho index " << PhotonIndex << endl;
          if ( higgs && debug ) {
            if ( PhotonIndex == LeadPhotonIndex ) cout << " Here we are with Lead " << endl; 
            if ( PhotonIndex == SubLeadPhotonIndex ) cout << " Here we are with Sub " << endl; 

          }
          if (higgs && PhotonIndex!=LeadPhotonIndex && PhotonIndex!=SubLeadPhotonIndex) continue;
          
          if (debug)  cout << " I am  passing from here " << endl;
          if (pho_isEBEEGap()[PhotonIndex]) continue;
          if ( highpt && Photonp4[PhotonIndex].Pt() < 20 ) continue;
          if ( presel && !PhotonPreSelection(PhotonIndex)) continue;
          //if (!nocuts && pho_cic4cutlevel_lead()->at(PhotonIndex).at(0)<=4) continue;


          //         if (pixel && pho_isEE()[PhotonIndex]) continue;
          //         if (pixel && (pho_isEE()[PhotonIndex] || pho_r9()[PhotonIndex]>1.0 || pho_hoe()[PhotonIndex]>0.05 || pho_sieie()[PhotonIndex]>0.011 ||
          //                       pho_ecalsumetconedr03()[PhotonIndex]+pho_hcalsumetconedr03()[PhotonIndex]+pho_trksumpthollowconedr03()[PhotonIndex]>6.0+.1211*rho())) continue;
          //if (pho_r9()[PhotonIndex]<0.3) continue;
          //if (pho_hoe()[PhotonIndex]>0.1) continue;
          //if (pho_isEB()[PhotonIndex] && pho_sieie()[PhotonIndex]>0.01) continue;
          //if (pho_isEE()[PhotonIndex] && pho_sieie()[PhotonIndex]>0.03) continue;
          //if (pho_trksumpthollowconedr03()[PhotonIndex]>2.0) continue;
          //if (pho_ecalsumetconedr03()[PhotonIndex]>4.0) continue;
          //if (pho_hcalsumetconedr03()[PhotonIndex]>5.0) continue;
          
          string PassValue = "Fail";
          if (!pixel && pho_isconv()[PhotonIndex]==1) PassValue = "Pass";
          if (pixel && pho_haspixseed()[PhotonIndex]==0) PassValue = "Pass";
        
          string Category = GetPhotonCat(PhotonIndex);

          TLorentzVector ZCandidate;
          if (!higgs) {
            double leadmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[LeadMuonIndex].Phi()));
            double leadmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[LeadMuonIndex].Eta());
            double leadmu_deltaR=leadmu_deltaeta*leadmu_deltaeta+leadmu_deltaphi*leadmu_deltaphi;
            histoContainer->Fill("LeadMuDeltaPhi",leadmu_deltaphi,weight);
            histoContainer->Fill("LeadMuDeltaEta",leadmu_deltaeta,weight);
            //if (pho_isEB()[PhotonIndex] && (leadmu_deltaeta<0.04 || leadmu_deltaphi<0.3) ) continue;
            //if (pho_isEE()[PhotonIndex] && (leadmu_deltaeta<0.08 || leadmu_deltaphi<0.3) ) continue;
            // if (leadmu_deltaR<0.3) continue;
          
            double subleadmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[SubLeadMuonIndex].Phi()));
            double subleadmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[SubLeadMuonIndex].Eta());
            double subleadmu_deltaR=subleadmu_deltaeta*subleadmu_deltaeta+subleadmu_deltaphi*subleadmu_deltaphi;
            histoContainer->Fill("SubLeadMuDeltaPhi",subleadmu_deltaphi,weight);
            histoContainer->Fill("SubLeadMuDeltaEta",subleadmu_deltaeta,weight);
            //if (pho_isEB()[PhotonIndex] && (subleadmu_deltaeta<0.04 || subleadmu_deltaphi<0.3) ) continue;
            //if (pho_isEE()[PhotonIndex] && (subleadmu_deltaeta<0.08 || subleadmu_deltaphi<0.3) ) continue;
            //if (subleadmu_deltaR<0.3) continue;
            if (leadmu_deltaR<subleadmu_deltaR) {
              NearMuonIndex=LeadMuonIndex;
              FarMuonIndex= SubLeadMuonIndex;
            } else {
              NearMuonIndex=SubLeadMuonIndex;
              FarMuonIndex = LeadMuonIndex;
            }
            if ( mu_glo_hasgsftrack()[NearMuonIndex]) continue;
            if ( MuonP4[FarMuonIndex].Pt()<30.0) continue;
       
            double nearmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[NearMuonIndex].Phi()));
            double nearmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[NearMuonIndex].Eta());
            double nearmu_deltaR=nearmu_deltaeta*nearmu_deltaeta+nearmu_deltaphi*nearmu_deltaphi;
            if ( nearMuonDRcut && nearmu_deltaR<0.2) continue; 

            //          if (leadmu_deltaR<subleadmu_deltaR && mu_glo_hasgsftrack()[LeadMuonIndex]) continue;
            // if (leadmu_deltaR>subleadmu_deltaR && mu_glo_hasgsftrack()[SubLeadMuonIndex]) continue;
            //          if (leadmu_deltaR<subleadmu_deltaR && MuonP4[SubLeadMuonIndex].Pt()<30.0) continue;
            //if (leadmu_deltaR>subleadmu_deltaR && MuonP4[LeadMuonIndex].Pt()<30.0) continue;
 
            ZCandidate = Photonp4[PhotonIndex]+MuonP4[LeadMuonIndex]+MuonP4[SubLeadMuonIndex];
            //if (ZCandidate.M()<75.0 && ZCandidate.M()>105.0) continue;
            if (debug && (ZCandidate.M()<87.2 || ZCandidate.M()>95.2)) cout << run() << ":" << lumis() << ":" << (unsigned int) event() << ":" << Photonp4[PhotonIndex].Pt() << endl;
          
            if (debug) cout << "Lead Muon Pt: " << MuonP4[LeadMuonIndex].Pt() << " SubLead Muon Pt: " << MuonP4[SubLeadMuonIndex].Pt() << endl;
            if (debug) cout << "Photon Pt: " << Photonp4[PhotonIndex].Pt() << " ZMass: " << ZCandidate.M() << endl;
          
            histoContainer->Fill("LeadMuPt",MuonP4[LeadMuonIndex].Pt(),weight);
            histoContainer->Fill("SubLeadMuPt",MuonP4[SubLeadMuonIndex].Pt(),weight);
            histoContainer->Fill("NearMuPt",MuonP4[NearMuonIndex].Pt(),weight);
            histoContainer->Fill("LeadMuDeltaR",leadmu_deltaR,weight);
            histoContainer->Fill("SubLeadMuDeltaR",subleadmu_deltaR,weight);
            histoContainer->Fill("NearMuDeltaR",nearmu_deltaR,weight);

          }
          

          string region=DetectorPosition(PhotonIndex);
          if (debug) cout << "Filling MVA Quantities" << endl;
          histoContainer->Fill("pho_pt",region,Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill("pho_idmva",region,pho_idmva()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_tiso1",region,pho_tmva_id_mit_tiso1()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_tiso2",region,pho_tmva_id_mit_tiso2()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_tiso3",region,pho_tmva_id_mit_tiso3()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_r9",region,pho_tmva_id_mit_r9()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_ecal",region,pho_tmva_id_mit_ecal()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_hcal",region,pho_tmva_id_mit_hcal()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_etawidth",region,pho_tmva_id_mit_etawidth()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_phiwidth",region,pho_tmva_id_mit_phiwidth()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_nvtx",region,pho_tmva_id_mit_nvtx()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_preshower",region,pho_tmva_id_mit_preshower()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_sceta",region,pho_tmva_id_mit_sceta()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_hoe",region,pho_tmva_id_mit_hoe()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_id_mit_sieie",region,pho_tmva_id_mit_sieie()[PhotonIndex],weight);




          histoContainer->Fill("Numvtx",PrimaryVertex.size(),weight);
          histoContainer->Fill("PhotonPt",Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill("PhotonEt",Photonp4[PhotonIndex].Et(),weight);
          histoContainer->Fill("PhotonEta",Photonp4[PhotonIndex].Eta(),weight);
          histoContainer->Fill("PhotonPhi",Photonp4[PhotonIndex].Phi(),weight);
          histoContainer->Fill("ZMass",ZCandidate.M(),weight);
          histoContainer->Fill("ZMassZoom",ZCandidate.M(),weight);

          histoContainer->Fill("Numvtx",Category,PrimaryVertex.size(),weight);
          histoContainer->Fill("PhotonPt",Category,Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill("PhotonEt",Category,Photonp4[PhotonIndex].Et(),weight);
          histoContainer->Fill("PhotonEta",Category,Photonp4[PhotonIndex].Eta(),weight);
          histoContainer->Fill("PhotonPhi",Category,Photonp4[PhotonIndex].Phi(),weight);
          histoContainer->Fill("ZMass",Category,ZCandidate.M(),weight);
          histoContainer->Fill("ZMassZoom",Category,ZCandidate.M(),weight);

          histoContainer->Fill(PassValue,"Numvtx",PrimaryVertex.size(),weight);
          histoContainer->Fill(PassValue,"PhotonPt",Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill(PassValue,"PhotonEt",Photonp4[PhotonIndex].Et(),weight);
          histoContainer->Fill(PassValue,"PhotonEta",Photonp4[PhotonIndex].Eta(),weight);
          histoContainer->Fill(PassValue,"PhotonPhi",Photonp4[PhotonIndex].Phi(),weight);
          histoContainer->Fill(PassValue,"ZMass",ZCandidate.M(),weight);
          histoContainer->Fill(PassValue,"ZMassZoom",ZCandidate.M(),weight);

          histoContainer->Fill(PassValue,"Numvtx",Category,PrimaryVertex.size(),weight);
          histoContainer->Fill(PassValue,"PhotonPt",Category,Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill(PassValue,"PhotonEt",Category,Photonp4[PhotonIndex].Et(),weight);
          histoContainer->Fill(PassValue,"PhotonEta",Category,Photonp4[PhotonIndex].Eta(),weight);
          histoContainer->Fill(PassValue,"PhotonPhi",Category,Photonp4[PhotonIndex].Phi(),weight);
          histoContainer->Fill(PassValue,"ZMass",Category,ZCandidate.M(),weight);
          histoContainer->Fill(PassValue,"ZMassZoom",Category,ZCandidate.M(),weight);
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

bool PhotonPreSelection(unsigned int PhotonIndex) {

  TLorentzVector PhotonP4 = *((TLorentzVector*) pho_p4()->At(PhotonIndex));
  double EtCorrEcalIso = pho_ecalsumetconedr03()[PhotonIndex] - 0.012*PhotonP4.Pt();
  double EtCorrHcalIso = pho_hcalsumetconedr03()[PhotonIndex] - 0.005*PhotonP4.Pt();
  double EtCorrTrkIso = pho_trksumpthollowconedr03()[PhotonIndex] - 0.002*PhotonP4.Pt();
  double PuCorrHcalEcal = pho_ecalsumetconedr03()[PhotonIndex] + pho_hcalsumetconedr03()[PhotonIndex] - rho()*0.17;
  double AbsTrkIsoCIC = pho_tkiso_recvtx_030_002_0000_10_01()->at(PhotonIndex).at(0);
  
  if (pho_r9()[PhotonIndex]<=0.9) {
    if (pho_isEB()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.014)) return false;
    if (pho_isEE()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.034)) return false;
    if (EtCorrEcalIso>4.0) return false;
    if (EtCorrHcalIso>4.0) return false;
    if (EtCorrTrkIso>4.0) return false;
    if (PuCorrHcalEcal>3.0) return false;
    if (AbsTrkIsoCIC>2.8) return false;
    if (pho_trksumpthollowconedr03()[PhotonIndex]>4.0) return false;
    if (pho_isconv()[PhotonIndex]!=1) return false;
    return true;
  } else {
    if (pho_isEB()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.082 || pho_sieie()[PhotonIndex]>0.014)) return false;
    if (pho_isEE()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.034)) return false;
    if (EtCorrEcalIso>50.0) return false;
    if (EtCorrHcalIso>50.0) return false;
    if (EtCorrTrkIso>50.0) return false;
    if (PuCorrHcalEcal>3.0) return false;
    if (AbsTrkIsoCIC>2.8) return false;
    if (pho_trksumpthollowconedr03()[PhotonIndex]>4.0) return false;
    if (pho_isconv()[PhotonIndex]!=1) return false;
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

  kFactor["ZMuMu"]=1.0;
  kFactor["TTJets"]=1.0;
  kFactor["Higgs120GeV"]=1.0;
  return kFactor;
  
}

map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight) {

  map<TString,double> WeightsMap;
  WeightsMap["None"]=1/globalweight;
  WeightsMap["ZMuMu_Summer11"]=kFactor["ZMuMu"]*1626.0/29743564.0;
  WeightsMap["ZMuMu_Fall11"]=kFactor["ZMuMu"]*1626.0/29743564.0;
  WeightsMap["TTJets"]=kFactor["TTJets"]*94.76/3701947.0;
  WeightsMap["HiggsS4"]=kFactor["Higgs120GeV"]*16.63/105132;
  WeightsMap["HiggsS6"]=kFactor["Higgs120GeV"]*16.63/105132;
  return WeightsMap;
  
}

template <class type> string makestring(type value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
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
    if (conv_ntracks()[i]!=2 && conv_ntracks()[i]!=1) continue;
    if (conv_ntracks()[i]==2 && (conv_chi2_probability()[i]<0.000001 || ConversionRefittedPairMomentum.Pt()<1)) continue;
    if (conv_ntracks()[i]==1 && ConversionRefittedPairMomentum.Pt()<1) continue;
    
    double deltaphi = DeltaPhi(Photonxyz.Phi(),ConversionVertex.Phi());
    double zfromconv = FindNewZConvLinear(ConversionVertex,Photonxyz,BeamSpot);
    //cout << "NTracks: " << conv_ntracks()[i] << " Conerion Refitted Pair Momentum Eta: " << ConversionRefittedPairMomentum.Eta() << endl;
    double deltaeta = fabs(Photonxyz.Eta() - etaTransformation(ConversionRefittedPairMomentum.Eta(),zfromconv));
    double deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    
    if (deltaR<MindeltaR) {
      MindeltaR=deltaR;
      ReturnIndex = i;
    }
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

void BookCutsAndCategories(HistoContainer *histoContainer, TString histname, TString histtitle, int bins, float lowerlimit, float upperlimit) {

  TString Cuts[3] = {"","Pass","Fail"};
  TString CutLabel[3] = {" "," Passing Electron Veto "," Failing Electron Veto "};
  TString Categories[5] = {"","_cat0","_cat1","_cat2","_cat3"};
  TString CatLabel[5] = {"","R9>0.94 Barrel","R9<0.94 Barrel","R9>0.94 Endcap","R9<0.94 Endcap"};
  
  for (int cut=0; cut<3; cut++) {
    for (int cat=0; cat<5; cat++) {
      TString histnametemp = Cuts[cut];
      histnametemp += histname;
      histnametemp += Categories[cat];
      TString histtitletemp = histtitle;
      histtitletemp.ReplaceAll(":Cuts",CutLabel[cut]);
      histtitletemp.ReplaceAll(":Cat",CatLabel[cat]);
      histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),bins,lowerlimit,upperlimit);
    }
  }

}

void BookHistograms(HistoContainer *histoContainer) {

  histoContainer->Add("NetCharge","Net Charged;Charge;Counts",5,-2.5,2.5);
  histoContainer->Add("MuMuMass","Mass of diMuon System;Mass (GeV);Counts",100,0,100);
  
  histoContainer->Add("LeadMuDeltaEta","#Delta#eta to lead #mu;#Delta#eta;Counts",25,0,5);
  histoContainer->Add("SubLeadMuDeltaEta","#Delta#eta to sublead #mu;#Delta#eta;Counts",25,0,5);
  histoContainer->Add("LeadMuDeltaPhi","#Delta#phi to lead #mu;#Delta#phi;Counts",34,0,3.4);
  histoContainer->Add("SubLeadMuDeltaPhi","#Delta#phi to sublead #mu;#Delta#phi;Counts",34,0,3.4);
  histoContainer->Add("LeadMuDeltaR","#DeltaR to lead #mu;#DeltaR;Counts",30,0,6);
  histoContainer->Add("SubLeadMuDeltaR","#DeltaR to sublead #mu;#DeltaR;Counts",30,0,6);
  histoContainer->Add("LeadMuPt","Pt of Lead Muon;Pt (GeV);Counts",100,0,100);
  histoContainer->Add("SubLeadMuPt","Pt of SubLead Muon;Pt (GeV);Counts",100,0,100);
  histoContainer->Add("NearMuPt","Pt of Near Muon;Pt (GeV);Counts",100,0,100);
  histoContainer->Add("NearMuDeltaR","Near Muon #DeltaR  #mu;#DeltaR;Counts",30,0,6);

  BookCutsAndCategories(histoContainer,"Numvtx","Number of Primary Verticies:Cuts:Cat;Number of Vertices;Counts",100,0,100);
  BookCutsAndCategories(histoContainer,"PhotonPt","Pt of Photon:Cuts:Cat;Pt (GeV);Counts",100,0,200);
  BookCutsAndCategories(histoContainer,"PhotonEt","Et of Photon:Cuts:Cat;Et (GeV);Counts",50,0,100);
  BookCutsAndCategories(histoContainer,"PhotonEta","#eta of Photon:Cuts:Cat;#eta;Counts",32,-3.2,3.2);
  BookCutsAndCategories(histoContainer,"PhotonPhi","#phi of Photon:Cuts:Cat;#phi;Counts",32,-3.2,3.2);
  BookCutsAndCategories(histoContainer,"ZMass","Invariant mass of #mu#mu#gamma system:Cuts:Cat;Mass (GeV);Counts",30,75,105);
  BookCutsAndCategories(histoContainer,"ZMassZoom","Invariant mass of #mu#mu#gamma system:Cuts:Cat;Mass (GeV);Counts",16,87.2,95.2);

  BookBarrelAndEndcap(histoContainer,"pho_pt","pho_pt: region; Pt (GeV) ;Counts",100,0.,100.);
  BookBarrelAndEndcap(histoContainer,"pho_idmva","pho_idmva: region;idmva value;Counts",150,-1,0.5);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_tiso1","pho_tmva_id_mit_tiso1: region;tmva_id_mit_tiso1;Counts",110,-10,100);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_tiso2","pho_tmva_id_mit_tiso2: region;tmva_id_mit_tiso2;Counts",110,-10,100);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_tiso3","pho_tmva_id_mit_tiso3: region;tmva_id_mit_tiso3;Counts",110,-10,100);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_r9","pho_tmva_id_mit_r9: region;tmva_id_mit_r9;Counts",150,0,1.5);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_ecal","pho_tmva_id_mit_ecal: region;tmva_id_mit_ecal;Counts",50,-10,40);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_hcal","pho_tmva_id_mit_hcal: region;tmva_id_mit_hcal;Counts",40,-10,30);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_etawidth","pho_tmva_id_mit_etawidth: region;tmva_id_mit_etawidth;Counts",100,0,0.1);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_phiwidth","pho_tmva_id_mit_phiwidth: region;tmva_id_mit_phiwidth;Counts",100,0,0.2);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_nvtx","pho_tmva_id_mit_nvtx: region;tmva_id_mit_nvtx;Counts",50,0,50);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_preshower","pho_tmva_id_mit_preshower: region;tmva_id_mit_preshower;Counts",50,0,0.5);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_sceta","pho_tmva_id_mit_sceta: region;tmva_id_mit_sceta;Counts",68,-3.4,3.4);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_hoe","pho_tmva_id_mit_hoe: region;tmva_id_mit_hoe;Counts",60,0,0.6);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_id_mit_sieie","pho_tmva_id_mit_sieie: region;tmva_id_mit_sieie;Counts",100,0,0.1);

}

void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  if (inputstring.Contains("Run2011A")) {
    inputfilelist.push_back(pair<string,int> ("Run2011A.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/Run2011A.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2011B")) {
    inputfilelist.push_back(pair<string,int> ("Run2011B.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/Run2011B.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Data")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/Run2011A.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/Run2011B.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Fall11")) {
    inputfilelist.push_back(pair<string,int> ("Fall11.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/DYToMuMu_M-20_CT10_TuneZ2_7TeV_Fall11.root",WeightsMap["ZMuMu_Fall11"]));
  }
  if (inputstring.Contains("Summer11")) {
    inputfilelist.push_back(pair<string,int> ("Summer11.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/DYToMuMu_M-20_CT10_TuneZ2_7TeV_Summer11.root",WeightsMap["ZMuMu_Summer11"]));
  }
  if (inputstring.Contains("TTJets")) {
    inputfilelist.push_back(pair<string,int> ("TTJets.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/TTJets.root",WeightsMap["TTJets"]));
  }
  if (inputstring.Contains("HiggsS4")) {
    inputfilelist.push_back(pair<string,int> ("HiggsS4.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/120GeVHiggs_S4.root",WeightsMap["HiggsS4"]));
  }
  if (inputstring.Contains("HiggsS6")) {
    inputfilelist.push_back(pair<string,int> ("HiggsS6.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/ZMuMuGamma/120GeVHiggs_S6.root",WeightsMap["HiggsS6"]));
  }

}

void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  string outfile;
  if (infile.rfind("/")!=string::npos) outfile=infile.substr(infile.rfind("/")+1);
  else outfile=infile;
  
  inputfilelist.push_back(pair<string,int> (outfile,1));
  if (inputstring.Contains("Run2011A") || inputstring.Contains("Run2011B") || inputstring.Contains("Data")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("TTJets")) inputvector.push_back(pair<string,float> (infile,WeightsMap["TTJets"]));
  if (inputstring.Contains("Fall11")) inputvector.push_back(pair<string,float> (infile,WeightsMap["ZMuMu_Fall11"]));
  if (inputstring.Contains("Summer11")) inputvector.push_back(pair<string,float> (infile,WeightsMap["ZMuMu_Summer11"]));
  if (inputstring.Contains("HiggsS4")) inputvector.push_back(pair<string,float> (infile,kFactor["Higgs120GeV"]*WeightsMap["HiggsS4"]));
  if (inputstring.Contains("HiigsS6")) inputvector.push_back(pair<string,float> (infile,kFactor["Higgs120GeV"]*WeightsMap["HiggsS6"]));

}

void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap) {

  if (inputstring.Contains("TTJets")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_TTJets.h"
  } else if (inputstring.Contains("Summer11")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_Summer11.h"
  } else if (inputstring.Contains("Fall11")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_Fall11.h"
  } else if (inputstring.Contains("HiggsS4")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_HiggsS4.h"
  } else if (inputstring.Contains("HiggsS6")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_HiggsS6.h"
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
