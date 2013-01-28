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
#include "TTree.h"

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

//bool GenMatch(TVector3 Photon);
bool PhotonPreSelection(unsigned int PhotonIndex);
bool PhotonPreSelectionPFbased(unsigned int PhotonIndex);
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
void BookBarrelAndEndcap2D(HistoContainer *histoContainer, TString histname, TString histtitle, int binsa, float lowerlimita, float upperlimita,int binsb, float lowerlimitb, float upperlimitb);
void BookBarrelAndEndcapProfiles(HistoContainer *histoContainer, TString histname, TString histtitle, int binsa, float lowerlimita, float upperlimita, float lowerlimitb, float upperlimitb);
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
  bool varCorrections2011  = false;
  bool varCorrections2012  = false;
  bool preselpf = false;
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
  float globalweight = 19092.8;

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
  if (InputArgs.Contains("PFiso")) preselpf=true;
  if (InputArgs.Contains("VarCorr2012")) varCorrections2012=true;
  if (InputArgs.Contains("VarCorr2011")) varCorrections2011=true;
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
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: Run2011A, Run2011B, Run2011, Run2012AB, Run2012All, ZMuMu, or TTJets." << endl;
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
      if (preselpf) outfilename.ReplaceAll(".root","_PhoPFPresel.root");
      if (nearMuonDRcut ) outfilename.ReplaceAll(".root","_NearMuonDRCut.root");
      if (varCorrections2012)  outfilename.ReplaceAll(".root","_Corr2012.root");
      if (varCorrections2011)  outfilename.ReplaceAll(".root","_Corr2011.root");
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

    unsigned int mypho_N=0;
    unsigned int mypho_det=0;
    float mypho_pt=-99.;
    float mypho_r9=-99.;
    float mypho_sieie=-99.;
    float mypho_sieip=-99.;
    float mypho_etawidth=-99.;
    float mypho_phiwidth=-99.;
    float mypho_s4ratio=-99.;
    float mypho_lambdaratio=-99.;
    float mypho_sceta=-99.;
    float mypho_eventrho=-99.;
    float mypho_ESEffSigmaRR=-99.;
    TTree* minitree;
    minitree = new TTree("ZMuMuGamMinitree","An example of ROOT tree with a few branches");
    minitree->Branch("pho_N", &mypho_N, "pho_N/I");
    minitree->Branch("pho_det",&mypho_det,"pho_det/I");
    minitree->Branch("pho_pt",&mypho_pt,"pho_pt/F");
    minitree->Branch("pho_r9",&mypho_r9,"pho_r9/F");
    minitree->Branch("pho_sieie",&mypho_sieie,"pho_sieie/F");
    minitree->Branch("pho_sieip",&mypho_sieip,"pho_sieip/F");
    minitree->Branch("pho_etawidth",&mypho_etawidth,"pho_etawidth/F");
    minitree->Branch("pho_phiwidth",&mypho_phiwidth,"pho_phiwidth/F");
    minitree->Branch("pho_s4ratio",&mypho_s4ratio,"pho_s4ratio/F");
    minitree->Branch("pho_lambdaratio",&mypho_lambdaratio,"pho_lambdaratio/F");
    minitree->Branch("pho_sceta",&mypho_sceta,"pho_sceta/F");
    minitree->Branch("pho_eventrho",&mypho_eventrho,"pho_eventrho/F");
    minitree->Branch("pho_ESEffSigmaRR",&mypho_ESEffSigmaRR,"pho_ESEffSigmaRR/F");

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
        if (ientry==0 && gv_n()==0) mc=false;

        if (ientry==0 && trigger) {
          cout << "Trigger Names:" << endl;
          for (unsigned int j=0; j<(unsigned int) hlt_path_names_HLT()->size(); j++) {
            TString TriggerName(hlt_path_names_HLT()->at(j));
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
          //if (debug) cout << "Has GSFkTrack: " << mu_glo_hasgsftrack()[j] << endl;
          //if (mu_glo_hasgsftrack()[j]) continue;
          if (Muon_p4.Pt()<10.0) continue;
          if (fabs(Muon_p4.Eta())>2.4) continue;
          MuonPtMap[Muon_p4.Pt()]=j;
        }
        map<double,unsigned int>::reverse_iterator PhotonIterator;
        unsigned int LeadMuonIndex=0;
        unsigned int SubLeadMuonIndex=0;
        unsigned int ZMuMuPhotonIndex=0;
        unsigned int LeadPhotonIndex=0;
        unsigned int SubLeadPhotonIndex = 0;
        unsigned int NearMuonIndex=0;
        unsigned int FarMuonIndex =0 ;
        if (!higgs) {
          if (MuonPtMap.size()<2) continue;
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
          if (debug) cout << "LeadMuonIndex " <<  LeadMuonIndex << " SubLeadMuonIndex " << SubLeadMuonIndex << endl;
          if (debug) cout << "Lead pt " << MuonP4[LeadMuonIndex].Pt() << " SubLead Pt: " << MuonP4[SubLeadMuonIndex].Pt() << endl;
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

        unsigned int nPho=0;
        for (unsigned int PhotonIndex=0; PhotonIndex<(unsigned int) pho_n(); PhotonIndex++) {
          unsigned int dete=0;
          if (!higgs && PhotonIndex!=ZMuMuPhotonIndex) continue;
          if (debug)  cout << " pho index " << PhotonIndex << endl;
          if ( higgs && debug ) {
            if ( PhotonIndex == LeadPhotonIndex ) cout << " Here we are with Lead " << endl;
            if ( PhotonIndex == SubLeadPhotonIndex ) cout << " Here we are with Sub " << endl;

          }
          if (higgs && PhotonIndex!=LeadPhotonIndex && PhotonIndex!=SubLeadPhotonIndex) continue;
          string region=DetectorPosition(PhotonIndex);

          if (debug)  cout << " I am  passing from here " << endl;
          if (pho_isEBEEGap()[PhotonIndex]) continue;

          histoContainer->Fill("allpho_pt",region,Photonp4[PhotonIndex].Pt(),weight);

          if (region == "Endcap") dete=1;

          double leadmu_deltaphi=-1;
          double leadmu_deltaeta=-1;
          double leadmu_deltaR=-1;
          double subleadmu_deltaphi=-1;
          double subleadmu_deltaeta=-1;
          double subleadmu_deltaR=-1;
          double nearmu_deltaphi=-1;
          double nearmu_deltaeta =-1;
          double nearmu_deltaR=-1;
          if (!higgs) { // define muon quantities if not an higgs sample
            leadmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[LeadMuonIndex].Phi()));
            leadmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[LeadMuonIndex].Eta());
            leadmu_deltaR=sqrt(leadmu_deltaeta*leadmu_deltaeta+leadmu_deltaphi*leadmu_deltaphi);

            subleadmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[SubLeadMuonIndex].Phi()));
            subleadmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[SubLeadMuonIndex].Eta());
            subleadmu_deltaR=sqrt(subleadmu_deltaeta*subleadmu_deltaeta+subleadmu_deltaphi*subleadmu_deltaphi);

            if (leadmu_deltaR<subleadmu_deltaR) {
              NearMuonIndex=LeadMuonIndex;
              FarMuonIndex= SubLeadMuonIndex;
            } else {
              NearMuonIndex=SubLeadMuonIndex;
              FarMuonIndex = LeadMuonIndex;
            }

            nearmu_deltaphi=fabs(DeltaPhi(Photonp4[PhotonIndex].Phi(),MuonP4[NearMuonIndex].Phi()));
            nearmu_deltaeta=fabs(Photonp4[PhotonIndex].Eta()-MuonP4[NearMuonIndex].Eta());
            nearmu_deltaR=sqrt(nearmu_deltaeta*nearmu_deltaeta+nearmu_deltaphi*nearmu_deltaphi);

          } // end of !higgs

          float pho_s4ratio = pho_e2x2()[PhotonIndex]/bc_s25()[sc_bcseedind()[pho_scind()[PhotonIndex]]];
          float rr2=pho_eseffsixix()[PhotonIndex]*pho_eseffsixix()[PhotonIndex]+pho_eseffsiyiy()[PhotonIndex]*pho_eseffsiyiy()[PhotonIndex];
          float pho_ESEffSigmaRR = 0.0;
          if(rr2>0. && rr2<999999.) {
            pho_ESEffSigmaRR = sqrt(rr2);
          }
          
          float pho_tmva_photonid_pfchargedisogood03 = pho_pfiso_mycharged03()->at(PhotonIndex).at(0);
          float pho_tmva_photonid_pfchargedisobad03  = 0.;
          for(int ivtx=0; ivtx<vtx_std_n(); ivtx++) {
            pho_tmva_photonid_pfchargedisobad03 = pho_pfiso_mycharged03()->at(PhotonIndex).at(ivtx)> pho_tmva_photonid_pfchargedisobad03 ? pho_pfiso_mycharged03()->at(PhotonIndex).at(ivtx): pho_tmva_photonid_pfchargedisobad03 ;
          }
          float  pho_tmva_photonid_pfphotoniso03  = pho_pfiso_myphoton03()[PhotonIndex];
          float  pho_tmva_photonid_pfneutraliso03 = pho_pfiso_myneutral03()[PhotonIndex];
          histoContainer->Fill("pho_pfchargedisogood03_vs_nearMuonDRBeforePhoPresel",region,nearmu_deltaR,pho_tmva_photonid_pfchargedisogood03,weight);
          histoContainer->Fill("pho_pfchargedisobad03_vs_nearMuonDRBeforePhoPresel",region,nearmu_deltaR,pho_tmva_photonid_pfchargedisobad03,weight);
          histoContainer->Fill("pho_pfphotoniso03_vs_nearMuonDRBeforePhoPresel",region,nearmu_deltaR,pho_tmva_photonid_pfphotoniso03,weight);
          histoContainer->Fill("pho_pfneutraliso03_vs_nearMuonDRBeforePhoPresel",region,nearmu_deltaR,pho_tmva_photonid_pfneutraliso03,weight);

          string Category = GetPhotonCat(PhotonIndex);

          TLorentzVector ZCandidate;
          if (!higgs) {
            histoContainer->Fill("LeadMuDeltaPhi",leadmu_deltaphi,weight);
            histoContainer->Fill("LeadMuDeltaEta",leadmu_deltaeta,weight);
            histoContainer->Fill("SubLeadMuDeltaPhi",subleadmu_deltaphi,weight);
            histoContainer->Fill("SubLeadMuDeltaEta",subleadmu_deltaeta,weight);

            if ( mu_glo_hasgsftrack()[NearMuonIndex]) continue;
            if ( MuonP4[FarMuonIndex].Pt()<30.0) continue;

            //histoContainer->Fill("pho_hcalIso_vs_nearMuonDRBeforeCut",region,nearmu_deltaR,pho_tmva_id_mit_hcal()[PhotonIndex],weight);
            //histoContainer->Fill("pho_hoe_vs_nearMuonDRBeforeCut",region,nearmu_deltaR,pho_tmva_id_mit_hoe()[PhotonIndex],weight);
            //histoContainer->Fill("pho_tiso2_vs_nearMuonDRBeforeCut",region,nearmu_deltaR,pho_tmva_id_mit_tiso2()[PhotonIndex],weight);
            histoContainer->Fill("NearMuDeltaRBeforeCut",nearmu_deltaR,weight);
            histoContainer->Fill("NearMuDeltaEtaDeltaPhiBeforeCut",nearmu_deltaeta,nearmu_deltaphi,weight);

            if ( nearMuonDRcut && nearmu_deltaR<0.2) continue;
            //if ( nearMuonDRcut && ( nearmu_deltaphi < 0.1 &&   nearmu_deltaeta < 0.3 ) ) continue;

            histoContainer->Fill("NearMuDeltaRAfterCut",nearmu_deltaR,weight);
            histoContainer->Fill("NearMuDeltaEtaDeltaPhiAfterCut",nearmu_deltaeta,nearmu_deltaphi,weight);

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

          }

          if ( highpt && Photonp4[PhotonIndex].Pt() < 20 ) continue;
          if ( presel && !PhotonPreSelection(PhotonIndex)) continue;
          if ( preselpf && !PhotonPreSelectionPFbased(PhotonIndex)) continue;

          float r9W=1.;
          float r9Slope=0.;
          float etawidthW=1.;
          float etawidthSlope=0.;
          float phiwidthW=1.;
          float phiwidthSlope=0.;
          float sietaietaW=1.;
          float sietaietaSlope=0.;
          float s4W=1.;
          float s4Slope=0.;
          float lambdarW=1.;
          float eseffSigW=1.;
          float eseffSigSlope=0.;

          if  ( varCorrections2011 && mc )   {
            r9W*=1.0035;
            phiwidthW*=0.99;
            etawidthW*=0.98;
            if ( region=="Barrel") {
              s4W=1.0055;
              sietaietaW*=0.87;
              sietaietaSlope=0.0011;

            } else {
              s4W=1.0085;
              sietaietaW*=0.99;
              eseffSigW=1.04;
            }
          }

          if  ( varCorrections2012 && mc )   {

            if ( region=="Barrel") {

              r9W=1.0045;
              r9Slope=0.0010;
              s4W=1.01894;
              s4Slope=-0.01034;
              etawidthW=1.04302;
              etawidthSlope=-0.000618;
              phiwidthW=1.00002;
              phiwidthSlope=-0.000371;

              if ( pho_sieie()[PhotonIndex]  < 0.0087 ) {
                sietaietaW=0.976591;
                sietaietaSlope = 0.000081;
              } else {
                sietaietaW=0.980055;
                sietaietaSlope =  0.0000508632;
              }

            }else {

              r9W=1.0086;
              r9Slope=-0.0007;
              s4W=1.04969;
              s4Slope=-0.03642;
              etawidthW=0.903254;
              etawidthSlope=0.001346;
              phiwidthW=0.99992;
              phiwidthSlope=4.8e-07;
              eseffSigW=1.00023;
              eseffSigSlope=0.0913;
              sietaietaW=0.99470;
              sietaietaSlope =0.00003;

            }

            /*
              phiwidthW*=0.99;
              if ( region=="Barrel") {
              if ( InputArgs.Contains("PFiso")) {
              r9W*=1.;
              sietaietaW*=0.9869;
              etawidthW*=0.9834;
              s4W=1.0007;
              } else {
              r9W*=1.0033;
              sietaietaW*=0.9869;
              etawidthW*=1.0143;
              s4W=1.0007;
              }
              lambdarW=0.99;
              } else {
              if ( InputArgs.Contains("PFiso")) {
              r9W*=1.0048;
              sietaietaW*=1.0115;
              etawidthW*=0.9886;
              s4W=0.9945;
              } else {
              r9W*=1.0021;
              sietaietaW*=1.0112;
              etawidthW*=0.9571;
              s4W=0.9945;
              }
              }

            */

          }

          if (debug) cout << " Corrections being applied " << " r9 " << r9W << " etawidth " << etawidthW << " phywidth" << phiwidthW << " sieie " <<  sietaietaW << " sieie slope " << sietaietaSlope << " s4 " << s4W << " lambdar " << lambdarW << " eseffSig " <<   eseffSigW << endl;
          if (debug) cout << "Filling MVA Quantities" << endl;

          histoContainer->Fill("pho_pt",region,Photonp4[PhotonIndex].Pt(),weight);
          histoContainer->Fill("pho_mitmva",region,pho_mitmva()->at(PhotonIndex).at(0),weight);
          histoContainer->Fill("pho_tmva_photonid_sieip",region,pho_sieip()[PhotonIndex],weight);
          histoContainer->Fill("pho_tmva_photonid_etawidth",region,pho_etawidth()[PhotonIndex]*etawidthW+etawidthSlope,weight);
          int scind=pho_scind()[PhotonIndex];
          histoContainer->Fill("pho_tmva_photonid_phiwidth",region,sc_sphi()[scind]*phiwidthW+phiwidthSlope,weight);
          histoContainer->Fill("pho_tmva_photonid_r9",region,pho_r9()[PhotonIndex]*r9W+r9Slope,weight);
          histoContainer->Fill("pho_tmva_photonid_s4ratio",region,pho_s4ratio*s4W+s4Slope,weight);
          histoContainer->Fill("pho_tmva_photonid_lambdaratio",region,pho_lambdaratio()[PhotonIndex]*lambdarW,weight);
          histoContainer->Fill("pho_tmva_photonid_sceta",region, ((TVector3*)sc_xyz()->At(scind))->Eta(),weight);
          histoContainer->Fill("pho_tmva_photonid_sceta", ((TVector3*)sc_xyz()->At(scind))->Eta(),weight);
          histoContainer->Fill("pho_tmva_photonid_eventrho",region,rho_algo1(),weight);

          histoContainer->Fill("pho_tmva_photonid_ESEffSigmaRR",region,pho_ESEffSigmaRR*eseffSigW+eseffSigSlope,weight);
          histoContainer->Fill("pho_tmva_photonid_sieie",region,pho_sieie()[PhotonIndex]*sietaietaW+sietaietaSlope,weight);
          histoContainer->Fill("pho_tmva_photonid_pfchargedisogood03",region,pho_tmva_photonid_pfchargedisogood03,weight);
          histoContainer->Fill("pho_tmva_photonid_pfchargedisobad03",region,pho_tmva_photonid_pfchargedisobad03,weight);
          histoContainer->Fill("pho_tmva_photonid_pfphotoniso03",region,pho_tmva_photonid_pfphotoniso03,weight);
          histoContainer->Fill("pho_tmva_photonid_pfneutraliso03",region,pho_tmva_photonid_pfneutraliso03,weight);

          histoContainer->Fill("Numvtx",PrimaryVertex.size(),weight);
          histoContainer->Fill("Rho",rho_algo1(),weight);
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

          string PassValue = "Fail";
          if (!pixel && pho_isconv()[PhotonIndex]==1) PassValue = "Pass";
          if (pixel && pho_haspixseed()[PhotonIndex]==0) PassValue = "Pass";

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

          float MVAcut=-0.013;
          if ( region == "Endcap" )
            MVAcut=-0.011;

          histoContainer->Fill("pho_idmva_vsPt",region,Photonp4[PhotonIndex].Pt(),pho_mitmva()->at(PhotonIndex).at(0), weight);

          if ( pho_mitmva()->at(PhotonIndex).at(0) < MVAcut )
            histoContainer->Fill("pho_pt_afterMVAcut",region,Photonp4[PhotonIndex].Pt(),weight);

          mypho_pt=-99.;
          mypho_r9=-99.;
          mypho_sieie=-99.;
          mypho_sieip=-99.;
          mypho_etawidth=-99.;
          mypho_phiwidth=-99.;
          mypho_s4ratio=-99.;
          mypho_lambdaratio=-99.;
          mypho_sceta=-99.;
          mypho_eventrho=-99.;
          mypho_ESEffSigmaRR=-99.;

          nPho++;
          mypho_N            =  nPho;
          mypho_det          =  dete;
          mypho_pt           =  Photonp4[PhotonIndex].Pt();
          mypho_r9           =  pho_r9()[PhotonIndex]*r9W+r9Slope;
          mypho_sieie        =  pho_sieie()[PhotonIndex]*sietaietaW+sietaietaSlope;
          mypho_sieip        =  pho_sieip()[PhotonIndex];
          mypho_etawidth     =  pho_etawidth()[PhotonIndex]*etawidthW+etawidthSlope;
          mypho_phiwidth     =  sc_sphi()[scind]*phiwidthW+phiwidthSlope;
          mypho_s4ratio      =  pho_s4ratio*s4W+s4Slope;
          mypho_lambdaratio  =  pho_lambdaratio()[PhotonIndex]*lambdarW;
          mypho_sceta        =   ((TVector3*)sc_xyz()->At(scind))->Eta();
          mypho_eventrho     =  rho_algo1();
          mypho_ESEffSigmaRR =  pho_ESEffSigmaRR*eseffSigW+eseffSigSlope;

          minitree->Fill();

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

/*
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
*/

bool PhotonPreSelection(unsigned int PhotonIndex) {

  TLorentzVector PhotonP4 = *((TLorentzVector*) pho_p4()->At(PhotonIndex));
  double EtCorrEcalIso = pho_ecalsumetconedr03()[PhotonIndex] - 0.012*PhotonP4.Pt();
  double EtCorrHcalIso = pho_hcalsumetconedr03()[PhotonIndex] - 0.005*PhotonP4.Pt();
  double EtCorrTrkIso = pho_trksumpthollowconedr03()[PhotonIndex] - 0.002*PhotonP4.Pt();
  double PuCorrHcalEcal = pho_ecalsumetconedr03()[PhotonIndex] + pho_hcalsumetconedr03()[PhotonIndex] - rho_algo1()*0.17;
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
    if ( pho_pfiso_mycharged02()->at(PhotonIndex).at(0)  > 4 )  return false;
    //    if (pho_isconv()[PhotonIndex]!=1) return false;
    return true;
  } else {
    if (pho_isEB()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.082 || pho_sieie()[PhotonIndex]>0.014)) return false;
    if (pho_isEE()[PhotonIndex] && (pho_hoe()[PhotonIndex]>0.075 || pho_sieie()[PhotonIndex]>0.034)) return false;
    if (EtCorrEcalIso>50.0) return false;
    if (EtCorrHcalIso>50.0) return false;
    if (EtCorrTrkIso>50.0) return false;
    if ( pho_pfiso_mycharged02()->at(PhotonIndex).at(0)  > 4 )  return false;
    // if (pho_isconv()[PhotonIndex]!=1) return false;
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

  kFactor["ZToMuMu"]=1.0;
  return kFactor;

}

map<TString,double> GetWeightsMap(map<TString,double> kFactor, double globalweight) {

  map<TString,double> WeightsMap;
  WeightsMap["None"]=1/globalweight;
  WeightsMap["ZToMuMu"]=kFactor["ZToMuMu"]*1/48819386.0*1871.0;
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

void BookBarrelAndEndcap2D(HistoContainer *histoContainer, TString histname, TString histtitle, int binsa, float lowerlimita, float upperlimita, int binsb, float lowerlimitb, float upperlimitb ) {
  TString Regions[2] = {"Barrel","Endcap"};
  for (int i=0; i<2; i++) {
    TString histnametemp = histname;
    TString histtitletemp = histtitle;
    histnametemp += Regions[i];
    histtitletemp.ReplaceAll("region",Regions[i]);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),binsa,lowerlimita,upperlimita,binsb,lowerlimitb,upperlimitb);
  }
}

void BookBarrelAndEndcapProfiles(HistoContainer *histoContainer, TString histname, TString histtitle, int binsa, float lowerlimita, float upperlimita,float lowerlimitb, float upperlimitb ) {
  TString Regions[2] = {"Barrel","Endcap"};
  for (int i=0; i<2; i++) {
    TString histnametemp = histname;
    TString histtitletemp = histtitle;
    histnametemp += Regions[i];
    histtitletemp.ReplaceAll("region",Regions[i]);
    histoContainer->Add(histnametemp.Data(),histtitletemp.Data(),binsa,lowerlimita,upperlimita,lowerlimitb,upperlimitb);
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
  histoContainer->Add("NearMuDeltaRBeforeCut","Near Muon #DeltaR  #mu;#DeltaR;Counts",30,0,1);
  histoContainer->Add("NearMuDeltaEtaDeltaPhiBeforeCut","Near Muon #Delta#phi vs #Delta#eta   #mu;#Delta#eta;#Delta#phi",30,0.,1.,30,0.,1. );
  histoContainer->Add("NearMuDeltaRAfterCut","Near Muon #DeltaR  #mu;#DeltaR;Counts",30,0,1);
  histoContainer->Add("NearMuDeltaEtaDeltaPhiAfterCut","Near Muon #Delta#phi vs #Delta#eta   #mu;#Delta#eta;#Delta#phi",30,0.,1.,30,0.,1. );

  BookBarrelAndEndcap2D(histoContainer,"pho_pfchargedisogood03_vs_nearMuonDRBeforePhoPresel","pho_pfchargedisogood03_vs_nearMuonDRBeforePhoPresel: region; ;nearMuondR;pfcharedisogood03;",30,0.,1.,20,-10,10);
  BookBarrelAndEndcap2D(histoContainer,"pho_pfchargedisobad03_vs_nearMuonDRBeforePhoPresel","pho_pfchargedisobad03_vs_nearMuonDRBeforePhoPresel: region; ;nearMuondR;pfcharedisobad03;",30,0.,1.,20,-10,10);
  BookBarrelAndEndcap2D(histoContainer,"pho_pfphotoniso03_vs_nearMuonDRBeforePhoPresel","pho_pfphotoniso03_vs_nearMuonDRBeforePhoPresel: region; ;nearMuondR;pfphotoniso03;",30,0.,1.,50,0.,10.);
  BookBarrelAndEndcap2D(histoContainer,"pho_pfneutraliso03_vs_nearMuonDRBeforePhoPresel","pho_pfneutraliso03_vs_nearMuonDRBeforePhoPresel: region; ;nearMuondR;pfneutraliso03;",30,0.,1.,20,-10,10);

  BookCutsAndCategories(histoContainer,"Numvtx","Number of Primary Verticies:Cuts:Cat;Number of Vertices;Counts",100,0,100);
  BookCutsAndCategories(histoContainer,"PhotonPt","Pt of Photon:Cuts:Cat;Pt (GeV);Counts",100,0,200);
  BookCutsAndCategories(histoContainer,"PhotonEt","Et of Photon:Cuts:Cat;Et (GeV);Counts",50,0,100);
  BookCutsAndCategories(histoContainer,"PhotonEta","#eta of Photon:Cuts:Cat;#eta;Counts",32,-3.2,3.2);
  BookCutsAndCategories(histoContainer,"PhotonPhi","#phi of Photon:Cuts:Cat;#phi;Counts",32,-3.2,3.2);
  BookCutsAndCategories(histoContainer,"ZMass","Invariant mass of #mu#mu#gamma system:Cuts:Cat;Mass (GeV);Counts",30,75,105);
  BookCutsAndCategories(histoContainer,"ZMassZoom","Invariant mass of #mu#mu#gamma system:Cuts:Cat;Mass (GeV);Counts",16,87.2,95.2);

  histoContainer->Add("Rho","Event Rho;Rho;Counts",100,0,50);
  
  BookBarrelAndEndcap(histoContainer,"allpho_pt","pho_pt: region; Pt (GeV) ;Counts",100,0.,100.);
  BookBarrelAndEndcap(histoContainer,"pho_pt","pho_pt: region; Pt (GeV) ;Counts",100,0.,100.);
  BookBarrelAndEndcap(histoContainer,"pho_pt_afterMVAcut","pho_pt: region; Pt (GeV) ;Counts",100,0.,100.);
  
  BookBarrelAndEndcap(histoContainer,"pho_pt_afternewMVAcut","pho_pt: region; Pt (GeV) ;Counts",100,0.,100.);
  BookBarrelAndEndcapProfiles(histoContainer,"pho_idmva_vsPt","pho_idmva: region;Pt;idmva value;",100,0.,100.,0.4,0.4);

  BookBarrelAndEndcap(histoContainer,"pho_mitmva","pho_mitmva: region;idmva value;Counts",50,-0.5,0.5);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_r9","pho_tmva_id_photonid_r9: region;tmva_id_photonid_r9;Counts",100,0.2,1.1);

  histoContainer->Add("pho_tmva_photonid_etawidthBarrel","pho_tmva_photonid_etawidth: Barrel;tmva_photonid_etawidth;Counts",50,0,0.02);
  histoContainer->Add("pho_tmva_photonid_etawidthEndcap","pho_tmva_photonid_etawidth: Endcap;tmva_photonid_etawidth;Counts",50,0,0.05);
  histoContainer->Add("pho_tmva_photonid_sieieBarrel","pho_tmva_id_photonid_sieie: Barrel;#sigma#_{ieie};Counts",100,0.004,0.014);
  histoContainer->Add("pho_tmva_photonid_sieieEndcap","pho_tmva_id_photonid_sieie: Endcap;#sigma#_{ieie};Counts",100,0.01,0.04);
  histoContainer->Add("pho_tmva_photonid_sieipBarrel","pho_tmva_id_photonid_sieip: Barrel;cov_{ieip};Counts",100,-0.0002,0.0002);
  histoContainer->Add("pho_tmva_photonid_sieipEndcap","pho_tmva_id_photonid_sieip: Endcap;cov_{ieip};Counts",100,-0.001,0.001);

  histoContainer->Add("pho_tmva_photonid_sceta","pho_tmva_id_photonid_sceta;tmva_id_photonid_sceta;Counts",68,-3.4,3.4);
  
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_phiwidth","pho_tmva_id_photonid_phiwidth: region;tmva_id_photonid_phiwidth;Counts",50,0,0.16);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_s4ratio","pho_tmva_id_photonid_s4ratio: region;tmva_id_photonid_s4ratio;Counts",70,0.3,1.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_lambdaratio","pho_tmva_id_photonid_lambdaratio: region;tmva_id_photonid_lambdaratio;Counts",50,0,1.5);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_sceta","pho_tmva_id_photonid_sceta: region;tmva_id_photonid_sceta;Counts",68,-3.4,3.4);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_eventrho","pho_tmva_id_photonid_eventrho: region;tmva_id_photonid_eventrho;Counts",100,0,60.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_ESEffSigmaRR","pho_tmva_id_photonid_ESEffSigmaRR: region;tmva_id_photonid_ESEffSigmaRR;Counts",100,0.,15.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_pfchargedisogood03","pho_tmva_id_photonid_pfchargedisogood03: region;tmva_id_photonid_pfchargedisogood03;Counts",80,0.,8.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_pfchargedisobad03","pho_tmva_id_photonid_pfchargedisobad03: region;tmva_id_photonid_pfchargedisobad03;Counts",80,0.,20.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_pfphotoniso03","pho_tmva_id_photonid_pfphotoniso03: region;tmva_id_photonid_pfphotoniso03;Counts",50,0,10.);
  BookBarrelAndEndcap(histoContainer,"pho_tmva_photonid_pfneutraliso03","pho_tmva_id_photonid_pfneutraliso03: region;tmva_id_photonid_pfneutraliso03;Counts",50,0,1.5);

}

void MakeFilesAndWeights(TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  if (inputstring.Contains("Data")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_AllData.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("Run2012ABC")) {
    inputfilelist.push_back(pair<string,int> ("Run2012ABC.root",6));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012A13Jul2012v1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012Arecover06Aug2012v1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012B13Jul2012v4.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012CPromptRecov1.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012CPromptRecov2.root",WeightsMap["None"]));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_Run2012PromptCMissing.root",WeightsMap["None"]));
    
  }
  if (inputstring.Contains("Run2012D")) {
    inputfilelist.push_back(pair<string,int> ("Run2012D.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_DoublePhoton_Run2012D-PromptReco-v1.root",WeightsMap["None"]));
  }
  if (inputstring.Contains("ZToMuMu")) {
    inputfilelist.push_back(pair<string,int> ("ZToMuMu.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc1/c/HiggsGammaGamma/ZMuMuGamma/Moriond2013/Reduced/Reduced_DYToMuMu2012.root",WeightsMap["ZToMuMu"]));
  }

}

void MakeFilesAndWeights(string infile, TString inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, map<TString,double> kFactor, map<TString,double> WeightsMap) {

  string outfile;
  if (infile.rfind("/")!=string::npos) outfile=infile.substr(infile.rfind("/")+1);
  else outfile=infile;

  inputfilelist.push_back(pair<string,int> (outfile,1));
  if (inputstring.Contains("Data") || inputstring.Contains("Run2012ABC") || inputstring.Contains("Run2012D")) inputvector.push_back(pair<string,float> (infile,WeightsMap["None"]));
  if (inputstring.Contains("ZToMuMu")) inputvector.push_back(pair<string,float> (infile,WeightsMap["ZToMuMu"]));

}

void MakePileUpWeights(TString inputstring, map<int,double> &PileUpMap) {

  if (inputstring.Contains("ZToMuMu")) {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/ZMuMu_ZToMuMu.h"
  } else {
    #include "ND_Hto2Photons/TreeReaders/interface/PileUpWeights/Dummy.h"
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
