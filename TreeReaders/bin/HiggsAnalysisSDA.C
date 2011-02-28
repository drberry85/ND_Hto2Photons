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
int DetectorPosition(sdaReader *currentTree, unsigned int index);
string MakeFileName(string filename, bool unweighted, bool dataweight);
void ProgressBar(int &percent);
void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData);
void PrintWeights();

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool dataweight = false;
  bool unweighted = false;

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
  if (filesAndWeights.size()==0) {
    cout << "Warning!!!! No valid inputs!!!! Please one of the following: 90GeV, 110GeV, 120GeV, 150GeV, PhotonPlusJet, EMEnriched, DoubleEMEnriched, QCDBEtoE, Born, or Box." << endl;
    cout << "Exiting Program!!!!" << endl;
    return 0;
  }
  
  for (vector<pair<string, int> >::iterator itFilePair = filelist.begin(); itFilePair != filelist.end(); ++itFilePair) {

    string outfilename = MakeFileName(itFilePair->first, unweighted, dataweight);
    
    TFile* outfile = new TFile(outfilename.c_str(),"RECREATE");
    outfile->cd();
    cout << outfilename << " created." << endl;

    TH1F* hNVertices;
    TH1F* hNSimVertices;
    TH1F* hZVertex;
    TH1F* hZSimVertex;
    TH1F* hZconv;
    TH1F* hZconvError;
    TH1F* hZError;
    TH1F* convr;
    TH2F* convrdZ;
    
    TH1F* hNPhotons[2];
    TH1F* hLeadEt[3][2];
    TH1F* hSubLeadEt[3][2];
    TH1F* hLeadEta[2];
    TH1F* hSubLeadEta[2];
    TH1F* hLeadPhi[2];
    TH1F* hSubLeadPhi[2];
    TH1F* hLeadR9[3][2];
    TH1F* hSubLeadR9[3][2];
    TH1F* hLeadHoE[3][2];
    TH1F* hSubLeadHoE[3][2];
    TH1F* hLeadTrkPtSumSolid03[3][2];
    TH1F* hSubLeadTrkPtSumSolid03[3][2];
    TH1F* hLeadEcalPtSumSolid03[3][2];
    TH1F* hSubLeadEcalPtSumSolid03[3][2];
    TH1F* hLeadHcalPtSumSolid03[3][2];
    TH1F* hSubLeadHcalPtSumSolid03[3][2];
    TH1F* hLeadSigmaIetaIeta[3][2];
    TH1F* hSubLeadSigmaIetaIeta[3][2];

    TH1F* hLeadZPV_[3][2];
    TH1F* hSubLeadZPV_[3][2];
    TH1F* hLeadDzPV_[3][2];
    TH1F* hSubLeadDzPV_[3][2];
    TH1F* hLeadDzPVconv[3][2];
    TH1F* hSubLeadDzPVconv[3][2];

    TH1D* h_mass_2gamma[3][2];
    TH1D* h_pt_2gamma[3][2];
    TH1D* h_pz_2gamma[3][2];
    TH1D* h_eta_2gamma[3][2];
    TH1D* h_phi_2gamma[3][2];
    TH1D* h_CosThetaStar[3][2];

    TH1D* h_mass_2gamma_2gold[3][2];
    TH1D* h_pt_2gamma_2gold[3][2];
    TH1D* h_pz_2gamma_2gold[3][2];
    TH1D* h_eta_2gamma_2gold[3][2];
    TH1D* h_phi_2gamma_2gold[3][2];
    TH1D* h_CosThetaStar_2gold[3][2];

    TH1D* h_mass_2gamma_1goodconv[3][2];
    TH1D* h_pt_2gamma_1goodconv[3][2];
    TH1D* h_pz_2gamma_1goodconv[3][2];
    TH1D* h_eta_2gamma_1goodconv[3][2];
    TH1D* h_phi_2gamma_1goodconv[3][2];
    TH1D* h_CosThetaStar_1goodconv[3][2];

    TH1D* h_mass_2gamma_1goodconv_newvertex[3][2];
    TH1D* h_pt_2gamma_1goodconv_newvertex[3][2];
    TH1D* h_pz_2gamma_1goodconv_newvertex[3][2];
    TH1D* h_eta_2gamma_1goodconv_newvertex[3][2];
    TH1D* h_phi_2gamma_1goodconv_newvertex[3][2];
    TH1D* h_CosThetaStar_1goodconv_newvertex[3][2];

    TH1D* h_mass_2gamma_1goodconv_40[3][2];
    TH1D* h_mass_2gamma_1goodconv_50[3][2];
    TH1D* h_mass_2gamma_1goodconv_60[3][2];
    TH1D* h_mass_2gamma_1goodconv_70[3][2];
    
    TH1D* h_mass_2gamma_1goodconv_simvertex[3][2];
    TH1D* h_pt_2gamma_1goodconv_simvertex[3][2];
    TH1D* h_pz_2gamma_1goodconv_simvertex[3][2];
    TH1D* h_eta_2gamma_1goodconv_simvertex[3][2];
    TH1D* h_phi_2gamma_1goodconv_simvertex[3][2];
    TH1D* h_CosThetaStar_1goodconv_simvertex[3][2];
    
    TH1D* h_mass_2gamma_1poorconv[3][2];
    TH1D* h_pt_2gamma_1poorconv[3][2];
    TH1D* h_pz_2gamma_1poorconv[3][2];
    TH1D* h_eta_2gamma_1poorconv[3][2];
    TH1D* h_phi_2gamma_1poorconv[3][2];
    TH1D* h_CosThetaStar_1poorconv[3][2];

    TH1D* h_mass_2gamma_2conv[3][2];
    TH1D* h_pt_2gamma_2conv[3][2];
    TH1D* h_pz_2gamma_2conv[3][2];
    TH1D* h_eta_2gamma_2conv[3][2];
    TH1D* h_phi_2gamma_2conv[3][2];
    TH1D* h_CosThetaStar_2conv[3][2];

    TH1D* h_mass_2gamma_leftover[3][2];
    TH1D* h_pt_2gamma_leftover[3][2];
    TH1D* h_pz_2gamma_leftover[3][2];
    TH1D* h_eta_2gamma_leftover[3][2];
    TH1D* h_phi_2gamma_leftover[3][2];
    TH1D* h_CosThetaStar_leftover[3][2];

    TH1D* h_lead_r9_cat0[3];
    TH1D* h_sublead_r9_cat0[3];
    TH1D* h_lead_r9_cat1[3];
    TH1D* h_sublead_r9_cat1[3];
    TH1D* h_lead_r9_cat2[3];
    TH1D* h_sublead_r9_cat2[3];
    TH1D* h_lead_r9_cat3[3];
    TH1D* h_sublead_r9_cat3[3];
    TH1D* h_lead_r9_cat4[3];
    TH1D* h_sublead_r9_cat4[3];

    TH1D* h_phi_conv[3][2];
    TH1D* h_eta_conv[3][2];
    TH1D* h_pt_conv[3][2];
    TH1D* h_z_conv[3][2];
    TH1D* h_r_conv[3][2];

    TH2F* h2_convVtxRvsZBarrel_[2];

    TH1F * hLeadEtMarco;
    TH1F * hSubLeadEtMarco;
    TH1F * h_mass_Marco;

    TH1F * hLeadEtMarcoCat[4];
    TH1F * hSubLeadEtMarcoCat[4];
    TH1F * h_mass_MarcoCat[4];

    hNVertices = new TH1F("NumberVertices","Number of Reconstructed Vertices;Number of Vertices; Counts",20,0,20);
    hNSimVertices = new TH1F("NumberSimVertices","Number of Simulated Vertices;Number of Sim Vertices; Counts",20,0,20);

    hZVertex = new TH1F("ZVertex","Z of Primary Vertex;Z (cm); Counts",100,-20,20);
    hZSimVertex = new TH1F("ZSimVertex","Z of Simulated Vertex;Z (com); Counts",100,-20,20);

    hZconv = new TH1F("Zconv","Z of Primary Vertex from Conversion;Z (cm); Counts",100,-20,20);
    hZconvError = new TH1F("hZconvError","Distance between the Z of the Primary Vertex from Conversion and Sim Vertex;Z (cm); Counts",100,-1,1);
    hZError = new TH1F("hZError","Distance between the Z of the Primary Vertex and Sim Vertex;Z (cm); Counts",100,-1,1);

    convr = new TH1F("convr","R of conversion; R (cm); Counts",100,0,100);
    convrdZ = new TH2F("convrdZ","Distance between the Z of the Primary Vertex from Conversion and Sim Vertex versus R of Conversion;Z of Primary Vertex from Conversion (cm);R of Conversion (cm)",100, -5, 5, 100, 0, 100);
    
    hLeadEtMarco = new TH1F("leadEtMarco_allEcal","Leading Photon Et with Marco's Cuts, Et (GeV); Counts",30,0.,150.);
    hSubLeadEtMarco = new TH1F("subleadEtMarco_allEcal","Subleading Photon Et with Marco's Cuts, Et (GeV); Counts",30,0.,150.);
    h_mass_Marco = new TH1F("h_mass_Marco_allEcal","Invariant Mass of Photons with Marco's Cuts; Mass (GeV); Counts",20,100.,200.);

    hLeadEtMarcoCat[0] = new TH1F("leadEtMarco_cat0","Leading Photon Et in Catagory 0, Et (GeV); Counts",30,0.,150.);
    hLeadEtMarcoCat[1] = new TH1F("leadEtMarco_cat1","Leading Photon Et in Catagory 1, Et (GeV); Counts",30,0.,150.);
    hLeadEtMarcoCat[2] = new TH1F("leadEtMarco_cat2","Leading Photon Et in Catagory 2, Et (GeV); Counts",30,0.,150.);
    hLeadEtMarcoCat[3] = new TH1F("leadEtMarco_cat3","Leading Photon Et in Catagory 3, Et (GeV); Counts",30,0.,150.);
    
    hSubLeadEtMarcoCat[0] = new TH1F("subleadEtMarco_cat0","SubLeading Photon Et in Catagory 0, Et (GeV); Counts",30,0.,150.);
    hSubLeadEtMarcoCat[1] = new TH1F("subleadEtMarco_cat1","SubLeading Photon Et in Catagory 1, Et (GeV); Counts",30,0.,150.);
    hSubLeadEtMarcoCat[2] = new TH1F("subleadEtMarco_cat2","SubLeading Photon Et in Catagory 2, Et (GeV); Counts",30,0.,150.);
    hSubLeadEtMarcoCat[3] = new TH1F("subleadEtMarco_cat3","SubLeading Photon Et in Catagory 3, Et (GeV); Counts",30,0.,150.);

    h_mass_MarcoCat[0] = new TH1F("h_mass_Marco_cat0","Invariant Mass of Photons in Catagory 0; Mass (GeV); Counts",20,100.,200.);
    h_mass_MarcoCat[1] = new TH1F("h_mass_Marco_cat1","Invariant Mass of Photons in Catagory 1; Mass (GeV); Counts",20,100.,200.);
    h_mass_MarcoCat[2] = new TH1F("h_mass_Marco_cat2","Invariant Mass of Photons in Catagory 2; Mass (GeV); Counts",20,100.,200.);
    h_mass_MarcoCat[3] = new TH1F("h_mass_Marco_cat3","Invariant Mass of Photons in Catagory 3; Mass (GeV); Counts",20,100.,200.);
    
    hNPhotons[0] = new TH1F("hNPhotonsAll","Num of photons in the event: all candidates",20,-0.5,19.5);
    hNPhotons[1] = new TH1F("hNPhotonsSel","Num of photons in the event: selected candidates",20,-0.5,19.5);

    hLeadEt[0][0] = new TH1F("leadPhoEtAll_allEcal","leading photon Et, all candidates: all ECAL",100,0.,500.);
    hLeadEt[1][0] = new TH1F("leadPhoEtAll_Barrel","leading photon Et, all candidates:  Barrel",100,0.,500.);
    hLeadEt[2][0] = new TH1F("leadPhoEtAll_Endcap","leading photon Et, all candidates:  Endcap",100,0.,500.);

    hSubLeadEt[0][0] = new TH1F("subleadPhoEtAll_allEcal","sub-leading photon Et, all candidates: all ECAL",100,0.,500.);
    hSubLeadEt[1][0] = new TH1F("subleadPhoEtAll_Barrel","sub-leading photon Et, all candidates:  Barrel",100,0.,500.);
    hSubLeadEt[2][0] = new TH1F("subleadPhoEtAll_Endcap","sub-leading photon Et, all candidates:  Endcap",100,0.,500.);

    hLeadEta[0] = new TH1F("leadPhoEtaAll_allEcal","leading photon Eta, all candidates: all ECAL",100,-3.,3.);
    hSubLeadEta[0] = new TH1F("subleadPhoEtaAll_allEcal","sub-leading photon Eta, all candidates: all ECAL",100,-3.,3.);

    hLeadPhi[0] = new TH1F("leadPhoPhiAll_allEcal","leading photon Phi, all candidates: all ECAL",32,-3.2,3.2);
    hSubLeadPhi[0] = new TH1F("subleadPhoPhiAll_allEcal","sub-leading photon Phi, all candidates: all ECAL",32,-3.2,3.2);
    
    hLeadR9[0][0] = new TH1F("leadPhoR9All_allEcal","leading photon R9, all candidates: all ECAL",100,0.,1.1);
    hLeadR9[1][0] = new TH1F("leadPhoR9All_Barrel","leading photon R9, all candidates: Barrel",100,0.,1.1);
    hLeadR9[2][0] = new TH1F("leadPhoR9All_Endcap","leading photon R9, all candidates: Endcap",100,0.,1.1);

    hSubLeadR9[0][0] = new TH1F("subleadPhoR9All_allEcal","sub-leading photon R9, all candidates: all ECAL",100,0.,1.1);
    hSubLeadR9[1][0] = new TH1F("subleadPhoR9All_Barrel","sub-leading photon R9, all candidates: Barrel",100,0.,1.1);
    hSubLeadR9[2][0] = new TH1F("subleadPhoR9All_Endcap","sub-leading photon R9, all candidates: Endcap",100,0.,1.1);

    hLeadHoE[0][0] = new TH1F("leadPhoHoEAll_allEcal","leading photon HoE, all candidates: all ECAL",100,0.,.1);
    hLeadHoE[1][0] = new TH1F("leadPhoHoEAll_Barrel","leading photon HoE, all candidates: Barrel",100,0.,.1);
    hLeadHoE[2][0] = new TH1F("leadPhoHoEAll_Endcap","leading photon HoE, all candidates: Endcap",100,0.,.1);

    hSubLeadHoE[0][0] = new TH1F("subleadPhoHoEAll_allEcal","sub-leading photon HoE, all candidates: all ECAL",100,0.,.1);
    hSubLeadHoE[1][0] = new TH1F("subleadPhoHoEAll_Barrel","sub-leading photon HoE, all candidates: Barrel",100,0.,.1);
    hSubLeadHoE[2][0] = new TH1F("subleadPhoHoEAll_Endcap","sub-leading photon HoE, all candidates: Endcap",100,0.,.1);

    hLeadTrkPtSumSolid03[0][0] = new TH1F("leadPhoTrkPtSumSolid03All_allEcal","leading photon trk pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hLeadTrkPtSumSolid03[1][0] = new TH1F("leadPhoTrkPtSumSolid03All_Barrel","leading photon trk pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hLeadTrkPtSumSolid03[2][0] = new TH1F("leadPhoTrkPtSumSolid03All_Endcap","leading photon trk pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hSubLeadTrkPtSumSolid03[0][0] = new TH1F("subleadPhoTrkPtSumSolid03All_allEcal","sub-leading photon trk pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hSubLeadTrkPtSumSolid03[1][0] = new TH1F("subleadPhoTrkPtSumSolid03All_Barrel","sub-leading photon trk pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hSubLeadTrkPtSumSolid03[2][0] = new TH1F("subleadPhoTrkPtSumSolid03All_Endcap","sub-leading photon trk pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hLeadEcalPtSumSolid03[0][0] = new TH1F("leadPhoEcalPtSumSolid03All_allEcal","leading photon ecal pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hLeadEcalPtSumSolid03[1][0] = new TH1F("leadPhoEcalPtSumSolid03All_Barrel","leading photon ecal pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hLeadEcalPtSumSolid03[2][0] = new TH1F("leadPhoEcalPtSumSolid03All_Endcap","leading photon ecal pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hSubLeadEcalPtSumSolid03[0][0] = new TH1F("subleadPhoEcalPtSumSolid03All_allEcal","sub-leading photon ecal pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hSubLeadEcalPtSumSolid03[1][0] = new TH1F("subleadPhoEcalPtSumSolid03All_Barrel","sub-leading photon ecal pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hSubLeadEcalPtSumSolid03[2][0] = new TH1F("subleadPhoEcalPtSumSolid03All_Endcap","sub-leading photon ecal pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hLeadHcalPtSumSolid03[0][0] = new TH1F("leadPhoHcalPtSumSolid03All_allEcal","leading photon hcal pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hLeadHcalPtSumSolid03[1][0] = new TH1F("leadPhoHcalPtSumSolid03All_Barrel","leading photon hcal pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hLeadHcalPtSumSolid03[2][0] = new TH1F("leadPhoHcalPtSumSolid03All_Endcap","leading photon hcal pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hSubLeadHcalPtSumSolid03[0][0] = new TH1F("subleadPhoHcalPtSumSolid03All_allEcal","sub-leading photon hcal pt sum dr=03, all candidates: all ECAL",100,0.,5.0);
    hSubLeadHcalPtSumSolid03[1][0] = new TH1F("subleadPhoHcalPtSumSolid03All_Barrel","sub-leading photon hcal pt sum dr=03, all candidates: Barrel",100,0.,5.0);
    hSubLeadHcalPtSumSolid03[2][0] = new TH1F("subleadPhoHcalPtSumSolid03All_Endcap","sub-leading photon hcal pt sum dr=03, all candidates: Endcap",100,0.,5.0);

    hLeadSigmaIetaIeta[0][0] = new TH1F("leadPhoSigmaIetaIetaAll_allEcal","leading photon SigmaIetaIeta, all candidates: all ECAL",100,0.,0.2);
    hLeadSigmaIetaIeta[1][0] = new TH1F("leadPhoSigmaIetaIetaAll_Barrel","leading photon SigmaIetaIeta, all candidates: Barrel",100,0.,0.2);
    hLeadSigmaIetaIeta[2][0] = new TH1F("leadPhoSigmaIetaIetaAll_Endcap","leading photon SigmaIetaIeta, all candidates: Endcap",100,0.,0.2);

    hSubLeadSigmaIetaIeta[0][0] = new TH1F("subleadPhoSigmaIetaIetaAll_allEcal","sub-leading photon SigmaIetaIeta, all candidates: all ECAL",100,0.,0.2);
    hSubLeadSigmaIetaIeta[1][0] = new TH1F("subleadPhoSigmaIetaIetaAll_Barrel","sub-leading photon SigmaIetaIeta, all candidates: Barrel",100,0.,0.2);
    hSubLeadSigmaIetaIeta[2][0] = new TH1F("subleadPhoSigmaIetaIetaAll_Endcap","sub-leading photon SigmaIetaIeta, all candidates: Endcap",100,0.,0.2);

    hLeadZPV_[0][0]= new TH1F("leadPhoZPVAll_allEcal"," Leading photon Z(PV) (cm), all candidates: all ecal", 100, -50.0, 50.0);
    hLeadZPV_[1][0]= new TH1F("leadPhoZPVAll_Barrel"," Leading photon Z(PV) (cm), all candidates: Barrel ", 100, -50.0, 50.0);
    hLeadZPV_[2][0]= new TH1F("leadPhoZPVAll_Endcap"," Leading photon Z(PV) (cm), all candidates: Endcap ", 100, -50.0, 50.0);

    hSubLeadZPV_[0][0]= new TH1F("subleadPhoZPVAll_allEcal","Subleading Z(PV) (cm), all candidates: all ecal", 100, -50.0, 50.0);
    hSubLeadZPV_[1][0]= new TH1F("subleadPhoZPVAll_Barrel","Subleading Z(PV) (cm), all candidates: Barrel", 100, -50.0, 50.0);
    hSubLeadZPV_[2][0]= new TH1F("subleadPhoZPVAll_Endcap","Subleading Z(PV) (cm), all candidates: Endcap", 100, -50.0, 50.0);

    hLeadDzPV_[0][0]= new TH1F("leadPhoDZPVAll_allEcal"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), all candidates: all ecal", 100, -20, 20);
    hLeadDzPV_[1][0]= new TH1F("leadPhoDZPVAll_Barrel"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), all candidates: Barrel", 100, -20, 20);
    hLeadDzPV_[2][0]= new TH1F("leadPhoDZPVAll_Endcap"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), all candidates", 100, -20, 20);

    hLeadDzPVconv[0][0]= new TH1F("leadPhoDZPVAllconv_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20, 20);
    hLeadDzPVconv[1][0]= new TH1F("leadPhoDZPVAllconv_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20, 20);
    hLeadDzPVconv[2][0]= new TH1F("leadPhoDZPVAllconv_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates", 100, -20, 20);
    
    hSubLeadDzPV_[0][0]= new TH1F("subleadPhoDZPVAll_allEcal","Subleading #Deltaz_{Zpho - Ztrue} (cm), all candidates: all ecal", 100, -20, 20);
    hSubLeadDzPV_[1][0]= new TH1F("subleadPhoDZPVAll_Barrel","Subleading #Deltaz_{Zpho - Ztrue} (cm), all candidates: Barrel", 100, -20, 20);
    hSubLeadDzPV_[2][0]= new TH1F("subleadPhoDZPVAll_Endcap","Subleading #Deltaz_{Zpho - Ztrue} (cm), all candidates: Endcap", 100, -20, 20);

    hSubLeadDzPVconv[0][0]= new TH1F("subleadPhoDZPVAllconv_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20, 20);
    hSubLeadDzPVconv[1][0]= new TH1F("subleadPhoDZPVAllconv_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20, 20);
    hSubLeadDzPVconv[2][0]= new TH1F("subleadPhoDZPVAllconv_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Endcap", 100, -20, 20);
    
    // candidates after pt and eta selection
    hLeadEt[0][1] = new TH1F("leadPhoEtSel_allEcal","leading photon Et, selected candidates: all ECAL",100,0.,500.);
    hLeadEt[1][1] = new TH1F("leadPhoEtSel_Barrel","leading photon Et, selected candidates:  Barrel",100,0.,500.);
    hLeadEt[2][1] = new TH1F("leadPhoEtSel_Endcap","leading photon Et, selected candidates:  Endcap",100,0.,500.);

    hSubLeadEt[0][1] = new TH1F("subleadPhoEtSel_allEcal","sub-leading photon Et, selected candidates: all ECAL",100,0.,500.);
    hSubLeadEt[1][1] = new TH1F("subleadPhoEtSel_Barrel","sub-leading photon Et, selected candidates:  Barrel",100,0.,500.);
    hSubLeadEt[2][1] = new TH1F("subleadPhoEtSel_Endcap","sub-leading photon Et, selected candidates:  Endcap",100,0.,500.);

    hLeadEta[1] = new TH1F("leadPhoEtaSel_allEcal","leading photon Eta, selected candidates: all ECAL",100,-3.,3.);
    hSubLeadEta[1] = new TH1F("subleadPhoEtaSel_allEcal","sub-leading photon Eta, selected candidates: all ECAL",100,-3.,3.);

    hLeadPhi[1] = new TH1F("leadPhoPhiSel_allEcal","leading photon Phi, selected candidates: all ECAL",32,-3.2,3.2);
    hSubLeadPhi[1] = new TH1F("subleadPhoPhiSel_allEcal","sub-leading photon Phi, selected candidates: all ECAL",32,-3.2,3.2);
    
    hLeadR9[0][1] = new TH1F("leadPhoR9Sel_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    hLeadR9[1][1] = new TH1F("leadPhoR9Sel_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    hLeadR9[2][1] = new TH1F("leadPhoR9Sel_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    hSubLeadR9[0][1] = new TH1F("subleadPhoR9Sel_allEcal","sub-leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    hSubLeadR9[1][1] = new TH1F("subleadPhoR9Sel_Barrel","sub-leading photon R9, selected candidates: Barrel",100,0.,1.1);
    hSubLeadR9[2][1] = new TH1F("subleadPhoR9Sel_Endcap","sub-leading photon R9, selected candidates: Endcap",100,0.,1.1);

    hLeadHoE[0][1] = new TH1F("leadPhoHoESel_allEcal","leading photon HoE, selected candidates: all ECAL",100,0.,.1);
    hLeadHoE[1][1] = new TH1F("leadPhoHoESel_Barrel","leading photon HoE, selected candidates: Barrel",100,0.,.1);
    hLeadHoE[2][1] = new TH1F("leadPhoHoESel_Endcap","leading photon HoE, selected candidates: Endcap",100,0.,.1);

    hSubLeadHoE[0][1] = new TH1F("subleadPhoHoESel_allEcal","sub-leading photon HoE, selected candidates: all ECAL",100,0.,.1);
    hSubLeadHoE[1][1] = new TH1F("subleadPhoHoESel_Barrel","sub-leading photon HoE, selected candidates: Barrel",100,0.,.1);
    hSubLeadHoE[2][1] = new TH1F("subleadPhoHoESel_Endcap","sub-leading photon HoE, selected candidates: Endcap",100,0.,.1);

    hLeadTrkPtSumSolid03[0][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_allEcal","leading photon trk pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hLeadTrkPtSumSolid03[1][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_Barrel","leading photon trk pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hLeadTrkPtSumSolid03[2][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_Endcap","leading photon trk pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hSubLeadTrkPtSumSolid03[0][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_allEcal","sub-leading photon trk pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hSubLeadTrkPtSumSolid03[1][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_Barrel","sub-leading photon trk pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hSubLeadTrkPtSumSolid03[2][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_Endcap","sub-leading photon trk pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hLeadEcalPtSumSolid03[0][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_allEcal","leading photon ecal pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hLeadEcalPtSumSolid03[1][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_Barrel","leading photon ecal pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hLeadEcalPtSumSolid03[2][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_Endcap","leading photon ecal pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hSubLeadEcalPtSumSolid03[0][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_allEcal","sub-leading photon ecal pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hSubLeadEcalPtSumSolid03[1][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_Barrel","sub-leading photon ecal pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hSubLeadEcalPtSumSolid03[2][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_Endcap","sub-leading photon ecal pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hLeadHcalPtSumSolid03[0][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_allEcal","leading photon hcal pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hLeadHcalPtSumSolid03[1][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_Barrel","leading photon hcal pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hLeadHcalPtSumSolid03[2][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_Endcap","leading photon hcal pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hSubLeadHcalPtSumSolid03[0][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_allEcal","sub-leading photon hcal pt sum dr=03, selected candidates: all ECAL",100,0.,5.0);
    hSubLeadHcalPtSumSolid03[1][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_Barrel","sub-leading photon hcal pt sum dr=03, selected candidates: Barrel",100,0.,5.0);
    hSubLeadHcalPtSumSolid03[2][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_Endcap","sub-leading photon hcal pt sum dr=03, selected candidates: Endcap",100,0.,5.0);

    hLeadSigmaIetaIeta[0][1] = new TH1F("leadPhoSigmaIetaIetaSel_allEcal","leading photon SigmaIetaIeta, selected candidates: all ECAL",100,0.,0.2);
    hLeadSigmaIetaIeta[1][1] = new TH1F("leadPhoSigmaIetaIetaSel_Barrel","leading photon SigmaIetaIeta, selected candidates: Barrel",100,0.,0.2);
    hLeadSigmaIetaIeta[2][1] = new TH1F("leadPhoSigmaIetaIetaSel_Endcap","leading photon SigmaIetaIeta, selected candidates: Endcap",100,0.,0.2);

    hSubLeadSigmaIetaIeta[0][1] = new TH1F("subleadPhoSigmaIetaIetaSel_allEcal","sub-leading photon SigmaIetaIeta, selected candidates: all ECAL",100,0.,0.2);
    hSubLeadSigmaIetaIeta[1][1] = new TH1F("subleadPhoSigmaIetaIetaSel_Barrel","sub-leading photon SigmaIetaIeta, selected candidates: Barrel",100,0.,0.2);
    hSubLeadSigmaIetaIeta[2][1] = new TH1F("subleadPhoSigmaIetaIetaSel_Endcap","sub-leading photon SigmaIetaIeta, selected candidates: Endcap",100,0.,0.2);

    hLeadZPV_[0][1]= new TH1F("leadPhoZPVSel_allEcal"," Leading photon Z(PV), selected candidates:  all ECAL", 100, -50.0, 50.0);
    hLeadZPV_[1][1]= new TH1F("leadPhoZPVSel_Barrel"," Leading photon Z(PV), selected candidates: Barrel", 100, -50.0, 50.0);
    hLeadZPV_[2][1]= new TH1F("leadPhoZPVSel_Endcap"," Leading photon Z(PV), selected candidates: Endcap", 100, -50.0, 50.0);

    hSubLeadZPV_[0][1]= new TH1F("subleadPhoZPVSel_allEcal","Subleading Z(PV) (cm), selected candidates: all ECAL ", 100, -50.0, 50.0);
    hSubLeadZPV_[1][1]= new TH1F("subleadPhoZPVSel_Barrel","Subleading Z(PV) (cm), selected candidates: Barrel", 100, -50.0, 50.0);
    hSubLeadZPV_[2][1]= new TH1F("subleadPhoZPVSel_Endcap","Subleading Z(PV) (cm), selected candidates: Endcap", 100, -50.0, 50.0);

    hLeadDzPV_[0][1]= new TH1F("leadPhoDZPVSel_allEcal"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), selected candidates:  all ECAL", 100, -20, 20);
    hLeadDzPV_[1][1]= new TH1F("leadPhoDZPVSel_Barrel"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), selected candidates: Barrel", 100, -20, 20);
    hLeadDzPV_[2][1]= new TH1F("leadPhoDZPVSel_Endcap"," Leading photon #Deltaz_{Zpho - Ztrue} (cm), selected candidates: Endcap", 100, -20, 20);

    hLeadDzPVconv[0][1]= new TH1F("leadPhoDZPVSelconv_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20, 20);
    hLeadDzPVconv[1][1]= new TH1F("leadPhoDZPVSelconv_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20, 20);
    hLeadDzPVconv[2][1]= new TH1F("leadPhoDZPVSelconv_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20, 20);
    
    hSubLeadDzPV_[0][1]= new TH1F("subleadPhoDZPVSel_allEcal","Subleading #Deltaz_{Zpho - Ztrue} (cm), selected candidates:  all ECAL", 100, -20, 20);
    hSubLeadDzPV_[1][1]= new TH1F("subleadPhoDZPVSel_Barrel","Subleading #Deltaz_{Zpho - Ztrue} (cm), selected candidates: Barrel", 100, -20, 20);
    hSubLeadDzPV_[2][1]= new TH1F("subleadPhoDZPVSel_Endcap","Subleading #Deltaz_{Zpho - Ztrue} (cm), selected candidates: Endcap", 100, -20, 20);

    hSubLeadDzPVconv[0][1]= new TH1F("subleadPhoDZPVSelconv_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20, 20);
    hSubLeadDzPVconv[1][1]= new TH1F("subleadPhoDZPVSelconv_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20, 20);
    hSubLeadDzPVconv[2][1]= new TH1F("subleadPhoDZPVSelconv_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20, 20);
    
    // diphoton system
    h_mass_2gamma[0][0]        = new TH1D("h_mass_2gammaAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma[0][0]          = new TH1D("h_pt_2gammaAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates", 200, 0., 200.0);
    h_pz_2gamma[0][0]          = new TH1D("h_pz_2gammaAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma[0][0]         = new TH1D("h_eta_2gammaAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma[0][0]         = new TH1D("h_phi_2gammaAllEB","Di-Photon #phi ;#phi(2#gamma) all barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar[0][0]       = new TH1D("h_CosThetaStarAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates", 60, 0., 1.);

    h_mass_2gamma[0][1]        = new TH1D("h_mass_2gammaAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma[0][1]          = new TH1D("h_pt_2gammaAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates", 200, 0., 200.0);
    h_pz_2gamma[0][1]          = new TH1D("h_pz_2gammaAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma[0][1]         = new TH1D("h_eta_2gammaAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma[0][1]         = new TH1D("h_phi_2gammaAllEE","Di-Photon #phi ;#phi(2#gamma) all endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar[0][1]       = new TH1D("h_CosThetaStarAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates", 60, 0., 1.);
    
    h_mass_2gamma[1][0]        = new TH1D("h_mass_2gammaSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma[1][0]          = new TH1D("h_pt_2gammaSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates", 200, 0., 200.0);
    h_pz_2gamma[1][0]          = new TH1D("h_pz_2gammaSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma[1][0]         = new TH1D("h_eta_2gammaSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma[1][0]         = new TH1D("h_phi_2gammaSelEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar[1][0]       = new TH1D("h_CosThetaStarSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates", 60, 0., 1.);

    h_mass_2gamma[1][1]        = new TH1D("h_mass_2gammaSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma[1][1]          = new TH1D("h_pt_2gammaSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates", 200, 0., 200.0);
    h_pz_2gamma[1][1]          = new TH1D("h_pz_2gammaSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma[1][1]         = new TH1D("h_eta_2gammaSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma[1][1]         = new TH1D("h_phi_2gammaSelEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar[1][1]       = new TH1D("h_CosThetaStarSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates", 60, 0., 1.);

    h_mass_2gamma[2][0]        = new TH1D("h_mass_2gammaMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma[2][0]          = new TH1D("h_pt_2gammaMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched barrel candidates", 200, 0., 200.0);
    h_pz_2gamma[2][0]          = new TH1D("h_pz_2gammaMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma[2][0]         = new TH1D("h_eta_2gammaMatchedEB","Di-Photon #eta ;#eta(2#gamma) matched barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma[2][0]         = new TH1D("h_hpi_2gammaMatchedEB","Di-Photon #hpi ;#hpi(2#gamma) matched barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar[2][0]       = new TH1D("h_CosThetaStarMatchedEB","cos#theta^{*};cos#theta^{*} matched barrel candidates", 60, 0., 1.);

    h_mass_2gamma[2][1]        = new TH1D("h_mass_2gammaMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma[2][1]          = new TH1D("h_pt_2gammaMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched endcap candidates", 200, 0., 200.0);
    h_pz_2gamma[2][1]          = new TH1D("h_pz_2gammaMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma[2][1]         = new TH1D("h_eta_2gammaMatchedEE","Di-Photon #eta ;#eta(2#gamma) matched endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma[2][1]         = new TH1D("h_phi_2gammaMatchedEE","Di-Photon #phi ;#phi(2#gamma) matched endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar[2][1]       = new TH1D("h_CosThetaStarMatchedEE","cos#theta^{*};cos#theta^{*} matched endcap candidates", 60, 0., 1.);
    
    h_mass_2gamma_2gold[0][0]        = new TH1D("h_mass_2gammaGoldenAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all golden barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[0][0]          = new TH1D("h_pt_2gammaGoldenAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all golden barrel candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[0][0]          = new TH1D("h_pz_2gammaGoldenAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all golden barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[0][0]         = new TH1D("h_eta_2gammaGoldenAllEB","Di-Photon #eta ;#eta(2#gamma) all golden barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[0][0]         = new TH1D("h_phi_2gammaGoldenAllEB","Di-Photon #phi ;#phi(2#gamma) all golden barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[0][0]       = new TH1D("h_CosThetaStarGoldenAllEB","cos#theta^{*};cos#theta^{*} all golden barrel candidates", 60, 0., 1.);

    h_mass_2gamma_2gold[0][1]        = new TH1D("h_mass_2gammaGoldenAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all golden endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[0][1]          = new TH1D("h_pt_2gammaGoldenAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all golden endcap candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[0][1]          = new TH1D("h_pz_2gammaGoldenAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all golden endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[0][1]         = new TH1D("h_eta_2gammaGoldenAllEE","Di-Photon #eta ;#eta(2#gamma) all golden endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[0][1]         = new TH1D("h_phi_2gammaGoldenAllEE","Di-Photon #phi ;#phi(2#gamma) all golden endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[0][1]       = new TH1D("h_CosThetaStarGoldenAllEE","cos#theta^{*};cos#theta^{*} all golden endcap candidates", 60, 0., 1.);

    h_mass_2gamma_2gold[1][0]        = new TH1D("h_mass_2gammaGoldenSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[1][0]          = new TH1D("h_pt_2gammaGoldenSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden selected barrel candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[1][0]          = new TH1D("h_pz_2gammaGoldenSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden selected barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[1][0]         = new TH1D("h_eta_2gammaGoldenSelEB","Di-Photon #eta ;#eta(2#gamma) golden selected barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[1][0]         = new TH1D("h_phi_2gammaGoldenSelEB","Di-Photon #phi ;#phi(2#gamma) golden selected barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[1][0]       = new TH1D("h_CosThetaStarGoldenSelEB","cos#theta^{*};cos#theta^{*} golden selected barrel candidates", 60, 0., 1.);

    h_mass_2gamma_2gold[1][1]        = new TH1D("h_mass_2gammaGoldenSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[1][1]          = new TH1D("h_pt_2gammaGoldenSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden selected endcap candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[1][1]          = new TH1D("h_pz_2gammaGoldenSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden selected endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[1][1]         = new TH1D("h_eta_2gammaGoldenSelEE","Di-Photon #eta ;#eta(2#gamma) golden selected endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[1][1]         = new TH1D("h_phi_2gammaGoldenSelEE","Di-Photon #phi ;#phi(2#gamma) golden selected endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[1][1]       = new TH1D("h_CosThetaStarGoldenSelEE","cos#theta^{*};cos#theta^{*} golden selected endcap candidates", 60, 0., 1.);

    h_mass_2gamma_2gold[2][0]        = new TH1D("h_mass_2gammaGoldenMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden matched barrel candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[2][0]          = new TH1D("h_pt_2gammaGoldenMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden matched barrel candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[2][0]          = new TH1D("h_pz_2gammaGoldenMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden matched barrel candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[2][0]         = new TH1D("h_eta_2gammaGoldenMatchedEB","Di-Photon #eta ;#eta(2#gamma) golden matched barrel candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[2][0]         = new TH1D("h_phi_2gammaGoldenMatchedEB","Di-Photon #phi ;#phi(2#gamma) golden matched barrel candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[2][0]       = new TH1D("h_CosThetaStarGoldenMatchedEB","cos#theta^{*};cos#theta^{*} golden matched barrel candidates", 60, 0., 1.);

    h_mass_2gamma_2gold[2][1]        = new TH1D("h_mass_2gammaGoldenMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden matched endcap candidates", 80, 80.0, 160.0);
    h_pt_2gamma_2gold[2][1]          = new TH1D("h_pt_2gammaGoldenMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden matched endcap candidates", 200, 0., 200.0);
    h_pz_2gamma_2gold[2][1]          = new TH1D("h_pz_2gammaGoldenMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden matched endcap candidates", 100, -1000., 1000.0);
    h_eta_2gamma_2gold[2][1]         = new TH1D("h_eta_2gammaGoldenMatchedEE","Di-Photon #eta ;#eta(2#gamma) golden matched endcap candidates", 160, -8.0, 8.0);
    h_phi_2gamma_2gold[2][1]         = new TH1D("h_phi_2gammaGoldenMatchedEE","Di-Photon #phi ;#phi(2#gamma) golden matched endcap candidates", 64, -3.2, 3.2);
    h_CosThetaStar_2gold[2][1]       = new TH1D("h_CosThetaStarGoldenMatchedEE","cos#theta^{*};cos#theta^{*} golden matched endcap candidates", 60, 0., 1.);
    
    h_mass_2gamma_1goodconv[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[0][0]          = new TH1D("h_pt_2gamma1goodconvAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[0][0]          = new TH1D("h_pz_2gamma1goodconvAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[0][0]         = new TH1D("h_eta_2gamma1goodconvAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[0][0]         = new TH1D("h_phi_2gamma1goodconvAllEB","Di-Photon #phi ;#phi(2#gamma) all barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[0][0]       = new TH1D("h_CosThetaStar1goodconvAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[0][1]          = new TH1D("h_pt_2gamma1goodconvAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[0][1]          = new TH1D("h_pz_2gamma1goodconvAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[0][1]         = new TH1D("h_eta_2gamma1goodconvAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[0][1]         = new TH1D("h_phi_2gamma1goodconvAllEE","Di-Photon #phi ;#phi(2#gamma) all endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[0][1]       = new TH1D("h_CosThetaStar1goodconvAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[1][0]          = new TH1D("h_pt_2gamma1goodconvSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[1][0]          = new TH1D("h_pz_2gamma1goodconvSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[1][0]         = new TH1D("h_eta_2gamma1goodconvSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[1][0]         = new TH1D("h_phi_2gamma1goodconvSelEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[1][0]       = new TH1D("h_CosThetaStar1goodconvSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[1][1]          = new TH1D("h_pt_2gamma1goodconvSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[1][1]          = new TH1D("h_pz_2gamma1goodconvSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[1][1]         = new TH1D("h_eta_2gamma1goodconvSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[1][1]         = new TH1D("h_phi_2gamma1goodconvSelEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[1][1]       = new TH1D("h_CosThetaStar1goodconvSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[2][0]          = new TH1D("h_pt_2gamma1goodconvMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[2][0]          = new TH1D("h_pz_2gamma1goodconvMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[2][0]         = new TH1D("h_eta_2gamma1goodconvMatchedEB","Di-Photon #eta ;#eta(2#gamma) matched barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[2][0]         = new TH1D("h_phi_2gamma1goodconvMatchedEB","Di-Photon #phi ;#phi(2#gamma) matched barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[2][0]       = new TH1D("h_CosThetaStar1goodconvMatchedEB","cos#theta^{*};cos#theta^{*} matched barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv[2][1]          = new TH1D("h_pt_2gamma1goodconvMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv[2][1]          = new TH1D("h_pz_2gamma1goodconvMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv[2][1]         = new TH1D("h_eta_2gamma1goodconvMatchedEE","Di-Photon #eta ;#eta(2#gamma) matched endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv[2][1]         = new TH1D("h_phi_2gamma1goodconvMatchedEE","Di-Photon #phi ;#phi(2#gamma) matched endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv[2][1]       = new TH1D("h_CosThetaStar1goodconvMatchedEE","cos#theta^{*};cos#theta^{*} matched endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEBnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[0][0]          = new TH1D("h_pt_2gamma1goodconvAllEBnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) all barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[0][0]          = new TH1D("h_pz_2gamma1goodconvAllEBnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) all barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[0][0]         = new TH1D("h_eta_2gamma1goodconvAllEBnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) all barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[0][0]         = new TH1D("h_phi_2gamma1goodconvAllEBnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) all barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[0][0]       = new TH1D("h_CosThetaStar1goodconvAllEBnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} all barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEEnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[0][1]          = new TH1D("h_pt_2gamma1goodconvAllEEnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) all endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[0][1]          = new TH1D("h_pz_2gamma1goodconvAllEEnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) all endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[0][1]         = new TH1D("h_eta_2gamma1goodconvAllEEnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) all endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[0][1]         = new TH1D("h_phi_2gamma1goodconvAllEEnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) all endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[0][1]       = new TH1D("h_CosThetaStar1goodconvAllEEnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} all endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEBnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[1][0]          = new TH1D("h_pt_2gamma1goodconvSelEBnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) selected barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[1][0]          = new TH1D("h_pz_2gamma1goodconvSelEBnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) selected barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[1][0]         = new TH1D("h_eta_2gamma1goodconvSelEBnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) selected barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[1][0]         = new TH1D("h_phi_2gamma1goodconvSelEBnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) selected barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[1][0]       = new TH1D("h_CosThetaStar1goodconvSelEBnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} selected barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEEnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[1][1]          = new TH1D("h_pt_2gamma1goodconvSelEEnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) selected endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[1][1]          = new TH1D("h_pz_2gamma1goodconvSelEEnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) selected endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[1][1]         = new TH1D("h_eta_2gamma1goodconvSelEEnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) selected endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[1][1]         = new TH1D("h_phi_2gamma1goodconvSelEEnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) selected endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[1][1]       = new TH1D("h_CosThetaStar1goodconvSelEEnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} selected endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEBnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[2][0]          = new TH1D("h_pt_2gamma1goodconvMatchedEBnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) matched barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[2][0]          = new TH1D("h_pz_2gamma1goodconvMatchedEBnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) matched barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[2][0]         = new TH1D("h_eta_2gamma1goodconvMatchedEBnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) matched barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[2][0]         = new TH1D("h_phi_2gamma1goodconvMatchedEBnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) matched barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[2][0]       = new TH1D("h_CosThetaStar1goodconvMatchedEBnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} matched barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_newvertex[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEEnewvertex", "Di-photon invariant mass with new vertex;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_newvertex[2][1]          = new TH1D("h_pt_2gamma1goodconvMatchedEEnewvertex","Di-photon P_{T} with new vertex;PT_{2#gamma} (GeV) matched endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_newvertex[2][1]          = new TH1D("h_pz_2gamma1goodconvMatchedEEnewvertex","Di-photon P_{T} with new vertex;Pz_{2#gamma} (GeV) matched endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_newvertex[2][1]         = new TH1D("h_eta_2gamma1goodconvMatchedEEnewvertex","Di-Photon #eta with new vertex;#eta(2#gamma) matched endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_newvertex[2][1]         = new TH1D("h_phi_2gamma1goodconvMatchedEEnewvertex","Di-Photon #phi with new vertex;#phi(2#gamma) matched endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_newvertex[2][1]       = new TH1D("h_CosThetaStar1goodconvMatchedEEnewvertex","cos#theta^{*} with new vertex;cos#theta^{*} matched endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEBsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[0][0]          = new TH1D("h_pt_2gamma1goodconvAllEBsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) all barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[0][0]          = new TH1D("h_pz_2gamma1goodconvAllEBsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) all barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[0][0]         = new TH1D("h_eta_2gamma1goodconvAllEBsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) all barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[0][0]         = new TH1D("h_phi_2gamma1goodconvAllEBsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) all barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[0][0]       = new TH1D("h_CosThetaStar1goodconvAllEBsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} all barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEEsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[0][1]          = new TH1D("h_pt_2gamma1goodconvAllEEsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) all endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[0][1]          = new TH1D("h_pz_2gamma1goodconvAllEEsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) all endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[0][1]         = new TH1D("h_eta_2gamma1goodconvAllEEsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) all endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[0][1]         = new TH1D("h_phi_2gamma1goodconvAllEEsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) all endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[0][1]       = new TH1D("h_CosThetaStar1goodconvAllEEsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} all endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEBsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[1][0]          = new TH1D("h_pt_2gamma1goodconvSelEBsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) selected barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[1][0]          = new TH1D("h_pz_2gamma1goodconvSelEBsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) selected barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[1][0]         = new TH1D("h_eta_2gamma1goodconvSelEBsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) selected barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[1][0]         = new TH1D("h_phi_2gamma1goodconvSelEBsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) selected barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[1][0]       = new TH1D("h_CosThetaStar1goodconvSelEBsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} selected barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEEsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[1][1]          = new TH1D("h_pt_2gamma1goodconvSelEEsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) selected endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[1][1]          = new TH1D("h_pz_2gamma1goodconvSelEEsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) selected endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[1][1]         = new TH1D("h_eta_2gamma1goodconvSelEEsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) selected endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[1][1]         = new TH1D("h_phi_2gamma1goodconvSelEEsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) selected endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[1][1]       = new TH1D("h_CosThetaStar1goodconvSelEEsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} selected endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEBsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[2][0]          = new TH1D("h_pt_2gamma1goodconvMatchedEBsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) matched barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[2][0]          = new TH1D("h_pz_2gamma1goodconvMatchedEBsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) matched barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[2][0]         = new TH1D("h_eta_2gamma1goodconvMatchedEBsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) matched barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[2][0]         = new TH1D("h_phi_2gamma1goodconvMatchedEBsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) matched barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[2][0]       = new TH1D("h_CosThetaStar1goodconvMatchedEBsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} matched barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_simvertex[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEEsimvertex", "Di-photon invariant mass with sim vertex;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1goodconv_simvertex[2][1]          = new TH1D("h_pt_2gamma1goodconvMatchedEEsimvertex","Di-photon P_{T} with sim vertex;PT_{2#gamma} (GeV) matched endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1goodconv_simvertex[2][1]          = new TH1D("h_pz_2gamma1goodconvMatchedEEsimvertex","Di-photon P_{T} with sim vertex;Pz_{2#gamma} (GeV) matched endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1goodconv_simvertex[2][1]         = new TH1D("h_eta_2gamma1goodconvMatchedEEsimvertex","Di-Photon #eta with sim vertex;#eta(2#gamma) matched endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1goodconv_simvertex[2][1]         = new TH1D("h_phi_2gamma1goodconvMatchedEEsimvertex","Di-Photon #phi with sim vertex;#phi(2#gamma) matched endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1goodconv_simvertex[2][1]       = new TH1D("h_CosThetaStar1goodconvMatchedEEsimvertex","cos#theta^{*} with sim vertex;cos#theta^{*} matched endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1goodconv_40[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEB40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_40[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEE40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_40[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEB40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_40[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEE40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_40[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEB40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_40[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEE40", "Di-photon invariant mass with an R cut of 40cm;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);

    h_mass_2gamma_1goodconv_50[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEB50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_50[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEE50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_50[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEB50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_50[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEE50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_50[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEB50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_50[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEE50", "Di-photon invariant mass with an R cut of 50cm;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);

    h_mass_2gamma_1goodconv_60[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEB60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_60[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEE60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_60[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEB60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_60[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEE60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_60[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEB60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_60[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEE60", "Di-photon invariant mass with an R cut of 60cm;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);

    h_mass_2gamma_1goodconv_70[0][0]        = new TH1D("h_mass_2gamma1goodconvAllEB70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_70[0][1]        = new TH1D("h_mass_2gamma1goodconvAllEE70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_70[1][0]        = new TH1D("h_mass_2gamma1goodconvSelEB70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_70[1][1]        = new TH1D("h_mass_2gamma1goodconvSelEE70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_70[2][0]        = new TH1D("h_mass_2gamma1goodconvMatchedEB70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_mass_2gamma_1goodconv_70[2][1]        = new TH1D("h_mass_2gamma1goodconvMatchedEE70", "Di-photon invariant mass with an R cut of 70cm;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);
    
    h_mass_2gamma_1poorconv[0][0]        = new TH1D("h_mass_2gamma1poorconvAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[0][0]          = new TH1D("h_pt_2gamma1poorconvAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[0][0]          = new TH1D("h_pz_2gamma1poorconvAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[0][0]         = new TH1D("h_eta_2gamma1poorconvAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[0][0]         = new TH1D("h_phi_2gamma1poorconvAllEB","Di-Photon #phi ;#phi(2#gamma) all barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[0][0]       = new TH1D("h_CosThetaStar1poorconvAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1poorconv[0][1]        = new TH1D("h_mass_2gamma1poorconvAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[0][1]          = new TH1D("h_pt_2gamma1poorconvAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[0][1]          = new TH1D("h_pz_2gamma1poorconvAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[0][1]         = new TH1D("h_eta_2gamma1poorconvAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[0][1]         = new TH1D("h_phi_2gamma1poorconvAllEE","Di-Photon #phi ;#phi(2#gamma) all endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[0][1]       = new TH1D("h_CosThetaStar1poorconvAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1poorconv[1][0]        = new TH1D("h_mass_2gamma1poorconvSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[1][0]          = new TH1D("h_pt_2gamma1poorconvSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[1][0]          = new TH1D("h_pz_2gamma1poorconvSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[1][0]         = new TH1D("h_eta_2gamma1poorconvSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[1][0]         = new TH1D("h_phi_2gamma1poorconvSelEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[1][0]       = new TH1D("h_CosThetaStar1poorconvSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1poorconv[1][1]        = new TH1D("h_mass_2gamma1poorconvSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[1][1]          = new TH1D("h_pt_2gamma1poorconvSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[1][1]          = new TH1D("h_pz_2gamma1poorconvSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[1][1]         = new TH1D("h_eta_2gamma1poorconvSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[1][1]         = new TH1D("h_phi_2gamma1poorconvSelEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[1][1]       = new TH1D("h_CosThetaStar1poorconvSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1poorconv[2][0]        = new TH1D("h_mass_2gamma1poorconvMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched barrel candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[2][0]          = new TH1D("h_pt_2gamma1poorconvMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched barrel candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[2][0]          = new TH1D("h_pz_2gamma1poorconvMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched barrel candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[2][0]         = new TH1D("h_eta_2gamma1poorconvMatchedEB","Di-Photon #eta ;#eta(2#gamma) matched barrel candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[2][0]         = new TH1D("h_phi_2gamma1poorconvMatchedEB","Di-Photon #phi ;#phi(2#gamma) matched barrel candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[2][0]       = new TH1D("h_CosThetaStar1poorconvMatchedEB","cos#theta^{*};cos#theta^{*} matched barrel candidates with one conversion", 60, 0., 1.);

    h_mass_2gamma_1poorconv[2][1]        = new TH1D("h_mass_2gamma1poorconvMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched endcap candidates with one conversion", 80, 80.0, 160.0);
    h_pt_2gamma_1poorconv[2][1]          = new TH1D("h_pt_2gamma1poorconvMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched endcap candidates with one conversion", 200, 0., 200.0);
    h_pz_2gamma_1poorconv[2][1]          = new TH1D("h_pz_2gamma1poorconvMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched endcap candidates with one conversion", 100, -1000., 1000.0);
    h_eta_2gamma_1poorconv[2][1]         = new TH1D("h_eta_2gamma1poorconvMatchedEE","Di-Photon #eta ;#eta(2#gamma) matched endcap candidates with one conversion", 160, -8.0, 8.0);
    h_phi_2gamma_1poorconv[2][1]         = new TH1D("h_phi_2gamma1poorconvMatchedEE","Di-Photon #phi ;#phi(2#gamma) matched endcap candidates with one conversion", 64, -3.2, 3.2);
    h_CosThetaStar_1poorconv[2][1]       = new TH1D("h_CosThetaStar1poorconvMatchedEE","cos#theta^{*};cos#theta^{*} matched endcap candidates with one conversion", 60, 0., 1.);
    
    h_mass_2gamma_2conv[0][0]        = new TH1D("h_mass_2gamma2convAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[0][0]          = new TH1D("h_pt_2gamma2convAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[0][0]          = new TH1D("h_pz_2gamma2convAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[0][0]         = new TH1D("h_eta_2gamma2convAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[0][0]         = new TH1D("h_phi_2gamma2convAllEB","Di-Photon #phi ;#phi(2#gamma) all barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[0][0]       = new TH1D("h_CosThetaStar2convAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_2conv[0][1]        = new TH1D("h_mass_2gamma2convAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[0][1]          = new TH1D("h_pt_2gamma2convAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[0][1]          = new TH1D("h_pz_2gamma2convAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[0][1]         = new TH1D("h_eta_2gamma2convAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[0][1]         = new TH1D("h_phi_2gamma2convAllEE","Di-Photon #phi ;#phi(2#gamma) all endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[0][1]       = new TH1D("h_CosThetaStar2convAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_2conv[1][0]        = new TH1D("h_mass_2gamma2convSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[1][0]          = new TH1D("h_pt_2gamma2convSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[1][0]          = new TH1D("h_pz_2gamma2convSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[1][0]         = new TH1D("h_eta_2gamma2convSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[1][0]         = new TH1D("h_phi_2gamma2convSelEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[1][0]       = new TH1D("h_CosThetaStar2convSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_2conv[1][1]        = new TH1D("h_mass_2gamma2convSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[1][1]          = new TH1D("h_pt_2gamma2convSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[1][1]          = new TH1D("h_pz_2gamma2convSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[1][1]         = new TH1D("h_eta_2gamma2convSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[1][1]         = new TH1D("h_phi_2gamma2convSelEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[1][1]       = new TH1D("h_CosThetaStar2convSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_2conv[2][0]        = new TH1D("h_mass_2gamma2convMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[2][0]          = new TH1D("h_pt_2gamma2convMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[2][0]          = new TH1D("h_pz_2gamma2convMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[2][0]         = new TH1D("h_eta_2gamma2convMatchedEB","Di-Photon #eta ;#eta(2#gamma) matched barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[2][0]         = new TH1D("h_phi_2gamma2convMatchedEB","Di-Photon #phi ;#phi(2#gamma) matched barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[2][0]       = new TH1D("h_CosThetaStar2convMatchedEB","cos#theta^{*};cos#theta^{*} matched barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_2conv[2][1]        = new TH1D("h_mass_2gamma2convMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_2conv[2][1]          = new TH1D("h_pt_2gamma2convMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_2conv[2][1]          = new TH1D("h_pz_2gamma2convMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_2conv[2][1]         = new TH1D("h_eta_2gamma2convMatchedEE","Di-Photon #eta ;#eta(2#gamma) matched endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_2conv[2][1]         = new TH1D("h_phi_2gamma2convMatchedEE","Di-Photon #phi ;#phi(2#gamma) matched endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_2conv[2][1]       = new TH1D("h_CosThetaStar2convMatchedEE","cos#theta^{*};cos#theta^{*} matched endcap candidates with two conversions", 60, 0., 1.);
    
    h_mass_2gamma_leftover[0][0]        = new TH1D("h_mass_2gammaleftoverAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[0][0]          = new TH1D("h_pt_2gammaleftoverAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[0][0]          = new TH1D("h_pz_2gammaleftoverAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[0][0]         = new TH1D("h_eta_2gammaleftoverAllEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[0][0]         = new TH1D("h_phi_2gammaleftoverAllEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[0][0]       = new TH1D("h_CosThetaStarleftoverAllEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_leftover[0][1]        = new TH1D("h_mass_2gammaleftoverAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[0][1]          = new TH1D("h_pt_2gammaleftoverAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[0][1]          = new TH1D("h_pz_2gammaleftoverAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[0][1]         = new TH1D("h_eta_2gammaleftoverAllEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[0][1]         = new TH1D("h_phi_2gammaleftoverAllEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[0][1]       = new TH1D("h_CosThetaStarleftoverAllEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_leftover[1][0]        = new TH1D("h_mass_2gammaleftoverSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[1][0]          = new TH1D("h_pt_2gammaleftoverSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[1][0]          = new TH1D("h_pz_2gammaleftoverSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[1][0]         = new TH1D("h_eta_2gammaleftoverSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[1][0]         = new TH1D("h_phi_2gammaleftoverSelEB","Di-Photon #phi ;#phi(2#gamma) selected barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[1][0]       = new TH1D("h_CosThetaStarleftoverSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_leftover[1][1]        = new TH1D("h_mass_2gammaleftoverSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[1][1]          = new TH1D("h_pt_2gammaleftoverSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[1][1]          = new TH1D("h_pz_2gammaleftoverSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[1][1]         = new TH1D("h_eta_2gammaleftoverSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[1][1]         = new TH1D("h_phi_2gammaleftoverSelEE","Di-Photon #phi ;#phi(2#gamma) selected endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[1][1]       = new TH1D("h_CosThetaStarleftoverSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_leftover[2][0]        = new TH1D("h_mass_2gammaleftoverMatchedEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched barrel candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[2][0]          = new TH1D("h_pt_2gammaleftoverMatchedEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched barrel candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[2][0]          = new TH1D("h_pz_2gammaleftoverMatchedEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched barrel candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[2][0]         = new TH1D("h_eta_2gammaleftoverMatchedEB","Di-Photon #eta ;#eta(2#gamma) matched barrel candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[2][0]         = new TH1D("h_phi_2gammaleftoverMatchedEB","Di-Photon #phi ;#phi(2#gamma) matched barrel candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[2][0]       = new TH1D("h_CosThetaStarleftoverMatchedEB","cos#theta^{*};cos#theta^{*} matched barrel candidates with two conversions", 60, 0., 1.);

    h_mass_2gamma_leftover[2][1]        = new TH1D("h_mass_2gammaleftoverMatchedEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) matched endcap candidates with two conversions", 80, 80.0, 160.0);
    h_pt_2gamma_leftover[2][1]          = new TH1D("h_pt_2gammaleftoverMatchedEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) matched endcap candidates with two conversions", 200, 0., 200.0);
    h_pz_2gamma_leftover[2][1]          = new TH1D("h_pz_2gammaleftoverMatchedEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) matched endcap candidates with two conversions", 100, -1000., 1000.0);
    h_eta_2gamma_leftover[2][1]         = new TH1D("h_eta_2gammaleftoverMatchedEE","Di-Photon #eta ;#eta(2#gamma) matched endcap candidates with two conversions", 160, -8.0, 8.0);
    h_phi_2gamma_leftover[2][1]         = new TH1D("h_phi_2gammaleftoverMatchedEE","Di-Photon #phi ;#phi(2#gamma) matched endcap candidates with two conversions", 64, -3.2, 3.2);
    h_CosThetaStar_leftover[2][1]       = new TH1D("h_CosThetaStarleftoverMatchedEE","cos#theta^{*};cos#theta^{*} matched endcap candidates with two conversions", 60, 0., 1.);
    
    h_lead_r9_cat0[0] = new TH1D("h_lead_r9_cat0_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_lead_r9_cat0[1] = new TH1D("h_lead_r9_cat0_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_lead_r9_cat0[2] = new TH1D("h_lead_r9_cat0_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_sublead_r9_cat0[0] = new TH1D("h_sublead_r9_cat0_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_sublead_r9_cat0[1] = new TH1D("h_sublead_r9_cat0_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_sublead_r9_cat0[2] = new TH1D("h_sublead_r9_cat0_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_lead_r9_cat1[0] = new TH1D("h_lead_r9_cat1_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_lead_r9_cat1[1] = new TH1D("h_lead_r9_cat1_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_lead_r9_cat1[2] = new TH1D("h_lead_r9_cat1_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_sublead_r9_cat1[0] = new TH1D("h_sublead_r9_cat1_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_sublead_r9_cat1[1] = new TH1D("h_sublead_r9_cat1_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_sublead_r9_cat1[2] = new TH1D("h_sublead_r9_cat1_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_lead_r9_cat2[0] = new TH1D("h_lead_r9_cat2_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_lead_r9_cat2[1] = new TH1D("h_lead_r9_cat2_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_lead_r9_cat2[2] = new TH1D("h_lead_r9_cat2_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_sublead_r9_cat2[0] = new TH1D("h_sublead_r9_cat2_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_sublead_r9_cat2[1] = new TH1D("h_sublead_r9_cat2_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_sublead_r9_cat2[2] = new TH1D("h_sublead_r9_cat2_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_lead_r9_cat3[0] = new TH1D("h_lead_r9_cat3_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_lead_r9_cat3[1] = new TH1D("h_lead_r9_cat3_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_lead_r9_cat3[2] = new TH1D("h_lead_r9_cat3_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_sublead_r9_cat3[0] = new TH1D("h_sublead_r9_cat3_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_sublead_r9_cat3[1] = new TH1D("h_sublead_r9_cat3_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_sublead_r9_cat3[2] = new TH1D("h_sublead_r9_cat3_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_lead_r9_cat4[0] = new TH1D("h_lead_r9_cat4_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_lead_r9_cat4[1] = new TH1D("h_lead_r9_cat4_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_lead_r9_cat4[2] = new TH1D("h_lead_r9_cat4_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_sublead_r9_cat4[0] = new TH1D("h_sublead_r9_cat4_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
    h_sublead_r9_cat4[1] = new TH1D("h_sublead_r9_cat4_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
    h_sublead_r9_cat4[2] = new TH1D("h_sublead_r9_cat4_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

    h_phi_conv[0][0] = new TH1D("h_phi_conv_All_AllECAL","#phi of Photon Conversion All ECAL; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[0][1] = new TH1D("h_phi_conv_Sel_AllECAL","#phi of Selected Photon Conversion All ECAL; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[1][0] = new TH1D("h_phi_conv_All_Barrel","#phi of Photon Conversion Barrel; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[1][1] = new TH1D("h_phi_conv_Sel_Barrel","#phi of Selected Photon Conversion Barrel; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[2][0] = new TH1D("h_phi_conv_All_Endcap","#phi of Photon Conversion Endcap; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[2][1] = new TH1D("h_phi_conv_Sel_Endcap","#phi of Selected Photon Conversion Endcap; #phi of Conversion", 64, -3.2, 3.2);

    h_eta_conv[0][0] = new TH1D("h_eta_conv_All_AllECAL","#eta of Photon Conversion All ECAL; #eta of Conversion", 60, -3.0, 3.0);
    h_eta_conv[0][1] = new TH1D("h_eta_conv_Sel_AllECAL","#eta of Selected Photon Conversion All ECAL; #eta of Conversion", 60, -3.0, 3.0);
    h_eta_conv[1][0] = new TH1D("h_eta_conv_All_Barrel","#eta of Photon Conversion Barrel; #eta of Conversion", 60, -3.0, 3.0);
    h_eta_conv[1][1] = new TH1D("h_eta_conv_Sel_Barrel","#eta of Selected Photon Conversion Barrel; #eta of Conversion", 60, -3.0, 3.0);
    h_eta_conv[2][0] = new TH1D("h_eta_conv_All_Endcap","#eta of Photon Conversion Endcap; #eta of Conversion", 60, -3.0, 3.0);
    h_eta_conv[2][1] = new TH1D("h_eta_conv_Sel_Endcap","#eta of Selected Photon Conversion Endcap; #eta of Conversion", 60, -3.0, 3.0);

    h_pt_conv[0][0] = new TH1D("h_pt_conv_All_AllECAL","pt of Photon Conversion All ECAL; pt of Conversion", 200, 0, 200);
    h_pt_conv[0][1] = new TH1D("h_pt_conv_Sel_AllECAL","pt of Selected Photon Conversion All ECAL; pt of Conversion", 200, 0, 200);
    h_pt_conv[1][0] = new TH1D("h_pt_conv_All_Barrel","pt of Photon Conversion Barrel; pt of Conversion", 200, 0, 200);
    h_pt_conv[1][1] = new TH1D("h_pt_conv_Sel_Barrel","pt of Selected Photon Conversion Barrel; pt of Conversion", 200, 0, 200);
    h_pt_conv[2][0] = new TH1D("h_pt_conv_All_Endcap","pt of Photon Conversion Endcap; pt of Conversion", 200, 0, 200);
    h_pt_conv[2][1] = new TH1D("h_pt_conv_Sel_Endcap","pt of Selected Photon Conversion Endcap; pt of Conversion", 200, 0, 200);

    h_z_conv[0][0] = new TH1D("h_z_conv_All_AllECAL","z of Photon Conversion All ECAL; z of Conversion", 200, -200, 200);
    h_z_conv[0][1] = new TH1D("h_z_conv_Sel_AllECAL","z of Selected Photon Conversion All ECAL; z of Conversion", 200, -200, 200);
    h_z_conv[1][0] = new TH1D("h_z_conv_All_Barrel","z of Photon Conversion Barrel; z of Conversion", 200, -200, 200);
    h_z_conv[1][1] = new TH1D("h_z_conv_Sel_Barrel","z of Selected Photon Conversion Barrel; z of Conversion", 200, -200, 200);
    h_z_conv[2][0] = new TH1D("h_z_conv_All_Endcap","z of Photon Conversion Endcap; z of Conversion", 200, -200, 200);
    h_z_conv[2][1] = new TH1D("h_z_conv_Sel_Endcap","z of Selected Photon Conversion Endcap; z of Conversion", 200, -200, 200);
  
    h_r_conv[0][0] = new TH1D("h_r_conv_All_AllECAL","r of Photon Conversion All ECAL; r of Conversion", 100, 0, 100);
    h_r_conv[0][1] = new TH1D("h_r_conv_Sel_AllECAL","r of Selected Photon Conversion All ECAL; r of Conversion", 100, 0, 100);
    h_r_conv[1][0] = new TH1D("h_r_conv_All_Barrel","r of Photon Conversion Barrel; r of Conversion", 100, 0, 100);
    h_r_conv[1][1] = new TH1D("h_r_conv_Sel_Barrel","r of Selected Photon Conversion Barrel; r of Conversion", 100, 0, 100);
    h_r_conv[2][0] = new TH1D("h_r_conv_All_Endcap","r of Photon Conversion Endcap; r of Conversion", 100, 0, 100);
    h_r_conv[2][1] = new TH1D("h_r_conv_Sel_Endcap","r of Selected Photon Conversion Endcap; r of Conversion", 100, 0, 100);
    
    h2_convVtxRvsZBarrel_[0] =   new TH2F("convVtxRvsZBarrelAll"," Photon  conversion vtx position all candidates Barrel",200, 0., 280., 200, 0., 80.);
    h2_convVtxRvsZBarrel_[1] =   new TH2F("convVtxRvsZBarrelSel"," Photon  conversion vtx position selected candidates Barrel",200, 0., 280., 200, 0., 80.);

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

        int iLeadDetector = 0;
        int iSubleadDetector = 0;

        TVector3 conversionVertex;

        if (i % (nentries/100) == 0 && bar) {
          if (percent == 100) percent = 99;
          ProgressBar(percent);
        }

        currentTree.GetEntry(i);

        TVector3 SimVertex = *((TVector3*) currentTree.simvtx->At(0));
        vector<TVector3> PrimaryVertex;
        vector<TVector3> ConversionVertex;
        vector<TVector3> ConversionPairMomentum;
        vector<TVector3> Photonxyz;
        vector<TLorentzVector> Photonp4;

        for (unsigned int j=0; j!=(unsigned int) currentTree.vtx_std_xyz->GetSize(); j++) PrimaryVertex.push_back(*((TVector3*) currentTree.vtx_std_xyz->At(j)));
        
        hZVertex->Fill(PrimaryVertex[0].Z());
        hZSimVertex->Fill(SimVertex.Z());

        if (currentTree.pho_n<1) continue;

        map <double, unsigned int> ptindex;
        for (unsigned int j=0; j!=(unsigned int) currentTree.pho_n; j++) {
          ConversionVertex.push_back(*((TVector3*) currentTree.pho_conv_vtx->At(j)));
          ConversionPairMomentum.push_back(*((TVector3*) currentTree.pho_conv_pair_momentum->At(j)));
          Photonxyz.push_back(*((TVector3*) currentTree.pho_calopos->At(j)));
          Photonp4.push_back(*((TLorentzVector*) currentTree.pho_p4->At(j)));
          ptindex[Photonp4[j].Pt()]=j;
        }

        double leadpt = -1;
        double subleadpt = -1;
        unsigned int leadindex = 0;
        unsigned int subleadindex = 0;

        bool sorted = sortpt(ptindex, currentTree.pho_n, leadpt, subleadpt, leadindex, subleadindex);

        if (!sorted) cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        //cout << "Final Lead Index: " << leadindex  << " (" << Photonp4[leadindex].Pt() << ") Sublead Index: " << subleadindex << " (" << Photonp4[subleadindex].Pt() << ") Number of Photons: " << currentTree.pho_n << endl;
        hNVertices->Fill(currentTree.vtx_std_xyz->GetSize());
        hNSimVertices->Fill(currentTree.simvtx->GetSize());

        //////////////// basic selection
        if (!preselection(Photonp4[leadindex].Pt(), Photonp4[subleadindex].Pt(), Photonxyz[leadindex].Eta(), Photonxyz[subleadindex].Eta(), currentTree.pho_isEBEEGap[leadindex], currentTree.pho_isEBEEGap[subleadindex])) continue;

        ////////////////////////////////////
        iLeadDetector = DetectorPosition(&currentTree, leadindex);
        iSubleadDetector = DetectorPosition(&currentTree, subleadindex);

        /////////////////////////
        int leadPhoCategory = photonCategory ( currentTree.pho_haspixseed[leadindex], currentTree.pho_r9[leadindex],  currentTree.pho_conv_ntracks[leadindex], currentTree.pho_conv_chi2_probability[leadindex] , Photonp4[leadindex].Pt()/ConversionPairMomentum[leadindex].Perp());
        int subleadPhoCategory = photonCategory (currentTree.pho_haspixseed[subleadindex], currentTree.pho_r9[subleadindex],  currentTree.pho_conv_ntracks[subleadindex], currentTree.pho_conv_chi2_probability[subleadindex] ,  Photonp4[subleadindex].Pt()/ConversionPairMomentum[subleadindex].Perp());
        int diPhoCategory = diPhotonCategory( leadPhoCategory, subleadPhoCategory );
        ////////////////////////////////////

        unsigned int convindex = 0;
        if (leadPhoCategory==2 || subleadPhoCategory==2) {
          if (!currentTree.pho_isEB[convindex]) continue;
          if (currentTree.pho_conv_validvtx[subleadindex]) convindex=subleadindex;
          if (currentTree.pho_conv_validvtx[leadindex]) convindex=leadindex;
          convr->Fill(ConversionVertex[convindex].Perp());
          convrdZ->Fill(currentTree.pho_conv_zofprimvtxfromtrks[convindex]-SimVertex.Z(),ConversionVertex[convindex].Perp());
          hZconv->Fill(currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          hZconvError->Fill(SimVertex.Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          hZError->Fill(SimVertex.Z()-PrimaryVertex[0].Z());
        }

        for (unsigned int j=0; j<(unsigned int) currentTree.pho_n; j++) {
          unsigned int detector = 1;
          if (!currentTree.pho_conv_validvtx[j]) continue;
          if (currentTree.pho_isEB[j]) detector=1;
          if (currentTree.pho_isEE[j]) detector=2;

          h_phi_conv[0][0]->Fill(Photonxyz[j].Phi());
          h_eta_conv[0][0]->Fill(Photonxyz[j].Eta());
          h_pt_conv[0][0]->Fill(Photonp4[j].Pt());
          h_z_conv[0][0]->Fill(ConversionVertex[j].z());
          h_r_conv[0][0]->Fill(ConversionVertex[j].Perp());

          h_phi_conv[detector][0]->Fill(Photonxyz[j].Phi());
          h_eta_conv[detector][0]->Fill(Photonxyz[j].Eta());
          h_pt_conv[detector][0]->Fill(Photonp4[j].Pt());
          h_z_conv[detector][0]->Fill(ConversionVertex[j].z());
          h_r_conv[detector][0]->Fill(ConversionVertex[j].Perp());
        }

        hLeadEt[0][0]->Fill(Photonp4[leadindex].Et(),weight);
        hLeadEta[0]->Fill(Photonxyz[leadindex].Eta(),weight);
        hLeadPhi[0]->Fill(Photonxyz[leadindex].Phi(),weight);
        hLeadR9[0][0]->Fill(currentTree.pho_r9[leadindex],weight);
        hLeadHoE[0][0]->Fill(currentTree.pho_hoe[leadindex],weight);
        hLeadTrkPtSumSolid03[0][0]->Fill(currentTree.pho_trksumptsolidconedr03[leadindex],weight);
        hLeadEcalPtSumSolid03[0][0]->Fill(currentTree.pho_ecalsumetconedr03[leadindex],weight);
        hLeadHcalPtSumSolid03[0][0]->Fill(currentTree.pho_hcalsumetconedr03[leadindex],weight);
        hLeadSigmaIetaIeta[0][0]->Fill(currentTree.pho_sieie[leadindex],weight);
        hLeadZPV_[0][0]->Fill(PrimaryVertex[0].Z(),weight);
        hLeadDzPV_[0][0]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==leadindex && leadPhoCategory==2) hLeadDzPVconv[0][0]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[leadindex]-SimVertex.Z(),weight);
        if (currentTree.pho_conv_ntracks[leadindex]==2 && currentTree.pho_conv_chi2_probability[leadindex]>0.0005 && (bool) currentTree.pho_isEB[leadindex]) {
          conversionVertex.SetXYZ(ConversionVertex[leadindex].X(),ConversionVertex[leadindex].Y(),ConversionVertex[leadindex].Z());
          h2_convVtxRvsZBarrel_[0]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }
        
        hNPhotons[0]->Fill(currentTree.pho_n,weight);

        hLeadEt[iLeadDetector][0]->Fill(Photonp4[leadindex].Et(),weight);
        hLeadR9[iLeadDetector][0]->Fill(currentTree.pho_r9[leadindex],weight);
        hLeadHoE[iLeadDetector][0]->Fill(currentTree.pho_hoe[leadindex],weight);
        hLeadTrkPtSumSolid03[iLeadDetector][0]->Fill(currentTree.pho_trksumptsolidconedr03[leadindex],weight);
        hLeadEcalPtSumSolid03[iLeadDetector][0]->Fill(currentTree.pho_ecalsumetconedr03[leadindex],weight);
        hLeadHcalPtSumSolid03[iLeadDetector][0]->Fill(currentTree.pho_hcalsumetconedr03[leadindex],weight);
        hLeadSigmaIetaIeta[iLeadDetector][0]->Fill(currentTree.pho_sieie[leadindex],weight);
        hLeadZPV_[iLeadDetector][0]->Fill(PrimaryVertex[0].Z(),weight);
        hLeadDzPV_[iLeadDetector][0]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==leadindex && leadPhoCategory==2) hLeadDzPVconv[iLeadDetector][0]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[leadindex]-SimVertex.Z(),weight);

        hSubLeadEt[0][0]->Fill(Photonp4[subleadindex].Et(),weight);
        hSubLeadEta[0]->Fill(Photonxyz[subleadindex].Eta(),weight);
        hSubLeadPhi[0]->Fill(Photonxyz[subleadindex].Phi(),weight);
        hSubLeadR9[0][0]->Fill(currentTree.pho_r9[subleadindex],weight);
        hSubLeadHoE[0][0]->Fill(currentTree.pho_hoe[subleadindex],weight);
        hSubLeadTrkPtSumSolid03[0][0]->Fill(currentTree.pho_trksumptsolidconedr03[subleadindex],weight);
        hSubLeadEcalPtSumSolid03[0][0]->Fill(currentTree.pho_ecalsumetconedr03[subleadindex],weight);
        hSubLeadHcalPtSumSolid03[0][0]->Fill(currentTree.pho_hcalsumetconedr03[subleadindex],weight);
        hSubLeadSigmaIetaIeta[0][0]->Fill(currentTree.pho_sieie[subleadindex],weight);
        hSubLeadZPV_[0][0]->Fill(PrimaryVertex[0].Z(),weight);
        hSubLeadDzPV_[0][0]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==subleadindex && subleadPhoCategory==2) hSubLeadDzPVconv[0][0]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[subleadindex]-SimVertex.Z(),weight);
        if (currentTree.pho_conv_ntracks[subleadindex]==2 && currentTree.pho_conv_chi2_probability[subleadindex]>0.0005 && (bool) currentTree.pho_isEB[subleadindex]) {
          conversionVertex.SetXYZ(ConversionVertex[leadindex].X(),ConversionVertex[leadindex].Y(),ConversionVertex[leadindex].Z());
          h2_convVtxRvsZBarrel_[0]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        /*if (Photonp4[subleadindex].Pt()>Photonp4[leadindex].Pt()) {
          cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
          cout << "LeadPt is " << Photonp4[leadindex].Pt() << endl;
          cout << "SubLeadPt is " << Photonp4[subleadindex].Pt() << endl;
        }*/

        hSubLeadEt[iSubleadDetector][0]->Fill(Photonp4[subleadindex].Et(),weight);
        hSubLeadR9[iSubleadDetector][0]->Fill(currentTree.pho_r9[subleadindex],weight);
        hSubLeadHoE[iSubleadDetector][0]->Fill(currentTree.pho_hoe[subleadindex],weight);
        hSubLeadTrkPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.pho_trksumptsolidconedr03[subleadindex],weight);
        hSubLeadEcalPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.pho_ecalsumetconedr03[subleadindex],weight);
        hSubLeadHcalPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.pho_hcalsumetconedr03[subleadindex],weight);
        hSubLeadSigmaIetaIeta[iSubleadDetector][0]->Fill(currentTree.pho_sieie[subleadindex],weight);
        hSubLeadZPV_[iSubleadDetector][0]->Fill(PrimaryVertex[0].Z(),weight);
        hSubLeadDzPV_[iSubleadDetector][0]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==subleadindex && subleadPhoCategory==2) hSubLeadDzPVconv[iSubleadDetector][0]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[subleadindex]-SimVertex.Z(),weight);

        //TLorentzVector VLead( currentTree.momentumX[leadindex],   currentTree.momentumY[leadindex],  currentTree.momentumZ[leadindex], currentTree.energy[leadindex]);
        //TLorentzVector VSubLead( currentTree.momentumX[subleadindex],   currentTree.momentumY[subleadindex],  currentTree.momentumZ[subleadindex], currentTree.energy[subleadindex]);
        TLorentzVector VLead(Photonp4[leadindex]);
        TLorentzVector VSubLead(Photonp4[subleadindex]);
        TLorentzVector VSum=VLead+VSubLead;
        double InvMass=fabs(VSum.M());

        // calculate Cos_theta_star 
        double beta_b = VSum.Beta();
        double gamma_b = VSum.Gamma();
        TVector3 directionV = VSum.Vect().Unit();

        // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
        TVector3 CrossVLead = VLead.Vect().Cross(directionV);
        double DotVLeadValue = VLead.Vect().Dot(directionV);
        double CrossVLeadValue = sqrt(CrossVLead.x()*CrossVLead.x()+CrossVLead.y()*CrossVLead.y()+CrossVLead.z()*CrossVLead.z());
        double sin_theta =  CrossVLeadValue/VLead.E();
        double cos_theta =  DotVLeadValue/VLead.E(); 
        double tg_thetas = sin_theta/(gamma_b*(cos_theta-beta_b));
        double cos_thetastar = 1.0/sqrt(1.0+tg_thetas*tg_thetas);

        double leadtheta_newvertex = 0;
        double leadeta_newvertex = 0;
        double leadpt_newvertex = 0;
        
        double subleadtheta_newvertex = 0;
        double subleadeta_newvertex = 0;
        double subleadpt_newvertex = 0;
        
        TLorentzVector VSum_newvertex(0,0,0,0);
        double InvMass_newvertex = 0;

        // calculate Cos_theta_star 
        double beta_b_newvertex = 0;
        double gamma_b_newvertex = 0;
        TVector3 directionV_newvertex(0,0,0);

        // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
        TVector3 CrossVLead_newvertex(0,0,0);
        double DotVLeadValue_newvertex = 0;
        double CrossVLeadValue_newvertex = 0;
        double sin_theta_newvertex = 0;
        double cos_theta_newvertex = 0;
        double tg_thetas_newvertex = 0;
        double cos_thetastar_newvertex = 0;
        
        //With New Vertex
        TLorentzVector VLead_newvertex(0,0,0,0);
        TLorentzVector VSubLead_newvertex(0,0,0,0);

        double leadtheta_simvertex = 0;
        double leadeta_simvertex = 0;
        double leadpt_simvertex = 0;
        
        double subleadtheta_simvertex = 0;
        double subleadeta_simvertex = 0;
        double subleadpt_simvertex = 0;
        
        TLorentzVector VSum_simvertex(0,0,0,0);
        double InvMass_simvertex = 0;

        // calculate Cos_theta_star 
        double beta_b_simvertex = 0;
        double gamma_b_simvertex = 0;
        TVector3 directionV_simvertex(0,0,0);

        // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
        TVector3 CrossVLead_simvertex(0,0,0);
        double DotVLeadValue_simvertex = 0;
        double CrossVLeadValue_simvertex = 0;
        double sin_theta_simvertex = 0;
        double cos_theta_simvertex = 0;
        double tg_thetas_simvertex = 0;
        double cos_thetastar_simvertex = 0;

        //With Sim Vertex
        TLorentzVector VLead_simvertex(0,0,0,0);
        TLorentzVector VSubLead_simvertex(0,0,0,0);
        
        if (diPhoCategory==2) {

          leadtheta_newvertex = atan2(Photonxyz[leadindex].Perp(),Photonxyz[leadindex].Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          leadeta_newvertex = -log(tan(leadtheta_newvertex/2));
          leadpt_newvertex = fabs(Photonp4[leadindex].E()*sin(leadtheta_newvertex));

          //cout "Lorentz Theta: " << Photonxyz[leadindex].Theta() << " Lorentz Eta: " << Photonxyz[leadindex].Eta() << " Lorentz Pt: " << Photonp4[leadindex].Pt() << endl;
          //cout "My Theta:      " << leadtheta_newvertex << " My Eta       " << leadeta_newvertex << " My Pt      " << leadpt_newvertex << endl;
        
          subleadtheta_newvertex =  atan2(Photonxyz[subleadindex].Perp(),Photonxyz[subleadindex].Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          subleadeta_newvertex = -log(tan(subleadtheta_newvertex/2));
          subleadpt_newvertex = fabs(Photonp4[subleadindex].E()*sin(subleadtheta_newvertex));
          //cout << "Lorentz Theta: " << Photonxyz[subleadindex].Theta() << " Lorentz Eta: " << Photonxyz[subleadindex].Eta() << " Lorentz Pt: " << Photonp4[subleadindex].Pt() << endl;
          //cout << "My Theta:      " << subleadtheta_newvertex << " My Eta       " << subleadeta_newvertex << " My Pt      " << subleadpt_newvertex << endl;
      
          VLead_newvertex.SetPtEtaPhiE(leadpt_newvertex,leadeta_newvertex,Photonxyz[leadindex].Phi(),Photonp4[leadindex].E());
          VSubLead_newvertex.SetPtEtaPhiE(subleadpt_newvertex,subleadeta_newvertex,Photonxyz[subleadindex].Phi(),Photonp4[subleadindex].E());
          //cout << "LeadPt: " << VLead.Pt() << " LeadEta: " << VLead.Eta() << " LeadPhi: " << VLead.Phi() << " LeadE: " << VLead.E() << endl;
          //cout << "MyLeadPt: " << VLead_newvertex.Pt() << " MyLeadEta: " << VLead_newvertex.Eta() << " MyLeadPhi: " << VLead_newvertex.Phi() << " MyLeadE: " << VLead_newvertex.E() << endl;
          //cout << "SubLeadPt: " << VSubLead.Pt() << " SubLeadEta: " << VSubLead.Eta() << " SubLeadPhi: " << VSubLead.Phi() << " SubLeadE: " << VSubLead.E() << endl;
          //cout << "MySubLeadPt: " << VSubLead_newvertex.Pt() << " MySubLeadEta: " << VSubLead_newvertex.Eta() << " MySubLeadPhi: " << VSubLead_newvertex.Phi() << " MySubLeadE: " << VSubLead_newvertex.E() << endl;
          VSum_newvertex=VLead_newvertex+VSubLead_newvertex;
          InvMass_newvertex=fabs(VSum_newvertex.M());

          // calculate Cos_theta_star 
          beta_b_newvertex  = VSum_newvertex.Beta();
          gamma_b_newvertex = VSum_newvertex.Gamma();
          directionV_newvertex= VSum_newvertex.Vect().Unit();

          // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
          CrossVLead_newvertex=VLead_newvertex.Vect().Cross(directionV);
          DotVLeadValue_newvertex=VLead_newvertex.Vect().Dot(directionV);
          CrossVLeadValue_newvertex=sqrt(CrossVLead_newvertex.x()*CrossVLead_newvertex.x()+CrossVLead_newvertex.y()*CrossVLead_newvertex.y()+CrossVLead_newvertex.z()*CrossVLead_newvertex.z());
          sin_theta_newvertex =  CrossVLeadValue_newvertex/VLead_newvertex.E();
          cos_theta_newvertex =  DotVLeadValue_newvertex/VLead_newvertex.E(); 
          tg_thetas_newvertex = sin_theta_newvertex/(gamma_b_newvertex*(cos_theta_newvertex-beta_b_newvertex));
          cos_thetastar_newvertex= 1.0/sqrt(1.0+tg_thetas_newvertex*tg_thetas_newvertex);

          leadtheta_simvertex = atan2(Photonxyz[leadindex].Perp(),Photonxyz[leadindex].Z()-SimVertex.Z());
          leadeta_simvertex = -log(tan(leadtheta_simvertex/2));
          leadpt_simvertex = fabs(Photonp4[leadindex].E()*sin(leadtheta_simvertex));

          //cout "Lorentz Theta: " << Photonxyz[leadindex].Theta() << " Lorentz Eta: " << Photonxyz[leadindex].Eta() << " Lorentz Pt: " << Photonp4[leadindex].Pt() << endl;
          //cout "My Theta:      " << leadtheta_simvertex << " My Eta       " << leadeta_simvertex << " My Pt      " << leadpt_simvertex << endl;
        
          subleadtheta_simvertex =  atan2(Photonxyz[subleadindex].Perp(),Photonxyz[subleadindex].Z()-SimVertex.Z());
          subleadeta_simvertex = -log(tan(subleadtheta_simvertex/2));
          subleadpt_simvertex = fabs(Photonp4[subleadindex].E()*sin(subleadtheta_simvertex));
          //cout << "Lorentz Theta: " << Photonxyz[subleadindex].Theta() << " Lorentz Eta: " << Photonxyz[subleadindex].Eta() << " Lorentz Pt: " << Photonp4[subleadindex].Pt() << endl;
          //cout << "My Theta:      " << subleadtheta_simvertex << " My Eta       " << subleadeta_simvertex << " My Pt      " << subleadpt_simvertex << endl;
      
          VLead_simvertex.SetPtEtaPhiE(leadpt_simvertex,leadeta_simvertex,Photonxyz[leadindex].Phi(),Photonp4[leadindex].E());
          VSubLead_simvertex.SetPtEtaPhiE(subleadpt_simvertex,subleadeta_simvertex,Photonxyz[subleadindex].Phi(),Photonp4[subleadindex].E());
          //cout << "LeadPt: " << VLead.Pt() << " LeadEta: " << VLead.Eta() << " LeadPhi: " << VLead.Phi() << " LeadE: " << VLead.E() << endl;
          //cout << "MyLeadPt: " << VLead_simvertex.Pt() << " MyLeadEta: " << VLead_simvertex.Eta() << " MyLeadPhi: " << VLead_simvertex.Phi() << " MyLeadE: " << VLead_simvertex.E() << endl;
          //cout << "SubLeadPt: " << VSubLead.Pt() << " SubLeadEta: " << VSubLead.Eta() << " SubLeadPhi: " << VSubLead.Phi() << " SubLeadE: " << VSubLead.E() << endl;
          //cout << "MySubLeadPt: " << VSubLead_simvertex.Pt() << " MySubLeadEta: " << VSubLead_simvertex.Eta() << " MySubLeadPhi: " << VSubLead_simvertex.Phi() << " MySubLeadE: " << VSubLead_simvertex.E() << endl;
          VSum_simvertex=VLead_simvertex+VSubLead_simvertex;
          InvMass_simvertex=fabs(VSum_simvertex.M());

          // calculate Cos_theta_star 
          beta_b_simvertex  = VSum_simvertex.Beta();
          gamma_b_simvertex = VSum_simvertex.Gamma();
          directionV_simvertex= VSum_simvertex.Vect().Unit();

          // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
          CrossVLead_simvertex=VLead_simvertex.Vect().Cross(directionV);
          DotVLeadValue_simvertex=VLead_simvertex.Vect().Dot(directionV);
          CrossVLeadValue_simvertex=sqrt(CrossVLead_simvertex.x()*CrossVLead_simvertex.x()+CrossVLead_simvertex.y()*CrossVLead_simvertex.y()+CrossVLead_simvertex.z()*CrossVLead_simvertex.z());
          sin_theta_simvertex =  CrossVLeadValue_simvertex/VLead_simvertex.E();
          cos_theta_simvertex =  DotVLeadValue_simvertex/VLead_simvertex.E(); 
          tg_thetas_simvertex = sin_theta_simvertex/(gamma_b_simvertex*(cos_theta_simvertex-beta_b_simvertex));
          cos_thetastar_simvertex= 1.0/sqrt(1.0+tg_thetas_simvertex*tg_thetas_simvertex);
        }
        
        if (Photonp4[leadindex].Pt()>40 && Photonp4[subleadindex].Pt()>30 && InvMass>90 && InvMass<250
            && MarcosCut(Photonp4[leadindex].Pt(), currentTree.pho_ecalsumetconedr04[leadindex], currentTree.pho_hcalsumetconedr04[leadindex], currentTree.pho_trksumpthollowconedr04[leadindex], currentTree.pho_haspixseed[leadindex], currentTree.pho_isEB[leadindex], currentTree.pho_isEE[leadindex], currentTree.pho_sieie[leadindex], currentTree.pho_hoe[leadindex])
            && MarcosCut(Photonp4[subleadindex].Pt(), currentTree.pho_ecalsumetconedr04[subleadindex], currentTree.pho_hcalsumetconedr04[subleadindex], currentTree.pho_trksumpthollowconedr04[subleadindex], currentTree.pho_haspixseed[subleadindex], currentTree.pho_isEB[subleadindex], currentTree.pho_isEE[subleadindex], currentTree.pho_sieie[subleadindex], currentTree.pho_hoe[subleadindex])
            ) {
          hLeadEtMarco->Fill(Photonp4[leadindex].Pt());
          hSubLeadEtMarco->Fill(Photonp4[subleadindex].Pt());
          h_mass_Marco->Fill(InvMass);

          int MarcosCategory = MarcosCutCategory(currentTree.pho_r9[leadindex], currentTree.pho_isEB[leadindex], currentTree.pho_isEE[leadindex], currentTree.pho_r9[subleadindex], currentTree.pho_isEB[subleadindex], currentTree.pho_isEE[subleadindex]);
          hLeadEtMarcoCat[MarcosCategory]->Fill(Photonp4[leadindex].Pt());
          hSubLeadEtMarcoCat[MarcosCategory]->Fill(Photonp4[subleadindex].Pt());
          h_mass_MarcoCat[MarcosCategory]->Fill(InvMass);

        }

        /// di-photon system before event selection
        int HiggsInWhichDetector = 0;
        if (currentTree.pho_isEB[leadindex] && currentTree.pho_isEB[subleadindex]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.pho_isEB[leadindex] && currentTree.pho_isEE[subleadindex]) || (currentTree.pho_isEE[leadindex] && currentTree.pho_isEB[subleadindex])) HiggsInWhichDetector=1; // at least one photon in endcap

        h_mass_2gamma[0][HiggsInWhichDetector]->Fill(InvMass,weight);
        h_pt_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
        h_pz_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
        h_eta_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
        h_phi_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
        h_CosThetaStar[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        // all photon categories together 
        h_mass_2gamma[0][HiggsInWhichDetector]->Fill(InvMass,weight);
        h_pt_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
        h_pz_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
        h_eta_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
        h_phi_2gamma[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
        h_CosThetaStar[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
        
        if (  diPhoCategory==1 ) {
          h_mass_2gamma_2gold[0][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_2gold[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_2gold[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_2gold[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_2gold[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_2gold[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        } else if ( diPhoCategory==2 ) {
          
          h_mass_2gamma_1goodconv[0][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_1goodconv[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_1goodconv[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_1goodconv[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_1goodconv[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_1goodconv[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          h_mass_2gamma_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          h_pt_2gamma_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(VSum_newvertex.Pt(),weight);
          h_pz_2gamma_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(VSum_newvertex.Pz(),weight);
          h_eta_2gamma_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(VSum_newvertex.Eta(),weight);
          h_phi_2gamma_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(VSum_newvertex.Phi(),weight);
          h_CosThetaStar_1goodconv_newvertex[0][HiggsInWhichDetector]->Fill(cos_thetastar_newvertex,weight);

          h_mass_2gamma_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(InvMass_simvertex,weight);
          h_pt_2gamma_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(VSum_simvertex.Pt(),weight);
          h_pz_2gamma_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(VSum_simvertex.Pz(),weight);
          h_eta_2gamma_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(VSum_simvertex.Eta(),weight);
          h_phi_2gamma_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(VSum_simvertex.Phi(),weight);
          h_CosThetaStar_1goodconv_simvertex[0][HiggsInWhichDetector]->Fill(cos_thetastar_simvertex,weight);

          if (ConversionVertex[convindex].Perp()<70) h_mass_2gamma_1goodconv_70[0][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<60) h_mass_2gamma_1goodconv_60[0][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<50) h_mass_2gamma_1goodconv_50[0][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<40) h_mass_2gamma_1goodconv_40[0][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          
        } else if ( diPhoCategory==3 ) {

          h_mass_2gamma_1poorconv[0][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_1poorconv[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_1poorconv[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_1poorconv[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_1poorconv[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_1poorconv[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        } else if ( diPhoCategory==4  ) {

          h_mass_2gamma_2conv[0][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_2conv[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_2conv[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_2conv[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_2conv[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_2conv[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        } else {

          h_mass_2gamma_leftover[0][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_leftover[0][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_leftover[0][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_leftover[0][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_leftover[0][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_leftover[0][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

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
				Photonp4[leadindex].Pt()/ConversionPairMomentum[leadindex].Perp());

        bool convsel2 = convSel(currentTree.pho_conv_ntracks[subleadindex],
				currentTree.pho_conv_validvtx[subleadindex] ,  
				currentTree.pho_conv_chi2_probability[subleadindex], 
				currentTree.pho_conv_dphitrksatvtx[subleadindex], 
				currentTree.pho_conv_paircotthetasep[subleadindex], 
				Photonp4[subleadindex].Pt()/ConversionPairMomentum[subleadindex].Perp());

        HiggsInWhichDetector = 0;
        if (currentTree.pho_isEB[leadindex] && currentTree.pho_isEB[subleadindex]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.pho_isEB[leadindex] && currentTree.pho_isEE[subleadindex]) || (currentTree.pho_isEE[leadindex] && currentTree.pho_isEB[subleadindex])) HiggsInWhichDetector=1; // at least one photon in endcap

        for (unsigned int j=0; j<(unsigned int) currentTree.pho_n; j++) {
          unsigned int detector = 1;
          if (!currentTree.pho_conv_validvtx[j]) continue;
          if (currentTree.pho_isEB[j]) detector=1;
          if (currentTree.pho_isEE[j]) detector=2;
          
          h_phi_conv[0][1]->Fill(Photonxyz[j].Phi());
          h_eta_conv[0][1]->Fill(Photonxyz[j].Eta());
          h_pt_conv[0][1]->Fill(Photonp4[j].Pt());
          h_z_conv[0][1]->Fill(ConversionVertex[j].z());
          h_r_conv[0][1]->Fill(ConversionVertex[j].Perp());

          h_phi_conv[detector][1]->Fill(Photonxyz[j].Phi());
          h_eta_conv[detector][1]->Fill(Photonxyz[j].Eta());
          h_pt_conv[detector][1]->Fill(Photonp4[j].Pt());
          h_z_conv[detector][1]->Fill(ConversionVertex[j].z());
          h_r_conv[detector][1]->Fill(ConversionVertex[j].Perp());
        }

        
        hLeadEt[0][1]->Fill(Photonp4[leadindex].Et(),weight);
        hLeadEta[1]->Fill(Photonxyz[leadindex].Eta(),weight);
        hLeadPhi[1]->Fill(Photonxyz[leadindex].Phi(),weight);
        hLeadR9[0][1]->Fill(currentTree.pho_r9[leadindex],weight);
        hLeadHoE[0][1]->Fill(currentTree.pho_hoe[leadindex],weight);
        hLeadTrkPtSumSolid03[0][1]->Fill(currentTree.pho_trksumptsolidconedr03[leadindex],weight);
        hLeadEcalPtSumSolid03[0][1]->Fill(currentTree.pho_ecalsumetconedr03[leadindex],weight);
        hLeadHcalPtSumSolid03[0][1]->Fill(currentTree.pho_hcalsumetconedr03[leadindex],weight);
        hLeadSigmaIetaIeta[0][1]->Fill(currentTree.pho_sieie[leadindex],weight);
        hLeadZPV_[0][1]->Fill(PrimaryVertex[0].Z(),weight);
        hLeadDzPV_[0][1]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==leadindex && leadPhoCategory==2) hLeadDzPVconv[0][1]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[leadindex]-SimVertex.Z(),weight);

        if (currentTree.pho_conv_ntracks[leadindex]==2 && currentTree.pho_conv_chi2_probability[leadindex]>0.0005 && (bool) currentTree.pho_isEB[leadindex]) {
          conversionVertex.SetXYZ(ConversionVertex[leadindex].X(),ConversionVertex[leadindex].Y(),ConversionVertex[leadindex].Z());
          h2_convVtxRvsZBarrel_[1]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        hSubLeadEt[0][1]->Fill(Photonp4[subleadindex].Et(),weight);
        hSubLeadEta[1]->Fill(Photonxyz[subleadindex].Eta(),weight);
        hSubLeadPhi[1]->Fill(Photonxyz[subleadindex].Phi(),weight);
        hSubLeadR9[0][1]->Fill(currentTree.pho_r9[subleadindex],weight);
        hSubLeadHoE[0][1]->Fill(currentTree.pho_hoe[subleadindex],weight);
        hSubLeadTrkPtSumSolid03[0][1]->Fill(currentTree.pho_trksumptsolidconedr03[subleadindex],weight);
        hSubLeadEcalPtSumSolid03[0][1]->Fill(currentTree.pho_ecalsumetconedr03[subleadindex],weight);
        hSubLeadHcalPtSumSolid03[0][1]->Fill(currentTree.pho_hcalsumetconedr03[subleadindex],weight);
        hSubLeadSigmaIetaIeta[0][1]->Fill(currentTree.pho_sieie[subleadindex],weight);
        hSubLeadZPV_[0][1]->Fill(PrimaryVertex[0].Z(),weight);
        hSubLeadDzPV_[0][1]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==subleadindex && subleadPhoCategory==2) hSubLeadDzPVconv[0][1]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[subleadindex]-SimVertex.Z(),weight);
        
        if (currentTree.pho_conv_ntracks[subleadindex]==2 && currentTree.pho_conv_chi2_probability[subleadindex]>0.0005 && (bool) currentTree.pho_isEB[subleadindex]) {
          conversionVertex.SetXYZ(ConversionVertex[leadindex].X(),ConversionVertex[leadindex].Y(),ConversionVertex[leadindex].Z());
          h2_convVtxRvsZBarrel_[1]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        hNPhotons[1]->Fill(currentTree.pho_n,weight);

        hLeadEt[iLeadDetector][1]->Fill(Photonp4[leadindex].Et(),weight);
        hLeadR9[iLeadDetector][1]->Fill(currentTree.pho_r9[leadindex],weight);
        hLeadHoE[iLeadDetector][1]->Fill(currentTree.pho_hoe[leadindex],weight);
        hLeadTrkPtSumSolid03[iLeadDetector][1]->Fill(currentTree.pho_trksumptsolidconedr03[leadindex],weight);
        hLeadEcalPtSumSolid03[iLeadDetector][1]->Fill(currentTree.pho_ecalsumetconedr03[leadindex],weight);
        hLeadHcalPtSumSolid03[iLeadDetector][1]->Fill(currentTree.pho_hcalsumetconedr03[leadindex],weight);
        hLeadSigmaIetaIeta[iLeadDetector][1]->Fill(currentTree.pho_sieie[leadindex],weight);
        hLeadZPV_[iLeadDetector][1]->Fill(PrimaryVertex[0].Z(),weight);
        hLeadDzPV_[iLeadDetector][1]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==leadindex && leadPhoCategory==2) hLeadDzPVconv[iLeadDetector][1]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[leadindex]-SimVertex.Z(),weight);
        
        hSubLeadEt[iSubleadDetector][1]->Fill(Photonp4[subleadindex].Et(),weight);
        hSubLeadR9[iSubleadDetector][1]->Fill(currentTree.pho_r9[subleadindex],weight);
        hSubLeadHoE[iSubleadDetector][1]->Fill(currentTree.pho_hoe[subleadindex],weight);
        hSubLeadTrkPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.pho_trksumptsolidconedr03[subleadindex],weight);
        hSubLeadEcalPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.pho_ecalsumetconedr03[subleadindex],weight);
        hSubLeadHcalPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.pho_hcalsumetconedr03[subleadindex],weight);
        hSubLeadSigmaIetaIeta[iSubleadDetector][1]->Fill(currentTree.pho_sieie[subleadindex],weight);
        hSubLeadZPV_[iSubleadDetector][1]->Fill(PrimaryVertex[0].Z(),weight);
        hSubLeadDzPV_[iSubleadDetector][1]->Fill(PrimaryVertex[0].Z()-SimVertex.Z(),weight);
        if (convindex==subleadindex && subleadPhoCategory==2) hSubLeadDzPVconv[iSubleadDetector][1]->Fill(currentTree.pho_conv_zofprimvtxfromtrks[subleadindex]-SimVertex.Z(),weight);
        
        if ( convsel1 )
          h2_convVtxRvsZBarrel_[1]->Fill(ConversionVertex[leadindex].Z(),sqrt(ConversionVertex[leadindex].X()*ConversionVertex[leadindex].X()+ConversionVertex[leadindex].Y()*ConversionVertex[leadindex].Y()));
        if ( convsel2 )
          h2_convVtxRvsZBarrel_[1]->Fill(ConversionVertex[leadindex].Z(),sqrt(ConversionVertex[leadindex].X()*ConversionVertex[leadindex].X()+ConversionVertex[leadindex].Y()*ConversionVertex[leadindex].Y()));


        if (  leadPhoCategory==0) 
          h_lead_r9_cat0[iLeadDetector]->Fill (currentTree.pho_r9[leadindex],weight);
        if (  subleadPhoCategory==0) 
          h_sublead_r9_cat0[iSubleadDetector]->Fill (currentTree.pho_r9[subleadindex],weight);


        if (  leadPhoCategory==1) 
          h_lead_r9_cat1[iLeadDetector]->Fill (currentTree.pho_r9[leadindex],weight);
        if (  subleadPhoCategory==1) 
          h_sublead_r9_cat1[iSubleadDetector]->Fill (currentTree.pho_r9[subleadindex],weight);


        if (  leadPhoCategory==2) 
          h_lead_r9_cat2[iLeadDetector]->Fill (currentTree.pho_r9[leadindex],weight);
        if (  subleadPhoCategory==2) 
          h_sublead_r9_cat2[iSubleadDetector]->Fill (currentTree.pho_r9[subleadindex],weight);

        if (  leadPhoCategory==3) 
          h_lead_r9_cat3[iLeadDetector]->Fill (currentTree.pho_r9[leadindex],weight);
        if (  subleadPhoCategory==3) 
          h_sublead_r9_cat3[iSubleadDetector]->Fill (currentTree.pho_r9[subleadindex],weight);


        if (  leadPhoCategory==4) 
          h_lead_r9_cat4[iLeadDetector]->Fill (currentTree.pho_r9[leadindex],weight);
        if (  subleadPhoCategory==4) 
          h_sublead_r9_cat4[iSubleadDetector]->Fill (currentTree.pho_r9[subleadindex],weight);
        
        // calculate invariant mass
        //VLead = TLorentzVector( currentTree.momentumX[leadindex],   currentTree.momentumY[leadindex],  currentTree.momentumZ[leadindex], currentTree.energy[leadindex]);
        //VSubLead = TLorentzVector( currentTree.momentumX[subleadindex],   currentTree.momentumY[subleadindex],  currentTree.momentumZ[subleadindex], currentTree.energy[subleadindex]);
        VLead = TLorentzVector(Photonp4[leadindex]);
        VSubLead = TLorentzVector(Photonp4[subleadindex]);
        VSum=VLead+VSubLead;
        InvMass=fabs(VSum.M());

        // calculate Cos_theta_star 
        beta_b  = VSum.Beta();
        gamma_b = VSum.Gamma();
        directionV= VSum.Vect().Unit();

        // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
        CrossVLead=VLead.Vect().Cross(directionV);
        DotVLeadValue=VLead.Vect().Dot(directionV);
        CrossVLeadValue=sqrt(CrossVLead.x()*CrossVLead.x()+CrossVLead.y()*CrossVLead.y()+CrossVLead.z()*CrossVLead.z());
        sin_theta =  CrossVLeadValue/VLead.E();
        cos_theta =  DotVLeadValue/VLead.E(); 
        tg_thetas = sin_theta/(gamma_b*(cos_theta-beta_b));
        cos_thetastar= 1.0/sqrt(1.0+tg_thetas*tg_thetas);

        if (diPhoCategory==2) {

          leadtheta_newvertex = atan2(Photonxyz[leadindex].Perp(),Photonxyz[leadindex].Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          leadeta_newvertex = -log(tan(leadtheta_newvertex/2));
          leadpt_newvertex = fabs(Photonp4[leadindex].E()*sin(leadtheta_newvertex));

          //cout "Lorentz Theta: " << Photonxyz[leadindex].Theta() << " Lorentz Eta: " << Photonxyz[leadindex].Eta() << " Lorentz Pt: " << Photonp4[leadindex].Pt() << endl;
          //cout "My Theta:      " << leadtheta_newvertex << " My Eta       " << leadeta_newvertex << " My Pt      " << leadpt_newvertex << endl;
        
          subleadtheta_newvertex =  atan2(Photonxyz[subleadindex].Perp(),Photonxyz[subleadindex].Z()-currentTree.pho_conv_zofprimvtxfromtrks[convindex]);
          subleadeta_newvertex = -log(tan(subleadtheta_newvertex/2));
          subleadpt_newvertex = fabs(Photonp4[subleadindex].E()*sin(subleadtheta_newvertex));
          //cout << "Lorentz Theta: " << Photonxyz[subleadindex].Theta() << " Lorentz Eta: " << Photonxyz[subleadindex].Eta() << " Lorentz Pt: " << Photonp4[subleadindex].Pt() << endl;
          //cout << "My Theta:      " << subleadtheta_newvertex << " My Eta       " << subleadeta_newvertex << " My Pt      " << subleadpt_newvertex << endl;
      
          VLead_newvertex.SetPtEtaPhiE(leadpt_newvertex,leadeta_newvertex,Photonxyz[leadindex].Phi(),Photonp4[leadindex].E());
          VSubLead_newvertex.SetPtEtaPhiE(subleadpt_newvertex,subleadeta_newvertex,Photonxyz[subleadindex].Phi(),Photonp4[subleadindex].E());
          //cout << "LeadPt: " << VLead.Pt() << " LeadEta: " << VLead.Eta() << " LeadPhi: " << VLead.Phi() << " LeadE: " << VLead.E() << endl;
          //cout << "MyLeadPt: " << VLead_newvertex.Pt() << " MyLeadEta: " << VLead_newvertex.Eta() << " MyLeadPhi: " << VLead_newvertex.Phi() << " MyLeadE: " << VLead_newvertex.E() << endl;
          //cout << "SubLeadPt: " << VSubLead.Pt() << " SubLeadEta: " << VSubLead.Eta() << " SubLeadPhi: " << VSubLead.Phi() << " SubLeadE: " << VSubLead.E() << endl;
          //cout << "MySubLeadPt: " << VSubLead_newvertex.Pt() << " MySubLeadEta: " << VSubLead_newvertex.Eta() << " MySubLeadPhi: " << VSubLead_newvertex.Phi() << " MySubLeadE: " << VSubLead_newvertex.E() << endl;
          VSum_newvertex=VLead_newvertex+VSubLead_newvertex;
          InvMass_newvertex=fabs(VSum_newvertex.M());

          // calculate Cos_theta_star 
          beta_b_newvertex  = VSum_newvertex.Beta();
          gamma_b_newvertex = VSum_newvertex.Gamma();
          directionV_newvertex= VSum_newvertex.Vect().Unit();

          // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
          CrossVLead_newvertex=VLead_newvertex.Vect().Cross(directionV);
          DotVLeadValue_newvertex=VLead_newvertex.Vect().Dot(directionV);
          CrossVLeadValue_newvertex=sqrt(CrossVLead_newvertex.x()*CrossVLead_newvertex.x()+CrossVLead_newvertex.y()*CrossVLead_newvertex.y()+CrossVLead_newvertex.z()*CrossVLead_newvertex.z());
          sin_theta_newvertex =  CrossVLeadValue_newvertex/VLead_newvertex.E();
          cos_theta_newvertex =  DotVLeadValue_newvertex/VLead_newvertex.E(); 
          tg_thetas_newvertex = sin_theta_newvertex/(gamma_b_newvertex*(cos_theta_newvertex-beta_b_newvertex));
          cos_thetastar_newvertex= 1.0/sqrt(1.0+tg_thetas_newvertex*tg_thetas_newvertex);

          leadtheta_simvertex = atan2(Photonxyz[leadindex].Perp(),Photonxyz[leadindex].Z()-SimVertex.Z());
          leadeta_simvertex = -log(tan(leadtheta_simvertex/2));
          leadpt_simvertex = fabs(Photonp4[leadindex].E()*sin(leadtheta_simvertex));

          //cout "Lorentz Theta: " << Photonxyz[leadindex].Theta() << " Lorentz Eta: " << Photonxyz[leadindex].Eta() << " Lorentz Pt: " << Photonp4[leadindex].Pt() << endl;
          //cout "My Theta:      " << leadtheta_simvertex << " My Eta       " << leadeta_simvertex << " My Pt      " << leadpt_simvertex << endl;
        
          subleadtheta_simvertex =  atan2(Photonxyz[subleadindex].Perp(),Photonxyz[subleadindex].Z()-SimVertex.Z());
          subleadeta_simvertex = -log(tan(subleadtheta_simvertex/2));
          subleadpt_simvertex = fabs(Photonp4[subleadindex].E()*sin(subleadtheta_simvertex));
          //cout << "Lorentz Theta: " << Photonxyz[subleadindex].Theta() << " Lorentz Eta: " << Photonxyz[subleadindex].Eta() << " Lorentz Pt: " << Photonp4[subleadindex].Pt() << endl;
          //cout << "My Theta:      " << subleadtheta_simvertex << " My Eta       " << subleadeta_simvertex << " My Pt      " << subleadpt_simvertex << endl;
      
          VLead_simvertex.SetPtEtaPhiE(leadpt_simvertex,leadeta_simvertex,Photonxyz[leadindex].Phi(),Photonp4[leadindex].E());
          VSubLead_simvertex.SetPtEtaPhiE(subleadpt_simvertex,subleadeta_simvertex,Photonxyz[subleadindex].Phi(),Photonp4[subleadindex].E());
          //cout << "LeadPt: " << VLead.Pt() << " LeadEta: " << VLead.Eta() << " LeadPhi: " << VLead.Phi() << " LeadE: " << VLead.E() << endl;
          //cout << "MyLeadPt: " << VLead_simvertex.Pt() << " MyLeadEta: " << VLead_simvertex.Eta() << " MyLeadPhi: " << VLead_simvertex.Phi() << " MyLeadE: " << VLead_simvertex.E() << endl;
          //cout << "SubLeadPt: " << VSubLead.Pt() << " SubLeadEta: " << VSubLead.Eta() << " SubLeadPhi: " << VSubLead.Phi() << " SubLeadE: " << VSubLead.E() << endl;
          //cout << "MySubLeadPt: " << VSubLead_simvertex.Pt() << " MySubLeadEta: " << VSubLead_simvertex.Eta() << " MySubLeadPhi: " << VSubLead_simvertex.Phi() << " MySubLeadE: " << VSubLead_simvertex.E() << endl;
          VSum_simvertex=VLead_simvertex+VSubLead_simvertex;
          InvMass_simvertex=fabs(VSum_simvertex.M());

          // calculate Cos_theta_star 
          beta_b_simvertex  = VSum_simvertex.Beta();
          gamma_b_simvertex = VSum_simvertex.Gamma();
          directionV_simvertex= VSum_simvertex.Vect().Unit();

          // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
          CrossVLead_simvertex=VLead_simvertex.Vect().Cross(directionV);
          DotVLeadValue_simvertex=VLead_simvertex.Vect().Dot(directionV);
          CrossVLeadValue_simvertex=sqrt(CrossVLead_simvertex.x()*CrossVLead_simvertex.x()+CrossVLead_simvertex.y()*CrossVLead_simvertex.y()+CrossVLead_simvertex.z()*CrossVLead_simvertex.z());
          sin_theta_simvertex =  CrossVLeadValue_simvertex/VLead_simvertex.E();
          cos_theta_simvertex =  DotVLeadValue_simvertex/VLead_simvertex.E(); 
          tg_thetas_simvertex = sin_theta_simvertex/(gamma_b_simvertex*(cos_theta_simvertex-beta_b_simvertex));
          cos_thetastar_simvertex= 1.0/sqrt(1.0+tg_thetas_simvertex*tg_thetas_simvertex);
        }
        
        HiggsInWhichDetector = 0;
        if (currentTree.pho_isEB[leadindex] && currentTree.pho_isEB[subleadindex]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.pho_isEB[leadindex] && currentTree.pho_isEE[subleadindex]) || (currentTree.pho_isEE[leadindex] && currentTree.pho_isEB[subleadindex])) HiggsInWhichDetector=1; // at least one photon in endcap

        // all photon categories together 
        h_mass_2gamma[1][HiggsInWhichDetector]->Fill(InvMass,weight);
        h_pt_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
        h_pz_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
        h_eta_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
        h_phi_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
        h_CosThetaStar[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        if (  diPhoCategory==1 ) {
          h_mass_2gamma_2gold[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_2gold[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        } else if ( diPhoCategory==2 ) {
          h_mass_2gamma_1goodconv[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_1goodconv[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_1goodconv[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_1goodconv[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_1goodconv[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_1goodconv[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          h_mass_2gamma_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          h_pt_2gamma_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(VSum_newvertex.Pt(),weight);
          h_pz_2gamma_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(VSum_newvertex.Pz(),weight);
          h_eta_2gamma_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(VSum_newvertex.Eta(),weight);
          h_phi_2gamma_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(VSum_newvertex.Phi(),weight);
          h_CosThetaStar_1goodconv_newvertex[1][HiggsInWhichDetector]->Fill(cos_thetastar_newvertex,weight);

          h_mass_2gamma_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(InvMass_simvertex,weight);
          h_pt_2gamma_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(VSum_simvertex.Pt(),weight);
          h_pz_2gamma_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(VSum_simvertex.Pz(),weight);
          h_eta_2gamma_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(VSum_simvertex.Eta(),weight);
          h_phi_2gamma_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(VSum_simvertex.Phi(),weight);
          h_CosThetaStar_1goodconv_simvertex[1][HiggsInWhichDetector]->Fill(cos_thetastar_simvertex,weight);

          if (ConversionVertex[convindex].Perp()<70) h_mass_2gamma_1goodconv_70[1][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<60) h_mass_2gamma_1goodconv_60[1][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<50) h_mass_2gamma_1goodconv_50[1][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
          if (ConversionVertex[convindex].Perp()<40) h_mass_2gamma_1goodconv_40[1][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);

        } else if ( diPhoCategory==3 ) {
          h_mass_2gamma_1poorconv[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_1poorconv[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_1poorconv[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_1poorconv[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_1poorconv[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_1poorconv[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);


        } else if ( diPhoCategory==4  ) {
          h_mass_2gamma_2conv[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_2conv[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
        } else {

          h_mass_2gamma_leftover[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar_leftover[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

        }
        //GenMatching
        if (!data) {

          // all photon categories together 
          h_mass_2gamma[2][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_phi_2gamma[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
          h_CosThetaStar[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          if (  diPhoCategory==1 ) {
            h_mass_2gamma_2gold[2][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_2gold[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2gold[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2gold[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_phi_2gamma_2gold[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
            h_CosThetaStar_2gold[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          } else if ( diPhoCategory==2 ) {
            h_mass_2gamma_1goodconv[2][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_1goodconv[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_1goodconv[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_1goodconv[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_phi_2gamma_1goodconv[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
            h_CosThetaStar_1goodconv[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

            h_mass_2gamma_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
            h_pt_2gamma_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(VSum_newvertex.Pt(),weight);
            h_pz_2gamma_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(VSum_newvertex.Pz(),weight);
            h_eta_2gamma_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(VSum_newvertex.Eta(),weight);
            h_phi_2gamma_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(VSum_newvertex.Phi(),weight);
            h_CosThetaStar_1goodconv_newvertex[2][HiggsInWhichDetector]->Fill(cos_thetastar_newvertex,weight);

            h_mass_2gamma_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(InvMass_simvertex,weight);
            h_pt_2gamma_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(VSum_simvertex.Pt(),weight);
            h_pz_2gamma_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(VSum_simvertex.Pz(),weight);
            h_eta_2gamma_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(VSum_simvertex.Eta(),weight);
            h_phi_2gamma_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(VSum_simvertex.Phi(),weight);
            h_CosThetaStar_1goodconv_simvertex[2][HiggsInWhichDetector]->Fill(cos_thetastar_simvertex,weight);

            if (ConversionVertex[convindex].Perp()<70) h_mass_2gamma_1goodconv_70[2][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
            if (ConversionVertex[convindex].Perp()<60) h_mass_2gamma_1goodconv_60[2][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
            if (ConversionVertex[convindex].Perp()<50) h_mass_2gamma_1goodconv_50[2][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);
            if (ConversionVertex[convindex].Perp()<40) h_mass_2gamma_1goodconv_40[2][HiggsInWhichDetector]->Fill(InvMass_newvertex,weight);

          } else if ( diPhoCategory==3 ) {
            h_mass_2gamma_1poorconv[2][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_1poorconv[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_1poorconv[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_1poorconv[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_phi_2gamma_1poorconv[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
            h_CosThetaStar_1poorconv[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);


          } else if ( diPhoCategory==4  ) {
            h_mass_2gamma_2conv[2][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_2conv[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2conv[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2conv[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_phi_2gamma_2conv[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
            h_CosThetaStar_2conv[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
          } else {

            h_mass_2gamma_leftover[2][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_leftover[2][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_leftover[2][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_leftover[2][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_phi_2gamma_leftover[2][HiggsInWhichDetector]->Fill(VSum.Phi(),weight);
            h_CosThetaStar_leftover[2][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          }
          
        }
        
      }
      currentFile->Close();
      delete currentFile;

    }
    
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


int DetectorPosition(sdaReader *currentTree, unsigned int index) {

  int ReturnValue=0;
  if (currentTree->pho_isEB[index]) ReturnValue=1;
  if (currentTree->pho_isEE[index]) ReturnValue=2;
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

string MakeFileName(string filename, bool unweighted, bool dataweight) {

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

  return outfilename;
  
}

void MakeFilesAndWeights(TString &inputstring, vector<pair<string, float> > &inputvector, vector<pair<string, int> > &inputfilelist, bool &isData) {

  float BranchingFraction = 0;

  if (inputstring.Contains("Data")) {
    inputfilelist.push_back(pair<string,int> ("Data.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc4/b/tkolberg/MPAntuples/Oct1_preselOriginalRefitMomentum/data.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc4/b/tkolberg/MPAntuples/Nov6_V00-00-13/Nov6_V00-00-13.root",1));
    isData=true;
  }
  if (inputstring.Contains("Yousidata")) {
    inputfilelist.push_back(pair<string,int> ("YousiData.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/YousiData/MPA_Run2010B_Nov5_12831nb.root",1));
    isData=true;
  }
  if (inputstring.Contains("90GeV") || inputstring.Contains("All")) {
    BranchingFraction = 0.000726;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis90GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon90.root",34.145*BranchingFraction/98996));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF90.root",1.7801*BranchingFraction/108813));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ90.root",3.1880*BranchingFraction/88000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
  }
  if (inputstring.Contains("95GeV") || inputstring.Contains("All")) {
    BranchingFraction = 0.00108;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis95GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon95.root",31.886*BranchingFraction/83986));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF95.root",1.6956*BranchingFraction/109579));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ95.root",2.9226*BranchingFraction/110000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
  }
  if (inputstring.Contains("100GeV") || inputstring.Contains("All")) {
    BranchingFraction = 0.00140;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis100GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF100.root",1.5929*BranchingFraction/109826));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ100.root",1.9366*BranchingFraction/110000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
    cout << "Warning no Gluon Fusion Samples!!!!!" << endl;
  }
  if (inputstring.Contains("105GeV") || inputstring.Contains("All")) {
    BranchingFraction = 0.001755;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis105GeV.root",2));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF105.root",1.5138*BranchingFraction/109835));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ105.root",1.6674*BranchingFraction/110000));
    cout << "Warning no Gluon Fusion Samples!!!!!" << endl;
  }
  if (inputstring.Contains("110GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.001939;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis110GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon110.root",20.493*BranchingFraction/109994));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF110.root",1.4405*BranchingFraction/105974));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ110.root",1.4421*BranchingFraction/110000));
  }
  if (inputstring.Contains("115GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002101;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis115GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon115.root",18.735*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF115.root",1.3712*BranchingFraction/109834));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ115.root",1.2524*BranchingFraction/110000));
  }
  if (inputstring.Contains("120GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002219;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis120GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon120.root",17.173*BranchingFraction/106151));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF120.root",1.3062*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ120.root",1.0921*BranchingFraction/110000));
  }
  if (inputstring.Contains("130GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.002240;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis130GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon130.root",14.579*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF130.root",1.1866*BranchingFraction/109848));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ130.root",0.8395*BranchingFraction/110000));
  }
  if (inputstring.Contains("140GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.001929;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis140GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon140.root",12.525*BranchingFraction/109991));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF140.root",1.0811*BranchingFraction/109842));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ140.root",0.6539*BranchingFraction/110000));
  }
  if (inputstring.Contains("150GeV") || inputstring.Contains("Signal") || inputstring.Contains("All")) {
    BranchingFraction = 0.001363;
    inputfilelist.push_back(pair<string,int> ("HiggsAnalysis150GeV.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon150.root",10.863*BranchingFraction/43500));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF150.root",0.9868*BranchingFraction/50000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ150.root",0.5155*BranchingFraction/48000));
  }
  if (inputstring.Contains("PhotonPlusJet") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("PhotonPlusJet.root",11));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet0to15.root",84200000.0/1057100));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet15to30.root",171700.0/1025840));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet30to50.root",16690.0/1025480));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet50to80.root",2722.0/1024608));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet80to120.root",447.2/1048215));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet120to170.root",84.17/1023361));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet170to300.root",22.64/1089000));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet300to470.root",1.493/1076926));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet470to800.root",0.1323/1093499));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet800to1400.root",0.003481/1092742));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet1400to1800.root",0.00001270/1097060));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet1800toInf.root",0.0000002936/1091360));
  }
  if (inputstring.Contains("EMEnriched") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("EMEnriched.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt20to30.root",236000000/(37169939/0.0104)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt30to80.root",59480000/(71845473/0.065)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt80to170.root",900000/(8073559/0.155)));
  }
  if (inputstring.Contains("Doubleemenriched") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("DoubleEMEnriched.root",1));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt10to20.root",20750000000/(31536145/0.0563)));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt20.root",293300000/(10912061/0.239)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt40.root",18700000/(21229315/0.00216)));
  }
  if (inputstring.Contains("Reweighteddoubleemenriched") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("ReweightedDoubleEMEnriched.root",1));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt10to20.root",1.15*20750000000/(31536145/0.0563)));
    //inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt20.root",1.15*293300000/(10912061/0.239)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt40.root",1.15*18700000/(21229315/0.00216)));
  }
  if (inputstring.Contains("QCDBCtoE") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("QCDBCtoE.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt20to30.root",236000000/(2243439/0.00056)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt30to80.root",59480000/(1995502/0.00230)));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt80to170.root",900000/(1043390/0.0104)));
  }
  if (inputstring.Contains("Born") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Born.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt10to25.root",236.4/523270));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt25to250.root",22.37/536230));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt250toInf.root",0.008072/541900));
  }
  if (inputstring.Contains("Box") || inputstring.Contains("Background") || inputstring.Contains("All")) {
    inputfilelist.push_back(pair<string,int> ("Box.root",3));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt10to25.root",358.2/792710));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt25to250.root",12.37/768815));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt250toInf.root",0.000208/790685));
  }
  if (inputstring.Contains("MPATest")) {
    inputfilelist.push_back(pair<string,int> ("MPATest.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon130.root",1));
  }
  if (inputstring.Contains("SDATest")) {
    inputfilelist.push_back(pair<string,int> ("SDATest.root",1));
    inputvector.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/hgg_reduced.root",1));
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

