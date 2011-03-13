///////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  1 15:30:00 2010 by ROOT version 5.22/00d
// from TTree Analysis/Analysis
// found on file: /data/ndpc1/b/tkolberg/MPAntuples/PhotonJet_Pt170_Spring10-START3X_V26_S09-v1_GEN-SIM-RECO/MultiPhotonAnalyzer.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <iostream>

using namespace std;

class sdaReader {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  TClonesArray    *sc_p4;
  Int_t           pho_n;
  Int_t           gp_n;
  Int_t           pho_isEB[100];   //[pho_n]
  Int_t           pho_isEE[100];   //[pho_n]
  Int_t           pho_isEBGap[100];   //[pho_n]
  Int_t           pho_isEEGap[100];   //[pho_n]
  Int_t           pho_isEBEEGap[100];   //[pho_n]
  Float_t         pho_see[100];   //[pho_n]
  Float_t         pho_sieie[100];   //[pho_n]
  Float_t         pho_e1x5[100];   //[pho_n]
  Float_t         pho_e2x5[100];   //[pho_n]
  Float_t         pho_e3x3[100];   //[pho_n]
  Float_t         pho_e5x5[100];   //[pho_n]
  Float_t         pho_emaxxtal[100];   //[pho_n]
  Float_t         pho_hoe[100];   //[pho_n]
  Float_t         pho_h1oe[100];   //[pho_n]
  Float_t         pho_h2oe[100];   //[pho_n]
  Float_t         pho_r1x5[100];   //[pho_n]
  Float_t         pho_r2x5[100];   //[pho_n]
  Float_t         pho_r9[100];   //[pho_n]
  Float_t         pho_ecalsumetconedr04[100];   //[pho_n]
  Float_t         pho_hcalsumetconedr04[100];   //[pho_n]
  Float_t         pho_hcal1sumetconedr04[100];   //[pho_n]
  Float_t         pho_hcal2sumetconedr04[100];   //[pho_n]
  Float_t         pho_trksumptsolidconedr04[100];   //[pho_n]
  Float_t         pho_trksumpthollowconedr04[100];   //[pho_n]
  Float_t         pho_ntrksolidconedr04[100];   //[pho_n]
  Float_t         pho_ntrkhollowconedr04[100];   //[pho_n]
  Float_t         pho_ecalsumetconedr03[100];   //[pho_n]
  Float_t         pho_hcalsumetconedr03[100];   //[pho_n]
  Float_t         pho_hcal1sumetconedr03[100];   //[pho_n]
  Float_t         pho_hcal2sumetconedr03[100];   //[pho_n]
  Float_t         pho_trksumptsolidconedr03[100];   //[pho_n]
  Float_t         pho_trksumpthollowconedr03[100];   //[pho_n]
  Float_t         pho_ntrksolidconedr03[100];   //[pho_n]
  Float_t         pho_ntrkhollowconedr03[100];   //[pho_n]
  TClonesArray    *pho_p4;
  TClonesArray    *pho_calopos;
  Int_t           pho_barrel[100];   //[pho_n]
  Int_t           pho_scind[100];   //[pho_n]
  Int_t           pho_haspixseed[100];   //[pho_n]
  Int_t           pho_hasconvtks[100];   //[pho_n]
  Int_t           pho_nconv[100];   //[pho_n]
  Int_t           pho_conv_ntracks[100];   //[pho_n]
  Float_t         pho_conv_pairinvmass[100];   //[pho_n]
  Float_t         pho_conv_paircotthetasep[100];   //[pho_n]
  Float_t         pho_conv_eoverp[100];   //[pho_n]
  Float_t         pho_conv_zofprimvtxfromtrks[100];   //[pho_n]
  Float_t         pho_conv_distofminapproach[100];   //[pho_n]
  Float_t         pho_conv_dphitrksatvtx[100];   //[pho_n]
  Float_t         pho_conv_dphitrksatecal[100];   //[pho_n]
  Float_t         pho_conv_detatrksatecal[100];   //[pho_n]
  Float_t         pho_conv_tk1_d0[100];   //[pho_n]
  Float_t         pho_conv_tk1_pout[100];   //[pho_n]
  Float_t         pho_conv_tk1_pin[100];   //[pho_n]
  Float_t         pho_conv_tk2_d0[100];   //[pho_n]
  Float_t         pho_conv_tk2_pout[100];   //[pho_n]
  Float_t         pho_conv_tk2_pin[100];   //[pho_n]
  Float_t         pho_conv_tk1_dz[100];   //[pho_n]
  Float_t         pho_conv_tk2_dz[100];   //[pho_n]
  Float_t         pho_conv_tk1_dzerr[100];   //[pho_n]
  Float_t         pho_conv_tk2_dzerr[100];   //[pho_n]
  Int_t           pho_conv_tk1_nh[100];   //[pho_n]
  Int_t           pho_conv_tk2_nh[100];   //[pho_n]
  Float_t         pho_conv_chi2[100];   //[pho_n]
  Float_t         pho_conv_chi2_probability[100];   //[pho_n]
  Int_t           pho_conv_ch1ch2[100];   //[pho_n]
  Int_t           pho_conv_validvtx[100];   //[pho_n]
  Int_t           pho_conv_MVALikelihood[100];   //[pho_n]
  Int_t           gp_pdgid[2000];   //[gp_n]
  Int_t           gp_status[2000];   //[gp_n]
  TClonesArray    *pho_conv_vtx;
  TClonesArray    *pho_conv_pair_momentum;
  TClonesArray    *pho_conv_refitted_momentum;
  TClonesArray    *pho_conv_vertexcorrected_p4;
  TClonesArray    *vtx_std_xyz;
  TClonesArray    *simvtx;
  
  // List of branches
  TBranch        *b_sc_p4;
  TBranch        *b_pho_n;
  TBranch        *b_gp_n;
  TBranch        *b_pho_isEB;   //!
  TBranch        *b_pho_isEE;   //!
  TBranch        *b_pho_isEBGap;   //!
  TBranch        *b_pho_isEEGap;   //!
  TBranch        *b_pho_isEBEEGap;   //!
  TBranch        *b_pho_see;   //!
  TBranch        *b_pho_sieie;   //!
  TBranch        *b_pho_e1x5;   //!
  TBranch        *b_pho_e2x5;   //!
  TBranch        *b_pho_e3x3;   //!
  TBranch        *b_pho_e5x5;   //!
  TBranch        *b_pho_emaxxtal;   //!
  TBranch        *b_pho_hoe;   //!
  TBranch        *b_pho_h1oe;   //!
  TBranch        *b_pho_h2oe;   //!
  TBranch        *b_pho_r1x5;   //!
  TBranch        *b_pho_r2x5;   //!
  TBranch        *b_pho_r9;   //!
  TBranch        *b_pho_ecalsumetconedr04;   //!
  TBranch        *b_pho_hcalsumetconedr04;   //!
  TBranch        *b_pho_hcal1sumetconedr04;   //!
  TBranch        *b_pho_hcal2sumetconedr04;   //!
  TBranch        *b_pho_trksumptsolidconedr04;   //!
  TBranch        *b_pho_trksumpthollowconedr04;   //!
  TBranch        *b_pho_ntrksolidconedr04;   //!
  TBranch        *b_pho_ntrkhollowconedr04;   //!
  TBranch        *b_pho_ecalsumetconedr03;   //!
  TBranch        *b_pho_hcalsumetconedr03;   //!
  TBranch        *b_pho_hcal1sumetconedr03;   //!
  TBranch        *b_pho_hcal2sumetconedr03;   //!
  TBranch        *b_pho_trksumptsolidconedr03;   //!
  TBranch        *b_pho_trksumpthollowconedr03;   //!
  TBranch        *b_pho_ntrksolidconedr03;   //!
  TBranch        *b_pho_ntrkhollowconedr03;   //!
  TBranch        *b_pho_p4;   //!
  TBranch        *b_pho_calopos;   //!
  TBranch        *b_pho_barrel;   //!
  TBranch        *b_pho_scind;   //!
  TBranch        *b_pho_haspixseed;   //!
  TBranch        *b_pho_hasconvtks;   //!
  TBranch        *b_pho_nconv;   //!
  TBranch        *b_pho_conv_ntracks;   //!
  TBranch        *b_pho_conv_pairinvmass;   //!
  TBranch        *b_pho_conv_paircotthetasep;   //!
  TBranch        *b_pho_conv_eoverp;   //!
  TBranch        *b_pho_conv_zofprimvtxfromtrks;   //!
  TBranch        *b_pho_conv_distofminapproach;   //!
  TBranch        *b_pho_conv_dphitrksatvtx;   //!
  TBranch        *b_pho_conv_dphitrksatecal;   //!
  TBranch        *b_pho_conv_detatrksatecal;   //!
  TBranch        *b_pho_conv_tk1_d0;   //!
  TBranch        *b_pho_conv_tk1_pout;   //!
  TBranch        *b_pho_conv_tk1_pin;   //!
  TBranch        *b_pho_conv_tk2_d0;   //!
  TBranch        *b_pho_conv_tk2_pout;   //!
  TBranch        *b_pho_conv_tk2_pin;   //!
  TBranch        *b_pho_conv_tk1_dz;   //!
  TBranch        *b_pho_conv_tk2_dz;   //!
  TBranch        *b_pho_conv_tk1_dzerr;   //!
  TBranch        *b_pho_conv_tk2_dzerr;   //!
  TBranch        *b_pho_conv_tk1_nh;   //!
  TBranch        *b_pho_conv_tk2_nh;   //!
  TBranch        *b_pho_conv_chi2;   //!
  TBranch        *b_pho_conv_chi2_probability;   //!
  TBranch        *b_pho_conv_ch1ch2;   //!
  TBranch        *b_pho_conv_validvtx;   //!
  TBranch        *b_pho_conv_MVALikelihood;   //!
  TBranch        *b_gp_pdgid;   //!
  TBranch        *b_gp_status;   //!
  TBranch        *b_pho_conv_vtx;   //!
  TBranch        *b_pho_conv_pair_momentum;   //!
  TBranch        *b_pho_conv_refitted_momentum;   //!
  TBranch        *b_pho_conv_vertexcorrected_p4;   //!
  TBranch        *b_vtx_std_xyz;   //!
  TBranch        *b_simvtx;   //!

  sdaReader(TFile *currentFile=0);
  virtual ~sdaReader();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init(TFile *currentFile);
  virtual Bool_t   Notify();
};

sdaReader::sdaReader(TFile *currentFile)
{
  if (currentFile == 0) cout << "Warning!!!!! Analysis tree is emtpy!!!!!" << endl;
  Init(currentFile);
}

sdaReader::~sdaReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t sdaReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void sdaReader::Init(TFile *currentFile)
{
  currentFile->cd();
  TTree * tree = (TTree *) currentFile->Get("event");
  
  // Set object pointer
  sc_p4 = 0;
  pho_p4 = 0;
  pho_calopos = 0;
  pho_conv_vtx = 0;
  pho_conv_pair_momentum = 0;
  pho_conv_refitted_momentum = 0;
  pho_conv_vertexcorrected_p4 = 0;
  vtx_std_xyz = 0;
  simvtx = 0;
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  if (fChain->FindLeaf("sc_p4")!=NULL) fChain->SetBranchAddress("sc_p4", &sc_p4, &b_sc_p4);
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  if (fChain->FindLeaf("gp_n")!=NULL) fChain->SetBranchAddress("gp_n", &gp_n, &b_gp_n);
  fChain->SetBranchAddress("pho_isEB", pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_isEBGap", pho_isEBGap, &b_pho_isEBGap);
  fChain->SetBranchAddress("pho_isEEGap", pho_isEEGap, &b_pho_isEEGap);
  fChain->SetBranchAddress("pho_isEBEEGap", pho_isEBEEGap, &b_pho_isEBEEGap);
  fChain->SetBranchAddress("pho_see", pho_see, &b_pho_see);
  fChain->SetBranchAddress("pho_sieie", pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_e1x5", pho_e1x5, &b_pho_e1x5);
  fChain->SetBranchAddress("pho_e2x5", pho_e2x5, &b_pho_e2x5);
  fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
  fChain->SetBranchAddress("pho_e5x5", pho_e5x5, &b_pho_e5x5);
  fChain->SetBranchAddress("pho_emaxxtal", pho_emaxxtal, &b_pho_emaxxtal);
  fChain->SetBranchAddress("pho_hoe", pho_hoe, &b_pho_hoe);
  fChain->SetBranchAddress("pho_h1oe", pho_h1oe, &b_pho_h1oe);
  fChain->SetBranchAddress("pho_h2oe", pho_h2oe, &b_pho_h2oe);
  fChain->SetBranchAddress("pho_r1x5", pho_r1x5, &b_pho_r1x5);
  fChain->SetBranchAddress("pho_r2x5", pho_r2x5, &b_pho_r2x5);
  fChain->SetBranchAddress("pho_r9", pho_r9, &b_pho_r9);
  fChain->SetBranchAddress("pho_ecalsumetconedr04", pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
  fChain->SetBranchAddress("pho_hcalsumetconedr04", pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
  fChain->SetBranchAddress("pho_hcal1sumetconedr04", pho_hcal1sumetconedr04, &b_pho_hcal1sumetconedr04);
  fChain->SetBranchAddress("pho_hcal2sumetconedr04", pho_hcal2sumetconedr04, &b_pho_hcal2sumetconedr04);
  fChain->SetBranchAddress("pho_trksumptsolidconedr04", pho_trksumptsolidconedr04, &b_pho_trksumptsolidconedr04);
  fChain->SetBranchAddress("pho_trksumpthollowconedr04", pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  fChain->SetBranchAddress("pho_ntrksolidconedr04", pho_ntrksolidconedr04, &b_pho_ntrksolidconedr04);
  fChain->SetBranchAddress("pho_ntrkhollowconedr04", pho_ntrkhollowconedr04, &b_pho_ntrkhollowconedr04);
  fChain->SetBranchAddress("pho_ecalsumetconedr03", pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
  fChain->SetBranchAddress("pho_hcalsumetconedr03", pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
  fChain->SetBranchAddress("pho_hcal1sumetconedr03", pho_hcal1sumetconedr03, &b_pho_hcal1sumetconedr03);
  fChain->SetBranchAddress("pho_hcal2sumetconedr03", pho_hcal2sumetconedr03, &b_pho_hcal2sumetconedr03);
  fChain->SetBranchAddress("pho_trksumptsolidconedr03", pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
  fChain->SetBranchAddress("pho_trksumpthollowconedr03", pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
  fChain->SetBranchAddress("pho_ntrksolidconedr03", pho_ntrksolidconedr03, &b_pho_ntrksolidconedr03);
  fChain->SetBranchAddress("pho_ntrkhollowconedr03", pho_ntrkhollowconedr03, &b_pho_ntrkhollowconedr03);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
  fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
  fChain->SetBranchAddress("pho_barrel", pho_barrel, &b_pho_barrel);
  fChain->SetBranchAddress("pho_scind", pho_scind, &b_pho_scind);
  fChain->SetBranchAddress("pho_haspixseed", pho_haspixseed, &b_pho_haspixseed);
  fChain->SetBranchAddress("pho_hasconvtks", pho_hasconvtks, &b_pho_hasconvtks);
  fChain->SetBranchAddress("pho_nconv", pho_nconv, &b_pho_nconv);
  fChain->SetBranchAddress("pho_conv_ntracks", pho_conv_ntracks, &b_pho_conv_ntracks);
  fChain->SetBranchAddress("pho_conv_pairinvmass", pho_conv_pairinvmass, &b_pho_conv_pairinvmass);
  fChain->SetBranchAddress("pho_conv_paircotthetasep", pho_conv_paircotthetasep, &b_pho_conv_paircotthetasep);
  fChain->SetBranchAddress("pho_conv_eoverp", pho_conv_eoverp, &b_pho_conv_eoverp);
  fChain->SetBranchAddress("pho_conv_zofprimvtxfromtrks", pho_conv_zofprimvtxfromtrks, &b_pho_conv_zofprimvtxfromtrks);
  fChain->SetBranchAddress("pho_conv_distofminapproach", pho_conv_distofminapproach, &b_pho_conv_distofminapproach);
  fChain->SetBranchAddress("pho_conv_dphitrksatvtx", pho_conv_dphitrksatvtx, &b_pho_conv_dphitrksatvtx);
  fChain->SetBranchAddress("pho_conv_dphitrksatecal", pho_conv_dphitrksatecal, &b_pho_conv_dphitrksatecal);
  fChain->SetBranchAddress("pho_conv_detatrksatecal", pho_conv_detatrksatecal, &b_pho_conv_detatrksatecal);
  fChain->SetBranchAddress("pho_conv_tk1_d0", pho_conv_tk1_d0, &b_pho_conv_tk1_d0);
  fChain->SetBranchAddress("pho_conv_tk1_pout", pho_conv_tk1_pout, &b_pho_conv_tk1_pout);
  fChain->SetBranchAddress("pho_conv_tk1_pin", pho_conv_tk1_pin, &b_pho_conv_tk1_pin);
  fChain->SetBranchAddress("pho_conv_tk2_d0", pho_conv_tk2_d0, &b_pho_conv_tk2_d0);
  fChain->SetBranchAddress("pho_conv_tk2_pout", pho_conv_tk2_pout, &b_pho_conv_tk2_pout);
  fChain->SetBranchAddress("pho_conv_tk2_pin", pho_conv_tk2_pin, &b_pho_conv_tk2_pin);
  fChain->SetBranchAddress("pho_conv_tk1_dz", pho_conv_tk1_dz, &b_pho_conv_tk1_dz);
  fChain->SetBranchAddress("pho_conv_tk2_dz", pho_conv_tk2_dz, &b_pho_conv_tk2_dz);
  fChain->SetBranchAddress("pho_conv_tk1_dzerr", pho_conv_tk1_dzerr, &b_pho_conv_tk1_dzerr);
  fChain->SetBranchAddress("pho_conv_tk2_dzerr", pho_conv_tk2_dzerr, &b_pho_conv_tk2_dzerr);
  fChain->SetBranchAddress("pho_conv_tk1_nh", pho_conv_tk1_nh, &b_pho_conv_tk1_nh);
  fChain->SetBranchAddress("pho_conv_tk2_nh", pho_conv_tk2_nh, &b_pho_conv_tk2_nh);
  fChain->SetBranchAddress("pho_conv_chi2", pho_conv_chi2, &b_pho_conv_chi2);
  fChain->SetBranchAddress("pho_conv_chi2_probability", pho_conv_chi2_probability, &b_pho_conv_chi2_probability);
  fChain->SetBranchAddress("pho_conv_ch1ch2", pho_conv_ch1ch2, &b_pho_conv_ch1ch2);
  fChain->SetBranchAddress("pho_conv_validvtx", pho_conv_validvtx, &b_pho_conv_validvtx);
  fChain->SetBranchAddress("pho_conv_MVALikelihood", pho_conv_MVALikelihood, &b_pho_conv_MVALikelihood);
  if (fChain->FindBranch("gp_pdgid")!=NULL) fChain->SetBranchAddress("gp_pdgid", gp_pdgid, &b_gp_pdgid);
  if (fChain->FindBranch("gp_status")!=NULL) fChain->SetBranchAddress("gp_status", gp_status, &b_gp_status);
  fChain->SetBranchAddress("pho_conv_vtx", &pho_conv_vtx, &b_pho_conv_vtx);
  fChain->SetBranchAddress("pho_conv_pair_momentum", &pho_conv_pair_momentum, &b_pho_conv_pair_momentum);
  fChain->SetBranchAddress("pho_conv_refitted_momentum", &pho_conv_refitted_momentum, &b_pho_conv_refitted_momentum);
  fChain->SetBranchAddress("pho_conv_vertexcorrected_p4", &pho_conv_vertexcorrected_p4, &b_pho_conv_vertexcorrected_p4);
  fChain->SetBranchAddress("vtx_std_xyz", &vtx_std_xyz, &b_vtx_std_xyz);
  if (fChain->FindLeaf("simvtx")!=NULL) fChain->SetBranchAddress("simvtx", &simvtx, &b_simvtx);
  if (fChain->FindLeaf("gp_vtx")!=NULL) fChain->SetBranchAddress("gp_vtx", &simvtx, &b_simvtx);
  
  Notify();
}

Bool_t sdaReader::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  return kTRUE;
}
