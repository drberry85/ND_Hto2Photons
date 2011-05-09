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
  Int_t           event;
  Int_t           lumis;
  Int_t           run;
  Int_t           bx;
  TClonesArray    *sc_p4;
  Int_t           pho_n;
  Float_t         pho_feta[100][100];
  Float_t         pho_crackcorr[100];
  Float_t         pho_localcorr[100];
  Int_t           pho_isEB[100];
  Int_t           pho_isEE[100];
  Int_t           pho_isEBGap[100];
  Int_t           pho_isEEGap[100];
  Int_t           pho_isEBEEGap[100];
  Float_t         pho_see[100];
  Float_t         pho_sieie[100];
  Float_t         pho_sipip[100];
  Float_t         pho_sieip[100];
  Float_t         pho_e1x5[100];
  Float_t         pho_e2x5[100];
  Float_t         pho_e3x3[100];
  Float_t         pho_e5x5[100];
  Float_t         pho_emaxxtal[100];
  Float_t         pho_hoe[100];
  Float_t         pho_h1oe[100];
  Float_t         pho_h2oe[100];
  Float_t         pho_r1x5[100];
  Float_t         pho_r2x5[100];
  Float_t         pho_r9[100];
  Float_t         pho_zernike20[100];
  Float_t         pho_zernike42[100];
  Float_t         pho_e2nd[100];
  Float_t         pho_e2x5right[100];
  Float_t         pho_e2x5left[100];
  Float_t         pho_e2x5Top[100];
  Float_t         pho_e2x5bottom[100];
  Float_t         pho_eright[100];
  Float_t         pho_eleft[100];
  Float_t         pho_etop[100];
  Float_t         pho_ebottom[100];
  Float_t         pho_e2overe9[100];
  Float_t         pho_seed_time[100];
  Float_t         pho_seed_outoftimechi2[100];
  Float_t         pho_seed_chi2[100];
  Float_t         pho_seed_recoflag[100];
  Float_t         pho_seed_severity[100];
  Float_t         pho_ecalsumetconedr04[100];
  Float_t         pho_hcalsumetconedr04[100];
  Float_t         pho_hcal1sumetconedr04[100];
  Float_t         pho_hcal2sumetconedr04[100];
  Float_t         pho_trksumptsolidconedr04[100];
  Float_t         pho_trksumpthollowconedr04[100];
  Float_t         pho_ntrksolidconedr04[100];
  Float_t         pho_ntrkhollowconedr04[100];
  Float_t         pho_ecalsumetconedr03[100];
  Float_t         pho_hcalsumetconedr03[100];
  Float_t         pho_hcal1sumetconedr03[100];
  Float_t         pho_hcal2sumetconedr03[100];
  Float_t         pho_trksumptsolidconedr03[100];
  Float_t         pho_trksumpthollowconedr03[100];
  Float_t         pho_ntrksolidconedr03[100];
  Float_t         pho_ntrkhollowconedr03[100];
  Int_t           pho_barrel[100];
  Int_t           pho_scind[100];
  Int_t           pho_haspixseed[100];
  Int_t           pho_hasconvtks[100];
  Int_t           pho_nconv[100];
  Int_t           pho_conv_ntracks[100];
  Float_t         pho_conv_pairinvmass[100];
  Float_t         pho_conv_paircotthetasep[100];
  Float_t         pho_conv_eoverp[100];
  Float_t         pho_conv_zofprimvtxfromtrks[100];
  Float_t         pho_conv_distofminapproach[100];
  Float_t         pho_conv_dphitrksatvtx[100];
  Float_t         pho_conv_dphitrksatecal[100];
  Float_t         pho_conv_detatrksatecal[100];
  Float_t         pho_conv_tk1_d0[100];
  Float_t         pho_conv_tk1_pout[100];
  Float_t         pho_conv_tk1_pin[100];
  Float_t         pho_conv_tk2_d0[100];
  Float_t         pho_conv_tk2_pout[100];
  Float_t         pho_conv_tk2_pin[100];
  Float_t         pho_conv_tk1_dz[100];
  Float_t         pho_conv_tk1_dzerr[100];
  Int_t           pho_conv_tk1_nh[100];
  Float_t         pho_conv_tk2_dz[100];
  Float_t         pho_conv_tk2_dzerr[100];
  Int_t           pho_conv_tk2_nh[100];
  Int_t           pho_conv_ch1ch2[100];
  Float_t         pho_conv_chi2[100];
  Float_t         pho_conv_chi2_probability[100];
  Int_t           pho_conv_validvtx[100];
  Int_t           pho_conv_MVALikelihood[100];
  TClonesArray    *pho_p4;
  TClonesArray    *pho_calopos;
  TClonesArray    *pho_conv_vtx;
  TClonesArray    *pho_conv_pair_momentum;
  TClonesArray    *pho_conv_refitted_momentum;
  TClonesArray    *pho_conv_vertexcorrected_p4;
  Int_t           conv_n;
  TClonesArray    *conv_p4;
  Int_t           conv_ntracks[1000];
  Float_t         conv_pairinvmass[1000];
  Float_t         conv_paircotthetasep[1000];
  Float_t         conv_eoverp[1000];
  Float_t         conv_distofminapproach[1000];
  Float_t         conv_dphitrksatvtx[1000];
  Float_t         conv_dphitrksatecal[1000];
  Float_t         conv_detatrksatecal[1000];
  Float_t         conv_dxy[1000];
  Float_t         conv_dz[1000];
  Float_t         conv_lxy[1000];
  Float_t         conv_lz[1000];
  Float_t         conv_zofprimvtxfromtrks[1000];
  vector<vector<unsigned short> > *conv_nHitsBeforeVtx;
  Int_t           conv_nSharedHits[1000];
  Int_t           conv_validvtx[1000];
  Int_t           conv_MVALikelihood[1000];
  Float_t         conv_chi2[1000];
  Float_t         conv_chi2_probability[1000];
  Float_t         conv_vtx_xErr[1000];
  Float_t         conv_vtx_yErr[1000];
  Float_t         conv_vtx_zErr[1000];
  Float_t         conv_tk1_dz[1000];
  Float_t         conv_tk2_dz[1000];
  Float_t         conv_tk1_dzerr[1000];
  Float_t         conv_tk2_dzerr[1000];
  Int_t           conv_tk1_nh[1000];
  Int_t           conv_tk2_nh[1000];
  Int_t           conv_ch1ch2[1000];
  Float_t         conv_tk1_d0[1000];
  Float_t         conv_tk1_pout[1000];
  Float_t         conv_tk1_pin[1000];
  Float_t         conv_tk2_d0[1000];
  Float_t         conv_tk2_pout[1000];
  Float_t         conv_tk2_pin[1000];
  TClonesArray    *conv_vtx;
  TClonesArray    *conv_pair_momentum;
  TClonesArray    *conv_refitted_momentum;
  Int_t           process_id;
  Float_t         weight;
  Float_t         pthat;
  Int_t           gp_n;
  TClonesArray    *gp_p4;
  Short_t         gp_status[2000];
  Short_t         gp_pdgid[2000];
  TClonesArray    *gp_vtx;
  Int_t           vtx_std_n;
  Float_t         vtx_std_x2dof[50];
  TClonesArray    *vtx_std_xyz;
  TClonesArray    *vtx_std_dxdydz;
  TClonesArray    *bs_xyz;
  Float_t         bs_sigmaZ;
  Float_t         bs_x0Error;
  Float_t         bs_y0Error;
  Float_t         bs_z0Error;
  Float_t         bs_sigmaZ0Error;
  
  // List of branches
  TBranch        *b_event;
  TBranch        *b_lumis;
  TBranch        *b_run;
  TBranch        *b_bx;
  TBranch        *b_sc_p4;
  TBranch        *b_pho_n;
  TBranch        *b_pho_feta;
  TBranch        *b_pho_crackcorr;
  TBranch        *b_pho_localcorr;
  TBranch        *b_pho_isEB;
  TBranch        *b_pho_isEE;
  TBranch        *b_pho_isEBGap;
  TBranch        *b_pho_isEEGap;
  TBranch        *b_pho_isEBEEGap;
  TBranch        *b_pho_see;
  TBranch        *b_pho_sieie;
  TBranch        *b_pho_sipip;
  TBranch        *b_pho_sieip;
  TBranch        *b_pho_e1x5;
  TBranch        *b_pho_e2x5;
  TBranch        *b_pho_e3x3;
  TBranch        *b_pho_e5x5;
  TBranch        *b_pho_emaxxtal;
  TBranch        *b_pho_hoe;
  TBranch        *b_pho_h1oe;
  TBranch        *b_pho_h2oe;
  TBranch        *b_pho_r1x5;
  TBranch        *b_pho_r2x5;
  TBranch        *b_pho_r9;
  TBranch        *b_pho_zernike20;
  TBranch        *b_pho_zernike42;
  TBranch        *b_pho_e2nd;
  TBranch        *b_pho_e2x5right;
  TBranch        *b_pho_e2x5left;
  TBranch        *b_pho_e2x5Top;
  TBranch        *b_pho_e2x5bottom;
  TBranch        *b_pho_eright;
  TBranch        *b_pho_eleft;
  TBranch        *b_pho_etop;
  TBranch        *b_pho_ebottom;
  TBranch        *b_pho_e2overe9;
  TBranch        *b_pho_seed_time;
  TBranch        *b_pho_seed_outoftimechi2;
  TBranch        *b_pho_seed_chi2;
  TBranch        *b_pho_seed_recoflag;
  TBranch        *b_pho_seed_severity;
  TBranch        *b_pho_ecalsumetconedr04;
  TBranch        *b_pho_hcalsumetconedr04;
  TBranch        *b_pho_hcal1sumetconedr04;
  TBranch        *b_pho_hcal2sumetconedr04;
  TBranch        *b_pho_trksumptsolidconedr04;
  TBranch        *b_pho_trksumpthollowconedr04;
  TBranch        *b_pho_ntrksolidconedr04;
  TBranch        *b_pho_ntrkhollowconedr04;
  TBranch        *b_pho_ecalsumetconedr03;
  TBranch        *b_pho_hcalsumetconedr03;
  TBranch        *b_pho_hcal1sumetconedr03;
  TBranch        *b_pho_hcal2sumetconedr03;
  TBranch        *b_pho_trksumptsolidconedr03;
  TBranch        *b_pho_trksumpthollowconedr03;
  TBranch        *b_pho_ntrksolidconedr03;
  TBranch        *b_pho_ntrkhollowconedr03;
  TBranch        *b_pho_barrel;
  TBranch        *b_pho_scind;
  TBranch        *b_pho_haspixseed;
  TBranch        *b_pho_hasconvtks;
  TBranch        *b_pho_nconv;
  TBranch        *b_pho_conv_ntracks;
  TBranch        *b_pho_conv_pairinvmass;
  TBranch        *b_pho_conv_paircotthetasep;
  TBranch        *b_pho_conv_eoverp;
  TBranch        *b_pho_conv_zofprimvtxfromtrks;
  TBranch        *b_pho_conv_distofminapproach;
  TBranch        *b_pho_conv_dphitrksatvtx;
  TBranch        *b_pho_conv_dphitrksatecal;
  TBranch        *b_pho_conv_detatrksatecal;
  TBranch        *b_pho_conv_tk1_d0;
  TBranch        *b_pho_conv_tk1_pout;
  TBranch        *b_pho_conv_tk1_pin;
  TBranch        *b_pho_conv_tk2_d0;
  TBranch        *b_pho_conv_tk2_pout;
  TBranch        *b_pho_conv_tk2_pin;
  TBranch        *b_pho_conv_tk1_dz;
  TBranch        *b_pho_conv_tk1_dzerr;
  TBranch        *b_pho_conv_tk1_nh;
  TBranch        *b_pho_conv_tk2_dz;
  TBranch        *b_pho_conv_tk2_dzerr;
  TBranch        *b_pho_conv_tk2_nh;
  TBranch        *b_pho_conv_ch1ch2;
  TBranch        *b_pho_conv_chi2;
  TBranch        *b_pho_conv_chi2_probability;
  TBranch        *b_pho_conv_validvtx;
  TBranch        *b_pho_conv_MVALikelihood;
  TBranch        *b_pho_p4;
  TBranch        *b_pho_calopos;
  TBranch        *b_pho_conv_vtx;
  TBranch        *b_pho_conv_pair_momentum;
  TBranch        *b_pho_conv_refitted_momentum;
  TBranch        *b_pho_conv_vertexcorrected_p4;
  TBranch        *b_conv_n;
  TBranch        *b_conv_p4;
  TBranch        *b_conv_ntracks;
  TBranch        *b_conv_pairinvmass;
  TBranch        *b_conv_paircotthetasep;
  TBranch        *b_conv_eoverp;
  TBranch        *b_conv_distofminapproach;
  TBranch        *b_conv_dphitrksatvtx;
  TBranch        *b_conv_dphitrksatecal;
  TBranch        *b_conv_detatrksatecal;
  TBranch        *b_conv_dxy;
  TBranch        *b_conv_dz;
  TBranch        *b_conv_lxy;
  TBranch        *b_conv_lz;
  TBranch        *b_conv_zofprimvtxfromtrks;
  TBranch        *b_conv_nHitsBeforeVtx;
  TBranch        *b_conv_nSharedHits;
  TBranch        *b_conv_validvtx;
  TBranch        *b_conv_MVALikelihood;
  TBranch        *b_conv_chi2;
  TBranch        *b_conv_chi2_probability;
  TBranch        *b_conv_vtx_xErr;
  TBranch        *b_conv_vtx_yErr;
  TBranch        *b_conv_vtx_zErr;
  TBranch        *b_conv_tk1_dz;
  TBranch        *b_conv_tk2_dz;
  TBranch        *b_conv_tk1_dzerr;
  TBranch        *b_conv_tk2_dzerr;
  TBranch        *b_conv_tk1_nh;
  TBranch        *b_conv_tk2_nh;
  TBranch        *b_conv_ch1ch2;
  TBranch        *b_conv_tk1_d0;
  TBranch        *b_conv_tk1_pout;
  TBranch        *b_conv_tk1_pin;
  TBranch        *b_conv_tk2_d0;
  TBranch        *b_conv_tk2_pout;
  TBranch        *b_conv_tk2_pin;
  TBranch        *b_conv_vtx;
  TBranch        *b_conv_pair_momentum;
  TBranch        *b_conv_refitted_momentum;
  TBranch        *b_process_id;
  TBranch        *b_weight;
  TBranch        *b_pthat;
  TBranch        *b_gp_n;
  TBranch        *b_gp_p4;
  TBranch        *b_gp_status;
  TBranch        *b_gp_pdgid;
  TBranch        *b_gp_vtx;
  TBranch        *b_vtx_std_n;
  TBranch        *b_vtx_std_x2dof;
  TBranch        *b_vtx_std_xyz;
  TBranch        *b_vtx_std_dxdydz;
  TBranch        *b_bs_xyz;
  TBranch        *b_bs_sigmaZ;
  TBranch        *b_bs_x0Error;
  TBranch        *b_bs_y0Error;
  TBranch        *b_bs_z0Error;
  TBranch        *b_bs_sigmaZ0Error;
  
  sdaReader(TFile *currentFile=0);
  virtual ~sdaReader();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init(TFile *currentFile);
  virtual Bool_t   Notify();
  virtual Bool_t   FindLeaf(const char* leafname);
  virtual Bool_t   FindBranch(const char* leafname);
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
  conv_p4 = 0;
  conv_nHitsBeforeVtx = 0;
  conv_vtx = 0;
  conv_pair_momentum = 0;
  conv_refitted_momentum = 0;
  gp_p4 = 0;
  gp_vtx = 0;
  vtx_std_xyz = 0;
  vtx_std_dxdydz = 0;
  bs_xyz = 0;  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("bx", &bx, &b_bx);
  fChain->SetBranchAddress("sc_p4", &sc_p4, &b_sc_p4);
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  fChain->SetBranchAddress("pho_feta", pho_feta, &b_pho_feta);
  fChain->SetBranchAddress("pho_crackcorr", pho_crackcorr, &b_pho_crackcorr);
  fChain->SetBranchAddress("pho_localcorr", pho_localcorr, &b_pho_localcorr);
  fChain->SetBranchAddress("pho_isEB", pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_isEBGap", pho_isEBGap, &b_pho_isEBGap);
  fChain->SetBranchAddress("pho_isEEGap", pho_isEEGap, &b_pho_isEEGap);
  fChain->SetBranchAddress("pho_isEBEEGap", pho_isEBEEGap, &b_pho_isEBEEGap);
  fChain->SetBranchAddress("pho_see", pho_see, &b_pho_see);
  fChain->SetBranchAddress("pho_sieie", pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_sipip", pho_sipip, &b_pho_sipip);
  fChain->SetBranchAddress("pho_sieip", pho_sieip, &b_pho_sieip);
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
  fChain->SetBranchAddress("pho_zernike20", pho_zernike20, &b_pho_zernike20);
  fChain->SetBranchAddress("pho_zernike42", pho_zernike42, &b_pho_zernike42);
  fChain->SetBranchAddress("pho_e2nd", pho_e2nd, &b_pho_e2nd);
  fChain->SetBranchAddress("pho_e2x5right", pho_e2x5right, &b_pho_e2x5right);
  fChain->SetBranchAddress("pho_e2x5left", pho_e2x5left, &b_pho_e2x5left);
  fChain->SetBranchAddress("pho_e2x5Top", pho_e2x5Top, &b_pho_e2x5Top);
  fChain->SetBranchAddress("pho_e2x5bottom", pho_e2x5bottom, &b_pho_e2x5bottom);
  fChain->SetBranchAddress("pho_eright", pho_eright, &b_pho_eright);
  fChain->SetBranchAddress("pho_eleft", pho_eleft, &b_pho_eleft);
  fChain->SetBranchAddress("pho_etop", pho_etop, &b_pho_etop);
  fChain->SetBranchAddress("pho_ebottom", pho_ebottom, &b_pho_ebottom);
  fChain->SetBranchAddress("pho_e2overe9", pho_e2overe9, &b_pho_e2overe9);
  fChain->SetBranchAddress("pho_seed_time", pho_seed_time, &b_pho_seed_time);
  fChain->SetBranchAddress("pho_seed_outoftimechi2", pho_seed_outoftimechi2, &b_pho_seed_outoftimechi2);
  fChain->SetBranchAddress("pho_seed_chi2", pho_seed_chi2, &b_pho_seed_chi2);
  fChain->SetBranchAddress("pho_seed_recoflag", pho_seed_recoflag, &b_pho_seed_recoflag);
  fChain->SetBranchAddress("pho_seed_severity", pho_seed_severity, &b_pho_seed_severity);
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
  fChain->SetBranchAddress("pho_conv_tk1_dzerr", pho_conv_tk1_dzerr, &b_pho_conv_tk1_dzerr);
  fChain->SetBranchAddress("pho_conv_tk1_nh", pho_conv_tk1_nh, &b_pho_conv_tk1_nh);
  fChain->SetBranchAddress("pho_conv_tk2_dz", pho_conv_tk2_dz, &b_pho_conv_tk2_dz);
  fChain->SetBranchAddress("pho_conv_tk2_dzerr", pho_conv_tk2_dzerr, &b_pho_conv_tk2_dzerr);
  fChain->SetBranchAddress("pho_conv_tk2_nh", pho_conv_tk2_nh, &b_pho_conv_tk2_nh);
  fChain->SetBranchAddress("pho_conv_ch1ch2", pho_conv_ch1ch2, &b_pho_conv_ch1ch2);
  fChain->SetBranchAddress("pho_conv_chi2", pho_conv_chi2, &b_pho_conv_chi2);
  fChain->SetBranchAddress("pho_conv_chi2_probability", pho_conv_chi2_probability, &b_pho_conv_chi2_probability);
  fChain->SetBranchAddress("pho_conv_validvtx", pho_conv_validvtx, &b_pho_conv_validvtx);
  fChain->SetBranchAddress("pho_conv_MVALikelihood", pho_conv_MVALikelihood, &b_pho_conv_MVALikelihood);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
  fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
  fChain->SetBranchAddress("pho_conv_vtx", &pho_conv_vtx, &b_pho_conv_vtx);
  fChain->SetBranchAddress("pho_conv_pair_momentum", &pho_conv_pair_momentum, &b_pho_conv_pair_momentum);
  fChain->SetBranchAddress("pho_conv_refitted_momentum", &pho_conv_refitted_momentum, &b_pho_conv_refitted_momentum);
  fChain->SetBranchAddress("pho_conv_vertexcorrected_p4", &pho_conv_vertexcorrected_p4, &b_pho_conv_vertexcorrected_p4);
  fChain->SetBranchAddress("conv_n", &conv_n, &b_conv_n);
  fChain->SetBranchAddress("conv_p4", &conv_p4, &b_conv_p4);
  fChain->SetBranchAddress("conv_ntracks", conv_ntracks, &b_conv_ntracks);
  fChain->SetBranchAddress("conv_pairinvmass", conv_pairinvmass, &b_conv_pairinvmass);
  fChain->SetBranchAddress("conv_paircotthetasep", conv_paircotthetasep, &b_conv_paircotthetasep);
  fChain->SetBranchAddress("conv_eoverp", conv_eoverp, &b_conv_eoverp);
  fChain->SetBranchAddress("conv_distofminapproach", conv_distofminapproach, &b_conv_distofminapproach);
  fChain->SetBranchAddress("conv_dphitrksatvtx", conv_dphitrksatvtx, &b_conv_dphitrksatvtx);
  fChain->SetBranchAddress("conv_dphitrksatecal", conv_dphitrksatecal, &b_conv_dphitrksatecal);
  fChain->SetBranchAddress("conv_detatrksatecal", conv_detatrksatecal, &b_conv_detatrksatecal);
  fChain->SetBranchAddress("conv_dxy", conv_dxy, &b_conv_dxy);
  fChain->SetBranchAddress("conv_dz", conv_dz, &b_conv_dz);
  fChain->SetBranchAddress("conv_lxy", conv_lxy, &b_conv_lxy);
  fChain->SetBranchAddress("conv_lz", conv_lz, &b_conv_lz);
  fChain->SetBranchAddress("conv_zofprimvtxfromtrks", conv_zofprimvtxfromtrks, &b_conv_zofprimvtxfromtrks);
  fChain->SetBranchAddress("conv_nHitsBeforeVtx", &conv_nHitsBeforeVtx, &b_conv_nHitsBeforeVtx);
  fChain->SetBranchAddress("conv_nSharedHits", conv_nSharedHits, &b_conv_nSharedHits);
  fChain->SetBranchAddress("conv_validvtx", conv_validvtx, &b_conv_validvtx);
  fChain->SetBranchAddress("conv_MVALikelihood", conv_MVALikelihood, &b_conv_MVALikelihood);
  fChain->SetBranchAddress("conv_chi2", conv_chi2, &b_conv_chi2);
  fChain->SetBranchAddress("conv_chi2_probability", conv_chi2_probability, &b_conv_chi2_probability);
  fChain->SetBranchAddress("conv_vtx_xErr", conv_vtx_xErr, &b_conv_vtx_xErr);
  fChain->SetBranchAddress("conv_vtx_yErr", conv_vtx_yErr, &b_conv_vtx_yErr);
  fChain->SetBranchAddress("conv_vtx_zErr", conv_vtx_zErr, &b_conv_vtx_zErr);
  fChain->SetBranchAddress("conv_tk1_dz", conv_tk1_dz, &b_conv_tk1_dz);
  fChain->SetBranchAddress("conv_tk2_dz", conv_tk2_dz, &b_conv_tk2_dz);
  fChain->SetBranchAddress("conv_tk1_dzerr", conv_tk1_dzerr, &b_conv_tk1_dzerr);
  fChain->SetBranchAddress("conv_tk2_dzerr", conv_tk2_dzerr, &b_conv_tk2_dzerr);
  fChain->SetBranchAddress("conv_tk1_nh", conv_tk1_nh, &b_conv_tk1_nh);
  fChain->SetBranchAddress("conv_tk2_nh", conv_tk2_nh, &b_conv_tk2_nh);
  fChain->SetBranchAddress("conv_ch1ch2", conv_ch1ch2, &b_conv_ch1ch2);
  fChain->SetBranchAddress("conv_tk1_d0", conv_tk1_d0, &b_conv_tk1_d0);
  fChain->SetBranchAddress("conv_tk1_pout", conv_tk1_pout, &b_conv_tk1_pout);
  fChain->SetBranchAddress("conv_tk1_pin", conv_tk1_pin, &b_conv_tk1_pin);
  fChain->SetBranchAddress("conv_tk2_d0", conv_tk2_d0, &b_conv_tk2_d0);
  fChain->SetBranchAddress("conv_tk2_pout", conv_tk2_pout, &b_conv_tk2_pout);
  fChain->SetBranchAddress("conv_tk2_pin", conv_tk2_pin, &b_conv_tk2_pin);
  fChain->SetBranchAddress("conv_vtx", &conv_vtx, &b_conv_vtx);
  fChain->SetBranchAddress("conv_pair_momentum", &conv_pair_momentum, &b_conv_pair_momentum);
  fChain->SetBranchAddress("conv_refitted_momentum", &conv_refitted_momentum, &b_conv_refitted_momentum);
  fChain->SetBranchAddress("process_id", &process_id, &b_process_id);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
  fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
  fChain->SetBranchAddress("gp_n", &gp_n, &b_gp_n);
  fChain->SetBranchAddress("gp_p4", &gp_p4, &b_gp_p4);
  fChain->SetBranchAddress("gp_status", gp_status, &b_gp_status);
  fChain->SetBranchAddress("gp_pdgid", gp_pdgid, &b_gp_pdgid);
  fChain->SetBranchAddress("gp_vtx", &gp_vtx, &b_gp_vtx);
  fChain->SetBranchAddress("vtx_std_n", &vtx_std_n, &b_vtx_std_n);
  fChain->SetBranchAddress("vtx_std_x2dof", vtx_std_x2dof, &b_vtx_std_x2dof);
  fChain->SetBranchAddress("vtx_std_xyz", &vtx_std_xyz, &b_vtx_std_xyz);
  fChain->SetBranchAddress("vtx_std_dxdydz", &vtx_std_dxdydz, &b_vtx_std_dxdydz);
  fChain->SetBranchAddress("bs_xyz", &bs_xyz, &b_bs_xyz);
  fChain->SetBranchAddress("bs_sigmaZ", &bs_sigmaZ, &b_bs_sigmaZ);
  fChain->SetBranchAddress("bs_x0Error", &bs_x0Error, &b_bs_x0Error);
  fChain->SetBranchAddress("bs_y0Error", &bs_y0Error, &b_bs_y0Error);
  fChain->SetBranchAddress("bs_z0Error", &bs_z0Error, &b_bs_z0Error);
  fChain->SetBranchAddress("bs_sigmaZ0Error", &bs_sigmaZ0Error, &b_bs_sigmaZ0Error);
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

Bool_t sdaReader::FindLeaf(const char* leafname) {
  if (fChain->FindLeaf(leafname)==NULL) return false;
  return true;
}

Bool_t sdaReader::FindBranch(const char* leafname) {
  if (fChain->FindBranch(leafname)==NULL) return false;
  return true;
}
