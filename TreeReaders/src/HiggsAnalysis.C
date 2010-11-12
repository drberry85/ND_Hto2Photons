#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"

#include <iostream>
#include <vector>
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

#include "ND_Hto2Photons/TreeReaders/interface/mpaReader.h"
#include "ND_Hto2Photons/TreeReaders/interface/Selections.h"

using namespace std;

int main(int argc, char * charmass[]) {

    float globalWeight = 35;
  
    vector<pair<string, float> > filesAndWeights;
    vector<string> filelist;
    
    string Mass=charmass[1];
    cout << "Mass is: " << Mass << endl;
    
    if (Mass=="90GeV") {
      filelist.push_back("HiggsAnalysis90GeV.root");
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon90MPA.root",20.493));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA90.root",1.4405));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ90MPA.root",0.9899));
      cout << "Warnging Weights Not Correct!!!!!" << endl;
    } else if (Mass=="110GeV") {
      filelist.push_back("HiggsAnalysis110GeV.root");
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon110MPA.root",20.493));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA110.root",1.4405));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ110MPA.root",0.9899));
    } else if (Mass=="120GeV") {
      filelist.push_back("HiggsAnalysis120GeV.root");
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon120MPA.root",17.173));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA120.root",1.3062));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ120MPA.root",0.4240));
    } else if (Mass=="150GeV") {
      filelist.push_back("HiggsAnalysis150GeV.root");
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon150MPA.root",10.863));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA150.root",1.0985));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ150MPA.root",0.2035));
    } else if (Mass=="") {
      cout << "Calculating all masses!" << endl;
      filelist.push_back("HiggsAnalysis90GeV.root");
      filelist.push_back("HiggsAnalysis110GeV.root");
      filelist.push_back("HiggsAnalysis120GeV.root");
      filelist.push_back("HiggsAnalysis150GeV.root");
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon90MPA.root",20.493));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA90.root",1.4405));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ90MPA.root",0.9899));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon110MPA.root",20.493));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA110.root",1.4405));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ110MPA.root",0.9899));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon120MPA.root",17.173));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA120.root",1.3062));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ120MPA.root",0.4240));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggGluon150MPA.root",10.863));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggVBFMPA150.root",1.0985));
      filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_6_3/MPA/HiggsSignal/HggQQ150MPA.root",0.2035));
    } else {
      cout << "Please Select a Higgs Mass!" << endl << "90GeV" << endl << "110GeV" << endl << "120GeV" << endl << "150GeV" << endl;
    }

    for (unsigned int itMasses=0; itMasses < filelist.size(); itMasses++) {
    
      TFile * outfile = new TFile(filelist[itMasses].c_str(),"RECREATE");

      TH1F * hNPhotons[2];
      TH1F * hLeadEt[3][2];
      TH1F * hSubLeadEt[3][2];
      TH1F * hLeadEta[2];
      TH1F * hSubLeadEta[2];
      TH1F * hLeadR9[3][2];
      TH1F * hSubLeadR9[3][2];
      TH1F * hLeadHoE[3][2];
      TH1F * hSubLeadHoE[3][2];
      TH1F * hLeadTrkPtSumSolid03[3][2];
      TH1F * hSubLeadTrkPtSumSolid03[3][2];
      TH1F * hLeadEcalPtSumSolid03[3][2];
      TH1F * hSubLeadEcalPtSumSolid03[3][2];
      TH1F * hLeadHcalPtSumSolid03[3][2];
      TH1F * hSubLeadHcalPtSumSolid03[3][2];
      TH1F * hLeadSigmaIetaIeta[3][2];
      TH1F * hSubLeadSigmaIetaIeta[3][2];

      TH1F*  hLeadZPV_[3][2];
      TH1F*  hSubLeadZPV_[3][2];
      TH1F*  hLeadDzPV_[3][2];
      TH1F*  hSubLeadDzPV_[3][2];

      TH1D* h_mass_2gamma[2][2];
      TH1D* h_pt_2gamma[2][2];
      TH1D* h_pz_2gamma[2][2];
      TH1D* h_eta_2gamma[2][2];
      TH1D* h_CosThetaStar[2][2];

      TH1D* h_mass_2gamma_2gold[2][2];
      TH1D* h_mass_2gamma_2goldConv[2][2];
      TH1D* h_pt_2gamma_2gold[2][2];
      TH1D* h_pz_2gamma_2gold[2][2];
      TH1D* h_eta_2gamma_2gold[2][2];
      TH1D* h_CosThetaStar_2gold[2][2];

      TH1D* h_mass_2gamma_1conv[2][2];
      TH1D* h_pt_2gamma_1conv[2][2];
      TH1D* h_pz_2gamma_1conv[2][2];
      TH1D* h_eta_2gamma_1conv[2][2];
      TH1D* h_CosThetaStar_1conv[2][2];

      TH1D* h_mass_2gamma_2conv[2][2];
      TH1D* h_pt_2gamma_2conv[2][2];
      TH1D* h_pz_2gamma_2conv[2][2];
      TH1D* h_eta_2gamma_2conv[2][2];
      TH1D* h_CosThetaStar_2conv[2][2];

      TH1D* h_mass_2gamma_leftover[2][2];
      TH1D* h_pt_2gamma_leftover[2][2];
      TH1D* h_pz_2gamma_leftover[2][2];
      TH1D* h_eta_2gamma_leftover[2][2];
      TH1D* h_CosThetaStar_leftover[2][2];


      TH2F*  h2_convVtxRvsZBarrel_[2];

      hNPhotons[0] = new TH1F("hNPhotonsAll","Num of photons in the event: all candidates",20,-0.5,19.5);
      hNPhotons[1] = new TH1F("hNPhotonsSel","Num of photons in the event: selected candidates",20,-0.5,19.5);

      hLeadEt[0][0] = new TH1F("leadPhoEtAll_allEcal","leading photon Et, all candidates: all ECAL",100,0.,100.);
      hLeadEt[1][0] = new TH1F("leadPhoEtAll_Barrel","leading photon Et, all candidates:  Barrel",100,0.,100.);
      hLeadEt[2][0] = new TH1F("leadPhoEtAll_Endcap","leading photon Et, all candidates:  Endcap",100,0.,100.);

      hSubLeadEt[0][0] = new TH1F("subleadPhoEtAll_allEcal","sub-leading photon Et, all candidates: all ECAL",100,0.,100.);
      hSubLeadEt[1][0] = new TH1F("subleadPhoEtAll_Barrel","sub-leading photon Et, all candidates:  Barrel",100,0.,100.);
      hSubLeadEt[2][0] = new TH1F("subleadPhoEtAll_Endcap","sub-leading photon Et, all candidates:  Endcap",100,0.,100.);

      hLeadEta[0] = new TH1F("leadPhoEtaAll_allEcal","leading photon Eta, all candidates: all ECAL",100,-3.,3.);
      hSubLeadEta[0] = new TH1F("subleadPhoEtaAll_allEcal","sub-leading photon Eta, all candidates: all ECAL",100,-3.,3.);

      hLeadR9[0][0] = new TH1F("leadPhoR9All_allEcal","leading photon R9, all candidates: all ECAL",100,0.,1.1);
      hLeadR9[1][0] = new TH1F("leadPhoR9All_Barrel","leading photon R9, all candidates: Barrel",100,0.,1.1);
      hLeadR9[2][0] = new TH1F("leadPhoR9All_Endcap","leading photon R9, all candidates: Endcap",100,0.,1.1);

      hSubLeadR9[0][0] = new TH1F("subleadPhoR9All_allEcal","sub-leading photon R9, all candidates: all ECAL",100,0.,1.1);
      hSubLeadR9[1][0] = new TH1F("subleadPhoR9All_Barrel","sub-leading photon R9, all candidates: Barrel",100,0.,1.1);
      hSubLeadR9[2][0] = new TH1F("subleadPhoR9All_Endcap","sub-leading photon R9, all candidates: Endcap",100,0.,1.1);

      hLeadHoE[0][0] = new TH1F("leadPhoHoEAll_allEcal","leading photon HoE, all candidates: all ECAL",100,0.,1.);
      hLeadHoE[1][0] = new TH1F("leadPhoHoEAll_Barrel","leading photon HoE, all candidates: Barrel",100,0.,1.);
      hLeadHoE[2][0] = new TH1F("leadPhoHoEAll_Endcap","leading photon HoE, all candidates: Endcap",100,0.,1.);

      hSubLeadHoE[0][0] = new TH1F("subleadPhoHoEAll_allEcal","sub-leading photon HoE, all candidates: all ECAL",100,0.,1.);
      hSubLeadHoE[1][0] = new TH1F("subleadPhoHoEAll_Barrel","sub-leading photon HoE, all candidates: Barrel",100,0.,1.);
      hSubLeadHoE[2][0] = new TH1F("subleadPhoHoEAll_Endcap","sub-leading photon HoE, all candidates: Endcap",100,0.,1.);

      hLeadTrkPtSumSolid03[0][0] = new TH1F("leadPhoTrkPtSumSolid03All_allEcal","leading photon trk pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hLeadTrkPtSumSolid03[1][0] = new TH1F("leadPhoTrkPtSumSolid03All_Barrel","leading photon trk pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hLeadTrkPtSumSolid03[2][0] = new TH1F("leadPhoTrkPtSumSolid03All_Endcap","leading photon trk pt sum dr=03, all candidates: Endcap",100,0.,50.);

      hSubLeadTrkPtSumSolid03[0][0] = new TH1F("subleadPhoTrkPtSumSolid03All_allEcal","sub-leading photon trk pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hSubLeadTrkPtSumSolid03[1][0] = new TH1F("subleadPhoTrkPtSumSolid03All_Barrel","sub-leading photon trk pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hSubLeadTrkPtSumSolid03[2][0] = new TH1F("subleadPhoTrkPtSumSolid03All_Endcap","sub-leading photon trk pt sum dr=03, all candidates: Endcap",100,0.,50.);

      hLeadEcalPtSumSolid03[0][0] = new TH1F("leadPhoEcalPtSumSolid03All_allEcal","leading photon ecal pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hLeadEcalPtSumSolid03[1][0] = new TH1F("leadPhoEcalPtSumSolid03All_Barrel","leading photon ecal pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hLeadEcalPtSumSolid03[2][0] = new TH1F("leadPhoEcalPtSumSolid03All_Endcap","leading photon ecal pt sum dr=03, all candidates: Endcap",100,0.,50.);

      hSubLeadEcalPtSumSolid03[0][0] = new TH1F("subleadPhoEcalPtSumSolid03All_allEcal","sub-leading photon ecal pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hSubLeadEcalPtSumSolid03[1][0] = new TH1F("subleadPhoEcalPtSumSolid03All_Barrel","sub-leading photon ecal pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hSubLeadEcalPtSumSolid03[2][0] = new TH1F("subleadPhoEcalPtSumSolid03All_Endcap","sub-leading photon ecal pt sum dr=03, all candidates: Endcap",100,0.,50.);

      hLeadHcalPtSumSolid03[0][0] = new TH1F("leadPhoHcalPtSumSolid03All_allEcal","leading photon hcal pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hLeadHcalPtSumSolid03[1][0] = new TH1F("leadPhoHcalPtSumSolid03All_Barrel","leading photon hcal pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hLeadHcalPtSumSolid03[2][0] = new TH1F("leadPhoHcalPtSumSolid03All_Endcap","leading photon hcal pt sum dr=03, all candidates: Endcap",100,0.,50.);

      hSubLeadHcalPtSumSolid03[0][0] = new TH1F("subleadPhoHcalPtSumSolid03All_allEcal","sub-leading photon hcal pt sum dr=03, all candidates: all ECAL",100,0.,50.);
      hSubLeadHcalPtSumSolid03[1][0] = new TH1F("subleadPhoHcalPtSumSolid03All_Barrel","sub-leading photon hcal pt sum dr=03, all candidates: Barrel",100,0.,50.);
      hSubLeadHcalPtSumSolid03[2][0] = new TH1F("subleadPhoHcalPtSumSolid03All_Endcap","sub-leading photon hcal pt sum dr=03, all candidates: Endcap",100,0.,50.);

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

      hLeadDzPV_[0][0]= new TH1F("leadPhoDZPVAll_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20.0, 20.0);
      hLeadDzPV_[1][0]= new TH1F("leadPhoDZPVAll_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20.0, 20.0);
      hLeadDzPV_[2][0]= new TH1F("leadPhoDZPVAll_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates", 100, -20.0, 20.0);

      hSubLeadDzPV_[0][0]= new TH1F("subleadPhoDZPVAll_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20.0, 20.0);
      hSubLeadDzPV_[1][0]= new TH1F("subleadPhoDZPVAll_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20.0, 20.0);
      hSubLeadDzPV_[2][0]= new TH1F("subleadPhoDZPVAll_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Endcap", 100, -20.0, 20.0);

      // candidates after pt and eta selection
      hLeadEt[0][1] = new TH1F("leadPhoEtSel_allEcal","leading photon Et, selected candidates: all ECAL",100,0.,100.);
      hLeadEt[1][1] = new TH1F("leadPhoEtSel_Barrel","leading photon Et, selected candidates:  Barrel",100,0.,100.);
      hLeadEt[2][1] = new TH1F("leadPhoEtSel_Endcap","leading photon Et, selected candidates:  Endcap",100,0.,100.);

      hSubLeadEt[0][1] = new TH1F("subleadPhoEtSel_allEcal","sub-leading photon Et, selected candidates: all ECAL",100,0.,100.);
      hSubLeadEt[1][1] = new TH1F("subleadPhoEtSel_Barrel","sub-leading photon Et, selected candidates:  Barrel",100,0.,100.);
      hSubLeadEt[2][1] = new TH1F("subleadPhoEtSel_Endcap","sub-leading photon Et, selected candidates:  Endcap",100,0.,100.);

      hLeadEta[1] = new TH1F("leadPhoEtaSel_allEcal","leading photon Eta, selected candidates: all ECAL",100,-3.,3.);
      hSubLeadEta[1] = new TH1F("subleadPhoEtaSel_allEcal","sub-leading photon Eta, selected candidates: all ECAL",100,-3.,3.);

      hLeadR9[0][1] = new TH1F("leadPhoR9Sel_allEcal","leading photon R9, selected candidates: all ECAL",100,0.,1.1);
      hLeadR9[1][1] = new TH1F("leadPhoR9Sel_Barrel","leading photon R9, selected candidates: Barrel",100,0.,1.1);
      hLeadR9[2][1] = new TH1F("leadPhoR9Sel_Endcap","leading photon R9, selected candidates: Endcap",100,0.,1.1);

      hSubLeadR9[0][1] = new TH1F("subleadPhoR9Sel_allEcal","sub-leading photon R9, selected candidates: all ECAL",100,0.,1.1);
      hSubLeadR9[1][1] = new TH1F("subleadPhoR9Sel_Barrel","sub-leading photon R9, selected candidates: Barrel",100,0.,1.1);
      hSubLeadR9[2][1] = new TH1F("subleadPhoR9Sel_Endcap","sub-leading photon R9, selected candidates: Endcap",100,0.,1.1);

      hLeadHoE[0][1] = new TH1F("leadPhoHoESel_allEcal","leading photon HoE, selected candidates: all ECAL",100,0.,1.);
      hLeadHoE[1][1] = new TH1F("leadPhoHoESel_Barrel","leading photon HoE, selected candidates: Barrel",100,0.,1.);
      hLeadHoE[2][1] = new TH1F("leadPhoHoESel_Endcap","leading photon HoE, selected candidates: Endcap",100,0.,1.);

      hSubLeadHoE[0][1] = new TH1F("subleadPhoHoESel_allEcal","sub-leading photon HoE, selected candidates: all ECAL",100,0.,1.);
      hSubLeadHoE[1][1] = new TH1F("subleadPhoHoESel_Barrel","sub-leading photon HoE, selected candidates: Barrel",100,0.,1.);
      hSubLeadHoE[2][1] = new TH1F("subleadPhoHoESel_Endcap","sub-leading photon HoE, selected candidates: Endcap",100,0.,1.);

      hLeadTrkPtSumSolid03[0][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_allEcal","leading photon trk pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hLeadTrkPtSumSolid03[1][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_Barrel","leading photon trk pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hLeadTrkPtSumSolid03[2][1] = new TH1F("leadPhoTrkPtSumSolid03Sel_Endcap","leading photon trk pt sum dr=03, selected candidates: Endcap",100,0.,50.);

      hSubLeadTrkPtSumSolid03[0][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_allEcal","sub-leading photon trk pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hSubLeadTrkPtSumSolid03[1][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_Barrel","sub-leading photon trk pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hSubLeadTrkPtSumSolid03[2][1] = new TH1F("subleadPhoTrkPtSumSolid03Sel_Endcap","sub-leading photon trk pt sum dr=03, selected candidates: Endcap",100,0.,50.);

      hLeadEcalPtSumSolid03[0][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_allEcal","leading photon ecal pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hLeadEcalPtSumSolid03[1][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_Barrel","leading photon ecal pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hLeadEcalPtSumSolid03[2][1] = new TH1F("leadPhoEcalPtSumSolid03Sel_Endcap","leading photon ecal pt sum dr=03, selected candidates: Endcap",100,0.,50.);

      hSubLeadEcalPtSumSolid03[0][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_allEcal","sub-leading photon ecal pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hSubLeadEcalPtSumSolid03[1][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_Barrel","sub-leading photon ecal pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hSubLeadEcalPtSumSolid03[2][1] = new TH1F("subleadPhoEcalPtSumSolid03Sel_Endcap","sub-leading photon ecal pt sum dr=03, selected candidates: Endcap",100,0.,50.);

      hLeadHcalPtSumSolid03[0][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_allEcal","leading photon hcal pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hLeadHcalPtSumSolid03[1][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_Barrel","leading photon hcal pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hLeadHcalPtSumSolid03[2][1] = new TH1F("leadPhoHcalPtSumSolid03Sel_Endcap","leading photon hcal pt sum dr=03, selected candidates: Endcap",100,0.,50.);

      hSubLeadHcalPtSumSolid03[0][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_allEcal","sub-leading photon hcal pt sum dr=03, selected candidates: all ECAL",100,0.,50.);
      hSubLeadHcalPtSumSolid03[1][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_Barrel","sub-leading photon hcal pt sum dr=03, selected candidates: Barrel",100,0.,50.);
      hSubLeadHcalPtSumSolid03[2][1] = new TH1F("subleadPhoHcalPtSumSolid03Sel_Endcap","sub-leading photon hcal pt sum dr=03, selected candidates: Endcap",100,0.,50.);

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

      hLeadDzPV_[0][1]= new TH1F("leadPhoDZPVSel_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20.0, 20.0);
      hLeadDzPV_[1][1]= new TH1F("leadPhoDZPVSel_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20.0, 20.0);
      hLeadDzPV_[2][1]= new TH1F("leadPhoDZPVSel_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20.0, 20.0);

      hSubLeadDzPV_[0][1]= new TH1F("subleadPhoDZPVSel_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20.0, 20.0);
      hSubLeadDzPV_[1][1]= new TH1F("subleadPhoDZPVSel_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20.0, 20.0);
      hSubLeadDzPV_[2][1]= new TH1F("subleadPhoDZPVSel_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20.0, 20.0);

      // diphoton system
      h_mass_2gamma[0][0]        = new TH1D("h_mass_2gammaAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates", 160, 80.0, 160.0);
      h_pt_2gamma[0][0]          = new TH1D("h_pt_2gammaAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates", 200, 0., 200.0);
      h_pz_2gamma[0][0]          = new TH1D("h_pz_2gammaAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates", 100, -1000., 1000.0);
      h_eta_2gamma[0][0]         = new TH1D("h_eta_2gammaAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates", 160, -8.0, 8.0);
      h_CosThetaStar[0][0]       = new TH1D("h_CosThetaStarAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates", 60, 0., 1.);

      h_mass_2gamma[0][1]        = new TH1D("h_mass_2gammaAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates", 160, 80.0, 160.0);
      h_pt_2gamma[0][1]          = new TH1D("h_pt_2gammaAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates", 200, 0., 200.0);
      h_pz_2gamma[0][1]          = new TH1D("h_pz_2gammaAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates", 100, -1000., 1000.0);
      h_eta_2gamma[0][1]         = new TH1D("h_eta_2gammaAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates", 160, -8.0, 8.0);
      h_CosThetaStar[0][1]       = new TH1D("h_CosThetaStarAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates", 60, 0., 1.);
    
      h_mass_2gamma[1][0]        = new TH1D("h_mass_2gammaSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates", 160, 80.0, 160.0);
      h_pt_2gamma[1][0]          = new TH1D("h_pt_2gammaSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates", 200, 0., 200.0);
      h_pz_2gamma[1][0]          = new TH1D("h_pz_2gammaSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates", 100, -1000., 1000.0);
      h_eta_2gamma[1][0]         = new TH1D("h_eta_2gammaSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates", 160, -8.0, 8.0);
      h_CosThetaStar[1][0]       = new TH1D("h_CosThetaStarSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates", 60, 0., 1.);

      h_mass_2gamma[1][1]        = new TH1D("h_mass_2gammaSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates", 160, 80.0, 160.0);
      h_pt_2gamma[1][1]          = new TH1D("h_pt_2gammaSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates", 200, 0., 200.0);
      h_pz_2gamma[1][1]          = new TH1D("h_pz_2gammaSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates", 100, -1000., 1000.0);
      h_eta_2gamma[1][1]         = new TH1D("h_eta_2gammaSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates", 160, -8.0, 8.0);
      h_CosThetaStar[1][1]       = new TH1D("h_CosThetaStarSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates", 60, 0., 1.);
    
      h_mass_2gamma_2gold[0][0]        = new TH1D("h_mass_2gammaGoldenAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all golden barrel candidates", 160, 80.0, 160.0);
      h_pt_2gamma_2gold[0][0]          = new TH1D("h_pt_2gammaGoldenAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all golden barrel candidates", 200, 0., 200.0);
      h_pz_2gamma_2gold[0][0]          = new TH1D("h_pz_2gammaGoldenAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all golden barrel candidates", 100, -1000., 1000.0);
      h_eta_2gamma_2gold[0][0]         = new TH1D("h_eta_2gammaGoldenAllEB","Di-Photon #eta ;#eta(2#gamma) all golden barrel candidates", 160, -8.0, 8.0);
      h_CosThetaStar_2gold[0][0]       = new TH1D("h_CosThetaStarGoldenAllEB","cos#theta^{*};cos#theta^{*} all golden barrel candidates", 60, 0., 1.);

      h_mass_2gamma_2gold[0][1]        = new TH1D("h_mass_2gammaGoldenAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all golden endcap candidates", 160, 80.0, 160.0);
      h_pt_2gamma_2gold[0][1]          = new TH1D("h_pt_2gammaGoldenAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all golden endcap candidates", 200, 0., 200.0);
      h_pz_2gamma_2gold[0][1]          = new TH1D("h_pz_2gammaGoldenAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all golden endcap candidates", 100, -1000., 1000.0);
      h_eta_2gamma_2gold[0][1]         = new TH1D("h_eta_2gammaGoldenAllEE","Di-Photon #eta ;#eta(2#gamma) all golden endcap candidates", 160, -8.0, 8.0);
      h_CosThetaStar_2gold[0][1]       = new TH1D("h_CosThetaStarGoldenAllEE","cos#theta^{*};cos#theta^{*} all golden endcap candidates", 60, 0., 1.);

      h_mass_2gamma_2gold[1][0]        = new TH1D("h_mass_2gammaGoldenSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected barrel candidates", 160, 80.0, 160.0);
      h_mass_2gamma_2goldConv[1][0]        = new TH1D("h_mass_2gammaGoldenConvSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected barrel candidates", 160, 80.0, 160.0);
      h_pt_2gamma_2gold[1][0]          = new TH1D("h_pt_2gammaGoldenSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden selected barrel candidates", 200, 0., 200.0);
      h_pz_2gamma_2gold[1][0]          = new TH1D("h_pz_2gammaGoldenSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden selected barrel candidates", 100, -1000., 1000.0);
      h_eta_2gamma_2gold[1][0]         = new TH1D("h_eta_2gammaGoldenSelEB","Di-Photon #eta ;#eta(2#gamma) golden selected barrel candidates", 160, -8.0, 8.0);
      h_CosThetaStar_2gold[1][0]       = new TH1D("h_CosThetaStarGoldenSelEB","cos#theta^{*};cos#theta^{*} golden selected barrel candidates", 60, 0., 1.);

      h_mass_2gamma_2gold[1][1]        = new TH1D("h_mass_2gammaGoldenSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected endcap candidates", 160, 80.0, 160.0);
      h_mass_2gamma_2goldConv[1][1]        = new TH1D("h_mass_2gammaGoldenConvSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) golden selected endcap candidates", 160, 80.0, 160.0);
      h_pt_2gamma_2gold[1][1]          = new TH1D("h_pt_2gammaGoldenSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) golden selected endcap candidates", 200, 0., 200.0);
      h_pz_2gamma_2gold[1][1]          = new TH1D("h_pz_2gammaGoldenSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) golden selected endcap candidates", 100, -1000., 1000.0);
      h_eta_2gamma_2gold[1][1]         = new TH1D("h_eta_2gammaGoldenSelEE","Di-Photon #eta ;#eta(2#gamma) golden selected endcap candidates", 160, -8.0, 8.0);
      h_CosThetaStar_2gold[1][1]       = new TH1D("h_CosThetaStarGoldenSelEE","cos#theta^{*};cos#theta^{*} golden selected endcap candidates", 60, 0., 1.);

      h_mass_2gamma_1conv[0][0]        = new TH1D("h_mass_2gamma1convAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates with one conversion", 160, 80.0, 160.0);
      h_pt_2gamma_1conv[0][0]          = new TH1D("h_pt_2gamma1convAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates with one conversion", 200, 0., 200.0);
      h_pz_2gamma_1conv[0][0]          = new TH1D("h_pz_2gamma1convAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates with one conversion", 100, -1000., 1000.0);
      h_eta_2gamma_1conv[0][0]         = new TH1D("h_eta_2gamma1convAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates with one conversion", 160, -8.0, 8.0);
      h_CosThetaStar_1conv[0][0]       = new TH1D("h_CosThetaStar1convAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates with one conversion", 60, 0., 1.);

      h_mass_2gamma_1conv[0][1]        = new TH1D("h_mass_2gamma1convAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates with one conversion", 160, 80.0, 160.0);
      h_pt_2gamma_1conv[0][1]          = new TH1D("h_pt_2gamma1convAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates with one conversion", 200, 0., 200.0);
      h_pz_2gamma_1conv[0][1]          = new TH1D("h_pz_2gamma1convAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates with one conversion", 100, -1000., 1000.0);
      h_eta_2gamma_1conv[0][1]         = new TH1D("h_eta_2gamma1convAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates with one conversion", 160, -8.0, 8.0);
      h_CosThetaStar_1conv[0][1]       = new TH1D("h_CosThetaStar1convAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates with one conversion", 60, 0., 1.);

      h_mass_2gamma_1conv[1][0]        = new TH1D("h_mass_2gamma1convSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with one conversion", 160, 80.0, 160.0);
      h_pt_2gamma_1conv[1][0]          = new TH1D("h_pt_2gamma1convSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with one conversion", 200, 0., 200.0);
      h_pz_2gamma_1conv[1][0]          = new TH1D("h_pz_2gamma1convSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with one conversion", 100, -1000., 1000.0);
      h_eta_2gamma_1conv[1][0]         = new TH1D("h_eta_2gamma1convSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with one conversion", 160, -8.0, 8.0);
      h_CosThetaStar_1conv[1][0]       = new TH1D("h_CosThetaStar1convSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with one conversion", 60, 0., 1.);

      h_mass_2gamma_1conv[1][1]        = new TH1D("h_mass_2gamma1convSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with one conversion", 160, 80.0, 160.0);
      h_pt_2gamma_1conv[1][1]          = new TH1D("h_pt_2gamma1convSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with one conversion", 200, 0., 200.0);
      h_pz_2gamma_1conv[1][1]          = new TH1D("h_pz_2gamma1convSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with one conversion", 100, -1000., 1000.0);
      h_eta_2gamma_1conv[1][1]         = new TH1D("h_eta_2gamma1convSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with one conversion", 160, -8.0, 8.0);
      h_CosThetaStar_1conv[1][1]       = new TH1D("h_CosThetaStar1convSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with one conversion", 60, 0., 1.);

      h_mass_2gamma_2conv[0][0]        = new TH1D("h_mass_2gamma2convAllEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all barrel candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_2conv[0][0]          = new TH1D("h_pt_2gamma2convAllEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) all barrel candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_2conv[0][0]          = new TH1D("h_pz_2gamma2convAllEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all barrel candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_2conv[0][0]         = new TH1D("h_eta_2gamma2convAllEB","Di-Photon #eta ;#eta(2#gamma) all barrel candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_2conv[0][0]       = new TH1D("h_CosThetaStar2convAllEB","cos#theta^{*};cos#theta^{*} all barrel candidates with two conversions", 60, 0., 1.);

      h_mass_2gamma_2conv[0][1]        = new TH1D("h_mass_2gamma2convAllEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) all endcap candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_2conv[0][1]          = new TH1D("h_pt_2gamma2convAllEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) all endcap candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_2conv[0][1]          = new TH1D("h_pz_2gamma2convAllEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) all endcap candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_2conv[0][1]         = new TH1D("h_eta_2gamma2convAllEE","Di-Photon #eta ;#eta(2#gamma) all endcap candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_2conv[0][1]       = new TH1D("h_CosThetaStar2convAllEE","cos#theta^{*};cos#theta^{*} all endcap candidates with two conversions", 60, 0., 1.);

      h_mass_2gamma_2conv[1][0]        = new TH1D("h_mass_2gamma2convSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_2conv[1][0]          = new TH1D("h_pt_2gamma2convSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_2conv[1][0]          = new TH1D("h_pz_2gamma2convSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_2conv[1][0]         = new TH1D("h_eta_2gamma2convSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_2conv[1][0]       = new TH1D("h_CosThetaStar2convSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with two conversions", 60, 0., 1.);

      h_mass_2gamma_2conv[1][1]        = new TH1D("h_mass_2gamma2convSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_2conv[1][1]          = new TH1D("h_pt_2gamma2convSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_2conv[1][1]          = new TH1D("h_pz_2gamma2convSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_2conv[1][1]         = new TH1D("h_eta_2gamma2convSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_2conv[1][1]       = new TH1D("h_CosThetaStar2convSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with two conversions", 60, 0., 1.);

      h_mass_2gamma_leftover[1][0]        = new TH1D("h_mass_2gammaleftoverSelEB", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected barrel candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_leftover[1][0]          = new TH1D("h_pt_2gammaleftoverSelEB","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected barrel candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_leftover[1][0]          = new TH1D("h_pz_2gammaleftoverSelEB","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected barrel candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_leftover[1][0]         = new TH1D("h_eta_2gammaleftoverSelEB","Di-Photon #eta ;#eta(2#gamma) selected barrel candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_leftover[1][0]       = new TH1D("h_CosThetaStarleftoverSelEB","cos#theta^{*};cos#theta^{*} selected barrel candidates with two conversions", 60, 0., 1.);

      h_mass_2gamma_leftover[1][1]        = new TH1D("h_mass_2gammaleftoverSelEE", "Di-photon invariant mass ;M_{#gamma#gamma} (GeV) selected endcap candidates with two conversions", 160, 80.0, 160.0);
      h_pt_2gamma_leftover[1][1]          = new TH1D("h_pt_2gammaleftoverSelEE","Di-photon P_{T} ;PT_{2#gamma} (GeV) selected endcap candidates with two conversions", 200, 0., 200.0);
      h_pz_2gamma_leftover[1][1]          = new TH1D("h_pz_2gammaleftoverSelEE","Di-photon P_{T} ;Pz_{2#gamma} (GeV) selected endcap candidates with two conversions", 100, -1000., 1000.0);
      h_eta_2gamma_leftover[1][1]         = new TH1D("h_eta_2gammaleftoverSelEE","Di-Photon #eta ;#eta(2#gamma) selected endcap candidates with two conversions", 160, -8.0, 8.0);
      h_CosThetaStar_leftover[1][1]       = new TH1D("h_CosThetaStarleftoverSelEE","cos#theta^{*};cos#theta^{*} selected endcap candidates with two conversions", 60, 0., 1.);





      h2_convVtxRvsZBarrel_[0] =   new TH2F("convVtxRvsZBarrelAll"," Photon  conversion vtx position all candidates Barrel",200, 0., 280., 200, 0., 80.);
      h2_convVtxRvsZBarrel_[1] =   new TH2F("convVtxRvsZBarrelSel"," Photon  conversion vtx position selected candidates Barrel",200, 0., 280., 200, 0., 80.);

      for (vector<pair<string, float> >::iterator itFile = filesAndWeights.begin(); itFile != filesAndWeights.end(); ++itFile) {

        string file = itFile->first;
        float weight = itFile->second * globalWeight;

        TFile * currentFile = new TFile(file.c_str());
        TTree * Analysis = (TTree *) currentFile->Get("NTuples/Analysis");

        weight *= 50000/Analysis->GetEntries();
        
        cout << "Reading the tree in file " << file << endl;
        mpaReader currentTree(Analysis);
      
        Long64_t nentries = currentTree.fChain->GetEntries();

        for ( Long64_t i = 0; i < nentries; ++i ) {
          int iDetector = 0;

          currentTree.GetEntry(i);

          if (currentTree.nPhotons<2) continue;

          /*if (currentTree.pt[1]>currentTree.pt[0]) {
            cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
            cout << "LeadPt is " << currentTree.pt[0] << endl;
            cout << "SubLeadPt is " << currentTree.pt[1] << endl;
            }*/
        
          hNPhotons[0]->Fill(currentTree.nPhotons,weight);
          if (currentTree.isEB[0]) iDetector=1;
          if (currentTree.isEE[0]) iDetector=2;

          hLeadEt[iDetector][0]->Fill(currentTree.et[0],weight);
          hLeadEta[0]->Fill(currentTree.eta[0],weight);
          hLeadR9[iDetector][0]->Fill(currentTree.r9[0],weight);

          if (currentTree.isEB[1]) iDetector=1;
          if (currentTree.isEE[1]) iDetector=2;

          hSubLeadEt[iDetector][0]->Fill(currentTree.et[1],weight);
          hSubLeadEta[0]->Fill(currentTree.eta[1],weight);
          hSubLeadR9[iDetector][0]->Fill(currentTree.r9[1],weight);

          TLorentzVector VLead( currentTree.momentumX[0],   currentTree.momentumY[0],  currentTree.momentumZ[0], currentTree.energy[0]);
          TLorentzVector VSubLead( currentTree.momentumX[1],   currentTree.momentumY[1],  currentTree.momentumZ[1], currentTree.energy[1]);
          TLorentzVector VSum=VLead+VSubLead;
          double InvMass=fabs(VSum.M());

          // calculate Cos_theta_star 
          double beta_b  = VSum.Beta();
          double gamma_b = VSum.Gamma();
          TVector3 directionV= VSum.Vect().Unit();

          // TVector3 directionV = VSum.Vect()*(1/VSum.Mag());
          TVector3 CrossVLead=VLead.Vect().Cross(directionV);
          double DotVLeadValue=VLead.Vect().Dot(directionV);
          double CrossVLeadValue=sqrt(CrossVLead.x()*CrossVLead.x()+CrossVLead.y()*CrossVLead.y()+CrossVLead.z()*CrossVLead.z());
          double sin_theta =  CrossVLeadValue/VLead.E();
          double cos_theta =  DotVLeadValue/VLead.E(); 
          double tg_thetas = sin_theta/(gamma_b*(cos_theta-beta_b));
          double cos_thetastar= 1.0/sqrt(1.0+tg_thetas*tg_thetas);

          int HiggsType = 0;
          if (currentTree.isEB[0] && currentTree.isEB[1]) HiggsType=0;
          if ((currentTree.isEB[0] && currentTree.isEE[1]) || (currentTree.isEE[0] && currentTree.isEB[1])) HiggsType=1;
          h_mass_2gamma[0][HiggsType]->Fill(InvMass,weight);
          h_pt_2gamma[0][HiggsType]->Fill(VSum.Pt(),weight);
          h_pz_2gamma[0][HiggsType]->Fill(VSum.Pz(),weight);
          h_eta_2gamma[0][HiggsType]->Fill(VSum.Eta(),weight);
          h_CosThetaStar[0][HiggsType]->Fill(cos_thetastar,weight);

          if (currentTree.r9[0]>.93 && currentTree.r9[1]>.93) {
            h_mass_2gamma_2gold[0][HiggsType]->Fill(InvMass,weight);
            h_pt_2gamma_2gold[0][HiggsType]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2gold[0][HiggsType]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2gold[0][HiggsType]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_2gold[0][HiggsType]->Fill(cos_thetastar,weight);
          }

          if (((currentTree.isConverted[0]==1 && currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005) || (currentTree.isConverted[1]==1 && currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>0.0005)) &&
              !(currentTree.isConverted[0]==1 && currentTree.isConverted[1]==1)) {
            h_mass_2gamma_1conv[0][HiggsType]->Fill(InvMass,weight);
            h_pt_2gamma_1conv[0][HiggsType]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_1conv[0][HiggsType]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_1conv[0][HiggsType]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_1conv[0][HiggsType]->Fill(cos_thetastar,weight);
          }

          if (currentTree.isConverted[0]==1 && currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005 && currentTree.isConverted[1]==1 && currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>.0005) {
            h_mass_2gamma_2conv[0][HiggsType]->Fill(InvMass,weight);
            h_pt_2gamma_2conv[0][HiggsType]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2conv[0][HiggsType]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2conv[0][HiggsType]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_2conv[0][HiggsType]->Fill(cos_thetastar,weight);
          }
          
          if (currentTree.pt[0]<30) continue;
          if (currentTree.pt[1]<25) continue;
          if (currentTree.eta[0]>2.5 || currentTree.eta[1]>2.5) continue;
          if (!(looseId(currentTree.pt[0],currentTree.ecalRecHitSumEtConeDR04[0],currentTree.hcalTowerSumEtConeDR04[0],currentTree.trkSumPtHollowConeDR04[0],(bool) currentTree.isEB[0],(bool) currentTree.isEE[0],currentTree.sigmaIetaIeta[0]))) continue;
          if (!(looseId(currentTree.pt[1],currentTree.ecalRecHitSumEtConeDR04[1],currentTree.hcalTowerSumEtConeDR04[1],currentTree.trkSumPtHollowConeDR04[1],(bool) currentTree.isEB[1],(bool) currentTree.isEE[1],currentTree.sigmaIetaIeta[1]))) continue;

          hNPhotons[1]->Fill(currentTree.nPhotons,weight);
          if (currentTree.isEB[0]) iDetector=1;
          if (currentTree.isEE[0]) iDetector=2;

          hLeadEt[iDetector][1]->Fill(currentTree.et[0],weight);
          hLeadEta[1]->Fill(currentTree.eta[0],weight);
          hLeadR9[iDetector][1]->Fill(currentTree.r9[0],weight);

          if (currentTree.isEB[1]) iDetector=1;
          if (currentTree.isEE[1]) iDetector=2;

          hSubLeadEt[iDetector][1]->Fill(currentTree.et[1],weight);
          hSubLeadEta[1]->Fill(currentTree.eta[1],weight);
          hSubLeadR9[iDetector][1]->Fill(currentTree.r9[1],weight);

          // calculate invariant mass

          VLead = TLorentzVector( currentTree.momentumX[0],   currentTree.momentumY[0],  currentTree.momentumZ[0], currentTree.energy[0]);
          VSubLead = TLorentzVector( currentTree.momentumX[1],   currentTree.momentumY[1],  currentTree.momentumZ[1], currentTree.energy[1]);
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


          int HiggsInWhichDetector = 0;
          if (currentTree.isEB[0] && currentTree.isEB[1]) HiggsInWhichDetector=0;  // both photons in barrel 
          if ((currentTree.isEB[0] && currentTree.isEE[1]) || (currentTree.isEE[0] && currentTree.isEB[1])) HiggsInWhichDetector=1; // at least one photon in endcap

          // all photon categories together 
          h_mass_2gamma[1][HiggsInWhichDetector]->Fill(InvMass,weight);
          h_pt_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
          h_pz_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
          h_eta_2gamma[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
          h_CosThetaStar[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

          // define the photon categories
          int category=0 ;         
          if (currentTree.r9[0]>.93 && currentTree.r9[1]>.93 && !( currentTree.nTracks[0]==2 && currentTree.nTracks[1]==2 ) ) {
	    category=1;	  // golden photons 
	  } else if ( currentTree.r9[0]>.93 && currentTree.r9[1]>.93 &&  currentTree.nTracks[0]==2 && currentTree.nTracks[1]==2 ) {
	    category=2;	  // golden photons with tracks
	  } else if ( ((currentTree.r9[0]>.93 && currentTree.isConverted[0]==0) && (currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>0.0005 && currentTree.r9[1]<=.93 ) ) ||
	       ((currentTree.r9[1]>.93 && currentTree.isConverted[0]==0) && (currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005 && currentTree.r9[0]<=.93  ) ) ) {
	    category=3; // 1 golden and 1 converted 
	  } 
	
	  //	  if (((currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005) || (currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>0.0005)) 
	  //   &&!(currentTree.isConverted[0]==1 && currentTree.isConverted[1]==1)) {
	  // category=2; // 1 golden and 1 converted 
	  //} 

	  else if (currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005 && currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>.0005) {
	    category=4; 	    // two conversions 
	  }


          if ( category==1 ) {
            h_mass_2gamma_2gold[1][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2gold[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_2gold[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
          } else if ( category==2 ) {

            h_mass_2gamma_2goldConv[1][HiggsInWhichDetector]->Fill(InvMass,weight);

          } else if ( category==3 ) {
	    
            h_mass_2gamma_1conv[1][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_1conv[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_1conv[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_1conv[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_1conv[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
          } else if ( category==4 ) {

            h_mass_2gamma_2conv[1][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_2conv[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_2conv[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);
          } else {

            h_mass_2gamma_leftover[1][HiggsInWhichDetector]->Fill(InvMass,weight);
            h_pt_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Pt(),weight);
            h_pz_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Pz(),weight);
            h_eta_2gamma_leftover[1][HiggsInWhichDetector]->Fill(VSum.Eta(),weight);
            h_CosThetaStar_leftover[1][HiggsInWhichDetector]->Fill(cos_thetastar,weight);

	  }
      
        }

      }
    
      outfile->Write();
      outfile->Close();

    }
}
