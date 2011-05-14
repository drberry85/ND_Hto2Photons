#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCut.h"
#include "TString.h"

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

int main(int argc, char * input[]) {

  bool bar = false;
  bool data = false;
  bool dataweight = false;
  bool unweighted = false;

  //float globalWeight = 29.19;
  float globalWeight = 1000;
  float BranchingFraction = 0;

  int FirstFileNum = 0;
  
  vector<pair<string, float> > filesAndWeights;
  vector<pair<string, int> > filelist;
    
  TString InputArgs(input[1]);
  //cout << "Mass is: " << Mass << endl;

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
  
  // Load Signal
  if (InputArgs.Contains("Data")) {
    filelist.push_back(pair<string,int> ("Data.root",2));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc4/b/tkolberg/MPAntuples/Oct1_preselOriginalRefitMomentum/data.root",1));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc4/b/tkolberg/MPAntuples/Nov6_V00-00-13/Nov6_V00-00-13.root",1));
    data=true;
  }
  if (InputArgs.Contains("Yousidata")) {
    filelist.push_back(pair<string,int> ("YousiData.root",1));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/YousiData/MPA_Run2010B_Nov5_12831nb.root",1));
    data=true;
  }
  if (InputArgs.Contains("90GeV") || InputArgs.Contains("All")) {
    BranchingFraction = 0.000726;
    filelist.push_back(pair<string,int> ("HiggsAnalysis90GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon90.root",34.145*BranchingFraction/98996));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF90.root",1.7801*BranchingFraction/108813));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ90.root",3.1880*BranchingFraction/88000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
  }
  if (InputArgs.Contains("95GeV") || InputArgs.Contains("All")) {
    BranchingFraction = 0.00108;
    filelist.push_back(pair<string,int> ("HiggsAnalysis95GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon95.root",31.886*BranchingFraction/83986));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF95.root",1.6956*BranchingFraction/109579));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ95.root",2.9226*BranchingFraction/110000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
  }
  if (InputArgs.Contains("100GeV") || InputArgs.Contains("All")) {
    BranchingFraction = 0.00140;
    filelist.push_back(pair<string,int> ("HiggsAnalysis100GeV.root",2));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF100.root",1.5929*BranchingFraction/109826));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ100.root",1.9366*BranchingFraction/110000));
    cout << "Warning Weights Not Correct!!!!!" << endl;
    cout << "Warning no Gluon Fusion Samples!!!!!" << endl;
  }
  if (InputArgs.Contains("105GeV") || InputArgs.Contains("All")) {
    BranchingFraction = 0.001755;
    filelist.push_back(pair<string,int> ("HiggsAnalysis105GeV.root",2));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF105.root",1.5138*BranchingFraction/109835));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ105.root",1.6674*BranchingFraction/110000));
    cout << "Warning no Gluon Fusion Samples!!!!!" << endl;
  }
  if (InputArgs.Contains("110GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.001939;
    filelist.push_back(pair<string,int> ("HiggsAnalysis110GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon110.root",20.493*BranchingFraction/109994));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF110.root",1.4405*BranchingFraction/105974));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ110.root",1.4421*BranchingFraction/110000));
  }
  if (InputArgs.Contains("115GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.002101;
    filelist.push_back(pair<string,int> ("HiggsAnalysis115GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon115.root",18.735*BranchingFraction/109991));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF115.root",1.3712*BranchingFraction/109834));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ115.root",1.2524*BranchingFraction/110000));
  }
  if (InputArgs.Contains("120GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.002219;
    filelist.push_back(pair<string,int> ("HiggsAnalysis120GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon120.root",17.173*BranchingFraction/106151));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF120.root",1.3062*BranchingFraction/109848));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ120.root",1.0921*BranchingFraction/110000));
  }
  if (InputArgs.Contains("130GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.002240;
    filelist.push_back(pair<string,int> ("HiggsAnalysis130GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon130.root",14.579*BranchingFraction/109991));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF130.root",1.1866*BranchingFraction/109848));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ130.root",0.8395*BranchingFraction/110000));
  }
  if (InputArgs.Contains("140GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.001929;
    filelist.push_back(pair<string,int> ("HiggsAnalysis140GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon140.root",12.525*BranchingFraction/109991));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF140.root",1.0811*BranchingFraction/109842));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ140.root",0.6539*BranchingFraction/110000));
  }
  if (InputArgs.Contains("150GeV") || InputArgs.Contains("Signal") || InputArgs.Contains("All")) {
    BranchingFraction = 0.001363;
    filelist.push_back(pair<string,int> ("HiggsAnalysis150GeV.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon150.root",10.863*BranchingFraction/43500));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsVBF150.root",0.9868*BranchingFraction/50000));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsQQ150.root",0.5155*BranchingFraction/48000));
  }
  if (InputArgs.Contains("PhotonPlusJet") || InputArgs.Contains("Background") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("PhotonPlusJet.root",11));
    //filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet0to15.root",84200000.0/1057100));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet15to30.root",171700.0/1025840));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet30to50.root",16690.0/1025480));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet50to80.root",2722.0/1024608));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet80to120.root",447.2/1048215));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet120to170.root",84.17/1023361));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet170to300.root",22.64/1089000));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet300to470.root",1.493/1076926));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet470to800.root",0.1323/1093499));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet800to1400.root",0.003481/1092742));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet1400to1800.root",0.00001270/1097060));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_PhotonPlusJet1800toInf.root",0.0000002936/1091360));
  }
  if (InputArgs.Contains("EMEnriched") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("EMEnriched.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt20to30.root",236000000/(37169939/0.0104)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt30to80.root",59480000/(71845473/0.065)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_EMEnrichedpt80to170.root",900000/(8073559/0.155)));
  }
  if (InputArgs.Contains("Doubleemenriched") || InputArgs.Contains("Background") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("DoubleEMEnriched.root",1));
    //filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt10to20.root",20750000000/(31536145/0.0563)));
    //filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt20.root",293300000/(10912061/0.239)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt40.root",18700000/(21229315/0.00216)));
  }
  if (InputArgs.Contains("Reweighteddoubleemenriched") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("ReweightedDoubleEMEnriched.root",1));
    //filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt10to20.root",1.15*20750000000/(31536145/0.0563)));
    //filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt20.root",1.15*293300000/(10912061/0.239)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDDoubleEMEnrichedpt40.root",1.15*18700000/(21229315/0.00216)));
  }
  if (InputArgs.Contains("QCDBCtoE") || InputArgs.Contains("Background") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("QCDBCtoE.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt20to30.root",236000000/(2243439/0.00056)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt30to80.root",59480000/(1995502/0.00230)));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_QCDBCtoEpt80to170.root",900000/(1043390/0.0104)));
  }
  if (InputArgs.Contains("Born") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("Born.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt10to25.root",236.4/523270));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt25to250.root",22.37/536230));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBorn_Pt250toInf.root",0.008072/541900));
  }
  if (InputArgs.Contains("Box") || InputArgs.Contains("Background") || InputArgs.Contains("All")) {
    filelist.push_back(pair<string,int> ("Box.root",3));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt10to25.root",358.2/792710));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt25to250.root",12.37/768815));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Background/MPA_DiPhotonBox_Pt250toInf.root",0.000208/790685));
  }
  if (InputArgs.Contains("MPATest")) {
    filelist.push_back(pair<string,int> ("MPATest.root",1));
    filesAndWeights.push_back(pair<string,float> ("/data/ndpc2/c/HiggsGammaGamma/MPA/Signal/MPA_HiggsGluon130.root",1));
  }

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

    string outfilename = "";

    if (unweighted) {
      outfilename = "Unweighted";
      outfilename += itFilePair->first;
    } else if (dataweight) {
      outfilename = "Dataweight";
      outfilename += itFilePair->first;
    } else {
      outfilename = itFilePair->first;
    }
    
    TFile* outfile = new TFile(outfilename.c_str(),"RECREATE");
    outfile->cd();
    cout << outfilename << " created." << endl;

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
    TH2F* h2_convVtxRvsZBarrel_[2];

    TH1F * hLeadEtMarco;
    TH1F * hSubLeadEtMarco;
    TH1F * h_mass_Marco;

    TH1F * hLeadEtMarcoCat[4];
    TH1F * hSubLeadEtMarcoCat[4];
    TH1F * h_mass_MarcoCat[4];

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

    hLeadDzPV_[0][0]= new TH1F("leadPhoDZPVAll_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20.0, 20.0);
    hLeadDzPV_[1][0]= new TH1F("leadPhoDZPVAll_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20.0, 20.0);
    hLeadDzPV_[2][0]= new TH1F("leadPhoDZPVAll_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), all candidates", 100, -20.0, 20.0);

    hSubLeadDzPV_[0][0]= new TH1F("subleadPhoDZPVAll_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: all ecal", 100, -20.0, 20.0);
    hSubLeadDzPV_[1][0]= new TH1F("subleadPhoDZPVAll_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Barrel", 100, -20.0, 20.0);
    hSubLeadDzPV_[2][0]= new TH1F("subleadPhoDZPVAll_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), all candidates: Endcap", 100, -20.0, 20.0);

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

    hLeadDzPV_[0][1]= new TH1F("leadPhoDZPVSel_allEcal"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20.0, 20.0);
    hLeadDzPV_[1][1]= new TH1F("leadPhoDZPVSel_Barrel"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20.0, 20.0);
    hLeadDzPV_[2][1]= new TH1F("leadPhoDZPVSel_Endcap"," Leading photon #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20.0, 20.0);

    hSubLeadDzPV_[0][1]= new TH1F("subleadPhoDZPVSel_allEcal","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates:  all ECAL", 100, -20.0, 20.0);
    hSubLeadDzPV_[1][1]= new TH1F("subleadPhoDZPVSel_Barrel","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Barrel", 100, -20.0, 20.0);
    hSubLeadDzPV_[2][1]= new TH1F("subleadPhoDZPVSel_Endcap","Subleading #Deltaz_{Zconv - Ztrue} (cm), selected candidates: Endcap", 100, -20.0, 20.0);

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
    h_phi_conv[0][1] = new TH1D("h_phi_conv_Sel_AllECAL","#phi of Slected Photon Conversion All ECAL; #phi of Conversion", 64, -3.2, 3.2);
    
    h_phi_conv[1][0] = new TH1D("h_phi_conv_All_Barrel","#phi of Photon Conversion Barrel; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[1][1] = new TH1D("h_phi_conv_Sel_Barrel","#phi of Selected Photon Conversion Barrel; #phi of Conversion", 64, -3.2, 3.2);
    
    h_phi_conv[2][0] = new TH1D("h_phi_conv_All_Endcap","#phi of Photon Conversion Endcap; #phi of Conversion", 64, -3.2, 3.2);
    h_phi_conv[2][1] = new TH1D("h_phi_conv_Sel_Endcap","#phi of Selected Photon Conversion Endcap; #phi of Conversion", 64, -3.2, 3.2);
    
    h2_convVtxRvsZBarrel_[0] =   new TH2F("convVtxRvsZBarrelAll"," Photon  conversion vtx position all candidates Barrel",200, 0., 280., 200, 0., 80.);
    h2_convVtxRvsZBarrel_[1] =   new TH2F("convVtxRvsZBarrelSel"," Photon  conversion vtx position selected candidates Barrel",200, 0., 280., 200, 0., 80.);

    for (int itFile = FirstFileNum; itFile<itFilePair->second+FirstFileNum; itFile++) {

      string file = filesAndWeights[itFile].first;
      float weight = filesAndWeights[itFile].second * globalWeight;

      if (unweighted) weight=1;
      if (itFilePair->first=="Data.root" || itFilePair->first=="YousiData.root") weight=1;
      
      TFile * currentFile = new TFile(file.c_str());
      currentFile->cd();

      TTree * Analysis = (TTree *) currentFile->Get("NTuples/Analysis");
      
      //cout << "FirstFileNum is " << FirstFileNum << " and itFile is: " << itFile << endl;
      cout << "\nReading the tree in file " << file << endl;

      mpaReader currentTree(Analysis);
      
      Long64_t nentries = currentTree.fChain->GetEntries();

      outfile->cd();

      int pct = 0;
      
      for ( Long64_t i = 0; i < nentries; ++i ) {
        int iLeadDetector = 0;
        int iSubleadDetector = 0;

        TVector3 conversionVertex;

        if (i % (nentries/100) == 0 && bar) {
          if (pct == 100) pct = 99;
          cout << "\033[100m";
          cout << "\r[";
          cout << "\033[42m";
          for (int ctr = 0; ctr <= pct / 2; ++ctr) cout << "-";
          cout << "\b>";
          cout << "\033[100m";
          for (int ctr = pct / 2; ctr < 49; ++ctr) cout << " ";
          cout << "]\033[0m";
          cout << " " << ++pct << "%";
          cout << flush;
        }
        
        currentTree.GetEntry(i);
        
        if (currentTree.nPhotons<1) continue;

        if (convSel(currentTree.nTracks[0], currentTree.convVtxValid[0], currentTree.convVtxChi2Prob[0], currentTree.convDPhiTracksAtVtx[0], currentTree.convpairCotThetaSeparation[0], currentTree.convEoverP[0], currentTree.convVtxR[0])) {
          h_phi_conv[0][0]->Fill(currentTree.phi[0],weight);
          if (currentTree.isEB[0]) h_phi_conv[1][0]->Fill(currentTree.phi[0],weight);
          if (currentTree.isEE[0]) h_phi_conv[2][0]->Fill(currentTree.phi[0],weight);
          if (looseId(currentTree.pt[0], currentTree.ecalRecHitSumEtConeDR04[0], currentTree.hcalTowerSumEtConeDR04[0], currentTree.trkSumPtHollowConeDR04[0], (bool) currentTree.isEB[0], (bool) currentTree.isEE[0], currentTree.sigmaIetaIeta[0], currentTree.hadronicOverEm[0])) {
            h_phi_conv[0][1]->Fill(currentTree.phi[0],weight);
            if (currentTree.isEB[0]) h_phi_conv[1][1]->Fill(currentTree.phi[0],weight);
            if (currentTree.isEE[0]) h_phi_conv[2][1]->Fill(currentTree.phi[0],weight);
          }
        }

        if (currentTree.nPhotons<2) continue;
        if (convSel(currentTree.nTracks[1], currentTree.convVtxValid[1], currentTree.convVtxChi2Prob[1], currentTree.convDPhiTracksAtVtx[1], currentTree.convpairCotThetaSeparation[1], currentTree.convEoverP[1], currentTree.convVtxR[1])) {
          h_phi_conv[0][0]->Fill(currentTree.phi[1],weight);
          if (currentTree.isEB[1]) h_phi_conv[1][0]->Fill(currentTree.phi[1],weight);
          if (currentTree.isEE[1]) h_phi_conv[2][0]->Fill(currentTree.phi[1],weight);
          if (looseId(currentTree.pt[1], currentTree.ecalRecHitSumEtConeDR04[1], currentTree.hcalTowerSumEtConeDR04[1], currentTree.trkSumPtHollowConeDR04[1], (bool) currentTree.isEB[1], (bool) currentTree.isEE[1], currentTree.sigmaIetaIeta[1], currentTree.hadronicOverEm[1])) {
            h_phi_conv[0][1]->Fill(currentTree.phi[1],weight);
            if (currentTree.isEB[1]) h_phi_conv[1][1]->Fill(currentTree.phi[1],weight);
            if (currentTree.isEE[1]) h_phi_conv[2][1]->Fill(currentTree.phi[1],weight);
          }
        }
        
        //////////////// basic selection
        if (currentTree.pt[0]<20) continue;
        if (currentTree.pt[1]<20) continue;
        if (fabs(currentTree.scEta[0])>2.5 || fabs(currentTree.scEta[1])>2.5) continue;
        if (currentTree.isEBEEGap[0] || currentTree.isEBEEGap[1] ) continue;
        ////////////////////////////////////
        if (currentTree.isEB[0]) iLeadDetector=1;
        if (currentTree.isEE[0]) iLeadDetector=2;
        if (currentTree.isEB[1]) iSubleadDetector=1;
        if (currentTree.isEE[1]) iSubleadDetector=2;
        /////////////////////////
        int leadPhoCategory = photonCategory ( currentTree.hasPixelSeed[0], currentTree.r9[0],  currentTree.nTracks[0], currentTree.convVtxChi2Prob[0] , currentTree.pt[0]/currentTree.convPairMomentumPerp[0]);
        int subleadPhoCategory = photonCategory (currentTree.hasPixelSeed[1], currentTree.r9[1],  currentTree.nTracks[1], currentTree.convVtxChi2Prob[1] ,  currentTree.pt[1]/currentTree.convPairMomentumPerp[1]);
        int diPhoCategory = diPhotonCategory( leadPhoCategory, subleadPhoCategory );
        ////////////////////////////////////
        
        hLeadEt[0][0]->Fill(currentTree.et[0],weight);
        hLeadEta[0]->Fill(currentTree.eta[0],weight);
        hLeadPhi[0]->Fill(currentTree.phi[0],weight);
        hLeadR9[0][0]->Fill(currentTree.r9[0],weight);
        hLeadHoE[0][0]->Fill(currentTree.hadronicOverEm[0],weight);
        hLeadTrkPtSumSolid03[0][0]->Fill(currentTree.trkSumPtSolidConeDR03[0],weight);
        hLeadEcalPtSumSolid03[0][0]->Fill(currentTree.ecalRecHitSumEtConeDR03[0],weight);
        hLeadHcalPtSumSolid03[0][0]->Fill(currentTree.hcalTowerSumEtConeDR03[0],weight);
        hLeadSigmaIetaIeta[0][0]->Fill(currentTree.sigmaIetaIeta[0],weight);
        hLeadZPV_[0][0]->Fill(currentTree.vtxZ,weight);
        hLeadDzPV_[0][0]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        if (currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005 && (bool) currentTree.isEB[0]) {
          conversionVertex.SetXYZ(currentTree.convVtxX[0],currentTree.convVtxY[0],currentTree.convVtxZ[0]);
          h2_convVtxRvsZBarrel_[0]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }
        
        hNPhotons[0]->Fill(currentTree.nPhotons,weight);

        hLeadEt[iLeadDetector][0]->Fill(currentTree.et[0],weight);
        hLeadR9[iLeadDetector][0]->Fill(currentTree.r9[0],weight);
        hLeadHoE[iLeadDetector][0]->Fill(currentTree.hadronicOverEm[0],weight);
        hLeadTrkPtSumSolid03[iLeadDetector][0]->Fill(currentTree.trkSumPtSolidConeDR03[0],weight);
        hLeadEcalPtSumSolid03[iLeadDetector][0]->Fill(currentTree.ecalRecHitSumEtConeDR03[0],weight);
        hLeadHcalPtSumSolid03[iLeadDetector][0]->Fill(currentTree.hcalTowerSumEtConeDR03[0],weight);
        hLeadSigmaIetaIeta[iLeadDetector][0]->Fill(currentTree.sigmaIetaIeta[0],weight);
        hLeadZPV_[iLeadDetector][0]->Fill(currentTree.vtxZ,weight);
        hLeadDzPV_[iLeadDetector][0]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        
        hSubLeadEt[0][0]->Fill(currentTree.et[1],weight);
        hSubLeadEta[0]->Fill(currentTree.eta[1],weight);
        hSubLeadPhi[0]->Fill(currentTree.phi[1],weight);
        hSubLeadR9[0][0]->Fill(currentTree.r9[1],weight);
        hSubLeadHoE[0][0]->Fill(currentTree.hadronicOverEm[1],weight);
        hSubLeadTrkPtSumSolid03[0][0]->Fill(currentTree.trkSumPtSolidConeDR03[1],weight);
        hSubLeadEcalPtSumSolid03[0][0]->Fill(currentTree.ecalRecHitSumEtConeDR03[1],weight);
        hSubLeadHcalPtSumSolid03[0][0]->Fill(currentTree.hcalTowerSumEtConeDR03[1],weight);
        hSubLeadSigmaIetaIeta[0][0]->Fill(currentTree.sigmaIetaIeta[1],weight);
        hSubLeadZPV_[0][0]->Fill(currentTree.vtxZ,weight);
        hSubLeadDzPV_[0][0]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        if (currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>0.0005 && (bool) currentTree.isEB[1]) {
          conversionVertex.SetXYZ(currentTree.convVtxX[1],currentTree.convVtxY[1],currentTree.convVtxZ[1]);
          h2_convVtxRvsZBarrel_[0]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        /*if (currentTree.pt[1]>currentTree.pt[0]) {
          cout << "WARNING! - Tree is not pt sorted!!!!!!!!!" << endl;
          cout << "LeadPt is " << currentTree.pt[0] << endl;
          cout << "SubLeadPt is " << currentTree.pt[1] << endl;
        }*/

        hSubLeadEt[iSubleadDetector][0]->Fill(currentTree.et[1],weight);
        hSubLeadR9[iSubleadDetector][0]->Fill(currentTree.r9[1],weight);
        hSubLeadHoE[iSubleadDetector][0]->Fill(currentTree.hadronicOverEm[1],weight);
        hSubLeadTrkPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.trkSumPtSolidConeDR03[1],weight);
        hSubLeadEcalPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.ecalRecHitSumEtConeDR03[1],weight);
        hSubLeadHcalPtSumSolid03[iSubleadDetector][0]->Fill(currentTree.hcalTowerSumEtConeDR03[1],weight);
        hSubLeadSigmaIetaIeta[iSubleadDetector][0]->Fill(currentTree.sigmaIetaIeta[1],weight);
        hSubLeadZPV_[iSubleadDetector][0]->Fill(currentTree.vtxZ,weight);
        hSubLeadDzPV_[iSubleadDetector][0]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        
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

        if (currentTree.pt[0]>40 && currentTree.pt[1]>30 && InvMass>90 && InvMass<250
            && MarcosCut(currentTree.pt[0], currentTree.ecalRecHitSumEtConeDR04[0], currentTree.hcalTowerSumEtConeDR04[0], currentTree.trkSumPtHollowConeDR04[0], currentTree.hasPixelSeed[0], currentTree.isEB[0], currentTree.isEE[0], currentTree.sigmaIetaIeta[0], currentTree.hadronicOverEm[0])
            && MarcosCut(currentTree.pt[1], currentTree.ecalRecHitSumEtConeDR04[1], currentTree.hcalTowerSumEtConeDR04[1], currentTree.trkSumPtHollowConeDR04[1], currentTree.hasPixelSeed[1], currentTree.isEB[1], currentTree.isEE[1], currentTree.sigmaIetaIeta[1], currentTree.hadronicOverEm[1])
            ) {
          hLeadEtMarco->Fill(currentTree.pt[0]);
          hSubLeadEtMarco->Fill(currentTree.pt[1]);
          h_mass_Marco->Fill(InvMass);

          int MarcosCategory = MarcosCutCategory(currentTree.r9[0], currentTree.isEB[0], currentTree.isEE[0], currentTree.r9[1], currentTree.isEB[1], currentTree.isEE[1]);
          hLeadEtMarcoCat[MarcosCategory]->Fill(currentTree.pt[0]);
          hSubLeadEtMarcoCat[MarcosCategory]->Fill(currentTree.pt[1]);
          h_mass_MarcoCat[MarcosCategory]->Fill(InvMass);

        }
        
        /// di-photon system before event selection
        int HiggsInWhichDetector = 0;
        if (currentTree.isEB[0] && currentTree.isEB[1]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.isEB[0] && currentTree.isEE[1]) || (currentTree.isEE[0] && currentTree.isEB[1])) HiggsInWhichDetector=1; // at least one photon in endcap

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
        if (currentTree.pt[0]<40) continue; // leading photon
        if (currentTree.pt[1]<30) continue; // subleading photon
        //isolation 
        
        if (fabs(currentTree.scEta[0])>2.5 || fabs(currentTree.scEta[1])>2.5) continue;
        if (!(looseId(currentTree.pt[0],
		      currentTree.ecalRecHitSumEtConeDR04[0],
		      currentTree.hcalTowerSumEtConeDR04[0],
		      currentTree.trkSumPtHollowConeDR04[0],
		      (bool) currentTree.isEB[0],
		      (bool) currentTree.isEE[0],
		      currentTree.sigmaIetaIeta[0],
		      currentTree.hadronicOverEm[0]))) continue;

        if (!(looseId(currentTree.pt[1],
		      currentTree.ecalRecHitSumEtConeDR04[1],
		      currentTree.hcalTowerSumEtConeDR04[1],
		      currentTree.trkSumPtHollowConeDR04[1],
		      (bool) currentTree.isEB[1],
		      (bool) currentTree.isEE[1],
		      currentTree.sigmaIetaIeta[1],
		      currentTree.hadronicOverEm[1]))) continue;

        bool convsel1 = convSel(currentTree.nTracks[0],
				currentTree.convVtxValid[0] ,  
				currentTree.convVtxChi2Prob[0], 
				currentTree.convDPhiTracksAtVtx[0], 
				currentTree.convpairCotThetaSeparation[0], 
                currentTree.pt[0]/currentTree.convPairMomentumPerp[0],
                currentTree.convVtxR[0]);

        bool convsel2 = convSel(currentTree.nTracks[1],
				currentTree.convVtxValid[1] ,  
				currentTree.convVtxChi2Prob[1], 
				currentTree.convDPhiTracksAtVtx[1], 
				currentTree.convpairCotThetaSeparation[1], 
                currentTree.pt[1]/currentTree.convPairMomentumPerp[1],
                currentTree.convVtxR[1]);

        HiggsInWhichDetector = 0;
        if (currentTree.isEB[0] && currentTree.isEB[1]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.isEB[0] && currentTree.isEE[1]) || (currentTree.isEE[0] && currentTree.isEB[1])) HiggsInWhichDetector=1; // at least one photon in endcap

        hLeadEt[0][1]->Fill(currentTree.et[0],weight);
        hLeadEta[1]->Fill(currentTree.eta[0],weight);
        hLeadPhi[1]->Fill(currentTree.phi[0],weight);
        hLeadR9[0][1]->Fill(currentTree.r9[0],weight);
        hLeadHoE[0][1]->Fill(currentTree.hadronicOverEm[0],weight);
        hLeadTrkPtSumSolid03[0][1]->Fill(currentTree.trkSumPtSolidConeDR03[0],weight);
        hLeadEcalPtSumSolid03[0][1]->Fill(currentTree.ecalRecHitSumEtConeDR03[0],weight);
        hLeadHcalPtSumSolid03[0][1]->Fill(currentTree.hcalTowerSumEtConeDR03[0],weight);
        hLeadSigmaIetaIeta[0][1]->Fill(currentTree.sigmaIetaIeta[0],weight);
        hLeadZPV_[0][1]->Fill(currentTree.vtxZ,weight);
        hLeadDzPV_[0][1]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        if (currentTree.nTracks[0]==2 && currentTree.convVtxChi2Prob[0]>0.0005 && (bool) currentTree.isEB[0]) {
          conversionVertex.SetXYZ(currentTree.convVtxX[0],currentTree.convVtxY[0],currentTree.convVtxZ[0]);
          h2_convVtxRvsZBarrel_[1]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        hSubLeadEt[0][1]->Fill(currentTree.et[1],weight);
        hSubLeadEta[1]->Fill(currentTree.eta[1],weight);
        hSubLeadPhi[1]->Fill(currentTree.phi[1],weight);
        hSubLeadR9[0][1]->Fill(currentTree.r9[1],weight);
        hSubLeadHoE[0][1]->Fill(currentTree.hadronicOverEm[1],weight);
        hSubLeadTrkPtSumSolid03[0][1]->Fill(currentTree.trkSumPtSolidConeDR03[1],weight);
        hSubLeadEcalPtSumSolid03[0][1]->Fill(currentTree.ecalRecHitSumEtConeDR03[1],weight);
        hSubLeadHcalPtSumSolid03[0][1]->Fill(currentTree.hcalTowerSumEtConeDR03[1],weight);
        hSubLeadSigmaIetaIeta[0][1]->Fill(currentTree.sigmaIetaIeta[1],weight);
        hSubLeadZPV_[0][1]->Fill(currentTree.vtxZ,weight);
        hSubLeadDzPV_[0][1]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);
        if (currentTree.nTracks[1]==2 && currentTree.convVtxChi2Prob[1]>0.0005 && (bool) currentTree.isEB[1]) {
          conversionVertex.SetXYZ(currentTree.convVtxX[1],currentTree.convVtxY[1],currentTree.convVtxZ[1]);
          h2_convVtxRvsZBarrel_[1]->Fill(conversionVertex.z(),conversionVertex.Perp(),weight);
        }

        hNPhotons[1]->Fill(currentTree.nPhotons,weight);

        hLeadEt[iLeadDetector][1]->Fill(currentTree.et[0],weight);
        hLeadR9[iLeadDetector][1]->Fill(currentTree.r9[0],weight);
        hLeadHoE[iLeadDetector][1]->Fill(currentTree.hadronicOverEm[0],weight);
        hLeadTrkPtSumSolid03[iLeadDetector][1]->Fill(currentTree.trkSumPtSolidConeDR03[0],weight);
        hLeadEcalPtSumSolid03[iLeadDetector][1]->Fill(currentTree.ecalRecHitSumEtConeDR03[0],weight);
        hLeadHcalPtSumSolid03[iLeadDetector][1]->Fill(currentTree.hcalTowerSumEtConeDR03[0],weight);
        hLeadSigmaIetaIeta[iLeadDetector][1]->Fill(currentTree.sigmaIetaIeta[0],weight);
        hLeadZPV_[iLeadDetector][1]->Fill(currentTree.vtxZ,weight);
        hLeadDzPV_[iLeadDetector][1]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);

        hSubLeadEt[iSubleadDetector][1]->Fill(currentTree.et[1],weight);
        hSubLeadR9[iSubleadDetector][1]->Fill(currentTree.r9[1],weight);
        hSubLeadHoE[iSubleadDetector][1]->Fill(currentTree.hadronicOverEm[1],weight);
        hSubLeadTrkPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.trkSumPtSolidConeDR03[1],weight);
        hSubLeadEcalPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.ecalRecHitSumEtConeDR03[1],weight);
        hSubLeadHcalPtSumSolid03[iSubleadDetector][1]->Fill(currentTree.hcalTowerSumEtConeDR03[1],weight);
        hSubLeadSigmaIetaIeta[iSubleadDetector][1]->Fill(currentTree.sigmaIetaIeta[1],weight);
        hSubLeadZPV_[iSubleadDetector][1]->Fill(currentTree.vtxZ,weight);
        hSubLeadDzPV_[iSubleadDetector][1]->Fill(currentTree.vtxZ-currentTree.simVertexZ,weight);

        if ( convsel1 )
          h2_convVtxRvsZBarrel_[1]->Fill(currentTree.convVtxZ[0],sqrt(currentTree.convVtxX[0]*currentTree.convVtxX[0]+currentTree.convVtxY[0]*currentTree.convVtxY[0]));
        if ( convsel2 )
          h2_convVtxRvsZBarrel_[1]->Fill(currentTree.convVtxZ[1],sqrt(currentTree.convVtxX[1]*currentTree.convVtxX[1]+currentTree.convVtxY[1]*currentTree.convVtxY[1]));


        if (  leadPhoCategory==0) 
          h_lead_r9_cat0[iLeadDetector]->Fill (currentTree.r9[0],weight);
        if (  subleadPhoCategory==0) 
          h_sublead_r9_cat0[iSubleadDetector]->Fill (currentTree.r9[1],weight);


        if (  leadPhoCategory==1) 
          h_lead_r9_cat1[iLeadDetector]->Fill (currentTree.r9[0],weight);
        if (  subleadPhoCategory==1) 
          h_sublead_r9_cat1[iSubleadDetector]->Fill (currentTree.r9[1],weight);


        if (  leadPhoCategory==2) 
          h_lead_r9_cat2[iLeadDetector]->Fill (currentTree.r9[0],weight);
        if (  subleadPhoCategory==2) 
          h_sublead_r9_cat2[iSubleadDetector]->Fill (currentTree.r9[1],weight);

        if (  leadPhoCategory==3) 
          h_lead_r9_cat3[iLeadDetector]->Fill (currentTree.r9[0],weight);
        if (  subleadPhoCategory==3) 
          h_sublead_r9_cat3[iSubleadDetector]->Fill (currentTree.r9[1],weight);


        if (  leadPhoCategory==4) 
          h_lead_r9_cat4[iLeadDetector]->Fill (currentTree.r9[0],weight);
        if (  subleadPhoCategory==4) 
          h_sublead_r9_cat4[iSubleadDetector]->Fill (currentTree.r9[1],weight);
        
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

        HiggsInWhichDetector = 0;
        if (currentTree.isEB[0] && currentTree.isEB[1]) HiggsInWhichDetector=0;  // both photons in barrel 
        if ((currentTree.isEB[0] && currentTree.isEE[1]) || (currentTree.isEE[0] && currentTree.isEB[1])) HiggsInWhichDetector=1; // at least one photon in endcap

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
        if (!data && (bool) currentTree.isGenMatched[0] && (bool) currentTree.isGenMatched[1]) {

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
