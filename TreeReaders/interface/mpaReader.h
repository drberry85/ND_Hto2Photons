//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  1 15:30:00 2010 by ROOT version 5.22/00d
// from TTree Analysis/Analysis
// found on file: /data/ndpc1/b/tkolberg/MPAntuples/PhotonJet_Pt170_Spring10-START3X_V26_S09-v1_GEN-SIM-RECO/MultiPhotonAnalyzer.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class mpaReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           orbit;
   Int_t           bunchCrossing;
   Int_t           luminosityBlock;
   Int_t           timesec;
   Int_t           nTTBits;
   Int_t           kMaxTTBits;
   UChar_t         TTBit[64];   //[nTTBits]
   UChar_t         TTBit34;
   UChar_t         TTBit40;
   UChar_t         TTBit41;
   UChar_t         TTBit0;
   Int_t           kMaxL1Bits;
   Int_t           nL1Bits;
   UChar_t         L1Bit[128];   //[nL1Bits]
//    UChar_t         HLT_L1SingleEG5;
//    UChar_t         HLT_Photon10_L1R;
//    UChar_t         HLT_Photon15_L1R;
//    UChar_t         HLT_Photon20_L1R;
//    UChar_t         HLT_Photon20_Cleaned_L1R;
//    UChar_t         HLT_Photon30_Cleaned_L1R;
   Int_t           kMaxTrigFlag;
   Int_t           nHLTBits;
   UChar_t         HLTBit[100];   //[nHLTBits]
   Int_t           nHfTowersP;
   Int_t           nHfTowersN;
   Int_t           nVtxGood;
   Float_t         simVertexX;
   Float_t         simVertexY;
   Float_t         simVertexZ;
   Float_t         ptHat;
   Double_t        beamSpotX;
   Double_t        beamSpotY;
   Double_t        beamSpotZ;
//   UChar_t         vtxIsFake;
   Double_t        vtxX;
   Double_t        vtxY;
   Double_t        vtxZ;
   Double_t        vtxXError;
   Double_t        vtxYError;
   Double_t        vtxZError;
   Int_t           vtxNTrk;
   Int_t           vtxNTrkWeight05;
   Double_t        vtxChi2;
   Double_t        vtxNdof;
   Double_t        vtxNormChi2;
   Int_t           isMETEmpty;
   Double_t        metEt;
   Double_t        metPx;
   Double_t        metPy;
   Double_t        metPz;
   Double_t        metE_longitudinal;
   Int_t           metNCorrections;
   Float_t         metCorEx;
   Float_t         metCorEy;
   Float_t         metCorSumEt;
   Float_t         metUncorrectedPt;
   Float_t         metUncorrectedPhi;
   UChar_t         metIsCaloMET;
   UChar_t         metIsRecoMET;
//    Double_t        metMaxEtInEmTowers;
//    Double_t        metMaxEtInHadTowers;
//    Double_t        metEtFractionHadronic;
//    Double_t        metEmEtFraction;
//    Double_t        metHadEtInHB;
//    Double_t        metHadEtInHO;
//    Double_t        metHadEtInHE;
//    Double_t        metHadEtInHF;
//    Double_t        metSignificance;
//    Double_t        metCaloSETInpHF;
//    Double_t        metCaloSETInmHF;
//    Double_t        metCaloMETInpHF;
//    Double_t        metCaloMETInmHF;
//    Double_t        metCaloMETPhiInpHF;
//    Double_t        metCaloMETPhiInmHF;
   Int_t           nJets;
   Int_t           kMaxJets;
   Float_t         jetPt[20];   //[nJets]
   Float_t         jetE[20];   //[nJets]
   Float_t         jetP[20];   //[nJets]
   Float_t         jetEta[20];   //[nJets]
   Float_t         jetPhi[20];   //[nJets]
   Float_t         jetCharge[20];   //[nJets]
   Int_t           jetNtrk[20];   //[nJets]
   UChar_t         jetIsCaloJet[20];   //[nJets]
   UChar_t         jetIsPFJet[20];   //[nJets]
   UChar_t         jetIsBasicJet[20];   //[nJets]
   Float_t         jetMaxEInEmTowers[20];   //[nJets]
   Float_t         jetMaxEInHadTowers[20];   //[nJets]
   Float_t         jetEnergyFractionHadronic[20];   //[nJets]
   Float_t         jetEmEnergyFraction[20];   //[nJets]
   Float_t         jetHadEnergyInHB[20];   //[nJets]
   Float_t         jetHadEnergyInHO[20];   //[nJets]
   Float_t         jetHadEnergyInHE[20];   //[nJets]
   Float_t         jetHadEnergyInHF[20];   //[nJets]
   Float_t         jetEmEnergyInEB[20];   //[nJets]
   Float_t         jetEmEnergyInEE[20];   //[nJets]
   Float_t         jetEmEnergyInHF[20];   //[nJets]
   Float_t         jetTowersArea[20];   //[nJets]
   Int_t           jetN90[20];   //[nJets]
   Int_t           jetN60[20];   //[nJets]
   Float_t         jetFHPD[20];   //[nJets]
   Float_t         jetFRBX[20];   //[nJets]
   Float_t         jetN90Hits[20];   //[nJets]
   Int_t           nPhotons;
   Int_t           kMaxPhotons;
   Float_t         p[10];   //[nPhotons]
   Float_t         et[10];   //[nPhotons]
   Float_t         energy[10];   //[nPhotons]
   Float_t         momentumX[10];   //[nPhotons]
   Float_t         momentumY[10];   //[nPhotons]
   Float_t         momentumZ[10];   //[nPhotons]
   Float_t         pt[10];   //[nPhotons]
   Float_t         eta[10];   //[nPhotons]
   Float_t         phi[10];   //[nPhotons]
   Float_t         r9[10];   //[nPhotons]
   UChar_t         isEBGap[10];   //[nPhotons]
   UChar_t         isEEGap[10];   //[nPhotons]
   UChar_t         isEBEEGap[10];   //[nPhotons]
   UChar_t         isTransGap[10];   //[nPhotons]
   UChar_t         isEB[10];   //[nPhotons]
   UChar_t         isEE[10];   //[nPhotons]
   Float_t         rawEnergy[10];   //[nPhotons]
   Float_t         preshowerEnergy[10];   //[nPhotons]
   Int_t           numOfPreshClusters[10];   //[nPhotons]
   Float_t         ESRatio[10];   //[nPhotons]
   Int_t           clustersSize[10];   //[nPhotons]
   Int_t           scSize[10];   //[nPhotons]
   Float_t         phiWidth[10];   //[nPhotons]
   Float_t         etaWidth[10];   //[nPhotons]
   Float_t         scEta[10];   //[nPhotons]
   Float_t         scPhi[10];   //[nPhotons]
   Float_t         maxEnergyXtal[10];   //[nPhotons]
   Float_t         sigmaEtaEta[10];   //[nPhotons]
   Float_t         sigmaIetaIeta[10];   //[nPhotons]
   Float_t         r1x5[10];   //[nPhotons]
   Float_t         r2x5[10];   //[nPhotons]
   Float_t         e1x5[10];   //[nPhotons]
   Float_t         e2x5[10];   //[nPhotons]
   Float_t         seedTime[10];   //[nPhotons]
   Float_t         seedChi2[10];   //[nPhotons]
   Float_t         seedOutOfTimeChi2[10];   //[nPhotons]
   Int_t           seedRecoFlag[10];   //[nPhotons]
   Int_t           seedSeverity[10];   //[nPhotons]
   Float_t         eMax[10];   //[nPhotons]
   Float_t         e2nd[10];   //[nPhotons]
   Float_t         e2x2[10];   //[nPhotons]
   Float_t         e3x2[10];   //[nPhotons]
   Float_t         e3x3[10];   //[nPhotons]
   Float_t         e4x4[10];   //[nPhotons]
   Float_t         e5x5[10];   //[nPhotons]
   Float_t         e2x5Right[10];   //[nPhotons]
   Float_t         e2x5Left[10];   //[nPhotons]
   Float_t         e2x5Top[10];   //[nPhotons]
   Float_t         e2x5Bottom[10];   //[nPhotons]
   Float_t         eRight[10];   //[nPhotons]
   Float_t         eLeft[10];   //[nPhotons]
   Float_t         eTop[10];   //[nPhotons]
   Float_t         eBottom[10];   //[nPhotons]
   Float_t         covPhiPhi[10];   //[nPhotons]
   Float_t         covEtaPhi[10];   //[nPhotons]
   Float_t         covEtaEta[10];   //[nPhotons]
   Float_t         hadronicOverEm[10];   //[nPhotons]
   Float_t         hadronicDepth1OverEm[10];   //[nPhotons]
   Float_t         hadronicDepth2OverEm[10];   //[nPhotons]
   Float_t         trackIso[10];   //[nPhotons]
   Float_t         caloIso[10];   //[nPhotons]
   Float_t         ecalIso[10];   //[nPhotons]
   Float_t         hcalIso[10];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR04[10];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR04[10];   //[nPhotons]
   Int_t           nTrkSolidConeDR04[10];   //[nPhotons]
   Int_t           nTrkHollowConeDR04[10];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR03[10];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR03[10];   //[nPhotons]
   Int_t           nTrkSolidConeDR03[10];   //[nPhotons]
   Int_t           nTrkHollowConeDR03[10];   //[nPhotons]
   UChar_t         hasConversionTracks[10];   //[nPhotons]
   UChar_t         hasPixelSeed[10];   //[nPhotons]
   UChar_t         isLoose[10];   //[nPhotons]
   UChar_t         isTight[10];   //[nPhotons]
   Int_t           nTracks[10];   //[nPhotons]
   UChar_t         isConverted[10];   //[nPhotons]
   Float_t         convPairInvariantMass[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumEta[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumEta[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumEta[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumEta[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumPhi[10];   //[nPhotons]
   Float_t         convpairCotThetaSeparation[10];   //[nPhotons]
   Float_t         convPairMomentumMag[10];   //[nPhotons]
   Float_t         convPairMomentumPerp[10];   //[nPhotons]
   Float_t         convPairMomentumPhi[10];   //[nPhotons]
   Float_t         convPairMomentumEta[10];   //[nPhotons]
   Float_t         convPairMomentumX[10];   //[nPhotons]
   Float_t         convPairMomentumY[10];   //[nPhotons]
   Float_t         convPairMomentumZ[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumMag[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumPerp[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumPhi[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumEta[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumX[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumY[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumZ[10];   //[nPhotons]
   Float_t         convDistOfMinimumApproach[10];   //[nPhotons]
   Float_t         convDPhiTracksAtVtx[10];   //[nPhotons]
   Float_t         convDPhiTracksAtEcal[10];   //[nPhotons]
   Float_t         convDEtaTracksAtEcal[10];   //[nPhotons]
   UChar_t         convVtxValid[10];   //[nPhotons]
   Float_t         convVtxEta[10];   //[nPhotons]
   Float_t         convVtxPhi[10];   //[nPhotons]
   Float_t         convVtxR[10];   //[nPhotons]
   Float_t         convVtxX[10];   //[nPhotons]
   Float_t         convVtxY[10];   //[nPhotons]
   Float_t         convVtxZ[10];   //[nPhotons]
  //   Float_t         convVtxChi2[10];   //[nPhotons]
  //   Float_t         convVtxNdof[10];   //[nPhotons]
   Float_t         convMVALikelihood[10];   //[nPhotons]
   Float_t         convVtxChi2Prob[10];   //[nPhotons]
   Float_t         convEoverP[10];   //[nPhotons]
   Float_t         convzOfPrimaryVertexFromTracks[10];   //[nPhotons]
   UChar_t         isGenMatched[10];   //[nPhotons]
   Float_t         genMatchedPt[10];   //[nPhotons]
   Float_t         genMatchedEta[10];   //[nPhotons]
   Float_t         genMatchedPhi[10];   //[nPhotons]
   Int_t           genMomId[10];   //[nPhotons]
   Int_t           genGrandMomId[10];   //[nPhotons]
   Int_t           genNSiblings[10];   //[nPhotons]
   Float_t         genCalIsoDR03[10];   //[nPhotons]
   Float_t         genTrkIsoDR03[10];   //[nPhotons]
   Float_t         genCalIsoDR04[10];   //[nPhotons]
   Float_t         genTrkIsoDR04[10];   //[nPhotons]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_timesec;   //!
   TBranch        *b_nTTBits;   //!
   TBranch        *b_kMaxTTBits;   //!
   TBranch        *b_TTBit;   //!
   TBranch        *b_TTBit34;   //!
   TBranch        *b_TTBit40;   //!
   TBranch        *b_TTBit41;   //!
   TBranch        *b_TTBit0;   //!
   TBranch        *b_kMaxL1Bits;   //!
   TBranch        *b_nL1Bits;   //!
   TBranch        *b_L1Bit;   //!
//    TBranch        *b_HLT_L1SingleEG5;   //!
//    TBranch        *b_HLT_Photon10_L1R;   //!
//    TBranch        *b_HLT_Photon15_L1R;   //!
//    TBranch        *b_HLT_Photon20_L1R;   //!
//    TBranch        *b_HLT_Photon20_Cleaned_L1R; //!
//    TBranch        *b_HLT_Photon30_Cleaned_L1R; //!
   TBranch        *b_kMaxTrigFlag;   //!
   TBranch        *b_nHLTBits;   //!
   TBranch        *b_HLTBit;   //!
   TBranch        *b_nHfTowersP;   //!
   TBranch        *b_nHfTowersN;   //!
   TBranch        *b_nVtxGood;   //!
   TBranch        *b_simVertexX;   //!
   TBranch        *b_simVertexY;   //!
   TBranch        *b_simVertexZ;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
  //   TBranch        *b_vtxIsFake;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_vtxXError;   //!
   TBranch        *b_vtxYError;   //!
   TBranch        *b_vtxZError;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNTrkWeight05;   //!
   TBranch        *b_vtxChi2;   //!
   TBranch        *b_vtxNdof;   //!
   TBranch        *b_vtxNormChi2;   //!
   TBranch        *b_isMETEmpty;   //!
   TBranch        *b_metEt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metPz;   //!
   TBranch        *b_metE_longitudinal;   //!
   TBranch        *b_metNCorrections;   //!
   TBranch        *b_metCorEx;   //!
   TBranch        *b_metCorEy;   //!
   TBranch        *b_metCorSumEt;   //!
   TBranch        *b_metUncorrectedPt;   //!
   TBranch        *b_metUncorrectedPhi;   //!
   TBranch        *b_metIsCaloMET;   //!
   TBranch        *b_metIsRecoMET;   //!
//    TBranch        *b_metMaxEtInEmTowers;   //!
//    TBranch        *b_metMaxEtInHadTowers;   //!
//    TBranch        *b_metEtFractionHadronic;   //!
//    TBranch        *b_metEmEtFraction;   //!
//    TBranch        *b_metHadEtInHB;   //!
//    TBranch        *b_metHadEtInHO;   //!
//    TBranch        *b_metHadEtInHE;   //!
//    TBranch        *b_metHadEtInHF;   //!
//    TBranch        *b_metSignificance;   //!
//    TBranch        *b_metCaloSETInpHF;   //!
//    TBranch        *b_metCaloSETInmHF;   //!
//    TBranch        *b_metCaloMETInpHF;   //!
//    TBranch        *b_metCaloMETInmHF;   //!
//    TBranch        *b_metCaloMETPhiInpHF;   //!
//    TBranch        *b_metCaloMETPhiInmHF;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_kMaxJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetP;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetNtrk;   //!
   TBranch        *b_jetIsCaloJet;   //!
   TBranch        *b_jetIsPFJet;   //!
   TBranch        *b_jetIsBasicJet;   //!
   TBranch        *b_jetMaxEInEmTowers;   //!
   TBranch        *b_jetMaxEInHadTowers;   //!
   TBranch        *b_jetEnergyFractionHadronic;   //!
   TBranch        *b_jetEmEnergyFraction;   //!
   TBranch        *b_jetHadEnergyInHB;   //!
   TBranch        *b_jetHadEnergyInHO;   //!
   TBranch        *b_jetHadEnergyInHE;   //!
   TBranch        *b_jetHadEnergyInHF;   //!
   TBranch        *b_jetEmEnergyInEB;   //!
   TBranch        *b_jetEmEnergyInEE;   //!
   TBranch        *b_jetEmEnergyInHF;   //!
   TBranch        *b_jetTowersArea;   //!
   TBranch        *b_jetN90;   //!
   TBranch        *b_jetN60;   //!
   TBranch        *b_jetFHPD;   //!
   TBranch        *b_jetFRBX;   //!
   TBranch        *b_jetN90Hits;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_kMaxPhotons;   //!
   TBranch        *b_p;   //!
   TBranch        *b_et;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_momentumX;   //!
   TBranch        *b_momentumY;   //!
   TBranch        *b_momentumZ;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_r9;   //!
   TBranch        *b_isEBGap;   //!
   TBranch        *b_isEEGap;   //!
   TBranch        *b_isEBEEGap;   //!
   TBranch        *b_isTransGap;   //!
   TBranch        *b_isEB;   //!
   TBranch        *b_isEE;   //!
   TBranch        *b_rawEnergy;   //!
   TBranch        *b_preshowerEnergy;   //!
   TBranch        *b_numOfPreshClusters;   //!
   TBranch        *b_ESRatio;   //!
   TBranch        *b_clustersSize;   //!
   TBranch        *b_scSize;   //!
   TBranch        *b_phiWidth;   //!
   TBranch        *b_etaWidth;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_maxEnergyXtal;   //!
   TBranch        *b_sigmaEtaEta;   //!
   TBranch        *b_sigmaIetaIeta;   //!
   TBranch        *b_r1x5;   //!
   TBranch        *b_r2x5;   //!
   TBranch        *b_e1x5;   //!
   TBranch        *b_e2x5;   //!
   TBranch        *b_seedTime;   //!
   TBranch        *b_seedChi2;   //!
   TBranch        *b_seedOutOfTimeChi2;   //!
   TBranch        *b_seedRecoFlag;   //!
   TBranch        *b_seedSeverity;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_e2nd;   //!
   TBranch        *b_e2x2;   //!
   TBranch        *b_e3x2;   //!
   TBranch        *b_e3x3;   //!
   TBranch        *b_e4x4;   //!
   TBranch        *b_e5x5;   //!
   TBranch        *b_e2x5Right;   //!
   TBranch        *b_e2x5Left;   //!
   TBranch        *b_e2x5Top;   //!
   TBranch        *b_e2x5Bottom;   //!
   TBranch        *b_eRight;   //!
   TBranch        *b_eLeft;   //!
   TBranch        *b_eTop;   //!
   TBranch        *b_eBottom;   //!
   TBranch        *b_covPhiPhi;   //!
   TBranch        *b_covEtaPhi;   //!
   TBranch        *b_covEtaEta;   //!
   TBranch        *b_hadronicOverEm;   //!
   TBranch        *b_hadronicDepth1OverEm;   //!
   TBranch        *b_hadronicDepth2OverEm;   //!
   TBranch        *b_trackIso;   //!
   TBranch        *b_caloIso;   //!
   TBranch        *b_ecalIso;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_trkSumPtSolidConeDR04;   //!
   TBranch        *b_trkSumPtHollowConeDR04;   //!
   TBranch        *b_nTrkSolidConeDR04;   //!
   TBranch        *b_nTrkHollowConeDR04;   //!
   TBranch        *b_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_trkSumPtSolidConeDR03;   //!
   TBranch        *b_trkSumPtHollowConeDR03;   //!
   TBranch        *b_nTrkSolidConeDR03;   //!
   TBranch        *b_nTrkHollowConeDR03;   //!
   TBranch        *b_hasConversionTracks;   //!
   TBranch        *b_hasPixelSeed;   //!
   TBranch        *b_isLoose;   //!
   TBranch        *b_isTight;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_isConverted;   //!
   TBranch        *b_convPairInvariantMass;   //!
   TBranch        *b_convTrack1InnerMomentumPerp;   //!
   TBranch        *b_convTrack1InnerMomentumEta;   //!
   TBranch        *b_convTrack1InnerMomentumPhi;   //!
   TBranch        *b_convTrack1OuterMomentumPerp;   //!
   TBranch        *b_convTrack1OuterMomentumEta;   //!
   TBranch        *b_convTrack1OuterMomentumPhi;   //!
   TBranch        *b_convTrack2InnerMomentumPerp;   //!
   TBranch        *b_convTrack2InnerMomentumEta;   //!
   TBranch        *b_convTrack2InnerMomentumPhi;   //!
   TBranch        *b_convTrack2OuterMomentumPerp;   //!
   TBranch        *b_convTrack2OuterMomentumEta;   //!
   TBranch        *b_convTrack2OuterMomentumPhi;   //!
   TBranch        *b_convpairCotThetaSeparation;   //!
   TBranch        *b_convPairMomentumMag;   //!
   TBranch        *b_convPairMomentumPerp;   //!
   TBranch        *b_convPairMomentumPhi;   //!
   TBranch        *b_convPairMomentumEta;   //!
   TBranch        *b_convPairMomentumX;   //!
   TBranch        *b_convPairMomentumY;   //!
   TBranch        *b_convPairMomentumZ;   //!
   TBranch        *b_convPairRefittedMomentumMag;   //!
   TBranch        *b_convPairRefittedMomentumPerp;   //!
   TBranch        *b_convPairRefittedMomentumPhi;   //!
   TBranch        *b_convPairRefittedMomentumEta;   //!
   TBranch        *b_convPairRefittedMomentumX;   //!
   TBranch        *b_convPairRefittedMomentumY;   //!
   TBranch        *b_convPairRefittedMomentumZ;   //!
   TBranch        *b_convDistOfMinimumApproach;   //!
   TBranch        *b_convDPhiTracksAtVtx;   //!
   TBranch        *b_convDPhiTracksAtEcal;   //!
   TBranch        *b_convDEtaTracksAtEcal;   //!
   TBranch        *b_convVtxValid;   //!
   TBranch        *b_convVtxEta;   //!
   TBranch        *b_convVtxPhi;   //!
   TBranch        *b_convVtxR;   //!
   TBranch        *b_convVtxX;   //!
   TBranch        *b_convVtxY;   //!
   TBranch        *b_convVtxZ;   //!
  //   TBranch        *b_convVtxChi2;   //!
  //   TBranch        *b_convVtxNdof;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convVtxChi2Prob;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convzOfPrimaryVertexFromTracks;   //!
   TBranch        *b_isGenMatched;   //!
   TBranch        *b_genMatchedPt;   //!
   TBranch        *b_genMatchedEta;   //!
   TBranch        *b_genMatchedPhi;   //!
   TBranch        *b_genMomId;   //!
   TBranch        *b_genGrandMomId;   //!
   TBranch        *b_genNSiblings;   //!
   TBranch        *b_genCalIsoDR03;   //!
   TBranch        *b_genTrkIsoDR03;   //!
   TBranch        *b_genCalIsoDR04;   //!
   TBranch        *b_genTrkIsoDR04;   //!

   mpaReader(TTree *tree=0);
   virtual ~mpaReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

mpaReader::mpaReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/ndpc1/b/tkolberg/MPAntuples/PhotonJet_Pt170_Spring10-START3X_V26_S09-v1_GEN-SIM-RECO/MultiPhotonAnalyzer.root");
      if (!f) {
         f = new TFile("/data/ndpc1/b/tkolberg/MPAntuples/PhotonJet_Pt170_Spring10-START3X_V26_S09-v1_GEN-SIM-RECO/MultiPhotonAnalyzer.root");
      }
      tree = (TTree*)gDirectory->Get("Analysis");

   }
   Init(tree);
}

mpaReader::~mpaReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mpaReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mpaReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mpaReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("timesec", &timesec, &b_timesec);
   fChain->SetBranchAddress("nTTBits", &nTTBits, &b_nTTBits);
   fChain->SetBranchAddress("kMaxTTBits", &kMaxTTBits, &b_kMaxTTBits);
   fChain->SetBranchAddress("TTBit", TTBit, &b_TTBit);
   fChain->SetBranchAddress("TTBit34", &TTBit34, &b_TTBit34);
   fChain->SetBranchAddress("TTBit40", &TTBit40, &b_TTBit40);
   fChain->SetBranchAddress("TTBit41", &TTBit41, &b_TTBit41);
   fChain->SetBranchAddress("TTBit0", &TTBit0, &b_TTBit0);
   fChain->SetBranchAddress("kMaxL1Bits", &kMaxL1Bits, &b_kMaxL1Bits);
   fChain->SetBranchAddress("nL1Bits", &nL1Bits, &b_nL1Bits);
   fChain->SetBranchAddress("L1Bit", L1Bit, &b_L1Bit);
//    fChain->SetBranchAddress("HLT_L1SingleEG5", &HLT_L1SingleEG5, &b_HLT_L1SingleEG5);
//    fChain->SetBranchAddress("HLT_Photon10_L1R", &HLT_Photon10_L1R, &b_HLT_Photon10_L1R);
//    fChain->SetBranchAddress("HLT_Photon15_L1R", &HLT_Photon15_L1R, &b_HLT_Photon15_L1R);
//    fChain->SetBranchAddress("HLT_Photon20_L1R", &HLT_Photon20_L1R, &b_HLT_Photon20_L1R);
//    fChain->SetBranchAddress("HLT_Photon20_Cleaned_L1R", &HLT_Photon20_Cleaned_L1R, &b_HLT_Photon20_Cleaned_L1R);
//    fChain->SetBranchAddress("HLT_Photon30_Cleaned_L1R", &HLT_Photon30_Cleaned_L1R, &b_HLT_Photon30_Cleaned_L1R);
   fChain->SetBranchAddress("kMaxTrigFlag", &kMaxTrigFlag, &b_kMaxTrigFlag);
   fChain->SetBranchAddress("nHLTBits", &nHLTBits, &b_nHLTBits);
   fChain->SetBranchAddress("HLTBit", HLTBit, &b_HLTBit);
   fChain->SetBranchAddress("nHfTowersP", &nHfTowersP, &b_nHfTowersP);
   fChain->SetBranchAddress("nHfTowersN", &nHfTowersN, &b_nHfTowersN);
   fChain->SetBranchAddress("nVtxGood", &nVtxGood, &b_nVtxGood);
   fChain->SetBranchAddress("simVertexX", &simVertexX, &b_simVertexX);
   fChain->SetBranchAddress("simVertexY", &simVertexY, &b_simVertexY);
   fChain->SetBranchAddress("simVertexZ", &simVertexZ, &b_simVertexZ);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   //   fChain->SetBranchAddress("vtxIsFake", &vtxIsFake, &b_vtxIsFake);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("vtxXError", &vtxXError, &b_vtxXError);
   fChain->SetBranchAddress("vtxYError", &vtxYError, &b_vtxYError);
   fChain->SetBranchAddress("vtxZError", &vtxZError, &b_vtxZError);
   fChain->SetBranchAddress("vtxNTrk", &vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNTrkWeight05", &vtxNTrkWeight05, &b_vtxNTrkWeight05);
   fChain->SetBranchAddress("vtxChi2", &vtxChi2, &b_vtxChi2);
   fChain->SetBranchAddress("vtxNdof", &vtxNdof, &b_vtxNdof);
   fChain->SetBranchAddress("vtxNormChi2", &vtxNormChi2, &b_vtxNormChi2);
   fChain->SetBranchAddress("isMETEmpty", &isMETEmpty, &b_isMETEmpty);
   fChain->SetBranchAddress("metEt", &metEt, &b_metEt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metPz", &metPz, &b_metPz);
   fChain->SetBranchAddress("metE_longitudinal", &metE_longitudinal, &b_metE_longitudinal);
   fChain->SetBranchAddress("metNCorrections", &metNCorrections, &b_metNCorrections);
   fChain->SetBranchAddress("metCorEx", &metCorEx, &b_metCorEx);
   fChain->SetBranchAddress("metCorEy", &metCorEy, &b_metCorEy);
   fChain->SetBranchAddress("metCorSumEt", &metCorSumEt, &b_metCorSumEt);
   fChain->SetBranchAddress("metUncorrectedPt", &metUncorrectedPt, &b_metUncorrectedPt);
   fChain->SetBranchAddress("metUncorrectedPhi", &metUncorrectedPhi, &b_metUncorrectedPhi);
   fChain->SetBranchAddress("metIsCaloMET", &metIsCaloMET, &b_metIsCaloMET);
   fChain->SetBranchAddress("metIsRecoMET", &metIsRecoMET, &b_metIsRecoMET);
//    fChain->SetBranchAddress("metMaxEtInEmTowers", &metMaxEtInEmTowers, &b_metMaxEtInEmTowers);
//    fChain->SetBranchAddress("metMaxEtInHadTowers", &metMaxEtInHadTowers, &b_metMaxEtInHadTowers);
//    fChain->SetBranchAddress("metEtFractionHadronic", &metEtFractionHadronic, &b_metEtFractionHadronic);
//    fChain->SetBranchAddress("metEmEtFraction", &metEmEtFraction, &b_metEmEtFraction);
//    fChain->SetBranchAddress("metHadEtInHB", &metHadEtInHB, &b_metHadEtInHB);
//    fChain->SetBranchAddress("metHadEtInHO", &metHadEtInHO, &b_metHadEtInHO);
//    fChain->SetBranchAddress("metHadEtInHE", &metHadEtInHE, &b_metHadEtInHE);
//    fChain->SetBranchAddress("metHadEtInHF", &metHadEtInHF, &b_metHadEtInHF);
//    fChain->SetBranchAddress("metSignificance", &metSignificance, &b_metSignificance);
//    fChain->SetBranchAddress("metCaloSETInpHF", &metCaloSETInpHF, &b_metCaloSETInpHF);
//    fChain->SetBranchAddress("metCaloSETInmHF", &metCaloSETInmHF, &b_metCaloSETInmHF);
//    fChain->SetBranchAddress("metCaloMETInpHF", &metCaloMETInpHF, &b_metCaloMETInpHF);
//    fChain->SetBranchAddress("metCaloMETInmHF", &metCaloMETInmHF, &b_metCaloMETInmHF);
//    fChain->SetBranchAddress("metCaloMETPhiInpHF", &metCaloMETPhiInpHF, &b_metCaloMETPhiInpHF);
//    fChain->SetBranchAddress("metCaloMETPhiInmHF", &metCaloMETPhiInmHF, &b_metCaloMETPhiInmHF);
//   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
//   fChain->SetBranchAddress("kMaxJets", &kMaxJets, &b_kMaxJets);
//   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
//   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
//   fChain->SetBranchAddress("jetP", jetP, &b_jetP);
//   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
//   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
//   fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge);
//   fChain->SetBranchAddress("jetNtrk", jetNtrk, &b_jetNtrk);
//   fChain->SetBranchAddress("jetIsCaloJet", jetIsCaloJet, &b_jetIsCaloJet);
//   fChain->SetBranchAddress("jetIsPFJet", jetIsPFJet, &b_jetIsPFJet);
//   fChain->SetBranchAddress("jetIsBasicJet", jetIsBasicJet, &b_jetIsBasicJet);
//   fChain->SetBranchAddress("jetMaxEInEmTowers", jetMaxEInEmTowers, &b_jetMaxEInEmTowers);
//   fChain->SetBranchAddress("jetMaxEInHadTowers", jetMaxEInHadTowers, &b_jetMaxEInHadTowers);
//   fChain->SetBranchAddress("jetEnergyFractionHadronic", jetEnergyFractionHadronic, &b_jetEnergyFractionHadronic);
//   fChain->SetBranchAddress("jetEmEnergyFraction", jetEmEnergyFraction, &b_jetEmEnergyFraction);
//   fChain->SetBranchAddress("jetHadEnergyInHB", jetHadEnergyInHB, &b_jetHadEnergyInHB);
//   fChain->SetBranchAddress("jetHadEnergyInHO", jetHadEnergyInHO, &b_jetHadEnergyInHO);
//   fChain->SetBranchAddress("jetHadEnergyInHE", jetHadEnergyInHE, &b_jetHadEnergyInHE);
//   fChain->SetBranchAddress("jetHadEnergyInHF", jetHadEnergyInHF, &b_jetHadEnergyInHF);
//   fChain->SetBranchAddress("jetEmEnergyInEB", jetEmEnergyInEB, &b_jetEmEnergyInEB);
//   fChain->SetBranchAddress("jetEmEnergyInEE", jetEmEnergyInEE, &b_jetEmEnergyInEE);
//   fChain->SetBranchAddress("jetEmEnergyInHF", jetEmEnergyInHF, &b_jetEmEnergyInHF);
//   fChain->SetBranchAddress("jetTowersArea", jetTowersArea, &b_jetTowersArea);
//   fChain->SetBranchAddress("jetN90", jetN90, &b_jetN90);
//   fChain->SetBranchAddress("jetN60", jetN60, &b_jetN60);
//   fChain->SetBranchAddress("jetFHPD", jetFHPD, &b_jetFHPD);
//   fChain->SetBranchAddress("jetFRBX", jetFRBX, &b_jetFRBX);
//   fChain->SetBranchAddress("jetN90Hits", jetN90Hits, &b_jetN90Hits);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("kMaxPhotons", &kMaxPhotons, &b_kMaxPhotons);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("et", et, &b_et);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("momentumX", momentumX, &b_momentumX);
   fChain->SetBranchAddress("momentumY", momentumY, &b_momentumY);
   fChain->SetBranchAddress("momentumZ", momentumZ, &b_momentumZ);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("r9", r9, &b_r9);
   fChain->SetBranchAddress("isEBGap", isEBGap, &b_isEBGap);
   fChain->SetBranchAddress("isEEGap", isEEGap, &b_isEEGap);
   fChain->SetBranchAddress("isEBEEGap", isEBEEGap, &b_isEBEEGap);
   fChain->SetBranchAddress("isTransGap", isTransGap, &b_isTransGap);
   fChain->SetBranchAddress("isEB", isEB, &b_isEB);
   fChain->SetBranchAddress("isEE", isEE, &b_isEE);
   fChain->SetBranchAddress("rawEnergy", rawEnergy, &b_rawEnergy);
   fChain->SetBranchAddress("preshowerEnergy", preshowerEnergy, &b_preshowerEnergy);
   fChain->SetBranchAddress("numOfPreshClusters", numOfPreshClusters, &b_numOfPreshClusters);
   fChain->SetBranchAddress("ESRatio", ESRatio, &b_ESRatio);
   fChain->SetBranchAddress("clustersSize", clustersSize, &b_clustersSize);
   fChain->SetBranchAddress("scSize", scSize, &b_scSize);
   fChain->SetBranchAddress("phiWidth", phiWidth, &b_phiWidth);
   fChain->SetBranchAddress("etaWidth", etaWidth, &b_etaWidth);
   fChain->SetBranchAddress("scEta", scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", scPhi, &b_scPhi);
   fChain->SetBranchAddress("maxEnergyXtal", maxEnergyXtal, &b_maxEnergyXtal);
   fChain->SetBranchAddress("sigmaEtaEta", sigmaEtaEta, &b_sigmaEtaEta);
   fChain->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta, &b_sigmaIetaIeta);
   fChain->SetBranchAddress("r1x5", r1x5, &b_r1x5);
   fChain->SetBranchAddress("r2x5", r2x5, &b_r2x5);
   fChain->SetBranchAddress("e1x5", e1x5, &b_e1x5);
   fChain->SetBranchAddress("e2x5", e2x5, &b_e2x5);
   fChain->SetBranchAddress("seedTime", seedTime, &b_seedTime);
   fChain->SetBranchAddress("seedChi2", seedChi2, &b_seedChi2);
   fChain->SetBranchAddress("seedOutOfTimeChi2", seedOutOfTimeChi2, &b_seedOutOfTimeChi2);
   fChain->SetBranchAddress("seedRecoFlag", seedRecoFlag, &b_seedRecoFlag);
   fChain->SetBranchAddress("seedSeverity", seedSeverity, &b_seedSeverity);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("e2nd", e2nd, &b_e2nd);
   fChain->SetBranchAddress("e2x2", e2x2, &b_e2x2);
   fChain->SetBranchAddress("e3x2", e3x2, &b_e3x2);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e4x4", e4x4, &b_e4x4);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("e2x5Right", e2x5Right, &b_e2x5Right);
   fChain->SetBranchAddress("e2x5Left", e2x5Left, &b_e2x5Left);
   fChain->SetBranchAddress("e2x5Top", e2x5Top, &b_e2x5Top);
   fChain->SetBranchAddress("e2x5Bottom", e2x5Bottom, &b_e2x5Bottom);
   fChain->SetBranchAddress("eRight", eRight, &b_eRight);
   fChain->SetBranchAddress("eLeft", eLeft, &b_eLeft);
   fChain->SetBranchAddress("eTop", eTop, &b_eTop);
   fChain->SetBranchAddress("eBottom", eBottom, &b_eBottom);
   fChain->SetBranchAddress("covPhiPhi", covPhiPhi, &b_covPhiPhi);
   fChain->SetBranchAddress("covEtaPhi", covEtaPhi, &b_covEtaPhi);
   fChain->SetBranchAddress("covEtaEta", covEtaEta, &b_covEtaEta);
   fChain->SetBranchAddress("hadronicOverEm", hadronicOverEm, &b_hadronicOverEm);
   fChain->SetBranchAddress("hadronicDepth1OverEm", hadronicDepth1OverEm, &b_hadronicDepth1OverEm);
   fChain->SetBranchAddress("hadronicDepth2OverEm", hadronicDepth2OverEm, &b_hadronicDepth2OverEm);
   fChain->SetBranchAddress("trackIso", trackIso, &b_trackIso);
   fChain->SetBranchAddress("caloIso", caloIso, &b_caloIso);
   fChain->SetBranchAddress("ecalIso", ecalIso, &b_ecalIso);
   fChain->SetBranchAddress("hcalIso", hcalIso, &b_hcalIso);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR04", ecalRecHitSumEtConeDR04, &b_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR04", hcalTowerSumEtConeDR04, &b_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR04", hcalDepth1TowerSumEtConeDR04, &b_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR04", hcalDepth2TowerSumEtConeDR04, &b_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("trkSumPtSolidConeDR04", trkSumPtSolidConeDR04, &b_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("trkSumPtHollowConeDR04", trkSumPtHollowConeDR04, &b_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("nTrkSolidConeDR04", nTrkSolidConeDR04, &b_nTrkSolidConeDR04);
   fChain->SetBranchAddress("nTrkHollowConeDR04", nTrkHollowConeDR04, &b_nTrkHollowConeDR04);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR03", ecalRecHitSumEtConeDR03, &b_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR03", hcalTowerSumEtConeDR03, &b_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR03", hcalDepth1TowerSumEtConeDR03, &b_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR03", hcalDepth2TowerSumEtConeDR03, &b_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("trkSumPtSolidConeDR03", trkSumPtSolidConeDR03, &b_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("trkSumPtHollowConeDR03", trkSumPtHollowConeDR03, &b_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("nTrkSolidConeDR03", nTrkSolidConeDR03, &b_nTrkSolidConeDR03);
   fChain->SetBranchAddress("nTrkHollowConeDR03", nTrkHollowConeDR03, &b_nTrkHollowConeDR03);
   fChain->SetBranchAddress("hasConversionTracks", hasConversionTracks, &b_hasConversionTracks);
   fChain->SetBranchAddress("hasPixelSeed", hasPixelSeed, &b_hasPixelSeed);
   fChain->SetBranchAddress("isLoose", isLoose, &b_isLoose);
   fChain->SetBranchAddress("isTight", isTight, &b_isTight);
   fChain->SetBranchAddress("nTracks", nTracks, &b_nTracks);
   fChain->SetBranchAddress("isConverted", isConverted, &b_isConverted);
   fChain->SetBranchAddress("convPairInvariantMass", convPairInvariantMass, &b_convPairInvariantMass);
   fChain->SetBranchAddress("convTrack1InnerMomentumPerp", convTrack1InnerMomentumPerp, &b_convTrack1InnerMomentumPerp);
   fChain->SetBranchAddress("convTrack1InnerMomentumEta", convTrack1InnerMomentumEta, &b_convTrack1InnerMomentumEta);
   fChain->SetBranchAddress("convTrack1InnerMomentumPhi", convTrack1InnerMomentumPhi, &b_convTrack1InnerMomentumPhi);
   fChain->SetBranchAddress("convTrack1OuterMomentumPerp", convTrack1OuterMomentumPerp, &b_convTrack1OuterMomentumPerp);
   fChain->SetBranchAddress("convTrack1OuterMomentumEta", convTrack1OuterMomentumEta, &b_convTrack1OuterMomentumEta);
   fChain->SetBranchAddress("convTrack1OuterMomentumPhi", convTrack1OuterMomentumPhi, &b_convTrack1OuterMomentumPhi);
   fChain->SetBranchAddress("convTrack2InnerMomentumPerp", convTrack2InnerMomentumPerp, &b_convTrack2InnerMomentumPerp);
   fChain->SetBranchAddress("convTrack2InnerMomentumEta", convTrack2InnerMomentumEta, &b_convTrack2InnerMomentumEta);
   fChain->SetBranchAddress("convTrack2InnerMomentumPhi", convTrack2InnerMomentumPhi, &b_convTrack2InnerMomentumPhi);
   fChain->SetBranchAddress("convTrack2OuterMomentumPerp", convTrack2OuterMomentumPerp, &b_convTrack2OuterMomentumPerp);
   fChain->SetBranchAddress("convTrack2OuterMomentumEta", convTrack2OuterMomentumEta, &b_convTrack2OuterMomentumEta);
   fChain->SetBranchAddress("convTrack2OuterMomentumPhi", convTrack2OuterMomentumPhi, &b_convTrack2OuterMomentumPhi);
   fChain->SetBranchAddress("convpairCotThetaSeparation", convpairCotThetaSeparation, &b_convpairCotThetaSeparation);
   fChain->SetBranchAddress("convPairMomentumMag", convPairMomentumMag, &b_convPairMomentumMag);
   fChain->SetBranchAddress("convPairMomentumPerp", convPairMomentumPerp, &b_convPairMomentumPerp);
   fChain->SetBranchAddress("convPairMomentumPhi", convPairMomentumPhi, &b_convPairMomentumPhi);
   fChain->SetBranchAddress("convPairMomentumEta", convPairMomentumEta, &b_convPairMomentumEta);
   fChain->SetBranchAddress("convPairMomentumX", convPairMomentumX, &b_convPairMomentumX);
   fChain->SetBranchAddress("convPairMomentumY", convPairMomentumY, &b_convPairMomentumY);
   fChain->SetBranchAddress("convPairMomentumZ", convPairMomentumZ, &b_convPairMomentumZ);
   fChain->SetBranchAddress("convPairRefittedMomentumMag", convPairRefittedMomentumMag, &b_convPairRefittedMomentumMag);
   fChain->SetBranchAddress("convPairRefittedMomentumPerp", convPairRefittedMomentumPerp, &b_convPairRefittedMomentumPerp);
   fChain->SetBranchAddress("convPairRefittedMomentumPhi", convPairRefittedMomentumPhi, &b_convPairRefittedMomentumPhi);
   fChain->SetBranchAddress("convPairRefittedMomentumEta", convPairRefittedMomentumEta, &b_convPairRefittedMomentumEta);
   fChain->SetBranchAddress("convPairRefittedMomentumX", convPairRefittedMomentumX, &b_convPairRefittedMomentumX);
   fChain->SetBranchAddress("convPairRefittedMomentumY", convPairRefittedMomentumY, &b_convPairRefittedMomentumY);
   fChain->SetBranchAddress("convPairRefittedMomentumZ", convPairRefittedMomentumZ, &b_convPairRefittedMomentumZ);
   fChain->SetBranchAddress("convDistOfMinimumApproach", convDistOfMinimumApproach, &b_convDistOfMinimumApproach);
   fChain->SetBranchAddress("convDPhiTracksAtVtx", convDPhiTracksAtVtx, &b_convDPhiTracksAtVtx);
   fChain->SetBranchAddress("convDPhiTracksAtEcal", convDPhiTracksAtEcal, &b_convDPhiTracksAtEcal);
   fChain->SetBranchAddress("convDEtaTracksAtEcal", convDEtaTracksAtEcal, &b_convDEtaTracksAtEcal);
   fChain->SetBranchAddress("convVtxValid", convVtxValid, &b_convVtxValid);
   fChain->SetBranchAddress("convVtxEta", convVtxEta, &b_convVtxEta);
   fChain->SetBranchAddress("convVtxPhi", convVtxPhi, &b_convVtxPhi);
   fChain->SetBranchAddress("convVtxR", convVtxR, &b_convVtxR);
   fChain->SetBranchAddress("convVtxX", convVtxX, &b_convVtxX);
   fChain->SetBranchAddress("convVtxY", convVtxY, &b_convVtxY);
   fChain->SetBranchAddress("convVtxZ", convVtxZ, &b_convVtxZ);
   //   fChain->SetBranchAddress("convVtxChi2", convVtxChi2, &b_convVtxChi2);
   //   fChain->SetBranchAddress("convVtxNdof", convVtxNdof, &b_convVtxNdof);
   fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convVtxChi2Prob", convVtxChi2Prob, &b_convVtxChi2Prob);
   fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convzOfPrimaryVertexFromTracks", convzOfPrimaryVertexFromTracks, &b_convzOfPrimaryVertexFromTracks);
   fChain->SetBranchAddress("isGenMatched", isGenMatched, &b_isGenMatched);
   fChain->SetBranchAddress("genMatchedPt", genMatchedPt, &b_genMatchedPt);
   fChain->SetBranchAddress("genMatchedEta", genMatchedEta, &b_genMatchedEta);
   fChain->SetBranchAddress("genMatchedPhi", genMatchedPhi, &b_genMatchedPhi);
   fChain->SetBranchAddress("genMomId", genMomId, &b_genMomId);
   fChain->SetBranchAddress("genGrandMomId", genGrandMomId, &b_genGrandMomId);
   fChain->SetBranchAddress("genNSiblings", genNSiblings, &b_genNSiblings);
   fChain->SetBranchAddress("genCalIsoDR03", genCalIsoDR03, &b_genCalIsoDR03);
   fChain->SetBranchAddress("genTrkIsoDR03", genTrkIsoDR03, &b_genTrkIsoDR03);
   fChain->SetBranchAddress("genCalIsoDR04", genCalIsoDR04, &b_genCalIsoDR04);
   fChain->SetBranchAddress("genTrkIsoDR04", genTrkIsoDR04, &b_genTrkIsoDR04);
   Notify();
}

Bool_t mpaReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mpaReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mpaReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void mpaReader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L mpaReader.C
//      Root > mpaReader t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
