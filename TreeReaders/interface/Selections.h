
bool spikeSel(
              bool isEB, bool isEE,
              int seedSeverity,
              int seedRecoFlag
              ) {
  return
    ( (isEB && (seedSeverity !=3 && seedSeverity !=4 ) && (seedRecoFlag != 2) ) || isEE );
}

bool preselection(
                  float leadpt,
                  float subleadpt,
                  float leadeta,
                  float subleadeta,
                  bool leadEBEEGap,
                  bool subleadEBEEGap
                  ) {
    
  bool FilterResult = true;

  if (leadpt<20 || subleadpt<20) FilterResult=false;
  if (abs(leadeta)>2.5 || abs(subleadeta)>2.5) FilterResult=false;
  if (leadEBEEGap || subleadEBEEGap) FilterResult=false;
  
  return FilterResult;
}

bool looseId(
             float pt,
             float ecalRecHitSumEtConeDR04,
             float hcalTowerSumEtConeDR04,
             float trkSumPtHollowConeDR04,
             bool isEB,
             bool isEE,
             float sigmaIetaIeta,
             float hadronicOverEm
             ) {

  if (ecalRecHitSumEtConeDR04 > 4.2 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.2 + 0.0025 * pt) return false;
  if (trkSumPtHollowConeDR04 > 2.0 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.0105) return false;
  if (isEE && sigmaIetaIeta > 0.030) return false;
  if (hadronicOverEm > 0.05) return false;



  return true;

}

bool looseControlId(
                    float pt,
                    float ecalRecHitSumEtConeDR04,
                    float hcalTowerSumEtConeDR04,
                    float trkSumPtHollowConeDR04,
                    bool isEB,
                    bool isEE,
                    float sigmaIetaIeta,
                    float hadronicOverEm
                    ) {

  if (ecalRecHitSumEtConeDR04 > 4.2 + 0.003 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.2 + 0.001 * pt) return false;
  if (trkSumPtHollowConeDR04 > 5.0 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.015) return false;
  if (isEE && sigmaIetaIeta > 0.045) return false;
  if (hadronicOverEm > 0.05) return false;



  return true;

}

bool tightId(
             float pt,
             float ecalRecHitSumEtConeDR04,
             float hcalTowerSumEtConeDR04,
             float trkSumPtHollowConeDR04,
             bool isEB,
             bool isEE,
             float sigmaIetaIeta,
             float hadronicOverEm
             ) {

  /*
  if (hadronicOverEm > 0.03) return false;
  if (ecalRecHitSumEtConeDR04 > 2.4 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 1.0 + 0.001 * pt) return false;
  if (trkSumPtHollowConeDR04 > 0.9 + 0.0025 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.010) return false;
  if (isEE && sigmaIetaIeta > 0.028) return false;
  */

  if (hadronicOverEm > 0.02) return false;
  if (ecalRecHitSumEtConeDR04 > 2  + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.0 + 0.0025 * pt) return false;
  if (trkSumPtHollowConeDR04 > 1.5 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.010) return false;
  if (isEE && sigmaIetaIeta > 0.028) return false;



  return true;

}

bool tightControlId(
                    float pt,
                    float ecalRecHitSumEtConeDR04,
                    float hcalTowerSumEtConeDR04,
                    float trkSumPtHollowConeDR04,
                    bool isEB,
                    bool isEE,
                    float sigmaIetaIeta,
                    float hadronicOverEm
                    ) {

  if (hadronicOverEm > 0.03) return false;
  if (ecalRecHitSumEtConeDR04 > 2.4 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 1.0 + 0.001 * pt) return false;
  if (trkSumPtHollowConeDR04 > 3.9 + 0.0025 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.015) return false;
  if (isEE && sigmaIetaIeta > 0.042) return false;



  return true;

}

bool convSel(
             int nTracks,
             bool convVtxValid,
             float convVtxChi2Prob,
             float convDPhiTracksAtVtx,
             float convpairCotThetaSeparation,
             float EoP,
             float R
              ) {

  if (nTracks != 2) return false;
  if (!convVtxValid) return false;
  if (convVtxChi2Prob < 0.0005) return false;

  //if (fabs(convDPhiTracksAtVtx) > 0.2) return false;
  //if (fabs(convpairCotThetaSeparation) > 0.3) return false;
  //if (EoP > 3.) return false;
  //    if (PoE > 3.) return false;

  return true;

}


int photonCategory (bool pixMatch, float r9, int nTracks,  float convVtxChi2Prob, float etOverPt, float R) {
  
  int cate=0;
  //  if (  !pixMatch ) { 
    if ( r9>0.93 && !pixMatch ) {
      cate=1; // golden
    } else if ( r9<=0.93 && nTracks==2 && convVtxChi2Prob >0.0005 && etOverPt< 3) {
      cate=2; // good reconstructed conversion
    } else if ( r9<=0.93 && nTracks==2 && convVtxChi2Prob <=0.0005 ) {
      cate=3; // poor reconstructed conversions
    } else if ( r9<=0.93 && nTracks<2 ) {
      cate=4;  // no tracks are reconstructed
    }
    //  }

  return cate;
  
}

int diPhotonCategory ( int c1, int c2) {
  
  int cate=0;
  if ( c1==1 && c2==1 ) {
    cate=1; // two golden photons

  } else if (  (c1 ==1 && c2 ==2) ||  (c1 ==2 && c2 ==1)  ) {
    cate=2; // 1 golden, 1 good conversion

  } else if (  (c1 ==1 && c2 ==3) ||  (c1 ==3 && c2 ==1)  ) {
    cate=3; // 1 golden, 1 poor conversion

  } else if (   c1 == 2 &&   c2 ==2 ) {
    cate=4; // 2 good conversions

  }  else if (  (c1 ==1 && c2==4) || (c1==4 && c2==1) ) {
    cate=5; // 1 golden, 1 no tracks reconstructed
    
  }  else if (  (c1==1 && c2==5) ||  (c1 ==5 && c2==1)  ) {
    cate=6; // 1 golden, rest
  }

  return cate;
  
}

bool MarcosCut(
               float pt,
             float ecalRecHitSumEtConeDR04,
             float hcalTowerSumEtConeDR04,
             float trkSumPtHollowConeDR04,
             bool hasPixelSeed,
             bool isEB,
             bool isEE,
             float sigmaIetaIeta,
             float hadronicOverEm
             ) {

  if (!isEB && !isEE) return false;
  if (hadronicOverEm > 0.02) return false;
  if (hasPixelSeed) return false;
  if (ecalRecHitSumEtConeDR04 > 2.0 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.0 + 0.0025 * pt) return false;
  if (trkSumPtHollowConeDR04 > 1.5 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.010) return false;
  if (isEE && sigmaIetaIeta > 0.028) return false;


  return true;

}

int MarcosCutCategory(float LeadR9, bool LeadisEB, bool LeadisEE, float SubLeadR9, bool SubLeadisEB, bool SubLeadisEE) {

  int ReturnValue = 5;

  if (LeadR9 > 0.93 && SubLeadR9 > 0.93 && LeadisEB && SubLeadisEB) {
    ReturnValue = 0;
  } else if ((LeadR9 < 0.93 || SubLeadR9 < 0.93) && LeadisEB && SubLeadisEB) {
    ReturnValue = 1;
  } else if (LeadR9 > 0.93 && SubLeadR9 > 0.93 && (LeadisEE || SubLeadisEE)) {
    ReturnValue = 2;
  } else if ((LeadR9 < 0.93 || SubLeadR9 < 0.93) && (LeadisEE || SubLeadisEE)) {
    ReturnValue = 3;
  }
  
  return ReturnValue;

}
