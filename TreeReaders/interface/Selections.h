bool spikeSel(
              bool isEB, bool isEE,
              int seedSeverity,
              int seedRecoFlag
              ) {
  return
    ( (isEB && (seedSeverity !=3 && seedSeverity !=4 ) && (seedRecoFlag != 2) ) || isEE );
}

bool preselection(
                  float pt,
                  float eta,
                  float scEta,
                  float hadronicOverEm
                  ) {
    
  if (pt < 21.) return false;
  if (fabs(eta) > 2.5) return false;
  if (fabs(scEta) > 1.4442 && fabs(scEta) < 1.566) return false;
  if (hadronicOverEm > 0.05) return false;

  return true;
}

bool looseId(
             float pt,
             float ecalRecHitSumEtConeDR04,
             float hcalTowerSumEtConeDR04,
             float trkSumPtHollowConeDR04,
             bool isEB,
             bool isEE,
             float sigmaIetaIeta
             ) {

  if (ecalRecHitSumEtConeDR04 > 4.2 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.2 + 0.0025 * pt) return false;
  if (trkSumPtHollowConeDR04 > 2.0 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.0105) return false;
  if (isEE && sigmaIetaIeta > 0.030) return false;



  return true;

}

bool looseControlId(
                    float pt,
                    float ecalRecHitSumEtConeDR04,
                    float hcalTowerSumEtConeDR04,
                    float trkSumPtHollowConeDR04,
                    bool isEB,
                    bool isEE,
                    float sigmaIetaIeta
                    ) {

  if (ecalRecHitSumEtConeDR04 > 4.2 + 0.003 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 2.2 + 0.001 * pt) return false;
  if (trkSumPtHollowConeDR04 > 5.0 + 0.001 * pt) return false;
  if (isEB && sigmaIetaIeta > 0.015) return false;
  if (isEE && sigmaIetaIeta > 0.045) return false;



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

  if (hadronicOverEm > 0.03) return false;
  if (ecalRecHitSumEtConeDR04 > 2.4 + 0.006 * pt) return false;
  if (hcalTowerSumEtConeDR04 > 1.0 + 0.001 * pt) return false;
  if (trkSumPtHollowConeDR04 > 0.9 + 0.0025 * pt) return false;
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
             float PoE
             ) {

  if (nTracks != 2) return false;
  if (!convVtxValid) return false;
  if (convVtxChi2Prob < 0.0005) return false;
  if (fabs(convDPhiTracksAtVtx) > 0.2) return false;
  if (fabs(convpairCotThetaSeparation) > 0.3) return false;
  if (EoP > 3.) return false;
  //    if (PoE > 3.) return false;

  return true;

}
