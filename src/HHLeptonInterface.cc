#include "Tools/Tools/interface/HHLeptonInterface.h"




// Constructors

HHLeptonInterface::HHLeptonInterface () {};

// Destructor
HHLeptonInterface::~HHLeptonInterface() {}

std::vector<int> HHLeptonInterface::get_dau_indexes(
    fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
    fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
    bRVec Muon_mediumId, bRVec Muon_tightId,
    fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
    bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
    bRVec Electron_mvaFall17V2Iso_WP90, fRVec Electron_pfRelIso03_all,
    fRVec Electron_dxy, fRVec Electron_dz,
    fRVec Tau_pt, fRVec Tau_eta, fRVec Tau_phi, fRVec Tau_mass,
    iRVec Tau_idDeepTau2017v2p1VSmu, iRVec Tau_idDeepTau2017v2p1VSe,
    iRVec Tau_idDeepTau2017v2p1VSjet, fRVec Tau_rawDeepTau2017v2p1VSjet,
    iRVec Tau_dz, iRVec Tau_decayMode    
  )
{
  // mutau channel
  std::vector<int> goodmuons;
  for (size_t imuon = 0; imuon < Muon_pt.size(); imuon ++) {
    if (fabs(Muon_eta[imuon]) > 2.1 || Muon_pfRelIso04_all[imuon] > 0.15
        || fabs(Muon_dxy[imuon]) > 0.045 || fabs(Muon_dz[imuon]) > 0.2
        || !Muon_tightId[imuon])
      continue;
    goodmuons.push_back(imuon);
  }  // loop over muons
  if (goodmuons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      if (Tau_idDeepTau2017v2p1VSmu[itau] < 15 || Tau_idDeepTau2017v2p1VSe[itau] < 7
          || Tau_idDeepTau2017v2p1VSjet[itau] < 1)
        continue;
      if (fabs(Tau_dz[itau]) > 0.2)
        continue;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        continue;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.)
        continue;
      goodtaus.push_back(itau);
    } // loop over taus

    std::vector<tau_pair> tau_pairs;
    
    for (auto & imuon: goodmuons) {
      auto muon_tlv = TLorentzVector();
      muon_tlv.SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon], Muon_phi[imuon], Muon_mass[imuon]);
      for (auto & itau: goodtaus) {
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(muon_tlv) < 0.5)
          continue;

        // ******************** //
        // MISSING: TRIGGER REQ //
        // ******************** //

        tau_pairs.push_back(tau_pair({imuon, (int) Muon_pfRelIso04_all[imuon], Muon_pt[imuon],
          itau, (int) Tau_idDeepTau2017v2p1VSjet[itau], Tau_pt[itau], 0}));

      } // loop over goodtaus
    } // loop over goodmuons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(tau_pairs[0].index1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {-1, -1, -1, -1};

      return {0, tau_pairs[0].index1, tau_pairs[0].index2, 0};
    }
  }  // goodmuons stuff


  // etau channel
  std::vector<int> goodelectrons;
  for (size_t iele = 0; iele < Electron_pt.size(); iele ++) {
    if (!Electron_mvaFall17V2Iso_WP80[iele]
        || fabs(Electron_dxy[iele]) > 0.045 || fabs(Electron_dz[iele]) > 0.2)
      continue;
    goodelectrons.push_back(iele);
  }  // loop over electrons
  if (goodelectrons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      if (Tau_idDeepTau2017v2p1VSmu[itau] < 15 || Tau_idDeepTau2017v2p1VSe[itau] < 7
          || Tau_idDeepTau2017v2p1VSjet[itau] < 1)
        continue;
      if (fabs(Tau_dz[itau]) > 0.2)
        continue;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        continue;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.)
        continue;
      goodtaus.push_back(itau);
    } // loop over taus

    std::vector<tau_pair> tau_pairs;
    
    for (auto & iele: goodelectrons) {
      auto electron_tlv = TLorentzVector();
      electron_tlv.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele],
        Electron_phi[iele], Electron_mass[iele]);
      for (auto & itau: goodtaus) {
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(electron_tlv) < 0.5)
          continue;

        // ******************** //
        // MISSING: TRIGGER REQ //
        // ******************** //

        tau_pairs.push_back(tau_pair({iele, (int) Electron_pfRelIso03_all[iele], Electron_pt[iele],
          itau, (int) Tau_idDeepTau2017v2p1VSjet[itau], Tau_pt[itau], 0}));

      } // loop over goodtaus
    } // loop over goodelectrons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(-1, tau_pairs[0].index1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {-1, -1, -1, -1};

      return {1, tau_pairs[0].index1, tau_pairs[0].index2, 0};
    }
  }  // goodmuons stuff

  std::vector<int> goodtaus;
  for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
    if (Tau_idDeepTau2017v2p1VSmu[itau] < 1 || Tau_idDeepTau2017v2p1VSe[itau] < 3
        || Tau_idDeepTau2017v2p1VSjet[itau] < 1)
      continue;
    if (fabs(Tau_dz[itau]) > 0.2)
      continue;
    if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
        && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
      continue;
    // common tau pt req for both single and cross triggers
    if (Tau_pt[itau] <= 20.)
      continue;
    goodtaus.push_back(itau);
  } // loop over taus
  
  if (goodtaus.size() >= 2) {
    std::vector<tau_pair> tau_pairs;
    for (auto & itau1 : goodtaus) {
      auto tau1_tlv = TLorentzVector();
      tau1_tlv.SetPtEtaPhiM(Tau_pt[itau1], Tau_eta[itau1], Tau_phi[itau1], Tau_mass[itau1]);
      for (auto & itau2 : goodtaus) {
        if (itau1 == itau2)
          continue;
        auto tau2_tlv = TLorentzVector();
        tau2_tlv.SetPtEtaPhiM(Tau_pt[itau2], Tau_eta[itau2], Tau_phi[itau2], Tau_mass[itau2]);
        if (tau1_tlv.DeltaR(tau2_tlv) < 0.5)
          continue;

        int pass_vbf = 0;
        // ********************************** //
        // MISSING: DITAU AND VBF TRIGGER REQ //
        // ********************************** //
        tau_pairs.push_back(tau_pair({itau1, (int) Tau_rawDeepTau2017v2p1VSjet[itau1], Tau_pt[itau1],
          itau2, (int) Tau_rawDeepTau2017v2p1VSjet[itau2], Tau_pt[itau2], pass_vbf}));
      }
    }
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      if (lepton_veto(-1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {-1, -1, -1, -1};
      return {2, tau_pairs[0].index1, tau_pairs[0].index2, tau_pairs[0].isVBFtrigger};
    }
  }

  return {-1, -1, -1, -1};
}


bool HHLeptonInterface::lepton_veto(int muon_index, int electron_index,
    fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_dz, fRVec Muon_dxy,
    fRVec Muon_pfRelIso04_all, bRVec Muon_mediumId, bRVec Muon_tightId,
    fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_dz, fRVec Electron_dxy,
    bRVec Electron_mvaFall17V2noIso_WP90, bRVec Electron_mvaFall17V2Iso_WP90,
    fRVec Electron_pfRelIso03_all)
{
  size_t nleps = 0;
  for (size_t imuon = 0; imuon < Muon_pt.size(); imuon++) {
    if ((int) imuon == muon_index)  // same muon than the chosen in the pair
      continue;
    if (fabs(Muon_eta[imuon]) > 2.4 || Muon_pt[imuon] < 10 || fabs(Muon_dz[imuon]) > 0.2
        || fabs(Muon_dxy[imuon]) > 0.045 || Muon_pfRelIso04_all[imuon] > 0.3)
      continue;
    if (!Muon_mediumId[imuon] and !Muon_tightId[imuon])
      continue;
    nleps++;
  }
  for (size_t iele = 0; iele < Electron_pt.size(); iele++) {
    if ((int) iele == electron_index)  // same electron than the chosen in the pair
      continue;
    if (fabs(Electron_eta[iele]) > 2.5 || Electron_pt[iele] < 10 || fabs(Electron_dz[iele]) > 0.2
        || fabs(Electron_dxy[iele]) > 0.045)
      continue;
    if (!((Electron_pfRelIso03_all[iele] < 0.3 && Electron_mvaFall17V2noIso_WP90[iele])
        || Electron_mvaFall17V2Iso_WP90[iele]))
      continue;
    nleps++;
  }
  return (nleps > 0);
}