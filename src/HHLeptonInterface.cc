#include "Tools/Tools/interface/HHLeptonInterface.h"

#include <algorithm>
#include <cassert>

// Constructors

HHLeptonInterface::HHLeptonInterface (
    int vvvl_vsjet, int vl_vse, int vvl_vse, int t_vsmu, int vl_vsmu,
    double BT_VsMu_threshold, double BT_VsE_threshold, double BT_VsJet_threshold) {
  vvvl_vsjet_ = vvvl_vsjet;
  vl_vse_ = vl_vse;
  vvl_vse_ = vvl_vse;
  t_vsmu_ = t_vsmu;
  vl_vsmu_ = vl_vsmu;
  BT_VsMu_threshold_ = BT_VsMu_threshold;
  BT_VsE_threshold_ = BT_VsE_threshold;
  BT_VsJet_threshold_ = BT_VsJet_threshold;
};

// Destructor
HHLeptonInterface::~HHLeptonInterface() {}

std::pair<lepton_output, cutflow_output> HHLeptonInterface::get_boosted_dau_indexes(
      bool doGenCutFlow,
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
      fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
      bRVec Muon_looseId, bRVec Muon_mediumId, bRVec Muon_tightId, iRVec Muon_charge,
      sRVec Muon_genPartIdx,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
      bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
      bRVec Electron_mvaFall17V2Iso_WP90, iRVec Electron_vidNestedWPBitmap, fRVec Electron_pfRelIso03_all,
      fRVec Electron_dxy, fRVec Electron_dz, iRVec Electron_charge,
      sRVec Electron_genPartIdx,
      fRVec boostedTau_pt, fRVec boostedTau_eta, fRVec boostedTau_phi, fRVec boostedTau_mass,
      iRVec boostedTau_idDeepTauVSmu, iRVec boostedTau_idDeepTauVSe,
      iRVec boostedTau_idDeepTauVSjet, fRVec boostedTau_rawDeepTauVSjet,
      iRVec boostedTau_decayMode, iRVec boostedTau_charge,
      sRVec boostedTau_genPartIdx, ROOT::VecOps::RVec<UChar_t> boostedTau_genPartFlav,
      iRVec boostedTau_muonCount, std::array<sRVec, 3> BT_muon_idx, 
      std::array<fRVec, 3> BT_muon_pt, std::array<fRVec, 3> BT_muon_correctedIso,
      iRVec boostedTau_electronCount, std::array<sRVec, 3> BT_electron_idx, 
      std::array<fRVec, 3> BT_electron_pt, std::array<fRVec, 3> BT_electron_correctedIso,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<trig_req> mutau_triggers, std::vector<trig_req> etau_triggers,
      std::vector<trig_req> tautau_triggers,
      int GenPairType, int genDau1_genPart_idx, int genDau2_genPart_idx
    )
{
  auto boostedTauGenMatch = [&] (int iBoostedTau_reco, int genVisTauToMatchTo) -> bool { 
    // Gen-matching boostedTau
    // genPartFlav==5 means true hadronic tau, in which case genPartIdx points to GenVisTau collection (otherwise points to GenParticle collection)
    // This needs to be used in case genVisTauToMatchTo is from the GenVisTau collection
    //std::cout << "GenMatch genPartFlav=" << (int)boostedTau_genPartFlav[iBoostedTau_reco] << " genPartIdx=" << (int)boostedTau_genPartIdx[iBoostedTau_reco] << std::endl;
    return boostedTau_genPartFlav[iBoostedTau_reco] == 5 && boostedTau_genPartIdx[iBoostedTau_reco] == genVisTauToMatchTo;
  };
  // Fail reasons of the gen daus, recording why the dau did not pass analysis selections
  // Reco is set by default, we remove the flag in case we find a gen-matched reco candidate of the right kind (depending on pairType)
  cutflow_output cutflow{};
  cutflow.dau1_fail.Reco = true;
  cutflow.dau2_fail.Reco = true;
  std::vector<int> goodBoostedTaus;
  int genDau2_boostedTau_idx = -1; // Index in boostedTau collection of genmatched boostedTau
  for (size_t itau = 0; itau < boostedTau_pt.size(); itau ++) {
    FailReason failReason;
    if (boostedTau_pt[itau] < 40) failReason.Pt = true; // Pt threshold arbitrary (Wisconsin uses >20, central Nano uses >40)
    //if (boostedTau_idDeepTauVSmu[itau] < BT_VsMu_threshold_) failReason.TauIdVsMu = true;
    //if (boostedTau_idDeepTauVSe[itau] < BT_VsE_threshold_) failReason.TauIdVsE = true;
    if (boostedTau_rawDeepTauVSjet[itau] < BT_VsJet_threshold_) failReason.TauIdVsJet = true;
    
    if (boostedTau_decayMode[itau] != 0 && boostedTau_decayMode[itau] != 1
        && boostedTau_decayMode[itau] != 10 && boostedTau_decayMode[itau] != 11)
      failReason.TauDM = true;
    if (failReason.pass())
      goodBoostedTaus.push_back(itau);

    if (doGenCutFlow) {
      if (boostedTauGenMatch(itau, genDau2_genPart_idx)) {
        cutflow.dau2_fail = failReason; // this will unset the Reco flag
        genDau2_boostedTau_idx = itau;
      } else if (GenPairType == 2 && boostedTauGenMatch(itau, genDau1_genPart_idx))
        cutflow.dau1_fail = failReason; // same
    }
  } // loop over boostedTaus

  // ---------------- MU-TAU channel  ----------------------
  std::vector<tau_pair> tau_pairs;
  bool foundMuonInBoostedTauCollection = false; // true if a genmatched muon has been found in the boostedTau muon collection
  for (auto & itau: goodBoostedTaus) {
    for (int imuonFromBT = 0; imuonFromBT < std::min(boostedTau_muonCount[itau], 3); imuonFromBT++) {
      FailReason failReason;
      int imuon = BT_muon_idx[imuonFromBT][itau];
      if (imuon < 0) {
        std::cerr << "Muon that was not selected in NanoAOD electron" << std::endl;// Should be very rare
        continue;
      }
      assert(imuon < (int)Muon_eta.size());

      if (fabs(Muon_eta[imuon]) >= 2.4) failReason.Eta = true;
      if (fabs(Muon_dxy[imuon]) > 0.045 || fabs(Muon_dz[imuon]) > 0.2) failReason.Vertex = true;
      if (!Muon_looseId[imuon]) failReason.LeptonID = true; // this cannot happen in theory, cut already at Nano level
      // No ISO requirement yet
      
      if (Muon_pt[imuon] < 20) failReason.Pt = true; // arbitrary threshold

      if (BT_muon_correctedIso[imuonFromBT][itau]/BT_muon_pt[imuonFromBT][itau] > 0.25) failReason.LeptonIso = true;

      auto muon_tlv = TLorentzVector();
      muon_tlv.SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon], Muon_phi[imuon], Muon_mass[imuon]);
      auto tau_tlv = TLorentzVector();
      tau_tlv.SetPtEtaPhiM(boostedTau_pt[itau], boostedTau_eta[itau], boostedTau_phi[itau], boostedTau_mass[itau]);
      // No trigger matching for MET trigger
      // if (!pass_trigger(
      //       muon_tlv.Pt(), muon_tlv.Eta(), muon_tlv.Phi(), 13,
      //       tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
      //       mutau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi))
      //     continue;

      if (failReason.pass())
        tau_pairs.push_back(tau_pair({imuon, BT_muon_correctedIso[imuonFromBT][itau]/BT_muon_pt[imuonFromBT][itau], Muon_pt[imuon],
            itau, boostedTau_rawDeepTauVSjet[itau], boostedTau_pt[itau], 0, 0}));
        
      if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[imuon] == genDau1_genPart_idx) {
        cutflow.dau1_fail = failReason;
        foundMuonInBoostedTauCollection = true;
      }
    }
  }

  if (tau_pairs.size() > 0) {
    std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSortHybrid);
    int ind1 = tau_pairs[0].index1;
    int ind2 = tau_pairs[0].index2;
    int isOS = (int) (Muon_charge[ind1] != boostedTau_charge[ind2]);
    
    if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[ind1] != genDau1_genPart_idx)
      cutflow.dau1_fail.WrongPair = true;
    if (doGenCutFlow && !boostedTauGenMatch(ind2, genDau2_genPart_idx))
      cutflow.dau2_fail.WrongPair = true;

    if (lepton_veto(tau_pairs[0].index1, -1,
        Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
        Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
        Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
        Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
        Electron_pfRelIso03_all))
      return {lepton_output({
        -1, -1, -1, -1, -1, -1,
        -1., -1., -1., -1, -1, -1, -1,
        -1., -1., -1, -1, -1, -1}), 
      cutflow.setLeptonVeto().setWrongChannel(GenPairType != 0)};


    return {lepton_output({0, ind1, ind2, 0, 0, isOS,
      Muon_eta[ind1], Muon_phi[ind1], Muon_pfRelIso04_all[ind1], -1, -1, -1, -1,
      boostedTau_eta[ind2], boostedTau_phi[ind2], boostedTau_decayMode[ind2],
      boostedTau_idDeepTauVSe[ind2], boostedTau_idDeepTauVSmu[ind2],
      boostedTau_idDeepTauVSjet[ind2]}), 
      cutflow.setWrongChannel(GenPairType != 0)};
  }

  // Finding muons that were not in the boostedTau collections
  // gen study on muons : because boostedTau muons branches have selection on ID, and we want to separate reco from ID in cutflow
  if (doGenCutFlow && GenPairType == 0 && !foundMuonInBoostedTauCollection) {
    for (size_t imuon = 0; imuon < Muon_pt.size(); imuon ++) {
      if (Muon_genPartIdx[imuon] == genDau1_genPart_idx) {
        // We found a reco muon matching truth : Could be because of deltaR failure or because of ID fail
        cutflow.dau1_fail.Reco = false;

        // Now we try to reproduce selections done at NanoAOD level

        if (genDau2_boostedTau_idx >= 0) {
          auto muon_tlv = TLorentzVector();
          muon_tlv.SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon],
              Muon_phi[imuon], Muon_mass[imuon]);
          auto tau_tlv = TLorentzVector();
          tau_tlv.SetPtEtaPhiM(boostedTau_pt[genDau2_boostedTau_idx], boostedTau_eta[genDau2_boostedTau_idx], boostedTau_phi[genDau2_boostedTau_idx], boostedTau_mass[genDau2_boostedTau_idx]);

          if (muon_tlv.DeltaR(tau_tlv) > 0.7 || muon_tlv.DeltaR(tau_tlv) <= 0.05)
            cutflow.deltaR = true;
        }
        cutflow.dau1_fail.LeptonID = !Muon_looseId[imuon];
      }
    }  // loop over electrons
  }

  // ---------------- E-TAU channel  ----------------------
  tau_pairs.clear();

  bool foundElectronInBoostedTauCollection = false;
  for (auto & itau: goodBoostedTaus) {
    // boostedTau_electronCount can be up to 4 which is a bug
    for (int ielectronFromBT = 0; ielectronFromBT < std::min(boostedTau_electronCount[itau], 3); ielectronFromBT++) {
      FailReason failReason;
      int iele = BT_electron_idx[ielectronFromBT][itau];
      if (iele < 0) {
        std::cerr << "Electron that was not selected in NanoAOD electron" << std::endl;
        continue;
      }
      assert(iele < (int)Electron_eta.size());

      // Isolation is done already when selecting electrons in boostedTau cone
      //if (!Electron_mvaFall17V2noIso_WP90[iele]) failReason.LeptonID = true;  // NoISO MVA 
      if (fabs(Electron_dxy[iele]) > 0.045 || fabs(Electron_dz[iele]) > 0.2) failReason.Vertex = true;
      if (fabs(Electron_eta[iele]) >= 2.5
        || (fabs(Electron_eta[iele]) > 1.44 && fabs(Electron_eta[iele]) < 1.57)) // exclude barrel/endcap transition region (no SFs available for EGamma)
          failReason.Eta = true;
      if (Electron_pt[iele] < 30) failReason.Pt = true; // arbitrary threshold

      // Isolation cut from https://indico.cern.ch/event/1420456/contributions/5991451/attachments/2878187/5041085/HeavyMassResonance_B2G_DIB_June_14_2024.pdf
      float isoOverPt = BT_electron_correctedIso[ielectronFromBT][itau]/BT_electron_pt[ielectronFromBT][itau];
      if ((
          fabs(Electron_eta[iele]) < 1.479 &&  isoOverPt >=  0.112 + 0.506 / BT_electron_pt[ielectronFromBT][itau]
        ) || (
          fabs(Electron_eta[iele]) >= 1.479 &&  isoOverPt >=  0.108 + 0.963 / BT_electron_pt[ielectronFromBT][itau]
      ))
        failReason.LeptonIso = true;
      
      auto electron_tlv = TLorentzVector();
      electron_tlv.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele],
          Electron_phi[iele], Electron_mass[iele]);
      auto tau_tlv = TLorentzVector();
      tau_tlv.SetPtEtaPhiM(boostedTau_pt[itau], boostedTau_eta[itau], boostedTau_phi[itau], boostedTau_mass[itau]);
      // no trigger matching for MET trigger
      // if (!pass_trigger(
      //       electron_tlv.Pt(), electron_tlv.Eta(), electron_tlv.Phi(), 11,
      //       tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
      //       mutau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi))
      //     continue;
      if (failReason.pass())
        tau_pairs.push_back(tau_pair({iele, BT_electron_correctedIso[ielectronFromBT][itau]/BT_electron_pt[ielectronFromBT][itau], Electron_pt[iele],
            itau, boostedTau_rawDeepTauVSjet[itau], boostedTau_pt[itau], 0, 0}));
      
      if (doGenCutFlow && GenPairType == 1 && Electron_genPartIdx[iele] == genDau1_genPart_idx) {
        cutflow.dau1_fail = failReason;
        foundElectronInBoostedTauCollection = true;
      }
        
    }
  }

  if (tau_pairs.size() > 0) {
    std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSortHybrid);

    int ind1 = tau_pairs[0].index1;
    int ind2 = tau_pairs[0].index2;
    int isOS = (int) (Electron_charge[ind1] != boostedTau_charge[ind2]);

    if (doGenCutFlow && Electron_genPartIdx[ind1] != genDau1_genPart_idx)
      cutflow.dau1_fail.WrongPair = true;
    if (doGenCutFlow && !boostedTauGenMatch(ind2, genDau2_genPart_idx))
      cutflow.dau2_fail.WrongPair = true;

    if (lepton_veto(-1, tau_pairs[0].index1,
        Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
        Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
        Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
        Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
        Electron_pfRelIso03_all))
      return {lepton_output({
        -1, -1, -1, -1, -1, -1,
        -1., -1., -1., -1, -1, -1, -1,
        -1., -1., -1, -1, -1, -1}), 
      cutflow.setLeptonVeto().setWrongChannel(GenPairType != 1)
      };

    return {lepton_output({1, ind1, ind2, 0, 0, isOS,
      Electron_eta[ind1], Electron_phi[ind1], Electron_pfRelIso03_all[ind1], -1, -1, -1, -1,
      boostedTau_eta[ind2], boostedTau_phi[ind2], boostedTau_decayMode[ind2],
      boostedTau_idDeepTauVSe[ind2], boostedTau_idDeepTauVSmu[ind2],
      boostedTau_idDeepTauVSjet[ind2]}), 
    cutflow.setWrongChannel(GenPairType != 1)
    };
  }

  // Finding electrons that were not in the boostedTau collections
  // gen study on electrons : because boostedTau electron branches have selection on ID, and we want to separate reco from ID in cutflow
  if (doGenCutFlow && GenPairType == 1 && !foundElectronInBoostedTauCollection) {
    for (size_t iele = 0; iele < Electron_pt.size(); iele ++) {
      if (Electron_genPartIdx[iele] == genDau1_genPart_idx) {
        // We found a reco electron matching truth : Could be because of deltaR failure or because of ID fail
        cutflow.dau1_fail.Reco = false;

        // Now we try to reproduce selections done at NanoAOD level

        if (genDau2_boostedTau_idx >= 0) {
          auto electron_tlv = TLorentzVector();
          electron_tlv.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele],
              Electron_phi[iele], Electron_mass[iele]);
          auto tau_tlv = TLorentzVector();
          tau_tlv.SetPtEtaPhiM(boostedTau_pt[genDau2_boostedTau_idx], boostedTau_eta[genDau2_boostedTau_idx], boostedTau_phi[genDau2_boostedTau_idx], boostedTau_mass[genDau2_boostedTau_idx]);

          if (electron_tlv.DeltaR(tau_tlv) > 0.6 || electron_tlv.DeltaR(tau_tlv) <= 0.05)
            cutflow.deltaR = true;
        }

        // Checking electron ID using Electron_vidNestedWPBitmap
        // Electron_vidNestedWPBitmap is split in 10 groups of 3bits (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2)
        // Each group consists of 3bits that has to be interpreted as a number (not a bitmask!), where 0=fail,1=veto,2=loose,3=medium,4=tight (NB 4->100 in binary) 
        // Other documenation on this : https://cms-talk.web.cern.ch/t/nanoaod-v9-bitmap/24831/3 https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolCERN2020EgammaShortExercise 
        bool failLooseIdWithoutIso = false;
        for (unsigned bit_nb = 0; bit_nb <= 9; bit_nb++) {
          if (bit_nb == 7) continue; // We ignore bit 7 ie GsfEleRelPFIsoScaledCut
          //                                                   start bit        bits per cut   
          int cut_value = (Electron_vidNestedWPBitmap[iele] >> (bit_nb*3)) & ((1 << 3    ) - 1);
          if (cut_value < 2) // 2 is loose
            failLooseIdWithoutIso = true;
        }
        cutflow.dau1_fail.LeptonID = failLooseIdWithoutIso;
      }
    }  // loop over electrons
  }

  // ---------------- TAU-TAU channel  ----------------------
  tau_pairs.clear();

  if (goodBoostedTaus.size() >= 2) {
    for (auto & itau1 : goodBoostedTaus) {
      auto tau1_tlv = TLorentzVector();
      tau1_tlv.SetPtEtaPhiM(boostedTau_pt[itau1], boostedTau_eta[itau1], boostedTau_phi[itau1], boostedTau_mass[itau1]);
      for (auto & itau2 : goodBoostedTaus) {
        if (itau1 == itau2)
          continue;
        FailReason failReason;
        auto tau2_tlv = TLorentzVector();
        tau2_tlv.SetPtEtaPhiM(boostedTau_pt[itau2], boostedTau_eta[itau2], boostedTau_phi[itau2], boostedTau_mass[itau2]);
        if (tau1_tlv.DeltaR(tau2_tlv) <= 0.8) {
          // trigger matching (not used for MET trigger)
          tau_pairs.push_back(tau_pair({itau1, boostedTau_rawDeepTauVSjet[itau1], boostedTau_pt[itau1],
            itau2, boostedTau_rawDeepTauVSjet[itau2], boostedTau_pt[itau2], 0, 0}));
        } else if (doGenCutFlow && GenPairType == 2 && // failed deltaR : log it in case boostedTaus are both genmatched
            ((boostedTauGenMatch(itau1, genDau1_genPart_idx) && boostedTauGenMatch(itau2, genDau2_genPart_idx))
            || (boostedTauGenMatch(itau1, genDau2_genPart_idx) && boostedTauGenMatch(itau2, genDau1_genPart_idx)))) {
              cutflow.deltaR = true;
        }
      }
    }
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);
      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (boostedTau_charge[ind1] != boostedTau_charge[ind2]);

      if (doGenCutFlow) {
        if (!(boostedTauGenMatch(ind1, genDau1_genPart_idx) || boostedTauGenMatch(ind2, genDau1_genPart_idx)))
          cutflow.dau1_fail.WrongPair = true;
        if (!(boostedTauGenMatch(ind1, genDau2_genPart_idx) || boostedTauGenMatch(ind2, genDau2_genPart_idx)))
          cutflow.dau2_fail.WrongPair = true;
      }

      if (lepton_veto(-1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1}), 
          cutflow.setLeptonVeto().setWrongChannel(GenPairType != 2)};

      return {lepton_output({2, ind1, ind2,
        0, 0, isOS,
        boostedTau_eta[ind1], boostedTau_phi[ind1], -1., boostedTau_decayMode[ind1],
        boostedTau_idDeepTauVSe[ind1], boostedTau_idDeepTauVSmu[ind1],
        boostedTau_idDeepTauVSjet[ind1],
        boostedTau_eta[ind2], boostedTau_phi[ind2], boostedTau_decayMode[ind2],
        boostedTau_idDeepTauVSe[ind2], boostedTau_idDeepTauVSmu[ind2],
        boostedTau_idDeepTauVSjet[ind2]}), 
        cutflow.setWrongChannel(GenPairType != 2)
        };
    }
  }

  return {lepton_output({
    -1, -1, -1, -1, -1, -1,
    -1., -1., -1., -1, -1, -1, -1,
    -1., -1., -1, -1, -1, -1}), 
    cutflow.setWrongChannel(!(GenPairType==0 || GenPairType==1 || GenPairType==2))};
}

std::pair<lepton_output, cutflow_output>  HHLeptonInterface::get_dau_indexes(
    bool doGenCutFlow,
    fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
    fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
    bRVec Muon_mediumId, bRVec Muon_tightId, iRVec Muon_charge,
    sRVec Muon_genPartIdx,
    fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
    bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
    bRVec Electron_mvaFall17V2Iso_WP90, fRVec Electron_pfRelIso03_all,
    fRVec Electron_dxy, fRVec Electron_dz, iRVec Electron_charge,
    sRVec Electron_genPartIdx,
    fRVec Tau_pt, fRVec Tau_eta, fRVec Tau_phi, fRVec Tau_mass,
    iRVec Tau_idDeepTauVSmu, iRVec Tau_idDeepTauVSe,
    iRVec Tau_idDeepTauVSjet, fRVec Tau_rawDeepTauVSjet,
    fRVec Tau_dz, iRVec Tau_decayMode, iRVec Tau_charge,
    sRVec Tau_genPartIdx, ROOT::VecOps::RVec<UChar_t> Tau_genPartFlav,
    iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
    std::vector<trig_req> mutau_triggers, std::vector<trig_req> etau_triggers,
    std::vector<trig_req> tautau_triggers, std::vector<trig_req> tautaujet_triggers,
    std::vector<trig_req> vbf_triggers,
    int GenPairType, int genDau1_genPart_idx, int genDau2_genPart_idx
  )
{
  auto HPSTauGenMatch = [&] (int iTau_reco, int genVisTauToMatchTo) -> bool { 
    // Gen-matching HPS Tau
    // genPartFlav==5 means true hadronic tau, in which case genPartIdx points to GenVisTau collection (otherwise points to GenParticle collection)
    // This needs to be used in case genVisTauToMatchTo is from the GenVisTau collection
    //std::cout << "GenMatch genPartFlav=" << (int)boostedTau_genPartFlav[iBoostedTau_reco] << " genPartIdx=" << (int)boostedTau_genPartIdx[iBoostedTau_reco] << std::endl;
    return Tau_genPartFlav[iTau_reco] == 5 && Tau_genPartIdx[iTau_reco] == genVisTauToMatchTo;
  };
  // Fail reasons of the gen daus, recording why the dau did not pass analysis selections
  // Reco is set by default, we remove the flag in case we find a gen-matched reco candidate of the right kind (depending on pairType)
  cutflow_output cutflow{};
  cutflow.dau1_fail.Reco = true;
  cutflow.dau2_fail.Reco = true;

  // mutau channel
  std::vector<int> goodmuons;
  for (size_t imuon = 0; imuon < Muon_pt.size(); imuon ++) {
    FailReason failReason;
    if (fabs(Muon_eta[imuon]) >= 2.4) failReason.Eta = true;
    if (Muon_pfRelIso04_all[imuon] > 0.15) failReason.LeptonIso = true;
    if (fabs(Muon_dxy[imuon]) > 0.045 || fabs(Muon_dz[imuon]) > 0.2) failReason.Vertex = true;
    if (!Muon_tightId[imuon]) failReason.LeptonID = true;
    
    if (failReason.pass())
      goodmuons.push_back(imuon);

    if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[imuon] == genDau1_genPart_idx)
      cutflow.dau1_fail = failReason;
  }  // loop over muons
  if (goodmuons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      FailReason tau_failReason;
      if (Tau_idDeepTauVSmu[itau] < t_vsmu_) tau_failReason.TauIdVsMu = true;
      if (Tau_idDeepTauVSe[itau] < vl_vse_) tau_failReason.TauIdVsE = true;
      if (Tau_idDeepTauVSjet[itau] < vvvl_vsjet_) tau_failReason.TauIdVsJet = true;
      if (fabs(Tau_dz[itau]) > 0.2) tau_failReason.Vertex = true;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        tau_failReason.TauDM = true;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.) tau_failReason.Pt = true;

      if (tau_failReason.pass())
        goodtaus.push_back(itau);

      if (doGenCutFlow && GenPairType == 0 && HPSTauGenMatch(itau, genDau2_genPart_idx)) {
        cutflow.dau2_fail = tau_failReason; // this will unset the Reco flag
      }
    } // loop over taus

    std::vector<tau_pair> tau_pairs;

    for (auto & imuon: goodmuons) {
      auto muon_tlv = TLorentzVector();
      muon_tlv.SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon], Muon_phi[imuon], Muon_mass[imuon]);
      for (auto & itau: goodtaus) {        
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(muon_tlv) < 0.5) {
          if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[imuon] == genDau1_genPart_idx && HPSTauGenMatch(itau, genDau2_genPart_idx)) 
            cutflow.deltaR = true;
          continue;
        }
          
        if (!pass_trigger(
            muon_tlv.Pt(), muon_tlv.Eta(), muon_tlv.Phi(), 13,
            tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
            mutau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi)) {
          if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[imuon] == genDau1_genPart_idx && HPSTauGenMatch(itau, genDau2_genPart_idx)) 
            cutflow.triggerFail = true;
          continue;
        }

        tau_pairs.push_back(tau_pair({imuon, Muon_pfRelIso04_all[imuon], Muon_pt[imuon],
          itau, Tau_rawDeepTauVSjet[itau], Tau_pt[itau], 0, 0}));

      } // loop over goodtaus
    } // loop over goodmuons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSortHybrid);

      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Muon_charge[ind1] != Tau_charge[ind2]);

      if (doGenCutFlow && GenPairType == 0 && Muon_genPartIdx[ind1] != genDau1_genPart_idx)
        cutflow.dau1_fail.WrongPair = true;
      if (doGenCutFlow && !HPSTauGenMatch(ind2, genDau2_genPart_idx))
        cutflow.dau2_fail.WrongPair = true;

      if (lepton_veto(tau_pairs[0].index1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
          return {lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1}), 
            cutflow.setLeptonVeto().setWrongChannel(GenPairType != 0)};

      return {lepton_output({0, ind1, ind2, 0, 0, isOS,
        Muon_eta[ind1], Muon_phi[ind1], Muon_pfRelIso04_all[ind1], -1, -1, -1, -1,
        Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
        Tau_idDeepTauVSe[ind2], Tau_idDeepTauVSmu[ind2],
        Tau_idDeepTauVSjet[ind2]}),
        cutflow.setWrongChannel(GenPairType != 0)};
    }
  }  // goodmuons stuff


  // etau channel
  std::vector<int> goodelectrons;
  for (size_t iele = 0; iele < Electron_pt.size(); iele ++) {
    FailReason failReason;
    if (!Electron_mvaFall17V2Iso_WP80[iele]) failReason.LeptonID = true;
    if (fabs(Electron_dxy[iele]) > 0.045 || fabs(Electron_dz[iele]) > 0.2) failReason.Vertex = true;
    if (fabs(Electron_eta[iele]) >= 2.5
        || (fabs(Electron_eta[iele]) > 1.44 && fabs(Electron_eta[iele]) < 1.57)) // exclude barrel/endcap transition region (no SFs available for EGamma)
      failReason.Eta = true;
    if (failReason.pass())
      goodelectrons.push_back(iele);

    if (doGenCutFlow && GenPairType == 1 && Electron_genPartIdx[iele] == genDau1_genPart_idx) {
      cutflow.dau1_fail = failReason;
    }
  }  // loop over electrons
  if (goodelectrons.size() >= 1) {
    std::vector<int> goodtaus;
    for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
      FailReason tau_failReason;
      if (Tau_idDeepTauVSmu[itau] < t_vsmu_) tau_failReason.TauIdVsMu = true;
      if (Tau_idDeepTauVSe[itau] < vl_vse_) tau_failReason.TauIdVsE = true;
      if (Tau_idDeepTauVSjet[itau] < vvvl_vsjet_) tau_failReason.TauIdVsJet = true;
      if (fabs(Tau_dz[itau]) > 0.2) tau_failReason.Vertex = true;
      if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
          && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
        tau_failReason.TauDM = true;
      // common tau pt req for both single and cross triggers
      if (Tau_pt[itau] <= 20.) tau_failReason.Pt = true;

      if (tau_failReason.pass())
        goodtaus.push_back(itau);
        
      if (doGenCutFlow && GenPairType == 1 && HPSTauGenMatch(itau, genDau2_genPart_idx)) {
        cutflow.dau1_fail = tau_failReason; // this will unset the Reco flag
      }
    } // loop over taus

    std::vector<tau_pair> tau_pairs;

    for (auto & iele: goodelectrons) {
      auto electron_tlv = TLorentzVector();
      electron_tlv.SetPtEtaPhiM(Electron_pt[iele], Electron_eta[iele],
        Electron_phi[iele], Electron_mass[iele]);
      for (auto & itau: goodtaus) {
        auto tau_tlv = TLorentzVector();
        tau_tlv.SetPtEtaPhiM(Tau_pt[itau], Tau_eta[itau], Tau_phi[itau], Tau_mass[itau]);
        if (tau_tlv.DeltaR(electron_tlv) < 0.5) {
          if (doGenCutFlow && GenPairType == 1 && Electron_genPartIdx[iele] == genDau1_genPart_idx && HPSTauGenMatch(itau, genDau2_genPart_idx)) 
            cutflow.deltaR = true;
          continue;
        }
          
        if (!pass_trigger(
            electron_tlv.Pt(), electron_tlv.Eta(), electron_tlv.Phi(), 11,
            tau_tlv.Pt(), tau_tlv.Eta(), tau_tlv.Phi(), 15,
            etau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi)) {
          if (doGenCutFlow && GenPairType == 1 && Electron_genPartIdx[iele] == genDau1_genPart_idx && HPSTauGenMatch(itau, genDau2_genPart_idx)) 
            cutflow.triggerFail = true;
          continue;
        }
        tau_pairs.push_back(tau_pair({iele, Electron_pfRelIso03_all[iele], Electron_pt[iele],
          itau, Tau_rawDeepTauVSjet[itau], Tau_pt[itau], 0, 0}));

      } // loop over goodtaus
    } // loop over goodelectrons
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSortHybrid);

      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Electron_charge[ind1] != Tau_charge[ind2]);

      if (doGenCutFlow && GenPairType == 1 && Electron_genPartIdx[ind1] != genDau1_genPart_idx)
        cutflow.dau1_fail.WrongPair = true;
      if (doGenCutFlow && !HPSTauGenMatch(ind2, genDau2_genPart_idx))
        cutflow.dau1_fail.WrongPair = true;

      if (lepton_veto(-1, tau_pairs[0].index1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1}),
          cutflow.setLeptonVeto().setWrongChannel(GenPairType != 1)};

      return {lepton_output({1, ind1, ind2, 0, 0, isOS,
          Electron_eta[ind1], Electron_phi[ind1], Electron_pfRelIso03_all[ind1], -1,
          -1, -1, -1,
          Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
          Tau_idDeepTauVSe[ind2], Tau_idDeepTauVSmu[ind2],
          Tau_idDeepTauVSjet[ind2]}),
        cutflow.setWrongChannel(GenPairType != 1)};
    }
  }  // goodmuons stuff

  std::vector<int> goodtaus;
  for (size_t itau = 0; itau < Tau_pt.size(); itau ++) {
    FailReason tau_failReason;
    if (Tau_idDeepTauVSmu[itau] < vl_vsmu_) tau_failReason.TauIdVsMu = true;
    if (Tau_idDeepTauVSe[itau] < vvl_vse_) tau_failReason.TauIdVsE = true;
    if (Tau_idDeepTauVSjet[itau] < vvvl_vsjet_) tau_failReason.TauIdVsJet = true;
    if (fabs(Tau_dz[itau]) > 0.2) tau_failReason.Vertex = true;
    if (Tau_decayMode[itau] != 0 && Tau_decayMode[itau] != 1
        && Tau_decayMode[itau] != 10 && Tau_decayMode[itau] != 11)
      tau_failReason.TauDM = true;
    // common tau pt req for both single and cross triggers
    if (Tau_pt[itau] <= 20.) tau_failReason.Pt = true;
    
    if (tau_failReason.pass())
      goodtaus.push_back(itau);

    if (doGenCutFlow) {
      if (HPSTauGenMatch(itau, genDau2_genPart_idx)) {
        cutflow.dau2_fail = tau_failReason; // this will unset the Reco flag
      } else if (GenPairType == 2 && HPSTauGenMatch(itau, genDau1_genPart_idx))
        cutflow.dau1_fail = tau_failReason; // same
    }
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

        bool isPairGenMatched = doGenCutFlow && GenPairType == 2 && 
              ((HPSTauGenMatch(itau1, genDau1_genPart_idx) && HPSTauGenMatch(itau2, genDau2_genPart_idx))
              || (HPSTauGenMatch(itau1, genDau2_genPart_idx) && HPSTauGenMatch(itau2, genDau1_genPart_idx)));
        if (tau1_tlv.DeltaR(tau2_tlv) < 0.5) {
          if (isPairGenMatched) { // failed deltaR : log it in case Taus are both genmatched
              cutflow.deltaR = true;
          }
          continue;
        }

        int pass_tautaujet = 0;
        int pass_vbf = 0;
        if (!pass_trigger(
            tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
            tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
            tautau_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi)) {
          if (pass_trigger(
              tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
              tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
              tautaujet_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi))
            pass_tautaujet = 1;
          else if (pass_trigger(
              tau1_tlv.Pt(), tau1_tlv.Eta(), tau1_tlv.Phi(), 15,
              tau2_tlv.Pt(), tau2_tlv.Eta(), tau2_tlv.Phi(), 15,
              vbf_triggers, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi))
            pass_vbf = 1;
          else {
            if (isPairGenMatched) 
              cutflow.triggerFail = true;
          }
        }
        tau_pairs.push_back(tau_pair({itau1, Tau_rawDeepTauVSjet[itau1], Tau_pt[itau1],
          itau2, Tau_rawDeepTauVSjet[itau2], Tau_pt[itau2], pass_tautaujet, pass_vbf}));
      }
    }
    if (tau_pairs.size() > 0) {
      std::stable_sort(tau_pairs.begin(), tau_pairs.end(), pairSort);

      int ind1 = tau_pairs[0].index1;
      int ind2 = tau_pairs[0].index2;
      int isOS = (int) (Tau_charge[ind1] != Tau_charge[ind2]);

      if (doGenCutFlow) {
        if (!(HPSTauGenMatch(ind1, genDau1_genPart_idx) || HPSTauGenMatch(ind2, genDau1_genPart_idx)))
          cutflow.dau1_fail.WrongPair = true;
        if (!(HPSTauGenMatch(ind1, genDau2_genPart_idx) || HPSTauGenMatch(ind2, genDau2_genPart_idx)))
          cutflow.dau2_fail.WrongPair = true;
      }

      if (lepton_veto(-1, -1,
          Muon_pt, Muon_eta, Muon_dz, Muon_dxy,
          Muon_pfRelIso04_all, Muon_mediumId, Muon_tightId,
          Electron_pt, Electron_eta, Electron_dz, Electron_dxy,
          Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2Iso_WP90,
          Electron_pfRelIso03_all))
        return {lepton_output({
          -1, -1, -1, -1, -1, -1,
          -1., -1., -1., -1, -1, -1, -1,
          -1., -1., -1, -1, -1, -1}),
          cutflow.setLeptonVeto().setWrongChannel(GenPairType != 2)};

      return {lepton_output({2, ind1, ind2,
        tau_pairs[0].isTauTauJetTrigger, tau_pairs[0].isVBFtrigger, isOS,
        Tau_eta[ind1], Tau_phi[ind1], -1., Tau_decayMode[ind1],
        Tau_idDeepTauVSe[ind1], Tau_idDeepTauVSmu[ind1],
        Tau_idDeepTauVSjet[ind1],
        Tau_eta[ind2], Tau_phi[ind2], Tau_decayMode[ind2],
        Tau_idDeepTauVSe[ind2], Tau_idDeepTauVSmu[ind2],
        Tau_idDeepTauVSjet[ind2]}),
        cutflow.setWrongChannel(GenPairType != 2)};
    }
  }

  return {lepton_output({
    -1, -1, -1, -1, -1, -1,
    -1., -1., -1., -1, -1, -1, -1,
    -1., -1., -1, -1, -1, -1}), cutflow.setWrongChannel(!(GenPairType==0 || GenPairType==1 || GenPairType==2))};
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

bool HHLeptonInterface::pass_trigger(
    float off_pt1, float off_eta1, float off_phi1, int obj_id1,
    float off_pt2, float off_eta2, float off_phi2, int obj_id2,
    std::vector<trig_req> triggers, 
    iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi)
{
  for (auto &trigger: triggers) {
    if (!trigger.pass)
      continue;
    if (off_pt1 < trigger.pt1_offline || off_pt2 < trigger.pt2_offline
        || abs(off_eta1) > trigger.eta1_offline || abs(off_eta2) > trigger.eta2_offline)
      continue;
    if (trigger.bits[0].size() > 0) {
      if (!match_trigger_object(off_eta1, off_phi1, obj_id1, trigger.pt1_online,
          TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, trigger.bits[0]))
        continue;
    }
    if (trigger.bits[1].size() > 0) {
      if (!match_trigger_object(off_eta2, off_phi2, obj_id2, trigger.pt2_online,
          TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, trigger.bits[1]))
        continue;
    }
    return true;
  }
  return false;
}

// trig_pt_threshold : threshold for TrigObj_pt (used for cross triggers, not sure if strictly needed)
bool HHLeptonInterface::match_trigger_object(float off_eta, float off_phi, int obj_id,
    float trig_pt_threshold, iRVec TrigObj_id, iRVec TrigObj_filterBits,
    fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
    std::vector<int> bits)
{
  for (size_t iobj = 0; iobj < TrigObj_id.size(); iobj++) {
    if (TrigObj_id[iobj] != obj_id) continue;
    if (TrigObj_pt[iobj] <= trig_pt_threshold) continue;
    auto const dPhi(std::abs(reco::deltaPhi(off_phi, TrigObj_phi[iobj])));
    auto const dEta(std::abs(off_eta - TrigObj_eta[iobj]));
    auto const delR2(dPhi * dPhi + dEta * dEta);
    if (delR2 > 0.5 * 0.5)
      continue;
    bool matched_bits = true;
    for (auto & bit : bits) {
      if ((TrigObj_filterBits[iobj] & bit) == 0) {
        matched_bits = false;
        break;
      }
    }
    if (!matched_bits)
      continue;
    return true;
  }
  return false;
}

