#include "Tools/Tools/interface/HHJetsInterface.h"

#include <stdexcept>

// Constructor
HHJetsInterface::HHJetsInterface (std::string model_0, std::string model_1, int year, bool isUL, float btag_wp, float fatjet_bbtag_wp):
  HHbtagger_(std::array<std::string, 2> { {model_0, model_1} }), btag_wp_(btag_wp), fatjet_bbtag_wp_(fatjet_bbtag_wp)
{
  year_ = year;
  if (isUL) max_bjet_eta = 2.5;
}


// Destructor
HHJetsInterface::~HHJetsInterface() {}



output HHJetsInterface::GetHHJets(
    unsigned long long int event, int pairType,
    rfRVec Jet_pt, rfRVec Jet_eta, rfRVec Jet_phi, rfRVec Jet_mass,
    rcRVec Jet_puId, rcRVec Jet_jetId, rfRVec Jet_btagDeepFlavB,
    rfRVec FatJet_pt, rfRVec FatJet_eta, rfRVec FatJet_phi, rfRVec FatJet_mass,
    rfRVec FatJet_msoftdrop, rcRVec FatJet_jetId, rfRVec FatJet_particleNet_XbbVsQCD,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi, 
    rfRVec GenPart_pt, rfRVec GenPart_eta, rfRVec GenPart_phi, rfRVec GenPart_mass,
    int genXbb_GenPartIdx, int gen_b1_GenPartIdx, int gen_b2_GenPartIdx,
    bool doGenCutFlow
    )
{
  auto deltaRGen = [&](TLorentzVector const& jet_tlv, int genPartIdx) {
    TLorentzVector gen_tlv;
    gen_tlv.SetPtEtaPhiM(GenPart_pt[genPartIdx], GenPart_eta[genPartIdx], GenPart_phi[genPartIdx], GenPart_mass[genPartIdx]);
    return jet_tlv.DeltaR(gen_tlv);
  };
  auto resolvedJetGenMatched = [&](TLorentzVector const& jet_tlv) {
    return deltaRGen(jet_tlv, gen_b1_GenPartIdx) < 0.4 || deltaRGen(jet_tlv, gen_b2_GenPartIdx) < 0.4;
  };
  auto resolvedJetGenMatched_idx = [&](int jet_idx) {
    auto jet_tlv = TLorentzVector();
    jet_tlv.SetPtEtaPhiM(Jet_pt[jet_idx], Jet_eta[jet_idx], Jet_phi[jet_idx], Jet_mass[jet_idx]);
    return resolvedJetGenMatched(jet_tlv);
  };
  auto boostedJetGenMatched = [&](TLorentzVector const& jet_tlv) {
    return deltaRGen(jet_tlv, genXbb_GenPartIdx) < 0.2; // Could even go down to 0.15
  };
  Jets_cutflow_output cutflow_output{}; // default-initialize so everything is filled with false (very important)
  std::vector<JetsFailReason> resolvedFailReasons; // List of failReason for genmatched AK4 jets. Will be put in cutflow_output later
  if (!(genXbb_GenPartIdx>=0 && gen_b1_GenPartIdx >= 0 && gen_b2_GenPartIdx >= 0)) doGenCutFlow = false;
  cutflow_output.fatjet_failReason.Reco = true; // Willl get overwritten if we find a FatJet genmatched

  std::vector<float> all_HHbtag_scores;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    all_HHbtag_scores.push_back(-999.);
  }

  int bjet1_idx = -1;
  int bjet2_idx = -1;
  int vbfjet1_idx = -1;
  int vbfjet2_idx = -1;
  int fatjet_idx = -1;
  std::vector <int> ctjet_indexes, fwjet_indexes;

  if (pairType < 0) {
    return output({all_HHbtag_scores, bjet1_idx, bjet2_idx, vbfjet1_idx, vbfjet2_idx,
      ctjet_indexes, fwjet_indexes, fatjet_idx, cutflow_output});
  }

  auto dau1_tlv = TLorentzVector();
  auto dau2_tlv = TLorentzVector();
  auto met_tlv = TLorentzVector();
  dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
  dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
  met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);
  
  // --- Resolved
  std::vector <jet_idx_btag> jet_indexes;
  std::vector <int> all_jet_indexes;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    JetsFailReason failReason;
    // JetID : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID    https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
    // PUId >=1 is loose
    if (Jet_puId[ijet] < 1 && Jet_pt[ijet] <= 50) failReason.JetPUID = true;
    if (Jet_jetId[ijet] < 6) // JetId == 6 (2017/2018) or 7 (2016) means pass tight and tightLepVeto IDs (tight jetId is mandatory, lepVeto is recommended)
      failReason.JetID = true;
    auto jet_tlv = TLorentzVector();
    jet_tlv.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
    if (jet_tlv.DeltaR(dau1_tlv) < 0.5 || jet_tlv.DeltaR(dau2_tlv) < 0.5)
      failReason.DeltaRDau = true;
    if (Jet_pt[ijet] <= 20) failReason.Pt = true;

    if (failReason.pass() && fabs(Jet_eta[ijet]) < 4.7)
      all_jet_indexes.push_back(ijet);
    
    if (fabs(Jet_eta[ijet]) >= max_bjet_eta) failReason.Eta = true;
    if (failReason.pass())
      jet_indexes.push_back(jet_idx_btag({(int) ijet, Jet_btagDeepFlavB[ijet]}));
    
    if (doGenCutFlow && resolvedJetGenMatched(jet_tlv))
      resolvedFailReasons.push_back(failReason);
  }
  if (jet_indexes.size() >= 2) {
    std::stable_sort(jet_indexes.begin(), jet_indexes.end(), jetSort);

    auto htt_tlv = dau1_tlv + dau2_tlv;
    std::vector <float> HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_;
    std::vector <float> HHbtag_jet_htt_deta_, HHbtag_jet_htt_dphi_, HHbtag_jet_deepFlavour_;

    for (auto &jet : jet_indexes) {
      auto jet_tlv = TLorentzVector();
      jet_tlv.SetPtEtaPhiM(Jet_pt[jet.idx], Jet_eta[jet.idx], Jet_phi[jet.idx], Jet_mass[jet.idx]);
      HHbtag_jet_pt_.push_back(jet_tlv.Pt());
      HHbtag_jet_eta_.push_back(jet_tlv.Eta());
      HHbtag_rel_jet_M_pt_.push_back(jet_tlv.M() / jet_tlv.Pt());
      HHbtag_rel_jet_E_pt_.push_back(jet_tlv.E() / jet_tlv.Pt());
      HHbtag_jet_htt_deta_.push_back(htt_tlv.Eta() - jet_tlv.Eta());
      HHbtag_jet_htt_dphi_.push_back(ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, jet_tlv));
      HHbtag_jet_deepFlavour_.push_back(jet.btag);
    }

    auto HHbtag_htt_met_dphi_ = (float) ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, met_tlv);
    auto HHbtag_htt_scalar_pt_ = (float) (dau1_tlv.Pt() + dau2_tlv.Pt());
    auto HHbtag_rel_met_pt_htt_pt_ = (float) met_tlv.Pt() / HHbtag_htt_scalar_pt_;
    auto HHbtag_htt_pt_ = (float) htt_tlv.Pt();
    auto HHbtag_htt_eta_ = (float) htt_tlv.Eta();
    int HHbtag_channel_ = -1;
    // Old HHbtag channel assignement. For Run3 this was changed to match analysis definition, but for run2 etau/mutau are swapped
    if (pairType == 0)
      HHbtag_channel_ = 1;
    else if (pairType == 1)
      HHbtag_channel_ = 0;
    else if (pairType == 2)
      HHbtag_channel_ = 2; // updated, used to be 1!

    auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_,
      HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_,
      HHbtag_jet_deepFlavour_, HHbtag_jet_htt_dphi_, year_, HHbtag_channel_,
      HHbtag_htt_pt_, HHbtag_htt_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, event);

    for (size_t ijet = 0; ijet < jet_indexes.size(); ijet++) {
      all_HHbtag_scores[jet_indexes[ijet].idx] = HHbtag_scores[ijet];
    }
    std::vector <jet_idx_btag> jet_indexes_hhbtag;
    for (size_t ijet = 0; ijet < jet_indexes.size(); ijet++) {
      jet_indexes_hhbtag.push_back(jet_idx_btag({jet_indexes[ijet].idx, HHbtag_scores[ijet]}));
    }
    std::stable_sort(jet_indexes_hhbtag.begin(), jet_indexes_hhbtag.end(), jetSort);
    bjet1_idx = jet_indexes_hhbtag[0].idx;
    bjet2_idx = jet_indexes_hhbtag[1].idx;

    if (Jet_pt[bjet1_idx] < Jet_pt[bjet2_idx]) {
      auto aux = bjet1_idx;
      bjet1_idx = bjet2_idx;
      bjet2_idx = aux;
    }
    if (doGenCutFlow) {
      if (!resolvedJetGenMatched_idx(bjet1_idx) || !resolvedJetGenMatched_idx(bjet2_idx)) {
        cutflow_output.wrongJet = true;
      }
    }
    /* // Disable VBF
    if (all_jet_indexes.size() >= 4) { // 2 bjets + 2 vbf jets
      std::vector <jet_pair_mass> vbfjet_indexes;
      for (size_t ijet = 0; ijet < all_jet_indexes.size(); ijet++) {
        auto jet1_index = all_jet_indexes[ijet];
        if (jet1_index == (int) bjet1_idx || jet1_index == (int) bjet2_idx)
          continue;
        for (size_t jjet = ijet + 1; jjet < all_jet_indexes.size(); jjet++) {
          auto jet2_index = all_jet_indexes[jjet];
          if (jet2_index == (int) bjet1_idx || jet2_index == (int) bjet2_idx)
            continue;
          if (Jet_pt[jet1_index] < 30 || Jet_pt[jet2_index] < 30)
            continue;
          auto jet1_tlv = TLorentzVector();
          auto jet2_tlv = TLorentzVector();
          jet1_tlv.SetPtEtaPhiM(Jet_pt[jet1_index], Jet_eta[jet1_index],
            Jet_phi[jet1_index], Jet_mass[jet1_index]);
          jet2_tlv.SetPtEtaPhiM(Jet_pt[jet2_index], Jet_eta[jet2_index],
            Jet_phi[jet2_index], Jet_mass[jet2_index]);
          auto jj_tlv = jet1_tlv + jet2_tlv;
          vbfjet_indexes.push_back(jet_pair_mass({jet1_index, jet2_index, (float) jj_tlv.M()}));
        }
      }
      if (vbfjet_indexes.size() > 0) {
        std::stable_sort(vbfjet_indexes.begin(), vbfjet_indexes.end(), jetPairSort);
        vbfjet1_idx = vbfjet_indexes[0].idx1;
        vbfjet2_idx = vbfjet_indexes[0].idx2;

        if (Jet_pt[vbfjet1_idx] < Jet_pt[vbfjet2_idx]) {
          auto aux = vbfjet1_idx;
          vbfjet1_idx = vbfjet2_idx;
          vbfjet2_idx = aux;
        }
      }
    }*/


    // additional central and forward jets
    for (auto & ijet : all_jet_indexes) {
      if ((int) ijet == bjet1_idx || (int) ijet == bjet2_idx
          || (int) ijet == vbfjet1_idx || (int) ijet == vbfjet2_idx)
        continue;
      if (fabs(Jet_eta[ijet]) < max_bjet_eta)
        ctjet_indexes.push_back(ijet);
      else if (fabs(Jet_eta[ijet]) < 4.7 && Jet_pt[ijet] > 30)
        fwjet_indexes.push_back(ijet);
    }
  } // jet_indexes.size() >= 2

  if (doGenCutFlow) { // Put the FailReason for the 2 jets. In case there are more than 2 gen-matched jets we choose the 2 most successful ones
    std::stable_sort(resolvedFailReasons.begin(), resolvedFailReasons.end(), [](JetsFailReason const& a, JetsFailReason const& b){ return a.countFails() < b.countFails();});
    if (resolvedFailReasons.size() >= 1)
      cutflow_output.jet1_failReason = resolvedFailReasons[0];
    if (resolvedFailReasons.size() >= 2)
      cutflow_output.jet2_failReason = resolvedFailReasons[1];
  }
  
  // -- Boosted :  Looking for AK8 jets for boosted
  // new definiton of boosted only requiring 1 AK8 jet (no subjets match)
  std::vector <jet_idx_btag> fatjet_indexes;
  for (size_t ifatjet = 0; ifatjet < FatJet_pt.size(); ifatjet++) {
    JetsFailReason failReason;
    if (FatJet_pt[ifatjet] < 250) failReason.Pt = true; // Probably this could be reduced to 200 ???
    if (FatJet_jetId[ifatjet] < 6) failReason.JetID = true; // FatJet should use the same JetId as AK4 jets
    if (FatJet_eta[ifatjet] >= 2.4) failReason.Eta = true;
    auto fatjet_tlv = TLorentzVector();
    fatjet_tlv.SetPtEtaPhiM(FatJet_pt[ifatjet], FatJet_eta[ifatjet],
      FatJet_phi[ifatjet], FatJet_mass[ifatjet]);
    if (fatjet_tlv.DeltaR(dau1_tlv) < 0.8) failReason.DeltaRDau = true;
    if (fatjet_tlv.DeltaR(dau2_tlv) < 0.8) failReason.DeltaRDau = true;
    if (FatJet_msoftdrop.at(ifatjet) < 30) failReason.SoftDrop = true;
    // failReason.print(std::cout);
    if (failReason.pass()) {
      fatjet_indexes.push_back(jet_idx_btag({(int) ifatjet, FatJet_particleNet_XbbVsQCD.at(ifatjet)}));
    }
    if (doGenCutFlow && boostedJetGenMatched(fatjet_tlv))
      cutflow_output.fatjet_failReason = failReason;
  }
  
  if (fatjet_indexes.size() != 0) {
    std::stable_sort(fatjet_indexes.begin(), fatjet_indexes.end(), jetSort);
    fatjet_idx = fatjet_indexes[0].idx;
    if (doGenCutFlow) {
      auto fatjet_tlv = TLorentzVector();
      fatjet_tlv.SetPtEtaPhiM(FatJet_pt[fatjet_idx], FatJet_eta[fatjet_idx],
        FatJet_phi[fatjet_idx], FatJet_mass[fatjet_idx]);
      if (!boostedJetGenMatched(fatjet_tlv))
        cutflow_output.wrongFatJet = true;
    } 
  }

  /*
  // Logic for removing overlap between boosted & resolved
  if (priorityMode == JetCategoryPriorityMode::Res2b_Boosted_Res1b_noPNetFail) {
    // Priority order res2b -> boosted -> res1b
    if (bjet1_idx >= 0 && bjet2_idx >= 0 && (Jet_btagDeepFlavB[bjet1_idx] >= btag_wp_ && Jet_btagDeepFlavB[bjet2_idx] >= btag_wp_)) {
      // resolved-2b
      jetCategory = JetCategory::Res_2b;
    } else if (fatjet_idx >= 0) {
      if (FatJet_particleNet_XbbVsQCD.at(fatjet_idx) >= fatjet_bbtag_wp_) 
        jetCategory = JetCategory::Boosted_bb; // boosted-bb
      else
        jetCategory = JetCategory::Boosted_failedPNet;
    } else if (bjet1_idx >= 0 && bjet2_idx >= 0 && (Jet_btagDeepFlavB[bjet1_idx] >= btag_wp_ || Jet_btagDeepFlavB[bjet2_idx] >= btag_wp_)) {
      // resolved-1b 
      jetCategory = JetCategory::Res_1b;
    } else {
      jetCategory = JetCategory::None;
    }
  }
  else if (priorityMode == JetCategoryPriorityMode::Boosted_Res2b_Res1b_noPNetFail) {
    if (fatjet_idx >= 0) {
      if (FatJet_particleNet_XbbVsQCD.at(fatjet_idx) >= fatjet_bbtag_wp_) 
        jetCategory = JetCategory::Boosted_bb; // boosted-bb
      else
        jetCategory = JetCategory::Boosted_failedPNet;
    }
    else if (bjet1_idx >= 0 && bjet2_idx >= 0 && (Jet_btagDeepFlavB[bjet1_idx] >= btag_wp_ && Jet_btagDeepFlavB[bjet2_idx] >= btag_wp_)) {
      // resolved-2b
      jetCategory = JetCategory::Res_2b;
    } else if (bjet1_idx >= 0 && bjet2_idx >= 0 && (Jet_btagDeepFlavB[bjet1_idx] >= btag_wp_ || Jet_btagDeepFlavB[bjet2_idx] >= btag_wp_)) {
      // resolved-1b 
      jetCategory = JetCategory::Res_1b;
    } else {
      jetCategory = JetCategory::None;
    }
  } else {
    throw std::invalid_argument("HHJets : Wrong JetCategoryPriorityMode");
  }
  */

  return output({all_HHbtag_scores, bjet1_idx, bjet2_idx, vbfjet1_idx, vbfjet2_idx,
    ctjet_indexes, fwjet_indexes, fatjet_idx, cutflow_output});
}


// GetScore
std::vector<float> HHJetsInterface::GetScore(
  std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
  std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_deepFlavour_,
  std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
  float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_)
{
  // Get HHbtag score
  auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_,
      HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_, HHbtag_jet_deepFlavour_, HHbtag_jet_htt_dphi_,
      HHbtag_year_, HHbtag_channel_, HHbtag_tauH_pt_, HHbtag_tauH_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, HHbtag_evt_);

  // Store HHbtag scores in a map<jet_idx,HHbtag_score>
  std::vector<float> jets_and_HHbtag;
  for (unsigned int i = 0; i < HHbtag_jet_pt_.size(); i++)
  {
    jets_and_HHbtag.push_back(HHbtag_scores.at(i));
  }

  return jets_and_HHbtag;
}
