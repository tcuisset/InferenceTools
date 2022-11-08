#include "Tools/Tools/interface/HHMulticlassInterface.h"

// Constructor
HHMulticlassInterface:: HHMulticlassInterface (int year, std::vector<std::pair<std::string, std::string>> modelSpecs)
 : mci_(year, modelSpecs)
{
  // Store target lambdas
  year_ = year;
}

// Destructor
HHMulticlassInterface::~HHMulticlassInterface() {}

std::vector<std::vector<float>> HHMulticlassInterface::GetPredictionsWithInputs(
  int EventNumber, int pairType,
  fRVec jet_pt, fRVec jet_eta, fRVec jet_phi, fRVec jet_mass,
  fRVec jet_btagDeepFlavB, fRVec jet_btagDeepFlavCvL, fRVec jet_btagDeepFlavCvB, fRVec jet_HHbtag,
  int dau1_index, int dau2_index,
  int bjet1_index, int bjet2_index, int vbfjet1_index, int vbfjet2_index,
  iRVec ctjet_indexes, iRVec fwjet_indexes,
  fRVec muon_pt, fRVec muon_eta, fRVec muon_phi, fRVec muon_mass,
  fRVec electron_pt, fRVec electron_eta, fRVec electron_phi, fRVec electron_mass,
  fRVec tau_pt, fRVec tau_eta, fRVec tau_phi, fRVec tau_mass,
  float met_pt, float met_phi,
  float htt_sv_pt, float htt_sv_eta, float htt_sv_phi, float htt_sv_mass
) 
{
  mci_.clearInputs();
  float dau1_pt = -1, dau1_eta = -1, dau1_phi = -1, dau1_mass = -1;

  if (pairType == 0) {
    dau1_pt = muon_pt.at(dau1_index);
    dau1_eta = muon_eta.at(dau1_index);
    dau1_phi = muon_phi.at(dau1_index);
    dau1_mass = muon_mass.at(dau1_index);
  } else if (pairType == 1) {
    dau1_pt = electron_pt.at(dau1_index);
    dau1_eta = electron_eta.at(dau1_index);
    dau1_phi = electron_phi.at(dau1_index);
    dau1_mass = electron_mass.at(dau1_index);
  } else if (pairType == 2) {
    dau1_pt = tau_pt.at(dau1_index);
    dau1_eta = tau_eta.at(dau1_index);
    dau1_phi = tau_phi.at(dau1_index);
    dau1_mass = tau_mass.at(dau1_index);
  }
  float dau2_pt = tau_pt.at(dau2_index);
  float dau2_eta = tau_eta.at(dau2_index);
  float dau2_phi = tau_phi.at(dau2_index);
  float dau2_mass = tau_mass.at(dau2_index);

  auto dau1_tlv = TLorentzVector();
  auto dau2_tlv = TLorentzVector();
  auto bjet1_tlv = TLorentzVector();
  auto bjet2_tlv = TLorentzVector();

  dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
  dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
  bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
    jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
  bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
    jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));

  float vbfjet1_pt = -999, vbfjet1_eta = -999, vbfjet1_phi = -999, vbfjet1_e = -999;
  float vbfjet2_pt = -999, vbfjet2_eta = -999, vbfjet2_phi = -999, vbfjet2_e = -999;
  if (vbfjet1_index >= 0) {
    auto vbfjet1_tlv = TLorentzVector();
    auto vbfjet2_tlv = TLorentzVector();
    vbfjet1_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet1_index), jet_eta.at(vbfjet1_index),
      jet_phi.at(vbfjet1_index), jet_mass.at(vbfjet1_index));
    vbfjet2_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet2_index), jet_eta.at(vbfjet2_index),
      jet_phi.at(vbfjet2_index), jet_mass.at(vbfjet2_index));

    vbfjet1_pt = vbfjet1_tlv.Pt();
    vbfjet1_eta = vbfjet1_tlv.Eta();
    vbfjet1_phi = vbfjet1_tlv.Phi();
    vbfjet1_e = vbfjet1_tlv.E();

    vbfjet2_pt = vbfjet2_tlv.Pt();
    vbfjet2_eta = vbfjet2_tlv.Eta();
    vbfjet2_phi = vbfjet2_tlv.Phi();
    vbfjet2_e = vbfjet2_tlv.E();
  }

  auto met_tlv = TLorentzVector();
  met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

  float tauh_sv_pt = -999, tauh_sv_eta = -999, tauh_sv_phi = -999, tauh_sv_e = -999, tauh_sv_ez = -999;
  if (htt_sv_mass > 0) {
    auto htt_svfit_tlv = TLorentzVector();
    htt_svfit_tlv.SetPtEtaPhiM(htt_sv_pt, htt_sv_eta, htt_sv_phi, htt_sv_mass);
    tauh_sv_pt = htt_svfit_tlv.Pt();
    tauh_sv_eta = htt_svfit_tlv.Eta();
    tauh_sv_phi = htt_svfit_tlv.Phi();
    tauh_sv_e = htt_svfit_tlv.E();
    tauh_sv_ez = std::pow(std::pow(htt_svfit_tlv.Pz(), 2) + std::pow(htt_svfit_tlv.M(), 2), 0.5);
  }

  float deepFlav1 = -1., deepFlav2 = -1., deepFlav_vbf1 = -1., deepFlav_vbf2 = -1.,
    CvsL_b1 = -1., CvsL_b2 = -1.,
    CvsL_vbf1 = -1., CvsL_vbf2 = -1., CvsB_b1 = -1., CvsB_b2 = -1.,
    CvsB_vbf1 = -1., CvsB_vbf2 = -1., HHbtag_b1 = -1., HHbtag_b2 = -1.,
    HHbtag_vbf1 = -1., HHbtag_vbf2 = -1.;
  deepFlav1 = jet_btagDeepFlavB.at(bjet1_index);
  deepFlav2 = jet_btagDeepFlavB.at(bjet2_index);
  CvsL_b1 = jet_btagDeepFlavCvL.at(bjet1_index);
  CvsL_b2 = jet_btagDeepFlavCvL.at(bjet2_index);
  CvsB_b1 = jet_btagDeepFlavCvB.at(bjet1_index);
  CvsB_b2 = jet_btagDeepFlavCvB.at(bjet2_index);
  HHbtag_b1 = jet_HHbtag.at(bjet1_index);
  HHbtag_b2 = jet_HHbtag.at(bjet2_index);
  if (vbfjet1_index >= 0) {
    CvsL_vbf1 = jet_btagDeepFlavCvL.at(vbfjet1_index);
    CvsL_vbf2 = jet_btagDeepFlavCvL.at(vbfjet2_index);
    CvsB_vbf1 = jet_btagDeepFlavCvB.at(vbfjet1_index);
    CvsB_vbf2 = jet_btagDeepFlavCvB.at(vbfjet2_index);
    HHbtag_vbf1 = jet_HHbtag.at(vbfjet1_index);
    HHbtag_vbf2 = jet_HHbtag.at(vbfjet2_index);
  }

  std::vector<float> ctjet_pt(3, -999), ctjet_eta(3, -999), ctjet_phi(3, -999), ctjet_e(3, -999),
    ctjet_deepflavor_b(3, -999), ctjet_hhbtag(3, -999);
  std::vector<float> fwjet_pt(2, -999), fwjet_eta(2, -999), fwjet_phi(2, -999), fwjet_e(2, -999);

  for (size_t ictjet = 0; ictjet < 3; ictjet++) {
    if (ctjet_indexes.size() <= ictjet)
      break;
    auto aux_tlv = TLorentzVector();
    aux_tlv.SetPtEtaPhiM(jet_pt.at(ctjet_indexes[ictjet]), jet_eta.at(ctjet_indexes[ictjet]),
      jet_phi.at(ctjet_indexes[ictjet]), jet_mass.at(ctjet_indexes[ictjet]));
    ctjet_pt[ictjet] = aux_tlv.Pt();
    ctjet_eta[ictjet] = aux_tlv.Eta();
    ctjet_phi[ictjet] = aux_tlv.Phi();
    ctjet_e[ictjet] = aux_tlv.E();
    ctjet_deepflavor_b[ictjet] = jet_btagDeepFlavB.at(ctjet_indexes[ictjet]);
    ctjet_hhbtag[ictjet] = jet_HHbtag.at(ctjet_indexes[ictjet]);
  }
  
  for (size_t ifwjet = 0; ifwjet < 2; ifwjet++) {
    if (fwjet_indexes.size() <= ifwjet)
      break;
    auto aux_tlv = TLorentzVector();
    aux_tlv.SetPtEtaPhiM(jet_pt.at(fwjet_indexes[ifwjet]), jet_eta.at(fwjet_indexes[ifwjet]),
      jet_phi.at(fwjet_indexes[ifwjet]), jet_mass.at(fwjet_indexes[ifwjet]));
    fwjet_pt[ifwjet] = aux_tlv.Pt();
    fwjet_eta[ifwjet] = aux_tlv.Eta();
    fwjet_phi[ifwjet] = aux_tlv.Phi();
    fwjet_e[ifwjet] = aux_tlv.E();
  }

  for (size_t j = 0; j < mci_.getNumberOfModels(); j++)
  {
    mci_.setInputs(j,
      {
        {"is_mutau", pairType == 0 ? 1 : 0},
        {"is_etau", pairType == 1 ? 1 : 0},
        {"is_tautau", pairType == 2 ? 1 : 0},

        {"is_2016", year_ == 2016 ? 1 : 0},
        {"is_2017", year_ == 2017 ? 1 : 0},
        {"is_2018", year_ == 2018 ? 1 : 0},

        {"bjet1_pt", bjet1_tlv.Pt()},
        {"bjet1_eta", bjet1_tlv.Eta()},
        {"bjet1_phi", bjet1_tlv.Phi()},
        {"bjet1_e", bjet1_tlv.E()},
        {"bjet1_deepflavor_b", deepFlav1},
        {"bjet1_hhbtag", HHbtag_b1},

        {"bjet2_pt", bjet2_tlv.Pt()},
        {"bjet2_eta", bjet2_tlv.Eta()},
        {"bjet2_phi", bjet2_tlv.Phi()},
        {"bjet2_e", bjet2_tlv.E()},
        {"bjet2_deepflavor_b", deepFlav2},
        {"bjet2_hhbtag", HHbtag_b2},

        {"ctjet1_pt", ctjet_pt[0]},
        {"ctjet1_eta", ctjet_eta[0]},
        {"ctjet1_phi", ctjet_phi[0]},
        {"ctjet1_e", ctjet_e[0]},
        {"ctjet1_deepflavor_b", ctjet_deepflavor_b[0]},
        {"ctjet1_hhbtag", ctjet_hhbtag[0]},

        {"ctjet2_pt", ctjet_pt[1]},
        {"ctjet2_eta", ctjet_eta[1]},
        {"ctjet2_phi", ctjet_phi[1]},
        {"ctjet2_e", ctjet_e[1]},
        {"ctjet2_deepflavor_b", ctjet_deepflavor_b[1]},
        {"ctjet2_hhbtag", ctjet_hhbtag[1]},

        {"ctjet3_pt", ctjet_pt[2]},
        {"ctjet3_eta", ctjet_eta[2]},
        {"ctjet3_phi", ctjet_phi[2]},
        {"ctjet3_e", ctjet_e[2]},
        {"ctjet3_deepflavor_b", ctjet_deepflavor_b[2]},
        {"ctjet3_hhbtag", ctjet_hhbtag[2]},

        {"fwjet1_pt", fwjet_pt[0]},
        {"fwjet1_eta", fwjet_eta[0]},
        {"fwjet1_phi", fwjet_phi[0]},
        {"fwjet1_e", fwjet_e[0]},

        {"fwjet2_pt", fwjet_pt[1]},
        {"fwjet2_eta", fwjet_eta[1]},
        {"fwjet2_phi", fwjet_phi[1]},
        {"fwjet2_e", fwjet_e[1]},

        {"vbfjet1_pt", vbfjet1_pt},
        {"vbfjet1_eta", vbfjet1_eta},
        {"vbfjet1_phi", vbfjet1_phi},
        {"vbfjet1_e", vbfjet1_e},
        {"vbfjet1_deepflavor_b", deepFlav_vbf1},
        {"vbfjet1_hhbtag", HHbtag_vbf1},

        {"vbfjet2_pt", vbfjet2_pt},
        {"vbfjet2_eta", vbfjet2_eta},
        {"vbfjet2_phi", vbfjet2_phi},
        {"vbfjet2_e", vbfjet2_e},
        {"vbfjet2_deepflavor_b", deepFlav_vbf2},
        {"vbfjet2_hhbtag", HHbtag_vbf2},

        {"bjet1_deepflavor_cvsb", CvsB_b1},
        {"bjet1_deepflavor_cvsl", CvsL_b1},

        {"bjet2_deepflavor_cvsb", CvsB_b2},
        {"bjet2_deepflavor_cvsl", CvsL_b2},

        {"vbfjet1_deepflavor_cvsb", CvsB_vbf1},
        {"vbfjet1_deepflavor_cvsl", CvsL_vbf1},

        {"vbfjet2_deepflavor_cvsb", CvsB_vbf2},
        {"vbfjet2_deepflavor_cvsl", CvsL_vbf2},

        {"lep1_pt", dau1_tlv.Pt()},
        {"lep1_eta", dau1_tlv.Eta()},
        {"lep1_phi", dau1_tlv.Phi()},
        {"lep1_e", dau1_tlv.E()},

        {"lep2_pt", dau2_tlv.Pt()},
        {"lep2_eta", dau2_tlv.Eta()},
        {"lep2_phi", dau2_tlv.Phi()},
        {"lep2_e", dau2_tlv.E()},

        {"met_pt", met_tlv.Pt()},
        {"met_phi", met_tlv.Phi()},

        {"bh_pt", (bjet1_tlv + bjet2_tlv).Pt()},
        {"bh_eta", (bjet1_tlv + bjet2_tlv).Eta()},
        {"bh_phi", (bjet1_tlv + bjet2_tlv).Phi()},
        {"bh_e", (bjet1_tlv + bjet2_tlv).E()},

        {"tauh_sv_pt", tauh_sv_pt},
        {"tauh_sv_eta", tauh_sv_eta},
        {"tauh_sv_phi", tauh_sv_phi},
        {"tauh_sv_e", tauh_sv_e},
        {"tauh_sv_ez", tauh_sv_ez}
      });
  }

  std::vector<std::vector<float>> output_scores;
  for (size_t j=0; j<mci_.getNumberOfModels(); j++) {
    auto model_scores = mci_.predict(EventNumber, j);
    std::vector<float> this_scores;
    for (size_t k = 0; k < model_scores.size(); k++) {
      this_scores.push_back(model_scores.at(k).second);
    }
    output_scores.push_back(this_scores);
  }
  return output_scores;
}

std::vector<std::string> HHMulticlassInterface::get_node_names(size_t imodel) {
    return mci_.getNodeNames(imodel);
}