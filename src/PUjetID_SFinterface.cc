#include "Tools/Tools/interface/PUjetID_SFinterface.h"

// Constructor
PUjetID_SFinterface::PUjetID_SFinterface (int year, std::string folder_path) {
  year_ = year;
  TFile* f_eff = TFile::Open(
    (folder_path + "h2_eff_mc_" + std::to_string(year) + "_L.root").c_str());
  TFile* f_eff_sf = TFile::Open(
    (folder_path + "h2_eff_sf_" + std::to_string(year) + "_L.root").c_str());
  TFile* f_mistag = TFile::Open(
    (folder_path + "h2_mistag_mc_" + std::to_string(year) + "_L.root").c_str());
  TFile* f_mistag_sf = TFile::Open(
    (folder_path + "h2_mistag_sf_" + std::to_string(year) + "_L.root").c_str());
  TFile* f_sf_err = TFile::Open(
    (folder_path + "scalefactorsPUID_81Xtraining.root").c_str());

  h_eff_ = (TH2F*) f_eff->Get(("h2_eff_mc" + std::to_string(year) + "_L").c_str());
  h_eff_sf_ = (TH2F*) f_eff_sf->Get(("h2_eff_sf" + std::to_string(year) + "_L").c_str());
  h_eff_sf_err_ = (TH2F*) f_sf_err->Get(("h2_eff_sf" + std::to_string(year) + "_L_Systuncty").c_str());
  h_mistag_ = (TH2F*) f_mistag->Get(("h2_mistag_mc" + std::to_string(year) + "_L").c_str());
  h_mistag_sf_ = (TH2F*) f_mistag_sf->Get(("h2_mistag_sf" + std::to_string(year) + "_L").c_str());
  h_mistag_sf_err_ = (TH2F*) f_sf_err->Get(
    ("h2_mistag_sf" + std::to_string(year) + "_L_Systuncty").c_str());
}


std::vector<double> PUjetID_SFinterface::get_pu_weights(
    fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass, iRVec Jet_jetId, iRVec Jet_puId,
    fRVec GenJet_pt, fRVec GenJet_eta, fRVec GenJet_phi, fRVec GenJet_mass,
    double dau1_pt, double dau1_eta, double dau1_phi, double dau1_mass,
    double dau2_pt, double dau2_eta, double dau2_phi, double dau2_mass)
{
  double P_MC = 1.;
  double P_DATA = 1.;
  double P_DATA_up = 1.;
  double P_DATA_down = 1.;
  double P_DATA_effic_up = 1.;
  double P_DATA_effic_down = 1.;
  double P_DATA_mistag_up = 1.;
  double P_DATA_mistag_down = 1.;
  double P_DATA_effic_eta_s2p5_up = 1.;
  double P_DATA_effic_eta_s2p5_down = 1.;
  double P_DATA_effic_eta_l2p5_up = 1.;
  double P_DATA_effic_eta_l2p5_down = 1.;
  double P_DATA_mistag_eta_s2p5_up = 1.;
  double P_DATA_mistag_eta_s2p5_down = 1.;
  double P_DATA_mistag_eta_l2p5_up = 1.;
  double P_DATA_mistag_eta_l2p5_down = 1.;

  auto dau1_tlv = TLorentzVector();
  auto dau2_tlv = TLorentzVector();
  dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
  dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);

  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    if (Jet_jetId[ijet] < 2)
      continue;
    auto jet_tlv = TLorentzVector();
    jet_tlv.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
    if (jet_tlv.Pt() < 20 || jet_tlv.Pt() > 50 || fabs(jet_tlv.Eta()) > 4.7)
      continue;
    if (jet_tlv.DeltaR(dau1_tlv) < 0.5 || jet_tlv.DeltaR(dau2_tlv) < 0.5)
      continue;

    // noisy jet removal for 2017
    if (year_ == 2017 && fabs(jet_tlv.Eta()) > 2.65 && fabs(jet_tlv.Eta()) < 3.139)
      continue;

    bool isRealJet = false;
    TLorentzVector genjet_tlv;
    for (unsigned int igj = 0; igj < GenJet_pt.size(); igj++) {
      genjet_tlv.SetPtEtaPhiM(GenJet_pt.at(igj), GenJet_eta.at(igj), GenJet_phi.at(igj), GenJet_mass.at(igj));
      if (jet_tlv.DeltaR(genjet_tlv) < 0.4) {
        isRealJet = true;
        break;
      }
    }

    auto eff_and_sf = get_eff_sf_and_error(isRealJet, (double) jet_tlv.Pt(), (double) jet_tlv.Eta());
    auto eff = eff_and_sf[0];
    auto sf = eff_and_sf[1];
    auto sf_err = eff_and_sf[2];
    auto sf_up = std::min(std::max(0., sf + sf_err), 5.);
    auto sf_down = std::min(std::max(0., sf - sf_err), 5.);
    if (Jet_puId[ijet] >= 4 || Jet_pt[ijet] > 50) {
      P_MC *= eff;
      P_DATA *= sf * eff;
      P_DATA_up *= sf_up * eff;
      P_DATA_down *= sf_down * eff;
      P_DATA_effic_up *= sf_up * eff;
      P_DATA_effic_down *= sf_down * eff;
      P_DATA_mistag_up *= sf * eff;  // true jet --> use nominal SF for mistag
      P_DATA_mistag_down *= sf * eff;  // true jet --> use nominal SF for mistag

      if (abs(jet_tlv.Eta()) <= 2.5) {
        P_DATA_effic_eta_s2p5_up *= sf_up * eff;
        P_DATA_effic_eta_s2p5_down *= sf_down * eff;
        P_DATA_effic_eta_l2p5_up *= sf * eff;
        P_DATA_effic_eta_l2p5_down *= sf * eff;
        P_DATA_mistag_eta_s2p5_up *= sf * eff;
        P_DATA_mistag_eta_s2p5_down *= sf * eff;
        P_DATA_mistag_eta_l2p5_up *= sf * eff;
        P_DATA_mistag_eta_l2p5_down *= sf * eff;
      } else {
        P_DATA_effic_eta_s2p5_up *= sf * eff;
        P_DATA_effic_eta_s2p5_down *= sf * eff;
        P_DATA_effic_eta_l2p5_up *= sf_up * eff;
        P_DATA_effic_eta_l2p5_down *= sf_down * eff;
        P_DATA_mistag_eta_s2p5_up *= sf * eff;
        P_DATA_mistag_eta_s2p5_down *= sf * eff;
        P_DATA_mistag_eta_l2p5_up *= sf * eff;
        P_DATA_mistag_eta_l2p5_down *= sf * eff;
      }
    } else {
      P_MC *= (1. - eff);
      P_DATA *= (1. - sf * eff);
      P_DATA_up *= (1. - sf_up * eff);
      P_DATA_down *= (1. - sf_down * eff);
      P_DATA_effic_up *= (1. - sf * eff);  // fake jet --> use nominal SF for effic
      P_DATA_effic_down *= (1. - sf * eff);  // fake jet --> use nominal SF for effic
      P_DATA_mistag_up *= (1. - sf_up * eff);
      P_DATA_mistag_down *= (1. - sf_down * eff);

      if (fabs(jet_tlv.Eta()) <= 2.5) {
        P_DATA_effic_eta_s2p5_up *= (1 - sf * eff);
        P_DATA_effic_eta_s2p5_down *= (1 - sf * eff);
        P_DATA_effic_eta_l2p5_up *= (1 - sf * eff);
        P_DATA_effic_eta_l2p5_down *= (1 - sf * eff);
        P_DATA_mistag_eta_s2p5_up *= (1 - sf_up * eff);
        P_DATA_mistag_eta_s2p5_down *= (1 - sf_down * eff);
        P_DATA_mistag_eta_l2p5_up *= (1 - sf * eff);
        P_DATA_mistag_eta_l2p5_down *= (1 - sf * eff);
      } else {
        P_DATA_effic_eta_s2p5_up *= (1 - sf * eff);
        P_DATA_effic_eta_s2p5_down *= (1 - sf * eff);
        P_DATA_effic_eta_l2p5_up *= (1 - sf * eff);
        P_DATA_effic_eta_l2p5_down *= (1 - sf * eff);
        P_DATA_mistag_eta_s2p5_up *= (1 - sf * eff);
        P_DATA_mistag_eta_s2p5_down *= (1 - sf * eff);
        P_DATA_mistag_eta_l2p5_up *= (1 - sf_up * eff);
        P_DATA_mistag_eta_l2p5_down *= (1 - sf_down * eff);
      }
    }
  }

  return {P_DATA / P_MC, P_DATA_up / P_MC, P_DATA_down / P_MC,
      P_DATA_effic_up / P_MC, P_DATA_effic_down / P_MC,
      P_DATA_mistag_up / P_MC, P_DATA_mistag_down / P_MC,
      P_DATA_effic_eta_s2p5_up / P_MC, P_DATA_effic_eta_s2p5_down / P_MC, 
      P_DATA_effic_eta_l2p5_up / P_MC, P_DATA_effic_eta_l2p5_down / P_MC,
      P_DATA_mistag_eta_s2p5_up / P_MC, P_DATA_mistag_eta_s2p5_down / P_MC,
      P_DATA_mistag_eta_l2p5_up / P_MC, P_DATA_mistag_eta_l2p5_down / P_MC};

}


// Destructor
PUjetID_SFinterface::~PUjetID_SFinterface() {}
