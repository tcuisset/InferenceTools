#include "Tools/Tools/interface/Htt_trigSFinterface.h"

// Constructor
Htt_trigSFinterface:: Htt_trigSFinterface (
    int year, float mutau_pt_th1, float mutau_pt_th2, float etau_pt_th1, float etau_pt_th2,
    std::string eTrgSF_file, std::string eTrgSF_name, bool eTrgSF_bool, 
    std::string eTauTrgSF_file, std::string eTauTrgSF_name, bool eTauTrgSF_bool,
    std::string muTrgSF_file, std::string muTrgSF_name, bool muTrgSF_bool,
    std::string muTauTrgSF_file, std::string muTauTrgSF_name, bool muTauTrgSF_bool,
    std::string tauTrgSF_ditau_file, std::string tauTrgSF_mutau_file,
    std::string tauTrgSF_etau_file, std::string tauTrgSF_vbf_file, std::string jetTrgSF_vbf_file
): 
  tauTrgSF_ditau(tauTrgSF_ditau_file, "ditau", "Medium"),
  tauTrgSF_mutau(tauTrgSF_mutau_file, "mutau", "Medium"),
  tauTrgSF_etau(tauTrgSF_etau_file, "etau", "Medium"),
  tauTrgSF_vbf(tauTrgSF_vbf_file, "ditauvbf", "Medium")
{
  year_ = year;
  mutau_pt_th1_ = mutau_pt_th1;
  mutau_pt_th2_ = mutau_pt_th2;
  etau_pt_th1_ = etau_pt_th1;
  etau_pt_th2_ = etau_pt_th2;
  
  if (eTrgSF_name != "") eTrgSF.init_ScaleFactor(eTrgSF_file, eTrgSF_name);
  else eTrgSF.init_EG_ScaleFactor(eTrgSF_file, eTrgSF_bool);

  if (eTauTrgSF_name != "") eTauTrgSF.init_ScaleFactor(eTauTrgSF_file, eTauTrgSF_name);
  else eTauTrgSF.init_EG_ScaleFactor(eTauTrgSF_file, eTauTrgSF_bool);

  if (muTrgSF_name != "") muTrgSF.init_ScaleFactor(muTrgSF_file, muTrgSF_name);
  else muTrgSF.init_EG_ScaleFactor(muTrgSF_file, muTrgSF_bool);

  if (muTauTrgSF_name != "") muTauTrgSF.init_ScaleFactor(muTauTrgSF_file, muTauTrgSF_name);
  else muTauTrgSF.init_EG_ScaleFactor(muTauTrgSF_file, muTauTrgSF_bool);
      
  TFile* tf = TFile::Open(jetTrgSF_vbf_file.c_str());
  jetTrgSF_vbf = (TH3F*) tf->Get("SF_mjj_pT1_pT2");
}

std::vector<double> Htt_trigSFinterface::get_scale_factors(int pairType, int isVBFtrigger,
    int dau1_decayMode, float dau1_pt, float dau1_eta,
    int dau2_decayMode, float dau2_pt, float dau2_eta,
    float vbfjet1_pt, float vbfjet1_eta, float vbfjet1_phi, float vbfjet1_mass,
    float vbfjet2_pt, float vbfjet2_eta, float vbfjet2_phi, float vbfjet2_mass) {

  std::vector <int> decayModes = {0, 1, 10, 11};

  /////////////////////////////////////////////////////////////////////////////////////////
  // mutau
  /////////////////////////////////////////////////////////////////////////////////////////
  
  if (pairType == 0) {
    std::vector<double> trigSF_mu, trigSF_tauup, trigSF_taudown;
    double trigSF_single, trigSF_cross;
    if (fabs(dau2_eta) < 2.1) {
      int passSingle = (dau1_pt > mutau_pt_th1_) ? 1 : 0;
      int passCross = (dau2_pt > mutau_pt_th2_) ? 1 : 0;

      // lepton trigger
      auto Eff_L_Data_nom = muTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_L_MC_nom = muTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_L_Data_Err = muTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_L_MC_Err = muTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // std::cout << " --> muTrgSF : pt " << dau1_pt << " ; eta " << dau1_eta << " --> " << Eff_L_Data_nom << std::endl;

      std::vector <double> Eff_L_Data = {
        Eff_L_Data_nom - Eff_L_Data_Err, Eff_L_Data_nom, Eff_L_Data_nom + Eff_L_Data_Err};
      std::vector <double> Eff_L_MC = {
        Eff_L_MC_nom - Eff_L_MC_Err, Eff_L_MC_nom, Eff_L_MC_nom + Eff_L_MC_Err};

      // cross trigger
      // mu leg

      // std::cout << " --> muTauTrgSF : pt " << dau1_pt << " ; eta " << dau1_eta << " --> " << Eff_l_Data_nom << std::endl;

      auto Eff_l_Data_nom = muTauTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_l_MC_nom = muTauTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_l_Data_Err = muTauTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_l_MC_Err = muTauTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      std::vector <double> Eff_l_Data = {
        Eff_l_Data_nom - Eff_l_Data_Err, Eff_l_Data_nom, Eff_l_Data_nom + Eff_l_Data_Err};
      std::vector <double> Eff_l_MC = {
        Eff_l_MC_nom - Eff_l_MC_Err, Eff_l_MC_nom, Eff_l_MC_nom + Eff_l_MC_Err};

      // tau leg
      auto Eff_tau_Data = tauTrgSF_mutau.getEfficiencyData(dau2_pt, dau2_decayMode, 0);
      auto Eff_tau_MC = tauTrgSF_mutau.getEfficiencyMC(dau2_pt, dau2_decayMode, 0);

      std::cout << " --> tauTrgSF_mutau " << Eff_tau_Data << std::endl;

      std::vector <double> Eff_Data_mu, Eff_MC_mu;
      for (size_t i = 0; i< Eff_l_Data.size(); i++) {
        Eff_Data_mu.push_back(passSingle * Eff_L_Data[i] - passCross * passSingle * std::min(Eff_l_Data[i], Eff_L_Data[i])
          * Eff_tau_Data + passCross * Eff_l_Data[i] * Eff_tau_Data);
        Eff_MC_mu.push_back(passSingle * Eff_L_MC[i] - passCross * passSingle * std::min(Eff_l_MC[i], Eff_L_MC[i])
          * Eff_tau_MC + passCross * Eff_l_MC[i] * Eff_tau_MC);
        trigSF_mu.push_back(Eff_Data_mu[i] / Eff_MC_mu[i]);
      }

      std::vector<double> Eff_tau_Data_tauup(4, Eff_tau_Data);
      std::vector<double> Eff_tau_Data_taudown(4, Eff_tau_Data);
      std::vector<double> Eff_tau_MC_tauup(4, Eff_tau_MC);
      std::vector<double> Eff_tau_MC_taudown(4, Eff_tau_MC);
      std::vector<double> Eff_Data_tauup(4);
      std::vector<double> Eff_Data_taudown(4);
      std::vector<double> Eff_MC_tauup(4);
      std::vector<double> Eff_MC_taudown(4);

      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        if (decayModes[idm] == dau2_decayMode) {
          Eff_tau_Data_tauup[idm] = tauTrgSF_mutau.getEfficiencyData(dau2_pt, dau2_decayMode, 1);
          Eff_tau_Data_taudown[idm] = tauTrgSF_mutau.getEfficiencyData(dau2_pt, dau2_decayMode, -1);
          Eff_tau_MC_tauup[idm] = tauTrgSF_mutau.getEfficiencyMC(dau2_pt, dau2_decayMode, 1);
          Eff_tau_MC_taudown[idm] = tauTrgSF_mutau.getEfficiencyMC(dau2_pt, dau2_decayMode, -1);
        }
      }

      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        Eff_Data_tauup[idm] = passSingle * Eff_L_Data[1] - passCross * passSingle * std::min(Eff_l_Data[1], Eff_L_Data[1])
          * Eff_tau_Data_tauup[idm] + passCross * Eff_l_Data[1] * Eff_tau_Data_tauup[idm];
        Eff_Data_taudown[idm] = passSingle * Eff_L_Data[1] - passCross * passSingle * std::min(Eff_l_Data[1], Eff_L_Data[1])
          * Eff_tau_Data_taudown[idm] + passCross * Eff_l_Data[1] * Eff_tau_Data_taudown[idm];
        Eff_MC_tauup[idm] = passSingle * Eff_L_MC[1] - passCross * passSingle * std::min(Eff_l_MC[1], Eff_L_MC[1])
          * Eff_tau_MC_tauup[idm] + passCross * Eff_l_MC[1] * Eff_tau_MC_tauup[idm];
        Eff_MC_taudown[idm] = passSingle * Eff_L_MC[1] - passCross * passSingle * std::min(Eff_l_MC[1], Eff_L_MC[1])
          * Eff_tau_MC_taudown[idm] + passCross * Eff_l_MC[1] * Eff_tau_MC_taudown[idm];

        trigSF_tauup.push_back(Eff_Data_tauup[idm] / Eff_MC_tauup[idm]);
        trigSF_taudown.push_back(Eff_Data_taudown[idm] / Eff_MC_taudown[idm]);
      }

      double SFl = muTauTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SFtau = tauTrgSF_mutau.getSF(dau2_pt, dau2_decayMode, 0);
      trigSF_cross = SFl * SFtau;   
    } 
    
    else {
      double SF = muTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SF_error = muTrgSF.get_ScaleFactorError(dau1_pt, dau1_eta);
      trigSF_mu = {SF - SF_error, SF, SF + SF_error};
      trigSF_cross = SF;
      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        trigSF_tauup.push_back(SF);
        trigSF_taudown.push_back(SF);
      }
    }

    trigSF_single = muTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);

    return {trigSF_mu[1], trigSF_single, trigSF_cross, trigSF_mu[2], trigSF_mu[0], trigSF_mu[1],
      trigSF_mu[1], trigSF_tauup[0], trigSF_tauup[1], trigSF_tauup[2], trigSF_tauup[3],
      trigSF_taudown[0], trigSF_taudown[1], trigSF_taudown[2], trigSF_taudown[3],
      trigSF_mu[1], trigSF_mu[1]};

  }
  /////////////////////////////////////////////////////////////////////////////////////////
  // etau
  ///////////////////////////////////////////////////////////////////////////////////////// 
  else if (pairType == 1) {

    std::vector<double> trigSF_e, trigSF_tauup, trigSF_taudown;
    double trigSF_single, trigSF_cross;
    if (fabs(dau2_eta) < 2.1 && year_ != 2016) {
      int passSingle = (dau1_pt > etau_pt_th1_) ? 1 : 0;
      int passCross = (dau2_pt > etau_pt_th2_) ? 1 : 0;

      // lepton trigger
      auto Eff_L_Data_nom = eTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_L_MC_nom = eTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_L_Data_Err = eTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_L_MC_Err = eTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // std::cout << " --> eTrgSF : pt " << dau1_pt << " ; eta " << dau1_eta << " --> " << Eff_L_Data_nom << std::endl;

      std::vector <double> Eff_L_Data = {
        Eff_L_Data_nom - Eff_L_Data_Err, Eff_L_Data_nom, Eff_L_Data_nom + Eff_L_Data_Err};
      std::vector <double> Eff_L_MC = {
        Eff_L_MC_nom - Eff_L_MC_Err, Eff_L_MC_nom, Eff_L_MC_nom + Eff_L_MC_Err};

      // cross trigger
      // mu leg
      auto Eff_l_Data_nom = eTauTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_l_MC_nom = eTauTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_l_Data_Err = eTauTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_l_MC_Err = eTauTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // std::cout << " --> eTauTrgSF : pt " << dau1_pt << " ; eta " << dau1_eta << " --> " << Eff_l_Data_nom << std::endl;

      std::vector <double> Eff_l_Data = {
        Eff_l_Data_nom - Eff_l_Data_Err, Eff_l_Data_nom, Eff_l_Data_nom + Eff_l_Data_Err};
      std::vector <double> Eff_l_MC = {Eff_l_MC_nom - Eff_l_MC_Err, Eff_l_MC_nom, Eff_l_MC_nom + Eff_l_MC_Err};

      // tau leg
      auto Eff_tau_Data = tauTrgSF_etau.getEfficiencyData(dau2_pt, dau2_decayMode, 0);
      auto Eff_tau_MC = tauTrgSF_etau.getEfficiencyMC(dau2_pt, dau2_decayMode, 0);

      std::vector <double> Eff_Data_e, Eff_MC_e;
      for (size_t i = 0; i< Eff_l_Data.size(); i++) {
        Eff_Data_e.push_back(passSingle * Eff_L_Data[i] - passCross * passSingle * std::min(Eff_l_Data[i], Eff_L_Data[i])
          * Eff_tau_Data + passCross * Eff_l_Data[i] * Eff_tau_Data);
        Eff_MC_e.push_back(passSingle * Eff_L_MC[i] - passCross * passSingle * std::min(Eff_l_MC[i], Eff_L_MC[i])
          * Eff_tau_MC + passCross * Eff_l_MC[i] * Eff_tau_MC);
        trigSF_e.push_back(Eff_Data_e[i] / Eff_MC_e[i]);
      }

      std::vector<double> Eff_tau_Data_tauup(4, Eff_tau_Data);
      std::vector<double> Eff_tau_Data_taudown(4, Eff_tau_Data);
      std::vector<double> Eff_tau_MC_tauup(4, Eff_tau_MC);
      std::vector<double> Eff_tau_MC_taudown(4, Eff_tau_MC);
      std::vector<double> Eff_Data_tauup(4);
      std::vector<double> Eff_Data_taudown(4);
      std::vector<double> Eff_MC_tauup(4);
      std::vector<double> Eff_MC_taudown(4);

      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        if (decayModes[idm] == dau2_decayMode) {
          Eff_tau_Data_tauup[idm] = tauTrgSF_etau.getEfficiencyData(dau2_pt, dau2_decayMode, 1);
          Eff_tau_Data_taudown[idm] = tauTrgSF_etau.getEfficiencyData(dau2_pt, dau2_decayMode, -1);
          Eff_tau_MC_tauup[idm] = tauTrgSF_etau.getEfficiencyMC(dau2_pt, dau2_decayMode, 1);
          Eff_tau_MC_taudown[idm] = tauTrgSF_etau.getEfficiencyMC(dau2_pt, dau2_decayMode, -1);
        }
      }

      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        Eff_Data_tauup[idm] = passSingle * Eff_L_Data[1] - passCross * passSingle * std::min(Eff_l_Data[1], Eff_L_Data[1])
          * Eff_tau_Data_tauup[idm] + passCross * Eff_l_Data[1] * Eff_tau_Data_tauup[idm];
        Eff_Data_taudown[idm] = passSingle * Eff_L_Data[1] - passCross * passSingle * std::min(Eff_l_Data[1], Eff_L_Data[1])
          * Eff_tau_Data_taudown[idm] + passCross * Eff_l_Data[1] * Eff_tau_Data_taudown[idm];
        Eff_MC_tauup[idm] = passSingle * Eff_L_MC[1] - passCross * passSingle * std::min(Eff_l_MC[1], Eff_L_MC[1])
          * Eff_tau_MC_tauup[idm] + passCross * Eff_l_MC[1] * Eff_tau_MC_tauup[idm];
        Eff_MC_taudown[idm] = passSingle * Eff_L_MC[1] - passCross * passSingle * std::min(Eff_l_MC[1], Eff_L_MC[1])
          * Eff_tau_MC_taudown[idm] + passCross * Eff_l_MC[1] * Eff_tau_MC_taudown[idm];

        trigSF_tauup.push_back(Eff_Data_tauup[idm] / Eff_MC_tauup[idm]);
        trigSF_taudown.push_back(Eff_Data_taudown[idm] / Eff_MC_taudown[idm]);
      }
      
      double SFl = eTauTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SFtau = tauTrgSF_etau.getSF(dau2_pt, dau2_decayMode, 0);
      trigSF_cross = SFl * SFtau;
    }

    else {
      double SF = eTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SF_error = eTrgSF.get_ScaleFactorError(dau1_pt, dau1_eta);
      trigSF_e = {SF - SF_error, SF, SF + SF_error};
      trigSF_cross = SF;
      for (size_t idm = 0; idm < decayModes.size(); idm++) {
        trigSF_tauup.push_back(SF);
        trigSF_taudown.push_back(SF);
      }
    }

    trigSF_single = eTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);

    return {trigSF_e[1], trigSF_single, trigSF_cross, trigSF_e[1], trigSF_e[1], trigSF_e[2],
      trigSF_e[0], trigSF_tauup[0], trigSF_tauup[1], trigSF_tauup[2], trigSF_tauup[3],
      trigSF_taudown[0], trigSF_taudown[1], trigSF_taudown[2], trigSF_taudown[3],
      trigSF_e[1], trigSF_e[1]};

  } 
  /////////////////////////////////////////////////////////////////////////////////////////
  // tautau
  /////////////////////////////////////////////////////////////////////////////////////////  
  else if (pairType == 2 && isVBFtrigger == 0) {

    auto SF1 = tauTrgSF_ditau.getSF(dau1_pt, dau1_decayMode, 0);
    auto SF2 = tauTrgSF_ditau.getSF(dau2_pt, dau2_decayMode, 0);

    std::vector<double> SF1_tauup(4, SF1);
    std::vector<double> SF1_taudown(4, SF1);
    std::vector<double> SF2_tauup(4, SF2);
    std::vector<double> SF2_taudown(4, SF2);

    for (size_t idm = 0; idm < decayModes.size(); idm++) {
      if (decayModes[idm] == dau1_decayMode) {
        SF1_tauup[idm] = tauTrgSF_ditau.getSF(dau1_pt, dau1_decayMode, 1);
        SF1_taudown[idm] = tauTrgSF_ditau.getSF(dau1_pt, dau1_decayMode, -1);
      }
      if (decayModes[idm] == dau2_decayMode) {
        SF2_tauup[idm] = tauTrgSF_ditau.getSF(dau2_pt, dau2_decayMode, 1);
        SF2_taudown[idm] = tauTrgSF_ditau.getSF(dau2_pt, dau2_decayMode, -1);
      }
    }
    double trigSF = SF1 * SF2;
    std::vector<double> trigSF_tauup, trigSF_taudown;
    for (size_t idm = 0; idm < decayModes.size(); idm++) {
      trigSF_tauup.push_back(SF1_tauup[idm] * SF2_tauup[idm]);
      trigSF_taudown.push_back(SF1_taudown[idm] * SF2_taudown[idm]);
    }
    return {trigSF, trigSF, trigSF, trigSF, trigSF, trigSF, trigSF,
      trigSF_tauup[0], trigSF_tauup[1], trigSF_tauup[2], trigSF_tauup[3],
      trigSF_taudown[0], trigSF_taudown[1], trigSF_taudown[2], trigSF_taudown[3],
      trigSF, trigSF};

  } 
  else if (pairType == 2 && year_ != 2016) {
    if (vbfjet1_pt >= 0 && vbfjet2_pt >= 0) {
      auto vbfjet1_tlv = TLorentzVector();
      auto vbfjet2_tlv = TLorentzVector();
      vbfjet1_tlv.SetPtEtaPhiM(vbfjet1_pt, vbfjet1_eta, vbfjet1_phi, vbfjet1_mass);
      vbfjet2_tlv.SetPtEtaPhiM(vbfjet2_pt, vbfjet2_eta, vbfjet2_phi, vbfjet2_mass);
      auto vbfjj_mass = (vbfjet1_tlv + vbfjet2_tlv).M();
      if (vbfjet1_tlv.Pt() > 140 && vbfjet2_tlv.Pt() > 60 && vbfjj_mass > 800
          && dau1_pt > 25 && dau2_pt > 25
          && (dau1_pt <= 40 || dau2_pt <= 40)) {
        auto jetSF = getContentHisto3D(
          jetTrgSF_vbf, vbfjj_mass, vbfjet1_tlv.Pt(), vbfjet2_tlv.Pt(), 0);
        auto jetSFerror = getContentHisto3D(
          jetTrgSF_vbf, vbfjj_mass, vbfjet1_tlv.Pt(), vbfjet2_tlv.Pt(), 0);
        auto SF1 = tauTrgSF_vbf.getSF(dau1_pt, dau1_decayMode, 0);
        auto SF2 = tauTrgSF_vbf.getSF(dau2_pt, dau2_decayMode, 0);

        std::vector<double> SF1_tauup(4, SF1);
        std::vector<double> SF1_taudown(4, SF1);
        std::vector<double> SF2_tauup(4, SF2);
        std::vector<double> SF2_taudown(4, SF2);
        for (size_t idm = 0; idm < decayModes.size(); idm++) {
          if (decayModes[idm] == dau1_decayMode) {
            SF1_tauup[idm] = tauTrgSF_ditau.getSF(dau1_pt, dau1_decayMode, 1);
            SF1_taudown[idm] = tauTrgSF_ditau.getSF(dau1_pt, dau1_decayMode, -1);
          }
          if (decayModes[idm] == dau2_decayMode) {
            SF2_tauup[idm] = tauTrgSF_ditau.getSF(dau2_pt, dau2_decayMode, 1);
            SF2_taudown[idm] = tauTrgSF_ditau.getSF(dau2_pt, dau2_decayMode, -1);
          }
        }

        std::vector<double> trigSF_vbfjet = {
          (jetSF - jetSFerror) * SF1 * SF2, jetSF * SF1 * SF2, (jetSF + jetSFerror) * SF1 * SF2
        };
        std::vector<double> trigSF_tauup, trigSF_taudown;
        for (size_t idm = 0; idm < decayModes.size(); idm++) {
          trigSF_tauup.push_back(jetSF * SF1_tauup[idm] * SF2_tauup[idm]);
          trigSF_taudown.push_back(jetSF * SF1_taudown[idm] * SF2_taudown[idm]);
        }

        return {trigSF_vbfjet[1], trigSF_vbfjet[1], trigSF_vbfjet[1], trigSF_vbfjet[1],
          trigSF_vbfjet[1], trigSF_vbfjet[1], trigSF_vbfjet[1],
          trigSF_tauup[0], trigSF_tauup[1], trigSF_tauup[2], trigSF_tauup[3],
          trigSF_taudown[0], trigSF_taudown[1], trigSF_taudown[2], trigSF_taudown[3],
          trigSF_vbfjet[2], trigSF_vbfjet[0]};
      }
    }
  }
  return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

// Destructor
Htt_trigSFinterface::~Htt_trigSFinterface() {}

