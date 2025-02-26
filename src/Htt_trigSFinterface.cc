#include "Tools/Tools/interface/Htt_trigSFinterface.h"

double GetScaleFactorError(double effData, double effMC, double errData, double errMC) {
  // Compute absolute error on SF=effData/effMC propagating uncertainties. Taken from KLUB
  double SF_error = 0.;
  if (effData==0. || effMC==0.) {
    std::cout<<"WARNING in GetScaleFactorError(double effData, double effMC, double errData, double errMC): efficiency in data or MC = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0. effData=" << effData << " effMC="<<effMC << std::endl;
    return 0.;
  }
  else {
    SF_error = pow((errData/effData),2) + pow((errMC/effMC),2);
    SF_error = pow(SF_error, 0.5)*(effData/effMC);
  }
  return SF_error;
}

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

Htt_trigSFresult Htt_trigSFinterface::get_scale_factors(int pairType, int isVBFtrigger,
    int dau1_decayMode, float dau1_pt, float dau1_eta,
    int dau2_decayMode, float dau2_pt, float dau2_eta) {

  std::vector <int> decayModes = {0, 1, 10, 11};

  Htt_trigSFresult res{};

  /////////////////////////////////////////////////////////////////////////////////////////
  // mutau
  /////////////////////////////////////////////////////////////////////////////////////////
  
  if (pairType == 0) {
    if (fabs(dau2_eta) < 2.1) { // region where cross tau trigger is defined
      // Use >= instead of > as it is possible that the pt is exactly equal to the threshold
      int passSingle = (dau1_pt >= mutau_pt_th1_) ? 1 : 0;
      int passCross = (dau2_pt >= mutau_pt_th2_) ? 1 : 0;
      assert(passSingle || passCross);
      res.passSingle = passSingle;
      res.passCross = passCross;

      // single-lepton trigger
      auto Eff_SL_Data = muTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_SL_MC = muTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_SL_Data_Err = muTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_SL_MC_Err = muTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // -- cross trigger
      // mu leg
      auto Eff_cross_mu_Data = muTauTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_cross_mu_MC = muTauTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_cross_mu_Data_Err = muTauTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_cross_mu_MC_Err = muTauTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // tau leg
      auto Eff_cross_tau_Data = tauTrgSF_mutau.getEfficiencyData(dau2_pt, dau2_decayMode, 0);
      auto Eff_cross_tau_MC = tauTrgSF_mutau.getEfficiencyMC(dau2_pt, dau2_decayMode, 0);

      // -- SF Nominal
      double Eff_Data = passSingle * Eff_SL_Data - passCross * passSingle * std::min(Eff_cross_mu_Data, Eff_SL_Data) * Eff_cross_tau_Data + passCross * Eff_cross_mu_Data * Eff_cross_tau_Data;
			double Eff_MC   = passSingle * Eff_SL_MC   - passCross * passSingle * std::min(Eff_cross_mu_MC  , Eff_SL_MC)   * Eff_cross_tau_MC   + passCross * Eff_cross_mu_MC   * Eff_cross_tau_MC;
      res.trigSF = Eff_Data / Eff_MC;

      // muon leg error
      /* Possible cases : 
      passSingle & passCross & Eff_cross_mu_Data <= Eff_SL_Data (most likely) : Err_Data_SL + Err_Data_cross_mu = Eff_SL_Data_Err 
      !passSingle & passCross : Err_Data_SL + Err_Data_cross_mu = Eff_cross_mu_Data_Err * Eff_cross_tau_Data
      others: ....
      */
      double Err_Data_SL = passSingle * Eff_SL_Data_Err - passCross * passSingle * (Eff_cross_mu_Data > Eff_SL_Data) * Eff_SL_Data_Err * Eff_cross_tau_Data;
			double Err_MC_SL   = passSingle * Eff_SL_MC_Err   - passCross * passSingle * (Eff_cross_mu_MC   > Eff_SL_MC)   * Eff_SL_MC_Err   * Eff_cross_tau_MC;

      double Err_Data_cross_mu = - passCross * passSingle * (Eff_cross_mu_Data <= Eff_SL_Data) * Eff_cross_mu_Data_Err * Eff_cross_tau_Data + passCross * Eff_cross_mu_Data_Err * Eff_cross_tau_Data;
			double Err_MC_cross_mu   = - passCross * passSingle * (Eff_cross_mu_MC   <= Eff_SL_MC)   * Eff_cross_mu_MC_Err   * Eff_cross_tau_MC   + passCross * Eff_cross_mu_MC_Err   * Eff_cross_tau_MC;

      double muErr = GetScaleFactorError(Eff_Data, Eff_MC, Err_Data_SL + Err_Data_cross_mu, Err_MC_SL + Err_MC_cross_mu);
      res.muUp = res.trigSF + muErr;
      res.muDown = res.trigSF - muErr;

      // tau leg errors
      // if single lepton trigger and cross trigger are passed and the lepton leg of the cross trigger is
      // less efficient than the single lepton trigger, the "combined" efficiency reduces to the single
      // lepton trigger efficiency and the uncertainties on the cross trigger legs are 0, if this is the
      // case in both Data and MC also the uncertainties on the trigger scale factors are 0. (If statement
      // to avoid unnecessary get_ScaleFactorError calls and warnings)
      if(passCross and !(passSingle and (Eff_cross_mu_Data <= Eff_SL_Data and Eff_cross_mu_MC <= Eff_SL_MC))) {
        for (size_t idm = 0; idm < decayModes.size(); idm++) {
          if (decayModes[idm] == dau2_decayMode) {
            double Eff_cross_tau_Data_Up = tauTrgSF_mutau.getEfficiencyData(dau2_pt, dau2_decayMode, 1);
            double Eff_cross_tau_MC_Up   = tauTrgSF_mutau.getEfficiencyMC(dau2_pt, dau2_decayMode, 1);
            double Eff_Data_Up = passSingle * Eff_SL_Data - passCross * passSingle * std::min(Eff_cross_mu_Data, Eff_SL_Data) * Eff_cross_tau_Data_Up   + passCross * Eff_cross_mu_Data * Eff_cross_tau_Data_Up;
            double Eff_MC_Up   = passSingle * Eff_SL_MC   - passCross * passSingle * std::min(Eff_cross_mu_MC  , Eff_SL_MC)   * Eff_cross_tau_MC_Up     + passCross * Eff_cross_mu_MC   * Eff_cross_tau_MC_Up;
            res.setDMError(decayModes[idm], GetScaleFactorError(Eff_Data, Eff_MC, Eff_Data_Up - Eff_Data, Eff_MC_Up - Eff_MC));
          }
        }
      }
    } 
    else { // single-lepton trigger only region 
      double SF = muTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SF_error = muTrgSF.get_ScaleFactorError(dau1_pt, dau1_eta);
      res.trigSF = SF;
      res.muUp = SF + SF_error;
      res.muDown = SF - SF_error;
    }

    res.setSystVariationsToNominal();
    return res;
  }
  /////////////////////////////////////////////////////////////////////////////////////////
  // etau
  ///////////////////////////////////////////////////////////////////////////////////////// 
  else if (pairType == 1) {
    if (fabs(dau2_eta) < 2.1 && year_ != 2016) {
      // Use >= instead of > as it is possible that the pt is exactly equal to the threshold
      int passSingle = (dau1_pt >= etau_pt_th1_) ? 1 : 0;
      int passCross = (dau2_pt >= etau_pt_th2_) ? 1 : 0;
      assert(passSingle || passCross);
      res.passSingle = passSingle;
      res.passCross = passCross;

      // lepton trigger
      auto Eff_SL_Data = eTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_SL_MC = eTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_SL_Data_Err = eTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_SL_MC_Err = eTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // cross trigger
      // e leg
      auto Eff_cross_ele_Data = eTauTrgSF.get_EfficiencyData(dau1_pt, dau1_eta);
      auto Eff_cross_ele_MC = eTauTrgSF.get_EfficiencyMC(dau1_pt, dau1_eta);
      auto Eff_cross_ele_Data_Err = eTauTrgSF.get_EfficiencyDataError(dau1_pt, dau1_eta);
      auto Eff_cross_ele_MC_Err = eTauTrgSF.get_EfficiencyMCError(dau1_pt, dau1_eta);

      // tau leg
      auto Eff_cross_tau_Data = tauTrgSF_etau.getEfficiencyData(dau2_pt, dau2_decayMode, 0);
      auto Eff_cross_tau_MC = tauTrgSF_etau.getEfficiencyMC(dau2_pt, dau2_decayMode, 0);


      // -- SF Nominal
      double Eff_Data = passSingle * Eff_SL_Data - passCross * passSingle * std::min(Eff_cross_ele_Data, Eff_SL_Data) * Eff_cross_tau_Data + passCross * Eff_cross_ele_Data * Eff_cross_tau_Data;
			double Eff_MC   = passSingle * Eff_SL_MC   - passCross * passSingle * std::min(Eff_cross_ele_MC  , Eff_SL_MC)   * Eff_cross_tau_MC   + passCross * Eff_cross_ele_MC   * Eff_cross_tau_MC;
      res.trigSF = Eff_Data / Eff_MC;

      // electron leg error
      /* Possible cases : 
      passSingle & passCross & Eff_cross_mu_Data <= Eff_SL_Data (most likely) : Err_Data_SL + Err_Data_cross_mu = Eff_SL_Data_Err 
      !passSingle & passCross : Err_Data_SL + Err_Data_cross_mu = Eff_cross_mu_Data_Err * Eff_cross_tau_Data
      others: ....
      */
      double Err_Data_SL = passSingle * Eff_SL_Data_Err - passCross * passSingle * (Eff_cross_ele_Data > Eff_SL_Data) * Eff_SL_Data_Err * Eff_cross_tau_Data;
			double Err_MC_SL   = passSingle * Eff_SL_MC_Err   - passCross * passSingle * (Eff_cross_ele_MC   > Eff_SL_MC)   * Eff_SL_MC_Err   * Eff_cross_tau_MC;

      double Err_Data_cross_ele = - passCross * passSingle * (Eff_cross_ele_Data <= Eff_SL_Data) * Eff_cross_ele_Data_Err * Eff_cross_tau_Data + passCross * Eff_cross_ele_Data_Err * Eff_cross_tau_Data;
			double Err_MC_cross_ele   = - passCross * passSingle * (Eff_cross_ele_MC   <= Eff_SL_MC)   * Eff_cross_ele_MC_Err   * Eff_cross_tau_MC   + passCross * Eff_cross_ele_MC_Err   * Eff_cross_tau_MC;

      double eleErr = GetScaleFactorError(Eff_Data, Eff_MC, Err_Data_SL + Err_Data_cross_ele, Err_MC_SL + Err_MC_cross_ele); 
      res.eleUp = res.trigSF + eleErr;
      res.eleDown = res.trigSF - eleErr;

      // tau leg errors
      // if single lepton trigger and cross trigger are passed and the lepton leg of the cross trigger is
      // less efficient than the single lepton trigger, the "combined" efficiency reduces to the single
      // lepton trigger efficiency and the uncertainties on the cross trigger legs are 0, if this is the
      // case in both Data and MC also the uncertainties on the trigger scale factors are 0. (If statement
      // to avoid unnecessary get_ScaleFactorError calls and warnings)
      if(passCross and !(passSingle and (Eff_cross_ele_Data <= Eff_SL_Data and Eff_cross_ele_MC <= Eff_SL_MC))) {
        for (size_t idm = 0; idm < decayModes.size(); idm++) {
          if (decayModes[idm] == dau2_decayMode) {
            double Eff_cross_tau_Data_Up = tauTrgSF_etau.getEfficiencyData(dau2_pt, dau2_decayMode, 1);
            double Eff_cross_tau_MC_Up   = tauTrgSF_etau.getEfficiencyMC(dau2_pt, dau2_decayMode, 1);
            // !passSingle & passCross : Eff_Data_Up =  Eff_cross_ele_Data * Eff_cross_tau_Data_Up
            double Eff_Data_Up = passSingle * Eff_SL_Data - passCross * passSingle * std::min(Eff_cross_ele_Data, Eff_SL_Data) * Eff_cross_tau_Data_Up   + passCross * Eff_cross_ele_Data * Eff_cross_tau_Data_Up;
            double Eff_MC_Up   = passSingle * Eff_SL_MC   - passCross * passSingle * std::min(Eff_cross_ele_MC  , Eff_SL_MC)   * Eff_cross_tau_MC_Up     + passCross * Eff_cross_ele_MC   * Eff_cross_tau_MC_Up;
            // !passSingle & passCross : error = GetScaleFactorError(Eff_data, Eff_MC, Eff_cross_ele_Data * (Eff_cross_tau_Data_Up-Eff_cross_tau_Data), ...)
            //std::cout << Eff_Data << " " << Eff_MC << " " << " Eff_Data_Up=" << Eff_Data_Up << " Eff_MC_Up=" << Eff_MC_Up << std::endl;
            res.setDMError(decayModes[idm], GetScaleFactorError(Eff_Data, Eff_MC, Eff_Data_Up - Eff_Data, Eff_MC_Up - Eff_MC));
          }
        }
      }
    }
    else {
      double SF = eTrgSF.get_ScaleFactor(dau1_pt, dau1_eta);
      double SF_error = eTrgSF.get_ScaleFactorError(dau1_pt, dau1_eta);
      res.trigSF = SF;
      res.eleUp = SF + SF_error;
      res.eleDown = SF - SF_error;
    }

    res.setSystVariationsToNominal();
    return res;
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
    res.trigSF = SF1 * SF2;
    res.DM0Up = SF1_tauup[0] * SF2_tauup[0];
    res.DM0Down = SF1_taudown[0] * SF2_taudown[0];
    res.DM1Up = SF1_tauup[1] * SF2_tauup[1];
    res.DM1Down = SF1_taudown[1] * SF2_taudown[1];
    res.DM10Up = SF1_tauup[2] * SF2_tauup[2];
    res.DM10Down = SF1_taudown[2] * SF2_taudown[2];
    res.DM11Up = SF1_tauup[3] * SF2_tauup[3];
    res.DM11Down = SF1_taudown[3] * SF2_taudown[3];
    res.setSystVariationsToNominal();
    return res;
  } 
  else {
    res.setSystVariationsToNominal();
    return res;
  }
}

MET_trigSF_interface::MET_trigSF_interface(std::string SF_file) {
  std::unique_ptr<TFile> file( TFile::Open(SF_file.c_str()) );
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not open MET trigger SF file " + SF_file);
  }
  trigSF_h.reset(file->Get<TH1F>("SF/MET_SF"));
  if (!trigSF_h)
    throw std::runtime_error("Failure reading MET trigger SF file " + SF_file);
  trigSF_h->SetDirectory(nullptr);
}


MET_trigSF_interface::trigSF_result MET_trigSF_interface::getSF(float MET_pt) {
  if (MET_pt < 120) {
    throw std::runtime_error("Attempting to get MET trigger SF with MET_pt<80GeV");
  }
  // if (MET_pt < 120) {
  //   std::cerr << "WARNING : attempting to get MET trigger SF outside of trigger plateau" << std::endl;
  //   MET_pt = 120;
  // }
  if (MET_pt >= 500)
    MET_pt = 500-0.1; // upper edge is excluded in histograms

  int binNumber = trigSF_h->GetXaxis()->FindFixBin(MET_pt);
  if (binNumber <= 0 || binNumber > trigSF_h->GetNbinsX())
    throw std::runtime_error(std::string("MET trigger SF read error. Tried to read MET_pt ") + std::to_string(MET_pt) + " which is outside bounds [120, 500]");

  MET_trigSF_interface::trigSF_result sf_res;
  sf_res.SF = trigSF_h->GetBinContent(binNumber);
  sf_res.SF_statup = sf_res.SF + trigSF_h->GetBinErrorUp(binNumber);
  sf_res.SF_statdown = sf_res.SF - trigSF_h->GetBinErrorLow(binNumber);
  return sf_res;
}
