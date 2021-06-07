#include "Tools/Tools/interface/HHbtagInterface.h"

// Constructor
HHbtagInterface::HHbtagInterface (std::string model_0, std::string model_1, int year):
  HHbtagger_(std::array<std::string, 2> { {model_0, model_1} })
{
  // Set the year just once
  // std::array<std::string, 2> models_;
  // std::ostringstream ss_model_0, ss_model_1;
  // ss_model_0 << model_0;
  // ss_model_1 << model_1;
  // models_[0] = ss_model_0.str();
  // models_[1] = ss_model_1.str();
  // HHbtagger_(models_);
}


// Destructor
HHbtagInterface::~HHbtagInterface() {}

// GetScore
std::vector<float> HHbtagInterface::GetScore(
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
