#ifndef HHJetsInterface_h
#define HHJetsInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHJetsInterface                                                                                        //
//                                                                                                                //
//   Class to compute HHbtag output.                                                                              //
//                                                                                                                //
//   Author: Jaime León Holgado                                                                                   //
//   Date  : June 2021                                                                                            //
//                                                                                                                //
//   modified from                                                                                                //
//                                                                                                                //
//   https://github.com/LLRCMS/KLUBAnalysis/blob/VBF_legacy/interface/HHbtagKLUBinterface.h                       //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

// HHbtag libraries
#include "HHTools/HHbtag/interface/HH_BTag.h"

// ROOT libraries
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>
#include <Math/VectorUtil.h>

// CMSSW
#include "DataFormats/Math/interface/deltaPhi.h"

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

struct jet_idx_btag {
  int idx;
  float btag;
};

struct jet_pair_mass {
  int idx1;
  int idx2;
  float inv_mass;
};

enum class JetCategory {
  Res_2b = 0,
  Res_1b = 1,
  Boosted_bb = 2,
  None = -1,
  Boosted_failedPNet = -2 // events which have a FatJet passing selections (pt, softdrop, etc) but failing ParticleNet cut. They are vetoed as we don't have SFs for jets failing PNet
};

enum class JetCategoryPriorityMode {
  Res2b_Boosted_Res1b_noPNetFail, // Priority to resolved2b, then boosted, then res1b but don't consider events with a fatjet passing all selections except PNet for res1b (because we don't have SFs for jets failing PNet). Used for HPSTaus
  Boosted_Res2b_Res1b_noPNetFail // Priority to boosted. Do not consider events with fatjet failing PNet
};

struct output {
  std::vector <float> hhbtag;
  int bjet_idx1;
  int bjet_idx2;
  int vbfjet_idx1;
  int vbfjet_idx2;
  std::vector<int> ctjet_indexes;
  std::vector<int> fwjet_indexes;
  JetCategory jetCategory;
  int fatjet_idx; // will be filled even when the FatJet fails PNet (so when category == 2 or -2)
};


bool jetSort (const jet_idx_btag& jA, const jet_idx_btag& jB)
{
  return (jA.btag > jB.btag);
}

bool jetPairSort (const jet_pair_mass& jA, const jet_pair_mass& jB)
{
  return (jA.inv_mass > jB.inv_mass);
}


// HHJetsInterface class
class HHJetsInterface {

  public:
    HHJetsInterface (std::string model_0, std::string model_1, int year, bool isUL, float btag_wp, float fatjet_bbtag_wp);
    ~HHJetsInterface ();
    
  output GetHHJets(unsigned long long int event, int pairType,
    fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass,
    iRVec Jet_puId, fRVec Jet_jetId, fRVec Jet_btagDeepFlavB,
    fRVec FatJet_pt, fRVec FatJet_eta, fRVec FatJet_phi, fRVec FatJet_mass,
    fRVec FatJet_msoftdrop, fRVec FatJet_jetId, fRVec FatJet_particleNet_XbbVsQCD,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi, JetCategoryPriorityMode priorityMode);

  std::vector<float> GetScore(
    std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
    std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_deepFlavour_,
    std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
    float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_);

  private:
    hh_btag::HH_BTag HHbtagger_;
    int year_;
    float max_bjet_eta = 2.4;
    float btag_wp_; // Working point for b-tagging score, used when giving priority to resolved category
    float fatjet_bbtag_wp_; // Working point for FatJet bb-tagging score, used when giving priority to boosted
};

#endif // HHJetsInterface
