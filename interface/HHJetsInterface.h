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
  size_t idx;
  float btag;
};

struct jet_pair_mass {
  int idx1;
  int idx2;
  float inv_mass;
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
    HHJetsInterface (std::string model_0, std::string model_1, int year);
    ~HHJetsInterface ();
    
  std::vector<float> GetHHJets(unsigned long long int event, int pairType,
    fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass,
    iRVec Jet_puId, fRVec Jet_jetId, fRVec Jet_btagDeepFlavB,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi);

    std::vector<float> GetScore(
      std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
      std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_deepFlavour_,
      std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
      float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_);
  
    void set_bjet_indexes(int idx1, int idx2) {
      bjet_indexes[0] = idx1;
      bjet_indexes[1] = idx2;
    }

    std::vector <int> get_bjet_indexes() {
      return bjet_indexes;
    }

    void set_vbfjet_indexes(int idx1, int idx2) {
      vbfjet_indexes[0] = idx1;
      vbfjet_indexes[1] = idx2;
    }

    std::vector <int> get_vbfjet_indexes() {
      return vbfjet_indexes;
    }

    void setBoosted(int isBoosted) {
      isBoosted_ = isBoosted;
    }

    int isBoosted() {
      return isBoosted_;
    }

  private:
    hh_btag::HH_BTag HHbtagger_;
    int year_;
    std::vector <int> bjet_indexes;
    std::vector <int> vbfjet_indexes;
    int isBoosted_;
};

#endif // HHJetsInterface