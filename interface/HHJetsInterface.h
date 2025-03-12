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
#include <ostream>

// HHbtag libraries
#include "HHTools/HHbtag/interface/HH_BTag.h"

// ROOT libraries
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>
#include <Math/VectorUtil.h>

// CMSSW
#include "DataFormats/Math/interface/deltaPhi.h"

// Using const& will lead to an error when using the wrong type, instead of silently converting, which can help catch mismatched branches
typedef ROOT::VecOps::RVec<float> const& rfRVec;
typedef ROOT::VecOps::RVec<UChar_t> const& rcRVec;

struct jet_idx_btag {
  int idx;
  float btag;
};

struct jet_pair_mass {
  int idx1;
  int idx2;
  float inv_mass;
};

/** Record for a jet of which cuts failed */
struct JetsFailReason {
  JetsFailReason() : Reco(false), Pt(false), Eta(false), JetID(false), JetPUID(false), SoftDrop(false), DeltaRDau(false)
  {}
  bool pass() const { return !(Reco || Pt || Eta || JetID || JetPUID || SoftDrop || DeltaRDau ); }
  int countFails() const { return Reco + Pt + Eta + JetID + JetPUID + SoftDrop + DeltaRDau; }
  void print(std::ostream& out) const {
    out << "JetsFailReason : ";
    if (Reco) out << "Reco ";
    if (Pt) out << "Pt ";
    if (Eta) out << "Eta ";
    if (JetID) out << "JetID ";
    if (JetPUID) out << "JetPUID ";
    if (SoftDrop) out << "SoftDrop ";
    if (DeltaRDau) out << "DeltaRDau ";
    out << std::endl;
  }

  bool Reco; // no gen-matched jet in collection

  bool Pt; bool Eta;
  bool JetID;
  bool JetPUID;
  bool SoftDrop;

  bool DeltaRDau; // Too close to a tau
};

struct Jets_cutflow_output {
    JetsFailReason jet1_failReason; // AK4 jets fail reason
    JetsFailReason jet2_failReason;
    JetsFailReason fatjet_failReason;

    bool wrongJet; // wrong AK4 jet (not genmatched)
    bool wrongFatJet; // wrong AK8 jet
};

struct output {
  std::vector <float> hhbtag;
  int bjet_idx1;
  int bjet_idx2;
  int vbfjet_idx1;
  int vbfjet_idx2;
  std::vector<int> ctjet_indexes;
  std::vector<int> fwjet_indexes;
  int fatjet_idx; // will be filled even when the FatJet fails PNet (so when category == 2 or -2)
  Jets_cutflow_output cutflow_output;
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
    
  output GetHHJets(unsigned long long int event, int pairType,
    rfRVec Jet_pt, rfRVec Jet_eta, rfRVec Jet_phi, rfRVec Jet_mass,
    rcRVec Jet_puId, rcRVec Jet_jetId, rfRVec Jet_btagDeepFlavB,
    rfRVec FatJet_pt, rfRVec FatJet_eta, rfRVec FatJet_phi, rfRVec FatJet_mass,
    rfRVec FatJet_msoftdrop, rcRVec FatJet_jetId, rfRVec FatJet_particleNet_XbbVsQCD,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi, 
    rfRVec GenPart_pt, rfRVec GenPart_eta, rfRVec GenPart_phi, rfRVec GenPart_mass,
    int genXbb_GenPartIdx, int gen_b1_GenPartIdx, int gen_b2_GenPartIdx,
    bool doGenCutFlow);

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


class HHJetsCategoryInterface {
public:
    HHJetsCategoryInterface (float btag_wp, float fatjet_bbtag_wp);
    JetCategory GetJetCategory(JetCategoryPriorityMode priorityMode, int bjet1_idx, int bjet2_idx, int fatjet_idx, float bjet1_btagDeepFlavB, float bjet2_btagDeepFlavB, float fatjet_particleNet_XbbVsQCD);

private:
    float btag_wp_; // Working point for b-tagging score, used when giving priority to resolved category
    float fatjet_bbtag_wp_; // Working point for FatJet bb-tagging score, used when giving priority to boosted
};

#endif // HHJetsInterface
