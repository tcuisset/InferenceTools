#ifndef HHLeptonInterface_h
#define HHLeptonInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHLeptonInterface                                                                                      //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <string>
#include <cmath>

// ROOT libraries
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>

// CMSSW
#include "DataFormats/Math/interface/deltaPhi.h"

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;
typedef ROOT::VecOps::RVec<int16_t> sRVec;
typedef ROOT::VecOps::RVec<UChar_t> uRVec;

struct tau_pair {
  int index1;
  float iso1;
  float pt1;
  int index2;
  float iso2;
  float pt2;
  int isTauTauJetTrigger;
  int isVBFtrigger;
};

struct trig_req {
  bool pass; // enable trigger
  float pt1_offline; // min offline pt of first leg (ele/mu/tauh)
  float eta1_offline; // max offline eta of first leg
  float pt2_offline; // min offline pt of second leg (hadronic tau)
  float eta2_offline; // max offline eta of second leg
  float pt1_online; // min TrigObj_pt of first leg
  float pt2_online; // min TrigObj_pt of second leg
  std::vector<std::vector<int>> bits; // List of trigger bits to match, for each trigger leg (put expanded powers of 2, ie 1, 2, 4, 8, etc).
};

struct lepton_output {
    int pairType;
    int dau1_index;
    int dau2_index;
    int isTauTauJetTrigger;
    int isVBFtrigger;
    int isOS;

    float dau1_eta;
    float dau1_phi;
    float dau1_iso;
    int dau1_decayMode;
    int dau1_idDeepTauVSe;
    int dau1_idDeepTauVSmu;
    int dau1_idDeepTauVSjet;

    float dau2_eta;
    float dau2_phi;
    int dau2_decayMode;
    int dau2_idDeepTauVSe;
    int dau2_idDeepTauVSmu;
    int dau2_idDeepTauVSjet;
};

struct FailReason {
  FailReason() : Reco(false), Pt(false), Eta(false), Vertex(false), LeptonID(false), LeptonIso(false), TauIdVsMu(false), TauIdVsE(false), TauIdVsJet(false), TauDM(false), WrongPair(false)
  {}
  bool pass() { return !(Reco || Pt || Eta || Vertex || LeptonID || LeptonIso || TauIdVsMu || TauIdVsE || TauIdVsJet || TauDM || WrongPair); }

  bool Reco;

  bool Pt; bool Eta;
  bool Vertex;
  bool LeptonID; bool LeptonIso;
  bool TauIdVsMu;
  bool TauIdVsE;
  bool TauIdVsJet;
  bool TauDM;

  bool WrongPair;
};

struct cutflow_output {
    // cutflow_output() : dau1_fail(), dau2_fail(), leptonVetoFail(false), deltaR(false)
    // {}
    FailReason dau1_fail;
    FailReason dau2_fail;
    bool leptonVetoFail;
    bool deltaR;
    bool wrongChannel;
};

// return true if pA > pB using the sorting criteria. For tautau (symmetrical)
bool pairSort (const tau_pair& pA, const tau_pair& pB)
{
  // first leg 1 iso
  float isoA = pA.iso1;
  float isoB = pB.iso1;
  if (isoA > isoB) return true; // NB: MVA iso ! Most iso --> highest MVA score
  else if (isoA < isoB) return false;

  // then leg 1 pt
  float ptA = pA.pt1;
  float ptB = pB.pt1;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // then leg 2 iso
  isoA = pA.iso2;
  isoB = pB.iso2;
  if (isoA > isoB) return true;
  else if (isoA < isoB) return false;

  // then leg 2 pt
  ptA = pA.pt2;
  ptB = pB.pt2;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // should be never here..
  return false;
}

// pairSorting for Lep+Tauh:
//  - lepton leg: lower score --> more isolated
//  - tauh leg: higher score --> more isolated
// Sorting strategy: iso leg1 -> pT leg1 -> iso leg2 -> pT leg2
bool pairSortHybrid (const tau_pair& pA, const tau_pair& pB)
{
  float isoA = pA.iso1;
  float isoB = pB.iso1;
  if (isoA < isoB) return true;
  else if (isoA < isoB) return false;

  // then leg 1 pt
  float ptA = pA.pt1;
  float ptB = pB.pt1;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // then leg 2 iso
  isoA = pA.iso2;
  isoB = pB.iso2;
  if (isoA > isoB) return true;
  else if (isoA < isoB) return false;

  // then leg 2 pt
  ptA = pA.pt2;
  ptB = pB.pt2;
  if (ptA > ptB) return true;
  else if (ptA < ptB) return false;

  // should be never here..
  return false;
}

// HHLeptonInterface class
class HHLeptonInterface {

  public:
    HHLeptonInterface (
      int vvvl_vsjet, int vl_vse, int vvl_vse, int t_vsmu, int vl_vsmu, // HPS DeepTau thresholds
      double BT_VsMu_threshold, double BT_VsE_threshold, double BT_VsJet_threshold // DeepBoostedTau thresholds
      );
    ~HHLeptonInterface ();
    // Trying to find suitable boosted taus (before looking at regular HPS taus)
    std::pair<lepton_output, cutflow_output> get_boosted_dau_indexes(
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
      sRVec boostedTau_genPartIdx,ROOT::VecOps::RVec<UChar_t> boostedTau_genPartFlav,
      iRVec boostedTau_muonCount, std::array<sRVec, 3> BT_muon_idx, 
      std::array<fRVec, 3> BT_muon_pt, std::array<fRVec, 3> BT_muon_correctedIso,
      iRVec boostedTau_electronCount, std::array<sRVec, 3> BT_electron_idx, 
      std::array<fRVec, 3> BT_electron_pt, std::array<fRVec, 3> BT_electron_correctedIso,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<trig_req> mutau_triggers, std::vector<trig_req> etau_triggers,
      std::vector<trig_req> tautau_triggers,
      int GenPairType, int genDau1_genPart_idx, int genDau2_genPart_idx
    );

    lepton_output get_dau_indexes(
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
      fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
      bRVec Muon_mediumId, bRVec Muon_tightId, iRVec Muon_charge,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
      bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
      bRVec Electron_mvaFall17V2Iso_WP90, fRVec Electron_pfRelIso03_all,
      fRVec Electron_dxy, fRVec Electron_dz, iRVec Electron_charge,
      fRVec Tau_pt, fRVec Tau_eta, fRVec Tau_phi, fRVec Tau_mass,
      iRVec Tau_idDeepTauVSmu, iRVec Tau_idDeepTauVSe,
      iRVec Tau_idDeepTauVSjet, fRVec Tau_rawDeepTauVSjet,
      fRVec Tau_dz, iRVec Tau_decayMode, iRVec Tau_charge,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<trig_req> mutau_triggers, std::vector<trig_req> etau_triggers,
      std::vector<trig_req> tautau_triggers, std::vector<trig_req> tautaujet_triggers, 
      std::vector<trig_req> vbf_triggers
    );

    bool lepton_veto(int muon_index, int electron_index,
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_dz, fRVec Muon_dxy,
      fRVec Muon_pfRelIso04_all, bRVec Muon_mediumId, bRVec Muon_tightId,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_dz, fRVec Electron_dxy,
      bRVec Electron_mvaFall17V2noIso_WP90, bRVec Electron_mvaFall17V2Iso_WP90,
      fRVec Electron_pfRelIso03_all);

    bool pass_trigger(
      float off_pt1, float off_eta1, float off_phi1, int obj_id1,
      float off_pt2, float off_eta2, float off_phi2, int obj_id2,
      std::vector<trig_req> triggers, 
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi);
    
    bool match_trigger_object(float off_eta, float off_phi, int obj_id, float trig_pt_threshold,
      iRVec TrigObj_id, iRVec TrigObj_filterBits, fRVec TrigObj_pt, fRVec TrigObj_eta, fRVec TrigObj_phi,
      std::vector<int> bits);

  private:
    // HPS taus DeepTau thresholds
    int vvvl_vsjet_;
    int vl_vse_;
    int vvl_vse_;
    int t_vsmu_;
    int vl_vsmu_;

    // boostedTaus DeepBoostedTaus threhols
    double BT_VsMu_threshold_;
    double BT_VsE_threshold_;
    double BT_VsJet_threshold_;
};


#endif // HHLeptonInterface_h