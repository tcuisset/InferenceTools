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

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

struct tau_pair {
  int index1;
  int iso1;
  float pt1;
  int index2;
  int iso2;
  float pt2;
  int isVBFtrigger;
};

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

// HHLeptonInterface class
class HHLeptonInterface {

  public:
    HHLeptonInterface ();
    ~HHLeptonInterface ();
    std::vector<int> get_dau_indexes(
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_phi, fRVec Muon_mass,
      fRVec Muon_pfRelIso04_all, fRVec Muon_dxy, fRVec Muon_dz,
      bRVec Muon_mediumId, bRVec Muon_tightId,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_phi, fRVec Electron_mass,
      bRVec Electron_mvaFall17V2Iso_WP80, bRVec Electron_mvaFall17V2noIso_WP90,
      bRVec Electron_mvaFall17V2Iso_WP90, fRVec Electron_pfRelIso03_all,
      fRVec Electron_dxy, fRVec Electron_dz,
      fRVec Tau_pt, fRVec Tau_eta, fRVec Tau_phi, fRVec Tau_mass,
      iRVec Tau_idDeepTau2017v2p1VSmu, iRVec Tau_idDeepTau2017v2p1VSe,
      iRVec Tau_idDeepTau2017v2p1VSjet, fRVec Tau_rawDeepTau2017v2p1VSjet,
      iRVec Tau_dz, iRVec Tau_decayMode    
    ); // triggers missing
    bool lepton_veto(int muon_index, int electron_index,
      fRVec Muon_pt, fRVec Muon_eta, fRVec Muon_dz, fRVec Muon_dxy,
      fRVec Muon_pfRelIso04_all, bRVec Muon_mediumId, bRVec Muon_tightId,
      fRVec Electron_pt, fRVec Electron_eta, fRVec Electron_dz, fRVec Electron_dxy,
      bRVec Electron_mvaFall17V2noIso_WP90, bRVec Electron_mvaFall17V2Iso_WP90,
      fRVec Electron_pfRelIso03_all);

  private:

};

#endif // HHLeptonInterface_h