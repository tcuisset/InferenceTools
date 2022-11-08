#ifndef HHMulticlassInterface_h
#define HHMulticlassInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHMulticlassInterface                                                                                  //
//                                                                                                                //
//   Class to compute Multi-class DNN outputs.                                                                    //
//                                                                                                                //
//   Author: Jaime León Holgado (CIEMAT)                                                                          //
//   Date  : November 2022                                                                                        //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>

#include "Tools/Tools/interface/MulticlassInterface.h"

// Math libraries
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>

// ROOT libraries
#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>

// Using names
using DNNVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

// HHMulticlassInterface class
class HHMulticlassInterface {

  public:
    HHMulticlassInterface (int year, std::vector<std::pair<std::string, std::string>> modelSpecs);
    ~HHMulticlassInterface ();

    std::vector<std::vector<float>> GetPredictionsWithInputs(
      int EventNumber, int pairType,
      fRVec jet_pt, fRVec jet_eta, fRVec jet_phi, fRVec jet_mass,
      fRVec jet_btagDeepFlavB, fRVec jet_btagDeepFlavCvL, fRVec jet_btagDeepFlavCvB, fRVec jet_HHbtag,
      int dau1_index, int dau2_index,
      int bjet1_JetIdx, int bjet2_JetIdx, int VBFjet1_JetIdx, int VBFjet2_JetIdx,
      iRVec ctjet_indexes, iRVec fwjet_indexes,
      fRVec muon_pt, fRVec muon_eta, fRVec muon_phi, fRVec muon_mass,
      fRVec electron_pt, fRVec electron_eta, fRVec electron_phi, fRVec electron_mass,
      fRVec tau_pt, fRVec tau_eta, fRVec tau_phi, fRVec tau_mass,
      float met_pt, float met_phi,
      float htt_sv_pt, float htt_sv_eta, float htt_sv_phi, float htt_sv_mass
    );

    std::vector<std::string> get_node_names(size_t imodel);

  private:
    MulticlassInterface mci_;
    int year_;
};

#endif // HHMulticlassInterface_h