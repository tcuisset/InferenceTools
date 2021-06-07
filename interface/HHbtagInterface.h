#ifndef HHbtagInterface_h
#define HHbtagInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHbtagInterface                                                                                        //
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


// HHbtagInterface class
class HHbtagInterface {

  public:
    HHbtagInterface (std::string model_0, std::string model_1, int year);
    ~HHbtagInterface ();

    std::vector<float> GetScore(
      std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
      std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_deepFlavour_,
      std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
      float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_);

  private:
    hh_btag::HH_BTag HHbtagger_;
};

#endif // HHbtagInterface