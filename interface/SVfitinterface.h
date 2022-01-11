#ifndef SVfitinterface_h
#define SVfitinterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class SVfitinterface                                                                                         //
//                                                                                                                //
//   Class to compute ClassicSVfit. Modified from                                                                 //
//    - https://github.com/LLRCMS/KLUBAnalysis/blob/VBF_legacy/src/SVfitKLUBinterface.cc                          //
//                                                                                                                //
//   Based on:                                                                                                    //
//    - https://github.com/LLRCMS/LLRHiggsTauTau/blob/102X_HH/NtupleProducer/plugins/ClassicSVfitInterface.cc     //
//                                                                                                                //
//   Original SVfit package from:                                                                                 //
//    - https://github.com/SVfit/ClassicSVfit                                                                     //
//                                                                                                                //
//   Author: Francesco Brivio (Milano-Bicocca)                                                                    //
//   Date  : April 2020                                                                                           //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <string>
#include <cmath>

// ROOT libraries
#include <TLorentzVector.h>
#include "TMatrixD.h"

// ClassicSVfit libraries
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

using namespace classic_svFit;

// SVfitinterface class
class SVfitinterface {

  public:
  SVfitinterface ();
    SVfitinterface (
      int verbosity, int pairType, int DM1, int DM2,
      Float_t tau1_pt, Float_t tau1_eta, Float_t tau1_phi, Float_t tau1_mass,
      Float_t tau2_pt, Float_t tau2_eta, Float_t tau2_phi, Float_t tau2_mass,
      Float_t met_pt, Float_t met_phi,
      Float_t met_covXX, Float_t met_covXY, Float_t met_covYY
    );
    ~SVfitinterface ();

    void SetInputs(int verbosity, int pairType, int DM1, int DM2,
      Float_t tau1_pt, Float_t tau1_eta, Float_t tau1_phi, Float_t tau1_mass,
      Float_t tau2_pt, Float_t tau2_eta, Float_t tau2_phi, Float_t tau2_mass,
      Float_t met_pt, Float_t met_phi,
      Float_t met_covXX, Float_t met_covXY, Float_t met_covYY);
    std::vector<double> FitAndGetResult();

  private:
    int verbosity_;
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
    double METx_;
    double METy_;
    TMatrixD covMET_;
    double kappa_;

};

#endif // SVfitinterface_h