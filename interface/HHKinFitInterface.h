#ifndef HHKinFitInterface_h
#define HHKinFitInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class HHKinFitInterface                                                                                         //
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
#include <TMatrixD.h>
#include <TVector2.h>

// HHKinFit2 libraries
#include "HHKinFit2/HHKinFit2Scenarios/interface/HHKinFitMasterHeavyHiggs.h"

// HHKinFitInterface class
class HHKinFitInterface {

  public:
    HHKinFitInterface (
      TLorentzVector const& bjet1,
      TLorentzVector const& bjet2,
      TLorentzVector const& tauvis1,
      TLorentzVector const& tauvis2,
      TVector2 const& met, TMatrixD const& met_cov,
      double sigmaEbjet1 = -1.0, double sigmaEbjet2 = -1.0,
      bool istruth=false,
      TLorentzVector const& higgsgen=TLorentzVector(0,0,0,0)
    );
    ~HHKinFitInterface ();

    std::vector<double> fit();
    void addHypo(int mh1, int mh2);

  private:
    HHKinFit2::HHKinFitMasterHeavyHiggs m_HHKinFit;
    int m_mh1; // store the values of the hypothesis mass (only one hypo supported) for the call to getMH and getChi2
    int m_mh2;
};

#endif // HHKinFitInterface_h