#include "Tools/Tools/interface/HHKinFitInterface.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHInvMConstraintException.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2Core/interface/exceptions/HHEnergyConstraintException.h"
#include <vector>

HHKinFitInterface::HHKinFitInterface(TLorentzVector const& bjet1, 
    TLorentzVector const& bjet2, 
    TLorentzVector const& tauvis1,
    TLorentzVector const& tauvis2,
    TVector2 const& met, 
    TMatrixD const& met_cov, 
    double sigmaEbjet1,
    double sigmaEbjet2,
    bool istruth, 
    TLorentzVector const&  heavyhiggsgen):
    m_HHKinFit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov, sigmaEbjet1, sigmaEbjet2, istruth, heavyhiggsgen)
{};

// Destructor
HHKinFitInterface::~HHKinFitInterface() {}

void HHKinFitInterface::addHypo(int mh1, int mh2) {
  m_HHKinFit.addHypo(mh1, mh2);
  return;
}

std::vector<double> HHKinFitInterface::fit(){
  bool good_fit = true;
  try {
    m_HHKinFit.fit();
  }
  catch (HHKinFit2::HHInvMConstraintException &e) {
    good_fit = false;
  }
  catch (HHKinFit2::HHEnergyRangeException &e) {
    good_fit = false;
  }
  catch (HHKinFit2::HHEnergyConstraintException &e) {
    good_fit = false;
  }
  
  if (good_fit) return {m_HHKinFit.getMH(), m_HHKinFit.getChi2()};
  else return {-999., -999.};
  
}