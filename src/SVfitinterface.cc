#include "Tools/Tools/interface/SVfitinterface.h"


// Constructors

SVfitinterface::SVfitinterface () {};

// Destructor
SVfitinterface::~SVfitinterface() {}

std::vector<double> SVfitinterface::FitAndGetResultWithInputs(
    int verbosity, int pairType, int DM1, int DM2,
    double tau1_pt, double tau1_eta, double tau1_phi, double tau1_mass,
    double tau2_pt, double tau2_eta, double tau2_phi, double tau2_mass,
    double met_pt, double met_phi,
    double met_covXX, double met_covXY, double met_covYY){

  int verbosity_;
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
  double METx_;
  double METy_;
  TMatrixD covMET_(2, 2);
  double kappa_;

  // verbosity
  verbosity_ = verbosity;
  measuredTauLeptons_.clear();
  // MET
  covMET_(0,0) = met_covXX;
  covMET_(0,1) = met_covXY;
  covMET_(1,0) = met_covXY;
  covMET_(1,1) = met_covYY;
  METx_ = met_pt * cos(met_phi);
  METy_ = met_pt * sin(met_phi);

  // Leptons
  double mass1, mass2;
  int decay1, decay2;
  classic_svFit::MeasuredTauLepton::kDecayType l1Type, l2Type;

  if (pairType == 0) // MuTau
  {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    mass1  = 105.658e-3;
    decay1 = -1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2;
    kappa_ = 4.;
  }

  else if (pairType == 1) // EleTau
  {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    mass1  = 0.51100e-3;
    decay1 = -1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2;
    kappa_ = 4.;
  }

  else // TauTau
  {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass1  = tau1_mass;
    decay1 = DM1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2; 
    kappa_ = 5.;
  }

  // Fill the measuredTauLeptons
  measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l1Type, tau1_pt, tau1_eta, tau1_phi, mass1, decay1));
  measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l2Type, tau2_pt, tau2_eta, tau2_phi, mass2, decay2));

  // Declare result: vector of SVfit (Pt,Eta,Phi,Mass)
  std::vector<double> result(4,-999.);

  // Declare algo
  ClassicSVfit algo(verbosity_);

  // Configure more options
  algo.addLogM_fixed(false, kappa_);
  algo.addLogM_dynamic(false);

  // Actually integrate
  algo.integrate(measuredTauLeptons_, METx_, METy_, covMET_);

  // Return SVfit quantities if the integration succeeded 
  // otherwise vector of -999 if the integration failed
  if (algo.isValidSolution())
  {
    result.at(0) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt();
    result.at(1) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
    result.at(2) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
    result.at(3) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass();
  }

  return result;

}
