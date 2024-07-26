#include "Tools/Tools/interface/HHDNNinterfaceNew.h"

// Constructor
HHDNNinterfaceNew::HHDNNinterfaceNew (int year_, std::vector<float> bjet_wps_, std::vector<float> cvsl_wps_, std::vector<float> cvsb_wps_)
 : year_((unsigned short) (year_ - 2016)), bjet_wps(bjet_wps_), cvsl_wps(cvsl_wps_), cvsb_wps(cvsb_wps_)
{
}

float HHDNNinterfaceNew::delta_r(const LorentzVector& v_0, const LorentzVector& v_1) {return ROOT::Math::VectorUtil::DeltaR(v_0, v_1);}

float HHDNNinterfaceNew::delta_r_boosted(const LorentzVector& v_0, const LorentzVector& v_1, const LorentzVector& ref){
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return DeltaR(boost(v_0, ref.BoostToCM()), boost(v_1, ref.BoostToCM()));
}

float HHDNNinterfaceNew::delta_phi(const LorentzVector& v_0, const LorentzVector& v_1) {return std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_0, v_1));}

float HHDNNinterfaceNew::calc_cos_delta(const LorentzVector& v, const LorentzVector& r) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return  CosTheta(boost(v, r.BoostToCM()), r);
}

float HHDNNinterfaceNew::calc_phi(const LorentzVector& l_1, const LorentzVector& l_2,
                                const LorentzVector& b_1, const LorentzVector& b_2, const LorentzVector& hh) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return Angle(boost(l_1, hh.BoostToCM()).Vect().Cross(boost(l_2, hh.BoostToCM()).Vect()),
                 boost(b_1, hh.BoostToCM()).Vect().Cross(boost(b_2, hh.BoostToCM()).Vect()));
}

float HHDNNinterfaceNew::calc_mt(const LorentzVector& v, const LorentzVector& met) {
    return std::sqrt(2.0*v.Pt()*met.Pt()*(1.0-std::cos(ROOT::Math::VectorUtil::DeltaPhi(v,met))));
}



void HHDNNinterfaceNew::computeIntermediateVariables() {
  h_bb = isBoosted ? fatjet : b1 + b2;
  h_tt_met = l1 + l2 + met;
  if (isBoosted || !KinFitConv) {  // HHKinFit didn't converge
      hh = SVfitConv ? h_bb+svfit : h_bb+h_tt_met;
  } else if (!SVfitConv) {  // HHKinFit converge but SVFit didn't
      hh = LorentzVector(h_bb.Px()+h_tt_met.Px(), h_bb.Py()+h_tt_met.Py(), h_bb.Pz()+h_tt_met.Pz(), KinFitMass);
  } else {
     hh = LorentzVector(h_bb.Px()+svfit.Px(), h_bb.Py()+svfit.Py(), h_bb.Pz()+svfit.Pz(), KinFitMass); 
  }
}

short HHDNNinterfaceNew::get_cvsb_flag(float score) const {
    short tag(0);
    for (float wp : cvsb_wps) if (score <= wp) tag++;
    return tag;
}

short HHDNNinterfaceNew::get_cvsl_flag(float score) const {
    short tag(0);
    for (float wp : cvsl_wps) if (score >= wp) tag++;
    return tag;
}
