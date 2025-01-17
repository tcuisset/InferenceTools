#include "Tools/Tools/interface/HHDNNinterfaceNew.h"

#include <Math/Vector4Dfwd.h>

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





void get_dnn_inputs(HHDNNinterfaceNew& dnnInt, int pairType, int isBoosted, int isBoostedTau, ULong64_t event,
    int bjet1_index, int bjet2_index,
    int fatjet_index, int vbfjet1_index, int vbfjet2_index,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    ROOT::RVec<float> const& jet_pt, ROOT::RVec<float> const& jet_eta, ROOT::RVec<float> const& jet_phi, ROOT::RVec<float> const& jet_mass,
    ROOT::RVec<float> const& fatjet_pt, ROOT::RVec<float> const& fatjet_eta, ROOT::RVec<float> const& fatjet_phi, ROOT::RVec<float> const& fatjet_mass, ROOT::RVec<float> const& fatjet_msoftdrop,
    float htt_sv_pt, float htt_sv_eta, float htt_sv_phi, float htt_sv_mass,
    float HHKinFit_mass, float HHKinFit_chi2, float met_pt, float met_phi,
    ROOT::RVec<float> const& Jet_btagDeepFlavB, ROOT::RVec<float> const& Jet_btagDeepFlavCvL, ROOT::RVec<float> const& Jet_btagDeepFlavCvB,
    ROOT::RVec<float> const& Jet_HHbtag
)
{          
    using LVector = ROOT::Math::PtEtaPhiMVector;

    dnnInt.pairType = pairType;
    dnnInt.isBoosted = isBoosted;
    dnnInt.isBoostedTau = isBoostedTau;

    dnnInt.l1 = LVector(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
    dnnInt.l2 = LVector(dau2_pt, dau2_eta, dau2_phi, dau2_mass);

    dnnInt.SVfitConv = htt_sv_mass >= 0;
    
    dnnInt.met.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

    if (htt_sv_mass > 0)
        dnnInt.svfit = LVector(htt_sv_pt, htt_sv_eta, htt_sv_phi, htt_sv_mass);
    else
        dnnInt.svfit = LVector(1., 1., 1., 1.);
    
    dnnInt.KinFitMass = HHKinFit_mass;
    dnnInt.KinFitChi2 = HHKinFit_chi2;
    if (!isBoosted) {
        assert (bjet1_index>=0 && bjet2_index >= 0);
        dnnInt.b1 = LVector(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
            jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
        dnnInt.b2 = LVector(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
            jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));
        
        dnnInt.deepFlav1 = Jet_btagDeepFlavB.at(bjet1_index);
        dnnInt.deepFlav2 = Jet_btagDeepFlavB.at(bjet2_index);
        dnnInt.CvsL_b1 = Jet_btagDeepFlavCvL.at(bjet1_index);
        dnnInt.CvsL_b2 = Jet_btagDeepFlavCvL.at(bjet2_index);
        dnnInt.CvsB_b1 = Jet_btagDeepFlavCvB.at(bjet1_index);
        dnnInt.CvsB_b2 = Jet_btagDeepFlavCvB.at(bjet2_index);
        assert (Jet_HHbtag.size() == jet_eta.size());
        dnnInt.HHbtag_b1 = Jet_HHbtag.at(bjet1_index);
        dnnInt.HHbtag_b2 = Jet_HHbtag.at(bjet2_index);

        dnnInt.KinFitConv = HHKinFit_chi2 >= 0;
                        
        dnnInt.fatjet = LVector(0., 0., 0., 0.);
        dnnInt.fatjet_softDropMass = 0.;

    }
    else {
        assert (fatjet_index>=0);
        dnnInt.fatjet = LVector(fatjet_pt.at(fatjet_index), fatjet_eta.at(fatjet_index),
            fatjet_phi.at(fatjet_index), fatjet_mass.at(fatjet_index));
        dnnInt.fatjet_softDropMass = fatjet_msoftdrop.at(fatjet_index);
                        
        dnnInt.b1 = LVector(0., 0., 0., 0.);
        dnnInt.b2 = LVector(0., 0., 0., 0.);
                        
        dnnInt.deepFlav1 = 0.;
        dnnInt.deepFlav2 = 0.;
        dnnInt.CvsL_b1 = 0.;
        dnnInt.CvsL_b2 = 0.;
        dnnInt.CvsB_b1 = 0.;
        dnnInt.CvsB_b2 = 0.;
        dnnInt.HHbtag_b1 = 0.;
        dnnInt.HHbtag_b2 = 0.;

        dnnInt.KinFitConv = false; // No KinFit for boosted bb
    }
        


    // MT2 computation
    /*
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    double MT2 = asymm_mt2_lester_bisect::get_mT2(
        bjet1_tlv.M(), bjet1_tlv.Px(), bjet1_tlv.Py(),
        bjet2_tlv.M(), bjet2_tlv.Px(), bjet2_tlv.Py(),
        dau1_tlv.Px() + dau2_tlv.Px() + met_tlv.Px(),
        dau1_tlv.Py() + dau2_tlv.Py() + met_tlv.Py(),
        dau1_tlv.M(), dau2_tlv.M(), 0.);
    */
    dnnInt.computeIntermediateVariables();
                        
    //return dnnInt;
}