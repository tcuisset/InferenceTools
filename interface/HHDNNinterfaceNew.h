#ifndef HHDNNinterfaceNew_h
#define HHDNNinterfaceNew_h

#include <vector>
#include <string>
#include <cmath>

#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>
#include <ROOT/RVec.hxx>


class HHDNNinterfaceNew {

  public:
    HHDNNinterfaceNew() = default;
    /* Parameters : 
     - year =  "regular" year ie 2016, ...
     - bjet_wps : working points of b-tagging algorithm
     - cvsl_wps : WPs of c-tagging for CvsL
     - cvsb_wps : CvsB
    */
    HHDNNinterfaceNew (int year, std::vector<float> bjet_wps, std::vector<float> cvsl_wps, std::vector<float> cvsb_wps);

    void computeIntermediateVariables();

    unsigned short year() const { return year_; }; // Year for DNN input (ie 0, 1, 2)
    bool boosted_bb() const { return isBoosted; }
    bool boostedTau() const { return isBoostedTau; }
    short channel() const { 
      // pairType | Our Def. | DNN
      // ------------------------
      // mutau   |    0     |  1
      // etau    |    1     |  2
      // tautau  |    2     |  0

      if(pairType == 0)
        return 1;
      else if (pairType == 1)
        return 2;
      else
        return 0;
    }
    bool is_vbf() const { return false; /* VBF is to be implemented */ }

    float hh_kinfit_m() const { return KinFitConv ? KinFitMass : 0.; }
    float hh_kinfit_chi2() const { return KinFitConv ? KinFitChi2 : 0.; }

    float sv_mass() const { return SVfitConv ? svfit.M() : std::nanf("1"); }

    float l_1_pT() const { return l1.Pt(); }
    float l_2_pT() const { return l2.Pt(); }
    float l_1_mt() const { return calc_mt(l1, met); }
    float l_2_mt() const { return calc_mt(l2, met); }

    float b_1_pT() const { return b1.Pt(); }
    float b_2_pT() const { return b2.Pt(); }

    float dR_l1_l2() const { return delta_r(l1, l2); }
    float dphi_sv_met() const { return SVfitConv ? delta_phi(svfit, met) : std::nanf("1"); }
    
    float h_bb_mass() const { return h_bb.M(); } /* resolved->(b1+b2).M() boosted->FatJet.M() */

    float b_2_hhbtag() const { return HHbtag_b2; }; /* 0 for boosted */
    short b_1_cvsb() const { return isBoosted ? 0 : get_cvsb_flag(CvsB_b1); }; 
    short b_1_cvsl() const { return isBoosted ? 0 : get_cvsl_flag(CvsL_b1); };
    short jet_1_quality() const {
      if (isBoosted) return 0;
      short tag_1(0);
      for (float wp : bjet_wps ) {
          if (deepFlav1 >= wp) tag_1++;
      }
      return tag_1;
    }
    short jet_2_quality() const {
      if (isBoosted) return 0;
      short tag_2(0);
      for (float wp : bjet_wps ) {
          if (deepFlav2 >= wp) tag_2++;
      }
      return tag_2;
    }

    float diH_mass_sv() const { return SVfitConv ? (h_bb+svfit).M() : std::nanf("1"); };
    float dphi_hbb_sv() const { return SVfitConv ? delta_phi(h_bb, svfit) : std::nanf("1"); };
    float h_bb_pT() const { return h_bb.Pt(); };

    float dR_l1_l2_x_sv_pT() const { return SVfitConv ? delta_r(l1, l2)*svfit.Pt() : std::nanf("1"); }
    float dR_l1_l2_boosted_htt_met() const { return delta_r_boosted(l1, l2, h_tt_met); };
    float phi() const { return isBoosted ? 0. : calc_phi(l1, l2, b1, b2, hh); }; /* The angle between the decay plane of the taus and the decay plane of the jets, expressed in the rest-frame of the ZZ/ZH system.*/
    float costheta_l2_httmet() const { return calc_cos_delta(l2, h_tt_met); };

  private:
    using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;
    short get_cvsb_flag(float score) const;
    short get_cvsl_flag(float score) const;

    static float delta_r(const LorentzVector& v_0, const LorentzVector& v_1);
    static float delta_r_boosted(const LorentzVector& v_0, const LorentzVector& v_1, const LorentzVector& ref);
    static float delta_phi(const LorentzVector& v_0, const LorentzVector& v_1);
    static float calc_cos_delta(const LorentzVector& v, const LorentzVector& r);
    static float calc_phi(const LorentzVector& l_1, const LorentzVector& l_2,
                                const LorentzVector& b_1, const LorentzVector& b_2, const LorentzVector& hh);
    static float calc_mt(const LorentzVector& v, const LorentzVector& met);


  public:
    unsigned short year_;
    short pairType;
    bool isBoosted; // boosted-bb (FatJet) versus resolved (2*AK4)
    bool isBoostedTau; // boostedTau vs HPS tau reconstruction
    unsigned long long int eventn;

    // 4-vectors of b-jets (for resolved b)
    LorentzVector b1;
    LorentzVector b2;
    LorentzVector fatjet; // FatJet with raw jet mass (for boosted bb)
    float fatjet_softDropMass; //Maybe could use ParticleNet regressed mass ?
    // float fatjet_nsubjettiness // possible future boosted DNN feature

    LorentzVector l1; // leptons
    LorentzVector l2;
    
    LorentzVector met;
    LorentzVector svfit; // tautau 4-vector from SVFit
    
    float KinFitMass;
    float KinFitChi2;
    bool KinFitConv; // did kinematic fit converge ?
    bool SVfitConv;
    float deepFlav1; // b-tagging scores of 2 AK4 jets (only for resolved category)
    float deepFlav2;
    float CvsL_b1;
    float CvsL_b2;
    float CvsB_b1;
    float CvsB_b2;
    float HHbtag_b1; // HHbtag score (resolved only)
    float HHbtag_b2;

    float DNN_res_mass; // Resonant mass for parametrized DNN
  
  private: // Variables computed as a function of the public attributes
    LorentzVector h_bb; // b1 + b2 (resolved) / FatJet (boosted)
    LorentzVector h_tt_met; // l1 + l2 + met
    LorentzVector hh; // mass of hh, using KinFit+SVFit if available (else fallbacks)

    std::vector<float> bjet_wps;
    std::vector<float> cvsl_wps;
    std::vector<float> cvsb_wps;
};

void get_dnn_inputs(HHDNNinterfaceNew& dnnInt, int pairType, int isBoosted, int isBoostedTau, ULong64_t event,
    int bjet1_index, int bjet2_index,
    int fatjet_index, int vbfjet1_index, int vbfjet2_index,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float bjet1_pt, float bjet1_eta, float bjet1_phi, float bjet1_mass,
    float bjet2_pt, float bjet2_eta, float bjet2_phi, float bjet2_mass,
    float fatjet_pt, float fatjet_eta, float fatjet_phi, float fatjet_mass, float fatjet_msoftdrop,
    float htt_sv_pt, float htt_sv_eta, float htt_sv_phi, float htt_sv_mass,
    float HHKinFit_mass, float HHKinFit_chi2, float met_pt, float met_phi,
    float bjet1_btagDeepFlavB, float bjet1_btagDeepFlavCvL, float bjet1_btagDeepFlavCvB,
    float bjet2_btagDeepFlavB, float bjet2_btagDeepFlavCvL, float bjet2_btagDeepFlavCvB,
    float bjet1_HHbtag, float bjet2_HHbtag
);

#endif // HHDNNinterface_h