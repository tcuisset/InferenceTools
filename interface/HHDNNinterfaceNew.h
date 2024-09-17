#ifndef HHDNNinterfaceNew_h
#define HHDNNinterfaceNew_h

#include <vector>
#include <string>
#include <cmath>

#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>


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
    bool boosted() const { return isBoosted; }
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

    float hh_kinfit_m() const { return KinFitConv ? KinFitMass : 0; }
    float hh_kinfit_chi2() const { return KinFitChi2; }

    float sv_mass() const { return SVfitConv ? svfit.M() : std::nanf("1"); }

    float l_1_pT() const { return l1.Pt(); }
    float l_2_pT() const { return l2.Pt(); }
    float l_1_mt() const { return calc_mt(l1, met); }
    float l_2_mt() const { return calc_mt(l2, met); }

    float b_1_pT() const { return b1.Pt(); }
    float b_2_pT() const { return b2.Pt(); }

    float dR_l1_l2() const { return delta_r(l1, l2); }
    float dphi_sv_met() const { return SVfitConv ? delta_phi(svfit, met) : std::nanf("1"); }
    
    float h_bb_mass() const { return h_bb.M(); }

    float b_2_hhbtag() const { return HHbtag_b2; };
    short b_1_cvsb() const { return get_cvsb_flag(CvsB_b1); };
    short b_1_cvsl() const { return get_cvsl_flag(CvsL_b1); };
    short jet_1_quality() const {
      short tag_1(0);
      for (float wp : bjet_wps ) {
          if (deepFlav1 >= wp) tag_1++;
      }
      return tag_1;
    }
    short jet_2_quality() const {
      short tag_2(0);
      for (float wp : bjet_wps ) {
          if (deepFlav2 >= wp) tag_2++;
      }
      return tag_2;
    }

    float diH_mass_sv() const { return SVfitConv ? (h_bb+svfit).M() : std::nanf("1"); };
    float dphi_hbb_sv() const { return SVfitConv ? delta_phi(h_bb, svfit) : std::nanf("1"); };
    float h_bb_pT() const { return delta_r(b1, b2)*h_bb.Pt(); };

    float dR_l1_l2_x_sv_pT() const { return SVfitConv ? delta_r(l1, l2)*svfit.Pt() : std::nanf("1"); }
    float dR_l1_l2_boosted_htt_met() const { return delta_r_boosted(l1, l2, h_tt_met); };
    float phi() const { return calc_phi(l1, l2, b1, b2, hh); };
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
    bool isBoosted;
    unsigned long long int eventn;

    // 4-vectors of b-jets and of leptons
    LorentzVector b1;
    LorentzVector b2;
    LorentzVector fatjet; // FatJet with raw jet mass
    float fatjet_softDropMass;
    // float fatjet_nsubjettiness // possible future boosted DNN feature

    LorentzVector l1;
    LorentzVector l2;
    
    LorentzVector met;
    LorentzVector svfit;
    
    float KinFitMass;
    float KinFitChi2;
    bool KinFitConv;
    bool SVfitConv;
    float deepFlav1;
    float deepFlav2;
    float CvsL_b1;
    float CvsL_b2;
    float CvsB_b1;
    float CvsB_b2;
    float HHbtag_b1;
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

#endif // HHDNNinterface_h