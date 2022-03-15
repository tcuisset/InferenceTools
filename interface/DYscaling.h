#ifndef DYscaling_h
#define DYscaling_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class DYscaling                                                                                          //
//                                                                                                                //
//                                                                                                                //
//   Author: Jaime Leon Holgado                                                                                   //
//   CIEMAT, January 2022                                                                                         //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

#include <vector>
#include <TLorentzVector.h>

#include <ROOT/RVec.hxx>

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

// DYscaling class
class DYscaling {

  public:
    DYscaling (int year);

    ~DYscaling ();

    std::vector<float> get_dy_scale(
      fRVec GenJet_pt,
      fRVec GenJet_eta,
      fRVec GenJet_phi,
      fRVec GenJet_mass,
      iRVec GenJet_hadronFlavour,
      int LHE_Nb,
      fRVec GenPart_pt,
      fRVec GenPart_eta,
      fRVec GenPart_phi,
      fRVec GenPart_mass,
      iRVec GenPart_statusFlags,
      iRVec GenPart_pdgId);

  private:
    // Computed in 2019 for 2018 data with deepFlavor - NOT USED for Legacy
    const float DYscale_LL[3] = {0.748154, 2.15445, 1.63619} ; // for now we use the same numbers computed with DY NLO sample
    const float DYscale_MM[3] = {0.862686, 1.08509, 1.10947} ; // for now we use the same numbers computed with DY NLO sample

    std::vector<float> DYscale_MH_vLowPt, DYscale_MH_LowPt, DYscale_MH_Med1Pt;
    std::vector<float> DYscale_MH_Med2Pt, DYscale_MH_HighPt, DYscale_MH_vHighPt;

    std::vector<float> DYscale_MH_vLowPt_2016 = {1.161, 0.515, 0.1};
    std::vector<float> DYscale_MH_LowPt_2016 = {1.151, 1.042, 1.150};
    std::vector<float> DYscale_MH_Med1Pt_2016 = {1.144, 1.152, 1.149};
    std::vector<float> DYscale_MH_Med2Pt_2016 = {1.151, 1.333, 1.218};
    std::vector<float> DYscale_MH_HighPt_2016 = {1.169, 1.458, 0.997};
    std::vector<float> DYscale_MH_vHighPt_2016 = {1.061, 1.963, 1.185};

    std::vector<float> DYscale_MH_vLowPt_2017 = {1.125, 0.01 , 0.01 };
    std::vector<float> DYscale_MH_LowPt_2017 = {1.326, 1.208, 1.016};
    std::vector<float> DYscale_MH_Med1Pt_2017 = {1.255, 1.317, 1.279};
    std::vector<float> DYscale_MH_Med2Pt_2017 = {1.198, 1.409, 1.374};
    std::vector<float> DYscale_MH_HighPt_2017 = {1.081, 1.687, 1.269};
    std::vector<float> DYscale_MH_vHighPt_2017 = {0.859, 1.595, 1.270};

    std::vector<float> DYscale_MH_vLowPt_2018 = {1.002, 0.01, 0.01};
    std::vector<float> DYscale_MH_LowPt_2018 = {1.209, 1.180, 1.162};
    std::vector<float> DYscale_MH_Med1Pt_2018 = {1.17 , 1.260, 1.288};
    std::vector<float> DYscale_MH_Med2Pt_2018 = {1.140, 1.357, 1.574};
    std::vector<float> DYscale_MH_HighPt_2018 = {1.037, 1.440, 1.603};
    std::vector<float> DYscale_MH_vHighPt_2018 = {0.835, 1.994, 1.037};

    std::vector<float> DYscale_MTT_vLowPt, DYscale_MTT_LowPt, DYscale_MTT_Med1Pt;
    std::vector<float> DYscale_MTT_Med2Pt, DYscale_MTT_HighPt, DYscale_MTT_vHighPt;

    std::vector<float> DYscale_MTT_vLowPt_2016 = {1.1630144, 0.010000393, 0.010000000};
    std::vector<float> DYscale_MTT_LowPt_2016 = {1.2194740, 1.1250249, 1.0609708};
    std::vector<float> DYscale_MTT_Med1Pt_2016 = {1.2536864, 1.2376837, 1.1901911};
    std::vector<float> DYscale_MTT_Med2Pt_2016 = {1.2763251, 1.2972053, 1.2731480};
    std::vector<float> DYscale_MTT_HighPt_2016 = {1.2785250, 1.4578434, 1.2241989};
    std::vector<float> DYscale_MTT_vHighPt_2016 = {1.1649714, 1.6778047, 1.1510545};

    std::vector<float> DYscale_MTT_vLowPt_2017 = {1.1092067, 0.010000025, 0.010000002};
    std::vector<float> DYscale_MTT_LowPt_2017 = {1.4229782, 0.025083862, 0.64151569};
    std::vector<float> DYscale_MTT_Med1Pt_2017 = {1.3143450, 1.0076030, 0.96584965};
    std::vector<float> DYscale_MTT_Med2Pt_2017 = {1.2485879, 1.2376612, 1.0745893};
    std::vector<float> DYscale_MTT_HighPt_2017 = {1.1273438, 1.5564765, 1.0578688};
    std::vector<float> DYscale_MTT_vHighPt_2017 = {0.87293424, 1.4469675, 1.1665250};
    
    std::vector<float> DYscale_MTT_vLowPt_2018 = {0.87720949, 0.010000006, 0.010000025};
    std::vector<float> DYscale_MTT_LowPt_2018 = {1.2191486, 0.010001064, 0.29051790 };
    std::vector<float> DYscale_MTT_Med1Pt_2018 = {1.1816037, 0.82760074, 0.84809836 };
    std::vector<float> DYscale_MTT_Med2Pt_2018 = {1.1579303, 1.1240148, 0.92974364 };
    std::vector<float> DYscale_MTT_HighPt_2018 = {1.0469869, 1.3690206, 1.0024774  };
    std::vector<float> DYscale_MTT_vHighPt_2018 = {0.80838089, 1.7465338, 0.73211715 };

    std::vector<float> DYscale_MTT_vLowPt_err, DYscale_MTT_LowPt_err, DYscale_MTT_Med1Pt_err;
    std::vector<float> DYscale_MTT_Med2Pt_err, DYscale_MTT_HighPt_err, DYscale_MTT_vHighPt_err;

    std::vector<float> DYscale_MTT_vLowPt_err_2016 = {0.0037446320, 0.011371985, 0.0071346649};
    std::vector<float> DYscale_MTT_LowPt_err_2016 = {0.0017404799, 0.023737485, 0.035623816};
    std::vector<float> DYscale_MTT_Med1Pt_err_2016 = {0.0027420531, 0.034357252, 0.043622931};
    std::vector<float> DYscale_MTT_Med2Pt_err_2016 = {0.0042018892, 0.055512921, 0.061654929};
    std::vector<float> DYscale_MTT_HighPt_err_2016 = {0.0045670912, 0.065499641, 0.064115483};
    std::vector<float> DYscale_MTT_vHighPt_err_2016 = {0.0084286985, 0.17367880, 0.13585326};
    
    std::vector<float> DYscale_MTT_vLowPt_err_2017 = {0.0020977215, 0.0014144677, 0.0015785401};
    std::vector<float> DYscale_MTT_LowPt_err_2017 = {0.0018789754, 0.028074073, 0.033416855};
    std::vector<float> DYscale_MTT_Med1Pt_err_2017 = {0.0022982702, 0.028221735, 0.032660541};
    std::vector<float> DYscale_MTT_Med2Pt_err_2017 = {0.0019811270, 0.024603319, 0.024254387};
    std::vector<float> DYscale_MTT_HighPt_err_2017 = {0.0030012172, 0.046576904, 0.035949714};
    std::vector<float> DYscale_MTT_vHighPt_err_2017 = {0.0061768066, 0.15213682, 0.091582069};

    std::vector<float> DYscale_MTT_vLowPt_err_2018 = {0.0014349486, 0.00085205781, 0.00094359815};
    std::vector<float> DYscale_MTT_LowPt_err_2018 = {0.0013417125, 0.0044914533, 0.025179988};
    std::vector<float> DYscale_MTT_Med1Pt_err_2018 = {0.0016262942, 0.021066553, 0.025609469};
    std::vector<float> DYscale_MTT_Med2Pt_err_2018 = {0.0014769770, 0.018984114, 0.019089746};
    std::vector<float> DYscale_MTT_HighPt_err_2018 = {0.0022691835, 0.035830965, 0.027803721};
    std::vector<float> DYscale_MTT_vHighPt_err_2018 = {0.0046517215, 0.11443309 , 0.064519398};
};

#endif // DYscaling_h