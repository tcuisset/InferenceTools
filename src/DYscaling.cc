#include "Tools/Tools/interface/DYscaling.h"

// Constructors

DYscaling::DYscaling (int year) {
  if (year == 2016) {
    DYscale_MH_vLowPt = DYscale_MH_vLowPt_2016;
    DYscale_MH_LowPt = DYscale_MH_LowPt_2016;
    DYscale_MH_Med1Pt = DYscale_MH_Med1Pt_2016;
    DYscale_MH_Med2Pt = DYscale_MH_Med2Pt_2016;
    DYscale_MH_HighPt = DYscale_MH_HighPt_2016;
    DYscale_MH_vHighPt = DYscale_MH_vHighPt_2016;

    DYscale_MTT_vLowPt = DYscale_MTT_vLowPt_2016;
    DYscale_MTT_LowPt = DYscale_MTT_LowPt_2016;
    DYscale_MTT_Med1Pt = DYscale_MTT_Med1Pt_2016;
    DYscale_MTT_Med2Pt = DYscale_MTT_Med2Pt_2016;
    DYscale_MTT_HighPt = DYscale_MTT_HighPt_2016;
    DYscale_MTT_vHighPt = DYscale_MTT_vHighPt_2016;

    DYscale_MTT_vLowPt_err = DYscale_MTT_vLowPt_err_2016;
    DYscale_MTT_LowPt_err = DYscale_MTT_LowPt_err_2016;
    DYscale_MTT_Med1Pt_err = DYscale_MTT_Med1Pt_err_2016;
    DYscale_MTT_Med2Pt_err = DYscale_MTT_Med2Pt_err_2016;
    DYscale_MTT_HighPt_err = DYscale_MTT_HighPt_err_2016;
    DYscale_MTT_vHighPt_err = DYscale_MTT_vHighPt_err_2016;
    
  } else if (year == 2017) {
    DYscale_MH_vLowPt = DYscale_MH_vLowPt_2017;
    DYscale_MH_LowPt = DYscale_MH_LowPt_2017;
    DYscale_MH_Med1Pt = DYscale_MH_Med1Pt_2017;
    DYscale_MH_Med2Pt = DYscale_MH_Med2Pt_2017;
    DYscale_MH_HighPt = DYscale_MH_HighPt_2017;
    DYscale_MH_vHighPt = DYscale_MH_vHighPt_2017;

    DYscale_MTT_vLowPt = DYscale_MTT_vLowPt_2017;
    DYscale_MTT_LowPt = DYscale_MTT_LowPt_2017;
    DYscale_MTT_Med1Pt = DYscale_MTT_Med1Pt_2017;
    DYscale_MTT_Med2Pt = DYscale_MTT_Med2Pt_2017;
    DYscale_MTT_HighPt = DYscale_MTT_HighPt_2017;
    DYscale_MTT_vHighPt = DYscale_MTT_vHighPt_2017;

    DYscale_MTT_vLowPt_err = DYscale_MTT_vLowPt_err_2017;
    DYscale_MTT_LowPt_err = DYscale_MTT_LowPt_err_2017;
    DYscale_MTT_Med1Pt_err = DYscale_MTT_Med1Pt_err_2017;
    DYscale_MTT_Med2Pt_err = DYscale_MTT_Med2Pt_err_2017;
    DYscale_MTT_HighPt_err = DYscale_MTT_HighPt_err_2017;
    DYscale_MTT_vHighPt_err = DYscale_MTT_vHighPt_err_2017;

  } else {
    DYscale_MH_vLowPt = DYscale_MH_vLowPt_2018;
    DYscale_MH_LowPt = DYscale_MH_LowPt_2018;
    DYscale_MH_Med1Pt = DYscale_MH_Med1Pt_2018;
    DYscale_MH_Med2Pt = DYscale_MH_Med2Pt_2018;
    DYscale_MH_HighPt = DYscale_MH_HighPt_2018;
    DYscale_MH_vHighPt = DYscale_MH_vHighPt_2018;

    DYscale_MTT_vLowPt = DYscale_MTT_vLowPt_2018;
    DYscale_MTT_LowPt = DYscale_MTT_LowPt_2018;
    DYscale_MTT_Med1Pt = DYscale_MTT_Med1Pt_2018;
    DYscale_MTT_Med2Pt = DYscale_MTT_Med2Pt_2018;
    DYscale_MTT_HighPt = DYscale_MTT_HighPt_2018;
    DYscale_MTT_vHighPt = DYscale_MTT_vHighPt_2018;

    DYscale_MTT_vLowPt_err = DYscale_MTT_vLowPt_err_2018;
    DYscale_MTT_LowPt_err = DYscale_MTT_LowPt_err_2018;
    DYscale_MTT_Med1Pt_err = DYscale_MTT_Med1Pt_err_2018;
    DYscale_MTT_Med2Pt_err = DYscale_MTT_Med2Pt_err_2018;
    DYscale_MTT_HighPt_err = DYscale_MTT_HighPt_err_2018;
    DYscale_MTT_vHighPt_err = DYscale_MTT_vHighPt_err_2018;
  }

};

// Destructor
DYscaling::~DYscaling() {}

std::vector<float> DYscaling::get_dy_scale(
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
    iRVec GenPart_pdgId) {
  TLorentzVector vgj;
  int nbs = 0;
  for (unsigned int igj = 0; igj < GenJet_pt.size(); igj++)
    {
      vgj.SetPtEtaPhiM(GenJet_pt.at(igj), GenJet_eta.at(igj), GenJet_phi.at(igj), GenJet_mass.at(igj));
      if (vgj.Pt() > 20 && fabs(vgj.Eta()) < 2.5)
        {
          int theFlav = GenJet_hadronFlavour.at(igj);
          if (abs(theFlav) == 5) nbs++;
        }
    }
  if (nbs > 2) nbs = 2;

  float DYscale_MH = 1.;
  float DYscale_MTT = 1.;
  float DYscale_MTT_up = 1.;
  float DYscale_MTT_down = 1.;

  int n_bJets = LHE_Nb;
  if (n_bJets > 2)
    n_bJets = 2;
  
  int idx = -1;
  for (unsigned int igen = 0; igen < GenPart_pt.size(); igen++) {
      bool isLast = (GenPart_statusFlags.at(igen) & (int) std::pow(2, 13)) != 0; // 13 = isLastCopy
      bool isPrompt = (GenPart_statusFlags.at(igen) & 1) != 0; //  0 = isPrompt
      if (GenPart_pdgId.at(igen) == 23 && isLast && isPrompt) { // Z0 + isLast + isPrompt
          idx = igen;
      }
  }
  
  if (idx >= 0) {
    TLorentzVector genZ;
    genZ.SetPtEtaPhiM(GenPart_pt.at(idx), GenPart_eta.at(idx), GenPart_phi.at(idx), GenPart_mass.at(idx));
    float genZ_pt = genZ.Pt();

    if (genZ_pt <= 10.) {
      DYscale_MH = DYscale_MH_vLowPt [n_bJets];
      DYscale_MTT = DYscale_MTT_vLowPt[n_bJets];
      DYscale_MTT_up = DYscale_MTT_vLowPt[n_bJets] + DYscale_MTT_vLowPt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_vLowPt[n_bJets] - DYscale_MTT_vLowPt_err[n_bJets];
    } else if (genZ_pt > 10. && genZ_pt <= 30.) {
      DYscale_MH = DYscale_MH_LowPt[n_bJets];
      DYscale_MTT = DYscale_MTT_LowPt[n_bJets];
      DYscale_MTT_up = DYscale_MTT_LowPt[n_bJets] + DYscale_MTT_LowPt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_LowPt[n_bJets] - DYscale_MTT_LowPt_err[n_bJets];
    } else if (genZ_pt > 30. && genZ_pt <= 50.) {
      DYscale_MH = DYscale_MH_Med1Pt[n_bJets];
      DYscale_MTT = DYscale_MTT_Med1Pt[n_bJets];
      DYscale_MTT_up  = DYscale_MTT_Med1Pt[n_bJets] + DYscale_MTT_Med1Pt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_Med1Pt[n_bJets] - DYscale_MTT_Med1Pt_err[n_bJets];
    } else if (genZ_pt > 50. && genZ_pt <= 100.) {
      DYscale_MH = DYscale_MH_Med2Pt[n_bJets];
      DYscale_MTT = DYscale_MTT_Med2Pt[n_bJets];
      DYscale_MTT_up = DYscale_MTT_Med2Pt[n_bJets] + DYscale_MTT_Med2Pt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_Med2Pt[n_bJets] - DYscale_MTT_Med2Pt_err[n_bJets];
    } else if (genZ_pt > 100. && genZ_pt <= 200.) {
      DYscale_MH = DYscale_MH_HighPt[n_bJets];
      DYscale_MTT = DYscale_MTT_HighPt[n_bJets];
      DYscale_MTT_up = DYscale_MTT_HighPt[n_bJets] + DYscale_MTT_HighPt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_HighPt[n_bJets] - DYscale_MTT_HighPt_err[n_bJets];
    } else {/* pT(Z)>=200. */ 
      DYscale_MH = DYscale_MH_vHighPt[n_bJets];
      DYscale_MTT = DYscale_MTT_vHighPt[n_bJets];
      DYscale_MTT_up = DYscale_MTT_vHighPt[n_bJets] + DYscale_MTT_vHighPt_err[n_bJets];
      DYscale_MTT_down = DYscale_MTT_vHighPt[n_bJets] - DYscale_MTT_vHighPt_err[n_bJets];
    }
  }
  return {DYscale_LL[nbs], DYscale_MM[nbs], DYscale_MH, DYscale_MTT, DYscale_MTT_up, DYscale_MTT_down};
}