#include "Tools/Tools/interface/DYreweighting.h"

// Constructors

DYreweighting::DYreweighting (int year) {
  if (year == 2016) {
    stitchWeights = stitchWeights_2016;
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
    stitchWeights = stitchWeights_2017;
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
    stitchWeights = stitchWeights_2018;
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
DYreweighting::~DYreweighting() {}

float DYreweighting::get_stitching_weight(int LHE_Nb, int LHE_Njets, float LHE_HT) {
  int njets = LHE_Njets;
  int nb = LHE_Nb;
  // these protections should be useless
  if (njets < 0) njets = 0;
  if (njets > 4) njets = 4;
  if (nb < 0)    nb = 0;
  if (nb > 4)    nb = 4;

  float ht = LHE_HT;
  int nht = 0;
  if      (ht  < 0                ) nht = 0;
  else if (ht >= 0    && ht < 70  ) nht = 0;
  else if (ht >= 70   && ht < 100 ) nht = 1;
  else if (ht >= 100  && ht < 200 ) nht = 2;
  else if (ht >= 200  && ht < 400 ) nht = 3;
  else if (ht >= 400  && ht < 600 ) nht = 4;
  else if (ht >= 600  && ht < 800 ) nht = 5;
  else if (ht >= 800  && ht < 1200) nht = 6;
  else if (ht >= 1200 && ht < 2500) nht = 7;
  else  /* ht >= 2500 */            nht = 8;

  return stitchWeights[njets][nb][nht];
}


std::vector<float> DYreweighting::get_dy_scale(std::vector<float> GenJet_pt, std::vector<float> GenJet_eta, std::vector<float> GenJet_phi, std::vector<float> GenJet_mass, std::vector<int> GenJet_hadronFlavour, int LHE_Nb,
    std::vector<float> GenPart_pt, std::vector<float> GenPart_eta, std::vector<float> GenPart_phi, std::vector<float> GenPart_mass, std::vector<int> GenPart_statusFlags, std::vector<int> GenPart_pdgId) {
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
    genZ.SetPtEtaPhiM(GenJet_pt.at(idx), GenJet_eta.at(idx), GenJet_phi.at(idx), GenJet_mass.at(idx));
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