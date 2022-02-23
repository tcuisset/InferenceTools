#ifndef Htt_trigSFInterface_h
#define Htt_trigSFInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class Htt_trigSFInterface                                                                                    //
//                                                                                                                //
//   Class to compute Htt trigger scale factors.                                                                  //
//                                                                                                                //
//   Author: Jaime León Holgado                                                                                   //
//   Date  : Feb 2022                                                                                             //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


// ROOT libraries
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>
#include <Math/VectorUtil.h>
#include <TH3.h>
#include <TFile.h>

// CMSSW
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/SFProvider.h"

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;


// Htt_trigSFInterface class
class Htt_trigSFinterface {
  
  public:
    Htt_trigSFinterface (
      int year, float mutau_pt_th1, float mutau_pt_th2, float etau_pt_th1, float etau_pt_th2,
      std::string eTrgSF_file, std::string eTauTrgSF_file, std::string muTrgSF_file,
      std::string muTauTrgSF_file, std::string tauTrgSF_ditau_file, std::string tauTrgSF_mutau_file,
      std::string tauTrgSF_etau_file, std::string tauTrgSF_vbf_file, std::string jetTrgSF_vbf_file);    
    ~Htt_trigSFinterface ();

    std::vector<double> get_scale_factors(int pairType, int isVBFtrigger,
      int dau1_decayMode, float dau1_pt, float dau1_eta,
      int dau2_decayMode, float dau2_pt, float dau2_eta,
      float vbfjet1_pt, float vbfjet1_eta, float vbfjet1_phi, float vbfjet1_mass,
      float vbfjet2_pt, float vbfjet2_eta, float vbfjet2_phi, float vbfjet2_mass);
  
    private:
    int year_;
    float mutau_pt_th1_;
    float mutau_pt_th2_;
    float etau_pt_th1_;
    float etau_pt_th2_;
    ScaleFactor eTrgSF = ScaleFactor();
    ScaleFactor eTauTrgSF = ScaleFactor();
    ScaleFactor muTrgSF = ScaleFactor();
    ScaleFactor muTauTrgSF;
    tau_trigger::SFProvider tauTrgSF_ditau;
    tau_trigger::SFProvider tauTrgSF_mutau;
    tau_trigger::SFProvider tauTrgSF_etau;
    tau_trigger::SFProvider tauTrgSF_vbf;
    TH3F* jetTrgSF_vbf;

    double getContentHisto3D(TH3F* TH3, double x, double y, double z, bool unc) {
      auto nbinsx = TH3->GetNbinsX();
      auto nbinsy = TH3->GetNbinsY();
      auto nbinsz = TH3->GetNbinsZ();

      auto ibinx = TH3->GetXaxis()->FindBin(x);
      auto ibiny = TH3->GetYaxis()->FindBin(y);
      auto ibinz = TH3->GetZaxis()->FindBin(z);

      if (ibinx == 0)
        ibinx = 1;
      else if (ibinx > nbinsx)
        ibinx = nbinsx;

      if (ibiny == 0)
        ibiny = 1;
      else if (ibiny > nbinsy)
        ibiny = nbinsy;

      if (ibinz == 0)
        ibinz = 1;
      else if (ibinz > nbinsz)
        ibinz = nbinsz;

      if (!unc)
          return TH3->GetBinContent(ibinx, ibiny, ibinz);
      else
          return TH3->GetBinError(ibinx, ibiny, ibinz);

    }
};

#endif // Htt_trigSFinterface
