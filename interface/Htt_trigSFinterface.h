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
#include <memory>


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

struct Htt_trigSFresult {
  Htt_trigSFresult() : passSingle(false), passCross(false), trigSF(1.), muUp(0.), muDown(0.), eleUp(0.), eleDown(0.), DM0Up(0.), DM0Down(0.), DM1Up(0.), DM1Down(0.), DM10Up(0.), DM10Down(0.), DM11Up(0.), DM11Down(0.) {}
  bool passSingle, passCross; // these only check if we are in a phase space where we look for single or cross trigger. It does not imply that the respective trigger has fired (if in overlap area, outside overlap area it does imply it)

  float trigSF; // nominal

  // uncertainties
  float muUp, muDown;
  float eleUp, eleDown;
  float DM0Up, DM0Down;
  float DM1Up, DM1Down;
  float DM10Up, DM10Down;
  float DM11Up, DM11Down;

  void setDMError(int dm, double error) {
    if (dm == 0) {
      DM0Up = trigSF + error;
      DM0Down = trigSF - error;
    } else if (dm == 1){
      DM1Up = trigSF + error;
      DM1Down = trigSF - error;
    } else if (dm == 10){
      DM10Up = trigSF + error;
      DM10Down = trigSF - error;
    } else if (dm == 11) {
      DM11Up = trigSF + error;
      DM11Down = trigSF - error;
    } else assert(false);
  }

  void setSystVariationsToNominal() {
    if (muUp==0.) muUp = trigSF;
    if (muDown==0.) muDown = trigSF;
    if (eleUp==0.) eleUp = trigSF;
    if (eleDown==0.) eleDown = trigSF;
    if (DM0Up==0.) DM0Up = trigSF;
    if (DM0Down==0.) DM0Down = trigSF;
    if (DM1Up==0.) DM1Up = trigSF;
    if (DM1Down==0.) DM1Down = trigSF;
    if (DM10Up==0.) DM10Up = trigSF;
    if (DM10Down==0.) DM10Down = trigSF;
    if (DM11Up==0.) DM11Up = trigSF;
    if (DM11Down==0.) DM11Down = trigSF;
 
  }
};


// Htt_trigSFInterface class
class Htt_trigSFinterface {
  
  public:
    Htt_trigSFinterface (
      int year, float mutau_pt_th1, float mutau_pt_th2, float etau_pt_th1, float etau_pt_th2,
      std::string eTrgSF_file, std::string eTrgSF_name, bool eTrgSF_bool, 
      std::string eTauTrgSF_file, std::string eTauTrgSF_name, bool eTauTrgSF_bool,
      std::string muTrgSF_file, std::string muTrgSF_name, bool muTrgSF_bool,
      std::string muTauTrgSF_file, std::string muTauTrgSF_name, bool muTauTrgSF_bool,
      std::string tauTrgSF_ditau_file, std::string tauTrgSF_mutau_file,
      std::string tauTrgSF_etau_file, std::string tauTrgSF_vbf_file, std::string jetTrgSF_vbf_file);    

    Htt_trigSFresult get_scale_factors(int pairType, int isVBFtrigger,
      int dau1_decayMode, float dau1_pt, float dau1_eta,
      int dau2_decayMode, float dau2_pt, float dau2_eta);
  
    //private:
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


class MET_trigSF_interface {
public:
  struct trigSF_result {
    float SF;
    float SF_statup;
    float SF_statdown;
  };

  MET_trigSF_interface(std::string SF_file);

  trigSF_result getSF(float MET_pt);

private:
  std::unique_ptr<TH1F> trigSF_h;
};


#endif // Htt_trigSFinterface
