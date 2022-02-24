#ifndef PUjetID_SFInterface_h
#define PUjetID_SFInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class PUjetID_SFInterface                                                                                    //
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
#include <TH2.h>
#include <TFile.h>

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;
typedef ROOT::VecOps::RVec<int> iRVec;


// PUjetID_SFInterface class
class PUjetID_SFinterface {

  public:
    PUjetID_SFinterface (int year, std::string folder_path);    
    ~PUjetID_SFinterface ();
    std::vector <double> get_pu_weights(
        fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass, iRVec Jet_jetId, iRVec Jet_puId,
        fRVec GenJet_pt, fRVec GenJet_eta, fRVec GenJet_phi, fRVec GenJet_mass,
        double dau1_pt, double dau1_eta, double dau1_phi, double dau1_mass,
        double dau2_pt, double dau2_eta, double dau2_phi, double dau2_mass);

  private:
    int year_;
    TH2F* h_eff_; 
    TH2F* h_eff_sf_;
    TH2F* h_eff_sf_err_;
    TH2F* h_mistag_;
    TH2F* h_mistag_sf_;
    TH2F* h_mistag_sf_err_;

    double getContentHisto2D(TH2F* histo, double x, double y) {
      auto nbinsx = histo->GetNbinsX();
      auto nbinsy = histo->GetNbinsY();
      auto ibinx = histo->GetXaxis()->FindBin(x);
      auto ibiny = histo->GetYaxis()->FindBin(y);

      if (ibinx == 0)
        ibinx = 1;
      else if (ibinx > nbinsx)
        ibinx = nbinsx;

      if (ibiny == 0)
        ibiny = 1;
      else if (ibiny > nbinsy)
        ibiny = nbinsy;

      return histo->GetBinContent(ibinx, ibiny);
    }

    std::vector <double> get_eff_sf_and_error(bool isReal, double pt, double eta) {
        if (pt < 20.)
          pt = 20.;
        else if (pt > 50.)
          pt = 50.;

        if (isReal)
          return {getContentHisto2D(h_eff_, (double) pt, (double) eta),
            getContentHisto2D(h_eff_sf_, (double) pt, (double) eta),
            getContentHisto2D(h_eff_sf_err_, (double) pt, (double) eta)};
        else
          return {getContentHisto2D(h_mistag_, (double) pt, (double) eta),
            getContentHisto2D(h_mistag_sf_, (double) pt, (double) eta),
            getContentHisto2D(h_mistag_sf_err_, (double) pt, (double) eta)};
    }

};

#endif // PUjetID_SFinterface
