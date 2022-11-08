// --------------------------------------------------------------------------------------------------------------
// //
//                                                                                                                //
//   class MulticlassInterface //
//                                                                                                                //
//   Class to compute DNN output from the Multiclass approach //
//                                                                                                                //
//   Author: Jaime Leon Holgado (CIEMAT), Marcel Rieger (CERN) // Date  : June 2020 //
//                                                                                                                //
// --------------------------------------------------------------------------------------------------------------
// //

#include "Tools/Tools/interface/FeatureProvider.h"

hmc::EventId FeatureProvider::getEventId() const {
  return hmc::EventId(ulong64Inputs_.at("EventNumber"));
}

void FeatureProvider::add(const std::string &featureName) {
  features_.emplace(featureName, 0.);
}

bool FeatureProvider::has(const std::string &featureName) {
  return features_.count(featureName) == 1;
}

float FeatureProvider::get(const std::string &featureName) const {
  const auto &it = features_.find(featureName);
  if (it == features_.end()) {
    throw std::runtime_error("FeatureProvider: unknown feature '" + featureName + "'");
  }
  return it->second;
}

FeatureProvider::FeatureProvider(int year, TTree *tree) : year_(year) {
  // define names of variables to read
  std::vector<std::string> boolNames = {};
  std::vector<std::string> intNames = {"pairType"};
  std::vector<std::string> ulong64Names = {"EventNumber"};
  std::vector<std::string> floatNames = {
    "bjet1_pt", "bjet1_eta", "bjet1_phi", "bjet1_e", "bjet1_bID_deepFlavor", "bjet1_CvsB",
    "bjet1_CvsL", "bjet1_HHbtag",
    "bjet2_pt", "bjet2_eta", "bjet2_phi", "bjet2_e", "bjet2_bID_deepFlavor", "bjet2_CvsB",
    "bjet2_CvsL", "bjet2_HHbtag",
    "addJetCentr1_pt", "addJetCentr1_eta", "addJetCentr1_phi", "addJetCentr1_e", "addJetCentr1_btag_deepFlavor", "addJetCentr1_HHbtag",
    "addJetCentr2_pt", "addJetCentr2_eta", "addJetCentr2_phi", "addJetCentr2_e", "addJetCentr2_btag_deepFlavor", "addJetCentr2_HHbtag",
    "addJetCentr3_pt", "addJetCentr3_eta", "addJetCentr3_phi", "addJetCentr3_e", "addJetCentr3_btag_deepFlavor", "addJetCentr3_HHbtag",
    "addJetForw1_pt", "addJetForw1_eta", "addJetForw1_phi", "addJetForw1_e",
    "addJetForw2_pt", "addJetForw2_eta", "addJetForw2_phi", "addJetForw2_e",
    "VBFjet1_pt", "VBFjet1_eta", "VBFjet1_phi", "VBFjet1_e", "VBFjet1_btag_deepFlavor",
    "VBFjet1_CvsB", "VBFjet1_CvsL", "VBFjet1_HHbtag",
    "VBFjet2_pt", "VBFjet2_eta", "VBFjet2_phi", "VBFjet2_e", "VBFjet2_btag_deepFlavor",
    "VBFjet2_CvsB", "VBFjet2_CvsL", "VBFjet2_HHbtag",
    "dau1_pt", "dau1_eta", "dau1_phi", "dau1_e", "dau2_pt", "dau2_eta", "dau2_phi", "dau2_e",
    "met_et", "met_phi",
    "bH_pt", "bH_eta", "bH_phi", "bH_e",
    "tauH_SVFIT_pt", "tauH_SVFIT_eta", "tauH_SVFIT_phi", "tauH_SVFIT_mass"
  };

  // register them in input maps and set branch addresses
  for (const auto &name : boolNames) {
    boolInputs_.emplace(name, 0.);
    tree->SetBranchAddress(name.c_str(), &boolInputs_.at(name));
  }
  for (const auto &name : intNames) {
    intInputs_.emplace(name, 0.);
    tree->SetBranchAddress(name.c_str(), &intInputs_.at(name));
  }
  for (const auto &name : ulong64Names) {
    ulong64Inputs_.emplace(name, 0);
    tree->SetBranchAddress(name.c_str(), &ulong64Inputs_.at(name));
  }
  for (const auto &name : floatNames) {
    floatInputs_.emplace(name, 0.);
    tree->SetBranchAddress(name.c_str(), &floatInputs_.at(name));
  }
}

void FeatureProvider::calculate() {
  // check if objects are set
  bool b1Set = floatInputs_.at("bjet1_pt") > 0;
  bool b2Set = floatInputs_.at("bjet2_pt") > 0;
  bool ctjet1Set = floatInputs_.at("addJetCentr1_pt") > 0;
  bool ctjet2Set = floatInputs_.at("addJetCentr2_pt") > 0;
  bool ctjet3Set = floatInputs_.at("addJetCentr3_pt") > 0;
  bool fwjet1Set = floatInputs_.at("addJetForw1_pt") > 0;
  bool fwjet2Set = floatInputs_.at("addJetForw2_pt") > 0;
  bool vbfj1Set = floatInputs_.at("VBFjet1_pt") > 0;
  bool vbfj2Set = floatInputs_.at("VBFjet2_pt") > 0;
  bool lep1Set = floatInputs_.at("dau1_pt") > 0;
  bool lep2Set = floatInputs_.at("dau2_pt") > 0;
  bool bHSet = floatInputs_.at("bH_pt") > 0;
  bool tauHSet = floatInputs_.at("tauH_SVFIT_pt") > 0;

  // define vectors for objects
  TLorentzVector tauH;
  tauH.SetPtEtaPhiM(floatInputs_.at("tauH_SVFIT_pt"), floatInputs_.at("tauH_SVFIT_eta"),
		    floatInputs_.at("tauH_SVFIT_phi"), floatInputs_.at("tauH_SVFIT_mass"));

  // loop through features and set values
  for (auto &it : features_) {
    if (it.first == "is_mutau") {
      it.second = float(intInputs_.at("pairType") == 0);
    } else if (it.first == "is_etau") {
      it.second = float(intInputs_.at("pairType") == 1);
    } else if (it.first == "is_tautau") {
      it.second = float(intInputs_.at("pairType") == 2);
    } else if (it.first == "is_2016") {
      it.second = (year_ == 2016);
    } else if (it.first == "is_2017") {
      it.second = (year_ == 2017);
    } else if (it.first == "is_2018") {
      it.second = (year_ == 2018);
    } else if (it.first == "bjet1_pt") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_pt"));
    } else if (it.first == "bjet1_eta") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_eta"));
    } else if (it.first == "bjet1_phi") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_phi"));
    } else if (it.first == "bjet1_e") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_e"));
    } else if (it.first == "bjet1_deepflavor_b") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_bID_deepFlavor"));
    } else if (it.first == "bjet1_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_CvsB"));
    } else if (it.first == "bjet1_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_CvsL"));
    } else if (it.first == "bjet1_hhbtag") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("bjet1_HHbtag"));
    } else if (it.first == "bjet2_pt") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_pt"));
    } else if (it.first == "bjet2_eta") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_eta"));
    } else if (it.first == "bjet2_phi") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_phi"));
    } else if (it.first == "bjet2_e") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_e"));
    } else if (it.first == "bjet2_deepflavor_b") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_bID_deepFlavor"));
    } else if (it.first == "bjet2_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_CvsB"));
    } else if (it.first == "bjet2_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_CvsL"));
    } else if (it.first == "bjet2_hhbtag") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("bjet2_HHbtag"));
    } else if (it.first == "ctjet1_pt") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_pt"));
    } else if (it.first == "ctjet1_eta") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_eta"));
    } else if (it.first == "ctjet1_phi") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_phi"));
    } else if (it.first == "ctjet1_e") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_e"));
    } else if (it.first == "ctjet1_deepflavor_b") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_btag_deepFlavor"));
    } else if (it.first == "ctjet1_hhbtag") {
      it.second = CHECK_EMPTY(ctjet1Set, floatInputs_.at("addJetCentr1_HHbtag"));
    } else if (it.first == "ctjet2_pt") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_pt"));
    } else if (it.first == "ctjet2_eta") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_eta"));
    } else if (it.first == "ctjet2_phi") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_phi"));
    } else if (it.first == "ctjet2_e") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_e"));
    } else if (it.first == "ctjet2_deepflavor_b") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_btag_deepFlavor"));
    } else if (it.first == "ctjet2_hhbtag") {
      it.second = CHECK_EMPTY(ctjet2Set, floatInputs_.at("addJetCentr2_HHbtag"));
    } else if (it.first == "ctjet3_pt") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_pt"));
    } else if (it.first == "ctjet3_eta") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_eta"));
    } else if (it.first == "ctjet3_phi") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_phi"));
    } else if (it.first == "ctjet3_e") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_e"));
    } else if (it.first == "ctjet3_deepflavor_b") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_btag_deepFlavor"));
    } else if (it.first == "ctjet3_hhbtag") {
      it.second = CHECK_EMPTY(ctjet3Set, floatInputs_.at("addJetCentr3_HHbtag"));
    } else if (it.first == "fwjet1_pt") {
      it.second = CHECK_EMPTY(fwjet1Set, floatInputs_.at("addJetForw1_pt"));
    } else if (it.first == "fwjet1_eta") {
      it.second = CHECK_EMPTY(fwjet1Set, floatInputs_.at("addJetForw1_eta"));
    } else if (it.first == "fwjet1_phi") {
      it.second = CHECK_EMPTY(fwjet1Set, floatInputs_.at("addJetForw1_phi"));
    } else if (it.first == "fwjet1_e") {
      it.second = CHECK_EMPTY(fwjet1Set, floatInputs_.at("addJetForw1_e"));
    } else if (it.first == "fwjet2_pt") {
      it.second = CHECK_EMPTY(fwjet2Set, floatInputs_.at("addJetForw2_pt"));
    } else if (it.first == "fwjet2_eta") {
      it.second = CHECK_EMPTY(fwjet2Set, floatInputs_.at("addJetForw2_eta"));
    } else if (it.first == "fwjet2_phi") {
      it.second = CHECK_EMPTY(fwjet2Set, floatInputs_.at("addJetForw2_phi"));
    } else if (it.first == "fwjet2_e") {
      it.second = CHECK_EMPTY(fwjet2Set, floatInputs_.at("addJetForw2_e"));
    } else if (it.first == "vbfjet1_pt") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_pt"));
    } else if (it.first == "vbfjet1_eta") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_eta"));
    } else if (it.first == "vbfjet1_phi") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_phi"));
    } else if (it.first == "vbfjet1_e") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_e"));
    } else if (it.first == "vbfjet1_deepflavor_b") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_btag_deepFlavor"));
    } else if (it.first == "vbfjet1_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_CvsB"));
    } else if (it.first == "vbfjet1_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_CvsL"));
    } else if (it.first == "vbfjet1_hhbtag") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBFjet1_HHbtag"));
    } else if (it.first == "vbfjet2_pt") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_pt"));
    } else if (it.first == "vbfjet2_eta") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_eta"));
    } else if (it.first == "vbfjet2_phi") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_phi"));
    } else if (it.first == "vbfjet2_e") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_e"));
    } else if (it.first == "vbfjet2_deepflavor_b") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_btag_deepFlavor"));
    } else if (it.first == "vbfjet2_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_CvsB"));
    } else if (it.first == "vbfjet2_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_CvsL"));
    } else if (it.first == "vbfjet2_hhbtag") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBFjet2_HHbtag"));
    } else if (it.first == "lep1_pt") {
      it.second = CHECK_EMPTY(lep1Set, floatInputs_.at("dau1_pt"));
    } else if (it.first == "lep1_eta") {
      it.second = CHECK_EMPTY(lep1Set, floatInputs_.at("dau1_eta"));
    } else if (it.first == "lep1_phi") {
      it.second = CHECK_EMPTY(lep1Set, floatInputs_.at("dau1_phi"));
    } else if (it.first == "lep1_e") {
      it.second = CHECK_EMPTY(lep1Set, floatInputs_.at("dau1_e"));
    } else if (it.first == "lep2_pt") {
      it.second = CHECK_EMPTY(lep2Set, floatInputs_.at("dau2_pt"));
    } else if (it.first == "lep2_eta") {
      it.second = CHECK_EMPTY(lep2Set, floatInputs_.at("dau2_eta"));
    } else if (it.first == "lep2_phi") {
      it.second = CHECK_EMPTY(lep2Set, floatInputs_.at("dau2_phi"));
    } else if (it.first == "lep2_e") {
      it.second = CHECK_EMPTY(lep2Set, floatInputs_.at("dau2_e"));
    } else if (it.first == "met_pt") {
      it.second = floatInputs_.at("met_et");
    } else if (it.first == "met_phi") {
      it.second = floatInputs_.at("met_phi");
    } else if (it.first == "bh_pt") {
      it.second = CHECK_EMPTY(bHSet, floatInputs_.at("bH_pt"));
    } else if (it.first == "bh_eta") {
      it.second = CHECK_EMPTY(bHSet, floatInputs_.at("bH_eta"));
    } else if (it.first == "bh_phi") {
      it.second = CHECK_EMPTY(bHSet, floatInputs_.at("bH_phi"));
    } else if (it.first == "bh_e") {
      it.second = CHECK_EMPTY(bHSet, floatInputs_.at("bH_e"));
    } else if (it.first == "tauh_sv_pt") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Pt());
    } else if (it.first == "tauh_sv_eta") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Eta());
    } else if (it.first == "tauh_sv_phi") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Phi());
    } else if (it.first == "tauh_sv_e") {  // used in v3b
      it.second = CHECK_EMPTY(tauHSet, tauH.E());
    } else if (it.first == "tauh_sv_ez") {  // used in v3 and before
      it.second = CHECK_EMPTY(tauHSet, pow(pow(tauH.Pz(), 2) + pow(tauH.M(), 2), 0.5));
    } else {
      throw std::runtime_error("MulticlassInference: unhandled feature '" + it.first + "'");
    }
  }
}
