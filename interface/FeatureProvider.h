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

// Standard libraries
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <exception>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/format.hpp>

// Multiclass DNN
#include "MulticlassInference/MulticlassInference/interface/hmc.h"

#define CHECK_EMPTY(COND, EXPR) ((COND) ? (EXPR) : (hmc::features::EMPTY))

class FeatureProvider {
public:
  FeatureProvider(int year, TTree *inTree);

  void calculate();  
  hmc::EventId getEventId() const;
  void add(const std::string &featureName);
  bool has(const std::string &featureName);
  float get(const std::string &featureName) const;
private:
  int year_;
  std::map<std::string, bool> boolInputs_;
  std::map<std::string, int> intInputs_;
  std::map<std::string, ULong64_t> ulong64Inputs_;
  std::map<std::string, float> floatInputs_;
  std::map<std::string, float> features_;

};
