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
#include "Tools/Tools/interface/FeatureProvider.h"

class MulticlassInterface {
public:
  MulticlassInterface(int year, const std::vector<std::pair<std::string, std::string>> &modelSpecs);
  ~MulticlassInterface();

  size_t getNumberOfModels() const;

  const std::vector<std::string>& getNodeNames(size_t modelIndex) const;

  void clearInputs();

  void setInput(size_t modelIndex, const std::string& featureName, float value, bool optional=true);

  void setInputs(size_t modelIndex, const std::vector<std::pair<std::string, float>>& inputs, bool optional=true);

  std::vector<std::pair<std::string, float>> predict(hmc::EventId eventId, size_t modelIndex);

  std::vector<std::pair<std::string, float>> predict(hmc::EventId eventId, size_t modelIndex,
    const std::vector<float>& inputs);

  void extendTree(TTree* tree, const std::string& branchPostfix = "");

private:
  int year_;
  std::vector<std::pair<std::string, std::string>> modelSpecs_;
  std::vector<hmc::Model *> models_;
  std::vector<std::vector<std::string>> nodeNames_;

  void checkModelIndex_(size_t modelIndex) const;
};
