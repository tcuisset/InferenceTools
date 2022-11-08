#include "Tools/Tools/interface/MulticlassInterface.h"

MulticlassInterface::MulticlassInterface(int year, const std::vector<std::pair<std::string, std::string>> &modelSpecs)
  : year_(year), modelSpecs_(modelSpecs) {
  // load models
  for (const auto &modelSpec : modelSpecs_) {
    const std::string &version = modelSpec.first;
    const std::string &tag = modelSpec.second;

    // load the model
    models_.push_back(hmc::loadModel(year_, version, tag));

    // save total node names
    nodeNames_.push_back(models_.back()->getAllNodeNames());
  }
}

MulticlassInterface::~MulticlassInterface() {}

size_t MulticlassInterface::getNumberOfModels() const {
  return models_.size();
}

void MulticlassInterface::checkModelIndex_(size_t modelIndex) const {
  if (modelIndex >= models_.size()) {
    throw std::runtime_error("modelIndex " + std::to_string(modelIndex)
           + "too large for number of models " + std::to_string(models_.size()));
  }
}

const std::vector<std::string>& MulticlassInterface::getNodeNames(size_t modelIndex) const {
  checkModelIndex_(modelIndex);
  return nodeNames_[modelIndex];
}

void MulticlassInterface::clearInputs() {
  for (auto& model : models_) {
    model->input.clear();
  }
}

void MulticlassInterface::setInput(size_t modelIndex, const std::string& featureName, float value, bool optional) {
  checkModelIndex_(modelIndex);
  if (!optional || models_[modelIndex]->hasFeatureSpec(featureName)) {
    models_[modelIndex]->input.setValue(featureName, value);
  }
}

void MulticlassInterface::setInputs(size_t modelIndex, const std::vector<std::pair<std::string, float>>& inputs, bool optional) {
  checkModelIndex_(modelIndex);
  for (const auto& pair : inputs) {
    if (!optional || models_[modelIndex]->hasFeatureSpec(pair.first)) {
      models_[modelIndex]->input.setValue(pair.first, pair.second);
    }
  }
}

std::vector<std::pair<std::string, float>> MulticlassInterface::predict(hmc::EventId eventId, size_t modelIndex) {
  checkModelIndex_(modelIndex);

  // check if the input spec is complete, i.e., if all values were set
  auto& inputSpec = models_[modelIndex]->input;
  if (!inputSpec.complete()) {
    throw std::runtime_error("model input spec is incomplete, only "
           + std::to_string(inputSpec.getNumberOfSetFeatures()) + " out of "
           + std::to_string(inputSpec.getNumberOfFeatures()) + " features set");
  }

  // run the prediction with input values
  return predict(eventId, modelIndex, inputSpec.getValues());
}

std::vector<std::pair<std::string, float>> MulticlassInterface::predict(hmc::EventId eventId, size_t modelIndex,
               const std::vector<float>& inputs) {
  checkModelIndex_(modelIndex);

  // run the model
  hmc::Model*& model = models_[modelIndex];
  std::vector<float> outputs;
  model->run(eventId, inputs, outputs);

  // extend outputs by merged nodes
  for (const auto& pair : model->getNodeMerging()) {
    float merged = 0.;
    for (const int& i : pair.second) {
      merged += outputs[i];
    }
    outputs.push_back(merged);
  }

  // the number of outputs after adding merged nodes must match the number of node names
  const std::vector<std::string>& nodeNames = nodeNames_[modelIndex];
  if (outputs.size() != nodeNames.size()) {
    throw std::runtime_error("number of outputs " + std::to_string(outputs.size())
           + " does not match number of all output node names " + std::to_string(nodeNames.size()));
  }

  // create the output structure
  std::vector<std::pair<std::string, float>> outputPairs;
  for (size_t i = 0; i < outputs.size(); i++) {
    outputPairs.emplace_back(nodeNames[i], outputs[i]);
  }

  return outputPairs;
}

void MulticlassInterface::extendTree(TTree* tree, const std::string& branchPostfix) {
  // create the feature provider
  FeatureProvider features(year_, tree);

  // loop through models, register features and create branches
  std::vector<std::vector<TBranch*>> branches;
  for (size_t i = 0; i < models_.size(); i++) {
    hmc::Model *&model = models_[i];
    const std::string &version = modelSpecs_[i].first;
    const std::string &tag = modelSpecs_[i].second;

    // register required features with the feature provider
    for (const auto &featureName : model->getFeatureNames()) {
      features.add(featureName);
    }

    // define branches per output node
    std::vector<TBranch*> modelBranches;
    for (const auto &nodeName : model->getAllNodeNames()) {
      auto branchName = "mdnn__" + version + "__" + tag + "__" + nodeName + branchPostfix;
      TBranch* b = tree->Branch(branchName.c_str(), model->output.getOutputAddress(nodeName),
        (branchName + "/F").c_str());
      modelBranches.push_back(b);
    }
    branches.push_back(modelBranches);
  }

  // start iterating
  int nEntries = tree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    // load the entry and calculate features
    tree->GetEntry(i);
    features.calculate();

    // fill features of all models and run them
    for (size_t i = 0; i < models_.size(); i++) {
      hmc::Model*& model = models_[i];
      model->input.clear();
      for (const auto &featureName : model->getFeatureNames()) {
        model->input.setValue(featureName, features.get(featureName));
      }

      model->run(features.getEventId());

      // fill branches
      for (TBranch *&b : branches[i]) {
        b->Fill();
      }
    }
  }
}
