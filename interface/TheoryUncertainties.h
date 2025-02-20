#ifndef TheoryUncertainties_h
#define TheoryUncertainties_h

// Author : Theo Cuisset (LLR)
// Implementation taken from columnflow

#include <array>
#include <ROOT/RVec.hxx>


// Returns nominal, up, down
std::array<float, 3> computePDFWeights(ROOT::VecOps::RVec<float> LHEPdfWeight);


// Returns nominal, up, down
std::array<float, 3> computeScaleUncertainties(ROOT::VecOps::RVec<float> LHEScaleWeight);

#endif