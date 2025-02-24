#include "Tools/Tools/interface/TheoryUncertainties.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

// Returns nominal, up, down
std::array<float, 3> computePDFWeights(ROOT::VecOps::RVec<float> LHEPdfWeight) {
    // Implementation from https://gitlab.cern.ch/cms-analysis/general/columnflow/-/blob/master/columnflow/production/cms/pdf.py
    if (LHEPdfWeight.size() != 101 && LHEPdfWeight.size() != 103) // 
        throw std::runtime_error(std::string("Invalid LHEPdfWeight with size ") + std::to_string(LHEPdfWeight.size()));
    float nominal = LHEPdfWeight[0];
    if (std::abs(nominal-1)>0.3)
        std::cout << "WARNING : Nominal LHEPdfWeight far from 1, value " << nominal << "\n";
        //throw std::runtime_error(std::string("Nominal LHEPdfWeight far from 1, value ") + std::to_string(nominal));
    // From https://arxiv.org/pdf/1510.03865v2 
    // suitable for MC sets like PDF4LHC15_mc
    auto sorted_values = ROOT::VecOps::Sort(ROOT::VecOps::Take(LHEPdfWeight/nominal, -(LHEPdfWeight.size()-1))); // Remove value 0 (nominal) then sort the weights
    float std_dev = (sorted_values[83]-sorted_values[15])/2.f;
    if (std_dev > 0.5f)
        std::cout << "WARNING : large pdf variation, stddev = " << std_dev << "\n";
    return {{1., 1.f+std_dev, 1.f-std_dev}};
}


// Returns nominal, up, down
std::array<float, 3> computeScaleUncertainties(ROOT::VecOps::RVec<float> LHEScaleWeight) {
    // Implementation from https://gitlab.cern.ch/cms-analysis/general/columnflow/-/blob/master/columnflow/production/cms/scale.py
    float murf_nominal;
    ROOT::RVec<ROOT::VecOps::RVec<float>::size_type> indices_to_take;
    if (LHEScaleWeight.size() == 9) {
        /*
        LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; 
        [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
        skipping : 4 (nominal), 2 (down/up crossed variation), 6 (up/down crossed variation) (crossed variations recommended to be ignored, slide 13 : https://indico.cern.ch/event/938672/contributions/3943718/attachments/2073936/3482265/MC_ContactReport_v3.pdf)
        */
        murf_nominal = LHEScaleWeight[4];
        indices_to_take = {{0, 1, 3, 5, 7, 8}};

    } else if (LHEScaleWeight.size() == 8) {
        /* LHE scale variation weights (w_var / w_nominal); [0] is MUF="0.5" MUR="0.5"; [1] is MUF="1.0" MUR="0.5"; [2] is MUF="2.0" MUR="0.5"; [3] is MUF="0.5" MUR="1.0"; [4] is MUF="2.0" MUR="1.0"; [5] is MUF="0.5" MUR="2.0"; [6] is MUF="1.0" MUR="2.0"; [7] is MUF="2.0" MUR="2.0"' */
        murf_nominal = 1;
        indices_to_take = {{0, 1, 3, 4, 6, 7}};
    } else
        throw std::runtime_error(std::string("Wrong nb of LHEScaleWeight, value ") + std::to_string(LHEScaleWeight.size()));

    if (std::abs(murf_nominal-1)>0.5)
        std::cout << "WARNING : Nominal mu renorm far from 1, value " << murf_nominal << "\n";
        //throw std::runtime_error(std::string("Nominal mu renorm far from 1, value ") + std::to_string(murf_nominal));
    
    auto murf_weights = ROOT::VecOps::Take(LHEScaleWeight, indices_to_take)/murf_nominal; 
    return {{1., ROOT::VecOps::Max(murf_weights), ROOT::VecOps::Min(murf_weights)}};
}
