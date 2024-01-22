// Script for testing the scale factor loading code
// USER_CXXFLAGS="-g -O0" scram b -j10 Tools/Tools HTT-utilities/LepEffInterface
#define private public
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#undef private
#include "Tools/Tools/interface/Htt_trigSFinterface.h"

std::string file = "/home/llr/cms/cuisset/bbtautau/ZHbbtautau/frameworkRepo/nanoaod_base_analysis/data/cmssw/CMSSW_12_3_0_pre6/src/HTT-utilities/trigSFs_UL_eleMu/sf_el_2018_HLTEle32.root"

float pt = 40
// 
ScaleFactor eTrgSF;
// UL2018 eTrgSF
eTrgSF.init_EG_ScaleFactor(file, true)

// eTrgSF.etaBinsH->GetXaxis()->GetBinLabel(3)
// eTrgSF.etaBinsH->GetXaxis()->GetBinLowEdge(3)
// eTrgSF.etaBinsH->GetXaxis()->GetBinUpEdge(3)


// std::string etaLabel = eTrgSF.FindEtaLabel(0.5, "data") // Eta0p0to0p8 : correct

// auto& graph = *(eTrgSF.eff_data["Eta0p0to0p8"])
// graph.Draw()

// int ptbin = eTrgSF.FindPtBin(eTrgSF.eff_data, etaLabel, pt); 


eTrgSF.get_EfficiencyData(pt, 0.5)




/// Testing muTauTrg which has absolute eta
#define private public
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#undef private
#include "Tools/Tools/interface/Htt_trigSFinterface.h"

std::string file = "/home/llr/cms/cuisset/bbtautau/ZHbbtautau/frameworkRepo/nanoaod_base_analysis/data/cmssw/CMSSW_12_3_0_pre6/src/HTT-utilities/trigSFs_UL_eleMu/sf_mu_2018_HLTMu20Tau27.root"

ScaleFactor muTauTrgSF;
muTauTrgSF.init_EG_ScaleFactor(file, true)

muTauTrgSF.etaIsAbsolute

muTauTrgSF.get_EfficiencyData(40, -1.3)
