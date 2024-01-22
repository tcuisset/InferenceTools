// Script for testing the scale factor loading code
// USER_CXXFLAGS="-g -O0" scram b -j10 Tools/Tools HTT-utilities/LepEffInterface
#define private public
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#undef private
#include "Tools/Tools/interface/Htt_trigSFinterface.h"

std::string file = "/home/llr/cms/cuisset/bbtautau/ZHbbtautau/frameworkRepo/nanoaod_base_analysis/data/cmssw/CMSSW_12_3_0_pre6/src/HTT-utilities/LepEffInterface/data/Muon/Run2018/Muon_Run2018_IsoMu24orIsoMu27.root"

float pt = 40
// 
ScaleFactor muTrgSF;
// UL2018 muTrgSF
muTrgSF.init_ScaleFactor(file, "ZMass")

// eTrgSF.etaBinsH->GetXaxis()->GetBinLabel(3)
// eTrgSF.etaBinsH->GetXaxis()->GetBinLowEdge(3)
// eTrgSF.etaBinsH->GetXaxis()->GetBinUpEdge(3)


muTrgSF.etaBinsH->GetXaxis()->FindFixBin(0.5)
muTrgSF.etaBinsH->GetXaxis()->GetBinLabel(muTrgSF.etaBinsH->GetXaxis()->FindFixBin(0.5))
// std::string etaLabel = muTrgSF.FindEtaLabel(0.5, "data") // Eta0p0to0p8 : correct

// auto& graph = *(eTrgSF.eff_data["Eta0p0to0p8"])
// graph.Draw()

// int ptbin = eTrgSF.FindPtBin(eTrgSF.eff_data, etaLabel, pt); 


muTrgSF.get_EfficiencyData(40, 0.5)