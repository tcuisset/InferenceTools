#ifndef AK8GenInterface_h
#define AK8GenInterface_h

// ------------------------------------------------------------ //
//                                                              //
//   class AK8GenInterface                                      //
//                                                              //
//                                                              //
//   Author: Francesco Brivio (INFN Milano-Bicocca)             //
//   Date  : June 2024                                          //
//                                                              //
// ------------------------------------------------------------ //

#include <vector>
#include <ROOT/RVec.hxx>

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

struct gen_match_output {
    std::vector<bool> Ak8_Zbb_matches;
    std::vector<bool> Ak8_Hbb_matches;
};

class AK8GenInterface {

  public:
    AK8GenInterface();
    ~AK8GenInterface();

    gen_match_output get_ak8_genmatch_info(
        iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
        iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
        fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass,
        fRVec FatJet_pt, fRVec FatJet_eta, fRVec FatJet_phi, fRVec FatJet_mass);
};

struct gen_jets_output {
    std::vector<int> genXbb_idxs; // index into GenParticle collection of the Z/H
    std::vector<int> genXbb_pdg; // PDG ID (ie 25 or 23)
    std::vector<std::vector<int>> genXbb_daughtersIdxs; // indices into GenParticle collection of the daughters of the X (ie b taus)
};

gen_jets_output get_jet_gen_info(
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass);


bool findGenZDecayToBB(iRVec GenPart_pdgId, iRVec GenPart_genPartIdxMother); // Check if we have a Z->bb decay from GenParticle collection. Returns true if such a decay exists

/** Info about a X->tautau candidate */
struct Xtautau_gen {
    int GenPartIdx; // index into GenParticle collection of the Z/H decaying to tautau
    int pdgId; // PDG ID (ie 25 or 23)
    std::vector<int> daughtersIdxs; // indices into GenParticle collection of the daughters taus of the X 
};
using gen_Xtautau_output = std::vector<Xtautau_gen>;
gen_Xtautau_output gen_Xtautau_info(iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother);

struct gen_daus_output {
    int genPairType; // -1=unmatched, -3=wrong GenVisTau match (deltaR very small between 2 taus)
    int genXtautau_GenPartIdx = -1;
    int genXtautau_pdgId = 0;

    int genTau1_GenPartIdx = -1; // index into GenPart collection of the tau decaying from Z/H
    int genTau2_GenPartIdx = -1;
    int genDau1_GenPartIdx = -1; // index into GenPart collection of the electron/muon decaying from the tau. Only for leptonic tau decays
    int genDau2_GenPartIdx = -1;
    int genDau1_GenVisTauIdx = -1; // index into GenVisTau collection (visible hadronic decay products of a tau) of the tau decaying from Z/H, only filled for hadronic taus
    int genDau2_GenVisTauIdx = -1;
};

gen_daus_output gen_daus_info(gen_Xtautau_output const& Xtautau_out, 
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    iRVec GenVisTau_genPartIdxMother);

#endif // AK8GenInterface_h
