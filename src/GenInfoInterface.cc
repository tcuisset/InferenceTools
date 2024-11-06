#include "Tools/Tools/interface/GenInfoInterface.h"
#include <iostream>
#include <algorithm>
#include <TLorentzVector.h>
#include <cassert>
#include <string>

// Constructor
AK8GenInterface::AK8GenInterface() {}

// Destructor
AK8GenInterface::~AK8GenInterface() {}

// Get ak8 match info
gen_match_output AK8GenInterface::get_ak8_genmatch_info(
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass,
    fRVec FatJet_pt, fRVec FatJet_eta, fRVec FatJet_phi, fRVec FatJet_mass)
{
    // Declare output
    gen_match_output genMatchOut;

    // Set all ouput to false
    for (size_t i = 0; i < FatJet_pt.size(); i++)
    {
        genMatchOut.Ak8_Zbb_matches.push_back(false);
        genMatchOut.Ak8_Hbb_matches.push_back(false);
    }

    // Declare vectors of gen H/Z->bb idxs
    std::vector<int> genHbb_idxs, genZbb_idxs;

    // Loop on gen particles to fill the genH and genZ idx vectors
    for (size_t idx = 0; idx < GenPart_pdgId.size(); idx++)
    {
        // Particle gen info
        int pdg    = GenPart_pdgId[idx];
        int mthIdx = GenPart_genPartIdxMother[idx];
        int mthPdg = (mthIdx > -1) ? GenPart_pdgId[mthIdx] : 0;

        bool isFirst       = (GenPart_statusFlags[idx] & (1<<12)) ? true : false;
        bool isHardProcess = (GenPart_statusFlags[idx] & (1<<7))  ? true : false;

        // If it's a b-quark
        if (abs(pdg) == 5 && isHardProcess && isFirst)
        {
            // Check whether it comes from H...
            if (abs(mthPdg) == 25)
            {
                // If motherH idx already saved (due to other b quark), continue
                if (std::find(genHbb_idxs.begin(), genHbb_idxs.end(), mthIdx) != genHbb_idxs.end())
                    continue;
                // else, save it
                else
                    genHbb_idxs.push_back(mthIdx);
            }
            // ...or Z
            else if (abs(mthPdg) == 23)
            {
                // If motherZ idx already saved (due to other b quark), continue
                if (std::find(genZbb_idxs.begin(), genZbb_idxs.end(), mthIdx) != genZbb_idxs.end())
                    continue;
                // else, save it
                else
                    genZbb_idxs.push_back(mthIdx);
            }
        }
    }

    // Now loop on FatJets and do the geometric matching (dR < 0.8) with the gen H/Z->bb
    for (size_t i = 0; i < FatJet_pt.size(); i++)
    {
        // Reco ak8 jet
        TLorentzVector fatjet_tlv;
        fatjet_tlv.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i], FatJet_mass[i]);

        // Matching with H->bb
        for (auto idx : genHbb_idxs)
        {
            TLorentzVector gen_Hbb;
            gen_Hbb.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
            if (fatjet_tlv.DeltaR(gen_Hbb) < 0.8)
                genMatchOut.Ak8_Hbb_matches[i] = true;
        }

        // Matching with Z->bb
        for (auto idx : genZbb_idxs)
        {
            TLorentzVector gen_Zbb;
            gen_Zbb.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
            if (fatjet_tlv.DeltaR(gen_Zbb) < 0.8)
                genMatchOut.Ak8_Zbb_matches[i] = true;
        }
    }

    return genMatchOut;
}



gen_jets_output get_jet_gen_info(
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass) {
    // Declare vectors of gen H/Z->bb idxs
    // X stands for Z or H
    gen_jets_output out;

    // Loop on gen particles 
    for (size_t idx = 0; idx < GenPart_pdgId.size(); idx++)
    {
        // Particle gen info
        int pdg    = GenPart_pdgId[idx];
        int mthIdx = GenPart_genPartIdxMother[idx];
        int mthPdg = (mthIdx > -1) ? GenPart_pdgId[mthIdx] : 0;

        bool isFirst       = (GenPart_statusFlags[idx] & (1<<12)) ? true : false;
        bool isHardProcess = (GenPart_statusFlags[idx] & (1<<7))  ? true : false;

        // If it's a b-quark
        if (abs(pdg) == 5 && isHardProcess && isFirst)
        {
            // Check whether it comes from H or Z
            if (abs(mthPdg) == 25 || abs(mthPdg) == 23) {
                auto mother_X = std::find(out.genXbb_idxs.begin(), out.genXbb_idxs.end(), mthIdx);
                if (mother_X == out.genXbb_idxs.end()) {
                    out.genXbb_idxs.push_back(mthIdx);
                    out.genXbb_pdg.push_back(mthPdg);
                    out.genXbb_daughtersIdxs.emplace_back();
                    mother_X = out.genXbb_idxs.end()-1;
                }
                out.genXbb_daughtersIdxs.at(mother_X-out.genXbb_idxs.begin()).push_back(idx);
            }
        }
    }

    if (out.genXbb_idxs.size() > 1) {
        std::cerr << "WARNING : more than one Z/H particle with b quark as daughter found in the event." << std::endl;
    }

    return out;
}


// Find all the X->tautau (where X=Z or H) at gen level
gen_Xtautau_output gen_Xtautau_info(iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother) {
    gen_Xtautau_output out;
    
    for (std::size_t i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
        int pdg    = GenPart_pdgId[i_gen];

        bool isLastCopy       = (GenPart_statusFlags[i_gen] & (1<<13)) ? true : false;
        bool isHardProcess = (GenPart_statusFlags[i_gen] & (1<<7))  ? true : false;
        bool isFromHardProcess = (GenPart_statusFlags[i_gen] & (1<<8))  ? true : false;
        // mthIdx == -1 /* it is the incoming parton */ || 
        // We want the last copy (ie after FSR) as that is the one used for GenVisTau_mother
        if (!(isHardProcess || isFromHardProcess) || !isLastCopy) continue;

        if (fabs(pdg) == 15 /*tau*/) {
            // Check whether this tau comes from H or Z
            // We might need to go up the chain in case of FSR (we can have H -> tau -> tau + gamma, and we want the last tau)
            int mthIdx = GenPart_genPartIdxMother[i_gen];
            while (mthIdx >= 0) {
                int mthPdg = GenPart_pdgId[mthIdx];
                if (abs(mthPdg) == 25 || abs(mthPdg) == 23) { // This tau indeed decays from a Z/H
                    // Check if we already encountered this X before (from another tau)
                    auto mother_X = std::find_if(out.begin(), out.end(), [mthIdx](Xtautau_gen const& x) { return x.GenPartIdx == mthIdx; });
                    if (mother_X == out.end()) {
                        // First time encountering this X : fill in the info
                        Xtautau_gen X;
                        X.GenPartIdx = mthIdx;
                        X.pdgId = mthPdg;
                        out.push_back(X);
                        mother_X = out.end()-1;
                    }
                    mother_X->daughtersIdxs.push_back(i_gen); // here we put the inital tau with isLastCopy
                    break;
                }
                mthIdx = GenPart_genPartIdxMother[mthIdx];
            }
        }
    }

    return out;
}

// https://github.com/cms-tau-pog/TauFW/blob/master/PicoProducer/python/analysis/GenFilterMuTau.py
/**
 Find a tau decay particle. 
 If tau->e or tau->mu, return the index into GenPart collection of the electron/muon. 
 If hadronic decay return -1
 Could not determine decay : -2 or -3
 Works even with a tau that is before FSR, will recurse down the daughters
*/
static int findTauDecay(int tau_idx, iRVec const& GenPart_pdgId, iRVec const& GenPart_genPartIdxMother) {
    int current_idx = tau_idx;
    for (int depth = 0; depth < 5; depth++) { // max depth of 5 (just in case...)
        bool foundDaughterToRecurse = false;
        for (std::size_t j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen ++) {
            if (GenPart_genPartIdxMother[j_gen] == current_idx) {
                // Found a daughter
                int pdg = GenPart_pdgId[j_gen];
                if (std::abs(pdg) == 11 || std::abs(pdg) == 13)
                    return j_gen; // e or mu
                else if (pdg == 22 || std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16)
                    continue; // photon / neutrino-e / neutrino-mu / neutrino-tau : ignore, keep looking at daughters
                else if (std::abs(pdg) == 15) {
                    current_idx = j_gen; // tau : recurse down, find daughters of this new tau
                    foundDaughterToRecurse = true;
                    break;
                }
                else if (std::abs(pdg) >= 100) {
                    return -1; // hadronic decay
                } else {
                    // Other stuff (quarks, ....) : ignore
                    std::cerr << "Found PDG="<< pdg << " decaying from " << current_idx << " (originally from " << tau_idx << ")" << std::endl;
                }
                
            }
        }
        if (!foundDaughterToRecurse)
            return -2; // Could not determine decay
    }
    return -3; // Could not determine decay (too deep)
}

gen_daus_output gen_daus_info(gen_Xtautau_output const& Xtautau_out, 
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    iRVec GenVisTau_genPartIdxMother) {
    
    gen_daus_output out = {}; // initialize to defaults

    // for (Xtautau_gen const& Xtautau : Xtautau_out) {
    //     std::cerr << "Xtautau " << Xtautau.GenPartIdx << "(PDG=" << Xtautau.pdgId << ") Children ";
    //     for (auto daughter : Xtautau.daughtersIdxs) {
    //         std::cerr << daughter << "(PDG=" << GenPart_pdgId.at(daughter) << ") ";
    //     }
    // }
    // std::cerr << std::endl;

    if (std::count_if(Xtautau_out.begin(), Xtautau_out.end(), [](Xtautau_gen const& Xtautau) { return Xtautau.daughtersIdxs.size() > 2;}) != 0)
        throw std::logic_error("More than 2 taus decayed from a Z/H");
    
    if (std::count_if(Xtautau_out.begin(), Xtautau_out.end(), [](Xtautau_gen const& Xtautau) { return Xtautau.daughtersIdxs.size() == 2;}) > 1)
        throw std::logic_error("More than 2 Z/H candidates decaying to tautau found");

    auto Xtautau_it = std::find_if(Xtautau_out.begin(), Xtautau_out.end(), [](Xtautau_gen const& Xtautau) { return Xtautau.daughtersIdxs.size() == 2;});
    if (Xtautau_it == Xtautau_out.end()) {
        std::cout << "WARNING no X->tautau candidate found" << std::endl;
        //throw std::runtime_error("no Xtautau");
        return out; // No X->tautau candidate found
    }
    Xtautau_gen const& Xtautau = *Xtautau_it;

    // std::cout << Xtautau.GenPartIdx << std::endl;
 
    out.genXtautau_GenPartIdx = Xtautau.GenPartIdx; // index into GenPart of the X decaying to tautau
    out.genXtautau_pdgId = Xtautau.pdgId;
    std::vector<int> const& taus_idx = Xtautau.daughtersIdxs; // index into GenPart of the taus
    if (taus_idx.size() != 2) {
        std::cerr << "WRONG NB OF TAUS" << std::endl;
        throw std::logic_error("WRONG NB OF TAUS");
    }
    std::vector<int> taus_daughter_GenPartIdx(taus_idx.size(), -1); // If we find a lepton daughter of the tau, we store its GenPart idx there. Otherwise keep it at -1
    std::vector<int> taus_daughter_pdg(taus_idx.size(), 0); // If we find a lepton daughter of the tau, we store its pdgId there. Otherwise keep it at 0
    std::vector<int> taus_genVisTauIdx(taus_idx.size(), -1); // If we find a GenVisTau daughter of the tau, store its GenVisTau idx otherwise -1

    // Sometimes there are weird things in GenParticle collection (FSR, etc) so we need a recusrive search
    // Looping over taus daughters from X to check if they are the mother of our j_gen
    for (std::size_t tau_index = 0; tau_index < taus_idx.size(); tau_index++) {
        taus_daughter_GenPartIdx.at(tau_index) = findTauDecay(taus_idx.at(tau_index), GenPart_pdgId, GenPart_genPartIdxMother);
        if (taus_daughter_GenPartIdx.at(tau_index) >= 0)
            taus_daughter_pdg.at(tau_index) = GenPart_pdgId[taus_daughter_GenPartIdx.at(tau_index)];
    }


    // looking for GenVisTau
    for (std::size_t i_genvistau = 0; i_genvistau < GenVisTau_genPartIdxMother.size(); i_genvistau ++) {
        // Looping over taus daughters from X to check if they are the mother of our GenVisTau
        for (std::size_t tau_index = 0; tau_index < taus_idx.size(); tau_index++) {
            int tau_GenPartIdx = taus_idx.at(tau_index);
            if (GenVisTau_genPartIdxMother[i_genvistau] == tau_GenPartIdx) {
                taus_genVisTauIdx.at(tau_index) = i_genvistau;
            }
        }
    }

    bool swap = false; // We want electron/muon to be dau1, so in case it is not the case set swap to true
    // NB : It is possible that no GenVisTau is found correspoding to a hadronic tau decay (usually due to out of eta range)
    if (std::abs(taus_daughter_pdg[0]) == 13 && taus_daughter_GenPartIdx[1] == -1 && taus_genVisTauIdx[0] == -1) {
        out.genPairType = 0; // mutau
    } else if (std::abs(taus_daughter_pdg[1]) == 13 && taus_daughter_GenPartIdx[0] == -1 && taus_genVisTauIdx[1] == -1) {
        out.genPairType = 0; // mutau
        swap = true;
    } else if (std::abs(taus_daughter_pdg[0]) == 11 && taus_daughter_GenPartIdx[1] == -1 && taus_genVisTauIdx[0] == -1) {
        out.genPairType = 1; // etau
    } else if (std::abs(taus_daughter_pdg[1]) == 11 && taus_daughter_GenPartIdx[0] == -1 && taus_genVisTauIdx[1] == -1) {
        out.genPairType = 1; // etau
        swap = true;
    } else if (taus_daughter_GenPartIdx[0] == -1 && taus_daughter_GenPartIdx[1] == -1) {
        out.genPairType = 2; // tautau
    } else if (taus_daughter_pdg[0] != 0 && taus_daughter_pdg[1] != 0 && taus_genVisTauIdx[0] == -1 && taus_genVisTauIdx[1] == -1) {
        out.genPairType = 3; // other decay mode
    } else {
        if ((taus_daughter_pdg[0] != 0 && taus_genVisTauIdx[0]>=0) || (taus_daughter_pdg[1] != 0 && taus_genVisTauIdx[1]>=0)) {
            // This case is when a tau decays leptonically according to GenParticle but there is a GenVisTau matched to it
            // Happens when deltaR(tau, tau)<1e-4 at gen level, because the matching is done based on gen-level deltaR !
            // https://github.com/cms-sw/cmssw/blob/7dc755530ad751991ea10bacdc0dbaecff716186/PhysicsTools/HepMCCandAlgos/plugins/GenVisTauProducer.cc#L78
            // We try to recover these cases by reassigning the GenVisTau to the dau2
            std::cout << "WARNING : tau decaying leptonically having a GenVisTau matched to them. Should be very rare, when deltaR<1e-4 between the two taus at gen level." << std::endl;
            if (std::abs(taus_daughter_pdg[0]) == 13 && taus_daughter_GenPartIdx[1] == -1 && taus_genVisTauIdx[0] >= 0 && taus_genVisTauIdx[1] == -1) {
                out.genPairType = 0; // mutau
            } else if (std::abs(taus_daughter_pdg[1]) == 13 && taus_daughter_GenPartIdx[0] == -1  && taus_genVisTauIdx[1] >= 0 && taus_genVisTauIdx[0] == -1) {
                out.genPairType = 0; // mutau
                swap = true;
            } else if (std::abs(taus_daughter_pdg[0]) == 11 && taus_daughter_GenPartIdx[1] == -1 && taus_genVisTauIdx[0] >= 0 && taus_genVisTauIdx[1] == -1) {
                out.genPairType = 1; // etau
            } else if (std::abs(taus_daughter_pdg[1]) == 11 && taus_daughter_GenPartIdx[0] == -1  && taus_genVisTauIdx[1] >= 0 && taus_genVisTauIdx[0] == -1) {
                out.genPairType = 1; // etau
                swap = true;
            }
            else
                throw std::runtime_error(std::string("Wrong tau decay modes : ") + "Xtautau GenPartIdx=" + std::to_string(Xtautau.GenPartIdx) +  " daughter PDGs = " + std::to_string(taus_daughter_pdg[0]) + "," + std::to_string(taus_daughter_pdg[1])
                + " ; genVisTaus idxs = " + std::to_string(taus_genVisTauIdx[0]) + "," + std::to_string(taus_genVisTauIdx[1]));
            if (swap) {
                out.genDau1_GenPartIdx = taus_daughter_GenPartIdx[1];
                out.genDau2_GenPartIdx = taus_daughter_GenPartIdx[0];
                out.genDau1_GenVisTauIdx = -1;
                out.genDau2_GenVisTauIdx = taus_genVisTauIdx[1];
            } else {
                out.genDau1_GenPartIdx = taus_daughter_GenPartIdx[0];
                out.genDau2_GenPartIdx = taus_daughter_GenPartIdx[1];
                out.genDau1_GenVisTauIdx = -1;
                out.genDau2_GenVisTauIdx = taus_genVisTauIdx[0];
            }
            return out;
        }
        else {
            throw std::runtime_error(std::string("Wrong tau decay modes : ") + "Xtautau GenPartIdx=" + std::to_string(Xtautau.GenPartIdx) +  " daughter PDGs = " + std::to_string(taus_daughter_pdg[0]) + "," + std::to_string(taus_daughter_pdg[1])
            + " ; genVisTaus idxs = " + std::to_string(taus_genVisTauIdx[0]) + "," + std::to_string(taus_genVisTauIdx[1]));
        }
    }

    if (swap) {
        out.genTau1_GenPartIdx = taus_idx[1];
        out.genTau1_GenPartIdx = taus_idx[0];
        out.genDau1_GenPartIdx = taus_daughter_GenPartIdx[1];
        out.genDau2_GenPartIdx = taus_daughter_GenPartIdx[0];
        out.genDau1_GenVisTauIdx = taus_genVisTauIdx[1];
        out.genDau2_GenVisTauIdx = taus_genVisTauIdx[0];
    } else {
        out.genTau1_GenPartIdx = taus_idx[0];
        out.genTau1_GenPartIdx = taus_idx[1];
        out.genDau1_GenPartIdx = taus_daughter_GenPartIdx[0];
        out.genDau2_GenPartIdx = taus_daughter_GenPartIdx[1];
        out.genDau1_GenVisTauIdx = taus_genVisTauIdx[0];
        out.genDau2_GenVisTauIdx = taus_genVisTauIdx[1];
    }

    return out;
}
