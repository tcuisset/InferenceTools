""" Modules using generator information : for gen studies, generator filter, etc """
import os
from analysis_tools.utils import import_root
from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()

class AK8GenRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            if not os.getenv("_GenInfoInterface"):
                os.environ["_GenInfoInterface"] = "_GenInfoInterface"
                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/GenInfoInterface.h".format(base))
            
            if not os.getenv("_AK8GenRDFProducer"):
                os.environ["_AK8GenRDFProducer"] = "_AK8GenRDFProducer"
                ROOT.gInterpreter.Declare("""
                    auto AK8Gen = AK8GenInterface();
                """)

    def run(self, df):
        if not self.isMC:
            return df, []

        variables = ["Ak8_Zbb_matches", "Ak8_Hbb_matches"]

        df = df.Define("ak8_gen_results", "AK8Gen.get_ak8_genmatch_info("
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, "
                            "GenPart_phi, GenPart_mass, "
                            "FatJet_pt{0}, FatJet_eta, FatJet_phi, FatJet_mass{0})"
                            .format(self.jet_syst))

        branches = []
        for var in variables:
            df = df.Define(f"gen{var}", f"ak8_gen_results.{var}")
            branches.append(f"gen{var}")

        return df, branches


def AK8GenRDF(**kwargs):
    """
    Module to store genMatch information between AK8 and  H/Z->bb

    :param isMC: flag of the dataset being MC or data
    :type : bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: AK8GenRDF
            path: Tools.Tools.AK8Gen
            parameters:
                isMC: self.dataset.process.isMC
    """
    return lambda: AK8GenRDFProducer(**kwargs)


class GenJetsRDFProducer:
    def __init__(self, *args, **kwargs):
        #super().__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            if not os.getenv("_GenInfoInterface"):
                os.environ["_GenInfoInterface"] = "_GenInfoInterface"
                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/GenInfoInterface.h".format(base))

    def run(self, df):
        if not self.isMC:
            return df, []

        df = df.Define("jet_gen_info", "get_jet_gen_info("
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, "
                            "GenPart_phi, GenPart_mass)")
        
        df = df.Define("genXbb_GenPartIdx", "jet_gen_info.genXbb_idxs.size()> 0 ? jet_gen_info.genXbb_idxs.at(0) : -1")
        df = df.Define("genXbb_pgdId", "jet_gen_info.genXbb_idxs.size()> 0 ? jet_gen_info.genXbb_pdg.at(0) : -1")
        df = df.Define("gen_b1_GenPartIdx", "jet_gen_info.genXbb_daughtersIdxs.size()> 0 ? jet_gen_info.genXbb_daughtersIdxs.at(0).at(0) : -1")
        df = df.Define("gen_b2_GenPartIdx", "(jet_gen_info.genXbb_daughtersIdxs.size()> 0 && jet_gen_info.genXbb_daughtersIdxs.at(0).size()>1) ? jet_gen_info.genXbb_daughtersIdxs.at(0).at(1) : -1")
        branches = ["genXbb_GenPartIdx", "genXbb_pgdId", "gen_b1_GenPartIdx", "gen_b2_GenPartIdx"]
        
        for x in ["Xbb", "_b1", "_b2"]:
            df = df.Define(f"gen{x}_pt", f"gen{x}_GenPartIdx >= 0 ? GenPart_pt[gen{x}_GenPartIdx] : -1")
            df = df.Define(f"gen{x}_eta", f"gen{x}_GenPartIdx >= 0 ? GenPart_eta[gen{x}_GenPartIdx] : -999")
            df = df.Define(f"gen{x}_phi", f"gen{x}_GenPartIdx >= 0 ? GenPart_phi[gen{x}_GenPartIdx] : -999")
            df = df.Define(f"gen{x}_mass", f"gen{x}_GenPartIdx >= 0 ? GenPart_mass[gen{x}_GenPartIdx] : -999")
            branches += [f"gen{x}_pt", f"gen{x}_eta", f"gen{x}_phi", f"gen{x}_mass"]
        
        df = df.Define("gen_deltaR_bb", "gen_b1_GenPartIdx>=0&&gen_b2_GenPartIdx>=0 ? deltaR(gen_b1_eta, gen_b1_phi, gen_b2_eta, gen_b2_phi) : -1")
        branches.append("gen_deltaR_bb")
        return df, branches


def GenJetsRDF(**kwargs):
    """
    Get gen information on Z/H->bb

    :param isMC: flag of the dataset being MC or data
    :type : bool

    Branches : 
     - genXbb_GenPartIdx : index into GenPart collection of the H/Z particle ddecaying into bb
     - genXbb_pgdId : pdg ID of the H/Z (ie 23 or 25)
     - gen_b1_GenPartIdx : index into GenPart collection of one of the b quarks decaying from H/Z
     - gen_b2_GenPartIdx : the other b
     - gen_deltaR_bb: deltaR at gen level between the 2 b quarks (-1 if no genmatch)
    """
    return lambda: GenJetsRDFProducer(**kwargs)


class GenLeptonRDFProducer():
    def __init__(self, *args, **kwargs):
        #super().__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            if not os.getenv("_GenInfoInterface"):
                os.environ["_GenInfoInterface"] = "_GenInfoInterface"
                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/GenInfoInterface.h".format(base))

    def run(self, df):
        if not self.isMC:
            return df, []
        branches = []

        df = df.Define("gen_Xtautau", "gen_Xtautau_info("
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother)")
        df = df.Define("gen_daus", "gen_daus_info(gen_Xtautau, "
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother, GenVisTau_genPartIdxMother)")
        
        import cppyy
        for branch in dir(cppyy.gbl.gen_daus_output): # this lists all attributes of the C++ object gen_daus_output
            if branch.startswith("_"):
                continue
            df = df.Define(branch, f"gen_daus.{branch}")
            branches.append(branch)
        
        def makeFourMomentumBranches(x):
            nonlocal df
            nonlocal branches
            df = df.Define(f"{x}_pt", f"{x}_GenPartIdx >= 0 ? GenPart_pt[{x}_GenPartIdx] : -1")
            df = df.Define(f"{x}_eta", f"{x}_GenPartIdx >= 0 ? GenPart_eta[{x}_GenPartIdx] : -999")
            df = df.Define(f"{x}_phi", f"{x}_GenPartIdx >= 0 ? GenPart_phi[{x}_GenPartIdx] : -999")
            df = df.Define(f"{x}_mass", f"{x}_GenPartIdx >= 0 ? GenPart_mass[{x}_GenPartIdx] : -999")
            branches += [f"{x}_pt", f"{x}_eta", f"{x}_phi", f"{x}_mass"]
        makeFourMomentumBranches("genXtautau")
        makeFourMomentumBranches("genTau1") # the tau itself (incl. neutrinos)
        makeFourMomentumBranches("genTau2")

        def makeFourMomentumBranchesGenVisTau(x):
            """ make branches choosing between GenVisTauIdx and GenPartIdx, whichever is available (ie hadronic or leptonic tau) """
            nonlocal df
            nonlocal branches
            df = df.Define(f"{x}_pt", f"{x}_GenVisTauIdx >= 0 ? GenVisTau_pt[{x}_GenVisTauIdx] : ({x}_GenPartIdx >= 0 ? GenPart_pt[{x}_GenPartIdx] : -1)")
            df = df.Define(f"{x}_eta", f"{x}_GenVisTauIdx >= 0 ? GenVisTau_eta[{x}_GenVisTauIdx] : ({x}_GenPartIdx >= 0 ? GenPart_eta[{x}_GenPartIdx] : -999)")
            df = df.Define(f"{x}_phi", f"{x}_GenVisTauIdx >= 0 ? GenVisTau_phi[{x}_GenVisTauIdx] : ({x}_GenPartIdx >= 0 ? GenPart_phi[{x}_GenPartIdx] : -999)")
            df = df.Define(f"{x}_mass", f"{x}_GenVisTauIdx >= 0 ? GenVisTau_mass[{x}_GenVisTauIdx] : ({x}_GenPartIdx >= 0 ? GenPart_mass[{x}_GenPartIdx] : -999)")
            df = df.Define(f"{x}_vis_status", f"{x}_GenVisTauIdx >= 0 ? GenVisTau_status[{x}_GenVisTauIdx] : -1") # Hadronic tau decay mode. 0=OneProng0PiZero, 1=OneProng1PiZero, 2=OneProng2PiZero, 10=ThreeProng0PiZero, 11=ThreeProng1PiZero, 15=Other
            branches += [f"{x}_pt", f"{x}_eta", f"{x}_phi", f"{x}_mass", f"{x}_vis_status"]
        makeFourMomentumBranchesGenVisTau("genDau1")
        makeFourMomentumBranchesGenVisTau("genDau2")
        
        # if genPairType==2, index into GenVisTau collection, otherwise index into GenPart collection (branch to be used for (boosted)Tau_genPartIdx checks)
        df = df.Define("genDau1_GenPartVisIdx", "genPairType == 2 ? genDau1_GenVisTauIdx : genDau1_GenPartIdx")
        df = df.Define("genDau2_GenPartVisIdx", "(genPairType >= 0 && genPairType <= 2) ? genDau2_GenVisTauIdx : genDau2_GenPartIdx")
        branches += ["genDau1_GenPartVisIdx", "genDau2_GenPartVisIdx"]
        
        df = df.Define("gen_deltaR_tautau", "genTau1_GenPartIdx>=0.&&genTau2_GenPartIdx>=0. ? deltaR(genTau1_eta, genTau1_phi, genTau2_eta, genTau2_phi) : -1")
        branches.append("gen_deltaR_tautau")
        df = df.Define("gen_deltaR_daudau", "genDau1_pt>-1.&&genDau2_pt>-1. ? deltaR(genDau1_eta, genDau1_phi, genDau2_eta, genDau2_phi) : -1")
        branches.append("gen_deltaR_daudau")

        return df, branches

def GenLeptonsRDF(**kwargs):
    """
    Get gen information on Z/H->tautau

    :param isMC: flag of the dataset being MC or data
    :type : bool

    General rules : genTau* include neutrinos from tau decay, genDau* only include visible energy
    Branches : 
     - genPairType : pairType using gen information : 0=mutau, 1=etau, 2=tautau, 3=other(ee,emu,muu), -1=no gen match found
     - genXtautau_GenPartIdx : index into GenPart collection of the H/Z particle ddecaying into tautau
     - genXtautau_pdgId : pdg ID of the H/Z (ie 23 or 25)
     - genXtautau_pt / eta / phi / mass
     
     - genTau1_GenPartIdx : index into GenPart collection of the tau decaying from Z/H
     - genTau1_pt / eta / phi / mass : generator tau information from GenParticle, including invisible tau energy
     - genDau1_GenPartIdx, genDau2_GenPartIdx : index into GenPart collection of the electron/muon decaying from the tau. Only for leptonic tau decays
     - genDau1_GenVisTauIdx, genDau2_GenVisTauIdx : index into GenVisTau collection (visible hadronic decay products of a tau) of the tau decaying from Z/H, only filled for hadronic taus
     - genDau1_pt / eta / phi / mass : gen info of the visible decay products of the tau (ie lepton 4-momentum for leptonic decay, GenVisTau info for hadronic taus)
     - genDau1_GenPartVisIdx : if genPairType==2, index into GenVisTau collection, otherwise index into GenPart collection (branch to be used for (boosted)Tau_genPartIdx checks)
     - genDau2_GenPartVisIdx : this is the same as genDau2_GenVisTauIdx, but is meant to be used for (boosted)Tau_genPartIdx checks (here in case of emu channels being implemented)
     
     - gen_deltaR_tautau : deltaR between the two taus at GenParticle level (incl. neutrinos)
     - gen_deltaR_daudau : deltaR between the visible taus (w/o neutrinos, using leptons/GenVisTau)
    """
    return lambda: GenLeptonRDFProducer(**kwargs)


class GenVariablesRDFProducer():
    def __init__(self, *args, **kwargs):
        #super().__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        # if self.isMC:
        #     if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
        #         ROOT.gSystem.Load("libToolsTools.so")

        #     if not os.getenv("_GenInfoInterface"):
        #         os.environ["AK8GenInterface"] = "AK8GenInterface"
        #         base = "{}/{}/src/Tools/Tools".format(
        #             os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        #         ROOT.gROOT.ProcessLine(".L {}/interface/GenInfoInterface.h".format(base))

    def run(self, df):
        if not self.isMC:
            return df, []
        branches = []

        df = df.Define("gen_deltaR_XX", "(genXtautau_GenPartIdx>=0 && genXbb_GenPartIdx>=0) ? deltaR(genXtautau_eta, genXtautau_phi, genXbb_eta, genXbb_phi) : -1.")
        branches.append("gen_deltaR_XX")
        return df, branches

def GenVariablesRDF(**kwargs):
    """
    Get gen information on Z/H

    :param isMC: flag of the dataset being MC or data
    :type : bool

    Branches : 
     - gen_deltaR_XX : deltaR at gen level between ZZ pair (or ZH or HH...)
    """
    return lambda: GenVariablesRDFProducer(**kwargs)


class BBTauTauFilterRDFProducer():
    """ Filter bbtautau events (ZZ/ZH) using LHE information
    !!!! Only works on NanoV12 or later (needs status=2 intermediate bosons in LHEPart collections, not available in NanoV9)
    It also depends on the way the sample is generated. There are 2 cases:
     - The Z/H is decayed by the matrix element (or by madspin), in which case the Z/H is LHEPart_status==2 and the decay products are LHEPart_status==1 (case of ZZTo2L2Q sample, madspin decay)
     - The H is decayed by Pythia, in which case the H is LHEPart_status==1 (on shell outgoing), and its decay products are not there. Case of ZH samples for the H only
    """
    def __init__(self, ProcType, isSigBBTT, isBkgBBTT, removeZH=False, *args, **kwargs):
        self.isSigBBTT = isSigBBTT
        self.isBkgBBTT = isBkgBBTT
        self.ProcType = ProcType
        self.removeZH = removeZH

        if self.isSigBBTT or self.isBkgBBTT:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            if not os.getenv("_GenInfoInterface"):
                os.environ["_GenInfoInterface"] = "_GenInfoInterface"
                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/GenInfoInterface.h".format(base))

    def run(self, df):
        if self.isSigBBTT or self.isBkgBBTT:
            df = df.Define("LHE_ZCount", "ROOT::VecOps::Sum(LHEPart_status == 2 && LHEPart_pdgId==23)")
            df = df.Define("LHE_HCount", "ROOT::VecOps::Sum(LHEPart_status >= 1 && LHEPart_pdgId==25)") # in case Higgs is decayed with Pythia, status==1. in case Matrix-element Higgs resonance, status==2
            df = df.Define("LHE_bCount", "ROOT::VecOps::Sum(LHEPart_status == 1 && ROOT::VecOps::abs(LHEPart_pdgId)==5 )")
            df = df.Define("LHE_tauCount", "ROOT::VecOps::Sum(LHEPart_status == 1 && ROOT::VecOps::abs(LHEPart_pdgId)==15 )")
            df = df.Define("LHE_emuCount", "ROOT::VecOps::Sum(LHEPart_status == 1 && (11 <= ROOT::VecOps::abs(LHEPart_pdgId) && ROOT::VecOps::abs(LHEPart_pdgId)<=14) )") # (also counts neutrinos)

            # define a new branch to check if it is or not a ZZ/ZH->bbtautau event
            if self.ProcType == "Zbb_Ztautau":
                print(" ### Running bbtautau filter for Zbb_Ztautau")
                # select : 
                # - 2 Z (we want ZZ on shell, reject continuum Z/gamma)
                # - at least 2 b (we can have more in case of ZZ+jet)
                # - exactly 2 taus
                # - no electron/muon (reject ZZ+jet with Z->ee,Z->tautau,jet=bb)
                genfilter = "LHE_ZCount == 2 && LHE_bCount >= 2 && LHE_tauCount == 2 && LHE_emuCount == 0"
            elif self.ProcType == "Zbb_Htautau" :
                print(" ### Running bb filter for Zbb_Htautau")
                # we cannot use LHE_isBBTauTau since that does not work for Pythia decays of the Higgs
                # because our dataset is Z-> anything & H->tautau (forced Pythia decay), we just have to filter for bb at LHE level.
                # but we cannot differentiate Z->bb,H->tautau,jet=uu from Z->uu,H->tautau,jet=bb just from LHE, need GenPart info
                # LHE_ZCount will always be 1 (with status 2), LHE_HCount will always be 1 (with status 1)
                # select
                # - at least 2 b (thus we can have either Z->bb,H->tautau or Z->anything,H->tautau,jet=bb)
                # - 2 b coming from a Z in GenParticle collection (reject Z->anything)
                genfilter = "LHE_bCount >= 2 && findGenZDecayToBB(GenPart_pdgId, GenPart_genPartIdxMother)"
                assert not self.removeZH
            elif self.ProcType == "Ztautau_Hbb":
                print(" ### Running tautau filter for Ztautau_Hbb")
                # same but dataset is Z->ll & H->bb
                # so find 2 LHE taus
                genfilter = "LHE_tauCount >= 2"
                assert not self.removeZH
            else:
                raise ValueError("BBTauTauFilterRDF not implemented for self.ProcType = ", self.ProcType)
            
            # filter the events with ZZ/ZH->bbtautau
            if self.isSigBBTT:
                # print(" ### DEBUG: isBBTT == 1")
                df = df.Filter(genfilter, "BBTauTauFilterRDF")
            # filter the events without ZZ/ZH->bbtautau
            elif self.isBkgBBTT:
                # print(" ### DEBUG: isBBTT == 0")
                df = df.Filter("!(" + genfilter + ")", "BBTauTauFilterRDF")
            
            if self.removeZH:
                print(" ### Removing ZH contribution")
                df = df.Filter("LHE_HCount == 0", "BBTauTauFilterRDF_removeZH")
                
            return df, ["LHE_ZCount", "LHE_HCount", "LHE_bCount", "LHE_tauCount", "LHE_emuCount"]
        
        else:
            return df, []

class BBTauTauFilterRDFProducerOld():
    """ Old bbtautau filter based on GenParticle collection """
    def __init__(self, ProcType, isSigBBTT, isBkgBBTT, *args, **kwargs):
        self.isSigBBTT = isSigBBTT
        self.isBkgBBTT = isBkgBBTT
        self.ProcType = ProcType
        # print(" ### DEBUG: isSigBBTT = {}".format(isSigBBTT))
        # print(" ### DEBUG: isBkgBBTT = {}".format(isBkgBBTT))

        if not os.getenv("_BBTauTauFilter"):
            os.environ["_BBTauTauFilter"] = "_BBTauTauFilter"

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint   = const ROOT::RVec<int>&;
                // Return {genPairType, genPart_idx for dau1, genPart_idx for  dau2}
                std::array<int, 3> GenPairType_Zbb_Ztautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother, Vint GenVisTau_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromZ = 0;
                    int n_tau_fromZ = 0;
                    int tau1_id = -1; // index into GenParticle collection of one of the tau
                    int tau2_id = -1;
                    
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5 /*b*/) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23/*Z*/)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15 /*tau*/) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23 /*Z*/)) {
                            n_tau_fromZ += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromZ > 1)) {
                        int dau1_gen_id = -1; // index into GenParticle collection of electron/muon decaying from gen tau nb 1
                        int dau1_pdg_id = -1;
                        int dau2_gen_id = -1;
                        int dau2_pdg_id = -1;
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            // Looking for children of taus
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = j_gen;
                                dau1_pdg_id = fabs(GenPart_pdgId.at(j_gen)); // electron or muon decaying from tau
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = j_gen;
                                dau2_pdg_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        int genDau1_genPartIdx = -1;
                        int genDau2_genPartIdx = -1;
                        if (((dau1_pdg_id == 11) && (dau2_pdg_id == 13)) || ((dau1_pdg_id == 13) && (dau2_pdg_id == 11))) 
                            GenPairType = 3; // decay modes not covered (e-mu)
                        else if (dau1_pdg_id == 11 || dau2_pdg_id == 11) {
                            GenPairType = 1; // electron
                            genDau1_genPartIdx = std::max(dau1_gen_id, dau2_gen_id); // One of dau1/2_gen_id will be filled, the other -1
                        } else if (dau1_pdg_id == 13 || dau2_pdg_id == 13) {
                            GenPairType = 0; // muon
                            genDau1_genPartIdx = std::max(dau1_gen_id, dau2_gen_id);
                        } else {
                            GenPairType = GenVisTau_genPartIdxMother.size() >= 2 ? 2 : -1; // di-hadronic
                        }
                        
                        for (int i_genvistau = 0; i_genvistau < GenVisTau_genPartIdxMother.size(); i_genvistau ++) {
                            if ((GenVisTau_genPartIdxMother[i_genvistau] == tau1_id || GenVisTau_genPartIdxMother[i_genvistau] == tau2_id)) {
                                // lepton-hadronic case
                                if (GenPairType >=0 && GenPairType <= 1)
                                    genDau2_genPartIdx = i_genvistau;
                                else if (GenPairType == 2) {
                                    // tau_h tau_h case
                                    if (genDau1_genPartIdx == -1)
                                        genDau1_genPartIdx = i_genvistau;
                                    else if (genDau2_genPartIdx == -1)
                                        genDau2_genPartIdx = i_genvistau;
                                    else
                                        std::cerr << "WARNING : Something weird is going on in tau decays" << std::endl;
                                }
                            }
                        }

                        return {GenPairType, genDau1_genPartIdx, genDau2_genPartIdx};
                    }
                    else {
                        return {-1, -1, -1};
                    }
                }
                
                int GenPairType_Zbb_Htautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromZ = 0;
                    int n_tau_fromH = 0;
                    int tau1_id = -1;
                    int tau2_id = -1;
                    int dau1_gen_id = -1;
                    int dau2_gen_id = -1;

                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_tau_fromH += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromH > 1)) {
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        if (((dau1_gen_id == 11) && (dau2_gen_id == 13)) || ((dau1_gen_id == 13) && (dau2_gen_id == 11))) GenPairType = 3; // decay modes not covered
                        else if ((dau1_gen_id == 11) || (dau2_gen_id == 11)) GenPairType = 0;
                        else if ((dau1_gen_id == 13) || (dau2_gen_id == 13)) GenPairType = 1;
                        else GenPairType = 2;
                        return GenPairType;
                    }
                    else {
                        return -1;
                    }
                }

                int GenPairType_Ztautau_Hbb (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromH = 0;
                    int n_tau_fromZ = 0;
                    int tau1_id = -1;
                    int tau2_id = -1;
                    int dau1_gen_id = -1;
                    int dau2_gen_id = -1;

                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_b_fromH += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_tau_fromZ += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromH > 1) && (n_tau_fromZ > 1)) {
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        if (((dau1_gen_id == 11) && (dau2_gen_id == 13)) || ((dau1_gen_id == 13) && (dau2_gen_id == 11))) GenPairType = 3; // decay modes not covered
                        else if ((dau1_gen_id == 11) || (dau2_gen_id == 11)) GenPairType = 0;
                        else if ((dau1_gen_id == 13) || (dau2_gen_id == 13)) GenPairType = 1;
                        else GenPairType = 2;
                        return GenPairType;
                    }
                    else {
                        return -1;
                    }
                }
            """)

    def run(self, df):
        if self.isSigBBTT or self.isBkgBBTT:

            # define a new branch to check if it is or not a ZZ/ZH->bbtautau event
            if self.ProcType == "Zbb_Ztautau":
                print(" ### Running bbtautau filter for Zbb_Ztautau")
                df = df.Define("GenPairType_idx_pair", """GenPairType_Zbb_Ztautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother,
                    GenVisTau_genPartIdxMother
                )""")
            elif self.ProcType == "Zbb_Htautau":
                print(" ### Running bbtautau filter for Zbb_Htautau")
                df = df.Define("GenPairType", """GenPairType_Zbb_Htautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            elif self.ProcType == "Ztautau_Hbb":
                print(" ### Running bbtautau filter for Ztautau_Hbb")
                df = df.Define("GenPairType", """GenPairType_Ztautau_Hbb(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            else:
                raise ValueError("BBTauTauFilterRDF not implemented for self.ProcType = ", self.ProcType)
            
            df = df.Define("GenPairType_old", "GenPairType_idx_pair[0]")
            df = df.Define("genDau1_genPartIdx_old", "GenPairType_idx_pair[1]")
            df = df.Define("genDau2_genPartIdx_old", "GenPairType_idx_pair[2]")
            
            # filter the events with ZZ/ZH->bbtautau
            if self.isSigBBTT:
                # print(" ### DEBUG: isBBTT == 1")
                df = df.Filter("GenPairType != -1", "BBTauTauFilterRDF")
            # filter the events without ZZ/ZH->bbtautau
            elif self.isBkgBBTT:
                # print(" ### DEBUG: isBBTT == 0")
                df = df.Filter("GenPairType == -1", "BBTauTauFilterRDF")
                
            return df, ["GenPairType_old", "genDau1_genPartIdx_old", "genDau2_genPartIdx_old"]
        
        else:
            return df, []

class BBTauTauFilterDummyRDFProducer():
    def run(self, df):
        df = df.Define("GenPairType", -1)
        return df, ["GenPairType"]

def BBTauTauFilterRDF(*args, **kwargs):
    """ removeZH : if True, remove Higgs contribution (meant to remove ZH from ZZTo2L2Q sample which is actually Z(Z/H) """

    ProcType = kwargs.pop("ProcType")
    isSigBBTT = kwargs.pop("isSigBBTT")
    isBkgBBTT = kwargs.pop("isBkgBBTT")

    return lambda: BBTauTauFilterRDFProducer(ProcType=ProcType, isSigBBTT=isSigBBTT, isBkgBBTT=isBkgBBTT, removeZH=kwargs.pop("removeZH", False), *args, **kwargs)

    if isSigBBTT or isBkgBBTT:
        return lambda: BBTauTauFilterRDFProducer(ProcType=ProcType, isSigBBTT=isSigBBTT, isBkgBBTT=isBkgBBTT, *args, **kwargs)
    else:
        return lambda: BBTauTauFilterDummyRDFProducer()



