import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst
from analysis_tools.utils import import_root

ROOT = import_root()


class HHProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(HHProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/HHKinFitInterface.h".format(base))
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("Hbb_pt%s" % self.systs, "F")
        self.out.branch("Hbb_eta%s" % self.systs, "F")
        self.out.branch("Hbb_phi%s" % self.systs, "F")
        self.out.branch("Hbb_mass%s" % self.systs, "F")

        self.out.branch("Htt_pt%s" % self.systs, "F")
        self.out.branch("Htt_eta%s" % self.systs, "F")
        self.out.branch("Htt_phi%s" % self.systs, "F")
        self.out.branch("Htt_mass%s" % self.systs, "F")

        self.out.branch("HH_pt%s" % self.systs, "F")
        self.out.branch("HH_eta%s" % self.systs, "F")
        self.out.branch("HH_phi%s" % self.systs, "F")
        self.out.branch("HH_mass%s" % self.systs, "F")

        self.out.branch("HH_svfit_pt%s" % self.systs, "F")
        self.out.branch("HH_svfit_eta%s" % self.systs, "F")
        self.out.branch("HH_svfit_phi%s" % self.systs, "F")
        self.out.branch("HH_svfit_mass%s" % self.systs, "F")

        self.out.branch("HHKinFit_mass%s" % self.systs, "F")
        self.out.branch("HHKinFit_chi2%s" % self.systs, "F")

        self.out.branch("VBFjj_mass%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaEta%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaPhi%s" % self.systs, "F")
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        htt_svfit = Object(event, "Htt_svfit")
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")
        jets = Collection(event, "Jet")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus) #FIXME FROM HERE
        bjet1, bjet2, bjet1_tlv, bjet2_tlv = self.get_bjets(event, jets)
        vbfjet1, vbfjet2, vbfjet1_tlv, vbfjet2_tlv = self.get_vbfjets(event, jets)
        met, met_tlv = self.get_met(event)

        met_tv = ROOT.TVector2(met_tlv.Px(), met_tlv.Py())

        cov = ROOT.TMatrixD(2, 2)
        cov[0][0] = met.covXX
        cov[0][1] = met.covXY
        cov[1][0] = met.covXY
        cov[1][1] = met.covYY

        # HH TLorentzVector
        htt_tlv = dau1_tlv + dau2_tlv
        hbb_tlv = bjet1_tlv + bjet2_tlv
        hh_tlv = htt_tlv + hbb_tlv

        # HH (with htt_svfit) TLorentzVector
        htt_svfit_tlv = ROOT.TLorentzVector()
        htt_svfit_tlv.SetPtEtaPhiM(
            eval("htt_svfit.pt%s" % self.systs), 
            eval("htt_svfit.eta%s" % self.systs), 
            eval("htt_svfit.phi%s" % self.systs), 
            eval("htt_svfit.mass%s" % self.systs),
        )
        hh_svfit_tlv = hbb_tlv + htt_svfit_tlv

        # HH Kin. Fit
        # constructor https://github.com/bvormwald/HHKinFit2/blob/CMSSWversion/HHKinFit2Scenarios/src/HHKinFitMasterHeavyHiggs.cpp#L27
        # kinFit = ROOT.HHKinFit2.HHKinFitMasterHeavyHiggs(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit = ROOT.HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit.addHypo(125, 125);
        results = kinFit.fit()
        HHKinFit_mass = results[0]
        HHKinFit_chi2 = results[1]

        # VBF jet composition
        if not vbfjet1:
            VBFjj_mass = -999.
            VBFjj_deltaEta = -999.
            VBFjj_deltaPhi = -999.
        else:
            VBFjj_tlv = vbfjet1_tlv + vbfjet2_tlv
            VBFjj_mass = VBFjj_tlv.M()
            VBFjj_deltaEta = vbfjet1_tlv.Eta() - vbfjet2_tlv.Eta()
            VBFjj_deltaPhi = vbfjet1_tlv.DeltaPhi(vbfjet2_tlv)

        self.out.fillBranch("Hbb_pt%s" % self.systs, hbb_tlv.Pt())
        self.out.fillBranch("Hbb_eta%s" % self.systs, hbb_tlv.Eta())
        self.out.fillBranch("Hbb_phi%s" % self.systs, hbb_tlv.Phi())
        self.out.fillBranch("Hbb_mass%s" % self.systs, hbb_tlv.M())

        self.out.fillBranch("Htt_pt%s" % self.systs, htt_tlv.Pt())
        self.out.fillBranch("Htt_eta%s" % self.systs, htt_tlv.Eta())
        self.out.fillBranch("Htt_phi%s" % self.systs, htt_tlv.Phi())
        self.out.fillBranch("Htt_mass%s" % self.systs, htt_tlv.M())

        self.out.fillBranch("HH_pt%s" % self.systs, hh_tlv.Pt())
        self.out.fillBranch("HH_eta%s" % self.systs, hh_tlv.Eta())
        self.out.fillBranch("HH_phi%s" % self.systs, hh_tlv.Phi())
        self.out.fillBranch("HH_mass%s" % self.systs, hh_tlv.M())

        self.out.fillBranch("HH_svfit_pt%s" % self.systs, hh_svfit_tlv.Pt())
        self.out.fillBranch("HH_svfit_eta%s" % self.systs, hh_svfit_tlv.Eta())
        self.out.fillBranch("HH_svfit_phi%s" % self.systs, hh_svfit_tlv.Phi())
        self.out.fillBranch("HH_svfit_mass%s" % self.systs, hh_svfit_tlv.M())

        self.out.fillBranch("HHKinFit_mass%s" % self.systs, HHKinFit_mass)
        self.out.fillBranch("HHKinFit_chi2%s" % self.systs, HHKinFit_chi2)

        self.out.fillBranch("VBFjj_mass%s" % self.systs, VBFjj_mass)
        self.out.fillBranch("VBFjj_deltaEta%s" % self.systs, VBFjj_deltaEta)
        self.out.fillBranch("VBFjj_deltaPhi%s" % self.systs, VBFjj_deltaPhi)

        return True


class HHKinFitRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.isMC = kwargs.pop("isMC")
        self.AnalysisType = AnalysisType
        super(HHKinFitRDFProducer, self).__init__(isMC=self.isMC, *args, **kwargs)

        if not os.getenv("_KINFIT"):
            os.environ["_KINFIT"] = "kinfit"

            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")
            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/HHKinFitInterface.h".format(base))
            ROOT.gInterpreter.Declare("""
                #include <TLorentzVector.h>
                #include "TVector.h"
                #include "TMatrixD.h"
                #include "TMath.h"

                using Vfloat = const ROOT::RVec<float>&;
                ROOT::RVec<double> compute_hhkinfit(
                        bool isMC, int pairType, int dau1_index, int dau2_index, int bjet1_index, int bjet2_index,
                        Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                        Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                        Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                        Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                        Vfloat jet_resolution, float met_pt, float met_phi,
                        float met_covXX, float met_covXY, float met_covYY, int target1, int target2) {
                    float dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;
                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_eta = muon_eta.at(dau1_index);
                        dau1_phi = muon_phi.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_eta = electron_eta.at(dau1_index);
                        dau1_phi = electron_phi.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_eta = tau_eta.at(dau1_index);
                        dau1_phi = tau_phi.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_eta = tau_eta.at(dau2_index);
                    dau2_phi = tau_phi.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);

                    auto bjet1_tlv = TLorentzVector();
                    auto bjet2_tlv = TLorentzVector();
                    auto dau1_tlv = TLorentzVector();
                    auto dau2_tlv = TLorentzVector();

                    dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
                    dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                    bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                        jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                    bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                        jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));

                    auto met_tv = TVector2(met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi));
                    auto cov = TMatrixD(2, 2);
                    cov[0][0] = met_covXX;
                    cov[0][1] = met_covXY;
                    cov[1][0] = met_covXY;
                    cov[1][1] = met_covYY;

                    double resolution1, resolution2;
                    if (isMC) {
                        resolution1 = bjet1_tlv.E() * jet_resolution.at(bjet1_index);
                        resolution2 = bjet2_tlv.E() * jet_resolution.at(bjet2_index);
                    } else {
                        resolution1 = -1;
                        resolution2 = -1;
                    }

                    auto kinFit = HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv,
                        met_tv, cov, resolution1, resolution2);
                    kinFit.addHypo(target1, target2);
                    return kinFit.fit();
                }
            """)

    def run(self, df):
        isMC = ("true" if self.isMC else "false")
        jet_resolution = "jet_pt_resolution"
        if not self.isMC:
            jet_resolution = "Jet_eta"  # placeholder

        # DEBUG
        if not self.AnalysisType:
            pp = "HH"; target1 = 125; target2 = 125
            print(" ### INFO: Running KinFit with default option for HH analysis")
        else:
            print(" ### INFO: Running KinFit with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      pp = "ZZ"; target1 = 91; target2 = 91
            elif self.AnalysisType == "Zbb_Htautau":    pp = "ZH"; target1 = 125; target2 = 91 # 1 is tautau, 2 is bb
            elif self.AnalysisType == "Ztautau_Hbb":    pp = "ZH"; target1 = 91; target2 = 125 

        branches = ["%sKinFit_mass%s" % (pp, self.systs), "%sKinFit_chi2%s" % (pp, self.systs)]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        df = df.Define("hhkinfit_result%s" % self.systs,
            "compute_hhkinfit({6}, pairType, dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx, "
                "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
                "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
                "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
                "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, {7},"
                "MET{5}_pt{4}, MET{5}_phi{4}, MET_covXX, MET_covXY, MET_covYY, {8}, {9})".format(
                    self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst, self.met_syst,
                    self.met_smear_tag, isMC, jet_resolution, target1, target2)
            ).Define("%sKinFit_mass%s" % (pp, self.systs), "hhkinfit_result%s[0]" % self.systs
            ).Define("%sKinFit_chi2%s" % (pp, self.systs), "hhkinfit_result%s[1]" % self.systs)
        return df, branches


class HHVarRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.AnalysisType = AnalysisType
        super(HHVarRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("_HHVAR_%s" % self.systs):
            os.environ["_HHVAR_%s" % self.systs] = "hhvar"
            ROOT.gInterpreter.Declare("""
                #include <TLorentzVector.h>
                using Vfloat = const ROOT::RVec<float>&;
                ROOT::RVec<double> get_hh_features%s(int pairType, int isBoosted,
                        int dau1_index, int dau2_index, int bjet1_index, int bjet2_index, int fatjet_index,
                        int vbfjet1_index, int vbfjet2_index,
                        Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                        Vfloat electron_pt, Vfloat electron_eta,
                        Vfloat electron_phi, Vfloat electron_mass,
                        Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                        Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                        Vfloat fatjet_pt, Vfloat fatjet_eta, Vfloat fatjet_phi, Vfloat fatjet_mass,
                        float met_pt, float met_phi,
                        double Htt_svfit_pt, double Htt_svfit_eta,
                        double Htt_svfit_phi, double Htt_svfit_mass) {

                    float dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;
                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_eta = muon_eta.at(dau1_index);
                        dau1_phi = muon_phi.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_eta = electron_eta.at(dau1_index);
                        dau1_phi = electron_phi.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_eta = tau_eta.at(dau1_index);
                        dau1_phi = tau_phi.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_eta = tau_eta.at(dau2_index);
                    dau2_phi = tau_phi.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);

                    auto dau1_tlv = TLorentzVector();
                    auto dau2_tlv = TLorentzVector();

                    dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
                    dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                    auto htt_tlv = dau1_tlv + dau2_tlv;
                    
                    
                    auto hbb_tlv = TLorentzVector();
                    if (!isBoosted) {
                        auto bjet1_tlv = TLorentzVector();
                        auto bjet2_tlv = TLorentzVector();
                        bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                            jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                        bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                            jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));
                        hbb_tlv = bjet1_tlv + bjet2_tlv;
                    }
                    else {
                        hbb_tlv.SetPtEtaPhiM(fatjet_pt.at(fatjet_index), fatjet_eta.at(fatjet_index),
                            fatjet_phi.at(fatjet_index), fatjet_mass.at(fatjet_index));
                    }
                    auto hh_tlv = htt_tlv + hbb_tlv;

                    auto met_tlv = TLorentzVector();
                    met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

                    auto htt_met_tlv = htt_tlv + met_tlv;

                    double hh_svfit_pt = -999., hh_svfit_eta = -999.;
                    double hh_svfit_phi = -999., hh_svfit_mass = -999.;
                    if (Htt_svfit_pt > 0) {
                        auto htt_svfit_tlv = TLorentzVector();
                        htt_svfit_tlv.SetPtEtaPhiM(Htt_svfit_pt, Htt_svfit_eta,
                            Htt_svfit_phi, Htt_svfit_mass);
                        auto hh_svfit_tlv = htt_svfit_tlv + hbb_tlv;
                        hh_svfit_pt = hh_svfit_tlv.Pt();
                        hh_svfit_eta = hh_svfit_tlv.Eta();
                        hh_svfit_phi = hh_svfit_tlv.Phi();
                        hh_svfit_mass = hh_svfit_tlv.M();
                    }
                    double vbfjj_mass = -999., vbfjj_deltaEta = -999., vbfjj_deltaPhi = -999.;
                    if (vbfjet1_index >= 0) {
                        auto vbfjet1_tlv = TLorentzVector();
                        auto vbfjet2_tlv = TLorentzVector();
                        vbfjet1_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet1_index), jet_eta.at(vbfjet1_index),
                            jet_phi.at(vbfjet1_index), jet_mass.at(vbfjet1_index));
                        vbfjet2_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet2_index), jet_eta.at(vbfjet2_index),
                            jet_phi.at(vbfjet2_index), jet_mass.at(vbfjet2_index));
                        auto vbfjj_tlv = vbfjet1_tlv + vbfjet2_tlv;
                        vbfjj_mass = vbfjj_tlv.M();
                        vbfjj_deltaEta = vbfjet1_tlv.Eta() - vbfjet2_tlv.Eta();
                        vbfjj_deltaPhi = vbfjet1_tlv.Phi() - vbfjet2_tlv.Phi();
                    }
                    return {
                        hbb_tlv.Pt(), hbb_tlv.Eta(), hbb_tlv.Phi(), hbb_tlv.M(),
                        htt_tlv.Pt(), htt_tlv.Eta(), htt_tlv.Phi(), htt_tlv.M(),
                        htt_met_tlv.Pt(), htt_met_tlv.Eta(), htt_met_tlv.Phi(), htt_met_tlv.M(),
                        hh_tlv.Pt(), hh_tlv.Eta(), hh_tlv.Phi(), hh_tlv.M(),
                        hh_svfit_pt, hh_svfit_eta, hh_svfit_phi, hh_svfit_mass,
                        vbfjj_mass, vbfjj_deltaEta, vbfjj_deltaPhi
                    };
                }
            """ % self.systs)

    def run(self, df):

        # DEBUG
        if not self.AnalysisType:
            p_b = "H"; p_t = "H"; pp = "HH"; p_sv = "H"
            print(" ### INFO: Running HHVar with default option for HH analysis")
        else:
            print(" ### INFO: Running HHVar with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      p_b = "Z"; p_t = "Z"; pp = "ZZ"; p_sv = "X"
            elif self.AnalysisType == "Zbb_Htautau":    p_b = "Z"; p_t = "H"; pp = "ZH"; p_sv = "X"
            elif self.AnalysisType == "Ztautau_Hbb":    p_b = "H"; p_t = "Z"; pp = "ZH"; p_sv = "X"

        features = ("{0}bb_pt{3},{0}bb_eta{3},{0}bb_phi{3},{0}bb_mass{3},"
            "{1}tt_pt{3},{1}tt_eta{3},{1}tt_phi{3},{1}tt_mass{3},"
            "{1}tt_met_pt{3},{1}tt_met_eta{3},{1}tt_met_phi{3},{1}tt_met_mass{3},"
            "{2}_pt{3},{2}_eta{3},{2}_phi{3},{2}_mass{3},"
            "{2}_svfit_pt{3},{2}_svfit_eta{3},{2}_svfit_phi{3},{2}_svfit_mass{3},"
            "VBFjj_mass{3},VBFjj_deltaEta{3},VBFjj_deltaPhi{3}".format(p_b, p_t, pp, self.systs))
        features = list(features.split(","))
        all_branches = df.GetColumnNames()
        if features[0] in all_branches:
            return df, []

        df = df.Define("hhfeatures%s" % self.systs, ("get_hh_features{6}(pairType, isBoosted, "
            "dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx, fatjet_JetIdx, VBFjet1_JetIdx, VBFjet2_JetIdx, "
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, "
            "FatJet_pt{3}, FatJet_eta, FatJet_phi, FatJet_mass{3}, "
            "MET{5}_pt{4}, MET{5}_phi{4}, "
            "{7}tt_svfit_pt{6}, {7}tt_svfit_eta{6}, {7}tt_svfit_phi{6}, {7}tt_svfit_mass{6})".format(
                self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst,
                self.met_syst, self.met_smear_tag, self.systs, p_sv)))

        for ifeat, feature in enumerate(features):
            df = df.Define(feature, "hhfeatures%s[%s]" % (self.systs, ifeat))
        return df, features


def HH(**kwargs):
    return lambda: HHProducer(**kwargs)


def HHKinFitRDF(*args, **kwargs):
    # The output of HHKinFitRDF changes when targeting H or Z, accornding to the 
    # mass hypothesis for the fit: at the moment there is a bug with the hypothesis 
    # To bypass it, it's necessary to manually change mh1 and mh2 values inside 
    # HHKinFit2/HHKinFit2Scenarios/interface/HHKinFitMasterHeavyHiggs.h
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 1 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHKinFitRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)


def HHVarRDF(*args, **kwargs):
    """ Computes variables relating to H/Z and to HH/ZZ/ZH (pt, eta, mass, phi) 
    For bb, for non-boosted sums the 2 AK4 jets 4-vectors. For boosted uses the AK8 jet directly
    """
    # The output of HHVarRDF is not affected by H or Z, but the output features
    # are called in different ways according to the ZZ or HH analysis
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 2 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHVarRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)
