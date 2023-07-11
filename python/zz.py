import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst
from analysis_tools.utils import import_root

ROOT = import_root()


class ZZProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(ZZProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/HHKinFitInterface.h".format(base))
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("Zbb_pt%s" % self.systs, "F")
        self.out.branch("Zbb_eta%s" % self.systs, "F")
        self.out.branch("Zbb_phi%s" % self.systs, "F")
        self.out.branch("Zbb_mass%s" % self.systs, "F")

        self.out.branch("Ztt_pt%s" % self.systs, "F")
        self.out.branch("Ztt_eta%s" % self.systs, "F")
        self.out.branch("Ztt_phi%s" % self.systs, "F")
        self.out.branch("Ztt_mass%s" % self.systs, "F")

        self.out.branch("ZZ_pt%s" % self.systs, "F")
        self.out.branch("ZZ_eta%s" % self.systs, "F")
        self.out.branch("ZZ_phi%s" % self.systs, "F")
        self.out.branch("ZZ_mass%s" % self.systs, "F")

        self.out.branch("ZZ_svfit_pt%s" % self.systs, "F")
        self.out.branch("ZZ_svfit_eta%s" % self.systs, "F")
        self.out.branch("ZZ_svfit_phi%s" % self.systs, "F")
        self.out.branch("ZZ_svfit_mass%s" % self.systs, "F")

        self.out.branch("ZZKinFit_mass%s" % self.systs, "F")
        self.out.branch("ZZKinFit_chi2%s" % self.systs, "F")

        self.out.branch("VBFjj_mass%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaEta%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaPhi%s" % self.systs, "F")
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ztt_svfit = Object(event, "Ztt_svfit")
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
        ztt_tlv = dau1_tlv + dau2_tlv
        zbb_tlv = bjet1_tlv + bjet2_tlv
        zz_tlv = ztt_tlv + zbb_tlv

        # HH (with ztt_svfit) TLorentzVector
        ztt_svfit_tlv = ROOT.TLorentzVector()
        ztt_svfit_tlv.SetPtEtaPhiM(
            eval("ztt_svfit.pt%s" % self.systs), 
            eval("ztt_svfit.eta%s" % self.systs), 
            eval("ztt_svfit.phi%s" % self.systs), 
            eval("ztt_svfit.mass%s" % self.systs),
        )
        zz_svfit_tlv = zbb_tlv + ztt_svfit_tlv

        # HH Kin. Fit
        # constructor https://github.com/bvormwald/HHKinFit2/blob/CMSSWversion/HHKinFit2Scenarios/src/HHKinFitMasterHeavyHiggs.cpp#L27
        # kinFit = ROOT.HHKinFit2.HHKinFitMasterHeavyHiggs(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit = ROOT.HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit.addHypo(91, 91);
        results = kinFit.fit()
        ZZKinFit_mass = results[0]
        ZZKinFit_chi2 = results[1]

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

        self.out.fillBranch("Zbb_pt%s" % self.systs, zbb_tlv.Pt())
        self.out.fillBranch("Zbb_eta%s" % self.systs, zbb_tlv.Eta())
        self.out.fillBranch("Zbb_phi%s" % self.systs, zbb_tlv.Phi())
        self.out.fillBranch("Zbb_mass%s" % self.systs, zbb_tlv.M())

        self.out.fillBranch("Ztt_pt%s" % self.systs, ztt_tlv.Pt())
        self.out.fillBranch("Ztt_eta%s" % self.systs, ztt_tlv.Eta())
        self.out.fillBranch("Ztt_phi%s" % self.systs, ztt_tlv.Phi())
        self.out.fillBranch("Ztt_mass%s" % self.systs, ztt_tlv.M())

        self.out.fillBranch("ZZ_pt%s" % self.systs, zz_tlv.Pt())
        self.out.fillBranch("ZZ_eta%s" % self.systs, zz_tlv.Eta())
        self.out.fillBranch("ZZ_phi%s" % self.systs, zz_tlv.Phi())
        self.out.fillBranch("ZZ_mass%s" % self.systs, zz_tlv.M())

        self.out.fillBranch("ZZ_svfit_pt%s" % self.systs, zz_svfit_tlv.Pt())
        self.out.fillBranch("ZZ_svfit_eta%s" % self.systs, zz_svfit_tlv.Eta())
        self.out.fillBranch("ZZ_svfit_phi%s" % self.systs, zz_svfit_tlv.Phi())
        self.out.fillBranch("ZZ_svfit_mass%s" % self.systs, zz_svfit_tlv.M())

        self.out.fillBranch("ZZKinFit_mass%s" % self.systs, ZZKinFit_mass)
        self.out.fillBranch("ZZKinFit_chi2%s" % self.systs, ZZKinFit_chi2)

        self.out.fillBranch("VBFjj_mass%s" % self.systs, VBFjj_mass)
        self.out.fillBranch("VBFjj_deltaEta%s" % self.systs, VBFjj_deltaEta)
        self.out.fillBranch("VBFjj_deltaPhi%s" % self.systs, VBFjj_deltaPhi)

        return True


class ZZKinFitRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        self.isMC = kwargs.pop("isMC")
        super(ZZKinFitRDFProducer, self).__init__(isMC=self.isMC, *args, **kwargs)

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
                        float met_covXX, float met_covXY, float met_covYY) {
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
                    kinFit.addHypo(91, 91);
                    return kinFit.fit();
                }
            """)

    def run(self, df):
        isMC = ("true" if self.isMC else "false")
        jet_resolution = "jet_pt_resolution"
        if not self.isMC:
            jet_resolution = "Jet_eta"  # placeholder
        branches = ["ZZKinFit_mass%s" % self.systs, "ZZKinFit_chi2%s" % self.systs]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        # print("### DEBUG: Running kin fitter")
        df = df.Define("zzkinfit_result%s" % self.systs,
            "compute_hhkinfit({6}, pairType, dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx, "
                "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
                "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
                "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
                "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, {7},"
                "MET{5}_pt{4}, MET{5}_phi{4}, MET_covXX, MET_covXY, MET_covYY)".format(
                    self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst, self.met_syst,
                    self.met_smear_tag, isMC, jet_resolution)
            ).Define("ZZKinFit_mass%s" % self.systs, "zzkinfit_result%s[0]" % self.systs
            ).Define("ZZKinFit_chi2%s" % self.systs, "zzkinfit_result%s[1]" % self.systs)
        return df, branches


class ZZVarRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(ZZVarRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("_HHVAR_%s" % self.systs):
            os.environ["_HHVAR_%s" % self.systs] = "hhvar"
            ROOT.gInterpreter.Declare("""
                #include <TLorentzVector.h>
                using Vfloat = const ROOT::RVec<float>&;
                ROOT::RVec<double> get_zz_features%s(int pairType,
                        int dau1_index, int dau2_index, int bjet1_index, int bjet2_index,
                        int vbfjet1_index, int vbfjet2_index,
                        Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                        Vfloat electron_pt, Vfloat electron_eta,
                        Vfloat electron_phi, Vfloat electron_mass,
                        Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                        Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                        float met_pt, float met_phi,
                        double Ztt_svfit_pt, double Ztt_svfit_eta,
                        double Ztt_svfit_phi, double Ztt_svfit_mass) {

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
                    auto ztt_tlv = dau1_tlv + dau2_tlv;
                    bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                        jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                    bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                        jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));
                    auto zbb_tlv = bjet1_tlv + bjet2_tlv;
                    auto zz_tlv = ztt_tlv + zbb_tlv;

                    auto met_tlv = TLorentzVector();
                    met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

                    auto ztt_met_tlv = ztt_tlv + met_tlv;

                    double zz_svfit_pt = -999., zz_svfit_eta = -999.;
                    double zz_svfit_phi = -999., zz_svfit_mass = -999.;
                    if (Ztt_svfit_pt > 0) {
                        auto ztt_svfit_tlv = TLorentzVector();
                        ztt_svfit_tlv.SetPtEtaPhiM(Ztt_svfit_pt, Ztt_svfit_eta,
                            Ztt_svfit_phi, Ztt_svfit_mass);
                        auto zz_svfit_tlv = ztt_svfit_tlv + zbb_tlv;
                        zz_svfit_pt = zz_svfit_tlv.Pt();
                        zz_svfit_eta = zz_svfit_tlv.Eta();
                        zz_svfit_phi = zz_svfit_tlv.Phi();
                        zz_svfit_mass = zz_svfit_tlv.M();
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
                        zbb_tlv.Pt(), zbb_tlv.Eta(), zbb_tlv.Phi(), zbb_tlv.M(),
                        ztt_tlv.Pt(), ztt_tlv.Eta(), ztt_tlv.Phi(), ztt_tlv.M(),
                        ztt_met_tlv.Pt(), ztt_met_tlv.Eta(), ztt_met_tlv.Phi(), ztt_met_tlv.M(),
                        zz_tlv.Pt(), zz_tlv.Eta(), zz_tlv.Phi(), zz_tlv.M(),
                        zz_svfit_pt, zz_svfit_eta, zz_svfit_phi, zz_svfit_mass,
                        vbfjj_mass, vbfjj_deltaEta, vbfjj_deltaPhi
                    };
                }
            """ % self.systs)

    def run(self, df):
        features = ("Zbb_pt{0},Zbb_eta{0},Zbb_phi{0},Zbb_mass{0},"
            "Ztt_pt{0},Ztt_eta{0},Ztt_phi{0},Ztt_mass{0},"
            "Ztt_met_pt{0},Ztt_met_eta{0},Ztt_met_phi{0},Ztt_met_mass{0},"
            "ZZ_pt{0},ZZ_eta{0},ZZ_phi{0},ZZ_mass{0},"
            "ZZ_svfit_pt{0},ZZ_svfit_eta{0},ZZ_svfit_phi{0},ZZ_svfit_mass{0},"
            "ZZ_VBFjj_mass{0},ZZ_VBFjj_deltaEta{0},ZZ_VBFjj_deltaPhi{0}".format(self.systs))
        features = list(features.split(","))
        all_branches = df.GetColumnNames()
        if features[0] in all_branches:
            return df, []

        df = df.Define("zzfeatures%s" % self.systs, ("get_zz_features{6}(pairType, "
            "dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx, VBFjet1_JetIdx, VBFjet2_JetIdx, "
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, "
            "MET{5}_pt{4}, MET{5}_phi{4}, "
            "Ztt_svfit_pt{6}, Ztt_svfit_eta{6}, Ztt_svfit_phi{6}, Ztt_svfit_mass{6})".format(
                self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst,
                self.met_syst, self.met_smear_tag, self.systs)))

        for ifeat, feature in enumerate(features):
            df = df.Define(feature, "zzfeatures%s[%s]" % (self.systs, ifeat))
        return df, features


def ZZ(**kwargs):
    return lambda: ZZProducer(**kwargs)


def ZZKinFitRDF(*args, **kwargs):
    return lambda: ZZKinFitRDFProducer(*args, **kwargs)


def ZZVarRDF(*args, **kwargs):
    return lambda: ZZVarRDFProducer(*args, **kwargs)
