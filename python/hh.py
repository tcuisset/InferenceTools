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
    def run(self, df):
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
                    int pairType, int dau1_index, int dau2_index, int bjet1_index, int bjet2_index,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                    Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                    float met_pt, float met_phi, float met_covXX, float met_covXY, float met_covYY) {
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

                auto kinFit = HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov);
                kinFit.addHypo(125, 125);
                return kinFit.fit();
            }
        """)

        df = df.Define("hhkinfit_result",
            "compute_hhkinfit(pairType, dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx, "
                "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
                "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
                "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
                "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, "
                "MET{5}_pt{4}, MET{5}_phi{4}, MET_covXX, MET_covXY, MET_covYY)".format(
                    self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst, self.met_syst,
                    self.met_smear_tag)).Define(
            "HHKinFit_mass%s" % self.systs, "hhkinfit_result[0]").Define(
            "HHKinFit_chi2%s" % self.systs, "hhkinfit_result[1]")
        return df, ["HHKinFit_mass%s" % self.systs, "HHKinFit_chi2%s" % self.systs]




def HH(**kwargs):
    return lambda: HHProducer(**kwargs)
