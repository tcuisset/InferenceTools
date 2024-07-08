import os
from array import array
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from Base.Modules.baseModules import JetLepMetSyst, JetLepMetModule
from analysis_tools.utils import import_root


ROOT = import_root()

class SVFitProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(SVFitProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/SVfitinterface.h".format(base))

        self.svfit = ROOT.SVfitinterface()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('Htt_svfit_pt%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_eta%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_phi%s' % self.systs, 'F')
        self.out.branch('Htt_svfit_mass%s' % self.systs, 'F')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        met, met_tlv = self.get_met(event)

        decayMode1 = (dau1.decayMode if dau1 in taus else -1)
        decayMode2 = (dau2.decayMode if dau2 in taus else -1)

        result = self.svfit.FitAndGetResultWithInputs(
            0, event.pairType, decayMode1, decayMode2,
            dau1_tlv.Pt(), dau1_tlv.Eta(), dau1_tlv.Phi(), dau1_tlv.M(),
            dau2_tlv.Pt(), dau2_tlv.Eta(), dau2_tlv.Phi(), dau2_tlv.M(),
            met_tlv.Pt(), met_tlv.Phi(), met.covXX, met.covXY, met.covYY
        )

        self.out.fillBranch("Htt_svfit_pt%s" % self.systs, result[0])
        self.out.fillBranch("Htt_svfit_eta%s" % self.systs, result[1])
        self.out.fillBranch("Htt_svfit_phi%s" % self.systs, result[2])
        self.out.fillBranch("Htt_svfit_mass%s" % self.systs, result[3])
        return True


class SVFitRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.AnalysisType = AnalysisType
        super(SVFitRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("_SVFIT"):
            os.environ["_SVFIT"] = "svfit"

            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")
            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/SVfitinterface.h".format(base))

            ROOT.gInterpreter.Declare("""
                auto svfit = SVfitinterface();
            """)

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint = const ROOT::RVec<int>&;
                ROOT::RVec<double> compute_svfit(
                        int pairType, int dau1_index, int dau2_index,
                        Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                        Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                        Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass, Vint Tau_decayMode,
                        float met_pt, float met_phi, float met_covXX, float met_covXY, float met_covYY) {
                    double dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;
                    int DM1=-1, DM2=-1;
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
                        DM1 = Tau_decayMode.at(dau1_index);
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_eta = tau_eta.at(dau2_index);
                    dau2_phi = tau_phi.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);
                    DM2 = Tau_decayMode.at(dau2_index);

                    return svfit.FitAndGetResultWithInputs(0, pairType, DM1, DM2,
                        dau1_pt, dau1_eta, dau1_phi, dau1_mass,
                        dau2_pt, dau2_eta, dau2_phi, dau2_mass,
                        met_pt, met_phi, met_covXX, met_covXY, met_covYY);
                }
            """)

    def run(self, df):

        # DEBUG
        if not self.AnalysisType:
            p = "H"
            print(" ### INFO: Running SVFit with default option for HH analysis")
        else:
            p = "X"
            # print(" ### INFO: Running SVFit with AnalysisType = {}".format(self.AnalysisType))
            # if self.AnalysisType == "Zbb_Ztautau" or self.AnalysisType == "Ztautau_Hbb":    p = "Z"
            # elif self.AnalysisType == "Zbb_Htautau":                                        p = "H"

        branches = ["%stt_svfit_pt%s"   % (p, self.systs), 
                    "%stt_svfit_eta%s"  % (p, self.systs),
                    "%stt_svfit_phi%s"  % (p, self.systs), 
                    "%stt_svfit_mass%s" % (p, self.systs)]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        df = df.Define("svfit_result%s" % self.systs,
            "compute_svfit(pairType, dau1_index, dau2_index, "
                f"Muon_pt{self.muon_syst}, Muon_eta, Muon_phi, Muon_mass{self.muon_syst}, "
                f"Electron_pt{self.electron_syst}, Electron_eta, Electron_phi, Electron_mass{self.electron_syst}, "
                f"Tau_pt{self.tau_syst}, Tau_eta, Tau_phi, Tau_mass{self.tau_syst}, Tau_decayMode, "
                f"MET{self.met_smear_tag}_pt{self.met_syst}, MET{self.met_smear_tag}_phi{self.met_syst}, MET_covXX, MET_covXY, MET_covYY)").Define(
            "%stt_svfit_pt%s" % (p, self.systs), "svfit_result%s[0]" % self.systs).Define(
            "%stt_svfit_eta%s" % (p, self.systs), "svfit_result%s[1]" % self.systs).Define(
            "%stt_svfit_phi%s" % (p, self.systs), "svfit_result%s[2]" % self.systs).Define(
            "%stt_svfit_mass%s" % (p, self.systs), "svfit_result%s[3]" % self.systs)
        return df, branches


def SVFit(**kwargs):
    return lambda: SVFitProducer(**kwargs)


def SVFitRDF(*args, **kwargs):
    # The output of SVFit is not affected by H or Z, but the output features
    # are called in different ways according to the analysis
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG : AnalysisType = {}".format(AnalysisType))
    return lambda: SVFitRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)
