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
        self.algo = kwargs.pop("algo", "ClassicSVFit")
        if self.algo == "ClassicSVFit":
            self.svfitFctName = "FitAndGetResultWithInputs"
        elif self.algo == "FastMTT":
            self.svfitFctName = "FitAndGetResultWithInputs_FastMTT"
        else:
            raise ValueError()
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
            f"svfit.{self.svfitFctName}(0, pairType, dau1_decayMode, dau2_decayMode, "
                f"dau1_pt{self.lep_syst}, dau1_eta, dau1_phi, dau1_mass{self.lep_syst}, "
                f"dau2_pt{self.lep_syst}, dau2_eta, dau2_phi, dau2_mass{self.lep_syst},"
                f"MET{self.met_smear_tag}_pt{self.met_syst}, MET{self.met_smear_tag}_phi{self.met_syst}, MET_covXX, MET_covXY, MET_covYY)").Define(
            "%stt_svfit_pt%s" % (p, self.systs), "svfit_result%s[0]" % self.systs).Define(
            "%stt_svfit_eta%s" % (p, self.systs), "svfit_result%s[1]" % self.systs).Define(
            "%stt_svfit_phi%s" % (p, self.systs), "svfit_result%s[2]" % self.systs).Define(
            "%stt_svfit_mass%s" % (p, self.systs), "svfit_result%s[3]" % self.systs)
        return df, branches


def SVFit(**kwargs):
    return lambda: SVFitProducer(**kwargs)


def SVFitRDF(*args, **kwargs):
    """ Compute SVFit or FastMTT
    Parameters : 
     - AnalysisType (not used anymore)
     - algo : can be SVFit or FastMTT. Make sure that the TauAnalysis/ClassicSVfit has been installed with the correct branch, and that Tools/Tools has been built with the FASTMTT_SUPPORT CppDefine for FastMTT
     """
    # The output of SVFit is not affected by H or Z, but the output features
    # are called in different ways according to the analysis
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG : AnalysisType = {}".format(AnalysisType))
    return lambda: SVFitRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)
