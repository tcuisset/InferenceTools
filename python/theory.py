""" PDF weights, alpha_s, etc 

Doc:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToPDF
https://indico.cern.ch/event/938672/contributions/3943718/attachments/2073936/3482265/MC_ContactReport_v3.pdf
https://gitlab.cern.ch/cms-analysis/general/columnflow/-/blob/master/columnflow/production/cms/pdf.py
https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/genWeightsTable_cfi.py
"""
import os

from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst
from analysis_tools.utils import import_root

ROOT = import_root()


class PdfWeightsRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.isMC = kwargs["isMC"]
        self.year = kwargs.pop("year")
        self.systematic_is_central = kwargs.pop("systematic_is_central")
        self.force_pdf_to_one = kwargs.pop("force_pdf_to_one", False)
        if self.isMC and self.systematic_is_central:
            if not os.getenv("_Htt_trigSF"):
                os.environ["_Htt_trigSF"] = "_Htt_trigSF"
                
                base = "{}/{}/src/Tools/Tools".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                    ROOT.gSystem.Load("libToolsTools.so")
                ROOT.gROOT.ProcessLine(".L {}/interface/TheoryUncertainties.h".format(base))               


    def run(self, df):
        if not self.isMC:
            return df, []
        
        if not self.systematic_is_central:
            df = df.Define("pdfWeight", "1.f")
            df = df.Define("scaleWeight", "LHEScaleWeight[4]")
            return df, ["pdfWeight", "scaleWeight"]
        
        if self.force_pdf_to_one:
            df = df.Define("pdfWeight", "1.f")
            df = df.Define("pdfWeight_up", "1.f")
            df = df.Define("pdfWeight_down", "1.f")
        else:
            df = df.Define("pdfWeights_values", "computePDFWeights(LHEPdfWeight)")
            df = df.Define("pdfWeight", "pdfWeights_values[0]")
            df = df.Define("pdfWeight_up", "pdfWeights_values[1]")
            df = df.Define("pdfWeight_down", "pdfWeights_values[2]")

        df = df.Define("scaleWeights_values", "computeScaleUncertainties(LHEScaleWeight)")
        df = df.Define("scaleWeight", "scaleWeights_values[0]")
        df = df.Define("scaleWeight_up", "scaleWeights_values[1]")
        df = df.Define("scaleWeight_down", "scaleWeights_values[2]")

        return df, ["pdfWeight", "pdfWeight_up", "pdfWeight_down", "scaleWeight", "scaleWeight_up", "scaleWeight_down"]


def PdfWeightsRDF(**kwargs):
    """ 
    Compute PDF & scale weights with uncertainties
    """
    return lambda: PdfWeightsRDFProducer(**kwargs)
