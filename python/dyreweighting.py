import os
from array import array
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from Base.Modules.baseModules import DummyModule
from analysis_tools.utils import import_root


ROOT = import_root()

class DYstitchingProducer(Module):
    def __init__(self, year, *args, **kwargs):
        super(DYstitchingProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/DYstitching.h".format(base))

        self.dy_reweighter = ROOT.DYstitching(year)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('DYstitchWeight', 'F')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        weight = self.dy_reweighter.get_stitching_weight(
            event.LHE_Nb, event.LHE_Njets, event.LHE_HT)
        self.out.fillBranch("DYstitchWeight", weight)
        return True


class DYstitchingRDFProducer():
    def __init__(self, year, isDY, *args, **kwargs):
        self.isDY = isDY
        if self.isDY:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/DYstitching.h".format(base))
            ROOT.gInterpreter.Declare('auto dy_stitcher = DYstitching(%s);' % year)

    def run(self, df):
        if self.isDY:
            df = df.Define("DYstitchWeight",
                "dy_stitcher.get_stitching_weight(LHE_Nb, LHE_Njets, LHE_HT)")
                # "50.")
        else:
            df = df.Define("DYstitchWeight", "1")
        return df, ["DYstitchWeight"]


def DYstitching(**kwargs):
    isDY = kwargs.pop("isDY", False)
    year = kwargs.pop("year")
    if isDY:
        return lambda: DYstitchingProducer(year=year, **kwargs)
    else:
        return lambda: DummyModule(**kwargs)


def DYstitchingRDF(*args, **kwargs):
    isDY = kwargs.pop("isDY", False)
    year = kwargs.pop("year")

    return lambda: DYstitchingRDFProducer(year=year, isDY=isDY, *args, **kwargs)


class DYscalingRDFProducer():
    def __init__(self, year, isDY, *args, **kwargs):
        self.isDY = isDY
        if self.isDY:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/DYscaling.h".format(base))
            ROOT.gInterpreter.Declare('auto dy_scaler = DYscaling(%s);' % year)

    def run(self, df):
        branches = ["DYscale_LL", "DYscale_MM", "DYscale_MH", "DYscale_MTT",
            "DYscale_MTT_up", "DYscale_MTT_down"]
        if self.isDY:
            df = df.Define("dyscalingweights",
                "dy_scaler.get_dy_scale(GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass,"
                    "GenJet_hadronFlavour, LHE_Nb, GenPart_pt, GenPart_eta, GenPart_phi,"
                    " GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")

            for ib, branch in enumerate(branches):
                df = df.Define(branch, "dyscalingweights[%s]" % ib)
        else:
            for branch in branches:
                df = df.Define(branch, "1")
        return df, branches


def DYscalingRDF(*args, **kwargs):
    isDY = kwargs.pop("isDY", False)
    year = kwargs.pop("year")

    return lambda: DYscalingRDFProducer(year=year, isDY=isDY, *args, **kwargs)
