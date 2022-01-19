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
        ROOT.gROOT.ProcessLine(".L {}/interface/DYreweighting.h".format(base))

        self.dy_reweighter = ROOT.DYreweighting(year)

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
            ROOT.gROOT.ProcessLine(".L {}/interface/DYreweighting.h".format(base))
            ROOT.gInterpreter.Declare('auto dy_reweighter = DYreweighting(%s);' % year)

    def run(self, df):
        if self.isDY:
            df = df.Define("DYstitchWeight",
                "dy_reweighter.get_stitching_weight(LHE_Nb, LHE_Njets, LHE_HT)")
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
