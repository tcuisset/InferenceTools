from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from PhysicsTools.NanoAODTools.postprocessing.modules.common.tauCorrProducer import (
    TauCorrectionsProducer
)
from Base.Modules.baseModules import DummyModule


class PrescaleWeightProducer(Module):
    def __init__(self, isMC, year, *args, **kwargs):
        super(PrescaleWeightProducer, self).__init__(*args, **kwargs)
        self.year = year

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("prescaleWeight", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if (self.year == 2017
                and event.HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg == 1):
            self.out.fillBranch("prescaleWeight", 0.65308574)
        elif (self.year == 2018
                and event.HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1 == 1):
            self.out.fillBranch("prescaleWeight", 0.990342)
        else:
            self.out.fillBranch("prescaleWeight", 1.)
        return True

def prescaleWeight(**kwargs):
    isMC = kwargs.pop("isMC")
    year = kwargs.pop("year")
    if isMC:
        return lambda: PrescaleWeightProducer(isMC, year, **kwargs)
    else:
        return lambda: DummyModule(**kwargs)


class PrescaleWeightRDFProducer():
    def __init__(self, isMC, year, *args, **kwargs):
        self.year = year
        self.isMC = isMC

    def run(self, df):
        if not self.isMC:
            return df, []
        
        if self.year == 2017:
            df = df.Define("prescaleWeight",
                "0.65308574 * (HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg == 1)"
                " + 1. * (HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg == 0)")
        elif self.year == 2018:
            df = df.Define("prescaleWeight",
                "0.990342 * (HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1 == 1)"
                " + 1. * (HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1 == 0)")
        else:
            df = df.Define("prescaleWeight", "1.")
        return df, ["prescaleWeight"]


def prescaleWeightRDF(**kwargs):
    isMC = kwargs.pop("isMC")
    year = int(kwargs.pop("year"))
    return lambda: PrescaleWeightRDFProducer(isMC, year, **kwargs)


def tauCorrections(**kwargs):
    isMC = kwargs.pop("isMC")
    year = int(kwargs.pop("year"))
    if not isMC:
        return lambda: DummyModule(**kwargs)
    if year == 2016:
        return lambda: TauCorrectionsProducer('2016Legacy', **kwargs)
    elif year == 2017:
        return lambda: TauCorrectionsProducer('2017ReReco', **kwargs)
    elif year == 2018:
        return lambda: TauCorrectionsProducer('2018ReReco', **kwargs)
    