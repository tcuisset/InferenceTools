""" Fix issues of Madgraph Run2 LO genweights having issues
See https://twiki.cern.ch/twiki/bin/view/CMS/MCKnownIssues#Different_weights_at_LHE_level_f
To be used only for Madgraph Run2 LO samples with no interference effects
"""
from analysis_tools.utils import import_root


ROOT = import_root()


class GenWeightLOFixRDFProducer():
    def __init__(self, setGenWeightToOne, isMC):
        self.setGenWeightToOne = setGenWeightToOne
        self.isMC = isMC

    def run(self, df):
        if not self.isMC:
            return df, []
        
        if self.setGenWeightToOne:
            df = df.Define("genWeightFixed", "1")
        else:
            df = df.Define("genWeightFixed", "genWeight")

        return df, ["genWeightFixed"]

def GenWeightLOFixRDF(*args, **kwargs):
    setGenWeightToOne = kwargs.pop("setGenWeightToOne", False)
    isMC = kwargs.pop("isMC")

    return lambda: GenWeightLOFixRDFProducer(setGenWeightToOne=setGenWeightToOne, isMC=isMC)
