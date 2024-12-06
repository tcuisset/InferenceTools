import os

from analysis_tools.utils import import_root
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()

class FatJetParticleNetSFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)

        if self.isMC:
            self.sampleType = kwargs.pop("fatjet_bb_type")
            if self.sampleType is None and self.isMC:
                print("WARNING : FatJetParticleNetSFProducer sampleType set to None for dataset " + kwargs["dataset"])

            if not os.getenv("_FatJetParticleNetSFProducer"):
                os.environ["_FatJetParticleNetSFProducer"] = "_FatJetParticleNetSFProducer"
                if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                    ROOT.gSystem.Load("libToolsTools.so")

                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/PNetSFInterface.h".format(base))

                ROOT.gInterpreter.Declare(f"""
                    auto PNetAK8SF = PNetSFInterface("{self.year}");
                """)

    def run(self, df):
        #fatjet_bb_tagging_branch = "fRVec(FatJet_particleNetLegacy_Xbb)/(fRVec(FatJet_particleNetLegacy_Xbb)+fRVec(FatJet_particleNetLegacy_QCD))"
        if not self.isMC:
            return df, []
        
        df = df.Define("fatjet_pNet_SF_vec", f""" 
            PNetAK8SF.getSFvec(FatJet_pt{self.jet_syst}, fatjet_JetIdx, genAk8_Zbb_matches, genAk8_Hbb_matches, "{self.sampleType}") 
        """)
        df = df.Define("fatjet_pNet_HP_SF",      "fatjet_pNet_SF_vec[0]")
        df = df.Define("fatjet_pNet_HP_SF_up",   "fatjet_pNet_SF_vec[1]")
        df = df.Define("fatjet_pNet_HP_SF_down", "fatjet_pNet_SF_vec[2]")
        df = df.Define("fatjet_pNet_MP_SF",      "fatjet_pNet_SF_vec[3]")
        df = df.Define("fatjet_pNet_MP_SF_up",   "fatjet_pNet_SF_vec[4]")
        df = df.Define("fatjet_pNet_MP_SF_down", "fatjet_pNet_SF_vec[5]")
        df = df.Define("fatjet_pNet_LP_SF",      "fatjet_pNet_SF_vec[6]")
        df = df.Define("fatjet_pNet_LP_SF_up",   "fatjet_pNet_SF_vec[7]")
        df = df.Define("fatjet_pNet_LP_SF_down", "fatjet_pNet_SF_vec[8]")


        return df, ["fatjet_pNet_HP_SF", "fatjet_pNet_HP_SF_up", "fatjet_pNet_HP_SF_down",
                    "fatjet_pNet_MP_SF", "fatjet_pNet_MP_SF_up", "fatjet_pNet_MP_SF_down",
                    "fatjet_pNet_LP_SF", "fatjet_pNet_LP_SF_up", "fatjet_pNet_LP_SF_down"]


def FatJetParticleNetSFProducerRDF(**kwargs):
    """
    Compute ParticleNet SFs for FatJet X->bb tagging
    Parameters : 
     - dataset : name of dataset, to determine what kind of SF to apply
    Output :
     - fatjet_pNet_HP/MP/LP_SF (_up/_down)
    """
    return lambda: FatJetParticleNetSFProducer(**kwargs)
