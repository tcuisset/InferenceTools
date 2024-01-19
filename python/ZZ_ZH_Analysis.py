import os

from analysis_tools.utils import import_root

ROOT = import_root()

class ZZEllipticalCutFilterRDFProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            bool apply_elliptical_cut_80(float Ztt_svfit_mass, float Zbb_mass) {
                float c_Ztt_svfit_mass = 105.0;
                float r_Ztt_svfit_mass = 51.0;
                float c_Zbb_mass = 118.0;
                float r_Zbb_mass = 113.0;
                float dist_Ztt_svfit_mass = (Ztt_svfit_mass - c_Ztt_svfit_mass) * (Ztt_svfit_mass - c_Ztt_svfit_mass)/(r_Ztt_svfit_mass * r_Ztt_svfit_mass);
                float dist_Zbb_mass = (Zbb_mass - c_Zbb_mass) * (Zbb_mass - c_Zbb_mass)/(r_Zbb_mass * r_Zbb_mass);
                return (dist_Ztt_svfit_mass + dist_Zbb_mass) < 1;
            }
        """)

    def run(self, df):
        # define a new branch to check if the event is inside or outside the ellipse
        df = df.Define("isInside", """apply_elliptical_cut_80(
            Ztt_svfit_mass,
            Zbb_mass
        )""")
        # filter the events outside the ellipse
        df = df.Filter("isInside == 1", "ZZEllipticalCutFilterRDF")
        return df, []

def ZZEllipticalCutFilterRDF(*args, **kwargs):
    return lambda: ZZEllipticalCutFilterRDFProducer(*args, **kwargs)

class BBTauTauFilterRDFProducer():
    def __init__(self, ProcType, isSigBBTT, isBkgBBTT, *args, **kwargs):
        self.isSigBBTT = isSigBBTT
        self.isBkgBBTT = isBkgBBTT
        self.ProcType = ProcType
        # print(" ### DEBUG: isSigBBTT = {}".format(isSigBBTT))
        # print(" ### DEBUG: isBkgBBTT = {}".format(isBkgBBTT))

        if not os.getenv("_BBTauTauFilter"):
            os.environ["_BBTauTauFilter"] = "_BBTauTauFilter"

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint   = const ROOT::RVec<int>&;
                bool find_Zbb_Ztautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                    bool FoundSignal = false;
                    int n_b_fromZ = 0;
                    int n_tau_fromZ = 0;
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_tau_fromZ += 1;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromZ > 1)) {
                        FoundSignal = true;
                    }
                    return FoundSignal;
                }
                bool find_Zbb_Htautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                    bool FoundSignal = false;
                    int n_b_fromZ = 0;
                    int n_tau_fromH = 0;
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_tau_fromH += 1;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromH > 1)) {
                        FoundSignal = true;
                    }
                    return FoundSignal;
                }
                bool find_Ztautau_Hbb (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                    bool FoundSignal = false;
                    int n_b_fromH = 0;
                    int n_tau_fromZ = 0;
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_b_fromH += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_tau_fromZ += 1;
                        }
                    }
                    if ((n_b_fromH > 1) && (n_tau_fromZ > 1)) {
                        FoundSignal = true;
                    }
                    return FoundSignal;
                }
            """)

    def run(self, df):
        if self.isSigBBTT or self.isBkgBBTT:

            # define a new branch to check if it is or not a ZZ/ZH->bbtautau event
            if self.ProcType == "Zbb_Ztautau":
                print(" ### Running bbtautau filter for Zbb_Ztautau")
                df = df.Define("isBBTT", """find_Zbb_Ztautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            elif self.ProcType == "Zbb_Htautau":
                print(" ### Running bbtautau filter for Zbb_Htautau")
                df = df.Define("isBBTT", """find_Zbb_Htautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            elif self.ProcType == "Ztautau_Hbb":
                print(" ### Running bbtautau filter for Ztautau_Hbb")
                df = df.Define("isBBTT", """find_Ztautau_Hbb(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            else:
                raise ValueError("BBTauTauFilterRDF not implemented for self.ProcType = ", self.ProcType)
            
            # filter the events with ZZ/ZH->bbtautau
            if self.isSigBBTT:
                # print(" ### DEBUG: isBBTT == 1")
                df = df.Filter("isBBTT == 1", "BBTauTauFilterRDF")
            # filter the events without ZZ/ZH->bbtautau
            elif self.isBkgBBTT:
                # print(" ### DEBUG: isBBTT == 0")
                df = df.Filter("isBBTT == 0", "BBTauTauFilterRDF")
                
        return df, []

class BBTauTauFilterDummyRDFProducer():
    def run(self, df):
        return df, []

def BBTauTauFilterRDF(*args, **kwargs):

    ProcType = kwargs.pop("ProcType")
    isSigBBTT = kwargs.pop("isSigBBTT")
    isBkgBBTT = kwargs.pop("isBkgBBTT")

    if isSigBBTT or isBkgBBTT:
        return lambda: BBTauTauFilterRDFProducer(ProcType=ProcType, isSigBBTT=isSigBBTT, isBkgBBTT=isBkgBBTT, *args, **kwargs)
    else:
        return lambda: BBTauTauFilterDummyRDFProducer()
