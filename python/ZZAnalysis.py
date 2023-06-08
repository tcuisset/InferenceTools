import os

from analysis_tools.utils import import_root

ROOT = import_root()

class ZZEllipticalCutFilterRDFProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            bool apply_elliptical_cut(float Htt_svfit_mass, float Hbb_mass) {
                float c_Htt_svfit_mass = 99.0;
                float r_Htt_svfit_mass = 39.0;
                float c_Hbb_mass = 101.0;
                float r_Hbb_mass = 101.0;
                float dist_Htt_svfit_mass = (Htt_svfit_mass - c_Htt_svfit_mass) * (Htt_svfit_mass - c_Htt_svfit_mass)/(r_Htt_svfit_mass * r_Htt_svfit_mass);
                float dist_Hbb_mass = (Hbb_mass - c_Hbb_mass) * (Hbb_mass - c_Hbb_mass)/(r_Hbb_mass * r_Hbb_mass);
                return (dist_Htt_svfit_mass + dist_Hbb_mass) < 1;
            }
        """)

    def run(self, df):
        # define a new branch to check if the event is inside or outside the ellipse
        df = df.Define("isInside", """apply_elliptical_cut(
            Htt_svfit_mass,
            Hbb_mass
        )""")
        # filter the events outside the ellipse
        df = df.Filter("isInside == 1")
        return df, []

def ZZEllipticalCutFilterRDF(*args, **kwargs):
    return lambda: ZZEllipticalCutFilterRDFProducer(*args, **kwargs)

class ZZBBTauTauFilterRDFProducer():
    def __init__(self, isZZsig, isZZbkg, *args, **kwargs):
        self.isZZsig = isZZsig
        self.isZZbkg = isZZbkg
        # print(" ### DEBUG: isZZsig = {}".format(isZZsig))
        # print(" ### DEBUG: isZZbkg = {}".format(isZZbkg))

        ROOT.gInterpreter.Declare("""
            using Vfloat = const ROOT::RVec<float>&;
            using Vint   = const ROOT::RVec<int>&;
            bool find_bb_tautau(Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                bool FoundSignal = false;
                int n_b_fromZ = 0;
                int n_tau_fromZ = 0;
                for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                    if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                    if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                        // std::cout << i_gen << std::endl;
                        n_b_fromZ += 1;
                    }
                    if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                        // std::cout << i_gen << std::endl;
                        n_tau_fromZ += 1;
                    }
                }
                if ((n_b_fromZ > 1) && (n_tau_fromZ > 1)) {
                    FoundSignal = true;
                }
                return FoundSignal;
            }
        """)

    def run(self, df):
        if self.isZZsig or self.isZZbkg:
            # define a new branch to check if it is or not a ZZ->bbtautau event
            df = df.Define("isZZTobbtautau", """find_bb_tautau(
                GenPart_pdgId,
                GenPart_genPartIdxMother
            )""")
            # filter the events with ZZ->bbtautau
            if self.isZZsig:
                # print(" ### DEBUG: isZZTobbtautau == 1")
                df = df.Filter("isZZTobbtautau == 1")
            # filter the events without ZZ->bbtautau
            elif self.isZZbkg:
                # print(" ### DEBUG: isZZTobbtautau == 0")
                df = df.Filter("isZZTobbtautau == 0")
        return df, []

def ZZBBTauTauFilterRDF(*args, **kwargs):
    isZZsig = kwargs.pop("isZZsig")
    isZZbkg = kwargs.pop("isZZbkg")

    # print("### DEBUG : isZZsig = {}, isZZbkg = {}".format(isZZsig, isZZbkg))
    return lambda: ZZBBTauTauFilterRDFProducer(isZZsig=isZZsig, isZZbkg=isZZbkg, *args, **kwargs)