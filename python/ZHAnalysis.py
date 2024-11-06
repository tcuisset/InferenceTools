# OUTDATED
# Tools (genFilters and elliptical cut) for ZH->bbtautau analysis (both non-resonant and resonant)

import os

from analysis_tools.utils import import_root

ROOT = import_root()


class ZHBBTauTauFilterRDFProducer():
    def __init__(self, isZHsig:list[str], isZHbkg:bool, *args, **kwargs):
        self.isZHsig = isZHsig
        self.isZHbkg = isZHbkg
        for val in isZHsig:
            assert (val in ["zbb_htt", "ztt_hbb"]), "isZHsignal should hold zbb_htt or ztt_hbb"
        # print(" ### DEBUG: isZZsig = {}".format(isZZsig))
        # print(" ### DEBUG: isZZbkg = {}".format(isZZbkg))

        if not os.getenv("_ZHBBTauTauFilter"):
            os.environ["_ZHBBTauTauFilter"] = "_ZHBBTauTauFilter"
            # PDG ids : 5=b, 15=tau-, 23=Z, 25=Higgs
            # 0 : not found signal, 1 : found Z->bb,H->tautau, 2 : found Z->tautau,H->bb
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint   = const ROOT::RVec<int>&;
                int find_zh_bb_tautau(Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {
                    int FoundSignal = 0;
                    int n_b_fromZ = 0;
                    int n_b_fromH = 0;
                    int n_tau_fromZ = 0;
                    int n_tau_fromH = 0;
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if (fabs(GenPart_pdgId.at(i_gen)) == 5) { // we have a b
                            if (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23) {
                                n_b_fromZ += 1;
                            }
                            else if (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25) {
                                n_b_fromH += 1;   
                            }
                        }
                        else if (fabs(GenPart_pdgId.at(i_gen)) == 15) { // we have a tau
                            if (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23) {
                                n_tau_fromZ += 1;
                            }
                            else if (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25) {
                                n_tau_fromH += 1;   
                            }
                        }
                    }
                    
                    if ((n_b_fromZ > 1) && (n_tau_fromH > 1)) {
                        FoundSignal = 1;
                    } else if ((n_b_fromH > 1) && (n_tau_fromZ > 1)) {
                        FoundSignal = 2;
                    }
                    return FoundSignal;
                }
            """)

    def run(self, df):
        if self.isZHsig or self.isZHbkg:
            # define a new branch to check if it is or not a ZH->bbtautau event
            df = df.Define("isZHTobbtautau", """find_zh_bb_tautau(
                GenPart_pdgId,
                GenPart_genPartIdxMother
            )""")
            foundSignalIds = []
            if "zbb_htt" in self.isZHsig:
                foundSignalIds.append(1)
            if "ztt_hbb" in self.isZHsig:
                foundSignalIds.append(2)
            if self.isZHbkg:
                assert len(foundSignalIds) == 0
                foundSignalIds.append(0)
            assert len(foundSignalIds) > 0
            df = df.Filter("||".join(f"(isZHTobbtautau == {foundSignalId})" for foundSignalId in foundSignalIds), "ZHBBTauTauFilterRDF_Sig")
        
        return df, []

class ZHBBTauTauFilterDummyRDFProducer():
    def run(self, df):
        return df, []

def ZHBBTauTauFilterRDF(*args, **kwargs):
    isZHsig = kwargs.pop("isZHsig")
    isZHsig = [isZHsig] if isinstance(isZHsig, str) else isZHsig # convert string to list if not already a list
    isZHbkg = kwargs.pop("isZHbkg")

    if isZHsig or isZHbkg:
        return lambda: ZHBBTauTauFilterRDFProducer(isZHsig=isZHsig, isZHbkg=isZHbkg, *args, **kwargs)
    else:
        return lambda: ZHBBTauTauFilterDummyRDFProducer()
