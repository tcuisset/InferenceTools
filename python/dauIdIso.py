from analysis_tools.utils import import_root

ROOT = import_root()

class dauIdIsoSFRDFProducer():
    def __init__(self, *args, **kwargs):
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                double get_dauIdIso_sf(
                        int pairType, int dau1_index, int dau2_index, 
                        Vfloat musf_id, Vfloat musf_reliso, Vfloat elesf, 
                        Vfloat Tau_sfDeepTau2017v2p1VSjet, Vfloat Tau_sfDeepTau2017v2p1VSe, 
                        Vfloat Tau_sfDeepTau2017v2p1VSmu) {
                    double idAndIsoSF_leg1 = 1.;
                    if (pairType == 0) {
                        idAndIsoSF_leg1 = musf_id.at(dau1_index) * musf_reliso.at(dau1_index);
                    } else if (pairType == 1) {
                        idAndIsoSF_leg1 = elesf.at(dau1_index);
                    } else if (pairType == 2) {
                        idAndIsoSF_leg1 = Tau_sfDeepTau2017v2p1VSjet.at(dau1_index) * 
                            Tau_sfDeepTau2017v2p1VSe.at(dau1_index) * 
                            Tau_sfDeepTau2017v2p1VSmu.at(dau1_index);
                    }
                     double idAndIsoSF_leg2 = Tau_sfDeepTau2017v2p1VSjet.at(dau2_index) * 
                        Tau_sfDeepTau2017v2p1VSe.at(dau2_index) * 
                        Tau_sfDeepTau2017v2p1VSmu.at(dau2_index);
                    return idAndIsoSF_leg1 * idAndIsoSF_leg2;
                }
            """)

    def run(self, df):
        if not self.isMC:
            return df, []

        df = df.Define("idAndIsoSF_leg1",
            "((pairType == 0) * musf_tight_id.at(dau1_index) * musf_tight_reliso.at(dau1_index))"
            "+ ((pairType == 1) * elesf_wp80iso.at(dau1_index))"
            "+ ((pairType == 2) * Tau_sfDeepTau2017v2p1VSjet_Medium.at(dau1_index)"
                " * Tau_sfDeepTau2017v2p1VSe_VVLoose.at(dau1_index)"
                " * Tau_sfDeepTau2017v2p1VSmu_VLoose.at(dau1_index))")
        df = df.Define("idAndIsoSF_leg2_deep_vsJet_pt",
            "Tau_sfDeepTau2017v2p1VSjet_Medium.at(dau2_index)")
        df = df.Define("idAndIsoSF_leg2_deep_vsEle",
            "Tau_sfDeepTau2017v2p1VSe_VVLoose.at(dau2_index)")
        df = df.Define("idAndIsoSF_leg2_deep_vsMu",
            "Tau_sfDeepTau2017v2p1VSmu_VLoose.at(dau2_index)")

        df = df.Define("idAndIsoAndFakeSF_deep_pt",
            "get_dauIdIso_sf(pairType, dau1_index, dau2_index, musf_tight_id, musf_tight_reliso, "
                "elesf_wp80iso, Tau_sfDeepTau2017v2p1VSjet_Medium, "
                "Tau_sfDeepTau2017v2p1VSe_VVLoose, Tau_sfDeepTau2017v2p1VSmu_VLoose)")

        return df, ["idAndIsoAndFakeSF_deep_pt"]


def dauIdIsoSFRDF(**kwargs):
    return lambda: dauIdIsoSFRDFProducer(**kwargs)
       