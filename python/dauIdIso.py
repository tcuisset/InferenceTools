from analysis_tools.utils import import_root

from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()

class dauIdIsoSFRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(dauIdIsoSFRDFProducer, self).__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                double get_dauIdIso_sf(
                        int pairType, int dau1_index, int dau2_index, Vfloat Tau_pt,
                        Vfloat musf_id, Vfloat musf_reliso, Vfloat elesf, 
                        Vfloat Tau_sfDeepTau2017v2p1VSjet_pt, Vfloat Tau_sfDeepTau2017v2p1VSjet_dm,
                        Vfloat Tau_sfDeepTau2017v2p1VSe, Vfloat Tau_sfDeepTau2017v2p1VSmu) {
                    double idAndIsoSF_leg1 = 1.;
                    if (pairType == 0) {
                        idAndIsoSF_leg1 = musf_id.at(dau1_index) * musf_reliso.at(dau1_index);
                    } else if (pairType == 1) {
                        idAndIsoSF_leg1 = elesf.at(dau1_index);
                    } else if (pairType == 2) {
                        if (Tau_pt[dau1_index] < 40)
                            idAndIsoSF_leg1 = Tau_sfDeepTau2017v2p1VSjet_pt.at(dau1_index);
                        else idAndIsoSF_leg1 = Tau_sfDeepTau2017v2p1VSjet_dm.at(dau1_index);

                        idAndIsoSF_leg1 *= Tau_sfDeepTau2017v2p1VSe.at(dau1_index) *
                            Tau_sfDeepTau2017v2p1VSmu.at(dau1_index);
                    }
                    double idAndIsoSF_leg2 = 1.;
                    if (pairType == 2 && Tau_pt[dau1_index] > 40)
                        idAndIsoSF_leg2 = Tau_sfDeepTau2017v2p1VSjet_dm.at(dau2_index);
                    else
                        idAndIsoSF_leg2 = Tau_sfDeepTau2017v2p1VSjet_pt.at(dau2_index);
                    idAndIsoSF_leg2 *= Tau_sfDeepTau2017v2p1VSe.at(dau2_index) *
                        Tau_sfDeepTau2017v2p1VSmu.at(dau2_index);
                    return idAndIsoSF_leg1 * idAndIsoSF_leg2;
                }
            """)

    def run(self, df):
        if not self.isMC:
            return df, []

        df = df.Define("idAndIsoAndFakeSF_deep_pt",
            "get_dauIdIso_sf(pairType, dau1_index, dau2_index, Tau_pt%s, "
                "musf_tight_id, musf_tight_reliso, "
                "elesf_wp80iso, Tau_sfDeepTau2017v2p1VSjet_pt_binned_Medium, "
                "Tau_sfDeepTau2017v2p1VSjet_dm_binned_Medium, "
                "Tau_sfDeepTau2017v2p1VSe_VVLoose, Tau_sfDeepTau2017v2p1VSmu_VLoose)"
                % self.tau_syst)

        return df, ["idAndIsoAndFakeSF_deep_pt"]


def dauIdIsoSFRDF(**kwargs):
    return lambda: dauIdIsoSFRDFProducer(**kwargs)
       