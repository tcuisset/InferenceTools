from copy import deepcopy as copy

from analysis_tools.utils import import_root

from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()

class dauIdIsoSFRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(dauIdIsoSFRDFProducer, self).__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            default_systs = ["muon_id", "muon_iso", "ele_iso", "tau_vsjet", "tau_vse", "tau_vsmu"]
            systnames = kwargs.pop("systs", ["central"] + default_systs)
            template = ["" for i in range(len(default_systs))]
            self.syst_names = []
            self.systs = []
            for name in systnames:
                if name == "central":
                    self.syst_names.append("")
                    self.systs.append(template)
                    continue
                try:
                    ind = default_systs.index(name)
                    for d in ["_up", "_down"]:
                        tmp = copy(template)
                        self.syst_names.append("_%s%s" % (name, d))
                        tmp[ind] = d
                        self.systs.append(tmp)
                except ValueError:
                    raise ValueError("Systematic %s not available" % name)

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
        branches = []
        for name, syst_dir in zip(self.syst_names, self.systs):
            df = df.Define("idAndIsoAndFakeSF%s" % name,
                "get_dauIdIso_sf(pairType, dau1_index, dau2_index, Tau_pt{0}, "
                    "musf_tight_id{1[0]}, musf_tight_reliso{1[1]}, "
                    "elesf_wp80iso{1[2]}, Tau_sfDeepTau2017v2p1VSjet_pt_binned_Medium{1[3]}, "
                    "Tau_sfDeepTau2017v2p1VSjet_dm_binned_Medium{1[3]}, "
                    "Tau_sfDeepTau2017v2p1VSe_VVLoose{1[4]},"
                    "Tau_sfDeepTau2017v2p1VSmu_VLoose{1[5]})".format(self.tau_syst, syst_dir))
            branches.append("idAndIsoAndFakeSF%s" % name)

        return df, branches


def dauIdIsoSFRDF(**kwargs):
    """
    Returns the ID and Isolation SF applied to the two leptons with the desired systematics.

    Required RDFModules: :ref:`HHLepton_HHLeptonRDF`, :ref:`Electron_eleSFRDF`,
        :ref:`Muon_muSFRDF`, :ref:`Tau_tauSFRDF`

    :param systs: Systematics to be considered. Default: [``central``, ``muon_id``, ``muon_iso``,
        ``ele_iso``, ``tau_vsjet``, ``tau_vse``, ``tau_vsmu``]. 
    :type systs: list of str

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: dauIdIsoSFRDF
            path: Tools.Tools.dauIdIso
            parameters:
                isMC: self.dataset.process.isMC
                systs: [central, ...]

    """
    return lambda: dauIdIsoSFRDFProducer(**kwargs)
       