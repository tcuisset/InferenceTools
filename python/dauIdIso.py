from copy import deepcopy as copy
import os

from analysis_tools.utils import import_root

from Base.Modules.baseModules import JetLepMetSyst
from Corrections.TAU.tauCorrections import listAllTauVSJetSystematics_pt, listAllTauVSJetSystematics_dm

ROOT = import_root()

class dauIdIsoSFRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(dauIdIsoSFRDFProducer, self).__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")
        self.runPeriod = kwargs.pop("runPeriod")
        self.computeSystematics = kwargs.pop("computeSystematics", True)

        if self.isMC:
            
            default_systs = ["muon_id", "muon_iso", "ele_reco", "ele_iso", "tau_vsjet_pt", "tau_vsjet_dm", "tau_vse", "tau_vsmu"]
            systnames = ["central"] 
            if self.computeSystematics and self.systematic_is_central:
                systnames += default_systs
            base_template = ["" for i in range(len(default_systs))]
            self.branch_names = []
            """ Suffixes of branch names to go with branch_templates """
            self.branch_templates = []
            """ List of the set of suffixes to apply to input branches. Each outer list is a 6-element list [mu SF ID suffix, mu SF ISO suffix, ....]
            for example ["", "syst_up", "", ...] for the mu SF ISO syst_up variation """
            for name in systnames:
                if name == "central":
                    self.branch_names.append("")
                    self.branch_templates.append(base_template)
                    continue
                try:
                    ind = default_systs.index(name)
                    dirs = kwargs.pop(name + "_syst_directions", ["_up", "_down"])
                    if name == "tau_vsjet_pt" and dirs == "all":
                        dirs = ["_"+x for x in listAllTauVSJetSystematics_pt(self.year, self.runPeriod, addCombined=False) if x != "nom"]
                    elif name == "tau_vsjet_dm" and dirs == "all":
                        dirs = ["_"+x for x in listAllTauVSJetSystematics_dm(self.year, self.runPeriod, addCombined=False) if x != "nom"]
                    for d in dirs:
                        tmp = copy(base_template)
                        self.branch_names.append("_%s%s" % (name, d))
                        tmp[ind] = d
                        self.branch_templates.append(tmp)
                except ValueError:
                    raise ValueError("Systematic %s not available" % name)
            
            if not os.getenv("_dauIdIsoSF"):
                os.environ["_dauIdIsoSF"] = "_dauIdIsoSF"
                ROOT.gInterpreter.Declare("""
                    using Vfloat = const ROOT::RVec<float>&;
                    double get_dauIdIso_sf(
                            int pairType, bool isBoostedTau, int dau1_index, int dau2_index, Vfloat Tau_pt,
                            Vfloat musf_id, Vfloat musf_reliso, Vfloat elesf_reco, Vfloat elesf_idiso, 
                            Vfloat Tau_sfDeepTau2017v2p1VSjet_pt, Vfloat Tau_sfDeepTau2017v2p1VSjet_dm,
                            Vfloat Tau_sfDeepTau2017v2p1VSe, Vfloat Tau_sfDeepTau2017v2p1VSmu_tautau, Vfloat Tau_sfDeepTau2017v2p1VSmu_lepton) {
                        if (pairType < 0) return 1.;
                        // We use separate VSmu WPs for tautau channel and for etau/mutau channel
                        double idAndIsoSF_leg1 = 1.;
                        if (pairType == 0) {
                            idAndIsoSF_leg1 = musf_id.at(dau1_index) * musf_reliso.at(dau1_index);
                        } else if (pairType == 1) {
                            idAndIsoSF_leg1 = elesf_reco.at(dau1_index) * elesf_idiso.at(dau1_index);
                        } else if (pairType == 2 && !isBoostedTau) {
                            if (Tau_pt[dau1_index] > 140)
                                idAndIsoSF_leg1 = Tau_sfDeepTau2017v2p1VSjet_pt.at(dau1_index);
                            else idAndIsoSF_leg1 = Tau_sfDeepTau2017v2p1VSjet_dm.at(dau1_index);

                            idAndIsoSF_leg1 *= Tau_sfDeepTau2017v2p1VSe.at(dau1_index) *
                                Tau_sfDeepTau2017v2p1VSmu_tautau.at(dau1_index);
                        }
                        double idAndIsoSF_leg2 = 1.;
                        if (!isBoostedTau) {
                            if (Tau_pt[dau2_index] > 140)
                                idAndIsoSF_leg2 = Tau_sfDeepTau2017v2p1VSjet_pt.at(dau2_index);
                            else
                                idAndIsoSF_leg2 = Tau_sfDeepTau2017v2p1VSjet_dm.at(dau2_index);
                            
                            if (pairType == 2)
                                idAndIsoSF_leg2 *= Tau_sfDeepTau2017v2p1VSmu_tautau.at(dau2_index);
                            else
                                idAndIsoSF_leg2 *= Tau_sfDeepTau2017v2p1VSmu_lepton.at(dau2_index);
        
                            idAndIsoSF_leg2 *= Tau_sfDeepTau2017v2p1VSe.at(dau2_index);
                        }
                        return idAndIsoSF_leg1 * idAndIsoSF_leg2;
                    }
                """)

    def run(self, df):
        if not self.isMC:
            return df, []
        branches = []
        for branch_name, branch_template in zip(self.branch_names, self.branch_templates):
            assert(len(branch_template) == 8)
            df = df.Define("idAndIsoAndFakeSF%s" % branch_name,
                "get_dauIdIso_sf(pairType, isBoostedTau, dau1_index, dau2_index, Tau_pt{0}, "
                    "musf_tight_id{1[0]}, musf_tight_reliso{1[1]}, "
                    "elesf_RecoAbove20{1[2]}, isBoostedTau ? elesf_Loose{1[3]} : elesf_wp80iso{1[3]}, Tau_sfDeepTau2017v2p1VSjet_pt_binned_Medium{1[4]}, "
                    "Tau_sfDeepTau2017v2p1VSjet_dm_binned_Medium{1[5]}, "
                    "Tau_sfDeepTau2017v2p1VSe_VVLoose{1[6]},"
                    "Tau_sfDeepTau2017v2p1VSmu_VLoose{1[7]}, Tau_sfDeepTau2017v2p1VSmu_Tight{1[7]})".format(self.tau_syst, branch_template))
            branches.append("idAndIsoAndFakeSF%s" % branch_name)

        return df, branches


def dauIdIsoSFRDF(**kwargs):
    """
    Returns the ID and Isolation SF applied to the two leptons with the desired systematics.

    Required RDFModules: :ref:`HHLepton_HHLeptonRDF`, :ref:`Electron_eleSFRDF`,
    :ref:`Muon_muSFRDF`, :ref:`Tau_tauSFRDF`
    Input branches : 
     - musf_tight_id, musf_tight_reliso
     - elesf_RecoAbove20, elesf_wp80iso, elesf_Loose
     - Tau_sfDeepTau2017v2p1VSjet_pt_binned_Medium & Tau_sfDeepTau2017v2p1VSjet_dm_binned_Medium
     - Tau_sfDeepTau2017v2p1VSmu_Tight & Tau_sfDeepTau2017v2p1VSmu_VLoose : the DeepTau SFS for the etau/mutau and the tautau channels respectively (different VsMu working points -> different SFs)

    Systematics considered : [``central``, ``muon_id``, ``muon_iso``,
        ```ele_reco```, ``ele_iso``, ``tau_vsjet``, ``tau_vse``, ``tau_vsmu``]. 
    
    For choosing, for one of the systematics considered, which variations to use (on top of central) :
    for example : muon_iso_syst_directions : default [_syst_up, _syst_down, _stat_up, _stat_down]
    default _syst_directions : [_up, _down]

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: dauIdIsoSFRDF
            path: Tools.Tools.dauIdIso
            parameters:
                isMC: self.dataset.process.isMC
                muon_iso_syst_directions: [_syst_up, _syst_down, _stat_up, _stat_down, ...]

    """
    return lambda: dauIdIsoSFRDFProducer(**kwargs)
       