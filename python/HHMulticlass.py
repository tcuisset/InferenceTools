import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from analysis_tools.utils import import_root
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()

class HHMulticlassRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        year = kwargs.pop("year")
        model_specs = kwargs.pop("model_specs")
        special_chars = "{}()[] "
        for char in special_chars:
            model_specs = model_specs.replace(char, "")
        model_specs = model_specs.split(",")
        self.model_specs = [(model_specs[i], model_specs[i + 1])
            for i in range(0, len(model_specs), 2)]
        model_specs = ['{"%s", "%s"}' % (model_specs[i], model_specs[i + 1])
            for i in range(0, len(model_specs), 2)]

        super(HHMulticlassRDFProducer, self).__init__(*args, **kwargs)

        if not os.getenv("_HHbbttMulticlass"):
            os.environ["_HHbbttMulticlass"] = "_HHbbttMulticlass"

            if os.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc10":
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/"
                    "external/eigen/d812f411c3f9-cms/include/")
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc10/"
                    "external/tensorflow/2.5.0/include/")
            elif os.path.expandvars("$CMT_SCRAM_ARCH") == "slc7_amd64_gcc820":
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                    "external/eigen/d812f411c3f9-bcolbf/include/eigen3")
                ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc7_amd64_gcc820/"
                    "external/tensorflow/2.1.0-bcolbf/include")
            else:
                raise ValueError("Architecture not considered")

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gSystem.Load("libToolsTools.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/HHMulticlassInterface.h".format(base))

            ROOT.gInterpreter.Declare("""
                auto hhmdnn = HHMulticlassInterface(%s, {%s});
            """ % (year, ", ".join(model_specs)))
        
            ROOT.gInterpreter.Declare("""
                std::vector<std::string> get_node_names(size_t imodel) {
                    return hhmdnn.get_node_names(imodel);
                }
            """)

    def run(self, df):
        branches = ["mdnn_hhbbtt_%s_{}%s" % ("_".join(model_spec), self.systs)
            for model_spec in self.model_specs]
        all_branches = df.GetColumnNames()

        df = df.Define("mdnn_output%s" % self.systs, "hhmdnn.GetPredictionsWithInputs("
            "event, pairType, "
            "Jet_pt{0}, Jet_eta, Jet_phi, Jet_mass{0}, "
            "Jet_btagDeepFlavB, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB, Jet_HHbtag, "
            "dau1_index, dau2_index, "
            "bjet1_JetIdx, bjet2_JetIdx, VBFjet1_JetIdx, VBFjet2_JetIdx, "
            "ctjet_indexes, fwjet_indexes, "
            "Muon_pt{1}, Muon_eta, Muon_phi, Muon_mass{1}, "
            "Electron_pt{2}, Electron_eta, Electron_phi, Electron_mass{2}, "
            "Tau_pt{3}, Tau_eta, Tau_phi, Tau_mass{3}, "
            "MET{5}_pt{6}, MET{5}_phi{6}, "
            "Htt_svfit_pt{4}, Htt_svfit_eta{4}, Htt_svfit_phi{4}, Htt_svfit_mass{4})".format(
                self.jet_syst, self.muon_syst, self.electron_syst, self.tau_syst, self.systs,
                self.met_smear_tag, self.met_syst)
            )

        branches_to_save = []
        for ib, branch_tmp in enumerate(branches):
            node_names = [elem for elem in ROOT.get_node_names(ib)]
            for inode, node in enumerate(node_names):
                df = df.Define(branch_tmp.format(node),
                    "mdnn_output%s[%s][%s]" % (self.systs, ib, inode))
                branches_to_save.append(branch_tmp.format(node))

        return df, branches_to_save


def HHMulticlassRDF(**kwargs):
    """
    Returns the HH->bbtt Muticlass DNN outputs.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    :param model_specs: versions and models to be used. Syntaxis: ``((version, model), ...)``.
    :type model_specs: tuple of tuple of str

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHMulticlassRDF
            path: Tools.Tools.HHMulticlass
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                model_specs: ((v5 , kl1_c2v1_c31_vbf))

    """
    return lambda: HHMulticlassRDFProducer(**kwargs)
