import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from analysis_tools.utils import import_root
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()

class ZZDNNRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        year = kwargs.pop("year")
        super(ZZDNNRDFProducer, self).__init__(*args, **kwargs)

        if not os.getenv("_HHbbttDNNNewout"):
            os.environ["_HHbbttDNNNewout"] = "HHbbttDNNNew"

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
            ROOT.gROOT.ProcessLine(".L {}/interface/HHDNNinterface.h".format(base))
            ROOT.gROOT.ProcessLine(".L {}/interface/lester_mt2_bisect.h".format(base))

            feature_file = kwargs.pop(
                "feature_file",
                "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2023-08-02-0/features.txt".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            )
            with open(feature_file) as f:
                lines = f.readlines()
            req_features = ', '.join(['"%s"' % line.strip() for line in lines])

            model_dir = kwargs.pop(
                "model_dir",
                "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2023-08-02-0/ensemble".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            )

            ROOT.gInterpreter.Declare("""
                auto hhdnn = HHDNNinterface("%s", {%s}, {1.}, %s);
            """ % (model_dir, req_features, year))

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<Float_t>&;
                ROOT::RVec<Float_t> get_dnn_outputs(int pairType, int isBoosted, int event,
                    int dau1_index, int dau2_index, int bjet1_index, int bjet2_index,
                    int vbfjet1_index, int vbfjet2_index,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                    Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                    float ztt_sv_pt, float ztt_sv_eta, float ztt_sv_phi, float ztt_sv_mass,
                    float ZZKinFit_mass, float ZZKinFit_chi2, float met_pt, float met_phi,
                    Vfloat Jet_btagDeepFlavB, Vfloat Jet_btagDeepFlavCvL, Vfloat Jet_btagDeepFlavCvB,
                    Vfloat Jet_HHbtag
                )
                {
                    float dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;

                    // pairType | Our Def. | DNN
                    // ------------------------
                    // mutau   |    0     |  1
                    // etau    |    1     |  2
                    // tautau  |    2     |  0

                    int channel;
                    if(pairType == 0)
                        channel = 1;
                    else if (pairType == 1)
                        channel = 2;
                    else
                        channel = 0;

                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_eta = muon_eta.at(dau1_index);
                        dau1_phi = muon_phi.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_eta = electron_eta.at(dau1_index);
                        dau1_phi = electron_phi.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_eta = tau_eta.at(dau1_index);
                        dau1_phi = tau_phi.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_eta = tau_eta.at(dau2_index);
                    dau2_phi = tau_phi.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);

                    auto dau1_tlv = TLorentzVector();
                    auto dau2_tlv = TLorentzVector();
                    auto bjet1_tlv = TLorentzVector();
                    auto bjet2_tlv = TLorentzVector();
                    auto vbfjet1_tlv = TLorentzVector();
                    auto vbfjet2_tlv = TLorentzVector();

                    dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
                    dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                    bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                        jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                    bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                        jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));

                    int nvbf = 0;
                    if (vbfjet1_index >= 0) {
                        vbfjet1_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet1_index), jet_eta.at(vbfjet1_index),
                            jet_phi.at(vbfjet1_index), jet_mass.at(vbfjet1_index));
                        vbfjet2_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet2_index), jet_eta.at(vbfjet2_index),
                            jet_phi.at(vbfjet2_index), jet_mass.at(vbfjet2_index));
                        nvbf = 2;
                    } else {
                        vbfjet1_tlv.SetPtEtaPhiM(1., 1., 1., 1.);
                        vbfjet2_tlv.SetPtEtaPhiM(1., 1., 1., 1.);
                    }
                    auto met_tlv = TLorentzVector();
                    met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

                    auto ztt_svfit_tlv = TLorentzVector();
                    if (ztt_sv_mass > 0)
                        ztt_svfit_tlv.SetPtEtaPhiM(ztt_sv_pt, ztt_sv_eta, ztt_sv_phi, ztt_sv_mass);
                    else
                        ztt_svfit_tlv.SetPtEtaPhiM(1., 1., 1., 1.);

                    // MT2 computation
                    asymm_mt2_lester_bisect::disableCopyrightMessage();
                    double MT2 = asymm_mt2_lester_bisect::get_mT2(
                        bjet1_tlv.M(), bjet1_tlv.Px(), bjet1_tlv.Py(),
                        bjet2_tlv.M(), bjet2_tlv.Px(), bjet2_tlv.Py(),
                        dau1_tlv.Px() + dau2_tlv.Px() + met_tlv.Px(),
                        dau1_tlv.Py() + dau2_tlv.Py() + met_tlv.Py(),
                        dau1_tlv.M(), dau2_tlv.M(), 0.);

                    float deepFlav1 = -1., deepFlav2 = -1., CvsL_b1 = -1., CvsL_b2 = -1.,
                        CvsL_vbf1 = -1., CvsL_vbf2 = -1., CvsB_b1 = -1., CvsB_b2 = -1.,
                        CvsB_vbf1 = -1., CvsB_vbf2 = -1., HHbtag_b1 = -1., HHbtag_b2 = -1.,
                        HHbtag_vbf1 = -1., HHbtag_vbf2 = -1.;
                    deepFlav1 = Jet_btagDeepFlavB.at(bjet1_index);
                    deepFlav2 = Jet_btagDeepFlavB.at(bjet2_index);
                    CvsL_b1 = Jet_btagDeepFlavCvL.at(bjet1_index);
                    CvsL_b2 = Jet_btagDeepFlavCvL.at(bjet2_index);
                    CvsB_b1 = Jet_btagDeepFlavCvB.at(bjet1_index);
                    CvsB_b2 = Jet_btagDeepFlavCvB.at(bjet2_index);
                    HHbtag_b1 = Jet_HHbtag.at(bjet1_index);
                    HHbtag_b2 = Jet_HHbtag.at(bjet2_index);
                    if (vbfjet1_index >= 0) {
                        CvsL_vbf1 = Jet_btagDeepFlavCvL.at(vbfjet1_index);
                        CvsL_vbf2 = Jet_btagDeepFlavCvL.at(vbfjet2_index);
                        CvsB_vbf1 = Jet_btagDeepFlavCvB.at(vbfjet1_index);
                        CvsB_vbf2 = Jet_btagDeepFlavCvB.at(vbfjet2_index);
                        HHbtag_vbf1 = Jet_HHbtag.at(vbfjet1_index);
                        HHbtag_vbf2 = Jet_HHbtag.at(vbfjet2_index);
                    }

                    return hhdnn.GetPredictionsWithInputs(
                      channel, isBoosted, nvbf, event,
                      bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, 
                      vbfjet1_tlv, vbfjet2_tlv, met_tlv, ztt_svfit_tlv, 
                      ZZKinFit_mass, ZZKinFit_chi2, (ZZKinFit_chi2 >= 0), (ztt_sv_mass >= 0), MT2,
                      deepFlav1, deepFlav2, CvsL_b1, CvsL_b2, CvsL_vbf1, CvsL_vbf2,
                      CvsB_b1, CvsB_b2, CvsB_vbf1, CvsB_vbf2,
                      HHbtag_b1, HHbtag_b2, HHbtag_vbf1, HHbtag_vbf2);
                }
            """)

    def run(self, df):
        branches = ["dnn_new_zzbbtt_kl_1%s" % self.systs]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        print(" ### DEBUG 1: Running New DNN")
        df = df.Define("dnn_output%s" % self.systs, "get_dnn_outputs("
            "pairType, isBoosted, event, "
            "dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx,VBFjet1_JetIdx, VBFjet2_JetIdx, "
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0},"
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, "
            "Ztt_svfit_pt{4}, Ztt_svfit_eta{4}, Ztt_svfit_phi{4}, Ztt_svfit_mass{4}, "
            "ZZKinFit_mass{4}, ZZKinFit_chi2{4}, MET{5}_pt{6}, MET{5}_phi{6}, "
            "Jet_btagDeepFlavB, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB, Jet_HHbtag)".format(
                self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst, self.systs,
                self.met_smear_tag, self.met_syst)
            ).Define(branches[0], "dnn_output%s[0]" % self.systs)
        return df, branches

def ZZDNNRDF(**kwargs):
    """
    Returns the HH->bbtt DNN output.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: ZZDNNRDF
            path: Tools.Tools.ZZDNN
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year

    """
    return lambda: ZZDNNRDFProducer(**kwargs)

