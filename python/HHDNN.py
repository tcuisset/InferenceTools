import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from analysis_tools.utils import import_root, randomize
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()


class HHDNNProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(HHDNNProducer, self).__init__(*args, **kwargs)
        base = "{}/{}/src/Tools/Tools".format(
            os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

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

        ROOT.gSystem.Load("libToolsTools.so")
        ROOT.gROOT.ProcessLine(".L {}/interface/HHDNNinterface.h".format(base))

        feature_file = kwargs.pop(
            "feature_file",
            "{}/{}/src/cms_runII_dnn_models/models/nonres_gluglu/2020-07-31-0/features.txt".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        )
        req_features = ROOT.vector(str)()
        with open(feature_file) as f:
            for feature in f.readlines():
                req_features.push_back(feature)       

        model_dir = kwargs.pop(
            "model_dir",
            "{}/{}/src/cms_runII_dnn_models/models/nonres_gluglu/2020-07-31-0/ensemble".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        )

        kls = ROOT.vector(float)()
        kls.push_back(1)  # FIXME, may not be needed

        self.hhdnn = ROOT.HHDNNinterface(model_dir, req_features, kls)
        self.hhdnn.SetGlobalInputs(kwargs.pop("year", 2018), 2)
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # self.out.branch("DNNoutSM_kl_1%s" % self.systs, "F")

        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        return True

        """process event, return True (go to next module) or False (fail, go to next event)"""
        htt_svfit = Object(event, "Htt_svfit")
        HHKinFit = Object(event, "HHKinFit")
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")
        jets = Collection(event, "Jet")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus) #FIXME FROM HERE
        bjet1, bjet2, bjet1_tlv, bjet2_tlv = self.get_bjets(event, jets)
        vbfjet1, vbfjet2, vbfjet1_tlv, vbfjet2_tlv = self.get_vbfjets(event, jets)
        if not vbfjet1:  # dummy tlv
            vbfjet1_tlv.SetPtEtaPhiM(1, 1, 1, 1)
            vbfjet2_tlv.SetPtEtaPhiM(1, 1, 1, 1)
        met, met_tlv = self.get_met(event)
        htt_svfit_tlv = ROOT.TLorentzVector()
        htt_svfit_tlv.SetPtEtaPhiM(
            eval("htt_svfit.pt%s" % self.systs), 
            eval("htt_svfit.eta%s" % self.systs), 
            eval("htt_svfit.phi%s" % self.systs), 
            eval("htt_svfit.mass%s" % self.systs),
        )

        # the definition of channels is a bit different between the DNN features
        # and our definition

        # Channel | Our Def. | DNN
        # ------------------------
        # mutau   |    0     |  1
        # etau    |    1     |  2
        # tautau  |    2     |  0
        
        if event.pairType == 0:
            channel = 1
        elif event.pairType == 1:
            channel = 2
        elif event.pairType == 2:
            channel = 0
        else:
            raise ValueError("pairType %s is not implemented" % event.pairType)

        self.hhdnn.SetEventInputs(channel, event.isBoosted, (2 if vbfjet1 else 0), event.event,
            bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, vbfjet1_tlv, vbfjet2_tlv, met_tlv,
            htt_svfit_tlv, HHKinFit.mass, HHKinFit.chi2, (HHKinFit.chi2 > 0),
            (eval("htt_svfit.mass%s" % self.systs) > 0),)

        # void SetEventInputs(Channel channel, int is_boosted, int nvbf, unsigned long long int eventn,
          # TLorentzVector b1, TLorentzVector b2, TLorentzVector l1, TLorentzVector l2, 
          # TLorentzVector vbf1, TLorentzVector vbf2, TLorentzVector met, TLorentzVector svfit, 
          # float KinFitMass, float KinFitChi2, bool KinFitConv, bool SVfitConv,
          # float HHbtag_b1, float HHbtag_b2, float HHbtag_vbf1, float HHbtag_vbf2,
          # float CvsL_b1, float CvsL_b2, float CvsL_vbf1, float CvsL_vbf2,
          # float CvsB_b1, float CvsB_b2, float CvsB_vbf1, float CvsB_vbf2,
          # float cv, float c2v, float c3, bool pass_massCut
        # );


        return True


def HH(**kwargs):
    return lambda: HHDNNProducer(**kwargs)


class HHDNNInputRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.AnalysisType = AnalysisType
        self.met_smear_tag_data = kwargs.pop("met_smear_tag_data")
        self.add_systs_to_output_branches = kwargs.pop("add_systs_to_output_branches", False) # set to True to not add self.systs as branch suffix (for DNN training with systematics)
        super(HHDNNInputRDFProducer, self).__init__(*args, **kwargs)

        if not os.getenv("_HHbbttDNNInterface"):
            os.environ["_HHbbttDNNInterface"] = "_HHbbttDNNInterface"
            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gSystem.Load("libToolsTools.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/HHDNNinterfaceNew.h".format(base))
            # ROOT.gROOT.ProcessLine(".L {}/interface/lester_mt2_bisect.h".format(base))

        if not self.AnalysisType:
            feature_file = "src/cms_runII_dnn_models/models/nonres_gluglu/2020-07-31-0/features.txt"
        else:
            feature_file = "src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-02-15/ZZbbtt-0/features.txt"
        feature_file = "{}/{}/{}".format(os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"), kwargs.pop("feature_file", feature_file)) 
        with open(feature_file) as f:
            self.features = [i.split('\n')[0] for i in f.readlines()]
        try:
            self.features.remove("res_mass")
        except ValueError: pass

        makeStr = lambda l: (", ".join([str(x) for x in l]))
        btag_wps = makeStr([0., kwargs["btag_wps"]["loose"], kwargs["btag_wps"]["medium"], kwargs["btag_wps"]["tight"]])
        cvl_wps = makeStr([0., kwargs["cvl_wps"]["loose"], kwargs["cvl_wps"]["medium"], kwargs["cvl_wps"]["tight"]])
        cvb_wps = makeStr([0.4,0.29,0.05,0.]) # hardcoded to the values used for Run2 legacy. The new WPs for cVSb don't have increasing thresholds so cannot be used with the current setup
        # UL 2018 DeepFlavour WPs for example : loose 0.246, medium 0.325, tight 0.267 (they are used in conjunction with cvsl), see https://indico.cern.ch/event/967689/contributions/4083041/attachments/2130779/3592013/CTagger_WP_UltraLegacy2018_DeepCSV_DeepJet.pdf

        self.variable_dnnInput = "hhdnninterface_" + "__".join(self.features)

        if not os.getenv("_HHbbttDNNInterface_" + self.variable_dnnInput):
            os.environ["_HHbbttDNNInterface_" + self.variable_dnnInput] = "_HHbbttDNNInterface_" + self.variable_dnnInput
            ROOT.gInterpreter.Declare("""
                auto %s = HHDNNinterfaceNew(%s, {%s}, {%s}, {%s});
            """ % (self.variable_dnnInput, self.year, btag_wps, cvl_wps, cvb_wps))
            # // ROOT::RVec<Float_t> const hhdnn_placeholder(%s, -1.); // vector to be used when we cannot compute the DNN inputs (boosted case)
            

    def run(self, df):

        # DEBUG
        if not self.AnalysisType:
            p_b = "H"; p_t = "H"; pp = "HH"
            print(" ### INFO: Running HHDNN with default option for HH analysis")
        else:
            print(" ### INFO: Running HHDNN with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      p_b = "Z"; p_t = "Z"; pp = "ZZ"; p_sv = "X"
            elif self.AnalysisType == "Zbb_Htautau":    p_b = "Z"; p_t = "H"; pp = "ZH"; p_sv = "X"
            elif self.AnalysisType == "Ztautau_Hbb":    p_b = "H"; p_t = "Z"; pp = "ZH"; p_sv = "X"

        met_smear_tag = self.met_smear_tag if self.isMC else self.met_smear_tag_data
        df = df.Define("filterForDNNInputs", "(jetCategory !=-1) && pairType>=0")
        # "fake" filter that relies on side effects. Will not filter anything but will run get_dnn_inputs
        # not very elegant but works
        df = df.Filter(f"""
            if (filterForDNNInputs) 
            get_dnn_inputs(
            {self.variable_dnnInput}, pairType, jetCategory==2||jetCategory==-2, isBoostedTau, event, 
            bjet1_JetIdx, bjet2_JetIdx,
            fatjet_JetIdx, VBFjet1_JetIdx, VBFjet2_JetIdx, 
            dau1_pt{self.lep_syst}, dau1_eta, dau1_phi, dau1_mass{self.lep_syst}, 
            dau2_pt{self.lep_syst}, dau2_eta, dau2_phi, dau2_mass{self.lep_syst}, 
            bjet1_pt{self.jet_syst}, bjet1_eta, bjet1_phi, bjet1_mass{self.jet_syst},
            bjet2_pt{self.jet_syst}, bjet2_eta, bjet2_phi, bjet2_mass{self.jet_syst},
            fatjet_pt{self.jet_syst}, fatjet_eta, fatjet_phi, fatjet_mass{self.jet_syst}, fatjet_msoftdrop,
            {p_sv}tt_svfit_pt{self.systs}, {p_sv}tt_svfit_eta{self.systs}, {p_sv}tt_svfit_phi{self.systs}, {p_sv}tt_svfit_mass{self.systs}, 
            {pp}KinFit_mass{self.systs}, {pp}KinFit_chi2{self.systs}, MET{met_smear_tag}_pt{self.met_syst}, MET{met_smear_tag}_phi{self.met_syst}, 
            bjet1_btagDeepFlavB, bjet1_btagDeepFlavCvL, bjet1_btagDeepFlavCvB,
            bjet2_btagDeepFlavB, bjet2_btagDeepFlavCvL, bjet2_btagDeepFlavCvB,
            bjet1_HHbtag, bjet2_HHbtag
            );
            return true;""")

        branches = []
        for feature in self.features:
            if self.add_systs_to_output_branches:
                branch = f"{feature}{self.systs}"
            else:
                branch = feature
            df = df.Define(branch, f"filterForDNNInputs ? {self.variable_dnnInput}.{feature}() : -1")
            branches.append(branch)
        return df, branches

class HHDNNRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, DNN_res_mass:float, *args, **kwargs):
        self.AnalysisType = AnalysisType
        self.DNN_res_mass = DNN_res_mass
        self.resonant_dnn = DNN_res_mass >= 0.
        # For the resonant DNN, we add the resonant mass to the branch name as we want to run different mass points in the same processing
        # We also add the suffix to the global C++ variable names
        if self.resonant_dnn:
            self.resonant_suffix = f"_{self.DNN_res_mass:.0f}"
        else:
            self.resonant_suffix = ""
        super(HHDNNRDFProducer, self).__init__(*args, **kwargs)

        if not os.getenv("_HHbbttDNNDefault"):
            os.environ["_HHbbttDNNDefault"] = "HHbbttDNNDefault"

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

        if not os.getenv("_HHbbttDNN"+self.resonant_suffix):
            os.environ["_HHbbttDNN"+self.resonant_suffix] = "HHbbttDNN"+self.resonant_suffix

            if not self.AnalysisType:
                model_dir = "{}/{}/src/cms_runII_dnn_models/models/nonres_gluglu/2020-07-31-0/".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            elif self.AnalysisType == "Zbb_Ztautau":
                if self.resonant_dnn:
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-04-29/ResZZbbtt-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                else:
                    # model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2023-08-02-0/".format( # old model
                    # model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-02-15/ZZbbtt-0/".format( # old model 2018
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-03-26/ZZbbtt-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            elif self.AnalysisType == "Zbb_Htautau": # or self.AnalysisType == "Ztautau_Hbb":
                if self.resonant_dnn:
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-04-29/ResZbbHtt-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                else:
                    # model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-02-15/ZbbHtt-0/".format( # old model 2018
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-03-26/ZbbHtt-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            elif self.AnalysisType == "Ztautau_Hbb":
                if self.resonant_dnn:
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-04-29/ResZttHbb-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                else:
                    # model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-02-15/ZttHbb-0/".format( # old model 2018
                    model_dir = "{}/{}/src/cms_runII_dnn_models/models/arc_checks/zz_bbtt/2024-03-26/ZttHbb-0/".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

            ensemble = model_dir + "ensemble"
            features = model_dir + "features.txt"

            feature_file = kwargs.pop("feature_file", features)
            with open(feature_file) as f:
                lines = f.readlines()
            req_features = ', '.join(['"%s"' % line.strip() for line in lines])

            model_dir = kwargs.pop("model_dir", ensemble)

            ROOT.gInterpreter.Declare("""
                auto hhdnn%s = HHDNNinterface("%s", {%s}, {1.}, %s);
            """ % (self.resonant_suffix, model_dir, req_features, self.year))

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<Float_t>&;
                ROOT::RVec<Float_t> get_dnn_outputs%s(int pairType, int isBoosted, int event,
                    int dau1_index, int dau2_index, int bjet1_index, int bjet2_index,
                    int vbfjet1_index, int vbfjet2_index,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass,
                    Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                    float htt_sv_pt, float htt_sv_eta, float htt_sv_phi, float htt_sv_mass,
                    float HHKinFit_mass, float HHKinFit_chi2, float met_pt, float met_phi,
                    Vfloat Jet_btagDeepFlavB, Vfloat Jet_btagDeepFlavCvL, Vfloat Jet_btagDeepFlavCvB,
                    Vfloat Jet_HHbtag, float DNN_res_mass
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
                    else
                        assert(false);
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

                    auto htt_svfit_tlv = TLorentzVector();
                    if (htt_sv_mass > 0)
                        htt_svfit_tlv.SetPtEtaPhiM(htt_sv_pt, htt_sv_eta, htt_sv_phi, htt_sv_mass);
                    else
                        htt_svfit_tlv.SetPtEtaPhiM(1., 1., 1., 1.);

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

                    return hhdnn%s.GetPredictionsWithInputs(
                      channel, isBoosted, nvbf, event,
                      bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, 
                      vbfjet1_tlv, vbfjet2_tlv, met_tlv, htt_svfit_tlv, 
                      HHKinFit_mass, HHKinFit_chi2, (HHKinFit_chi2 >= 0), (htt_sv_mass >= 0), MT2,
                      deepFlav1, deepFlav2, CvsL_b1, CvsL_b2, CvsL_vbf1, CvsL_vbf2,
                      CvsB_b1, CvsB_b2, CvsB_vbf1, CvsB_vbf2,
                      HHbtag_b1, HHbtag_b2, HHbtag_vbf1, HHbtag_vbf2, DNN_res_mass);
                }
            """ % (self.resonant_suffix, self.resonant_suffix))

    def run(self, df):

        # DEBUG
        if not self.AnalysisType:
            p_b = "H"; p_t = "H"; pp = "HH"
            print(" ### INFO: Running HHDNN with default option for HH analysis")
        else:
            print(" ### INFO: Running HHDNN with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      p_b = "Z"; p_t = "Z"; pp = "ZZ"; p_sv = "X"
            elif self.AnalysisType == "Zbb_Htautau":    p_b = "Z"; p_t = "H"; pp = "ZH"; p_sv = "X"
            elif self.AnalysisType == "Ztautau_Hbb":    p_b = "H"; p_t = "Z"; pp = "ZH"; p_sv = "X"
        
        branches = ["dnn_%sbbtt_kl_1%s%s" % (pp, self.resonant_suffix, self.systs)]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        if self.jet_category == "boosted":
            isBoostedJet = "true"
        elif self.jet_category == "resolved":
            isBoostedJet = "false"
        else:
            raise ValueError()
        df = df.Define("dnn_output%s%s" % (self.resonant_suffix, self.systs), f"get_dnn_outputs{self.resonant_suffix}("
            f"pairType, {isBoostedJet}, event, "
            f"dau1_index, dau2_index, bjet1_JetIdx, bjet2_JetIdx,VBFjet1_JetIdx, VBFjet2_JetIdx, "
            f"Muon_pt{self.muon_syst}, Muon_eta, Muon_phi, Muon_mass{self.muon_syst},"
            f"Electron_pt{self.electron_syst}, Electron_eta, Electron_phi, Electron_mass{self.electron_syst}, "
            f"Tau_pt{self.tau_syst}, Tau_eta, Tau_phi, Tau_mass{self.tau_syst}, "
            f"Jet_pt{self.jet_syst}, Jet_eta, Jet_phi, Jet_mass{self.jet_syst}, "
            f"{p_sv}tt_svfit_pt{self.systs}, {p_sv}tt_svfit_eta{self.systs}, {p_sv}tt_svfit_phi{self.systs}, {p_sv}tt_svfit_mass{self.systs}, "
            f"{pp}KinFit_mass{self.systs}, {pp}KinFit_chi2{self.systs}, MET{self.met_smear_tag}_pt{self.met_syst}, MET{self.met_smear_tag}_phi{self.met_syst}, "
            f"Jet_btagDeepFlavB, Jet_btagDeepFlavCvL, Jet_btagDeepFlavCvB, Jet_HHbtag, {self.DNN_res_mass})"
            ).Define(branches[0], "dnn_output%s%s[0]" % (self.resonant_suffix, self.systs))
        return df, branches

def HHDNNInputRDF(**kwargs):
    """
    Returns the DNN input.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHDNNInputRDF
            path: Tools.Tools.HHDNN
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                AnalysisType: self.config.get_aux('AnalysisType', False)
                btag_wps: self.config.btag_algo_wps
                cvl_wps: self.config.ctag_algo_wps.CvL
                cvb_wps: self.config.ctag_algo_wps.CvB

    """
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 1 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHDNNInputRDFProducer(AnalysisType=AnalysisType, **kwargs)

def HHDNNRDF(**kwargs):
    """
    Returns the HH->bbtt or ZZ/ZH->bbtt DNN output.

    Lepton and jet systematics (used for pt and mass variables) can be modified using the parameters
    from :ref:`BaseModules_JetLepMetSyst`.

    Parameter DNN_res_mass : set to -1 for non-resonant DNN. set to a positive mass value, to use parametrized resonant DNN.
    The mass is encoded into the output branch name.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHDNNRDF
            path: Tools.Tools.HHDNN
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                AnalysisType: self.config.get_aux('AnalysisType', False)
                DNN_res_mass: -1

    """
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 2 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHDNNRDFProducer(AnalysisType=AnalysisType, DNN_res_mass=kwargs.pop("DNN_res_mass", -1.), **kwargs)

