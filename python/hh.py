import os
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst
from analysis_tools.utils import import_root

ROOT = import_root()


class HHProducer(JetLepMetModule):
    def __init__(self, *args, **kwargs):
        super(HHProducer, self).__init__(*args, **kwargs)
        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        ROOT.gROOT.ProcessLine(".L {}/interface/HHKinFitInterface.h".format(base))
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("Hbb_pt%s" % self.systs, "F")
        self.out.branch("Hbb_eta%s" % self.systs, "F")
        self.out.branch("Hbb_phi%s" % self.systs, "F")
        self.out.branch("Hbb_mass%s" % self.systs, "F")

        self.out.branch("Htt_pt%s" % self.systs, "F")
        self.out.branch("Htt_eta%s" % self.systs, "F")
        self.out.branch("Htt_phi%s" % self.systs, "F")
        self.out.branch("Htt_mass%s" % self.systs, "F")

        self.out.branch("HH_pt%s" % self.systs, "F")
        self.out.branch("HH_eta%s" % self.systs, "F")
        self.out.branch("HH_phi%s" % self.systs, "F")
        self.out.branch("HH_mass%s" % self.systs, "F")

        self.out.branch("HH_svfit_pt%s" % self.systs, "F")
        self.out.branch("HH_svfit_eta%s" % self.systs, "F")
        self.out.branch("HH_svfit_phi%s" % self.systs, "F")
        self.out.branch("HH_svfit_mass%s" % self.systs, "F")

        self.out.branch("HHKinFit_mass%s" % self.systs, "F")
        self.out.branch("HHKinFit_chi2%s" % self.systs, "F")

        self.out.branch("VBFjj_mass%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaEta%s" % self.systs, "F")
        self.out.branch("VBFjj_deltaPhi%s" % self.systs, "F")
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        htt_svfit = Object(event, "Htt_svfit")
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")
        jets = Collection(event, "Jet")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus) #FIXME FROM HERE
        bjet1, bjet2, bjet1_tlv, bjet2_tlv = self.get_bjets(event, jets)
        vbfjet1, vbfjet2, vbfjet1_tlv, vbfjet2_tlv = self.get_vbfjets(event, jets)
        met, met_tlv = self.get_met(event)

        met_tv = ROOT.TVector2(met_tlv.Px(), met_tlv.Py())

        cov = ROOT.TMatrixD(2, 2)
        cov[0][0] = met.covXX
        cov[0][1] = met.covXY
        cov[1][0] = met.covXY
        cov[1][1] = met.covYY

        # HH TLorentzVector
        htt_tlv = dau1_tlv + dau2_tlv
        hbb_tlv = bjet1_tlv + bjet2_tlv
        hh_tlv = htt_tlv + hbb_tlv

        # HH (with htt_svfit) TLorentzVector
        htt_svfit_tlv = ROOT.TLorentzVector()
        htt_svfit_tlv.SetPtEtaPhiM(
            eval("htt_svfit.pt%s" % self.systs), 
            eval("htt_svfit.eta%s" % self.systs), 
            eval("htt_svfit.phi%s" % self.systs), 
            eval("htt_svfit.mass%s" % self.systs),
        )
        hh_svfit_tlv = hbb_tlv + htt_svfit_tlv

        # HH Kin. Fit
        # constructor https://github.com/bvormwald/HHKinFit2/blob/CMSSWversion/HHKinFit2Scenarios/src/HHKinFitMasterHeavyHiggs.cpp#L27
        # kinFit = ROOT.HHKinFit2.HHKinFitMasterHeavyHiggs(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit = ROOT.HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv, met_tv, cov)
        kinFit.addHypo(125, 125);
        results = kinFit.fit()
        HHKinFit_mass = results[0]
        HHKinFit_chi2 = results[1]

        # VBF jet composition
        if not vbfjet1:
            VBFjj_mass = -999.
            VBFjj_deltaEta = -999.
            VBFjj_deltaPhi = -999.
        else:
            VBFjj_tlv = vbfjet1_tlv + vbfjet2_tlv
            VBFjj_mass = VBFjj_tlv.M()
            VBFjj_deltaEta = vbfjet1_tlv.Eta() - vbfjet2_tlv.Eta()
            VBFjj_deltaPhi = vbfjet1_tlv.DeltaPhi(vbfjet2_tlv)

        self.out.fillBranch("Hbb_pt%s" % self.systs, hbb_tlv.Pt())
        self.out.fillBranch("Hbb_eta%s" % self.systs, hbb_tlv.Eta())
        self.out.fillBranch("Hbb_phi%s" % self.systs, hbb_tlv.Phi())
        self.out.fillBranch("Hbb_mass%s" % self.systs, hbb_tlv.M())

        self.out.fillBranch("Htt_pt%s" % self.systs, htt_tlv.Pt())
        self.out.fillBranch("Htt_eta%s" % self.systs, htt_tlv.Eta())
        self.out.fillBranch("Htt_phi%s" % self.systs, htt_tlv.Phi())
        self.out.fillBranch("Htt_mass%s" % self.systs, htt_tlv.M())

        self.out.fillBranch("HH_pt%s" % self.systs, hh_tlv.Pt())
        self.out.fillBranch("HH_eta%s" % self.systs, hh_tlv.Eta())
        self.out.fillBranch("HH_phi%s" % self.systs, hh_tlv.Phi())
        self.out.fillBranch("HH_mass%s" % self.systs, hh_tlv.M())

        self.out.fillBranch("HH_svfit_pt%s" % self.systs, hh_svfit_tlv.Pt())
        self.out.fillBranch("HH_svfit_eta%s" % self.systs, hh_svfit_tlv.Eta())
        self.out.fillBranch("HH_svfit_phi%s" % self.systs, hh_svfit_tlv.Phi())
        self.out.fillBranch("HH_svfit_mass%s" % self.systs, hh_svfit_tlv.M())

        self.out.fillBranch("HHKinFit_mass%s" % self.systs, HHKinFit_mass)
        self.out.fillBranch("HHKinFit_chi2%s" % self.systs, HHKinFit_chi2)

        self.out.fillBranch("VBFjj_mass%s" % self.systs, VBFjj_mass)
        self.out.fillBranch("VBFjj_deltaEta%s" % self.systs, VBFjj_deltaEta)
        self.out.fillBranch("VBFjj_deltaPhi%s" % self.systs, VBFjj_deltaPhi)

        return True


class HHKinFitRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.isMC = kwargs.pop("isMC")
        self.AnalysisType = AnalysisType
        super(HHKinFitRDFProducer, self).__init__(isMC=self.isMC, *args, **kwargs)

        if not os.getenv("_KINFIT"):
            os.environ["_KINFIT"] = "kinfit"

            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")
            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/HHKinFitInterface.h".format(base))
            ROOT.gInterpreter.Declare("""
                #include <TLorentzVector.h>
                #include "TVector.h"
                #include "TMatrixD.h"
                #include "TMath.h"

                using Vfloat = const ROOT::RVec<float>&;
                ROOT::RVec<double> compute_hhkinfit(
                        bool isMC, int pairType, int bjet1_index, int bjet2_index,
                        float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
                        float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
                        Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                        Vfloat jet_resolution, float met_pt, float met_phi,
                        float met_covXX, float met_covXY, float met_covYY, int target1, int target2) {
                    if (bjet1_index < 0 || bjet2_index < 0)
                        return {-1., -1.}; // KinFit needs 2 b jets

                    auto bjet1_tlv = TLorentzVector();
                    auto bjet2_tlv = TLorentzVector();
                    auto dau1_tlv = TLorentzVector();
                    auto dau2_tlv = TLorentzVector();

                    dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
                    dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                    bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                        jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                    bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                        jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));

                    auto met_tv = TVector2(met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi));
                    auto cov = TMatrixD(2, 2);
                    cov[0][0] = met_covXX;
                    cov[0][1] = met_covXY;
                    cov[1][0] = met_covXY;
                    cov[1][1] = met_covYY;

                    double resolution1, resolution2;
                    if (isMC) {
                        resolution1 = bjet1_tlv.E() * jet_resolution.at(bjet1_index);
                        resolution2 = bjet2_tlv.E() * jet_resolution.at(bjet2_index);
                    } else {
                        resolution1 = -1;
                        resolution2 = -1;
                    }

                    auto kinFit = HHKinFitInterface(bjet1_tlv, bjet2_tlv, dau1_tlv, dau2_tlv,
                        met_tv, cov, resolution1, resolution2);
                    kinFit.addHypo(target1, target2);
                    return kinFit.fit();
                }
            """)

    def run(self, df):
        isMC = ("true" if self.isMC else "false")
        jet_resolution = "jet_pt_resolution"
        if not self.isMC:
            jet_resolution = "Jet_eta"  # placeholder

        # DEBUG
        if not self.AnalysisType:
            pp = "HH"; target1 = 125; target2 = 125
            print(" ### INFO: Running KinFit with default option for HH analysis")
        else:
            print(" ### INFO: Running KinFit with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      pp = "ZZ"; target1 = 91; target2 = 91
            elif self.AnalysisType == "Zbb_Htautau":    pp = "ZH"; target1 = 125; target2 = 91 # 1 is tautau, 2 is bb
            elif self.AnalysisType == "Ztautau_Hbb":    pp = "ZH"; target1 = 91; target2 = 125 

        branches = ["%sKinFit_mass%s" % (pp, self.systs), "%sKinFit_chi2%s" % (pp, self.systs)]
        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        df = df.Define(f"hhkinfit_result{self.systs}",
            f"compute_hhkinfit({isMC}, pairType, bjet1_JetIdx, bjet2_JetIdx, "
                f"dau1_pt{self.lep_syst}, dau1_eta, dau1_phi, dau1_mass{self.lep_syst}, "
                f"dau2_pt{self.lep_syst}, dau2_eta, dau2_phi, dau2_mass{self.lep_syst}, "
                f"Jet_pt{self.jet_syst}, Jet_eta, Jet_phi, Jet_mass{self.jet_syst}, {jet_resolution},"
                f"MET{self.met_smear_tag}_pt{self.muon_syst}, MET{self.met_smear_tag}_phi{self.muon_syst}, MET_covXX, MET_covXY, MET_covYY, "
                f"{target1}, {target2})"
            ).Define("%sKinFit_mass%s" % (pp, self.systs), "hhkinfit_result%s[0]" % self.systs
            ).Define("%sKinFit_chi2%s" % (pp, self.systs), "hhkinfit_result%s[1]" % self.systs)
        return df, branches


class HHVarRDFProducer(JetLepMetSyst):
    def __init__(self, AnalysisType, *args, **kwargs):
        self.AnalysisType = AnalysisType
        super(HHVarRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("_HHVAR_%s" % self.systs):
            os.environ["_HHVAR_%s" % self.systs] = "hhvar"
            ROOT.gInterpreter.Declare("""
                #include <TLorentzVector.h>
                using Vfloat = const ROOT::RVec<float>&;
                ROOT::RVec<double> get_hh_features%s(int pairType, int isBoosted,
                        int bjet1_index, int bjet2_index, int fatjet_index,
                        int vbfjet1_index, int vbfjet2_index,
                        float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
                        float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
                        Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass,
                        Vfloat fatjet_pt, Vfloat fatjet_eta, Vfloat fatjet_phi, Vfloat fatjet_mass,
                        float met_pt, float met_phi,
                        double Htt_svfit_pt, double Htt_svfit_eta,
                        double Htt_svfit_phi, double Htt_svfit_mass) {

                    auto dau1_tlv = TLorentzVector();
                    auto dau2_tlv = TLorentzVector();

                    dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
                    dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                    auto htt_tlv = dau1_tlv + dau2_tlv;
                    
                    
                    auto hbb_tlv = TLorentzVector();
                    if (!isBoosted && bjet1_index >= 0 && bjet2_index >= 0) {
                        auto bjet1_tlv = TLorentzVector();
                        auto bjet2_tlv = TLorentzVector();
                        bjet1_tlv.SetPtEtaPhiM(jet_pt.at(bjet1_index), jet_eta.at(bjet1_index),
                            jet_phi.at(bjet1_index), jet_mass.at(bjet1_index));
                        bjet2_tlv.SetPtEtaPhiM(jet_pt.at(bjet2_index), jet_eta.at(bjet2_index),
                            jet_phi.at(bjet2_index), jet_mass.at(bjet2_index));
                        hbb_tlv = bjet1_tlv + bjet2_tlv;
                    }
                    else if (isBoosted && fatjet_index >= 0) {
                        hbb_tlv.SetPtEtaPhiM(fatjet_pt.at(fatjet_index), fatjet_eta.at(fatjet_index),
                            fatjet_phi.at(fatjet_index), fatjet_mass.at(fatjet_index));
                    }
                    auto hh_tlv = htt_tlv + hbb_tlv;

                    auto met_tlv = TLorentzVector();
                    met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

                    auto htt_met_tlv = htt_tlv + met_tlv;

                    double hh_svfit_pt = -999., hh_svfit_eta = -999.;
                    double hh_svfit_phi = -999., hh_svfit_mass = -999.;
                    if (Htt_svfit_pt > 0) {
                        auto htt_svfit_tlv = TLorentzVector();
                        htt_svfit_tlv.SetPtEtaPhiM(Htt_svfit_pt, Htt_svfit_eta,
                            Htt_svfit_phi, Htt_svfit_mass);
                        auto hh_svfit_tlv = htt_svfit_tlv + hbb_tlv;
                        hh_svfit_pt = hh_svfit_tlv.Pt();
                        hh_svfit_eta = hh_svfit_tlv.Eta();
                        hh_svfit_phi = hh_svfit_tlv.Phi();
                        hh_svfit_mass = hh_svfit_tlv.M();
                    }
                    double vbfjj_mass = -999., vbfjj_deltaEta = -999., vbfjj_deltaPhi = -999.;
                    if (vbfjet1_index >= 0) {
                        auto vbfjet1_tlv = TLorentzVector();
                        auto vbfjet2_tlv = TLorentzVector();
                        vbfjet1_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet1_index), jet_eta.at(vbfjet1_index),
                            jet_phi.at(vbfjet1_index), jet_mass.at(vbfjet1_index));
                        vbfjet2_tlv.SetPtEtaPhiM(jet_pt.at(vbfjet2_index), jet_eta.at(vbfjet2_index),
                            jet_phi.at(vbfjet2_index), jet_mass.at(vbfjet2_index));
                        auto vbfjj_tlv = vbfjet1_tlv + vbfjet2_tlv;
                        vbfjj_mass = vbfjj_tlv.M();
                        vbfjj_deltaEta = vbfjet1_tlv.Eta() - vbfjet2_tlv.Eta();
                        vbfjj_deltaPhi = vbfjet1_tlv.Phi() - vbfjet2_tlv.Phi();
                    }
                    return {
                        hbb_tlv.Pt(), hbb_tlv.Eta(), hbb_tlv.Phi(), hbb_tlv.M(),
                        htt_tlv.Pt(), htt_tlv.Eta(), htt_tlv.Phi(), htt_tlv.M(),
                        htt_met_tlv.Pt(), htt_met_tlv.Eta(), htt_met_tlv.Phi(), htt_met_tlv.M(),
                        hh_tlv.Pt(), hh_tlv.Eta(), hh_tlv.Phi(), hh_tlv.M(),
                        hh_svfit_pt, hh_svfit_eta, hh_svfit_phi, hh_svfit_mass,
                        vbfjj_mass, vbfjj_deltaEta, vbfjj_deltaPhi
                    };
                }
            """ % self.systs)

    def run(self, df):

        # DEBUG
        if not self.AnalysisType:
            p_b = "H"; p_t = "H"; pp = "HH"; p_sv = "H"
            print(" ### INFO: Running HHVar with default option for HH analysis")
        else:
            print(" ### INFO: Running HHVar with AnalysisType = {}".format(self.AnalysisType))
            if self.AnalysisType == "Zbb_Ztautau":      p_b = "Z"; p_t = "Z"; pp = "ZZ"; p_sv = "X"
            elif self.AnalysisType == "Zbb_Htautau":    p_b = "Z"; p_t = "H"; pp = "ZH"; p_sv = "X"
            elif self.AnalysisType == "Ztautau_Hbb":    p_b = "H"; p_t = "Z"; pp = "ZH"; p_sv = "X"

        features = (f"{p_b}bb_pt{self.jet_syst},{p_b}bb_eta,{p_b}bb_phi,{p_b}bb_mass{self.jet_syst},"
            f"{p_t}tt_pt{self.lep_syst},{p_t}tt_eta,{p_t}tt_phi,{p_t}tt_mass{self.lep_syst},"
            # basically every systematic is propagated to MET, thus to SVFit & KinFit. also they will change eta & phi probably
            f"{p_t}tt_met_pt{self.systs},{p_t}tt_met_eta{self.systs},{p_t}tt_met_phi{self.systs},{p_t}tt_met_mass{self.systs},"
            f"{pp}_pt{self.systs},{pp}_eta{self.systs},{pp}_phi{self.systs},{pp}_mass{self.systs},"
            f"{pp}_svfit_pt{self.systs},{pp}_svfit_eta{self.systs},{pp}_svfit_phi{self.systs},{pp}_svfit_mass{self.systs},"
            f"VBFjj_mass{self.systs},VBFjj_deltaEta,VBFjj_deltaPhi".format(p_b, p_t, pp, self.systs))
        features = list(features.split(","))
        all_branches = df.GetColumnNames()
        if features[0] in all_branches:
            return df, []

        df = df.Define("hhfeatures%s" % self.systs, (f"get_hh_features{self.systs}(pairType, isBoosted, "
            f"bjet1_JetIdx, bjet2_JetIdx, fatjet_JetIdx, VBFjet1_JetIdx, VBFjet2_JetIdx, "
            f"dau1_pt{self.lep_syst}, dau1_eta, dau1_phi, dau1_mass{self.lep_syst}, "
            f"dau2_pt{self.lep_syst}, dau2_eta, dau2_phi, dau2_mass{self.lep_syst}, "
            f"Jet_pt{self.jet_syst}, Jet_eta, Jet_phi, Jet_mass{self.jet_syst}, "
            f"FatJet_pt{self.jet_syst}, FatJet_eta, FatJet_phi, FatJet_mass{self.jet_syst}, "
            f"MET{self.met_smear_tag}_pt{self.met_syst}, MET{self.met_smear_tag}_phi{self.met_syst}, "
            f"{p_sv}tt_svfit_pt{self.systs}, {p_sv}tt_svfit_eta{self.systs}, {p_sv}tt_svfit_phi{self.systs}, {p_sv}tt_svfit_mass{self.systs})"
        ))

        for ifeat, feature in enumerate(features):
            df = df.Define(feature, "hhfeatures%s[%s]" % (self.systs, ifeat))
        return df, features


def HH(**kwargs):
    return lambda: HHProducer(**kwargs)


def HHKinFitRDF(*args, **kwargs):
    # The output of HHKinFitRDF changes when targeting H or Z, accornding to the 
    # mass hypothesis for the fit: at the moment there is a bug with the hypothesis 
    # To bypass it, it's necessary to manually change mh1 and mh2 values inside 
    # HHKinFit2/HHKinFit2Scenarios/interface/HHKinFitMasterHeavyHiggs.h
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 1 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHKinFitRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)


def HHVarRDF(*args, **kwargs):
    """ Computes variables relating to H/Z and to HH/ZZ/ZH (pt, eta, mass, phi) 
    For bb, for non-boosted sums the 2 AK4 jets 4-vectors. For boosted uses the AK8 jet directly
    """
    # The output of HHVarRDF is not affected by H or Z, but the output features
    # are called in different ways according to the ZZ or HH analysis
    AnalysisType = kwargs.pop("AnalysisType", False)

    # print("### DEBUG 2 : AnalysisType = {}".format(AnalysisType))
    return lambda: HHVarRDFProducer(AnalysisType=AnalysisType, *args, **kwargs)
