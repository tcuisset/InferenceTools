# based on https://github.com/LLRCMS/KLUBAnalysis/blob/VBF_legacy/src/PuJetIdSF.cc

import os
from copy import deepcopy as copy

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from analysis_tools.utils import import_root, getContentHisto2D
from Base.Modules.baseModules import JetLepMetModule, DummyModule, JetLepMetSyst

ROOT = import_root()

class PUjetID_SFProducer(JetLepMetModule):
    def __init__(self, year, *args, **kwargs):
        super(PUjetID_SFProducer, self).__init__(*args, **kwargs)
        self.year = year

        f_eff = ROOT.TFile.Open(os.path.expandvars(
            "$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/h2_eff_mc_%s_L.root" % self.year))
        f_eff_sf = ROOT.TFile.Open(os.path.expandvars(
            "$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/h2_eff_sf_%s_L.root" % self.year))
        f_mistag = ROOT.TFile.Open(os.path.expandvars(
            "$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/h2_mistag_mc_%s_L.root" % self.year))
        f_mistag_sf = ROOT.TFile.Open(os.path.expandvars(
            "$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/h2_mistag_sf_%s_L.root" % self.year))
        f_sf_err = ROOT.TFile.Open(os.path.expandvars(
            "$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/scalefactorsPUID_81Xtraining.root"))

        self.h_eff_ = copy(f_eff.Get("h2_eff_mc%s_L" % self.year))
        self.h_eff_sf_ = copy(f_eff_sf.Get("h2_eff_sf%s_L" % self.year))
        self.h_eff_sf_err_ = copy(f_sf_err.Get("h2_eff_sf%s_L_Systuncty" % self.year))
        self.h_mistag_ = copy(f_mistag.Get("h2_mistag_mc%s_L" % self.year))
        self.h_mistag_sf_ = copy(f_mistag_sf.Get("h2_mistag_sf%s_L" % self.year))
        self.h_mistag_sf_err_ = copy(f_sf_err.Get("h2_mistag_sf%s_L_Systuncty" % self.year))

        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('PUjetID_SF', 'F')
        self.out.branch('PUjetID_SF_up', 'F')
        self.out.branch('PUjetID_SF_down', 'F')
        self.out.branch('PUjetID_SF_eff_up', 'F')
        self.out.branch('PUjetID_SF_eff_down', 'F')
        self.out.branch('PUjetID_SF_mistag_up', 'F')
        self.out.branch('PUjetID_SF_mistag_down', 'F')
        self.out.branch('PUjetID_SF_eff_eta_s2p5_up', 'F')
        self.out.branch('PUjetID_SF_eff_eta_s2p5_down', 'F')
        self.out.branch('PUjetID_SF_mistag_eta_s2p5_up', 'F')
        self.out.branch('PUjetID_SF_mistag_eta_s2p5_down', 'F')
        self.out.branch('PUjetID_SF_eff_eta_l2p5_up', 'F')
        self.out.branch('PUjetID_SF_eff_eta_l2p5_down', 'F')
        self.out.branch('PUjetID_SF_mistag_eta_l2p5_up', 'F')
        self.out.branch('PUjetID_SF_mistag_eta_l2p5_down', 'F')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def get_eff_sf_and_error(self, isReal, pt, eta):
        if pt < 20.:
            pt = 20.
        elif pt > 50.:
            pt = 50.

        if isReal:
            return (getContentHisto2D(self.h_eff_, pt, eta),
                getContentHisto2D(self.h_eff_sf_, pt, eta),
                getContentHisto2D(self.h_eff_sf_err_, pt, eta))
        else:
            return (getContentHisto2D(self.h_mistag_, pt, eta),
                getContentHisto2D(self.h_mistag_sf_, pt, eta),
                getContentHisto2D(self.h_mistag_sf_err_, pt, eta))

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        P_MC = 1.
        P_DATA = 1.
        P_DATA_up = 1.
        P_DATA_down = 1.
        P_DATA_effic_up = 1.
        P_DATA_effic_down = 1.
        P_DATA_mistag_up = 1.
        P_DATA_mistag_down = 1.
        P_DATA_effic_eta_s2p5_up = 1.
        P_DATA_effic_eta_s2p5_down = 1.
        P_DATA_effic_eta_l2p5_up = 1.
        P_DATA_effic_eta_l2p5_down = 1.
        P_DATA_mistag_eta_s2p5_up = 1.
        P_DATA_mistag_eta_s2p5_down = 1.
        P_DATA_mistag_eta_l2p5_up = 1.
        P_DATA_mistag_eta_l2p5_down = 1.

        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")
        jets = Collection(event, "Jet")
        genjets = Collection(event, "GenJet")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)

        for jet in jets:
            if jet.jetId < 2:   
                continue
            jet_tlv = ROOT.TLorentzVector()
            jet_tlv.SetPtEtaPhiM(
                eval("jet.pt%s" % self.jet_syst),
                jet.eta,
                jet.phi,
                eval("jet.mass%s" % self.jet_syst)
            )
            if jet_tlv.Pt() < 20 or jet_tlv.Pt() > 50 or abs(jet_tlv.Eta()) > 4.7:
                continue
            if jet_tlv.DeltaR(dau1_tlv) < 0.5 or jet_tlv.DeltaR(dau2_tlv) < 0.5:
                continue

            # noisy jet removal for 2017
            # https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorkingLegacyRun2#Jets
            if self.year == 2017 and abs(tlv_jet.Eta()) > 2.65 and abs(tlv_jet.Eta()) < 3.139:
                continue

            isRealJet = False
            for genjet in genjets:
                genjet_tlv = ROOT.TLorentzVector()
                genjet_tlv.SetPtEtaPhiM(genjet.pt, genjet.eta, genjet.phi, genjet.mass)
                if jet_tlv.DeltaR(genjet_tlv) < 0.4:
                    isRealJet = True
                    break

            passPUjetIDLoose = jet.puId >= 4 or eval("jet.pt%s" % self.jet_syst) > 50
            eff, sf, sf_err = self.get_eff_sf_and_error(isRealJet, jet_tlv.Pt(), jet_tlv.Eta())
            sf_up = min([max([0, sf + sf_err]), 5])
            sf_down = min([max([0, sf - sf_err]), 5])

            if passPUjetIDLoose:
                P_MC *= eff;
                P_DATA *= sf*eff;
                P_DATA_up *= sf_up*eff;
                P_DATA_down *= sf_down*eff;
                P_DATA_effic_up *= sf_up*eff;
                P_DATA_effic_down *= sf_down*eff;
                P_DATA_mistag_up *= sf*eff;  # true jet --> use nominal SF for mistag
                P_DATA_mistag_down *= sf*eff;  # true jet --> use nominal SF for mistag

                if abs(jet_tlv.Eta()) <= 2.5:
                    P_DATA_effic_eta_s2p5_up *= sf_up * eff
                    P_DATA_effic_eta_s2p5_down *= sf_down * eff
                    P_DATA_effic_eta_l2p5_up *= sf * eff
                    P_DATA_effic_eta_l2p5_down *= sf * eff
                    P_DATA_mistag_eta_s2p5_up *= sf * eff
                    P_DATA_mistag_eta_s2p5_down *= sf * eff
                    P_DATA_mistag_eta_l2p5_up *= sf * eff
                    P_DATA_mistag_eta_l2p5_down *= sf * eff
                else:
                    P_DATA_effic_eta_s2p5_up *= sf * eff
                    P_DATA_effic_eta_s2p5_down *= sf * eff
                    P_DATA_effic_eta_l2p5_up *= sf_up * eff
                    P_DATA_effic_eta_l2p5_down *= sf_down * eff
                    P_DATA_mistag_eta_s2p5_up *= sf * eff
                    P_DATA_mistag_eta_s2p5_down *= sf * eff
                    P_DATA_mistag_eta_l2p5_up *= sf * eff
                    P_DATA_mistag_eta_l2p5_down *= sf * eff
            else:
                P_MC *= (1. - eff);
                P_DATA *= (1. - sf * eff);
                P_DATA_up *= (1. - sf_up * eff);
                P_DATA_down *= (1. - sf_down * eff);
                P_DATA_effic_up *= (1. - sf * eff);  # fake jet --> use nominal SF for effic
                P_DATA_effic_down *= (1. - sf * eff);  # fake jet --> use nominal SF for effic
                P_DATA_mistag_up *= (1. - sf_up * eff);
                P_DATA_mistag_down *= (1. - sf_down * eff);

                if abs(jet_tlv.Eta()) <= 2.5:
                    P_DATA_effic_eta_s2p5_up *= (1 - sf * eff)
                    P_DATA_effic_eta_s2p5_down *= (1 - sf * eff)
                    P_DATA_effic_eta_l2p5_up *= (1 - sf * eff)
                    P_DATA_effic_eta_l2p5_down *= (1 - sf * eff)
                    P_DATA_mistag_eta_s2p5_up *= (1 - sf_up * eff)
                    P_DATA_mistag_eta_s2p5_down *= (1 - sf_down * eff)
                    P_DATA_mistag_eta_l2p5_up *= (1 - sf * eff)
                    P_DATA_mistag_eta_l2p5_down *= (1 - sf * eff)
                else:
                    P_DATA_effic_eta_s2p5_up *= (1 - sf * eff)
                    P_DATA_effic_eta_s2p5_down *= (1 - sf * eff)
                    P_DATA_effic_eta_l2p5_up *= (1 - sf * eff)
                    P_DATA_effic_eta_l2p5_down *= (1 - sf * eff)
                    P_DATA_mistag_eta_s2p5_up *= (1 - sf * eff)
                    P_DATA_mistag_eta_s2p5_down *= (1 - sf * eff)
                    P_DATA_mistag_eta_l2p5_up *= (1 - sf_up * eff)
                    P_DATA_mistag_eta_l2p5_down *= (1 - sf_down * eff)

        self.out.fillBranch('PUjetID_SF', P_DATA / P_MC)
        self.out.fillBranch('PUjetID_SF_up', P_DATA_up / P_MC)
        self.out.fillBranch('PUjetID_SF_down', P_DATA_down / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_up', P_DATA_effic_up / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_down', P_DATA_effic_down / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_up', P_DATA_mistag_up / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_down', P_DATA_mistag_down / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_eta_s2p5_up', P_DATA_effic_eta_s2p5_up / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_eta_s2p5_down', P_DATA_effic_eta_s2p5_down / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_eta_s2p5_up', P_DATA_effic_eta_l2p5_up / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_eta_s2p5_down', P_DATA_effic_eta_l2p5_down / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_eta_l2p5_up', P_DATA_mistag_eta_s2p5_up / P_MC)
        self.out.fillBranch('PUjetID_SF_eff_eta_l2p5_down', P_DATA_mistag_eta_s2p5_down / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_eta_l2p5_up', P_DATA_mistag_eta_l2p5_up / P_MC)
        self.out.fillBranch('PUjetID_SF_mistag_eta_l2p5_down', P_DATA_mistag_eta_l2p5_down / P_MC)
        return True


def PUjetID_SF(**kwargs):
    isMC = kwargs.pop("isMC")
    year = kwargs.pop("year")
    if not isMC:
        return lambda: DummyModule(**kwargs)
    return lambda: PUjetID_SFProducer(year, **kwargs)


class PUjetID_SFRDFProducer(JetLepMetSyst):
    def __init__(self, year, *args, **kwargs):
        super(PUjetID_SFRDFProducer, self).__init__(*args, **kwargs)
        self.year = year
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")
            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/PUjetID_SFinterface.h".format(base))

            folder_path = os.path.expandvars("$CMSSW_BASE/src/Tools/Tools/python/pujetid_sf/")
            ROOT.gInterpreter.Declare("""
                auto PUjetID_SF = PUjetID_SFinterface(%s, "%s");
            """ % (int(self.year), folder_path))
        
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using VInt = const ROOT::RVec<int>&;
                std::vector<double> get_pujetid_sf (
                    Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, Vfloat Jet_mass,
                    VInt Jet_puId, Vfloat Jet_jetId,
                    Vfloat GenJet_pt, Vfloat GenJet_eta, Vfloat GenJet_phi, Vfloat GenJet_mass,
                    int pairType, int dau1_index, int dau2_index,
                    Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_eta, Vfloat electron_phi, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_eta, Vfloat tau_phi, Vfloat tau_mass
                )
                {
                    float dau1_pt, dau1_eta, dau1_phi, dau1_mass, dau2_pt, dau2_eta, dau2_phi, dau2_mass;
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
                    } else {
                        dau1_pt = -999.;
                        dau1_eta = -999.;
                        dau1_phi = -999.;
                        dau1_mass = -999.;
                    }
                    dau2_pt = tau_pt.at(dau2_index);
                    dau2_eta = tau_eta.at(dau2_index);
                    dau2_phi = tau_phi.at(dau2_index);
                    dau2_mass = tau_mass.at(dau2_index);

                    return PUjetID_SF.get_pu_weights(
                        Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_jetId, Jet_puId,
                        GenJet_pt, GenJet_eta,GenJet_phi, GenJet_mass,
                        dau1_pt, dau1_eta, dau1_phi, dau1_mass,
                        dau2_pt, dau2_eta, dau2_phi, dau2_mass);
                }
            """)

    def run(self, df):
        if not self.isMC:
            return df, []
        branches = ['PUjetID_SF2', 'PUjetID_SF_up2', 'PUjetID_SF_down2',
            'PUjetID_SF_eff_up2', 'PUjetID_SF_eff_down2',
            'PUjetID_SF_mistag_up2', 'PUjetID_SF_mistag_down2',
            'PUjetID_SF_eff_eta_s2p5_up2', 'PUjetID_SF_eff_eta_s2p5_down2',
            'PUjetID_SF_mistag_eta_s2p5_up2', 'PUjetID_SF_mistag_eta_s2p5_down2',
            'PUjetID_SF_eff_eta_l2p5_up2', 'PUjetID_SF_eff_eta_l2p5_down2',
            'PUjetID_SF_mistag_eta_l2p5_up2', 'PUjetID_SF_mistag_eta_l2p5_down2']

        df = df.Define("pujetid_sf_results", "get_pujetid_sf("
            "Jet_pt{3}, Jet_eta, Jet_phi, Jet_mass{3}, "
            "Jet_puId, Jet_jetId, "
            "GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, "
            "pairType, dau1_index, dau2_index, "
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}"
        ")".format(
            self.muon_syst, self.electron_syst, self.tau_syst, self.jet_syst))
        for ib, branch in enumerate(branches):
            df = df.Define(branch, "pujetid_sf_results[%s]" % ib)
        return df, branches


def PUjetID_SFRDF(**kwargs):
    year = kwargs.pop("year")
    return lambda: PUjetID_SFRDFProducer(year, **kwargs)

