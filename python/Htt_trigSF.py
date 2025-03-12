import os
from copy import deepcopy as copy

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from analysis_tools.utils import import_root, getContentHisto3D
from Base.Modules.baseModules import JetLepMetModule, JetLepMetSyst

ROOT = import_root()

class Htt_trigSFProducer(JetLepMetModule):
    def __init__(self, isMC, year, *args, **kwargs):
        super(Htt_trigSFProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC
        self.year = year

        if year == 2016:
            self.mutau_pt_th1 = 23.
            self.mutau_pt_th2 = 25.
        elif year in [2017, 2018]:
            self.mutau_pt_th1 = 25.
            self.mutau_pt_th2 = 32.
            self.etau_pt_th1 = 33.
            self.etau_pt_th2 = 35.

        base = "{}/{}/src/HTT-utilities/LepEffInterface".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
        base_tau = "{}/{}/src/TauAnalysisTools/TauTriggerSFs".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

        if not os.getenv("_Htt_trigSF"):
            os.environ["_Htt_trigSF"] = "_Htt_trigSF"
            if "/libHTT-utilitiesLepEffInterface.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libHTT-utilitiesLepEffInterface.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/ScaleFactor.h".format(base))

            if "/libTauAnalysisToolsTauTriggerSFs.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libTauAnalysisToolsTauTriggerSFs.so")
            ROOT.gROOT.ProcessLine(".L {}/interface/SFProvider.h".format(base_tau))

        if self.isMC:
            self.eTrgSF = ROOT.ScaleFactor()
            self.eTauTrgSF = ROOT.ScaleFactor()
            self.muTrgSF = ROOT.ScaleFactor()
            self.muTauTrgSF = ROOT.ScaleFactor()

            if self.year == 2016:
                self.eTrgSF.init_ScaleFactor(
                    "{}/data/Electron/Run2016/Electron_Run2016_legacy_Ele25.root".format(base))
                self.eTauTrgSF = None  # No ele cross triggers in 2016
                self.muTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2016/Muon_Run2016_legacy_IsoMu22.root".format(base))
                self.muTauTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2016/Muon_Mu19leg_2016BtoH_eff.root".format(base))
                self.tauTrgSF_ditau = ROOT.SFProvider(
                    "{}/data/2016_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "ditau", "Medium")
                self.tauTrgSF_mutau = ROOT.SFProvider(
                    "{}/data/2016_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "mutau", "Medium")
                self.tauTrgSF_etau = None
                self.tauTrgSF_vbf = None
                self.jetTrgSF_vbf = None
            elif self.year == 2017:
                self.eTrgSF.init_ScaleFactor(
                    "{}/data/Electron/Run2017/Electron_Ele32orEle35_fix.root".format(base))
                self.eTauTrgSF.init_ScaleFactor(
                    "{}/data/Electron/Run2017/Electron_EleTau_Ele24_fix.root".format(base))
                self.muTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2017/Muon_IsoMu24orIsoMu27.root".format(base))
                self.muTauTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2017/Muon_MuTau_IsoMu20.root".format(base))
                self.tauTrgSF_ditau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "ditau", "Medium")
                self.tauTrgSF_mutau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "mutau", "Medium")
                self.tauTrgSF_etau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "etau", "Medium")
                self.tauTrgSF_vbf = ROOT.tau_trigger.SFProvider(
                    "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "ditauvbf", "Medium")
                f = "{}/data/2017_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                tf = ROOT.TFile.Open(f)
                self.jetTrgSF_vbf = copy(tf.Get("SF_mjj_pT1_pT2"))
            elif self.year == 2018:
                self.eTrgSF.init_ScaleFactor(
                    "{}/data/Electron/Run2018/Electron_Run2018_Ele32orEle35.root".format(base))
                self.eTauTrgSF.init_ScaleFactor(
                    "{}/data/Electron/Run2018/Electron_Run2018_Ele24.root".format(base))
                self.muTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2018/Muon_Run2018_IsoMu24orIsoMu27.root".format(base))
                self.muTauTrgSF.init_ScaleFactor(
                    "{}/data/Muon/Run2018/Muon_Run2018_IsoMu20.root".format(base))
                self.tauTrgSF_ditau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "ditau", "Medium")
                self.tauTrgSF_mutau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "mutau", "Medium")
                self.tauTrgSF_etau = ROOT.tau_trigger.SFProvider(
                    "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "etau", "Medium")
                self.tauTrgSF_vbf = ROOT.tau_trigger.SFProvider(
                    "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau),
                    "ditauvbf", "Medium")
                f = "{}/data/2018_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                tf = ROOT.TFile.Open(f)
                self.jetTrgSF_vbf = copy(tf.Get("SF_mjj_pT1_pT2"))
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('trigSF', 'F')
        self.out.branch('trigSF_single', 'F')
        self.out.branch('trigSF_cross', 'F')
        self.out.branch('trigSF_muUp', 'F')
        self.out.branch('trigSF_muDown', 'F')
        self.out.branch('trigSF_eleUp', 'F')
        self.out.branch('trigSF_eleDown', 'F')
        self.out.branch('trigSF_DM0Up', 'F')
        self.out.branch('trigSF_DM1Up', 'F')
        self.out.branch('trigSF_DM10Up', 'F')
        self.out.branch('trigSF_DM11Up', 'F')
        self.out.branch('trigSF_DM0Down', 'F')
        self.out.branch('trigSF_DM1Down', 'F')
        self.out.branch('trigSF_DM10Down', 'F')
        self.out.branch('trigSF_DM11Down', 'F')
        self.out.branch('trigSF_vbfjetUp', 'F')
        self.out.branch('trigSF_vbfjetDown', 'F')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        dau1, dau2, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        decaymodes = [0, 1, 10, 11]

        if event.pairType == 0 and self.isMC:
            # eta region covered both by cross-trigger and single lepton trigger
            if abs(dau2_tlv.Eta()) < 2.1:
                passSingle = (1 if dau1_tlv.Pt() >= self.mutau_pt_th1 else 0)
                passCross = (1 if dau2_tlv.Pt() >= self.mutau_pt_th2 else 0)

                # lepton trigger
                SFL_Data = self.muTrgSF.get_EfficiencyData(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_MC = self.muTrgSF.get_EfficiencyMC(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_Data_Err = self.muTrgSF.get_EfficiencyDataError(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_MC_Err = self.muTrgSF.get_EfficiencyMCError(dau1_tlv.Pt(), dau1_tlv.Eta())

                SFL_Data = [SFL_Data - 1. * SFL_Data_Err, SFL_Data, SFL_Data + 1. * SFL_Data_Err]
                SFL_MC = [SFL_MC - 1. * SFL_MC_Err, SFL_MC, SFL_MC + 1. * SFL_MC_Err]

                # cross trigger
                # mu leg
                SFl_Data = self.muTauTrgSF.get_EfficiencyData(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_MC = self.muTauTrgSF.get_EfficiencyMC(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_Data_Err = self.muTauTrgSF.get_EfficiencyDataError(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_MC_Err = self.muTauTrgSF.get_EfficiencyMCError(dau1_tlv.Pt(), dau1_tlv.Eta())

                SFl_Data = [SFl_Data - 1. * SFl_Data_Err, SFl_Data, SFl_Data + 1. * SFl_Data_Err]
                SFl_MC = [SFl_MC - 1. * SFl_MC_Err, SFl_MC, SFl_MC + 1. * SFl_MC_Err]

                # tau leg
                SFtau_Data = self.tauTrgSF_mutau.getEfficiencyData(dau2_tlv.Pt(), dau2.decayMode, 0)
                SFtau_MC = self.tauTrgSF_mutau.getEfficiencyMC(dau2_tlv.Pt(), dau2.decayMode, 0)

                # efficiencies
                Eff_Data_mu =\
                    [passSingle * SFL_Data[i] - passCross * passSingle * min([SFl_Data[i], SFL_Data[i]])
                    * SFtau_Data + passCross * SFl_Data[i] * SFtau_Data for i in range(3)]
                Eff_MC_mu =\
                    [passSingle * SFL_MC[i] - passCross * passSingle * min([SFl_MC[i], SFL_MC[i]])
                    * SFtau_MC + passCross * SFl_MC[i] * SFtau_MC for i in range(3)]

                SFtau_Data_tauup = [SFtau_Data for i in range(len(decaymodes))]
                SFtau_Data_taudown = [SFtau_Data for i in range(len(decaymodes))]
                SFtau_MC_tauup = [SFtau_MC for i in range(len(decaymodes))]
                SFtau_MC_taudown = [SFtau_MC for i in range(len(decaymodes))]
                Eff_Data_tauup = [Eff_Data_mu for i in range(len(decaymodes))]
                Eff_Data_taudown = [Eff_Data_mu for i in range(len(decaymodes))]
                Eff_MC_tauup = [Eff_MC_mu for i in range(len(decaymodes))]
                Eff_MC_taudown = [Eff_MC_mu for i in range(len(decaymodes))]

                for idm, dm in enumerate(decaymodes):
                    if dm == int(dau2.decayMode):
                        SFtau_Data_tauup[idm] = self.tauTrgSF_mutau.getEfficiencyData(
                            dau2_tlv.Pt(), dau2.decayMode, 1)
                        SFtau_Data_taudown[idm] = self.tauTrgSF_mutau.getEfficiencyData(
                            dau2_tlv.Pt(), dau2.decayMode, -1)
                        SFtau_MC_tauup[idm] = self.tauTrgSF_mutau.getEfficiencyMC(
                            dau2_tlv.Pt(), dau2.decayMode, 1)
                        SFtau_MC_taudown[idm] = self.tauTrgSF_mutau.getEfficiencyMC(
                            dau2_tlv.Pt(), dau2.decayMode, -1)

                Eff_Data_tauup = [
                    passSingle * SFL_Data[1] - passCross * passSingle * min([SFl_Data[1], SFL_Data[1]])
                        * SFtau_Data_tauup[idm] + passCross * SFl_Data[1] * SFtau_Data_tauup[idm]
                    for idm in range(len(decaymodes))]
                Eff_Data_taudown = [
                    passSingle * SFL_Data[1] - passCross * passSingle * min([SFl_Data[1], SFL_Data[1]])
                        * SFtau_Data_taudown[idm] + passCross * SFl_Data[1] * SFtau_Data_taudown[idm]
                    for idm in range(len(decaymodes))]
                Eff_MC_tauup = [
                    passSingle * SFL_MC[1] - passCross * passSingle * min([SFl_MC[1], SFL_MC[1]])
                        * SFtau_MC_tauup[idm] + passCross * SFl_MC[1] * SFtau_MC_tauup[idm]
                    for idm in range(len(decaymodes))]
                Eff_MC_taudown = [
                    passSingle * SFL_MC[1] - passCross * passSingle * min([SFl_MC[1], SFL_MC[1]])
                        * SFtau_MC_taudown[idm] + passCross * SFl_MC[1] * SFtau_MC_taudown[idm]
                    for idm in range(len(decaymodes))]

                trigSF_mu = [Eff_Data_mu[i] / Eff_MC_mu[i] for i in range(len(Eff_Data_mu))]
                trigSF_tauup = [Eff_Data_tauup[i] / Eff_MC_tauup[i] for i in range(len(Eff_Data_tauup))]
                trigSF_taudown = [Eff_Data_taudown[i] / Eff_MC_taudown[i] for i in range(len(Eff_Data_taudown))]
                # trig SF for analysis only with cross-trigger
                SFl = self.muTauTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFtau = self.tauTrgSF_mutau.getSF(dau2_tlv.Pt(), dau2.decayMode, 0)
                trigSF_cross = SFl * SFtau

            else:
                SF = self.muTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())
                SF_error = self.muTrgSF.get_ScaleFactorError(dau1_tlv.Pt(), dau1_tlv.Eta())
                trigSF_mu = [SF - 1. * SF_error, SF, SF + 1. * SF_error]
                trigSF_cross = SF
                trigSF_tauup = [SF for i in range(len(decaymodes))]
                trigSF_taudown = [SF for i in range(len(decaymodes))]

            # trig SF for analysis only with single-mu trigger
            trigSF_single = self.muTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())

            self.out.fillBranch('trigSF', trigSF_mu[1])
            self.out.fillBranch('trigSF_single', trigSF_single)
            self.out.fillBranch('trigSF_cross', trigSF_cross)
            self.out.fillBranch('trigSF_muUp', trigSF_mu[2])
            self.out.fillBranch('trigSF_muDown', trigSF_mu[0])
            self.out.fillBranch('trigSF_eleUp', trigSF_mu[1])
            self.out.fillBranch('trigSF_eleDown', trigSF_mu[1])
            self.out.fillBranch('trigSF_DM0Up', trigSF_tauup[0])
            self.out.fillBranch('trigSF_DM1Up', trigSF_tauup[1])
            self.out.fillBranch('trigSF_DM10Up', trigSF_tauup[2])
            self.out.fillBranch('trigSF_DM11Up', trigSF_tauup[3])
            self.out.fillBranch('trigSF_DM0Down', trigSF_taudown[0])
            self.out.fillBranch('trigSF_DM1Down', trigSF_taudown[1])
            self.out.fillBranch('trigSF_DM10Down', trigSF_taudown[2])
            self.out.fillBranch('trigSF_DM11Down', trigSF_taudown[3])
            self.out.fillBranch('trigSF_vbfjetUp', trigSF_mu[1])
            self.out.fillBranch('trigSF_vbfjetDown', trigSF_mu[1])
            return True

        elif event.pairType == 1 and self.isMC:
            # eta region covered both by cross-trigger and single lepton trigger
            # in 2016 there is no cross etau trigger
            if abs(dau2_tlv.Eta()) < 2.1 and self.year != 2016:
                passSingle = (1 if dau1_tlv.Pt() >= self.etau_pt_th1 else 0)
                passCross = (1 if dau2_tlv.Pt() >= self.etau_pt_th2 else 0)

                # lepton trigger
                SFL_Data = self.eTrgSF.get_EfficiencyData(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_MC = self.eTrgSF.get_EfficiencyMC(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_Data_Err = self.eTrgSF.get_EfficiencyDataError(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFL_MC_Err = self.eTrgSF.get_EfficiencyMCError(dau1_tlv.Pt(), dau1_tlv.Eta())

                SFL_Data = [SFL_Data - 1. * SFL_Data_Err, SFL_Data, SFL_Data + 1. * SFL_Data_Err]
                SFL_MC = [SFL_MC - 1. * SFL_MC_Err, SFL_MC, SFL_MC + 1. * SFL_MC_Err]

                # cross trigger
                # e leg
                SFl_Data = self.eTauTrgSF.get_EfficiencyData(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_MC = self.eTauTrgSF.get_EfficiencyMC(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_Data_Err = self.eTauTrgSF.get_EfficiencyDataError(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFl_MC_Err = self.eTauTrgSF.get_EfficiencyMCError(dau1_tlv.Pt(), dau1_tlv.Eta())

                SFl_Data = [SFl_Data - 1. * SFl_Data_Err, SFl_Data, SFl_Data + 1. * SFl_Data_Err]
                SFl_MC = [SFl_MC - 1. * SFl_MC_Err, SFl_MC, SFl_MC + 1. * SFl_MC_Err]

                # tau leg
                SFtau_Data = self.tauTrgSF_etau.getEfficiencyData(dau2_tlv.Pt(), dau2.decayMode, 0)
                SFtau_MC = self.tauTrgSF_etau.getEfficiencyMC(dau2_tlv.Pt(), dau2.decayMode, 0)

                # efficiencies            
                Eff_Data_e =\
                    [passSingle * SFL_Data[i]
                        - passCross * passSingle * min([SFl_Data[i], SFL_Data[i]]) * SFtau_Data
                        + passCross * SFl_Data[i] * SFtau_Data for i in range(3)]
                Eff_MC_e =\
                    [passSingle * SFL_MC[i]
                        - passCross * passSingle * min([SFl_MC[i], SFL_MC[i]]) * SFtau_MC
                        + passCross * SFl_MC[i] * SFtau_MC for i in range(3)]

                SFtau_Data_tauup = [SFtau_Data for i in range(len(decaymodes))]
                SFtau_Data_taudown = [SFtau_Data for i in range(len(decaymodes))]
                SFtau_MC_tauup = [SFtau_MC for i in range(len(decaymodes))]
                SFtau_MC_taudown = [SFtau_MC for i in range(len(decaymodes))]
                Eff_Data_tauup = [Eff_Data_e for i in range(len(decaymodes))]
                Eff_Data_taudown = [Eff_Data_e for i in range(len(decaymodes))]
                Eff_MC_tauup = [Eff_MC_e for i in range(len(decaymodes))]
                Eff_MC_taudown = [Eff_MC_e for i in range(len(decaymodes))]

                for idm, dm in enumerate(decaymodes):
                    if dm == int(dau2.decayMode):
                        SFtau_Data_tauup[idm] = self.tauTrgSF_etau.getEfficiencyData(
                            dau2_tlv.Pt(), dau2.decayMode, 1)
                        SFtau_Data_taudown[idm] = self.tauTrgSF_etau.getEfficiencyData(
                            dau2_tlv.Pt(), dau2.decayMode, -1)
                        SFtau_MC_tauup[idm] = self.tauTrgSF_etau.getEfficiencyMC(
                            dau2_tlv.Pt(), dau2.decayMode, 1)
                        SFtau_MC_taudown[idm] = self.tauTrgSF_etau.getEfficiencyMC(
                            dau2_tlv.Pt(), dau2.decayMode, -1)

                Eff_Data_tauup = [
                    passSingle * SFL_Data[1]
                        - passCross * passSingle * min([SFl_Data[1], SFL_Data[1]]) * SFtau_Data_tauup[idm]
                        + passCross * SFl_Data[1] * SFtau_Data_tauup[idm]
                    for idm in range(len(decaymodes))]
                Eff_Data_taudown = [
                    passSingle * SFL_Data[1]
                        - passCross * passSingle * min([SFl_Data[1], SFL_Data[1]]) * SFtau_Data_taudown[idm]
                        + passCross * SFl_Data[1] * SFtau_Data_taudown[idm]
                    for idm in range(len(decaymodes))]
                Eff_MC_tauup = [
                    passSingle * SFL_MC[1]
                        - passCross * passSingle * min([SFl_MC[1], SFL_MC[1]]) * SFtau_MC_tauup[idm]
                        + passCross * SFl_MC[1] * SFtau_MC_tauup[idm]
                    for idm in range(len(decaymodes))]
                Eff_MC_taudown = [
                    passSingle * SFL_MC[1]
                        - passCross * passSingle * min([SFl_MC[1], SFL_MC[1]]) * SFtau_MC_taudown[idm]
                        + passCross * SFl_MC[1] * SFtau_MC_taudown[idm]
                    for idm in range(len(decaymodes))]

                trigSF_e = [Eff_Data_e[i] / Eff_MC_e[i] for i in range(len(Eff_Data_e))]
                trigSF_tauup = [Eff_Data_tauup[i] / Eff_MC_tauup[i] for i in range(len(Eff_Data_tauup))]
                trigSF_taudown = [Eff_Data_taudown[i] / Eff_MC_taudown[i]
                    for i in range(len(Eff_Data_taudown))]
                # trig SF for analysis only with cross-trigger
                SFl = self.eTauTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())
                SFtau = self.tauTrgSF_etau.getSF(dau2_tlv.Pt(), dau2.decayMode, 0)
                trigSF_cross = SFl * SFtau

            else:
                SF = self.eTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())
                SF_error = self.eTrgSF.get_ScaleFactorError(dau1_tlv.Pt(), dau1_tlv.Eta())
                trigSF_e = [SF - 1. * SF_error, SF, SF + 1. * SF_error]
                trigSF_cross = SF
                trigSF_tauup = [SF for i in range(len(decaymodes))]
                trigSF_taudown = [SF for i in range(len(decaymodes))]

            # trig SF for analysis only with single-e trigger
            trigSF_single = self.eTrgSF.get_ScaleFactor(dau1_tlv.Pt(), dau1_tlv.Eta())

            self.out.fillBranch('trigSF', trigSF_e[1])
            self.out.fillBranch('trigSF_single', trigSF_single)
            self.out.fillBranch('trigSF_cross', trigSF_cross)
            self.out.fillBranch('trigSF_muUp', trigSF_e[1])
            self.out.fillBranch('trigSF_muDown', trigSF_e[1])
            self.out.fillBranch('trigSF_eleUp', trigSF_e[2])
            self.out.fillBranch('trigSF_eleDown', trigSF_e[0])
            self.out.fillBranch('trigSF_DM0Up', trigSF_tauup[0])
            self.out.fillBranch('trigSF_DM1Up', trigSF_tauup[1])
            self.out.fillBranch('trigSF_DM10Up', trigSF_tauup[2])
            self.out.fillBranch('trigSF_DM11Up', trigSF_tauup[3])
            self.out.fillBranch('trigSF_DM0Down', trigSF_taudown[0])
            self.out.fillBranch('trigSF_DM1Down', trigSF_taudown[1])
            self.out.fillBranch('trigSF_DM10Down', trigSF_taudown[2])
            self.out.fillBranch('trigSF_DM11Down', trigSF_taudown[3])
            self.out.fillBranch('trigSF_vbfjetUp', trigSF_e[1])
            self.out.fillBranch('trigSF_vbfjetDown', trigSF_e[1])
            return True

        elif event.pairType == 2 and self.isMC and event.isVBFtrigger == 0:
            
            SF1 = self.tauTrgSF_ditau.getSF(dau1_tlv.Pt(), dau1.decayMode, 0)
            SF2 = self.tauTrgSF_ditau.getSF(dau2_tlv.Pt(), dau2.decayMode, 0)

            SF1_tauup = [SF1 for i in range(len(decaymodes))]
            SF1_taudown = [SF1 for i in range(len(decaymodes))]
            SF2_tauup = [SF2 for i in range(len(decaymodes))]
            SF2_taudown = [SF2 for i in range(len(decaymodes))]

            for idm, dm in enumerate(decaymodes):
                if dm == int(dau1.decayMode):
                    SF1_tauup[idm] = self.tauTrgSF_ditau.getSF(dau1_tlv.Pt(), dau1.decayMode, 1)
                    SF1_taudown[idm] = self.tauTrgSF_ditau.getSF(dau1_tlv.Pt(), dau1.decayMode, -1)
                if dm == int(dau2.decayMode):
                    SF2_tauup[idm] = self.tauTrgSF_ditau.getSF(dau2_tlv.Pt(), dau2.decayMode, 1)
                    SF2_taudown[idm] = self.tauTrgSF_ditau.getSF(dau2_tlv.Pt(), dau2.decayMode, -1)

            trigSF = SF1 * SF2
            trigSF_tauup = [SF1_tauup[i] * SF2_tauup[i] for i in range(len(SF1_tauup))]
            trigSF_taudown = [SF1_taudown[i] * SF2_taudown[i] for i in range(len(SF1_taudown))]

            self.out.fillBranch('trigSF', trigSF)
            self.out.fillBranch('trigSF_single', trigSF)
            self.out.fillBranch('trigSF_cross', trigSF)
            self.out.fillBranch('trigSF_muUp', trigSF)
            self.out.fillBranch('trigSF_muDown', trigSF)
            self.out.fillBranch('trigSF_eleUp', trigSF)
            self.out.fillBranch('trigSF_eleDown', trigSF)
            self.out.fillBranch('trigSF_DM0Up', trigSF_tauup[0])
            self.out.fillBranch('trigSF_DM1Up', trigSF_tauup[1])
            self.out.fillBranch('trigSF_DM10Up', trigSF_tauup[2])
            self.out.fillBranch('trigSF_DM11Up', trigSF_tauup[3])
            self.out.fillBranch('trigSF_DM0Down', trigSF_taudown[0])
            self.out.fillBranch('trigSF_DM1Down', trigSF_taudown[1])
            self.out.fillBranch('trigSF_DM10Down', trigSF_taudown[2])
            self.out.fillBranch('trigSF_DM11Down', trigSF_taudown[3])
            self.out.fillBranch('trigSF_vbfjetUp', trigSF)
            self.out.fillBranch('trigSF_vbfjetDown', trigSF)
            return True

        elif event.pairType == 2 and self.isMC and self.year != 2016:
            jets = Collection(event, "Jet")
            vbfjet1, vbfjet2, vbfjet1_tlv, vbfjet2_tlv = self.get_vbfjets(event, jets)
            
            
            if vbfjet1 and vbfjet2:
                vbfjj_mass = (vbfjet1_tlv + vbfjet2_tlv).M()
                if (vbfjet1_tlv.Pt() > 140 and vbfjet2_tlv.Pt() > 60 and vbfjj_mass > 800
                        and dau1_tlv.Pt() > 25 and dau2_tlv.Pt() > 25
                        and (dau1_tlv.Pt() <= 40 or dau2_tlv.Pt() <= 40)):
                    # jet leg sf
                    jetSF = getContentHisto3D(
                        self.jetTrgSF_vbf, vbfjj_mass, vbfjet1_tlv.Pt(), vbfjet2_tlv.Pt(), 0)
                    jetSFerror = getContentHisto3D(
                        self.jetTrgSF_vbf, vbfjj_mass, vbfjet1_tlv.Pt(), vbfjet2_tlv.Pt(), 1)
                    # tau leg sf
                    SF1 = self.tauTrgSF_vbf.getSF(dau1_tlv.Pt(), dau1.decayMode, 0)
                    SF2 = self.tauTrgSF_vbf.getSF(dau2_tlv.Pt(), dau2.decayMode, 0)

                    SF1_tauup = [SF1 for i in range(len(decaymodes))]
                    SF1_taudown = [SF1 for i in range(len(decaymodes))]
                    SF2_tauup = [SF2 for i in range(len(decaymodes))]
                    SF2_taudown = [SF2 for i in range(len(decaymodes))]

                    for idm, dm in enumerate(decaymodes):
                        if dm == int(dau1.decayMode):
                            SF1_tauup[idm] = self.tauTrgSF_vbf.getSF(dau1_tlv.Pt(), dau1.decayMode, 1)
                            SF1_taudown[idm] = self.tauTrgSF_vbf.getSF(dau1_tlv.Pt(), dau1.decayMode, -1)
                        if dm == int(dau2.decayMode):
                            SF2_tauup[idm] = self.tauTrgSF_vbf.getSF(dau2_tlv.Pt(), dau2.decayMode, 1)
                            SF2_taudown[idm] = self.tauTrgSF_vbf.getSF(dau2_tlv.Pt(), dau2.decayMode, -1)

                    trigSF_vbfjet = [
                        (jetSF - jetSFerror) * SF1 * SF2,
                        jetSF * SF1 * SF2,
                        (jetSF + jetSFerror) * SF1 * SF2
                    ]

                    trigSF_tauup = [jetSF * SF1_tauup[i] * SF2_tauup[i] for i in range(len(SF1_tauup))]
                    trigSF_taudown = [jetSF * SF1_taudown[i] * SF2_taudown[i] for i in range(len(SF1_taudown))]

                    self.out.fillBranch('trigSF', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_single', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_cross', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_muUp', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_muDown', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_eleUp', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_eleDown', trigSF_vbfjet[1])
                    self.out.fillBranch('trigSF_DM0Up', trigSF_tauup[0])
                    self.out.fillBranch('trigSF_DM1Up', trigSF_tauup[1])
                    self.out.fillBranch('trigSF_DM10Up', trigSF_tauup[2])
                    self.out.fillBranch('trigSF_DM11Up', trigSF_tauup[3])
                    self.out.fillBranch('trigSF_DM0Down', trigSF_taudown[0])
                    self.out.fillBranch('trigSF_DM1Down', trigSF_taudown[1])
                    self.out.fillBranch('trigSF_DM10Down', trigSF_taudown[2])
                    self.out.fillBranch('trigSF_DM11Down', trigSF_taudown[3])
                    self.out.fillBranch('trigSF_vbfjetUp', trigSF_vbfjet[2])
                    self.out.fillBranch('trigSF_vbfjetDown', trigSF_vbfjet[0])
                    return True

        elif not self.isMC:
            self.out.fillBranch('trigSF', 1.)
            self.out.fillBranch('trigSF_single', 1.)
            self.out.fillBranch('trigSF_cross', 1.)
            self.out.fillBranch('trigSF_muUp', 1.)
            self.out.fillBranch('trigSF_muDown', 1.)
            self.out.fillBranch('trigSF_eleUp', 1.)
            self.out.fillBranch('trigSF_eleDown', 1.)
            self.out.fillBranch('trigSF_DM0Up', 1.)
            self.out.fillBranch('trigSF_DM1Up', 1.)
            self.out.fillBranch('trigSF_DM10Up', 1.)
            self.out.fillBranch('trigSF_DM11Up', 1.)
            self.out.fillBranch('trigSF_DM0Down', 1.)
            self.out.fillBranch('trigSF_DM1Down', 1.)
            self.out.fillBranch('trigSF_DM10Down', 1.)
            self.out.fillBranch('trigSF_DM11Down', 1.)
            self.out.fillBranch('trigSF_vbfjetUp', 1.)
            self.out.fillBranch('trigSF_vbfjetDown', 1.)
            return True
        else:
            self.out.fillBranch('trigSF', 0.)
            self.out.fillBranch('trigSF_single', 0.)
            self.out.fillBranch('trigSF_cross', 0.)
            self.out.fillBranch('trigSF_muUp', 0.)
            self.out.fillBranch('trigSF_muDown', 0.)
            self.out.fillBranch('trigSF_eleUp', 0.)
            self.out.fillBranch('trigSF_eleDown', 0.)
            self.out.fillBranch('trigSF_DM0Up', 0.)
            self.out.fillBranch('trigSF_DM1Up', 0.)
            self.out.fillBranch('trigSF_DM10Up', 0.)
            self.out.fillBranch('trigSF_DM11Up', 0.)
            self.out.fillBranch('trigSF_DM0Down', 0.)
            self.out.fillBranch('trigSF_DM1Down', 0.)
            self.out.fillBranch('trigSF_DM10Down', 0.)
            self.out.fillBranch('trigSF_DM11Down', 0.)
            self.out.fillBranch('trigSF_vbfjetUp', 0.)
            self.out.fillBranch('trigSF_vbfjetDown', 0.)
            return True


def Htt_trigSF(**kwargs):
    return lambda: Htt_trigSFProducer(**kwargs)


class Htt_trigSFRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        self.isUL = kwargs.pop("isUL")
        self.ispreVFP = kwargs.pop("ispreVFP")
        self.computeSystematics = kwargs.pop("computeSystematics", True)
        self.met_smear_tag_data = kwargs.pop("met_smear_tag_data")
        super(Htt_trigSFRDFProducer, self).__init__(*args, **kwargs)

        # Thresholds to determine which trigger fired (cross or single lepton) 
        if self.year == 2016:
            mutau_pt_th1 = 23. #singleMuon trigger threshold (used with offline). 
            mutau_pt_th2 = 25. #cross mutau trigger tau leg threshold (used with offline)
            etau_pt_th1 = -1. # no cross etau trigger for 2016
            etau_pt_th2 = -1.
        elif self.year in [2017, 2018]:
            mutau_pt_th1 = 25. 
            mutau_pt_th2 = 32.
            etau_pt_th1 = 33.
            etau_pt_th2 = 35. # cross etau trigger tau leg threshold (offline)

        if self.isMC:

            if not os.getenv("_Htt_trigSF"):
                os.environ["_Htt_trigSF"] = "_Htt_trigSF"
                
                base = "{}/{}/src/Tools/Tools".format(
                        os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                    ROOT.gSystem.Load("libToolsTools.so")
                ROOT.gROOT.ProcessLine(".L {}/interface/Htt_trigSFinterface.h".format(base))

                base = "{}/{}/src/HTT-utilities/LepEffInterface".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                base_tau = "{}/{}/src/TauAnalysisTools/TauTriggerSFs".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                base_egamma = "{}/{}/src/HTT-utilities/EgammaPogSF_UL".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                base_mu = "{}/{}/src/HTT-utilities/MuPogSF_UL".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                # SFs from gammagamma->tautau analysis (Cecile Caillol), same as HH->bbtt resonant. Includes singleElectron and cross-e, cross-mu trigger SFs
                base_cross = "{}/{}/src/HTT-utilities/trigSFs_UL_eleMu".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))

                if not self.isUL:
                    if self.year == 2016:
                        eTrgSF = "{}/data/Electron/Run2016/Electron_Run2016_legacy_Ele25.root".format(base)
                        # using 2017 as dummy
                        eTauTrgSF = "{}/data/Electron/Run2017/Electron_EleTau_Ele24_fix.root".format(base)
                        muTrgSF = "{}/data/Muon/Run2016/Muon_Run2016_legacy_IsoMu22.root".format(base)
                        muTauTrgSF = "{}/data/Muon/Run2016/Muon_Mu19leg_2016BtoH_eff.root".format(base)
                        tauTrgSF_ditau = "{}/data/2016_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_mutau = "{}/data/2016_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # using 2017 as dummy
                        tauTrgSF_etau = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # using 2017 as dummy
                        tauTrgSF_vbf = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # using 2017 as dummy
                        jetTrgSF_vbf = "{}/data/2017_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                    elif self.year == 2017:
                        eTrgSF = "{}/data/Electron/Run2017/Electron_Ele32orEle35_fix.root".format(base)
                        eTauTrgSF = "{}/data/Electron/Run2017/Electron_EleTau_Ele24_fix.root".format(base)
                        muTrgSF = "{}/data/Muon/Run2017/Muon_IsoMu24orIsoMu27.root".format(base)
                        muTauTrgSF = "{}/data/Muon/Run2017/Muon_MuTau_IsoMu20.root".format(base)
                        tauTrgSF_ditau = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_mutau = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_etau = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_vbf = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        jetTrgSF_vbf = "{}/data/2017_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                    elif self.year == 2018:
                        eTrgSF = "{}/data/Electron/Run2018/Electron_Run2018_Ele32orEle35.root".format(base)
                        eTauTrgSF = "{}/data/Electron/Run2018/Electron_Run2018_Ele24.root".format(base)
                        muTrgSF = "{}/data/Muon/Run2018/Muon_Run2018_IsoMu24orIsoMu27.root".format(base)
                        muTauTrgSF = "{}/data/Muon/Run2018/Muon_Run2018_IsoMu20.root".format(base)
                        tauTrgSF_ditau = "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_mutau = "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_etau = "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_vbf = "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        jetTrgSF_vbf = "{}/data/2018_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                else:
                    if self.year == 2016:
                        if self.ispreVFP:   
                            prefix = "pre"
                            eTrgSF = "{}/sf_el_2016{}_HLTEle25.root".format(base_cross, prefix); eTrgName = ""; eTrgBool = "true"
                        else:               
                            prefix = "post"
                            eTrgSF = "{}/sf_el_2016{}_HLTEle25_fix.root".format(base_cross, prefix); eTrgName = ""; eTrgBool = "true"
                        eTauTrgSF = "{}/sf_el_2017_HLTEle24Tau30.root".format(base_cross); eTauTrgName = ""; eTauTrgBool = "true" # using 2017 as dummy
                        muTauTrgSF = "{}/sf_mu_2016{}_HLTMu20Tau27.root".format(base_cross, prefix); muTauTrgName = ""; muTauTrgBool = "true"
                        if self.ispreVFP:
                            muTrgSF = "{}/Run2/UL/2016_preVFP/2016_preVFP_trigger/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root".format(base_mu); 
                            muTrgName = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt"; 
                            muTrgBool = ""
                            tauTrgSF_ditau = "{}/data/2016ULpreVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                            tauTrgSF_mutau = "{}/data/2016ULpreVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                            tauTrgSF_etau  = "{}/data/2016ULpreVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        else:
                            muTrgSF = "{}/Run2/UL/2016_postVFP/2016_postVFP_trigger/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root".format(base_mu); 
                            muTrgName = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt"; 
                            muTrgBool = ""
                            tauTrgSF_ditau = "{}/data/2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                            tauTrgSF_mutau = "{}/data/2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                            tauTrgSF_etau  = "{}/data/2016ULpostVFP_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_vbf = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau) # using 2017 as dummy
                        jetTrgSF_vbf = "{}/data/2017_VBFHTauTauTrigger_JetLegs.root".format(base_tau) # using 2017 as dummy
                    elif self.year == 2017:
                        eTrgSF = "{}/sf_el_2017_HLTEle32.root".format(base_cross); eTrgName = ""; eTrgBool = "true"
                        eTauTrgSF = "{}/sf_el_2017_HLTEle24Tau30.root".format(base_cross); eTauTrgName = ""; eTauTrgBool = "true"
                        muTrgSF = "{}/Run2/UL/2017/2017_trigger/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root".format(base_mu); muTrgName = "NUM_IsoMu27_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt"; muTrgBool = ""
                        muTauTrgSF = "{}/sf_mu_2017_HLTMu20Tau27.root".format(base_cross); muTauTrgName = ""; muTauTrgBool = "true"
                        tauTrgSF_ditau = "{}/data/2017UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_mutau = "{}/data/2017UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_etau = "{}/data/2017UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # use Legacy file, since Sfs for DiTau OR VBF are not yet available
                        tauTrgSF_vbf = "{}/data/2017_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # for vbf using legacy while being computed
                        jetTrgSF_vbf = "{}/data/2017_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                    elif self.year == 2018:
                        # TH2D x=pT, y=eta -> init_EG_ScaleFactor. Checked, works
                        eTrgSF = "{}/sf_el_2018_HLTEle32.root".format(base_cross); eTrgName = ""; eTrgBool = "true"
                        # same format
                        eTauTrgSF = "{}/sf_el_2018_HLTEle24Tau30.root".format(base_cross); eTauTrgName = ""; eTauTrgBool = "true"
                        # TH2D x=eta, y=pT -> init_ScaleFactor(file, HistoName) (absolute eta). Checked, works
                        muTrgSF = "{}/Run2/UL/2018/2018_trigger/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root".format(base_mu); muTrgName = "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt"; muTrgBool = ""
                        # TH2D x=pT, y=eta -> init_EG_ScaleFactor
                        muTauTrgSF = "{}/sf_mu_2018_HLTMu20Tau27.root".format(base_cross); muTauTrgName = ""; muTauTrgBool = "true"
                        tauTrgSF_ditau = "{}/data/2018UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_mutau = "{}/data/2018UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        tauTrgSF_etau = "{}/data/2018UL_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # use Legacy file, since Sfs for DiTau OR VBF are not yet available
                        tauTrgSF_vbf = "{}/data/2018_tauTriggerEff_DeepTau2017v2p1.root".format(base_tau)
                        # for jet using legacy while being computed
                        jetTrgSF_vbf = "{}/data/2018_VBFHTauTauTrigger_JetLegs.root".format(base_tau)
                    else:
                        raise ValueError("TauTriggerSFs not implemented yet for 2017")

                # TrgSF: path to the TriggerSF file
                # TrgName: if empty string "" the contructor is init_EG_ScaleFactor, else init_ScaleFactor (default legacy is ZMass)
                # TrgBool: used inside init_EG_ScaleFactor
                ROOT.gInterpreter.Declare("""
                    auto Htt_trigSF = Htt_trigSFinterface(%s, %s, %s, %s, %s,
                        "%s", "%s", "%s", 
                        "%s", "%s", "%s", 
                        "%s", "%s", "%s", 
                        "%s", "%s", "%s", 
                        "%s", "%s", "%s", "%s", "%s");
                """ % (int(self.year), mutau_pt_th1, mutau_pt_th2, etau_pt_th1, etau_pt_th2,
                    eTrgSF, eTrgName, eTrgBool, 
                    eTauTrgSF, eTauTrgName, eTauTrgBool,
                    muTrgSF, muTrgName, muTrgBool, 
                    muTauTrgSF, muTauTrgName, muTauTrgBool,
                    tauTrgSF_ditau, tauTrgSF_mutau, tauTrgSF_etau, tauTrgSF_vbf, jetTrgSF_vbf))

               
                if self.year == 2016 and self.ispreVFP:
                    met_tag = "2016APV"
                else:
                    met_tag = str(self.year) 
                MET_trig_sf_file = f"{os.getenv('CMT_CMSSW_BASE')}/{os.getenv('CMT_CMSSW_VERSION')}/src/Tools/Tools/data/met_trig_sf/{met_tag}_MetTriggerSFs.root"
                ROOT.gInterpreter.Declare(f"""
                    auto MET_trigSF_interface_obj = MET_trigSF_interface("{MET_trig_sf_file}");
                    """"""
                    MET_trigSF_interface::trigSF_result computeMETTriggerSF(float MET_pt, bool isBoostedTau, int pairType) {
                        if (isBoostedTau) {
                            if (MET_pt < 180 && pairType >= 0) {
                                throw std::runtime_error("Trying to get MET trigger SF for boostedTau with MET_pt<180GeV whilst pairType>=0&isBoostedTau");
                                
                            } else if (MET_pt >= 140) {
                                return  MET_trigSF_interface_obj.getSF(MET_pt);
                            }
                            // The case (MET_pt<180 && pairType<0) is when there is no pairType filter and we have boostedTaus that pass MET trigger but not MET pt cut
                            // in which case we can use SFs if MET_pt >140 (it will be outside of trigger plateau though)
                        } 
                        return MET_trigSF_interface::trigSF_result{1., 1., 1.};
                        
                    }
                """ )


    def run(self, df):
        if not self.isMC:
            return df, []
        
        df = df.Define("htt_trigsf", (
            "isBoostedTau ? Htt_trigSFresult() : "
            "Htt_trigSF.get_scale_factors(pairType, isVBFtrigger, "
            f"dau1_decayMode, dau1_pt{self.lep_syst}, dau1_eta, "
            f"dau2_decayMode, dau2_pt{self.lep_syst}, dau2_eta)"
            ))
        met_smear_tag = self.met_smear_tag if self.isMC else self.met_smear_tag_data
        df = df.Define("MET_trigsf_res", f"computeMETTriggerSF(MET{met_smear_tag}_pt{self.met_syst}, isBoostedTau, pairType)")

        # nominal
        df = df.Define("trigSF", "isBoostedTau ? MET_trigsf_res.SF : htt_trigsf.trigSF")

        df = df.Define("passSingle", "htt_trigsf.passSingle")
        df = df.Define("passCross", "htt_trigsf.passCross")

        if not self.computeSystematics or not self.systematic_is_central:
            return df, ["trigSF", "passSingle", "passCross"]
        
        # lepton trigger variations
        branches_lepTriggers = [
            'muUp', 'muDown', 'eleUp', 'eleDown',
            'DM0Up', 'DM1Up', 'DM10Up', 'DM11Up', 
            'DM0Down', 'DM1Down', 'DM10Down', 'DM11Down']
        for syst_branch in branches_lepTriggers:
            # syst variations of lepton triggers
            # if boostedTau, just set to the nominal SF from MET trigger
            df = df.Define(f"trigSF_{syst_branch}", f"isBoostedTau ? trigSF :  htt_trigsf.{syst_branch}")

        # MET trigger variations
        for met_trig_syst in ["statup", "statdown"]:
            # if not boostedTau, set to nominal SF from lepton triggers
            df = df.Define(f"trigSF_met_{met_trig_syst}", f"isBoostedTau ? MET_trigsf_res.SF_{met_trig_syst} : trigSF")

        return df, ["trigSF", "passSingle", "passCross"] + [f"trigSF_{b}" for b in branches_lepTriggers] + ["trigSF_met_statup", "trigSF_met_statdown"]


def Htt_trigSFRDF(**kwargs):
    """ Computes trigger scale factors for lepton trigger and MET triggers
    Depends on lepton systematics (electron S&S, tau energy scale) for lepton triggers.
    Depends on MET variations for MET triggers (so basically every systematic)
    Chooses MET or lepton trigger depending on isBoostedTau

    Branches : 
     - trigSF : nominal trigger SF (either lepton trigger or MET trigger depending on isBoostedTau)
     - trigSF_single : SF of the single lepton trigger
     - trigSF_cross : SF of the cross lepton trigger 
     - trigSF_muUp/muDown/eleUp/eleDown/DM0Up/DM1Up/DM10Up/DM11Up/DM0Down/DM1Down/DM10Down/DM11Down : systematic variations of lepton triggers. Set to "trigSF" in case of MET trigger
     - trigSF_met_statup/statdown : syst variations of MET trigger. Set to "trigSF" in case of lepton trigger
    """
    return lambda: Htt_trigSFRDFProducer(**kwargs)
