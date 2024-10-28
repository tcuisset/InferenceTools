import os

from analysis_tools.utils import import_root

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from Tools.Tools.tau_utils import LeptonTauPair, TriggerChecker, lepton_veto
from Base.Modules.baseModules import JetLepMetSyst, JetLepMetModule

ROOT = import_root()

class HHLeptonProducer(JetLepMetModule):
    def __init__(self, isMC, year, runEra, *args, **kwargs):
        super(HHLeptonProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC
        self.year = year
        self.runEra = runEra
        self.trigger_checker = TriggerChecker(year)

        if self.year == 2016:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
                "HLT_IsoTkMu22", "HLT_IsoTkMu22_eta2p1"]
            self.trigger_checker.mutau_crosstriggers = ["HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
                "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1"]
            self.trigger_checker.etau_triggers = ["HLT_Ele25_eta2p1_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = []
            if not self.isMC:
                if self.runEra != "H":
                    self.trigger_checker.tautau_triggers = [
                        "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg"]
                else:
                    self.trigger_checker.tautau_triggers = [
                        "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg"]
            else:
                self.trigger_checker.tautau_triggers = [
                    "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
                    "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg"]
            self.trigger_checker.vbf_triggers = []

        elif self.year == 2017:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu24", "HLT_IsoMu27"]
            self.trigger_checker.mutau_crosstriggers = [
                "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1"]
            self.trigger_checker.etau_triggers = ["HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = [
                "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1"]
            self.trigger_checker.tautau_triggers = [
                "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg"]
            if self.runEra in ["D", "E", "F"] or self.isMC:
                self.trigger_checker.vbf_triggers = [
                    "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg"]
            else:
                self.trigger_checker.vbf_triggers = []

        elif self.year == 2018:
            self.trigger_checker.mutau_triggers = ["HLT_IsoMu24", "HLT_IsoMu27"]
            # lista = lambda e: ([1, 1] if e == 0 else [2, 2])
            self.trigger_checker.mutau_crosstriggers = lambda e: (
                ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"])
            self.trigger_checker.etau_triggers = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf"]
            self.trigger_checker.etau_crosstriggers = lambda e: (
                ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"])
            self.trigger_checker.tautau_triggers = lambda e: (
                ["HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg"])
            self.trigger_checker.vbf_triggers = lambda e: (
                ["HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1"]
                    if (e.run < 317509 and not self.isMC)
                    else ["HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"])
                
        pass

    #def beginJob(self):
    #    pass

    #def endJob(self):
    #    pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('pairType', 'I')
        self.out.branch('dau1_index', 'I')
        self.out.branch('dau2_index', 'I')
        self.out.branch('isVBFtrigger', 'I')
        self.out.branch('isOS', 'I')

        self.out.branch('dau1_eta', 'F')
        self.out.branch('dau1_phi', 'F')
        self.out.branch('dau1_dxy', 'F')
        self.out.branch('dau1_dz', 'F')
        self.out.branch('dau1_q', 'I')
        self.out.branch('dau1_iso', 'F')
        self.out.branch('dau1_decayMode', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSe', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSmu', 'I')
        self.out.branch('dau1_idDeepTau2017v2p1VSjet', 'I')

        self.out.branch('dau2_eta', 'F')
        self.out.branch('dau2_phi', 'F')
        self.out.branch('dau2_dxy', 'F')
        self.out.branch('dau2_dz', 'F')
        self.out.branch('dau2_q', 'I')
        self.out.branch('dau2_iso', 'F')
        self.out.branch('dau2_decayMode', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSe', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSmu', 'I')
        self.out.branch('dau2_idDeepTau2017v2p1VSjet', 'I')
        
        self.histo = ROOT.TH1D("InsideHHLepton", "", 21, -1, 20)
        bins = [
            "all", "goodmuon events", "muon-tau pairs", "deltaR > 0.5",
            "pass muon trigger", "pass lepton veto",
            "goodele events", "ele-tau pairs", "deltaR > 0.5",
            "pass ele trigger", "pass lepton veto",
            "goodtau events", "tau-tau pairs",
            "pass tau trigger", "pass lepton veto",
        ]
        for ibin, binname in enumerate(bins):
            self.histo.GetXaxis().SetBinLabel(ibin + 1, binname)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        if "histos" not in [key.GetName() for key in outputFile.GetListOfKeys()]:
            outputFile.mkdir("histos")
        outputFile.cd("histos")
        self.histo.Write()
        prevdir.cd()
    
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        self.histo.Fill(-1)
        # muon-tau channels
        goodmuons = []
        for imuon, muon in enumerate(muons):
            if (abs(muon.eta) > 2.1 or muon.pfRelIso04_all > 0.15 or abs(muon.dxy) > 0.045
                    or abs(muon.dz) > 0.2 or not muon.tightId):
                continue
            goodmuons.append((imuon, muon))
        if goodmuons:
            self.histo.Fill(0)
            goodtaus = []
            for itau, tau in enumerate(taus):
                if (tau.idDeepTau2017v2p1VSmu < 15 or tau.idDeepTau2017v2p1VSe < 7
                        or tau.idDeepTau2017v2p1VSjet < 1):
                    continue
                if abs(tau.dz) > 0.2:
                    continue
                if tau.decayMode not in [0, 1, 10, 11]:
                    continue
                # common tau pt req for both single and cross triggers
                if eval("tau.pt%s" % self.tau_syst) <= 20.:
                    continue
                goodtaus.append((itau, tau))
            muontaupairs = []
            for (imuon, muon) in goodmuons:
                for (itau, tau) in goodtaus:
                    self.histo.Fill(1)
                    if tau.DeltaR(muon) < 0.5: continue
                    self.histo.Fill(2)
                    if not self.trigger_checker.check_mutau(event,
                            eval("muon.pt%s" % self.muon_syst), muon.eta, muon.phi,
                            eval("tau.pt%s" % self.tau_syst), tau.eta, tau.phi, th1=1, th2=5):
                        continue
                    self.histo.Fill(3)
                    muontaupair = LeptonTauPair(
                        muon, eval("muon.pt%s" % self.muon_syst), muon.pfRelIso04_all,
                        tau, eval("tau.pt%s" % self.tau_syst), tau.rawDeepTau2017v2p1VSjet)
                    # if muontaupair.check_charge():
                    muontaupairs.append((imuon, itau, muontaupair))

            if len(muontaupairs) != 0:
                muontaupairs.sort(key=lambda x: x[2], reverse=True)
                muon, tau = muontaupairs[0][2].pair

                fail_lepton_veto, _ = lepton_veto(electrons, muons, taus, muon)
                if fail_lepton_veto:
                    return False
                self.histo.Fill(4)

                self.out.fillBranch("pairType", 0)
                self.out.fillBranch("isVBFtrigger", 0)
                self.out.fillBranch("isOS", int(muontaupairs[0][2].check_charge()))

                self.out.fillBranch("dau1_index", muontaupairs[0][0])
                self.out.fillBranch("dau1_eta", muon.eta)
                self.out.fillBranch("dau1_phi", muon.phi)
                self.out.fillBranch("dau1_dxy", muon.dxy)
                self.out.fillBranch("dau1_dz", muon.dz)
                self.out.fillBranch("dau1_q", muon.charge)
                self.out.fillBranch("dau1_iso", muon.pfRelIso04_all)
                self.out.fillBranch("dau1_decayMode", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", -1)

                self.out.fillBranch("dau2_index", muontaupairs[0][1])
                self.out.fillBranch("dau2_eta", tau.eta)
                self.out.fillBranch("dau2_phi", tau.phi)
                self.out.fillBranch("dau2_dxy", tau.dxy)
                self.out.fillBranch("dau2_dz", tau.dz)
                self.out.fillBranch("dau2_q", tau.charge)
                self.out.fillBranch("dau2_iso", tau.rawIso)
                self.out.fillBranch("dau2_decayMode", tau.decayMode)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau.idDeepTau2017v2p1VSe)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau.idDeepTau2017v2p1VSmu)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau.idDeepTau2017v2p1VSjet)
                return True

        # electron-tau channels
        goodelectrons = []
        for ielectron, electron in enumerate(electrons):
            # if (not (electron.mvaFall17V2Iso_WP80 or electron.mvaFall17V2noIso_WP80)
            if ((not electron.mvaFall17V2Iso_WP80)
                    or abs(electron.dxy) > 0.045 or abs(electron.dz) > 0.2):
                continue
            goodelectrons.append((ielectron, electron))
        if goodelectrons:
            self.histo.Fill(5)
            goodtaus = []
            for itau, tau in enumerate(taus):
                if (tau.idDeepTau2017v2p1VSmu < 15 or tau.idDeepTau2017v2p1VSe < 7
                        or tau.idDeepTau2017v2p1VSjet < 1):
                    continue
                if abs(tau.dz) > 0.2:
                    continue
                if tau.decayMode not in [0, 1, 10, 11]:
                    continue
                # common tau pt req for both single and cross triggers
                if eval("tau.pt%s" % self.tau_syst) <= 20.:
                    continue
                goodtaus.append((itau, tau))

            electrontaupairs = []
            for (ielectron, electron) in goodelectrons:
                for (itau, tau) in goodtaus:
                    self.histo.Fill(6)
                    if tau.DeltaR(electron) < 0.5: continue
                    self.histo.Fill(7)
                    if not self.trigger_checker.check_etau(event,
                            eval("electron.pt%s" % self.electron_syst), electron.eta, electron.phi,
                            eval("tau.pt%s" % self.tau_syst), tau.eta, tau.phi, th1=1, th2=5):
                        continue
                    self.histo.Fill(8)
                    electrontaupair = LeptonTauPair(
                        electron, eval("electron.pt%s" % self.electron_syst), electron.pfRelIso03_all,
                        tau, eval("tau.pt%s" % self.tau_syst), tau.rawDeepTau2017v2p1VSjet)
                    # if electrontaupair.check_charge():
                    electrontaupairs.append((ielectron, itau, electrontaupair))

            if len(electrontaupairs) != 0:
                electrontaupairs.sort(key=lambda x: x[2], reverse=True)
                electron, tau = electrontaupairs[0][2].pair

                fail_lepton_veto, _ = lepton_veto(electrons, muons, taus, electron)
                if fail_lepton_veto:
                    return False
                self.histo.Fill(9)
                self.out.fillBranch("pairType", 1)
                self.out.fillBranch("isVBFtrigger", 0)
                self.out.fillBranch("isOS", int(electrontaupairs[0][2].check_charge()))

                self.out.fillBranch("dau1_index", electrontaupairs[0][0])
                self.out.fillBranch("dau1_eta", electron.eta)
                self.out.fillBranch("dau1_phi", electron.phi)
                self.out.fillBranch("dau1_dxy", electron.dxy)
                self.out.fillBranch("dau1_dz", electron.dz)
                self.out.fillBranch("dau1_q", electron.charge)
                self.out.fillBranch("dau1_iso", electron.pfRelIso03_all)
                self.out.fillBranch("dau1_decayMode", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", -1)
                self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", -1)

                self.out.fillBranch("dau2_index", electrontaupairs[0][1])
                self.out.fillBranch("dau2_eta", tau.eta)
                self.out.fillBranch("dau2_phi", tau.phi)
                self.out.fillBranch("dau2_dxy", tau.dxy)
                self.out.fillBranch("dau2_dz", tau.dz)
                self.out.fillBranch("dau2_q", tau.charge)
                self.out.fillBranch("dau2_iso", tau.rawIso)
                self.out.fillBranch("dau2_decayMode", tau.decayMode)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau.idDeepTau2017v2p1VSe)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau.idDeepTau2017v2p1VSmu)
                self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau.idDeepTau2017v2p1VSjet)

                return True

        goodtaus = []
        for itau, tau in enumerate(taus):
            if (tau.idDeepTau2017v2p1VSmu < 1 or tau.idDeepTau2017v2p1VSe < 3
                    or tau.idDeepTau2017v2p1VSjet < 1):
                continue
            if abs(tau.dz) > 0.2:
                continue
            if tau.decayMode not in [0, 1, 10, 11]:
                continue
            goodtaus.append((itau, tau))

        if goodtaus:
            self.histo.Fill(10)
        tautaupairs = []
        for i in range(len(goodtaus)):
            for j in range(len(goodtaus)):
                if i == j:
                    continue
                self.histo.Fill(11)
                tau1_index = goodtaus[i][0]
                tau1 = goodtaus[i][1]
                tau2_index = goodtaus[j][0]
                tau2 = goodtaus[j][1]

                if tau1.DeltaR(tau2) < 0.5: continue

                pass_ditau = self.trigger_checker.check_tautau(event,
                    eval("tau1.pt%s" % self.tau_syst), tau1.eta, tau1.phi,
                    eval("tau2.pt%s" % self.tau_syst), tau2.eta, tau2.phi, abs_th1=40, abs_th2=40)
                # passing vbf trigger ONLY
                pass_vbf = (not pass_ditau) and self.trigger_checker.check_vbftautau(event,
                    eval("tau1.pt%s" % self.tau_syst), tau1.eta, tau1.phi,
                    eval("tau2.pt%s" % self.tau_syst), tau2.eta, tau2.phi, abs_th1=25, abs_th2=25)

                if not (pass_ditau or pass_vbf):
                    continue
                self.histo.Fill(12)
                pass_vbf = int(pass_vbf)
                tautaupair = LeptonTauPair(
                    tau1, eval("tau1.pt%s" % self.tau_syst), tau1.rawDeepTau2017v2p1VSjet,
                    tau2, eval("tau2.pt%s" % self.tau_syst), tau2.rawDeepTau2017v2p1VSjet)
                # if tautaupair.check_charge():
                tautaupairs.append((tau1_index, tau2_index, tautaupair, pass_vbf))

        if len(tautaupairs) != 0:
            tautaupairs.sort(key=lambda x: x[2], reverse=True)
            tau1, tau2 = tautaupairs[0][2].pair

            fail_lepton_veto, _ = lepton_veto(electrons, muons, taus)
            if fail_lepton_veto:
                return False
            self.histo.Fill(13)

            self.out.fillBranch("pairType", 2)
            self.out.fillBranch("isVBFtrigger", tautaupairs[0][3])
            self.out.fillBranch("isOS", int(tautaupairs[0][2].check_charge()))

            self.out.fillBranch("dau1_index", tautaupairs[0][0])
            self.out.fillBranch("dau1_eta", tau1.eta)
            self.out.fillBranch("dau1_phi", tau1.phi)
            self.out.fillBranch("dau1_dxy", tau1.dxy)
            self.out.fillBranch("dau1_dz", tau1.dz)
            self.out.fillBranch("dau1_q", tau1.charge)
            self.out.fillBranch("dau1_iso", tau1.rawIso)
            self.out.fillBranch("dau1_decayMode", tau1.decayMode)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSe", tau1.idDeepTau2017v2p1VSe)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSmu", tau1.idDeepTau2017v2p1VSmu)
            self.out.fillBranch("dau1_idDeepTau2017v2p1VSjet", tau1.idDeepTau2017v2p1VSjet)

            self.out.fillBranch("dau2_index", tautaupairs[0][1])
            self.out.fillBranch("dau2_eta", tau2.eta)
            self.out.fillBranch("dau2_phi", tau2.phi)
            self.out.fillBranch("dau2_dxy", tau2.dxy)
            self.out.fillBranch("dau2_dz", tau2.dz)
            self.out.fillBranch("dau2_q", tau2.charge)
            self.out.fillBranch("dau2_iso", tau2.rawIso)
            self.out.fillBranch("dau2_decayMode", tau2.decayMode)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSe", tau2.idDeepTau2017v2p1VSe)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSmu", tau2.idDeepTau2017v2p1VSmu)
            self.out.fillBranch("dau2_idDeepTau2017v2p1VSjet", tau2.idDeepTau2017v2p1VSjet)

            return True
        return False


class HHLeptonVariableProducer(JetLepMetModule):
    def __init__(self, isMC, *args, **kwargs):
        super(HHLeptonVariableProducer, self).__init__(*args, **kwargs)
        self.isMC = isMC

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch('dau1_pt%s' % self.lep_syst, 'F')
        self.out.branch('dau1_mass%s' % self.lep_syst, 'F')
        self.out.branch('dau2_pt%s' % self.lep_syst, 'F')
        self.out.branch('dau2_mass%s' % self.lep_syst, 'F')
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
        taus = Collection(event, "Tau")

        _, _, dau1_tlv, dau2_tlv = self.get_daus(event, muons, electrons, taus)
        self.out.fillBranch("dau1_pt%s" % self.lep_syst, dau1_tlv.Pt())
        self.out.fillBranch("dau1_mass%s" % self.lep_syst, dau1_tlv.M())
        self.out.fillBranch("dau2_pt%s" % self.lep_syst, dau2_tlv.Pt())
        self.out.fillBranch("dau2_mass%s" % self.lep_syst, dau2_tlv.M())
        return True


def HHLepton(**kwargs):
    return lambda: HHLeptonProducer(**kwargs)


def HHLeptonVariable(**kwargs):
    return lambda: HHLeptonVariableProducer(**kwargs)


class HHLeptonRDFProducer(JetLepMetSyst):
    def __init__(self, isMC, year, runEra, pairType_filter, *args, **kwargs):
        super(HHLeptonRDFProducer, self).__init__(isMC=isMC, year=year, *args, **kwargs)
        self.isMC = isMC
        self.doGenCutFlow = kwargs.pop("doGenCutFlow", False)
        self.year = year
        self.runEra = runEra
        self.isRun3 = kwargs.pop("isRun3", False)
        self.pairType_filter = pairType_filter
        self.deeptau_version = kwargs.pop("deeptau_version", "2017v2p1")
        self.deepboostedtau_version = kwargs.pop("deepboostedtau_version", "2018v2p7")
        self.isV10 = kwargs.pop("isV10", False)
        self.useBoostedTaus = kwargs.pop("useBoostedTaus")
        self.tau_priority = kwargs.pop("tau_priority")
        # HPS tau deepTau working points
        vvvl_vsjet = kwargs.pop("vvvl_vsjet")
        vl_vse = kwargs.pop("vl_vse")
        vvl_vse = kwargs.pop("vvl_vse")
        t_vsmu = kwargs.pop("t_vsmu")
        vl_vsmu = kwargs.pop("vl_vsmu")
        # deepBoostedTau WPs
        if self.useBoostedTaus:
            # VsMu & VsE are not used
            boostedTau_VsMu_threshold = kwargs.pop("boostedTau_VsMu_threshold")
            boostedTau_VsE_threshold = kwargs.pop("boostedTau_VsE_threshold")
            boostedTau_VsJet_threshold = kwargs.pop("boostedTau_VsJet_threshold") # !!!!raw threshold
        else:
            boostedTau_VsMu_threshold, boostedTau_VsE_threshold, boostedTau_VsJet_threshold = 0, 0, 0

        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        # DOCUMENTATION
        # The list of triggers here is checked against the list of branches in NanoAOD, non-existent ones are replaced by false.
        # A "Vbool triggers" array is built with the HLT_* trigger paths and passed to get_*tau_triggers functions
        # Then a trig_req object is built. pass = HLT_* branch (or false in case it does not exist). 
        # struct trig_req {bool pass; float pt1; float eta1; float pt2; float eta2; std::vector<std::vector<int>> bits; };
        # the values here are used to make cuts on *offline* objects, nb1 is electron/muon/hadTau, nb2 is hadTau

        # boostedTau triggers are not actually used if we use MET trigger

        self.mutau_triggers = ["HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
            "HLT_IsoTkMu24", "HLT_IsoTkMu22_eta2p1", "HLT_IsoMu24", "HLT_IsoMu27",
            "HLT_IsoMu19_eta2p1_LooseIsoPFTau20", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"]
        #self.mu_boostedTau_triggers = ["HLT_Mu50", "HLT_OldMu100", "HLT_TkMu100"]
        self.etau_triggers = ["HLT_Ele25_eta2p1_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
            "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"]
        #self.e_boostedTau_triggers = ["HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"] # TODO see if we should include HLT_Photon200
        self.tautau_triggers = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1"]
        #self.boostedTau_boostedTau_triggers = ["HLT_AK8PFHT800_TrimMass50"]
        self.tautaujet_triggers = ["HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60"]
        self.vbf_triggers = ["HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg",
            "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1",
            "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"]
        
        # 2018 only
        assert self.year == 2018
        self.boostedTau_MET_triggers = "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight"

        if not os.getenv("_HHLepton"):
            os.environ["_HHLepton"] = "_HHLepton"

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/HHLeptonInterface.h".format(base))
            ROOT.gInterpreter.Declare(f"""
                auto HHLepton = HHLeptonInterface(
                    {vvvl_vsjet}, {vl_vse}, {vvl_vse}, {t_vsmu}, {vl_vsmu},
                    {boostedTau_VsMu_threshold}, {boostedTau_VsE_threshold}, {boostedTau_VsJet_threshold});
            """)

            if self.year == 2016:
                if not self.isV10:
                    raise ValueError("TauTriggerSFs not implemented yet for 2016 v9")
                else:
                    ROOT.gInterpreter.Declare("""
                        using Vbool = const ROOT::RVec<Bool_t>&;
                        std::vector<trig_req> get_mutau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
                            // 2016 is special for trigger : see https://github.com/cms-sw/cmssw/blob/0ec1f22895570d81284e01d8189d86d5622e52ca/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L267
                            // Trig filterBits : 2 is Iso, 8 is IsoTkMu
                            trigger_reqs.push_back(trig_req({triggers[2], 26, 2.4, 0, 0, 24, 0, {{2}, {}}})); // HLT_IsoMu24
                            trigger_reqs.push_back(trig_req({triggers[8], 26, 2.4, 0, 0, 24, 0, {{8}, {}}})); // HLT_IsoTkMu24
                            // https://twiki.cern.ch/twiki/bin/view/CMS/TauTrigger
                            // For filter bits, TWiki seems wrong : says filterBits&32 but this bit does not exist 
                            // We use 2 (Iso hltL3cr*IsoFiltered0p09) & 4 (OverlapFilter PFTau *OverlapFilter*IsoMu*PFTau*)
                            // trigger_reqs.push_back(trig_req({triggers[6], 20, 2.1, 25, 2.1, {{2, 4}, {1, 32}}})); // HLT_IsoMu19_eta2p1_LooseIsoPFTau20 : not sure what it is, but not in list of recommended triggers by tau pog
                            trigger_reqs.push_back(trig_req({triggers[7], 21, 2.1, 25, 2.1, 19, 0, {{2, 4}, {1, 32}}})); // HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_etau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // Bit 1 (->2^1=2) is 1e (WPTight).
                            // eta is 2.1 (pre-pixel upgrade for 2016)
                            trigger_reqs.push_back(trig_req({triggers[0], 26, 2.1, 0, 0, 25, 0, {{2}, {}}})); // HLT_Ele25_eta2p1_WPTight_Gsf (from EG twiki : unprescaled but L1 turn on limited, still can be used)
                            // no cross (single ele has similar pt threshold)
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautau_triggers(
                                Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            if (!isMC) {
                                if ((runEra >= 2) && (runEra <= 7)) { // B to G inclusive (tau twiki says F but is wrong)
                                    trigger_reqs.push_back(trig_req({triggers[0], 40, 2.1, 40, 2.1, 0, 0, {{2, 256}, {2, 256}}})); // HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg
                                }
                                else if (runEra == 8) { // H
                                    trigger_reqs.push_back(trig_req({triggers[1], 40, 2.1, 40, 2.1, 0, 0, {{2, 256}, {2, 256}}})); // HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg
                                }               
                            }
                            else {
                                trigger_reqs.push_back(trig_req({triggers[0], 40, 2.1, 40, 2.1, 0, 0, {{2, 256}, {2, 256}}})); // HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg
                                // trigger_reqs.push_back(trig_req({triggers[1], 40, 2.1, 40, 2.1, 0, 0, {{2, 256}, {2, 256}}})); // HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg (data only)
                            }
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautaujet_triggers(Vbool triggers, bool isRun3) {
                            std::vector<trig_req> trigger_reqs;
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_vbf_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            return trigger_reqs;
                        }
                    """)

            if self.year == 2017:
                if not self.isV10:
                    raise ValueError("TauTriggerSFs not implemented yet for 2017 v9")
                else:
                    ROOT.gInterpreter.Declare("""
                        using Vbool = const ROOT::RVec<Bool_t>&;
                        std::vector<trig_req> get_mutau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // trigger_reqs.push_back(trig_req({triggers[4], 26, 2.4, 20, 2.3, 24, 0, {{2, 8}, {}}})); // HLT_IsoMu24 prescaled in 2017 -> do not use (we don't have SFs)
                            trigger_reqs.push_back(trig_req({triggers[5], 29, 2.4, 0, 0, 27, 0, {{2, 8}, {}}})); // HLT_IsoMu27
                            // No doc on this trigger in Tau twiki
                            // Filter bits muon leg : 4 : Iso (not sure if needed but does not hurt), 64 : 1mu-1tau
                            // Tau leg : bit 0 (-> 1) : LooseChargedIso, bit 9 (->512) : mu-tau
                            trigger_reqs.push_back(trig_req({triggers[8], 22, 2.1, 32, 2.1, 20, 27, {{4, 64}, {1, 512}}})); // HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_etau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // 2017 is special for EG (new pixels) -> big mess (https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary and https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers )
                            // Our SFs are for Ele32 & use the special emulation of HLT_Ele32_WPTight_Gsf for data runs where this does not exist (see twiki)
                            // bits : hltEle32L1DoubleEGWPTightGsfTrackIsoFilter -> bit 1 in Nano -> filter 2 =(2^1)
                            // hltEGL1SingleEGOrFilter -> bit 10  in Nano -> filter 1024
                            trigger_reqs.push_back(trig_req({triggers[1], 33, 2.5, 20, 2.3, 32, 0, {{2, 1024}, {}}})); // HLT_Ele32_WPTight_Gsf_L1DoubleEG with bits for HLT_Ele32_WPTight_Gsf emulation
                            // trigger_reqs.push_back(trig_req({triggers[2], 33, 2.3, 20, 2.3, 32, 0, {{2}, {}}})); // HLT_Ele32_WPTight_Gsf 
                            // trigger_reqs.push_back(trig_req({triggers[3], 36, 2.3, 20, 2.3, 35, 0, {{2}, {}}})); // HLT_Ele35_WPTight_Gsf
                            // Filter bits are difference nanoV9 vs v12
                            // Tau leg : bit 0 (->1) : LooseChargedIso, bit 8 (->256) ditau
                            trigger_reqs.push_back(trig_req({triggers[4], 26, 2.1, 35, 2.1, 0, 0, {{64}, {1, 256}}})); // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautau_triggers(
                                Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // Tau POG doc for ditau : 
                            // TrigObj_id == 15 && (TrigObj_filterBits&64)!=0 && (
                            //        ( (TrigObj_filterBits&4)!=0 && (TrigObj_filterBits&16)!=0 ) ||
                            //        ( TrigObj_pt > 40 && ( (TrigObj_filterBits&2)!=0 && (TrigObj_filterBits&16)!=0 ) 
                            //    || (TrigObj_filterBits&4)!=0 ) ) 
                            // Summary : (4&16&64) || (2&16&64&Trig_pt>40) || (4&64)
                            // Filter bits :; 1(->2):MediumChargedIso, 2(->4):TightChargedIso, 4(->16):TightID OOSC photons (*TightOOSCPhotons*), 6(->64):charged iso di-tau (hlt*DoublePFTau*TrackPt1*ChargedIsolation*Dz02*)
                            trigger_reqs.push_back(trig_req({triggers[2], 40, 2.1, 40, 2.1, 0, 0, {{4, 16, 64}, {4, 16, 64}}})); // HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg
                            // the offline req is indeed 40 GeV according to tau pog twiki, + need to put online 40 GeV req.
                            trigger_reqs.push_back(trig_req({triggers[3], 40, 2.1, 40, 2.1, 40, 40, {{2, 16, 64}, {2, 16, 64}}})); // HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg
                            trigger_reqs.push_back(trig_req({triggers[4], 40, 2.1, 40, 2.1, 0, 0, {{4, 64}, {4, 64}}})); // HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautaujet_triggers(Vbool triggers, bool isRun3) {
                            std::vector<trig_req> trigger_reqs;
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_vbf_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            return trigger_reqs;
                        }
                    """)

            if self.year == 2018 or self.year == 2022:
                if not self.isV10:
                    assert False
                    ROOT.gInterpreter.Declare("""
                        using Vbool = const ROOT::RVec<Bool_t>&;
                        std::vector<trig_req> get_mutau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            trigger_reqs.push_back(trig_req({triggers[4], 25, 2.3, 20, 2.3, {{2, 8}, {}}}));
                            trigger_reqs.push_back(trig_req({triggers[5], 28, 2.3, 20, 2.3, {{2, 8}, {}}}));
                            if (!isMC && run < 317509)
                                trigger_reqs.push_back(trig_req({triggers[8], 21, 2.3, 32, 2.1, {{64}, {1, 256}}}));
                            else
                                trigger_reqs.push_back(trig_req({triggers[9], 21, 2.3, 32, 2.1, {{4}, {1, 16}}}));
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_etau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            trigger_reqs.push_back(trig_req({triggers[2], 33, 2.3, 20, 2.3, {{2}, {}}}));
                            trigger_reqs.push_back(trig_req({triggers[3], 36, 2.3, 20, 2.3, {{2}, {}}}));
                            if (!isMC && run < 317509)
                                trigger_reqs.push_back(trig_req({triggers[4], 25, 2.3, 35, 2.1, {{64}, {1, 128}}}));
                            else
                                trigger_reqs.push_back(trig_req({triggers[5], 25, 2.3, 35, 2.1, {{8}, {1, 16}}}));
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautau_triggers(
                                Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            trigger_reqs.push_back(trig_req({triggers[2], 40, 2.1, 40, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                            trigger_reqs.push_back(trig_req({triggers[3], 40, 2.1, 40, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                            if (!isMC && run < 317509)
                                trigger_reqs.push_back(trig_req({triggers[4], 40, 2.1, 40, 2.1, {{64, 4}, {64, 4}}}));
                            else
                                trigger_reqs.push_back(trig_req({triggers[5], 40, 2.1, 40, 2.1, {{2, 16}, {2, 16}}}));
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_tautaujet_triggers(Vbool triggers, bool isRun3) { // DUMMY FUNCTION
                            std::vector<trig_req> trigger_reqs;
                            if (isRun3)
                                trigger_reqs.push_back(trig_req({triggers[0], 35, 2.1, 35, 2.1, {{8, 32, 128, 16384}, {8, 32, 128, 16384}}}));
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_vbf_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            if (!isMC && run < 317509)
                                trigger_reqs.push_back(trig_req({triggers[1], 25, 2.1, 25, 2.1, {{64, 4, 8}, {64, 4, 8}}}));
                            else
                                trigger_reqs.push_back(trig_req({triggers[2], 25, 2.1, 25, 2.1, {{512, 1, 16}, {512, 1, 16}}}));
                            return trigger_reqs;
                        }
                    """)
                else:
                    ROOT.gInterpreter.Declare("""
                        using Vbool = const ROOT::RVec<Bool_t>&;
                        std::vector<trig_req> get_mutau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            trigger_reqs.push_back(trig_req({triggers[4], 26, 2.4, 20, 2.3, 24, 0, {{2, 8}, {}}})); // HLT_IsoMu24 (not prescaled for 2018)
                            // trigger_reqs.push_back(trig_req({triggers[5], 29, 2.4, 20, 2.3, 27, 0, {{2, 8}, {}}})); // HLT_IsoMu27 -> not to be used for 2018
                            // https://twiki.cern.ch/twiki/bin/view/CMS/TauTrigger
                            // The TWiki suggests using HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 rather than HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 but our SFs are presumably for the latter (trigger used for gamma gamma ->tautau)
                            // Muon filter bits : 2(->4): *OverlapFilterIsoMu*PFTau*, 6(->64): hlt*OverlapFilterIsoMu*PFTau* (the 2 seem identical ?)
                            // Tau filter bits : 0(->1): LooseChargedIso, 5(->32): *Hps*, 9(->512): hlt*OverlapFilterIsoMu*PFTau* (mutau)
                            // Requiring HPS is not in  Tau POG wiki and probably is not needed as all paths matching hlt*OverlapFilterIsoMu*PFTau* are HPS when that is available
                            // Requiring LooseChargedIso is also not in Tau POG twiki
                            if (!isMC && run < 317509) {
                                trigger_reqs.push_back(trig_req({triggers[8], 22, 2.1, 32, 2.1, 0, 27, {{64}, {1, 512}}})); // HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1
                            }
                            else {
                                // used to be {1, 32} for tau leg which is wrong I think ? that would require every tau LooseIso HPS trigger
                                trigger_reqs.push_back(trig_req({triggers[9], 22, 2.1, 32, 2.1, 0, 27, {{64}, {1, 512}}})); // HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1
                            }
                            return trigger_reqs;
                        }
                        /*std::vector<trig_req> get_mu_boostedTau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs; 
                            // TODO checked for 2018 but 2016 is different !!!!!!!!!!!! TODO check 2016
                            // https://indico.cern.ch/event/1080036/contributions/4542924/attachments/2318322/3947065/210928_ULTriggerSF_kHwang.pdf
                            // TrigObj filter bits : 
                            // 10: ["hltL3fL1sMu*L3Filtered50*","hltL3fL1sMu*TkFiltered50*"],"1mu (Mu50)"
                            // 11: ["hltL3fL1sMu*L3Filtered100*","hltL3fL1sMu*TkFiltered100*"],"1mu (Mu100)"
                            trigger_reqs.push_back(trig_req({triggers[0], 52, 2.4, 20, 2.3, 50, 0, {{10}, {}}})); // HLT_Mu50 ie hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q
                            trigger_reqs.push_back(trig_req({triggers[1], 52, 2.4, 20, 2.3, 100, 0, {{11}, {}}})); // HLT_OldMu100 ie hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q
                            trigger_reqs.push_back(trig_req({triggers[2], 52, 2.4, 20, 2.3, 100, 0, {{11}, {}}})); // HLT_TkMu100 ie hltL3fL1sMu25f0TkFiltered100Q
                            return trigger_reqs;
                        }*/
                        std::vector<trig_req> get_etau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            trigger_reqs.push_back(trig_req({triggers[2], 33, 2.5, 20, 2.3, 32, 0, {{2}, {}}})); // HLT_Ele32_WPTight_Gsf
                            // trigger_reqs.push_back(trig_req({triggers[3], 36, 2.1, 20, 2.3, {{2}, {}}})); // HLT_Ele35_WPTight_Gsf
                            // gg->tautau analysis does OR with HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 but extremly few events pass this and not our trigger (<0.01% in ttbar)
                            if (!isMC && run < 317509) {
                                trigger_reqs.push_back(trig_req({triggers[4], 25, 2.1, 35, 2.1, 0, 0, {{64}, {1, 256}}})); // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
                            }
                            else {
                                trigger_reqs.push_back(trig_req({triggers[5], 25, 2.1, 35, 2.1, 0, 0, {{8}, {1, 256}}})); // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1
                            }
                            return trigger_reqs;
                        }
                        /*std::vector<trig_req> get_e_boostedTau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs; // filter bits checked for 2018 but TODO check for 2016 !! + TODO check thresholds for SFs when we have SFs
                            // Filter bits (2018)
                            // 11: "filter('hltEle*CaloIdVTGsfTrkIdTGsfDphiFilter')","1e (CaloIdVT_GsfTrkIdT)"
                            // 12: path('HLT_Ele*PFJet*')","1e (PFJet)"
                            trigger_reqs.push_back(trig_req({triggers[0], 120, 2.5, 20, 2.3, 115, 0, {{11}, {}}})); // HLT_Ele115_CaloIdVT_GsfTrkIdT ie hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter
                            trigger_reqs.push_back(trig_req({triggers[1], 55, 2.5, 20, 2.3, 50, 0, {{11, 12}, {}}})); // HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 ie hltEle50CaloIdVTGsfTrkIdTCentralPFJet165EleCleaned
                            return trigger_reqs;
                        }*/
                        std::vector<trig_req> get_tautau_triggers(
                                Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // Tau POG doc for ditau : (same as 2017)
                            // TrigObj_id == 15 && (TrigObj_filterBits&64)!=0 && (
                            //        ( (TrigObj_filterBits&4)!=0 && (TrigObj_filterBits&16)!=0 ) ||
                            //        ( TrigObj_pt > 40 && ( (TrigObj_filterBits&2)!=0 && (TrigObj_filterBits&16)!=0 ) 
                            //    || (TrigObj_filterBits&4)!=0 ) ) 
                            // Summary : (4&16&64) || (2&16&64&Trig_pt>40) || (4&64)
                            // Filter bits : 1(->2):MediumChargedIso, 2(->4):TightChargedIso, 4(->16):TightID OOSC photons (*TightOOSCPhotons*), 6(->64):charged iso di-tau (hlt*DoublePFTau*TrackPt1*ChargedIsolation*Dz02*)
                            if (!isMC && run < 317509) {
                                trigger_reqs.push_back(trig_req({triggers[2], 40, 2.1, 40, 2.1, 0, 0, {{4, 16, 64}, {4, 16, 64}}})); // HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg
                                // the offline req is indeed 40 GeV according to tau pog twiki, + need to put online 40 GeV req.
                                trigger_reqs.push_back(trig_req({triggers[3], 40, 2.1, 40, 2.1, 40, 40, {{2, 16, 64}, {2, 16, 64}}})); // HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg
                                trigger_reqs.push_back(trig_req({triggers[4], 40, 2.1, 40, 2.1, 0, 0, {{4, 64}, {4, 64}}})); // HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg
                            }
                            else {
                                trigger_reqs.push_back(trig_req({triggers[5], 40, 2.1, 40, 2.1, 0, 0, {{64}, {64}}})); // HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg
                            }
                            return trigger_reqs;
                        }
                        /*std::vector<trig_req> get_boostedTau_boostedTau_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // TrigObj bits for FatJet : filterBits done, TODO trigger matching & thresholds etc
                            // 0: mksel("coll('hltAK8PFJetsCorrected')"),    #1, always present
                            // 1: mksel(["hltAK8SingleCaloJet200"]),         #2, always present
                            // 2: mksel("coll('hltAK8PFSoftDropJets230')"),  #4, present if nothing else below is fired, otherwise 12, 20, 28, 52, 60
                            // 3: mksel(["hltAK8SinglePFJets230SoftDropMass40BTagParticleNetBB0p35",
                            // 4: "hltAK8SinglePFJets250SoftDropMass40BTagParticleNetBB0p35",
                            // 5: "hltAK8SinglePFJets275SoftDropMass40BTagParticleNetBB0p35"]), # 12 if nothing below is fired, #28 if also "hltAK8DoublePFJetSDModMass30", #60 if also "hltAK8DoublePFJetSDModMass50" 
                            // 6: mksel(["hltAK8DoublePFJetSDModMass30"]), # 16 if onthing else (except #1), 20 if also #4, 28 if also #12
                            // 7: mksel(["hltAK8DoublePFJetSDModMass50"]), # 48 if also (obviously) "hltAK8DoublePFJetSDModMass30", 52 if also #4, #60 if all above
                            
                            trigger_reqs.push_back(trig_req({triggers[0], 0, 0, 0, 0, 0, 0, {{0}, {}}})); // HLT_AK8PFHT800_TrimMass50
                            return trigger_reqs;
                        }*/
                        std::vector<trig_req> get_tautaujet_triggers(Vbool triggers, bool isRun3) {
                            std::vector<trig_req> trigger_reqs;
                            return trigger_reqs;
                        }
                        std::vector<trig_req> get_vbf_triggers(
                                Vbool triggers, bool isMC, int run, int runEra) {
                            std::vector<trig_req> trigger_reqs;
                            // if (!isMC && run < 317509)
                            //     trigger_reqs.push_back(trig_req({triggers[1], 25, 2.1, 25, 2.1, {{1}, {1}}})); // HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1
                            // else
                            //     trigger_reqs.push_back(trig_req({triggers[2], 25, 2.1, 25, 2.1, {{1, 32}, {1, 32}}})); // HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1
                            return trigger_reqs;
                        }
                    """)

    def run(self, df):
        variables = ["pairType", "dau1_index", "dau2_index",
            "isTauTauJetTrigger", "isVBFtrigger", "isOS",
            "dau1_eta", "dau1_phi", "dau1_iso", "dau1_decayMode",
            "dau1_idDeepTauVSe", "dau1_idDeepTauVSmu",
            "dau1_idDeepTauVSjet",
            "dau2_eta", "dau2_phi", "dau2_decayMode",
            "dau2_idDeepTauVSe", "dau2_idDeepTauVSmu",
            "dau2_idDeepTauVSjet"
        ]                  

        all_branches = df.GetColumnNames()
        for ib, branch in enumerate(self.mutau_triggers):
            if branch not in all_branches:
                self.mutau_triggers[ib] = "false"
        for ib, branch in enumerate(self.etau_triggers):
            if branch not in all_branches:
                self.etau_triggers[ib] = "false"
        for ib, branch in enumerate(self.tautau_triggers):
            if branch not in all_branches:
                self.tautau_triggers[ib] = "false"
        for ib, branch in enumerate(self.tautaujet_triggers):
            if branch not in all_branches:
                self.tautaujet_triggers[ib] = "false"
        for ib, branch in enumerate(self.vbf_triggers):
            if branch not in all_branches:
                self.vbf_triggers[ib] = "false"

        runEras = ["dum", "A", "B", "C", "D", "E", "F", "G", "H"]
        runEra = None
        for _irun, _runEra in enumerate(runEras):
            if self.runEra == _runEra:
                runEra = _irun
                break
        assert runEra != None

        df = df.Define("mutau_triggers", "get_mutau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.mutau_triggers), ("true" if self.isMC else "false"), runEra))
        # df = df.Define("mu_boostedTau_triggers", "get_mu_boostedTau_triggers({%s}, %s, run, %s)" % (
        #     ", ".join(self.mu_boostedTau_triggers), ("true" if self.isMC else "false"), runEra))
        df = df.Define("etau_triggers", "get_etau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.etau_triggers), ("true" if self.isMC else "false"), runEra))
        # df = df.Define("e_boostedTau_triggers", "get_e_boostedTau_triggers({%s}, %s, run, %s)" % (
        #     ", ".join(self.e_boostedTau_triggers), ("true" if self.isMC else "false"), runEra))
        df = df.Define("tautau_triggers", "get_tautau_triggers({%s}, %s, %s, run, %s)" % (
            ", ".join(self.tautau_triggers), ("true" if self.isMC else "false"),
            ("true" if self.isRun3 else "false"), runEra))
        # df = df.Define("boostedTau_boostedTau_triggers", "get_boostedTau_boostedTau_triggers({%s}, %s, run, %s)" % (
        #     ", ".join(self.boostedTau_boostedTau_triggers), ("true" if self.isMC else "false"), runEra))
        df = df.Define("tautaujet_triggers", "get_tautaujet_triggers({%s}, %s)" % (
            ", ".join(self.tautaujet_triggers), ("true" if self.isRun3 else "false")))
        df = df.Define("vbf_triggers", "get_vbf_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.vbf_triggers), ("true" if self.isMC else "false"), runEra))

        Electron_mvaIso_WP80 = "Electron_mvaIso_WP80"
        if Electron_mvaIso_WP80 not in all_branches:
            Electron_mvaIso_WP80 = "Electron_mvaFall17V2Iso_WP80"
        Electron_mvaIso_WP90 = "Electron_mvaIso_WP90"
        if Electron_mvaIso_WP90 not in all_branches:
            Electron_mvaIso_WP90 = "Electron_mvaFall17V2Iso_WP90"
        Electron_mvaNoIso_WP90 = "Electron_mvaNoIso_WP90"
        if Electron_mvaNoIso_WP90 not in all_branches:
            Electron_mvaNoIso_WP90 = "Electron_mvaFall17V2noIso_WP90"

        branches = []
        
        muonGenPartIdx_branch = "Muon_genPartIdx" if self.doGenCutFlow else "{}"
        electronGenPartIdx_branch = "Electron_genPartIdx" if self.doGenCutFlow else "{}"
        boostedTau_genPartIdx_branch = "boostedTau_genPartIdx, boostedTau_genPartFlav" if self.doGenCutFlow else "{}, {}"
        genPairTypeEtc_branch = "GenPairType, genDau1_genPartIdx, genDau2_genPartIdx" if self.doGenCutFlow else "-1, -1, -1"
        if self.useBoostedTaus:
            # boosted taus
            df = df.Define("hh_lepton_results_boostedTaus", "HHLepton.get_boosted_dau_indexes("
                f"{str(self.doGenCutFlow).lower()}, "
                f"Muon_pt{self.muon_syst}, Muon_eta, Muon_phi, Muon_mass{self.muon_syst}, "
                f"Muon_pfRelIso04_all, Muon_dxy, Muon_dz, "
                f"Muon_looseId, Muon_mediumId, Muon_tightId, Muon_charge, "
                f"{muonGenPartIdx_branch}, "
                f"Electron_pt{self.electron_syst}, Electron_eta, Electron_phi, Electron_mass{self.electron_syst}, "
                f"{Electron_mvaIso_WP80}, {Electron_mvaNoIso_WP90}, {Electron_mvaIso_WP90}, Electron_vidNestedWPBitmap, Electron_pfRelIso03_all, "
                f"Electron_dxy, Electron_dz, Electron_charge, "
                f"{electronGenPartIdx_branch}, "
                # TODO boostedTau systematic variations on pt & mass
                f"boostedTau_pt, boostedTau_eta, boostedTau_phi, boostedTau_mass, "
                f"boostedTau_idDeepTau{self.deepboostedtau_version}VSmu, boostedTau_idDeepTau{self.deepboostedtau_version}VSe, "
                f"boostedTau_idDeepTau{self.deepboostedtau_version}VSjet, boostedTau_rawDeepTau{self.deepboostedtau_version}VSjet, "
                "boostedTau_decayMode, boostedTau_charge, "
                f"{boostedTau_genPartIdx_branch}, "
                "boostedTau_Mcounter, { {boostedTau_LeadingMuon_muonIdx, boostedTau_SubLeadingMuon_muonIdx, boostedTau_SubSubLeadingMuon_muonIdx} }, "
                "{ {boostedTau_LeadingMuonPt, boostedTau_SubLeadingMuonPt, boostedTau_SubSubLeadingMuonPt} }, { {boostedTau_LeadingMuonCorrIso, boostedTau_SubLeadingMuonCorrIso, boostedTau_SubSubLeadingMuonCorrIso} }, "
                "boostedTau_Ecounter, { {boostedTau_LeadingElectron_electronIdx, boostedTau_SubLeadingElectron_electronIdx, boostedTau_SubSubLeadingElectron_electronIdx} }, "
                "{ {boostedTau_LeadingElectronPt, boostedTau_SubLeadingElectronPt, boostedTau_SubSubLeadingElectronPt} }, { {boostedTau_LeadingElectronCorrIso, boostedTau_SubLeadingElectronCorrIso, boostedTau_SubSubLeadingElectronCorrIso} }, "
                "TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, "
                "{}, {}, {}, "
                #"mu_boostedTau_triggers, e_boostedTau_triggers, boostedTau_boostedTau_triggers, "
                f"{genPairTypeEtc_branch}"
            ")"
            )
            df = df.Define("pairType_boostedTaus", "hh_lepton_results_boostedTaus.first.pairType")
            branches.append("pairType_boostedTaus")

            if self.doGenCutFlow:
                import cppyy
                # save all FailReason attributes to a branch
                for dauId in [1, 2]:
                    for failReason in dir(cppyy.gbl.FailReason): # this lists all attributes of the C++ object FailReason
                        if failReason.startswith("_") or failReason == "pass":
                            continue
                        df = df.Define(f"cutflow_boostedTaus_dau{dauId}_{failReason}", f"hh_lepton_results_boostedTaus.second.dau{dauId}_fail.{failReason}")
                        branches.append(f"cutflow_boostedTaus_dau{dauId}_{failReason}")
                df = df.Define("cutflow_boostedTaus_leptonVetoFail", "hh_lepton_results_boostedTaus.second.leptonVetoFail")
                branches.append("cutflow_boostedTaus_leptonVetoFail")
                df = df.Define("cutflow_boostedTaus_deltaR", "hh_lepton_results_boostedTaus.second.deltaR")
                branches.append("cutflow_boostedTaus_deltaR")

        # HPS taus
        df = df.Define("hh_lepton_results_HPStaus", "HHLepton.get_dau_indexes("
            f"Muon_pt{self.muon_syst}, Muon_eta, Muon_phi, Muon_mass{self.muon_syst}, "
            f"Muon_pfRelIso04_all, Muon_dxy, Muon_dz, Muon_mediumId, Muon_tightId, Muon_charge, "
            f"Electron_pt{self.electron_syst}, Electron_eta, Electron_phi, Electron_mass{self.electron_syst}, "
            f"{Electron_mvaIso_WP80}, {Electron_mvaNoIso_WP90}, {Electron_mvaIso_WP90}, Electron_pfRelIso03_all, "
            f"Electron_dxy, Electron_dz, Electron_charge, "
            f"Tau_pt{self.tau_syst}, Tau_eta, Tau_phi, Tau_mass{self.tau_syst}, "
            f"Tau_idDeepTau{self.deeptau_version}VSmu, Tau_idDeepTau{self.deeptau_version}VSe, "
            f"Tau_idDeepTau{self.deeptau_version}VSjet, Tau_rawDeepTau{self.deeptau_version}VSjet, "
            "Tau_dz, Tau_decayMode, Tau_charge, "
            "TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, "
            "mutau_triggers, etau_triggers, tautau_triggers, tautaujet_triggers, vbf_triggers"
        ")"
        )
        df = df.Define("pairType_HPSTaus", "hh_lepton_results_HPStaus.pairType")
        branches.append("pairType_HPSTaus")

        if self.useBoostedTaus:
            # boostedTau category need MET trigger fired && offline MET cut to avoid MET turn on (offline cut from Wisconsin analysis)
            boostedTau_trigger_req = f"({self.boostedTau_MET_triggers}) && MET_pt > 180"
            if self.tau_priority == "HPS":
                df = df.Define("isBoostedTau", f"hh_lepton_results_HPStaus.pairType < 0 && {boostedTau_trigger_req}")
            elif self.tau_priority == "boosted":
                df = df.Define("isBoostedTau", f"hh_lepton_results_boostedTaus.first.pairType >= 0 && ({boostedTau_trigger_req})")
            else:
                raise ValueError("tau_priority should be 'HPS' or 'boosted'")
            
            df = df.Define("hh_lepton_results", "isBoostedTau ? hh_lepton_results_boostedTaus.first : hh_lepton_results_HPStaus")
        else:
            df = df.Alias("hh_lepton_results", "hh_lepton_results_HPStaus")
            df = df.Define("isBoostedTau", "False")
        branches.append("isBoostedTau")

        for var in variables:
            branchName = var
            if "DeepTau" in branchName:
                branchName = var[:var.index("VS")] + self.deeptau_version + var[var.index("VS"):]
            df = df.Define(branchName, "hh_lepton_results.%s" % var)
            branches.append(branchName)
        
        # add raw DeepBoostedTau score
        df = df.Define("dau1_rawIdDeepTauVSjet", """
            if (pairType == 2) { """
                  f"return isBoostedTau ? boostedTau_rawDeepTau{self.deepboostedtau_version}VSjet[dau1_index] : Tau_rawDeepTau{self.deeptau_version}VSjet[dau1_index];"
            """
            } else {
                return -1.f; 
            }"""      
        )
        branches.append("dau1_rawIdDeepTauVSjet")
        df = df.Define("dau2_rawIdDeepTauVSjet", f"if (pairType >=0) return isBoostedTau ? boostedTau_rawDeepTau{self.deepboostedtau_version}VSjet[dau2_index] : Tau_rawDeepTau{self.deeptau_version}VSjet[dau2_index]; else return -1.f;")
        branches.append("dau2_rawIdDeepTauVSjet")

        if self.pairType_filter:
            df = df.Filter("pairType >= 0", "HHLeptonRDF")

        return df, branches
        # return df, []


def HHLeptonRDF(**kwargs):
    """
    Returns the index of the two selected taus + several of their variables not affected by
    systematics.

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    :param runEra: run period in caps (data only)
    :type runEra: str

    :param isV10: whether the input sample is from nanoaodV10 (default: ``False``)
    :type isV10: bool

    :param deeptau_version: version of the DeepTau discriminator (default: ``2017v2p1``)
    :type deeptau_version: str

    :param vvvl_vsjet: VVVLoose DeepTauVSjet WP value
    :type vvvl_vsjet: int

    :param vl_vse: VLoose DeepTauVSe WP value
    :type vl_vse: int

    :param vvl_vse: VVLoose DeepTaVSe WP value
    :type vvl_vse: int

    :param t_vsmu: Tight DeepTauVSmu WP value
    :type t_vsmu: int

    :param vl_vsmu: VLoose DeepTauVSmu WP value
    :type vl_vsmu: int

    DeepBoostedTau thresholds : 
    boostedTau_VsMu_threshold, boostedTau_VsE_threshold -> not used
    boostedTau_VsJet_threshold -> *raw* threshold for VsJet score

    :param tau_priority: Give priority to "HPS" taus or to "boosted"Taus

    :param filter: whether to filter out output events if they don't have 2 lepton candidates
    :type filter: bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHLeptonRDF
            path: Tools.Tools.HHLepton
            parameters:
                isMC: self.dataset.process.isMC
                isV10: self.dataset.has_tag("nanoV10")
                year: self.config.year
                runEra: self.dataset.runEra
                runEra: self.dataset.runEra
                vvvl_vsjet: self.config.deeptau.vsjet.VVVLoose
                vl_vse: self.config.deeptau.vse.VLoose
                vvl_vse: self.config.deeptau.vse.VVLoose
                t_vsmu: self.config.deeptau.vsmu.Tight
                vl_vsmu: self.config.deeptau.vsmu.VLoose
                pairType_filter: True

    """
    pairType_filter = kwargs.pop("pairType_filter")
    return lambda: HHLeptonRDFProducer(pairType_filter=pairType_filter, **kwargs)


class HHLeptonVarRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(HHLeptonVarRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("DAU_VAR"):
            os.environ["DAU_VAR"] = "true"
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using VInt = const ROOT::RVec<int>&;
                std::vector<float> get_lepton_values (
                    int pairType, bool isBoostedTau, int dau1_index, int dau2_index,
                    Vfloat muon_pt, Vfloat muon_mass,
                    Vfloat electron_pt, Vfloat electron_mass,
                    Vfloat tau_pt, Vfloat tau_mass,
                    Vfloat boostedtau_pt, Vfloat boostedtau_mass
                )
                {
                    float dau1_pt, dau1_mass, dau2_pt, dau2_mass;
                    if (pairType == 0) {
                        dau1_pt = muon_pt.at(dau1_index);
                        dau1_mass = muon_mass.at(dau1_index);
                    } else if (pairType == 1) {
                        dau1_pt = electron_pt.at(dau1_index);
                        dau1_mass = electron_mass.at(dau1_index);
                    } else if (pairType == 2 && !isBoostedTau) {
                        dau1_pt = tau_pt.at(dau1_index);
                        dau1_mass = tau_mass.at(dau1_index);
                    } else if (pairType == 2 && isBoostedTau) {
                        dau1_pt = boostedtau_pt.at(dau1_index);
                        dau1_mass = boostedtau_mass.at(dau1_index);
                    } else {
                        return {-999., -999., -999., -999.};
                    }
                    
                    if (isBoostedTau) {
                        dau2_pt = boostedtau_pt.at(dau2_index);
                        dau2_mass = boostedtau_mass.at(dau2_index);
                    } else {
                        dau2_pt = tau_pt.at(dau2_index);
                        dau2_mass = tau_mass.at(dau2_index);
                    }
                    return {dau1_pt, dau1_mass, dau2_pt, dau2_mass};
                }
            """)

    def run(self, df):
        branches = [f"dau1_pt{self.lep_syst}", f"dau1_mass{self.lep_syst}", f"dau2_pt{self.lep_syst}", f"dau2_mass{self.lep_syst}"]

        all_branches = df.GetColumnNames()
        if branches[0] in all_branches:
            return df, []

        df = df.Define(f"lepton_values{self.lep_syst}{self.tau_syst}", "get_lepton_values("
            "pairType, isBoostedTau, dau1_index, dau2_index, "
            "Muon_pt{0}, Muon_mass{0}, Electron_pt{1}, Electron_mass{1}, "
            "Tau_pt{2}, Tau_mass{2}, boostedTau_pt{3}, boostedTau_mass{3})".format(self.muon_syst, self.electron_syst, self.tau_syst, "")) # TODO boostedTau systematics

        for ib, branch in enumerate(branches):
            df = df.Define(branch, f"lepton_values{self.lep_syst}{self.tau_syst}[{ib}]")

        return df, branches


def HHLeptonVarRDF(**kwargs):
    """
    Returns the pt and mass of the selected taus possibly including systematics

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    Output branches: dau1/2_pt/mass_{lep_syst}

    Required RDFModules: :ref:`HHLepton_HHLeptonRDF`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHLeptonVarRDF
            path: Tools.Tools.HHLepton
            parameters:
                isMC: self.dataset.process.isMC

    """
    return lambda: HHLeptonVarRDFProducer(**kwargs)


class HHDiTauJetRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(HHDiTauJetRDFProducer, self).__init__(*args, **kwargs)
        if not os.getenv("HH_DITAUJET"):
            os.environ["HH_DITAUJET"] = "true"
            
            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using VInt = const ROOT::RVec<int>&;
                float deltaR(float eta_1, float eta_2, float phi_1, float phi_2) {
                   const float deta = eta_1 - eta_2;
                   const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
                   const float dRsq = std::pow(deta,2) + std::pow(dphi,2);
                   return sqrt(dRsq);
                }
                int pass_ditaujet(
                    int pairType, int isTauTauJetTrigger,
                    Vfloat tau_eta, Vfloat tau_phi, int dau1_index, int dau2_index,
                    Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi)
                {
                    if (pairType != 2 || isTauTauJetTrigger != 1)  // tautau channel, ditau+jet trig
                        return 1;
                    float dau1_eta = tau_eta[dau1_index];
                    float dau2_eta = tau_eta[dau2_index];
                    float dau1_phi = tau_phi[dau1_index];
                    float dau2_phi = tau_phi[dau2_index];
                    for (int ijet = 0; ijet < jet_pt.size(); ijet++) {
                        if (jet_pt[ijet] < 65)  // 60 from the HLT Path + 5
                            continue;
                        // Overlap removal criteria
                        if ((deltaR(dau1_eta, jet_eta[ijet], dau1_phi, jet_phi[ijet]) > 0.5) &&
                                (deltaR(dau2_eta, jet_eta[ijet], dau2_phi, jet_phi[ijet]) > 0.5))
                            return 1;
                    }
                    return 0;
                }
            """)

    def run(self, df):
        df = df.Define("passTauTauJet",
            "pass_ditaujet(pairType, isTauTauJetTrigger, "
                "Tau_eta, Tau_phi, dau1_index, dau2_index, "
                "Jet_pt{0}, Jet_eta, Jet_phi)".format(self.jet_syst))
        return df.Filter("passTauTauJet == 1", "HHDiTauJetRDF"), []


def HHDiTauJetRDF(**kwargs):
    """
    Filters events passing the ditau+jet trigger but not the offline criteria

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    Required RDFModules: :ref:`HHLepton_HHLeptonRDF`.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHDiTauJetRDF
            path: Tools.Tools.HHLepton
            parameters:
                isMC: self.dataset.process.isMC

    """
    return lambda: HHDiTauJetRDFProducer(**kwargs)
        
        

class HHLeptonGenProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            if not os.getenv("_HHLepton"):
                os.environ["_HHLepton"] = "_HHLepton"

                base = "{}/{}/src/Tools/Tools".format(
                    os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
                ROOT.gROOT.ProcessLine(".L {}/interface/HHLeptonInterface.h".format(base))
                ROOT.gInterpreter.Declare("""
                    auto AK8Gen = AK8GenInterface();
                """)

    def run(self, df):
        if not self.isMC:
            return df, []
        
        df = df.Define("")

        variables = ["Ak8_Zbb_matches", "Ak8_Hbb_matches"]

        df = df.Define("ak8_gen_results", "AK8Gen.get_ak8_genmatch_info("
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, "
                            "GenPart_phi, GenPart_mass, "
                            "FatJet_pt{0}, FatJet_eta, FatJet_phi, FatJet_mass{0})"
                            .format(self.jet_syst))

        branches = []
        for var in variables:
            df = df.Define(f"gen{var}", f"ak8_gen_results.{var}")
            branches.append(f"gen{var}")

        return df, branches


def AK8GenRDF(**kwargs):
    """
    Module to store genMatch information between AK8 and  H/Z->bb

    :param isMC: flag of the dataset being MC or data
    :type : bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: AK8GenRDF
            path: Tools.Tools.AK8Gen
            parameters:
                isMC: self.dataset.process.isMC
    """
    return lambda: AK8GenRDFProducer(**kwargs)

        