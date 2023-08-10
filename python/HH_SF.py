from Corrections.JME.PUjetID_SF import PUjetID_SFRDFProducer
from Corrections.BTV.btag_SF import btag_SFRDFProducer


class HHPUjetID_SFRDFProducer(PUjetID_SFRDFProducer):
    def __init__(self, year, *args, **kwargs):
        super(HHPUjetID_SFRDFProducer, self).__init__(year, *args, **kwargs)
        self.lep_pt = "{dau1_pt%s, dau2_pt%s}" % (self.systs, self.systs)
        self.lep_eta = "{dau1_eta, dau2_eta}"
        self.lep_phi = "{dau1_phi, dau2_phi}"
        self.lep_mass = "{dau1_mass%s, dau2_mass%s}" % (self.systs, self.systs)


def HHPUjetID_SFRDF(**kwargs):
    """
    Module to compute PU Jet Id scale factors for the HH->bbtautau analysis.
    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHPUjetID_SFRDF
            path: Tools.Tools.HH_SF
            parameters:
                year: self.config.year
                isMC: self.dataset.process.isMC
                isUL: self.dataset.has_tag('ul')
                ispreVFP: self.config.get_aux("isPreVFP", False)

    """
    year = kwargs.pop("year")
    return lambda: HHPUjetID_SFRDFProducer(year, **kwargs)


class HHbtag_SFRDFProducer(btag_SFRDFProducer):
    def __init__(self, year, *args, **kwargs):
        super(HHbtag_SFRDFProducer, self).__init__(year, *args, **kwargs)
        self.lep_pt = "{dau1_pt%s, dau2_pt%s}" % (self.systs, self.systs)
        self.lep_eta = "{dau1_eta, dau2_eta}"
        self.lep_phi = "{dau1_phi, dau2_phi}"
        self.lep_mass = "{dau1_mass%s, dau2_mass%s}" % (self.systs, self.systs)


def HHbtag_SFRDF(**kwargs):
    """
    Module to obtain btagging deepJet SFs with their uncertainties
    for the HH->bbtautau analysis.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHbtag_SFRDF
            path: Tools.Tools.HH_SF
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                isUL: self.dataset.has_tag('ul')
                reshape_uncertainties: [central, ...]

    """
    year = kwargs.pop("year")
    return lambda: HHbtag_SFRDFProducer(year, **kwargs)