import os

from analysis_tools.utils import import_root

ROOT = import_root()

class ZZEllipticalCutFilterRDFProducer():
    def __init__(self, *args, **kwargs):
        ROOT.gInterpreter.Declare("""
            bool apply_elliptical_cut_80(float Ztt_svfit_mass, float Zbb_mass) {
                float c_Ztt_svfit_mass = 105.0;
                float r_Ztt_svfit_mass = 51.0;
                float c_Zbb_mass = 118.0;
                float r_Zbb_mass = 113.0;
                float dist_Ztt_svfit_mass = (Ztt_svfit_mass - c_Ztt_svfit_mass) * (Ztt_svfit_mass - c_Ztt_svfit_mass)/(r_Ztt_svfit_mass * r_Ztt_svfit_mass);
                float dist_Zbb_mass = (Zbb_mass - c_Zbb_mass) * (Zbb_mass - c_Zbb_mass)/(r_Zbb_mass * r_Zbb_mass);
                return (dist_Ztt_svfit_mass + dist_Zbb_mass) < 1;
            }
        """)

    def run(self, df):
        # define a new branch to check if the event is inside or outside the ellipse
        df = df.Define("isInside", """apply_elliptical_cut_80(
            Ztt_svfit_mass,
            Zbb_mass
        )""")
        # filter the events outside the ellipse
        df = df.Filter("isInside == 1", "ZZEllipticalCutFilterRDF")
        return df, []

def ZZEllipticalCutFilterRDF(*args, **kwargs):
    return lambda: ZZEllipticalCutFilterRDFProducer(*args, **kwargs)

class BBTauTauFilterRDFProducer():
    def __init__(self, ProcType, isSigBBTT, isBkgBBTT, *args, **kwargs):
        self.isSigBBTT = isSigBBTT
        self.isBkgBBTT = isBkgBBTT
        self.ProcType = ProcType
        # print(" ### DEBUG: isSigBBTT = {}".format(isSigBBTT))
        # print(" ### DEBUG: isBkgBBTT = {}".format(isBkgBBTT))

        if not os.getenv("_BBTauTauFilter"):
            os.environ["_BBTauTauFilter"] = "_BBTauTauFilter"

            ROOT.gInterpreter.Declare("""
                using Vfloat = const ROOT::RVec<float>&;
                using Vint   = const ROOT::RVec<int>&;
                // Return {genPairType, genPart_idx for dau1, genPart_idx for  dau2}
                std::array<int, 3> GenPairType_Zbb_Ztautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother, Vint GenVisTau_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromZ = 0;
                    int n_tau_fromZ = 0;
                    int tau1_id = -1; // index into GenParticle collection of one of the tau
                    int tau2_id = -1;
                    
                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5 /*b*/) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23/*Z*/)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15 /*tau*/) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23 /*Z*/)) {
                            n_tau_fromZ += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromZ > 1)) {
                        int dau1_gen_id = -1; // index into GenParticle collection of electron/muon decaying from gen tau nb 1
                        int dau1_pdg_id = -1;
                        int dau2_gen_id = -1;
                        int dau2_pdg_id = -1;
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            // Looking for children of taus
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = j_gen;
                                dau1_pdg_id = fabs(GenPart_pdgId.at(j_gen)); // electron or muon decaying from tau
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = j_gen;
                                dau2_pdg_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        int genDau1_genPartIdx = -1;
                        int genDau2_genPartIdx = -1;
                        if (((dau1_pdg_id == 11) && (dau2_pdg_id == 13)) || ((dau1_pdg_id == 13) && (dau2_pdg_id == 11))) 
                            GenPairType = 3; // decay modes not covered (e-mu)
                        else if (dau1_pdg_id == 11 || dau2_pdg_id == 11) {
                            GenPairType = 1; // electron
                            genDau1_genPartIdx = std::max(dau1_gen_id, dau2_gen_id); // One of dau1/2_gen_id will be filled, the other -1
                        } else if (dau1_pdg_id == 13 || dau2_pdg_id == 13) {
                            GenPairType = 0; // muon
                            genDau1_genPartIdx = std::max(dau1_gen_id, dau2_gen_id);
                        } else {
                            GenPairType = GenVisTau_genPartIdxMother.size() >= 2 ? 2 : -1; // di-hadronic
                        }
                        
                        for (int i_genvistau = 0; i_genvistau < GenVisTau_genPartIdxMother.size(); i_genvistau ++) {
                            if ((GenVisTau_genPartIdxMother[i_genvistau] == tau1_id || GenVisTau_genPartIdxMother[i_genvistau] == tau2_id)) {
                                // lepton-hadronic case
                                if (GenPairType >=0 && GenPairType <= 1)
                                    genDau2_genPartIdx = i_genvistau;
                                else if (GenPairType == 2) {
                                    // tau_h tau_h case
                                    if (genDau1_genPartIdx == -1)
                                        genDau1_genPartIdx = i_genvistau;
                                    else if (genDau2_genPartIdx == -1)
                                        genDau2_genPartIdx = i_genvistau;
                                    else
                                        std::cerr << "WARNING : Something weird is going on in tau decays" << std::endl;
                                }
                            }
                        }

                        return {GenPairType, genDau1_genPartIdx, genDau2_genPartIdx};
                    }
                    else {
                        return {-1, -1, -1};
                    }
                }
                
                int GenPairType_Zbb_Htautau (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromZ = 0;
                    int n_tau_fromH = 0;
                    int tau1_id = -1;
                    int tau2_id = -1;
                    int dau1_gen_id = -1;
                    int dau2_gen_id = -1;

                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_b_fromZ += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_tau_fromH += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromZ > 1) && (n_tau_fromH > 1)) {
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        if (((dau1_gen_id == 11) && (dau2_gen_id == 13)) || ((dau1_gen_id == 13) && (dau2_gen_id == 11))) GenPairType = 3; // decay modes not covered
                        else if ((dau1_gen_id == 11) || (dau2_gen_id == 11)) GenPairType = 0;
                        else if ((dau1_gen_id == 13) || (dau2_gen_id == 13)) GenPairType = 1;
                        else GenPairType = 2;
                        return GenPairType;
                    }
                    else {
                        return -1;
                    }
                }

                int GenPairType_Ztautau_Hbb (Vint GenPart_pdgId, Vint GenPart_genPartIdxMother) {

                    int GenPairType = -1;
                    int n_b_fromH = 0;
                    int n_tau_fromZ = 0;
                    int tau1_id = -1;
                    int tau2_id = -1;
                    int dau1_gen_id = -1;
                    int dau2_gen_id = -1;

                    for (int i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
                        if (GenPart_genPartIdxMother.at(i_gen) == -1) continue; // it is the incoming parton
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 5) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 25)) {
                            n_b_fromH += 1;
                        }
                        if ((fabs(GenPart_pdgId.at(i_gen)) == 15) && (GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen)) == 23)) {
                            n_tau_fromZ += 1;
                            if (tau1_id == -1) tau1_id = i_gen;
                            else if (tau2_id == -1) tau2_id = i_gen;
                        }
                    }
                    if ((n_b_fromH > 1) && (n_tau_fromZ > 1)) {
                        for (int j_gen = 0; j_gen < GenPart_pdgId.size(); j_gen++) {
                            if ((GenPart_genPartIdxMother.at(j_gen) == tau1_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau1_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                            else if ((GenPart_genPartIdxMother.at(j_gen) == tau2_id) && ((fabs(GenPart_pdgId.at(j_gen)) == 11) || (fabs(GenPart_pdgId.at(j_gen)) == 13))) {
                                dau2_gen_id = fabs(GenPart_pdgId.at(j_gen));
                            }
                        }
                        if (((dau1_gen_id == 11) && (dau2_gen_id == 13)) || ((dau1_gen_id == 13) && (dau2_gen_id == 11))) GenPairType = 3; // decay modes not covered
                        else if ((dau1_gen_id == 11) || (dau2_gen_id == 11)) GenPairType = 0;
                        else if ((dau1_gen_id == 13) || (dau2_gen_id == 13)) GenPairType = 1;
                        else GenPairType = 2;
                        return GenPairType;
                    }
                    else {
                        return -1;
                    }
                }
            """)

    def run(self, df):
        if self.isSigBBTT or self.isBkgBBTT:

            # define a new branch to check if it is or not a ZZ/ZH->bbtautau event
            if self.ProcType == "Zbb_Ztautau":
                print(" ### Running bbtautau filter for Zbb_Ztautau")
                df = df.Define("GenPairType_idx_pair", """GenPairType_Zbb_Ztautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother,
                    GenVisTau_genPartIdxMother
                )""")
            elif self.ProcType == "Zbb_Htautau":
                print(" ### Running bbtautau filter for Zbb_Htautau")
                df = df.Define("GenPairType", """GenPairType_Zbb_Htautau(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            elif self.ProcType == "Ztautau_Hbb":
                print(" ### Running bbtautau filter for Ztautau_Hbb")
                df = df.Define("GenPairType", """GenPairType_Ztautau_Hbb(
                    GenPart_pdgId,
                    GenPart_genPartIdxMother
                )""")
            else:
                raise ValueError("BBTauTauFilterRDF not implemented for self.ProcType = ", self.ProcType)
            
            df = df.Define("GenPairType", "GenPairType_idx_pair[0]")
            df = df.Define("genDau1_genPartIdx", "GenPairType_idx_pair[1]")
            df = df.Define("genDau2_genPartIdx", "GenPairType_idx_pair[2]")
            
            # filter the events with ZZ/ZH->bbtautau
            if self.isSigBBTT:
                # print(" ### DEBUG: isBBTT == 1")
                df = df.Filter("GenPairType != -1", "BBTauTauFilterRDF")
            # filter the events without ZZ/ZH->bbtautau
            elif self.isBkgBBTT:
                # print(" ### DEBUG: isBBTT == 0")
                df = df.Filter("GenPairType == -1", "BBTauTauFilterRDF")
                
            return df, ["GenPairType", "genDau1_genPartIdx", "genDau2_genPartIdx"]
        
        else:
            return df, []

class BBTauTauFilterDummyRDFProducer():
    def run(self, df):
        df = df.Define("GenPairType", -1)
        return df, ["GenPairType"]

def BBTauTauFilterRDF(*args, **kwargs):

    ProcType = kwargs.pop("ProcType")
    isSigBBTT = kwargs.pop("isSigBBTT")
    isBkgBBTT = kwargs.pop("isBkgBBTT")

    return lambda: BBTauTauFilterRDFProducer(ProcType=ProcType, isSigBBTT=isSigBBTT, isBkgBBTT=isBkgBBTT, *args, **kwargs)

    if isSigBBTT or isBkgBBTT:
        return lambda: BBTauTauFilterRDFProducer(ProcType=ProcType, isSigBBTT=isSigBBTT, isBkgBBTT=isBkgBBTT, *args, **kwargs)
    else:
        return lambda: BBTauTauFilterDummyRDFProducer()
