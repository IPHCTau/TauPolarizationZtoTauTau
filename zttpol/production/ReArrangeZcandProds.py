import os

from typing import Optional
from columnflow.util import maybe_import
from columnflow.production import Producer, producer
#from httcp.production.ReconstructPi0 import reconstructPi0

from zttpol.util import getGenTauDecayMode
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



def convert_to_coffea_p4(zipped_item, typetag : Optional[str]="PtEtaPhiMLorentzVector"):
    return ak.zip(
        zipped_item,
        with_name = typetag,
        behavior  = coffea.nanoevents.methods.vector.behavior,
    )

def getMaxEtaTauStrip(pt):
    temp = 0.20 * np.power(pt, -0.66)
    ref1 = ak.where(temp > 0.15, 0.15, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2

def getMaxPhiTauStrip(pt):
    temp = 0.35 * np.power(pt, -0.71)
    ref1 = ak.where(temp > 0.30, 0.30, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2


def reconstructPi0(
        zcandp4,
        photons,
        method: Optional[str] = "simpleIC"
):
    photons = ak.with_name(photons, "PtEtaPhiMLorentzVector")
    photons_sorted_pt_indices = ak.argsort(photons.pt, ascending=False)
    photons = photons[photons_sorted_pt_indices]

    p4_pi0 = None


    if method == "simpleIC":
        photons_px = ak.sum(photons.px, axis=1)
        photons_py = ak.sum(photons.py, axis=1)
        photons_pt = np.sqrt(photons_px ** 2 + photons_py ** 2)

        #pt_pi0    = photons[:, 0:1].pt
        pt_pi0     = photons_pt
        eta_pi0    = photons[:, 0:1].eta
        phi_pi0    = photons[:, 0:1].phi
        pdgid_pi0  = ak.values_astype(111 * ak.ones_like(eta_pi0), "int64")
        mass_pi0   = 0.135 * ak.ones_like(eta_pi0)
        charge_pi0 = ak.values_astype(ak.zeros_like(eta_pi0), "int32")
        tauidx_pi0 = photons.tauIdx[:,:1]
        
        
        p4_pi0 = convert_to_coffea_p4({
            "pt"    : pt_pi0,
            "eta"   : eta_pi0,
            "phi"   : phi_pi0,
            "mass"  : mass_pi0,
            "pdgId" : pdgid_pi0,
            "charge": charge_pi0,
            "tauIdx": tauidx_pi0,
        })
        
    elif method == "simpleMB":
        # https://indico.cern.ch/event/1289627/contributions/5444322/attachments/2674971/4638385/mbluj_pi0-direction_CPHTTworkshop_June2023.pdf
        # I will try Michals XGB model later
        
        pi0RecoM = 0.136 #approximate pi0 peak from fits in PF paper
        pi0RecoW = 0.013 

        pdgid_pi0  = ak.values_astype(111 * ak.ones_like(photons.pt), "int64")
        mass_pi0   = 0.135 * ak.ones_like(photons.pt)
        charge_pi0 = ak.values_astype(ak.zeros_like(photons.pt), "int32")
        tauidx_pi0 = photons.tauIdx

        #has_atleast_one_photon = ak.num(photons.pt, axis=1) > 0
        #zcandp4 = ak.where(has_atleast_one_photon, zcandp4, zcandp4[:,:0])
        
        #photons_p4 = ak.where(has_atleast_one_photon,
        #                      photons[:,0:1],
        #                      photons[:,:0])
        
        #deta_photons_zcand = (photons_p4).metric_table(zcandp4, metric = lambda a,b: np.abs(a.eta - b.eta))
        #dphi_photons_zcand = (photons_p4).metric_table(zcandp4, metric = lambda a,b: np.abs(a.delta_phi(b)))

        deta_photons_zcand = ak.firsts((photons).metric_table(zcandp4, metric = lambda a,b: np.abs(a.eta - b.eta)), axis=-1)
        dphi_photons_zcand = ak.firsts((photons).metric_table(zcandp4, metric = lambda a,b: np.abs(a.delta_phi(b))), axis=-1)
        
        #maxeta_photons = getMaxEtaTauStrip(photons_p4.pt)
        #maxphi_photons = getMaxPhiTauStrip(photons_p4.pt)

        maxeta_photons = getMaxEtaTauStrip(photons.pt)
        maxphi_photons = getMaxPhiTauStrip(photons.pt)

        mask_photons = ((np.abs(deta_photons_zcand) < maxeta_photons)
                        & (np.abs(dphi_photons_zcand) < maxphi_photons))
        mask_photons = ak.fill_none(mask_photons, False)
        
        mass_pi0 = mass_pi0[mask_photons][:,:1]
        charge_pi0 = charge_pi0[mask_photons][:,:1]
        tauidx_pi0 = tauidx_pi0[mask_photons][:,:1]
        pdgid_pi0  = pdgid_pi0[mask_photons][:,:1]
        
        strip_photons_p4 = convert_to_coffea_p4(
            {
                "pt"   : photons.pt[mask_photons],
                "eta"  : photons.eta[mask_photons],
                "phi"  : photons.phi[mask_photons],
                "mass" : photons.mass[mask_photons],
            }
        )


        has_one_photon = ak.num(strip_photons_p4.pt, axis=1) == 1
        
        strip_photons_p4_pair = ak.combinations(strip_photons_p4, 2, axis=1)
        strip_photons_p4_pair_0, strip_photons_p4_pair_1 = ak.unzip(strip_photons_p4_pair)
        strip_photons_mass = (strip_photons_p4_pair_0 + strip_photons_p4_pair_1).mass
        mass_val = np.abs(strip_photons_mass - pi0RecoM)
        mass_sorted_idx = ak.argsort(mass_val, axis=1)
        strip_photons_p4_pair_sorted = strip_photons_p4_pair[mass_sorted_idx]
        mass_mask = mass_val < 2 * pi0RecoW
        evt_mask_no_pair = ak.sum(mass_mask, axis=1) == 0
        strip_photons_p4_pair_sorted_pass_mass = strip_photons_p4_pair_sorted[mass_mask]
        strip_photons_p4_mass_selected = ak.concatenate([strip_photons_p4_pair_sorted_pass_mass["0"][:,:1],
                                                         strip_photons_p4_pair_sorted_pass_mass["1"][:,:1]], axis=1)

        sel_strip_photons_p4 = ak.where(has_one_photon,
                                        strip_photons_p4,
                                        ak.where(evt_mask_no_pair,
                                                 strip_photons_p4[:,:1],
                                                 strip_photons_p4_mass_selected)
                                    )
        
        #from IPython import embed; embed()
        #sel_strip_pizero_p4 = ak.from_regular(ak.sum(sel_strip_photons_p4, axis=-1)[:,None])

        
        dummy = sel_strip_photons_p4.px[:,:0]

        #dummy = photons[ak.argsort(photons.pt, ascending=False)][:,:1] # CHANGE
        _mask = ak.num(sel_strip_photons_p4.px, axis=1) > 0
        
        px = ak.from_regular(ak.fill_none(ak.sum(sel_strip_photons_p4.px, axis=1), 0.0)[:,None])
        px = ak.where(_mask, px, dummy) # CHANGE 
        py = ak.from_regular(ak.fill_none(ak.sum(sel_strip_photons_p4.py, axis=1), 0.0)[:,None])
        py = ak.where(_mask, py, dummy) # CHANGE 
        pz = ak.from_regular(ak.fill_none(ak.sum(sel_strip_photons_p4.pz, axis=1), 0.0)[:,None])
        pz = ak.where(_mask, pz, dummy) # CHANGE 
        energy = ak.from_regular(ak.fill_none(ak.sum(sel_strip_photons_p4.energy, axis=1), 0.0)[:,None])
        energy = ak.where(_mask, energy, dummy) # CHANGE 
        
        p4 = convert_to_coffea_p4(
            {
                "x": px, "y": py, "z": pz, "t": energy,
            },
            typetag = "LorentzVector",
        )
        
        p4_pi0 = convert_to_coffea_p4(
            {
                "pt": p4.pt,
                "eta": p4.eta,
                "phi": p4.phi,
                "mass": ak.enforce_type(mass_pi0, "var * float32"),
                "pdgId": ak.enforce_type(pdgid_pi0, "var * int64"),
                "charge": ak.enforce_type(charge_pi0, "var * int32"),
                "tauIdx": ak.enforce_type(tauidx_pi0, "var * int16")
            }
        )

        
    return p4_pi0


def getpions(decay_gentau: ak.Array) -> ak.Array :
    """
    ispion_pos = lambda prod: ((prod.pdgId ==  211) | (prod.pdgId ==  321)) # | (prod.pdgId ==  323) | (prod.pdgId ==  325) | (prod.pdgId ==  327) | (prod.pdgId ==  329) )
    ispion_neg = lambda prod: ((prod.pdgId == -211) | (prod.pdgId == -321)) # | (prod.pdgId == -323) | (prod.pdgId == -325) | (prod.pdgId == -327) | (prod.pdgId == -329) )
    """
    ispion_pos = lambda prod: ((prod.pdgId ==  211) | (prod.pdgId ==  321) | (prod.pdgId ==  323) | (prod.pdgId ==  10321) | (prod.pdgId ==  10211))
    ispion_neg = lambda prod: ((prod.pdgId == -211) | (prod.pdgId == -321) | (prod.pdgId == -323) | (prod.pdgId == -10321) | (prod.pdgId == -10211))

    pions_tau  = decay_gentau[(ispion_pos(decay_gentau) | ispion_neg(decay_gentau))]

    has_three_pions = lambda prod : ak.sum((ispion_pos(prod) | ispion_neg(prod)), axis=1) == 3
    has_two_pos_one_neg_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 2) & (ak.sum(ispion_neg(prod), axis=1) == 1)
    has_two_neg_one_pos_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 1) & (ak.sum(ispion_neg(prod), axis=1) == 2)
    get_sorted_pion_indices   = lambda pions: ak.where(has_two_pos_one_neg_pions(pions),
                                                       ak.argsort(pions.pdgId, ascending=True),
                                                       ak.where(has_two_neg_one_pos_pions(pions),
                                                                ak.argsort(pions.pdgId, ascending=False),
                                                                ak.local_index(pions.pdgId)))
    
    pions_tau_sorted_indices = get_sorted_pion_indices(pions_tau)
    sorted_pions_tau = pions_tau[pions_tau_sorted_indices]

    return sorted_pions_tau


def getphotons(decay_tau: ak.Array) -> ak.Array :
    isphoton   = lambda prod: (prod.pdgId == 22)
    photons = decay_tau[isphoton(decay_tau)]
    return photons


def getgenpizeros(decay_gentau: ak.Array) -> ak.Array :
    #ispizero = lambda col: (np.abs(col.pdgId) == 111) | (np.abs(col.pdgId) == 311) | (np.abs(col.pdgId) == 130) | (np.abs(col.pdgId) == 310)
    ispizero = lambda col: ((col.pdgId == 111) 
                            | (col.pdgId == 311) 
                            | (col.pdgId == 130) 
                            | (col.pdgId == 310))
    #| (col.pdgId == 313)
    #| (col.pdgId == 315)
    #| (col.pdgId == 317)
    #| (col.pdgId == 319))
    pizeros_tau = decay_gentau[ispizero(decay_gentau)]

    return pizeros_tau


def presel_decay_pis(zcand, zcand_pi):
    #from IPython import embed; embed()
    dummy = zcand_pi[:,:0]
    mask02 = ak.fill_none(ak.firsts((zcand.decayMode >= 0) & (zcand.decayMode <= 2), axis=1), False)
    mask10 = ak.fill_none(ak.firsts(zcand.decayMode >= 10, axis=1), False)
    zcand_pi = ak.where(mask02,
                        zcand_pi[:,0:1],
                        ak.where(mask10,
                                 zcand_pi[:,0:3],
                                 dummy))
    #from IPython import embed; embed()
    return zcand_pi

def presel_decay_pi0s(zcand, zcand_pi0):
    dummy = zcand_pi0[:,:0]
    mask12 = ak.fill_none(ak.firsts(((zcand.decayMode == 1) | (zcand.decayMode == 2) | (zcand.decayMode == 11)), axis=1), False)
    zcand_pi0 = ak.where(mask12,
                         zcand_pi0[:,0:1],
                         dummy)
    return zcand_pi0


@producer(
    uses={
        "channel_id", 
        "zcand.{pt,eta,phi,mass,decayMode,charge,IPx,IPy,IPz}",
        "zcandprod.{pt,eta,phi,mass,pdgId}",
        optional("zcand.pt_fastMTT"),
        optional("zcand.eta_fastMTT"),
        optional("zcand.phi_fastMTT"),
        optional("zcand.mass_fastMTT"),
    },
)
def reArrangeDecayProducts(
        self: Producer,
        events: ak.Array,
        **kwargs
) :
    zcand      = events.zcand
    zcandprod  = events.zcandprod

    zcand1     = zcand[:, 0:1]
    zcand2     = zcand[:, 1:2]
    #zcand1prod = ak.firsts(zcandprod[:,0:1], axis=1)
    #zcand2prod = ak.firsts(zcandprod[:,1:2], axis=1)

    #from IPython import embed; embed()
    
    zcand1prod = zcandprod[:,0]
    zcand2prod = zcandprod[:,1]

    #zcand1prod_photons = getphotons(zcand1prod)
    #zcand2prod_photons = getphotons(zcand2prod)
    
    zcand1prod_pions = getpions(zcand1prod)
    zcand2prod_pions = getpions(zcand2prod)

    zcand1prod_pizeros = getgenpizeros(zcand1prod)
    zcand2prod_pizeros = getgenpizeros(zcand2prod)

    
    # zcand1 and its decay products
    p4_zcand1     = ak.with_name(zcand1, "PtEtaPhiMLorentzVector")
    p4_zcand1_pi  = ak.with_name(zcand1prod_pions, "PtEtaPhiMLorentzVector")
    #p4_zcand1_pi  = presel_decay_pis(p4_zcand1, p4_zcand1_pi) # safe
    #p4_zcand1_pi0 = reconstructPi0(p4_zcand1, zcand1prod_photons)
    p4_zcand1_pi0 = ak.with_name(zcand1prod_pizeros, "PtEtaPhiMLorentzVector")
    #p4_zcand1_pi0 = presel_decay_pi0s(p4_zcand1, p4_zcand1_pi0) # safe

    # zcand2 and its decay products
    p4_zcand2     = ak.with_name(zcand2, "PtEtaPhiMLorentzVector")
    p4_zcand2_pi  = ak.with_name(zcand2prod_pions, "PtEtaPhiMLorentzVector")
    #p4_zcand2_pi  = presel_decay_pis(p4_zcand2, p4_zcand2_pi)	# safe 
    #p4_zcand2_pi0 = reconstructPi0(p4_zcand2, zcand2prod_photons)
    p4_zcand2_pi0 = ak.with_name(zcand2prod_pizeros, "PtEtaPhiMLorentzVector")
    #p4_zcand2_pi0 = presel_decay_pi0s(p4_zcand2, p4_zcand2_pi0)	# safe  
    
    zcand1AndProds = ak.concatenate([p4_zcand1, p4_zcand1_pi, p4_zcand1_pi0], axis=1)
    zcand2AndProds = ak.concatenate([p4_zcand2, p4_zcand2_pi, p4_zcand2_pi0], axis=1)

    #from IPython import embed; embed()

    return events, {"p4h1"       : p4_zcand1, 
                    "p4h1pi"     : p4_zcand1_pi, 
                    "p4h1pi0"    : p4_zcand1_pi0, 
                    "p4h2"       : p4_zcand2, 
                    "p4h2pi"     : p4_zcand2_pi, 
                    "p4h2pi0"    : p4_zcand2_pi0}

@producer(
    uses={
        "channel_id",
        "GenTau.*",
        "GenTauProd.*"
    },
    #produces={
    #    "GenTau.decayMode",
    #},
    mc_only=True,
)
def reArrangeGenDecayProducts(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, dict] :
    gzcand       = events.GenTau
    gzcandprod   = events.GenTauProd

    #gzcandprod_indices = gzcand.distinctChildrenIdxG
    #gzcandprod         = events.GenPart._apply_global_index(gzcandprod_indices)

    #zcandprod_dm = getGenTauDecayMode(gzcandprod)
    #events = set_ak_column(events, "GenTau.decayMode", zcandprod_dm)
    
    
    zcand1     = gzcand[:, 0:1]
    zcand2     = gzcand[:, 1:2]
    #zcand1prod = ak.firsts(gzcandprod[:,0:1], axis=1)
    #zcand2prod = ak.firsts(gzcandprod[:,1:2], axis=1)
    zcand1prod = gzcandprod[:,0]
    zcand2prod = gzcandprod[:,1]

    
    zcand1prod_pions = getpions(zcand1prod)
    zcand2prod_pions = getpions(zcand2prod)

    zcand1prod_pizeros = getgenpizeros(zcand1prod)
    zcand2prod_pizeros = getgenpizeros(zcand2prod)

    # zcand1 and its decay products
    p4_zcand1     = ak.with_name(zcand1, "PtEtaPhiMLorentzVector")
    p4_zcand1_pi  = ak.with_name(zcand1prod_pions, "PtEtaPhiMLorentzVector")
    p4_zcand1_pi  = presel_decay_pis(p4_zcand1, p4_zcand1_pi) # safe      
    p4_zcand1_pi0 = ak.with_name(zcand1prod_pizeros, "PtEtaPhiMLorentzVector")
    p4_zcand1_pi0 = presel_decay_pi0s(p4_zcand1, p4_zcand1_pi0) # safe  
    
    # zcand2 and its decay products
    p4_zcand2     = ak.with_name(zcand2, "PtEtaPhiMLorentzVector")
    p4_zcand2_pi  = ak.with_name(zcand2prod_pions, "PtEtaPhiMLorentzVector")
    p4_zcand2_pi  = presel_decay_pis(p4_zcand2, p4_zcand2_pi)	# safe     
    p4_zcand2_pi0 = ak.with_name(zcand2prod_pizeros, "PtEtaPhiMLorentzVector")
    p4_zcand2_pi0 = presel_decay_pi0s(p4_zcand2, p4_zcand2_pi0)	# safe  
    
    zcand1AndProds = ak.concatenate([p4_zcand1, p4_zcand1_pi, p4_zcand1_pi0], axis=1)
    zcand2AndProds = ak.concatenate([p4_zcand2, p4_zcand2_pi, p4_zcand2_pi0], axis=1)

    #from IPython import embed; embed()

    return events, {"p4h1"        : p4_zcand1, 
                    "p4h1pi"      : p4_zcand1_pi, 
                    "p4h1pi0"     : p4_zcand1_pi0, 
                    "p4h2"        : p4_zcand2, 
                    "p4h2pi"      : p4_zcand2_pi, 
                    "p4h2pi0"     : p4_zcand2_pi0}
