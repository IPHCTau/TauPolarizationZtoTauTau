# coding: utf-8

"""
Process ID producer relevant for the stitching of the DY samples.
"""

import functools

import law

from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

from zttpol.util import transverse_mass, IF_RUN2, IF_RUN3, IF_DATASET_IS_DY,IF_DATASET_IS_W

np = maybe_import("numpy")
ak = maybe_import("awkward")
sp = maybe_import("scipy")
maybe_import("scipy.sparse")


logger = law.logger.get_logger(__name__)

set_ak_column_i64 = functools.partial(set_ak_column, value_type=np.int64)


# ################################# #
#            ABCD masks             #
# ################################# #
@producer(
    uses={
        "channel_id",
        "zcand.*",
        "Jet.pt", "bJet.pt",
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
        #"classifier_score",
    },
    produces={
        "is_os",
        "is_b_veto",
        "is_low_mt",
        "is_high_mt",
        "is_lep_1",
        "is_iso_1", "is_iso_2",
        "is_real_1", "is_real_2",
        "is_fake_1", "is_fake_2",
        "is_pi_1", "is_pi_2",
        "is_rho_1", "is_rho_2",
        "is_a1_1pr_2pi0_1", "is_a1_1pr_2pi0_2",
        "is_a1_3pr_0pi0_1", "is_a1_3pr_0pi0_2",
        "is_a1_3pr_1pi0_1", "is_a1_3pr_1pi0_2",
        "is_ipsig_0to1_1",
        "has_0jet","has_1jet","has_2jet",
        ## classifier based columns
        #"is_tautau_dy",
        #"is_tautau_fake",
        #"is_tautau_higgs",
        ## futher binning of higgs node
        #"is_tautau_higgs_bin_1",
        #"is_tautau_higgs_bin_2",
        #"is_tautau_higgs_bin_3",
        #"is_tautau_higgs_bin_4",
        #"is_tautau_higgs_bin_5",
    },
    exposed=False,
)
def build_abcd_masks(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:

    # get channels from the config
    ch_emu    = self.config_inst.get_channel("emu")
    ch_etau   = self.config_inst.get_channel("etau")
    ch_mutau  = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    zcand = ak.with_name(events.zcand, "PtEtaPhiMLorentzVector")
    z1 = zcand[:,0:1]
    z2 = zcand[:,1:2]

    # OS --> opposite sign
    is_os = (z1.charge * z2.charge) < 0
    is_os = ak.fill_none(ak.any(is_os, axis=1), False)

    # BVETO --> events with / without bjets
    is_b_veto = ak.num(events.bJet.pt, axis=1) == 0

    # LOWMT --> Required for leptonic channels
    met = events.MET if self.config_inst.campaign.x.run == 2 else events.PuppiMET
    met = ak.with_name(met, "PtEtaPhiMLorentzVector")
    is_low_mt = transverse_mass(z1, met) < 65.0
    is_high_mt = ~is_low_mt & (transverse_mass(z1, met) < 200.0)
    
    is_low_mt  = ak.fill_none(ak.any(is_low_mt, axis=1), False)
    is_high_mt = ak.fill_none(ak.any(is_high_mt, axis=1), False)

    is_lep_1 = z1.decayMode < 0
    is_lep_1 = ak.fill_none(ak.any(is_lep_1, axis=1), False)

    # IPSig
    is_ipsig_0to1_1 = np.abs(z1.IPsig) < 1.0
    is_ipsig_0to1_1 = ak.fill_none(ak.any(is_ipsig_0to1_1, axis=1), False)

    is_pi_1 = z1.decayMode[:,0] == 0            # z1 to pion
    is_pi_2 = z2.decayMode[:,0] == 0            # z2 to pion

    is_rho_1 = z1.decayMode[:,0] == 1           # z1 to rho
    is_rho_2 = z2.decayMode[:,0] == 1           # z2 to rho

    is_a1_1pr_2pi0_1 = z1.decayMode[:,0] == 2   # z1 to a1 (DM 2)
    is_a1_1pr_2pi0_2 = z2.decayMode[:,0] == 2   # z2 to a1 (DM 2)

    is_a1_3pr_0pi0_1 = z1.decayMode[:,0] == 10  # z1 to a1 (DM 10)
    is_a1_3pr_0pi0_2 = z2.decayMode[:,0] == 10  # z2 to a1 (DM 10)

    is_a1_3pr_1pi0_1 = z1.decayMode[:,0] == 11  # z1 to a1 (DM 11)
    is_a1_3pr_1pi0_2 = z2.decayMode[:,0] == 11  # z2 to a1 (DM 11)

    # njet categories
    has_0jet = ak.num(events.Jet.pt, axis=1) == 0
    has_1jet = ak.num(events.Jet.pt, axis=1) == 1
    has_2jet = ak.num(events.Jet.pt, axis=1) >= 2


    # tau tagger wp
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_info = self.config_inst.x.deep_tau_info[tau_tagger]

    vs_jet_wp       = lambda tau_tagger_info, ch : tau_tagger_info.wp.vs_j[tau_tagger_info.vs_j[ch]]

    # ISO1 --> required for tau-tau channel only to categorise events on the basis of
    # leading tau isolation
    is_iso_1_dummy = z1.rawIdx < 0
    is_iso_1 = ak.where(events.channel_id == ch_tautau.id,
                        ak.values_astype(z1.isolation, np.int32) >= vs_jet_wp(tau_tagger_info, ch_tautau.name),
                        ak.where(events.channel_id == ch_mutau.id,
                                 z1.isolation <= 0.15,
                                 ak.where(events.channel_id == ch_etau.id,
                                          z1.isolation <= 0.3,
                                          is_iso_1_dummy)
                                 )
                        )
    is_iso_1 = ak.fill_none(ak.any(is_iso_1, axis=1), False)

    # ISO2 --> required for all channels
    id_etau_pass   = z2.isolation >= vs_jet_wp(tau_tagger_info, ch_etau.name)
    id_mutau_pass  = z2.isolation >= vs_jet_wp(tau_tagger_info, ch_mutau.name)
    id_tautau_pass = z2.isolation >= vs_jet_wp(tau_tagger_info, ch_tautau.name)
    
    is_iso_2 = ak.where(events.channel_id == ch_tautau.id,
                        id_tautau_pass,
                        ak.where(events.channel_id == ch_mutau.id,
                                 id_mutau_pass,
                                 ak.where(events.channel_id == ch_etau.id,
                                          id_etau_pass,
                                          is_iso_1_dummy)
                                 )
                        )
    is_iso_2 = ak.fill_none(ak.any(is_iso_2, axis=1), False)

    # REAL1 --> to get the contribution of real MC taus only with genPartFlav > 0
    # only required for tau-tau channel
    # maybe redundant: true info for e and mu. Probably not gonna used
    is_real_1 = is_fake_1 = events.event >= 0

    if self.dataset_inst.is_mc:
        is_real_1 = ak.where(events.channel_id == ch_etau.id,
                             ((z1.genPartFlav == 1) | (z1.genPartFlav == 15) | (z1.genPartFlav == 22)), # true ele
                             ak.where(events.channel_id == ch_mutau.id,
                                      ((z1.genPartFlav == 1) | (z1.genPartFlav == 15)), # true mu
                                      ak.where(events.channel_id == ch_tautau.id,
                                               ((z1.genPartFlav > 0) & (z1.genPartFlav < 6)), # true tau
                                               is_iso_1_dummy)))
        is_real_1 = ak.fill_none(ak.any(is_real_1, axis=1), False)

        is_fake_1 = ak.where(events.channel_id == ch_etau.id,
                             ((z1.genPartFlav == 0) | (z1.genPartFlav == 3) | (z1.genPartFlav == 4) | (z1.genPartFlav == 5)), # fake ele
                             ak.where(events.channel_id == ch_mutau.id,
                                      ((z1.genPartFlav == 0) | (z1.genPartFlav == 3) | (z1.genPartFlav == 4) | (z1.genPartFlav == 5)), # true mu
                                      ak.where(events.channel_id == ch_tautau.id,
                                               ((z1.genPartFlav == 0) | (z1.genPartFlav == 6)), # true tau
                                               is_iso_1_dummy)))
        is_fake_1 = ak.fill_none(ak.any(is_fake_1, axis=1), False)
        
    # REAL2 --> the same for the sub-leading tau
    # valid for all channels
    is_real_2 = is_fake_2 = events.event >= 0    

    if self.dataset_inst.is_mc:
        is_real_2 = (z2.genPartFlav > 0) & (z2.genPartFlav < 6)
        if self.dataset_inst.has_tag("is_w") or self.dataset_inst.has_tag("no_lhe_weights") or self.dataset_inst.has_tag("is_st"):
            # temporary, fake contribution from WJets Simulation, not data driven
            is_real_2 = ak.where(events.channel_id == ch_etau.id,
                                 z2.genPartFlav >= 0,
                                 ak.where(events.channel_id == ch_mutau.id,
                                          z2.genPartFlav >= 0,
                                          is_real_2))
        is_real_2 = ak.fill_none(ak.any(is_real_2, axis=1), False)
        
        is_fake_2 = ((z2.genPartFlav == 0) | (z2.genPartFlav == 6))
        is_fake_2 = ak.fill_none(ak.any(is_fake_2, axis=1), False)
        


    # set columns
    events = set_ak_column(events, "is_os",     is_os)
    events = set_ak_column(events, "is_iso_1",  is_iso_1)
    events = set_ak_column(events, "is_iso_2",  is_iso_2)
    events = set_ak_column(events, "is_real_1", is_real_1)
    events = set_ak_column(events, "is_real_2", is_real_2)
    # real and fake are the same for the data
    events = set_ak_column(events, "is_fake_1", is_fake_1)
    events = set_ak_column(events, "is_fake_2", is_fake_2)
    # mt
    events = set_ak_column(events, "is_low_mt", is_low_mt)
    events = set_ak_column(events, "is_high_mt", is_high_mt)
    # bveto
    events = set_ak_column(events, "is_b_veto", is_b_veto)
    # for CP categories
    events = set_ak_column(events, "is_lep_1",   is_lep_1)
    events = set_ak_column(events, "is_pi_1",   is_pi_1)
    events = set_ak_column(events, "is_pi_2",   is_pi_2)
    events = set_ak_column(events, "is_rho_1",  is_rho_1)
    events = set_ak_column(events, "is_rho_2",  is_rho_2)
    events = set_ak_column(events, "is_a1_1pr_2pi0_1",  is_a1_1pr_2pi0_1)
    events = set_ak_column(events, "is_a1_1pr_2pi0_2",  is_a1_1pr_2pi0_2)
    events = set_ak_column(events, "is_a1_3pr_0pi0_1",  is_a1_3pr_0pi0_1)
    events = set_ak_column(events, "is_a1_3pr_0pi0_2",  is_a1_3pr_0pi0_2)
    events = set_ak_column(events, "is_a1_3pr_1pi0_1",  is_a1_3pr_1pi0_1)
    events = set_ak_column(events, "is_a1_3pr_1pi0_2",  is_a1_3pr_1pi0_2)

    events = set_ak_column(events, "is_ipsig_0to1_1", is_ipsig_0to1_1)

    events = set_ak_column(events, "has_0jet", has_0jet)
    events = set_ak_column(events, "has_1jet", has_1jet)
    events = set_ak_column(events, "has_2jet", has_2jet)


    """
    max_idx = ak.fill_none(ak.argmax(events.classifier_score, axis=1), -1)

    is_tautau_dy    = max_idx == 0
    is_tautau_higgs = max_idx == 1
    is_tautau_fake  = max_idx == 2
    
    events = set_ak_column(events, "is_tautau_dy",    is_tautau_dy)
    events = set_ak_column(events, "is_tautau_fake",  is_tautau_fake)    
    events = set_ak_column(events, "is_tautau_higgs", is_tautau_higgs)

    # further categories with Higgs BDT score
    tautau_higgs_score = ak.fill_none(ak.firsts(events.classifier_score[:,1:2], axis=1), -99.9)
    
    is_tautau_higgs_bin_1 = (tautau_higgs_score >= 0.0)  & (tautau_higgs_score < 0.55)
    is_tautau_higgs_bin_2 = (tautau_higgs_score >= 0.55) & (tautau_higgs_score < 0.65)
    is_tautau_higgs_bin_3 = (tautau_higgs_score >= 0.65) & (tautau_higgs_score < 0.80)
    is_tautau_higgs_bin_4 = (tautau_higgs_score >= 0.80) & (tautau_higgs_score < 0.90)
    is_tautau_higgs_bin_5 = (tautau_higgs_score >= 0.90) & (tautau_higgs_score <= 1.0)

    events = set_ak_column(events, "is_tautau_higgs_bin_1", is_tautau_higgs_bin_1)
    events = set_ak_column(events, "is_tautau_higgs_bin_2", is_tautau_higgs_bin_2)
    events = set_ak_column(events, "is_tautau_higgs_bin_3", is_tautau_higgs_bin_3)
    events = set_ak_column(events, "is_tautau_higgs_bin_4", is_tautau_higgs_bin_4)
    events = set_ak_column(events, "is_tautau_higgs_bin_5", is_tautau_higgs_bin_5)
    """

    return events
