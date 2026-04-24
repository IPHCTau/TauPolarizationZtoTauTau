import law
import functools

from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, layout_ak_array
from columnflow.columnar_util import optional_column as optional

from columnflow.util import maybe_import, safe_div

from zttpol.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)



@producer(
    uses={
        # custom columns created upstream, probably by a selector
        "Tau.{pt,eta,phi,mass,genPartFlav,decayModeHPS,decayMode}",
    },
    produces={
        "tau_weight",
    } | {
        f"tau_weight_{direction}"
        for direction in ["up", "down"]
    },
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_tau_file=(lambda self, external_files: external_files.tau_sf),
    # for genuine tauSF from IC
    get_genuine_tau_file=(lambda self, external_files: external_files.gen_tau_sf), 
    # for jet leg
    #get_jetleg_file=(lambda self, external_files: external_files.ditau_jet_trig_sf.path),
    # trigger SF
    #get_tau_trig_file=(lambda self, external_files: external_files.tau_trig_sf),
    # function to determine the tau tagger name
    get_tau_tagger=(lambda self: self.config_inst.x.deep_tau_tagger),
)
def tau_id_weights(self: Producer, events: ak.Array, do_syst: bool, **kwargs) -> ak.Array:
    """
    Producer for tau ID weights. Requires an external file in the config under ``tau_sf``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "tau_sf": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/TAU/2017_UL/tau.json.gz",  # noqa
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files.

    The name of the tagger should be given as an auxiliary entry in the config:

    .. code-block:: python

        cfg.x.tau_tagger = "DeepTau2017v2p1"

    It is used to extract correction set names such as "DeepTau2017v2p1VSjet". *get_tau_tagger* can
    be adapted in a subclass in case it is stored differently in the config.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """
    wp_config = self.config_inst.x.deep_tau_info[self.config_inst.x.deep_tau_tagger]

    # get channels from the config
    #ch_etau = self.config_inst.get_channel("etau")
    #ch_mutau = self.config_inst.get_channel("mutau")
    #ch_tautau = self.config_inst.get_channel("tautau")

    # helper to bring a flat sf array into the shape of taus, and multiply across the tau axis
    reduce_mul = lambda sf: ak.prod(layout_ak_array(sf, events.Tau.pt), axis=1, mask_identity=False)
    #shape_sf = lambda sf: ak.prod(ak.unflatten(sf, ak.num(events.Tau.pt, axis=1)), axis=1, mask_identity=False)
    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayModeHPS, axis=1)
    dmPNet = flat_np_view(events.Tau.decayMode, axis=1)
    genmatch = flat_np_view(events.Tau.genPartFlav, axis=1)


    channel_id_brdcast, _  = ak.broadcast_arrays(events.channel_id[:,None], events.Tau.genPartFlav)
    channel_id_flat        = flat_np_view(channel_id_brdcast, axis=1)

    
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }

    dm_mask = (
        (events.Tau.decayModeHPS == 0) |
        (events.Tau.decayModeHPS == 1) |
        (events.Tau.decayModeHPS == 10)|
        (events.Tau.decayModeHPS == 11)
    )
    dmPNet_mask = (
        (events.Tau.decayMode == 0) |
        (events.Tau.decayMode == 1) |
        (events.Tau.decayMode == 2) |
        (events.Tau.decayMode == 10)|
        (events.Tau.decayMode == 11)
    )
  
    # start with ones
    sf_nom = np.ones_like(pt, dtype=np.float32)
    #sf_nom_evt = np.ones_like(flat_jet_1_pt, dtype=np.float32)

    
    # helpers to create corrector arguments
    if self.id_vs_e_corrector.version == 0:
        args_vs_e   = lambda mask, id_vs_e_wp, syst  : (abseta[mask], genmatch[mask], id_vs_e_wp, syst)
    elif self.id_vs_e_corrector.version in (1,):
        args_vs_e   = lambda mask, id_vs_e_wp, syst  : (abseta[mask], dm[mask], genmatch[mask], id_vs_e_wp, syst)
    else:
        raise NotImplementedError
        
    args_vs_mu  = lambda mask, id_vs_m_wp, syst  : (abseta[mask], genmatch[mask], id_vs_m_wp, syst)    
    #args_vs_jet = lambda mask, id_vs_j_wp, id_vs_e_wp, syst : (pt[mask], dm[mask], genmatch[mask], id_vs_j_wp, id_vs_e_wp, syst, "dm")
    args_vs_jet = lambda mask, id_vs_j_wp, id_vs_e_wp, syst : (pt[mask], dmPNet[mask], genmatch[mask], id_vs_j_wp, id_vs_e_wp, syst, "dm")
    

    shifts = ["nom"]
    if  do_syst:
        shifts=[*shifts,"up", "down"]

    sf_values = {}
    sf_values_dm = {}
    sf_values_e = {}
    sf_values_mu = {}

    
    for the_shift in shifts:
        
        sf_values[the_shift] = sf_nom.copy()
        
        ##############################################################
        # ------------------- for ele -> tau fake ------------------ #
        ##############################################################

        e_mask = ((events.Tau.genPartFlav == tau_part_flav["prompt_e"]) | (events.Tau.genPartFlav == tau_part_flav["tau->e"]))
        e_mask = e_mask & (events.Tau.decayModeHPS != 5) & (events.Tau.decayModeHPS != 6)
        e_mask = flat_np_view(e_mask, axis=1)
        sf_values[the_shift][e_mask] = self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask, wp_config["vs_e"][self.config_inst.x.channel], the_shift))
        

        ##############################################################
        # ------------------ for muon -> tau fake ------------------ #
        ##############################################################

        #from IPython import embed; embed()
        
        mu_mask = ((events.Tau.genPartFlav == tau_part_flav["prompt_mu"]) | (events.Tau.genPartFlav == tau_part_flav["tau->mu"]))
        mu_mask = flat_np_view(mu_mask, axis=1)
        sf_values[the_shift][mu_mask] = self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask, wp_config["vs_m"][self.config_inst.x.channel], the_shift))

        ##############################################################
        # -------------------- for genuine taus -------------------- #
        ##############################################################

        tau_mask = (events.Tau.genPartFlav == tau_part_flav["tau_had"]) & dmPNet_mask
        tau_mask = flat_np_view(tau_mask, axis=1)

        if the_shift == "nom":
            sf_values[the_shift][tau_mask] = self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,
                                                                                            wp_config["vs_j"][self.config_inst.x.channel],
                                                                                            wp_config["vs_e"][self.config_inst.x.channel],
                                                                                            the_shift))

        else:
            for idm in [0, 1, 2, 10, 11]:
                tau_dmX_mask = tau_mask & (dmPNet == idm)
                ch_id_tau_dmX_mask = channel_id_flat[tau_dmX_mask]
                shift_name = f"syst_TES_{self.config_inst.campaign.x.year}_{self.config_inst.campaign.x.postfix}_dm{idm}_{the_shift}"
                if idm == 2:
                    shift_name = f"syst_TES_{self.config_inst.campaign.x.year}_{self.config_inst.campaign.x.postfix}_dm1_{the_shift}"

                sf_values[the_shift][tau_dmX_mask] = self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_dmX_mask,
                                                                                                    wp_config["vs_j"][self.config_inst.x.channel],
                                                                                                    wp_config["vs_e"][self.config_inst.x.channel],
                                                                                                    shift_name))
                
        wt_name = "tau_weight" if the_shift == "nom" else f"tau_weight_{the_shift}"
        events = set_ak_column(events, wt_name, reduce_mul(sf_values[the_shift]), value_type=np.float32)

        
    return events



@tau_id_weights.requires
def tau_id_weights_requires(self: Producer,
                            task: law.Task,
                            reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)


@tau_id_weights.setup
def tau_id_weights_setup(self: Producer,
                         task: law.Task,
                         reqs: dict, inputs: dict, reader_targets: law.util.InsertableDict) -> None:
    bundle = reqs["external_files"]

    # create the trigger and id correctors
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_tau_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    genuine_tau_correction_set = correctionlib.CorrectionSet.from_string(
        self.get_genuine_tau_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    #trig_correction_set = correctionlib.CorrectionSet.from_string(
    #    self.get_tau_trig_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    #)
    tagger_name = self.get_tau_tagger()
    # id
    #self.id_vs_jet_corrector = correction_set[f"{tagger_name}VSjet"]
    self.id_vs_jet_corrector = genuine_tau_correction_set[f"tau_sf_pt-dm_{tagger_name}VSjet_{self.config_inst.campaign.x.year}_{self.config_inst.campaign.x.postfix}"]
    self.id_vs_e_corrector = correction_set[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector = correction_set[f"{tagger_name}VSmu"]

    # check versions
    assert self.id_vs_jet_corrector.version in (0, 1, 2, 3)
    assert self.id_vs_e_corrector.version in (0, 1)
    assert self.id_vs_mu_corrector.version in (0, 1)



        

@producer(
    uses={
        f"TauSpinner*" 
    },
    produces={
        # Version with _alt are duplicated weights computed with different setup of
        # neutral currents parameterization and can be used to estimate uncertainty
        # of the weighting. However it was neglected in the Run-2 analysis
        "tauspinner_weight",
        "tauspinner_weight_up",   # same as cpeven
        "tauspinner_weight_down", # same as cpodd 
        "tauspinner_weight_cpeven",    # 0
        "tauspinner_weight_cpeven_alt",# 0
        "tauspinner_weight_cpmix",     # 0.25
        "tauspinner_weight_cpmix_alt", # 0.25_alt
        "tauspinner_weight_cpmixm",    # -0.25
        "tauspinner_weight_cpmixm_alt",# -0.25_alt
        "tauspinner_weight_cpalpha0p375",     # 0.375
        "tauspinner_weight_cpalpha0p375_alt", # 0.375_alt
        "tauspinner_weight_cpodd",     # 0.5
        "tauspinner_weight_cpodd_alt", # 0.5
    },
    mc_only=True,
)
def tauspinner_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    # https://github.com/hephysicist/CPinHToTauTau/blob/desy_dev/httcp/production/weights.py#L411
    """
    wt_names     = ['nom',
                    'cpeven',
                    'cpeven_alt',
                    'cpmix',
                    'cpmix_alt',
                    'cpmixm',
                    'cpmixm_alt',
                    'cpalpha0p375',
                    'cpalpha0p375_alt',
                    'cpodd',
                    'cpodd_alt']

    branch_names = ['',
                    'weight_cp_0',
                    'weight_cp_0_alt',
                    'weight_cp_0p25',
                    'weight_cp_0p25_alt',
                    'weight_cp_minus0p25',
                    'weight_cp_minus0p25_alt',
                    'weight_cp_0p375',
                    'weight_cp_0p375_alt',
                    'weight_cp_0p5',
                    'weight_cp_0p5_alt']
    
    weight_map   = zip(wt_names, branch_names)
    
    for (wt_name, branch) in weight_map:
        _name = ""
        if wt_name == "nom":
            #weight = (events.TauSpinner.weight_cp_0p5 + events.TauSpinner.weight_cp_0)/2.
            weight = events.TauSpinner.weight_cp_0p25
        else:
            _name = f"_{wt_name}"
            weight = events.TauSpinner[branch]

        #buf = ak.to_numpy(weight)
        #if any(np.isnan(buf)):
        #    warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
        #    buf[np.isnan(buf)] = 0
        #    weight = buf

        weight = ak.nan_to_num(weight, nan=0.0) # probably creating problems while producing datacards
        #weight = ak.nan_to_num(weight, nan=0.001) # Setting very low value to get same rate in all spinner shgfts
        nan_wts = ak.sum((weight == 0.0))
        if nan_wts > 0:
            logger.critical(f"{nan_wts} events with NaN values of {wt_name} out of {len(events)} events, and those NaNs are replaced by 0.0")

        #filter_effs = self.config_inst.x.signal_filter_efficiency
        #if self.dataset_inst.name in filter_effs.keys():
        #    weight *= filter_effs[self.dataset_inst.name]
        #else:
        #    logger.warning(f"No filter efficiency is registered for {self.dataset_inst.name}")
            
        events = set_ak_column_f32(events, f"tauspinner_weight{_name}", weight)
        if _name == "_cpeven": # redundant, needs to be resolved later
            events = set_ak_column_f32(events, f"tauspinner_weight_up", weight)
        elif _name == "_cpodd": # redundant, needs to be resolved later
            events = set_ak_column_f32(events, f"tauspinner_weight_down", weight)

    return events

