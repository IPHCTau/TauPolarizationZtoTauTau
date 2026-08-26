# coding: utf-8

"""
Tau scale factor production.
"""

from __future__ import annotations

import functools

import law

import dataclasses

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, load_correction_set, DotDict
from columnflow.columnar_util import set_ak_column, flat_np_view, layout_ak_array
from columnflow.types import Any

ak = maybe_import("awkward")
np = maybe_import("numpy")


# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@dataclasses.dataclass
class TauSFConfig:
    # vsJet
    tagger_vsJet: str
    correction_vsJet: str = ""
    syst_vsJet: list[str] = dataclasses.field(default_factory=list)

    # vsEle
    tagger_vsEle: str = ""
    correction_vsEle: str = ""
    syst_vsEle: list[str] = dataclasses.field(default_factory=list)

    # vsMu
    tagger_vsMu: str = ""
    correction_vsMu: str = ""
    syst_vsMu: list[str] = dataclasses.field(default_factory=list)

    @classmethod
    def new(cls, obj):
        if isinstance(obj, cls):
            return obj
        if isinstance(obj, dict):
            return cls(**obj)
        if isinstance(obj, (list, tuple)):
            return cls(*obj)

        raise ValueError(f"cannot convert {obj} to {cls.__name__}")


@producer(
    uses={
        "channel_id",
        "Tau.{mass,pt,eta,phi,decayMode,genPartFlav}",
    },
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_tau_file=(lambda self, external_files: external_files.tau_sf),
    # to get the tau config
    get_tau_config=(lambda self: TauSFConfig.new(self.config_inst.x("tau_sf", self.config_inst.x("tau_sf_names", None)))),  # noqa: E501
    weight_name="tau_weight",
    supported_versions={1, 2, 3},
)
def tau_weights(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """

    """
    analysis_channel = self.config_inst.x.channel
    
    tagger_for_gen_tau    = self.tau_config.tagger_vsJet
    tagger_for_ele_to_tau = self.tau_config.tagger_vsEle
    tagger_for_mu_to_tau  = self.tau_config.tagger_vsMu

    # get channels
    ch_etau = self.config_inst.get_channel("etau")
    ch_mutau = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    # channel name and id
    channel_name_id = {
        ch_etau.name : ch_etau.id,
        ch_mutau.name : ch_mutau.id,
        ch_tautau.name : ch_tautau.id
    }
    
    # get taus: one for e/mutau, two for tautau
    etau_mask = events.channel_id == ch_etau.id
    mutau_mask = events.channel_id == ch_mutau.id
    tautau_mask = events.channel_id == ch_tautau.id

    taus = ak.where(
        (etau_mask | mutau_mask),
        events.Tau[:, :1],
        ak.where(
            tautau_mask,
            events.Tau[:, :2],
            events.Tau[:, :0],
        ),
    )
    taus_flat = ak.flatten(taus, axis=1)

    # create a channel id array in the same shape of flat taus
    ch_flat = ak.where(
        etau_mask,
        [[ch_etau.id]],
        ak.where(
            mutau_mask,
            [[ch_mutau.id]],
            ak.where(
                tautau_mask,
                [2 * [ch_tautau.id]],
                [[]],
            ),
        ),
    )
    ch_flat = ak.flatten(ak.values_astype(ch_flat, np.uint8))

    # store some common values
    abseta_flat = abs(taus_flat.eta)
    dm_mask = (
        (taus_flat.decayMode == 0) |
        (taus_flat.decayMode == 1) |
        (taus_flat.decayMode == 10) |
        (taus_flat.decayMode == 11)
    )

    vs_jet_wp = self.config_inst.x.tauIDWPs_config[tagger_for_gen_tau].vs_j[analysis_channel]
    vs_ele_wp = self.config_inst.x.tauIDWPs_config[tagger_for_gen_tau].vs_e[analysis_channel]
    vs_mu_wp  = self.config_inst.x.tauIDWPs_config[tagger_for_gen_tau].vs_m[analysis_channel]


    # helpers to compute scale factors for various tau sources (genuine, fakes) and decay modes
    # genuine taus (separately for decay modes)
    def fill_genuine_tau(sfs_flat: np.array, syst: str, mask: np.array | ak.Array | None = None) -> None:
        genuine_mask = dm_mask & (taus_flat.genPartFlav == 5)
        if mask is not None:
            genuine_mask = genuine_mask & mask

        inputs = {
            "pt": taus_flat.pt[genuine_mask],
            "dm": taus_flat.decayMode[genuine_mask],
            "genmatch": 5,
            "wp": vs_jet_wp,
            "wp_VSe": vs_ele_wp,
            "syst": syst,
            "flag": "dm",
        }
        sfs_flat[genuine_mask] = self.id_vs_jet_corrector.evaluate(
            *(inputs[inp.name] for inp in self.id_vs_jet_corrector.inputs),
        )

        
    # electrons faking taus (separately for decay modes)
    def fill_e_fakes(sfs_flat: np.array, syst: str, mask: np.array | ak.Array | None = None) -> None:
        fake_mask = dm_mask & ((taus_flat.genPartFlav == 1) | (taus_flat.genPartFlav == 3))
        if mask is not None:
            fake_mask = fake_mask & mask
        inputs = {
            "eta": abseta_flat[fake_mask],
            "dm": taus_flat.decayMode[fake_mask],
            "genmatch": taus_flat.genPartFlav[fake_mask],
            "wp": vs_ele_wp,
            "syst": syst,
        }
        sfs_flat[fake_mask] = self.id_vs_e_corrector.evaluate(
            *(inputs[inp.name] for inp in self.id_vs_e_corrector.inputs),
        )

    # muons faking taus
    def fill_mu_fakes(sfs_flat: np.array, syst: str, mask: np.array | ak.Array | None = None) -> None:
        fake_mask = ((taus_flat.genPartFlav == 2) | (taus_flat.genPartFlav == 4))
        if mask is not None:
            fake_mask = fake_mask & mask
        inputs = {
            "eta": abseta_flat[fake_mask],
            "genmatch": taus_flat.genPartFlav[fake_mask],
            "wp": vs_mu_wp,
            "wp_VSe": vs_ele_wp,
            "wp_VSjet": vs_jet_wp,
            "syst": syst,
        }
        sfs_flat[fake_mask] = self.id_vs_mu_corrector.evaluate(
            *(inputs[inp.name] for inp in self.id_vs_mu_corrector.inputs),
        )

    
    # helper to reshape sfs_flat to the shape of taus, multiply across tau axis and store the results
    def add_weight(events: ak.Array, weight_name: str, sfs_flat: np.array) -> ak.Array:
        sfs = layout_ak_array(sfs_flat, taus.pt)
        events = set_ak_column_f32(events, weight_name, ak.prod(sfs, axis=1, mask_identity=False))
        return events

    # prepare per-tau scale factors, starting with ones, then fill values for specific tau sources
    sfs_flat = np.ones(len(taus_flat), dtype=np.float32)
    fill_genuine_tau(sfs_flat, "nom")
    fill_e_fakes(sfs_flat, "nom")
    fill_mu_fakes(sfs_flat, "nom")
    events = add_weight(events, f"{self.weight_name}", sfs_flat)


    # variations
    for direction in ["up", "down"]:
        syst_vsJet = self.tau_config.syst_vsJet
        for idm_syst in syst_vsJet:
            _sfs_flat = sfs_flat.copy()
            fill_genuine_tau(_sfs_flat, f"{idm_syst[-1]}_{direction}", mask=(taus_flat.decayMode == idm_syst[0]))
            events = add_weight(events, f"{self.weight_name}_tau_dm{idm_syst[0]}_{direction}", _sfs_flat)

            
        # electron fakes
        syst_vsEle = self.tau_config.syst_vsEle
        if len(syst_vsEle) == 0:
            for dm in [0,1,10,11]:
                #for region_name, region_mask in [
                #        ("barrel", (abseta_flat < 1.5)),
                #        ("endcap", (abseta_flat >= 1.5)),
                #]:
                _sfs_flat = sfs_flat.copy()
                #fill_e_fakes(_sfs_flat, direction, mask=(region_mask & taus_flat.decayMode == dm))
                #events = add_weight(events, f"{self.weight_name}_e_dm{dm}_{region_name}_{direction}", _sfs_flat)
                fill_e_fakes(_sfs_flat, direction, mask=(taus_flat.decayMode == dm))
                events = add_weight(events, f"{self.weight_name}_e_dm{dm}_{direction}", _sfs_flat)
                
        # muon fakes
        syst_vsMu = self.tau_config.syst_vsMu
        if len(syst_vsMu) == 0:
            for region_name, region_mask in [
                    ("0p0To0p4", (abseta_flat < 0.4)),
                    ("0p4To0p8", ((abseta_flat >= 0.4) & (abseta_flat < 0.8))),
                    ("0p8To1p2", ((abseta_flat >= 0.8) & (abseta_flat < 1.2))),
                    ("1p2To1p7", ((abseta_flat >= 1.2) & (abseta_flat < 1.7))),
                    ("1p7To2p3", (abseta_flat >= 1.7)),
            ]:
                _sfs_flat = sfs_flat.copy()
                fill_mu_fakes(_sfs_flat, direction, mask=region_mask)
                events = add_weight(events, f"{self.weight_name}_mu_{region_name}_{direction}", _sfs_flat)


    return events


@tau_weights.init
def tau_weights_init(self: Producer, **kwargs) -> None:
    super(tau_weights, self).init_func(**kwargs)

    # add the product of nominal and up/down variations to produced columns
    self.produces.add(f"{self.weight_name}")
    for unc in ["tau_dm{0,1,10,11}",
                #"e_dm{0,1,10,11}_{barrel,endcap}",
                "e_dm{0,1,10,11}",
                "mu_{0p0To0p4,0p4To0p8,0p8To1p2,1p2To1p7,1p7To2p3}"]:
        self.produces.add(f"{self.weight_name}_{unc}_{{up,down}}")
    

@tau_weights.requires
def tau_weights_requires(self: Producer, task: law.Task, reqs: dict, **kwargs) -> None:
    super(tau_weights, self).requires_func(task=task, reqs=reqs, **kwargs)

    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)


    
@tau_weights.setup
def tau_weights_setup(
        self: Producer,
        task: law.Task,
        reqs: dict[str, DotDict[str, Any]],
        inputs: dict[str, Any],
        reader_targets: law.util.InsertableDict,
        **kwargs,
) -> None:
    super(tau_weights, self).setup_func(task=task, reqs=reqs, inputs=inputs, reader_targets=reader_targets, **kwargs)

    bundle = reqs["external_files"]

    # load the corrector
    correction_set = load_correction_set(self.get_tau_file(bundle.files))

    # load the config
    self.tau_config = self.get_tau_config()

    self.id_vs_jet_corrector = correction_set[f"{self.tau_config.tagger_vsJet}VSjet"]
    self.id_vs_e_corrector = correction_set[f"{self.tau_config.tagger_vsEle}VSe"]
    self.id_vs_mu_corrector = correction_set[f"{self.tau_config.tagger_vsMu}VSmu"]

    # check versions
    assert self.id_vs_jet_corrector.version in {1, 2, 3, 4}, f"unsupported vs_jet: {self.id_vs_jet_corrector.version}"
    assert self.id_vs_e_corrector.version in {1}, f"unsupported vs_e: {self.id_vs_e_corrector.version}"
    assert self.id_vs_mu_corrector.version in {1}, f"unsupported vs_mu: {self.id_vs_mu_corrector.version}"
