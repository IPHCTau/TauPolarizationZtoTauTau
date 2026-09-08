# coding: utf-8

"""
Lepton pair selection as possible Z-Candidate (τ_µ - τ_h)
Requirements:
  select_base
Workflow:
  call base selector, get leptons, jets from there
  mutau_selection : to get the µτ pair
  assign channel id
  save trigger info as new column
  extra lepton veto
  build Zcand i.e. assigning some properties
  build ZandProds
  save process ids for to be stitched datasets, if enabled
  return events and save the electron and tau indices    
"""

import law
import order as od

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
#from columnflow.selection.stats import increment_stats

from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import (
    EMPTY_FLOAT, Route, set_ak_column, IF_DATA, IF_MC
)
from columnflow.columnar_util import optional_column as optional

from zttpol.selection.physics_objects import (
    gentau_selection, jet_cleaning
)
from zttpol.selection.select_base import select_base #, custom_increment_stats
from zttpol.selection.lepton_veto import extra_lepton_veto
from zttpol.selection.match_trigobj import match_trigobjs_dilep_diff
from zttpol.selection.zcand import selzcand, selzcandprod
from zttpol.selection.debug import debug_main
from zttpol.selection.stats import custom_increment_stats

#from zttpol.production.columnvalid import make_column_valid

from zttpol.util import (
    transverse_mass,
    transverse_mass_emu,
    D_zeta,
    IF_DATASET_IS_DY,
    IF_DATASET_IS_W,
    IF_DATASET_IS_SIGNAL,
    IF_RUN2, IF_RUN3,
    IF_DATASET_HAS_LHE_WEIGHTS,
    make_finite,
)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
hist = maybe_import("hist")


#from zttpol.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep


def sort_pairs(dtrpairs: ak.Array,)->ak.Array:
    # Just to get the indices
    # Redundatnt as already sorted by their isolation
    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso04_all, ascending=True)
    # Sort the pairs based on pfRelIso04_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso04_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pfRelIso04_all[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pfRelIso04_all[:,1:2], axis=1),
        False)

    # Sort the pairs based on pt if pfRelIso04_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)

    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False)
    
    # if so, sort the pairs with electrons pfRelIso04_all
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pfRelIso04_all, ascending=True),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    
    # check if the first two pairs have electrons with same pfRelIso04_all
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].pfRelIso04_all[:,:1], axis=1) == ak.firsts(dtrpairs["1"].pfRelIso04_all[:,1:2], axis=1),
        False)
    
    # Sort the pairs based on pt if pfRelIso04_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    return dtrpairs


@selector(
    uses={
        # Muon
        "Muon.{pt,eta,phi,mass,charge,pfRelIso04_all}",
        # Electron
        "Electron.{pt,eta,phi,mass,charge,pfRelIso04_all}",
        # MET
        "PuppiMET.{pt,phi}",
    },
    exposed=False,
)
def emu_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    results = SelectionResult()

    # mus and eles e.g.
    # mus:  [ [m1], [m1,m2], [m1],    [m1], [m1,m2] ]
    # eles: [ [e1], [e1],    [e1,e2], [],   [e1,e2] ]
    
    mus  = events.Muon[lep1_indices]
    eles  = events.Electron[lep2_indices]
    met = events.PuppiMET

    # Sorting leptons by isolation [ascending]
    eles_sort_idxs = ak.argsort(eles.pfRelIso04_all, axis=-1, ascending=True)
    eles = eles[eles_sort_idxs]
    mus_sort_idxs = ak.argsort(mus.pfRelIso04_all, axis=-1, ascending=True)
    mus = mus[mus_sort_idxs]

    # pair of leptons: probable higgs candidate -> leps_pair
    # e.g. [ [(e1,m1)],
    #        [(e1,m1),(e1,m2)],
    #        [(e1,m1),(e2,m1)],
    #        [],
    #        [(e1,m1),(e1,m2),(e2,m1),(e2,m2)]
    #      ]
    
    leps_pair  = ak.cartesian([mus, eles], axis=1)
    
    # unzip to get individuals
    # e.g.
    # lep1 -> lep_pair["0"] -> [ [m1],
    #                            [m1,m1],
    #                            [m1,m2],
    #                            [],
    #                            [m1,m1,m2,m2]
    #                          ]
    # lep2 -> lep_pair["1"] -> [ [e1],
    #                            [e1,e2],
    #                            [e1,e1],
    #                            [],
    #                            [e1,e2,e1,e2]
    #                          ]
    lep1, lep2 = ak.unzip(leps_pair)

    npair = ak.num(leps_pair, axis=1)
    results += SelectionResult(steps = {'emu pairs pre pre-selection': npair >= 1})

    invm = (1*lep1 + 1*lep2).mass
    preselection = {
        #"emu_is_os"         : (lep1.charge * lep2.charge) < 0,
        "emu_dr_0p3"        : (1*lep1).delta_r(1*lep2) > 0.3,
        "emu_mT_60"         : transverse_mass_emu(lep1, lep2,  met) < 60,
        "emu_dzeta_35"      : D_zeta(lep1, lep2, met) > -35,
        "emu_invmass_40"    : (1*lep1 + 1*lep2).mass > 40,
        "emu_invmass_85"    : (1*lep1 + 1*lep2).mass < 85
    }

    # get preselected pairs
    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["emu_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    #good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair = leps_pair[good_pair_mask]
    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["emu_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_dilep_diff(leps_pair, trigger_results)

    pair_selection_steps["emu_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

    lep1, lep2 = ak.unzip(leps_pair)

    # take the 1st pair
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    results += SelectionResult(aux=pair_selection_steps)
    
    return results, leps_pair, trigIds, trigTypes


@selector(
    uses={
        select_base,
        emu_selection,
        extra_lepton_veto,
        selzcand,
        custom_increment_stats,
        jet_cleaning,
        IF_MC("mc_weight"),
    },
    produces={
        select_base,
        emu_selection,
        selzcand,
        "channel_id",
        "cross_e_mu_triggered",
    },
    exposed=True,
)
def select_emu(
        self: Selector,
        events: ak.Array,
        stats: defaultdict,
        hists: DotDict[str, hist.Hist],
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    events,\
        no_sel,\
        results,\
        electron_indices,\
        veto_electron_indices,\
        muon_indices,\
        veto_muon_indices,\
        tau_indices,\
        jet_indices = self[select_base](events)
    
    emu_results, emu_pair, emu_trig_ids, emu_trig_types = self[emu_selection](events,
                                                                              muon_indices,
                                                                              electron_indices,
                                                                              results,
                                                                              call_force=True)
    results += emu_results

    has_one_emu_pair = ak.num(emu_pair.rawIdx, axis=1) == 2
    results += SelectionResult(steps={"one emu pair (sanity check)": has_one_emu_pair})

    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_emu_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    results += SelectionResult(steps={"fall into emu channel (sanity check)": events.channel_id == self.config_inst.get_channel(self.config_inst.x.channel).id})
    
    # muon
    cross_e_mu_triggered = (
        ak.zeros_like(events.event, dtype=bool)
        if ak.sum(ak.num(emu_trig_types, axis=1)) == 0
        else ak.fill_none(ak.any(emu_trig_types == "cross_e_mu", axis=1), False)
    )

    events = set_ak_column(events, "cross_e_mu_triggered", cross_e_mu_triggered)

    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_electron_indices,
                                                                veto_muon_indices,
                                                                emu_pair)
    results += extra_lepton_veto_results

    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, emu_pair)
    results += zcand_results

    events, jet_clean_result, _ = self[jet_cleaning](events,
                                                     jet_indices,
                                                     results,
                                                     **kwargs)
    results += jet_clean_result
    
    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    
    events, results = self[custom_increment_stats](events=events,
                                                   task=kwargs["task"],
                                                   results=results,
                                                   stats=stats,
                                                   no_sel=no_sel,
                                                   event_sel=event_sel)
    
    results += SelectionResult(
        objects = {
            "Electron" : {
                "Electron": emu_pair.rawIdx[:,1:2],
            },
            "Muon" : {
                "Muon": emu_pair.rawIdx[:,0:1],
            },
            "Tau" : {
                "Tau": tau_indices,
            },
        },
    )    

    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

    #events = make_finite(events, debug=self.config_inst.x.verbose.selection.main)
    #events = self[make_column_valid](events, **kwargs)
    
    return events, results
