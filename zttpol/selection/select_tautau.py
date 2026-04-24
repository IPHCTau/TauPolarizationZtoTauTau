# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

import law

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

#from columnflow.production.processes import process_ids
from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
#from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

from zttpol.selection.physics_objects import (
    gentau_selection, jet_cleaning
)
from zttpol.selection.select_base import select_base
from zttpol.selection.lepton_veto import extra_lepton_veto
from zttpol.selection.custom_stats import custom_increment_stats
from zttpol.selection.match_trigobj import match_trigobjs_fullhad
from zttpol.selection.zcand import selzcand, selzcandprod
from zttpol.selection.debug import debug_main

#from zttpol.production.stitching_NLO import process_ids_dy
#from zttpol.production.stitching_LO import process_ids_w
#from zttpol.production.extra_weights import scale_mc_weight

from zttpol.util import transverse_mass
from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


#from zttpol.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep



def sort_pairs(dtrpairs: ak.Array)->ak.Array:

    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    
    return dtrpairs



@selector(
    uses={
        "Tau.{pt,eta,phi,mass,rawIdx,charge,rawDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSe,idDeepTau2018v2p5VSmu}",
        optional("Tau.genPartFlav"),
        "Jet.{pt,eta,phi,mass}",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        trigger_results: SelectionResult,
        jet_indices : ak.Array,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    taus            = events.Tau[lep_indices]
    # Extra channel specific selections on tau
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["tautau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["tautau"]

    is_good_tau     = (
        (taus.pt > 20.0)
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    
    # Sorting leps [Tau] by deeptau [descending]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]

    leps_pair  = ak.combinations(taus, 2, axis=1)    
    lep1, lep2 = ak.unzip(leps_pair)

    #from IPython import embed; embed()
    
    preselection = {
        #"tautau_tau1_iso"      : (lep1.idDeepTau2018v2p5VSjet >= tau_tagger_wps.vs_j[vs_jet_wp]),
        "tautau_is_pt_35"      : (lep1.pt > 35.0) & (lep2.pt > 35.0), # just changed 40.0 to 35.0 (19.12.2024)
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        #"tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        "tautau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40.0, # invariant_mass(lep1, lep2) > 40
    }

    
    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    category_selections = {}
    pair_selection_steps["tautau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair = leps_pair[good_pair_mask]

    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["tautau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes, jet_idx = match_trigobjs_fullhad(leps_pair,
                                                                    trigger_results,
                                                                    events.Jet[jet_indices])

    
    pair_selection_steps["tautau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    lep1, lep2 = ak.unzip(leps_pair)

    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
                                    
    #has_jet_ = ak.num(jet_idx, axis=1) > 0
    ditau_ids = trigIds[trigIds == 15151]
    ditaujet_ids = trigIds[trigIds == 15152]

    has_jet = jet_idx >= 0
    a = ak.concatenate([(ditaujet_ids == 15152), has_jet], axis=1)
    b = ak.sum(a, axis=1) == 2
    ditaujet_ids_matched = ak.where(b, ditaujet_ids, ditaujet_ids[:,:0])

    trig_ids_matched = ak.concatenate([ditau_ids,ditaujet_ids_matched], axis=1)
    trig_types_matched = trigTypes[trig_ids_matched > 0]
    
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    sort_idx = ak.argsort(leps_pair.pt, ascending=False)
    leps_pair = leps_pair[sort_idx]

    leps_pair_dummy = leps_pair[:,:0]
    leps_pair_matched = ak.where(ak.num(trig_ids_matched) > 0, leps_pair, leps_pair_dummy)

    #from IPython import embed; embed()
    

    #return SelectionResult(
    #    aux = pair_selection_steps,
    #), leps_pair, trigIds, trigTypes

    return SelectionResult(
        aux = pair_selection_steps,
    ), leps_pair_matched, trig_ids_matched, trig_types_matched, jet_idx





@selector(
    uses={
        select_base,
        tautau_selection,
        extra_lepton_veto,
        selzcand,
        selzcandprod,
        gentau_selection,
        increment_stats,
        custom_increment_stats,
        jet_cleaning,
    },
    produces={
        select_base,
        tautau_selection,
        selzcandprod,
        gentau_selection,
        "channel_id",
        "cross_triggered",
        "cross_jet_triggered"
    },
    exposed=True,
)
def select_tautau(
        self: Selector,
        events: ak.Array,
        stats: defaultdict,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    events,results,electron_indices,veto_electron_indices,muon_indices,veto_muon_indices,tau_indices,jet_indices = self[select_base](events)
    
    tautau_results, tautau_pair, tautau_trig_ids, tautau_trig_types, ditaujet_jet_indices = self[tautau_selection](events,
                                                                                                                   tau_indices,
                                                                                                                   results,
                                                                                                                   jet_indices,
                                                                                                                   call_force=True)
    
    results += tautau_results

    has_one_tautau_pair = ak.num(tautau_pair.rawIdx, axis=1) == 2
    match_tautau_pair_result = SelectionResult(
        steps = {
            "has_one_tautau_pair"  : has_one_tautau_pair,
        },
    )
    results += match_tautau_pair_result

    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_tautau_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)
    
    
    # tau
    cross_tau_triggered = ak.any(tautau_trig_types == 'cross_tau_tau', axis=1)
    cross_tau_jet_triggered  = ak.any(tautau_trig_types == 'cross_tau_tau_jet', axis=1)

    
    events = set_ak_column(events, "cross_triggered", cross_tau_triggered)
    events = set_ak_column(events, "cross_jet_triggered",  cross_tau_jet_triggered)

    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_electron_indices,
                                                                veto_muon_indices,
                                                                tautau_pair)
    results += extra_lepton_veto_results
    
    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, tautau_pair)
    results += zcand_results

    
    # Zcand prod results
    # -------------------------------------------- #
    # Here, inside hcandprod, the replacement of raw
    # taus by calibrated taus takes place
    # -------------------------------------------- #
    events, zcandprod_results = self[selzcandprod](events, zcand_array)
    results += zcandprod_results

    events, jet_clean_result, _ = self[jet_cleaning](events,
                                                     jet_indices,
                                                     results,
                                                     ditaujet_jet_indices,
                                                     #call_force=True,
                                                     **kwargs)
    results += jet_clean_result

    

    # gen particles info
    # ############################################ #
    # After building the higgs candidates, one can
    # switch on the production of GenTau. Those gentaus
    # will be selected which are matched to the hcand
    # hcand-gentau match = True/False (via config)
    # ############################################ #
    if self.config_inst.x.extra_tags.genmatch:
        #if "is_signal" in list(self.dataset_inst.aux.keys()):
        print(" --->>> Zcand-gentau matching")
        events, gentau_results = self[gentau_selection](events, True)
        results += gentau_results


    
    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    
    events, results = self[custom_increment_stats]( 
        events,
        results,
        stats,
    )

    weight_map = {
        "num_filtered_events": Ellipsis,
        "num_events_selected": event_sel,
    }
    group_map = {}
    group_combinations = []
    if self.dataset_inst.is_mc:
        weight_map["sum_filtered_mc_weight"] = events.mc_weight
        weight_map["sum_mc_weight_selected"] = (events.mc_weight, event_sel)
        # groups
        group_map = {
            **group_map,
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
        }
        # combinations
        #group_combinations.append(("process"))

    events, results = self[increment_stats](
        events=events,
        results=results,
        stats=stats,
        weight_map=weight_map,
        group_map=group_map,
        #group_combinations=group_combinations,
        **kwargs,
    )


    results += SelectionResult(
        objects = {
            "Tau" : {
                "Tau": tautau_pair.rawIdx,
            }
        },
    )
    
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)


    #from IPython import embed; embed()


    
    return events, results
