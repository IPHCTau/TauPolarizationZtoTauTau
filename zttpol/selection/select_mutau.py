# coding: utf-8

"""
Lepton pair selection as possible Z-Candidate (τ_µ - τ_h)
Requirements:
  select_base
Workflow:
  etau_selection : to get the µτ pair
  assign channel id
  save trigger info as new column
  extra lepton veto
  build Zcand i.e. assigning some properties
  build ZandProds
  save process ids for to be stitched datasets, if enabled
  return events and save the electron and tau indices    
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
from zttpol.selection.match_trigobj import match_trigobjs_semilep
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

    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso04_all, ascending=True)

    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pfRelIso04_all[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pfRelIso04_all[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1),
        False
    )

    # if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if rawDeepTau2017v2p1VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    return dtrpairs




@selector(
    uses={
        # muon
        "Muon.{pt,eta,phi,mass,charge,pfRelIso04_all}",
        # tau
        "Tau.{pt,eta,phi,mass,charge,rawDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSe,idDeepTau2018v2p5VSmu}",
        # met
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
    },
    exposed=False,
)
def mutau_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    # lep1 and lep2 e.g.
    # lep1: [ [m1], [m1],    [m1,m2], [],   [m1,m2] ]
    # lep2: [ [t1], [t1,t2], [t1],    [t1], [t1,t2] ]
    muons = events.Muon[lep1_indices] 
    taus  = events.Tau[lep2_indices]

    # Extra channel specific selections on m or tau
    # -------------------- #
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["mutau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["mutau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["mutau"]
    
    is_good_tau     = (
        (taus.pt > 20.0)
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp]) # will be used in categorization
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]
    # -------------------- # 
    
    met = events.MET if self.config_inst.campaign.x.run == 2 else events.PuppiMET

    # Sorting lep1 [Electron] by isolation [ascending]
    muons_sort_idxs = ak.argsort(muons.pfRelIso04_all, axis=-1, ascending=True)
    muons = muons[muons_sort_idxs]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
        
    leps_pair        = ak.cartesian([muons, taus], axis=1)
    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)

    preselection = {
        #"mutau_is_os"         : (lep1.charge * lep2.charge) < 0, # will be used in catgorization
        "mutau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,   # deltaR(lep1, lep2) > 0.5,
        #"mutau_mT_50"         : transverse_mass(lep1, met) < 50, # will be used in catgorization
        "mutau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40,      # invariant_mass(lep1, lep2) > 40
    }

    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["mutau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask
        
    leps_pair = leps_pair[good_pair_mask]

    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["mutau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)
    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_semilep(leps_pair, trigger_results, channel_1='mu')


    #from IPython import embed; embed()

    pair_selection_steps["mutau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

    lep1, lep2 = ak.unzip(leps_pair)


    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    #trigId = trigIds[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)
    

    return SelectionResult(
        aux = pair_selection_steps,
    ), leps_pair, trigIds, trigTypes




@selector(
    uses={
        select_base,
        mutau_selection,
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
        mutau_selection,
        selzcandprod,
        gentau_selection,
        "channel_id",
        "single_triggered",
        "cross_triggered",
    },
    exposed=True,
)
def select_mutau(
        self: Selector,
        events: ak.Array,
        stats: defaultdict,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    events,results,electron_indices,veto_electron_indices,muon_indices,veto_muon_indices,tau_indices,jet_indices = self[select_base](events)
    
    mutau_results, mutau_pair, mutau_trig_ids, mutau_trig_types = self[mutau_selection](events,
                                                                                        muon_indices,
                                                                                        tau_indices,
                                                                                        results,
                                                                                        call_force=True)
    results += mutau_results

    has_one_mutau_pair = ak.num(mutau_pair.rawIdx, axis=1) == 2
    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_mutau_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    
    # muon
    single_mu_triggered = ak.any(mutau_trig_types == 'single_mu', axis=1)
    cross_mu_triggered  = ak.any(mutau_trig_types == 'cross_mu_tau', axis=1)

    
    events = set_ak_column(events, "single_triggered", single_mu_triggered)
    events = set_ak_column(events, "cross_triggered",  cross_mu_triggered)

    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_electron_indices,
                                                                veto_muon_indices,
                                                                mutau_pair)
    results += extra_lepton_veto_results


    print("Extra lepton veto done")
    
    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, mutau_pair)
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
                                                     #ditaujet_jet_indices,
                                                     #call_force=True,
                                                     **kwargs)
    results += jet_clean_result

    #from IPython import embed; embed(); exit()    

    # gen particles info
    # ############################################ #
    # After building the higgs candidates, one can
    # switch on the production of GenTau. Those gentaus
    # will be selected which are matched to the hcand
    # hcand-gentau match = True/False (via config)
    # ############################################ #
    if self.config_inst.x.extra_tags.genmatch:
        #if "is_signal" in list(self.dataset_inst.aux.keys()):
        print(" --->>> zcand-gentau matching")
        events, gentau_results = self[gentau_selection](events, True)
        results += gentau_results


    
    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel

    #events = self[process_ids](events, **kwargs)

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
    

    #from IPython import embed; embed()

    results += SelectionResult(
        objects = {
            "Muon" : {
                "Muon": mutau_pair.rawIdx[:,0:1],
            },
            "Tau" : {
                "Tau": mutau_pair.rawIdx[:,1:2],
            }
        },
    )
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

    
    return events, results

