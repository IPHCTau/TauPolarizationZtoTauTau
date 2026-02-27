# coding: utf-8

"""
Lepton pair selection as possible Z-Candidate (τ_e - τ_h)
Requirements:
  select_base
Workflow:
  etau_selection : to get the eτ pair
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

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
#from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

from zttpol.selection.lepton_veto import *
from zttpol.selection.zcand import zcand, zcandprod
from zttpol.selection.custom_stats import custom_increment_stats
from zttpol.selection.debug import debug_main
from zttpol.selection.match_trigobj import match_trigobjs_semilep

from zttpol.production.stitching_NLO import process_ids_dy
from zttpol.production.stitching_LO import process_ids_w
from zttpol.production.extra_weights import scale_mc_weight


from zttpol.util import transverse_mass
from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


from zttpol.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep





def sort_pairs(dtrpairs: ak.Array)->ak.Array:
    # Just to get the indices
    # Redundatnt as already sorted by their isolation
    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso03_all, ascending=True)
    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"].pfRelIso03_all[:,:1], axis=1) == ak.firsts(dtrpairs["0"].pfRelIso03_all[:,1:2], axis=1),
        False)

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
    
    # check if the first two pairs have taus with same rawDeepTau2018v2p5VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1) == ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if rawDeepTau2018v2p5VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]

    return dtrpairs



@selector(
    uses={
        # Electron
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Electron.charge", "Electron.pfRelIso03_all", "Electron.rawIdx",
        # Tau
        optional("Tau.pt"),
        optional("Tau.pt_etau"),
        "Tau.eta", "Tau.phi",
        optional("Tau.mass"),
        optional("Tau.mass_etau"),
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet", "Tau.rawIdx",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
        optional("Tau.genPartFlav"),
        # MET
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
    },
    exposed=False,
)
def etau_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    #events, electron_indices, muon_indices, tau_indices, base_sel_results = events[select_base](events)
    
    eles  = events.Electron[electron_indices]
    taus  = events.Tau[tau_indices]
    
    # Extra channel specific selections on e or tau
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_e_wp         = self.config_inst.x.deep_tau_info[tau_tagger].vs_e["etau"]
    vs_mu_wp        = self.config_inst.x.deep_tau_info[tau_tagger].vs_m["etau"]
    vs_jet_wp       = self.config_inst.x.deep_tau_info[tau_tagger].vs_j["etau"]

    
    is_good_tau     = (
        (taus.pt > 20.0)
        #(taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m[vs_mu_wp])
    )

    taus = taus[is_good_tau]

    # puppi for Run3
    met = events.MET if self.config_inst.campaign.x.year < 2022 else events.PuppiMET

    # Sorting lep1 [Electron] by isolation [ascending]
    eles_sort_idxs = ak.argsort(eles.pfRelIso03_all, axis=-1, ascending=True)
    eles = eles[eles_sort_idxs]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
    
    leps_pair  = ak.cartesian([eles, taus], axis=1)
    
    lep1, lep2         = ak.unzip(leps_pair)

        
    preselection = {
        "etau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "etau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,
        "etau_mT_50"         : transverse_mass(lep1, met) < 50
    }

    # get preselected pairs
    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["etau_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask

    good_pair_mask = ak.fill_none(good_pair_mask, False)
    leps_pair  = leps_pair[good_pair_mask]
    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    pair_selection_steps["etau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    #from IPython import embed; embed()
    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_semilep(leps_pair, trigger_results, channel_1='e')

    pair_selection_steps["etau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

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
        etau_selection,
    },
    produces={
        select_base,
        etau_selection,
        "channel_id",
        "single_triggered",
        "cross_triggered",
    },
    exposed=True,
)
def select_etau(
        self: Selector,
        events: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    events,results,electron_indices,veto_electron_indices,muon_indices,veto_muon_indices,tau_indices,jet_indices = events[select_base](events)
    
    etau_results, etau_pair, etau_trig_ids, etau_trig_types = self[etau_selection](events,
                                                                                   electron_indices,
                                                                                   tau_indices,
                                                                                   results,
                                                                                   call_force=True)
    results += etau_results

    has_one_etau_pair = ak.num(etau_pair.rawIdx, axis=1) == 2
    match_etau_pair_result = SelectionResult(
        steps = {
            "has_one_etau_pair"  : has_one_etau_pair,
        },
    )
    results += match_etau_pair_result

    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_etau_pair, self.config_inst.get_channel("etau"), 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)
    
    
    # ele
    single_e_triggered = ak.any(etau_trigger_types == 'single_e', axis=1)
    cross_e_triggered  = ak.any(etau_trigger_types == 'cross_e_tau', axis=1)

    
    events = set_ak_column(events, "single_triggered", single_e_triggered)
    events = set_ak_column(events, "cross_triggered",  cross_e_triggered)

    Zcands = etau_pair[:,None]
    
    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_ele_indices,
                                                                veto_muon_indices,
                                                                Zcands)
    results += extra_lepton_veto_results
    
    # Zcand results
    events, Zcand_array, Zcand_results = self[zcand](events, Zcands)
    results += Zcand_results

    
    # Zcand prod results
    # -------------------------------------------- #
    # Here, inside hcandprod, the replacement of raw
    # taus by calibrated taus takes place
    # -------------------------------------------- #
    events, Zcandprod_results = self[zcandprod](events, Zcand_array)
    results += Zcandprod_results


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


    outliers_mask_for_stitching = None
    
    if self.process_ids_dy is not None:
        logger.warning("DY stitching is ON : AGAIN")
        events, outliers_mask_for_stitching = self[self.process_ids_dy](events, **kwargs)
    elif self.process_ids_w is not None:
        logger.warning("W stitching is ON : AGAIN")
        events, outliers_mask_for_stitching = self[self.process_ids_w](events, **kwargs)
    else:
        events = self[process_ids](events, **kwargs)

    if outliers_mask_for_stitching is not None:
        n_outliers = ak.sum(outliers_mask_for_stitching)
        if n_outliers > 0:
            logger.warning(f"{self.dataset_inst} has {ak.sum(outliers_mask_for_stitching)} outlier events. Removing those events only")
        outliers_result_for_stitching = SelectionResult(
            steps = {"reject_stitching_outliers": ~outliers_mask_for_stitching}
        )
        results += outliers_result_for_stitching

    
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
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        #group_combinations=group_combinations,
        **kwargs,
    )
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

        return events, SelectionResult(steps = results.steps,
                                       objects = {
                                           "Electron" : {
                                               "Electron": etau_pair.rawIdx[:,0:1],
                                           },
                                           "Tau" : {
                                               "Tau": etau_pair.rawIdx[:,1:2],
                                           }
                                       },
                                       )
    


@select_etau.init
def select_etau_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) is None:
        return

    self.process_ids_dy: process_ids_dy | None = None
    #from IPython import embed; embed()
    if self.dataset_inst.has_tag("is_dy_m50"):
        if self.config_inst.x.allow_dy_stitching:
            #logger.warning("DY stitching is ON")
            #print(f"stitching: {self.config_inst.x.dy_stitching.items()}")
            # check if this dataset is covered by any dy id producer
            for name, dy_cfg in self.config_inst.x.dy_stitching.items():
                #print(f"dataset : {name}, {dy_cfg}")
                dataset_inst = dy_cfg["inclusive_dataset"]
                # the dataset is "covered" if its process is a subprocess of that of the dy dataset
                if dataset_inst.has_process(self.dataset_inst.processes.get_first()):
                    self.process_ids_dy = process_ids_dy.derive(f"process_ids_dy_{name}", cls_dict={
                        "dy_inclusive_dataset": dataset_inst,
                        "dy_leaf_processes": dy_cfg["leaf_processes"],
                    })

                    # add it as a dependency
                    self.uses.add(self.process_ids_dy)
                    self.produces.add(self.process_ids_dy)
                    
                    # stop after the first match
                    break
            
    self.process_ids_w: process_ids_w | None = None
    if self.dataset_inst.has_tag("is_w"):
        if self.config_inst.x.allow_w_stitching:
            #logger.warning("W stitching is ON")
            # check if this dataset is covered by any dy id producer
            for name, w_cfg in self.config_inst.x.w_stitching.items():
                dataset_inst = w_cfg["inclusive_dataset"]
                # the dataset is "covered" if its process is a subprocess of that of the dy dataset
                if dataset_inst.has_process(self.dataset_inst.processes.get_first()):
                    self.process_ids_w = process_ids_w.derive(f"process_ids_w_{name}", cls_dict={
                        "w_inclusive_dataset": dataset_inst,
                        "w_leaf_processes": w_cfg["leaf_processes"],
                    })

                    # add it as a dependency
                    self.uses.add(self.process_ids_w)
                    self.produces.add(self.process_ids_w)

                    # stop after the first match
                    break
