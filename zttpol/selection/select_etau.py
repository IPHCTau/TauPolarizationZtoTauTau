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

from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
#from columnflow.selection.stats import increment_stats

from columnflow.util import maybe_import
from columnflow.columnar_util import (
    EMPTY_FLOAT, Route, set_ak_column, IF_DATA, IF_MC
)
from columnflow.columnar_util import optional_column as optional

from zttpol.selection.physics_objects import (
    gentau_selection, jet_cleaning
)
from zttpol.selection.select_base import select_base
from zttpol.selection.lepton_veto import extra_lepton_veto
from zttpol.selection.match_trigobj import match_trigobjs_semilep
from zttpol.selection.zcand import selzcand, selzcandprod
from zttpol.selection.debug import debug_main
from zttpol.selection.stats import custom_increment_stats

from zttpol.util import (
    transverse_mass,
    IF_DATASET_IS_DY,
    IF_DATASET_IS_W,
    IF_DATASET_IS_SIGNAL,
    IF_RUN2, IF_RUN3,
    IF_DATASET_HAS_LHE_WEIGHTS
)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")




def sort_pairs(
        dtrpairs: ak.Array,
        ele_iso = "pfRelIso03_all",
        tau_iso = "rawDeepTau2018v2p5VSjet"
) -> ak.Array:
    # Just to get the indices
    # Redundatnt as already sorted by their isolation
    sorted_idx = ak.argsort(dtrpairs["0"][ele_iso], ascending=True)
    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"][ele_iso][:,:1], axis=1) == ak.firsts(dtrpairs["0"][ele_iso][:,1:2], axis=1),
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
                          ak.argsort(dtrpairs["1"][tau_iso], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    
    # check if the first two pairs have taus with same rawDeepTau2018v2p5VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"][tau_iso][:,:1], axis=1) == ak.firsts(dtrpairs["1"][tau_iso][:,1:2], axis=1),
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
        "Electron.{pt,eta,phi,mass,charge,pfRelIso03_all,rawIdx}",
        # Tau
        "Tau.{pt,eta,phi,mass,charge,rawDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSe,idDeepTau2018v2p5VSmu,rawIdx}",
        # MET
        "PuppiMET.{pt,phi}",
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

    results = SelectionResult()
    
    eles  = events.Electron[lep1_indices]
    taus  = events.Tau[lep2_indices]
    
    # Extra channel specific selections on m or tau
    vs_e_wp         = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_e["etau"]
    vs_mu_wp        = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_m["etau"]
    vs_jet_wp       = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_j["etau"]
    
    is_good_tau     = (
        (taus.pt > 20.0)
        & (taus.idDeepTau2018v2p5VSe   >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_m[vs_mu_wp])
    )
    # idDeepTau2018v2p5VSjet will be used in the categorization

    taus = taus[is_good_tau]


    # Sorting lep1 [Electron] by isolation [ascending]
    eles_sort_idxs = ak.argsort(eles.pfRelIso03_all, axis=-1, ascending=True)
    eles = eles[eles_sort_idxs]
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
    
    leps_pair  = ak.cartesian([eles, taus], axis=1)    
    lep1, lep2 = ak.unzip(leps_pair)

    npair = ak.num(lep1, axis=1)
    results += SelectionResult(steps = {'etau pairs pre pre-selection': npair >= 1})    
        
    preselection = {
        #"etau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "etau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,
        "etau_mT_60"         : transverse_mass(lep1, events.PuppiMET) < 60
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
    results += SelectionResult(steps = {'etau pairs post pre-selection': npair >= 1})
    pair_selection_steps["etau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)

    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_semilep(leps_pair, trigger_results, channel_1='e')

    npair = ak.num(leps_pair["0"], axis=1)
    results += SelectionResult(steps = {'etau pairs post trigobj match': npair >= 1})    
    pair_selection_steps["etau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

    lep1, lep2 = ak.unzip(leps_pair)


    # take the 1st pair and 1st trigger id
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)

    results += SelectionResult(aux=pair_selection_steps)

    return results, leps_pair, trigIds, trigTypes



@selector(
    uses={
        select_base,
        etau_selection,
        extra_lepton_veto,
        selzcand,
        selzcandprod,
        gentau_selection,
        custom_increment_stats,
        jet_cleaning,
        IF_MC("mc_weight"),
    },
    produces={
        select_base,
        etau_selection,
        selzcandprod,
        gentau_selection,
        "channel_id",
        "single_triggered",
        "cross_triggered",
    },
    exposed=True,
)
def select_etau(
        self: Selector,
        events: ak.Array,
        stats: defaultdict,
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
    
    etau_results, etau_pair, etau_trig_ids, etau_trig_types = self[etau_selection](events,
                                                                                   electron_indices,
                                                                                   tau_indices,
                                                                                   results,
                                                                                   call_force=True)
    results += etau_results

    has_one_etau_pair = ak.num(etau_pair.rawIdx, axis=1) == 2
    results += SelectionResult(steps={"one etau pair (sanity check)": has_one_etau_pair})

    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_etau_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    results += SelectionResult(steps={"fall into etau channel (sanity check)": events.channel_id == self.config_inst.get_channel(self.config_inst.x.channel).id})
    
    single_e_triggered = ak.any(etau_trig_types == 'single_e', axis=1)
    cross_e_triggered  = ak.any(etau_trig_types == 'cross_e_tau', axis=1)
    
    events = set_ak_column(events, "single_triggered", single_e_triggered)
    events = set_ak_column(events, "cross_triggered",  cross_e_triggered)

    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_electron_indices,
                                                                veto_muon_indices,
                                                                etau_pair)
    results += extra_lepton_veto_results

    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, etau_pair)
    results += zcand_results
    
    
    # Zcand prod results
    events, zcandprod_results = self[selzcandprod](events, zcand_array)
    results += zcandprod_results


    events, jet_clean_result, _ = self[jet_cleaning](events,
                                                     jet_indices,
                                                     results,
                                                     #ditaujet_jet_indices,
                                                     #call_force=True,
                                                     **kwargs)
    results += jet_clean_result


    # gen particles info
    if self.config_inst.x.extra_tags.genmatch:
        #if "is_signal" in list(self.dataset_inst.aux.keys()):
        print(" --->>> zcand-gentau matching")
        events, gentau_results = self[gentau_selection](events, True)
        results += gentau_results


    
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
                "Electron": etau_pair.rawIdx[:,0:1],
            },
            "Muon" : {
                "Muon": muon_indices,
            },            
            "Tau" : {
                "Tau": etau_pair.rawIdx[:,1:2],
            },
        },
    )
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

    
    return events, results
    

