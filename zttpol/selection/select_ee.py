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
from zttpol.selection.match_trigobj import match_trigobjs_dilep
from zttpol.selection.zcand import selzcand, selzcandprod
from zttpol.selection.debug import debug_main
from zttpol.selection.stats import custom_increment_stats

#from zttpol.production.columnvalid import make_column_valid

from zttpol.util import (
    transverse_mass,
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



def sort_pairs(
        dtrpairs: ak.Array,
        ele_iso = "pfRelIso03_all",
) -> ak.Array:
    sorted_idx = ak.argsort(dtrpairs["0"][ele_iso], ascending=True)

    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"][ele_iso][:,:1], axis=1) == ak.firsts(dtrpairs["0"][ele_iso][:,1:2], axis=1),
        False
    )
    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    return dtrpairs




@selector(
    uses={
        # Electron
        "Electron.{pt,eta,phi,mass,charge,pfRelIso03_all,rawIdx}",        
        # met
        "PuppiMET.{pt,phi}",
    },
    exposed=False,
)
def ee_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        tau_indices: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:

    results = SelectionResult()
    
    taus            = events.Tau[tau_indices]

    # Extra channel specific selections on m or tau
    vs_e_wp         = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_e["ee"]
    vs_mu_wp        = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_m["ee"]
    vs_jet_wp       = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_j["ee"]
    
    is_good_tau     = (
        (taus.pt > 20.0)
        & (taus.idDeepTau2018v2p5VSjet >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_j[vs_jet_wp])
        & (taus.idDeepTau2018v2p5VSe   >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_m[vs_mu_wp])
    )
    taus = taus[is_good_tau] # to be used to veto events
    results += SelectionResult(steps = {'no tauh': ak.num(taus, axis=1) == 0})
    
    eles = events.Electron[lep_indices] 

    # Sorting muons by isolation [ascending]
    eles_sort_idxs = ak.argsort(eles.pfRelIso03_all, axis=-1, ascending=True)
    eles = eles[eles_sort_idxs]
        
    leps_pair = ak.combinations(eles, 2, axis=1)
    lep1, lep2 = ak.unzip(leps_pair)

    npair = ak.num(lep1, axis=1)
    results += SelectionResult(steps = {'ee pairs pre pre-selection': npair >= 1})

    invm = (1*lep1 + 1*lep2).mass
    preselection = {
        #"ee_is_os"          : (lep1.charge * lep2.charge) < 0,              # will be used in catgorization
        "ee_dr_0p5"          : (1*lep1).delta_r(1*lep2) > 0.5,               # deltaR(lep1, lep2) > 0.5,
        "e1met_mT_60"        : transverse_mass(lep1, events.PuppiMET) < 60,  # 
        "e2met_mT_60"        : transverse_mass(lep2, events.PuppiMET) < 60,  # 
        "ee_invm_40_140"     : (invm > 40.0) & (invm < 140.0),               # 40 < invm < 140
    }

    good_pair_mask = lep1.rawIdx >= 0
    pair_selection_steps = {}
    pair_selection_steps["ee_starts_with"] = good_pair_mask
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask
        
    leps_pair = leps_pair[good_pair_mask]

    # check nPairs
    npair = ak.num(leps_pair["0"], axis=1)
    results += SelectionResult(steps = {'ee pairs post pre-selection': npair >= 1})
    pair_selection_steps["ee_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)
    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_dilep(leps_pair, trigger_results, channel_1='e')

    npair = ak.num(leps_pair["0"], axis=1)
    results += SelectionResult(steps = {'ee pairs post trigobj match': npair >= 1})
    pair_selection_steps["ee_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

    lep1, lep2 = ak.unzip(leps_pair)


    # take the 1st pair
    lep1 = lep1[:,:1]
    lep2 = lep2[:,:1]
    
    # rebuild the pair with the 1st one only
    leps_pair = ak.concatenate([lep1, lep2], axis=1)
    sort_idx = ak.argsort(leps_pair.pt, ascending=False)
    leps_pair = leps_pair[sort_idx]
    
    results += SelectionResult(aux=pair_selection_steps)
    
    return results, leps_pair, trigIds, trigTypes



@selector(
    uses={
        select_base,
        ee_selection,
        extra_lepton_veto,
        selzcand,
        custom_increment_stats,
        jet_cleaning,
        IF_MC("mc_weight"),
        #make_column_valid,
    },
    produces={
        select_base,
        ee_selection,
        selzcand,
        "channel_id",
        "single_triggered",
        #make_column_valid,
    },
    exposed=True,
)
def select_ee(
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
    
    ee_results, ee_pair, ee_trig_ids, ee_trig_types = self[ee_selection](events,
                                                                         electron_indices,
                                                                         tau_indices,
                                                                         results,
                                                                         call_force=True)
    results += ee_results

    has_one_ee_pair = ak.num(ee_pair.rawIdx, axis=1) == 2
    results += SelectionResult(steps={"one ee pair (sanity check)": has_one_ee_pair})
    
    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_ee_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    results += SelectionResult(steps={"fall into ee channel (sanity check)": events.channel_id == self.config_inst.get_channel(self.config_inst.x.channel).id})
    
    # muon
    single_e_triggered = ak.any(ee_trig_types == 'single_e', axis=1)
    
    events = set_ak_column(events, "single_triggered", single_e_triggered)

    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_electron_indices,
                                                                veto_muon_indices,
                                                                ee_pair)
    results += extra_lepton_veto_results

    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, ee_pair)
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
                "Electron": ee_pair.rawIdx,
            },
            "Muon" : {
                "Muon": muon_indices,
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

