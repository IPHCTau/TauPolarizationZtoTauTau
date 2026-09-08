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
hist = maybe_import("hist")



def sort_pairs(
        dtrpairs: ak.Array,
        muon_iso = "pfRelIso04_all",
        tau_iso = "rawDeepTau2018v2p5VSjet"
) -> ak.Array:
    """
    Sorting mu-tau if multiple pairs exist
    Choose the most probable one in the end if trigeer objects are matched with
    Else, check the next one
    """
    sorted_idx = ak.argsort(dtrpairs["0"][muon_iso], ascending=True)

    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]

    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        ak.firsts(dtrpairs["0"][muon_iso][:,:1], axis=1) == ak.firsts(dtrpairs["0"][muon_iso][:,1:2], axis=1),
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
                          ak.argsort(dtrpairs["1"][tau_iso], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]

    # check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
    where_same_iso_2 = ak.fill_none(
        ak.firsts(dtrpairs["1"][tau_iso][:,:1], axis=1) == ak.firsts(dtrpairs["1"][tau_iso][:,1:2], axis=1),
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
        "Muon.{pt,eta,phi,mass,charge,pfRelIso04_all,rawIdx}",
        # tau
        "Tau.{pt,eta,phi,mass,charge,decayMode,decayModePNet,rawDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSjet,idDeepTau2018v2p5VSe,idDeepTau2018v2p5VSmu,rawIdx}",
        # met
        "PuppiMET.{pt,phi}",
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

    results = SelectionResult()
    
    # muons and taus e.g.
    # lep1: [ [m1], [m1],    [m1,m2], [],   [m1,m2] ]
    # lep2: [ [t1], [t1,t2], [t1],    [t1], [t1,t2] ]
    muons = events.Muon[lep1_indices] 
    taus  = events.Tau[lep2_indices]

    # Extra channel specific selections on m or tau
    vs_e_wp         = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_e["mutau"]
    vs_mu_wp        = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_m["mutau"]
    vs_jet_wp       = self.config_inst.x.tauIDWPs_config["DeepTau2018v2p5"].vs_j["mutau"]
    
    is_good_tau     = (
        (taus.pt > 20.0)
        & (taus.idDeepTau2018v2p5VSe   >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_e[vs_e_wp])
        & (taus.idDeepTau2018v2p5VSmu  >= self.config_inst.x.tauIDWPs["DeepTau2018v2p5"].vs_m[vs_mu_wp])
    )
    # idDeepTau2018v2p5VSjet will be used in the categorization

    taus = taus[is_good_tau]

    # Sorting muons by isolation [ascending]
    muons_sort_idxs = ak.argsort(muons.pfRelIso04_all, axis=-1, ascending=True)
    muons = muons[muons_sort_idxs]
    # sorting taus by deepTau score
    # WARNING : Don't forget to change it if any other tau tagger to be used
    taus_sort_idx = ak.argsort(taus.rawDeepTau2018v2p5VSjet, axis=-1, ascending=False)
    taus = taus[taus_sort_idx]
        
    leps_pair = ak.cartesian([muons, taus], axis=1)
    lep1, lep2 = ak.unzip(leps_pair)

    npair = ak.num(lep1, axis=1)
    results += SelectionResult(steps = {'mutau pairs pre pre-selection': npair >= 1})

    
    preselection = {
        #"mutau_is_os"         : (lep1.charge * lep2.charge) < 0,             # will be used in catgorization
        "mutau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,               # deltaR(lep1, lep2) > 0.5,
        "mutau_mT_60"         : transverse_mass(lep1, events.PuppiMET) < 60,  # will be used in catgorization
        "mutau_invmass_40"    : (1*lep1 + 1*lep2).mass > 40,                  # invariant_mass(lep1, lep2) > 40
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
    results += SelectionResult(steps = {'mutau pairs post pre-selection': npair >= 1})
    pair_selection_steps["mutau_before_trigger_matching"] = leps_pair["0"].pt >= 0.0
    
    # sort the pairs if many
    leps_pair = ak.where(npair > 1, sort_pairs(leps_pair), leps_pair)
    
    # match trigger objects for all pairs
    leps_pair, trigIds, trigTypes = match_trigobjs_semilep(leps_pair, trigger_results, channel_1='mu')

    npair = ak.num(leps_pair["0"], axis=1)
    results += SelectionResult(steps = {'mutau pairs post trigobj match': npair >= 1})
    pair_selection_steps["mutau_after_trigger_matching"] = leps_pair["0"].pt >= 0.0

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
        mutau_selection,
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
        mutau_selection,
        selzcand,
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
    
    mutau_results, mutau_pair, mutau_trig_ids, mutau_trig_types = self[mutau_selection](events,
                                                                                        muon_indices,
                                                                                        tau_indices,
                                                                                        results,
                                                                                        call_force=True)
    results += mutau_results

    has_one_mutau_pair = ak.num(mutau_pair.rawIdx, axis=1) == 2
    results += SelectionResult(steps={"one mutau pair (sanity check)": has_one_mutau_pair})
    
    # define channel ID
    channel_id = ak.values_astype(ak.where(has_one_mutau_pair, self.config_inst.get_channel(self.config_inst.x.channel).id, 0), np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    results += SelectionResult(steps={"fall into mutau channel (sanity check)": events.channel_id == self.config_inst.get_channel(self.config_inst.x.channel).id})
    
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

    # Zcand results
    events, zcand_array, zcand_results = self[selzcand](events, mutau_pair)
    results += zcand_results


    #from IPython import embed; embed()
    
    
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
                                                   #hists=hists,
                                                   no_sel=no_sel,
                                                   event_sel=event_sel)
    
    results += SelectionResult(
        objects = {
            "Electron" : {
                "Electron": electron_indices,
            },
            "Muon" : {
                "Muon": mutau_pair.rawIdx[:,0:1],
            },
            "Tau" : {
                "Tau": mutau_pair.rawIdx[:,1:2],
            },
        },
    )
    
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        debug_main(events,
                   results,
                   self.config_inst.x.triggers)

    
    return events, results

