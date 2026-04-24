# coding: utf-8

"""
Base Selections common for all channels
"""

from __future__ import annotations

from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

import law
import order as od

from typing import Optional

from columnflow.selection import Selector, SelectionResult, selector

from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters
from columnflow.selection.cms.jets import jet_veto_map  

from columnflow.production.processes import process_ids
from columnflow.production.util import attach_coffea_behavior
##from columnflow.production.cms.top_pt_weight import gen_parton_top

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import IF_DATA, IF_MC

from zttpol.selection.physics_objects import (
    muon_selection, electron_selection, tau_selection, jet_selection, genZ_selection
)

from zttpol.selection.lepton_veto import tau_veto_from_dy
from zttpol.selection.trigger import trigger_selection
#from zttpol.selection.event_category import get_categories

from zttpol.production.columnvalid import make_column_valid

from zttpol.util import filter_by_triggers, get_objs_p4, trigger_object_matching_deep, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL


logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")




# exposed selector (can be invoked from the command line)
@selector(
    uses={
        #"event",
        # selectors / producers called within _this_ selector
        attach_coffea_behavior,
        json_filter, 
        met_filters, 
        process_ids,
        trigger_selection,
        ##IF_DATASET_IS_DY(genZ_selection),
        IF_MC(genZ_selection),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        ##jet_cleaning,
        #get_categories,
        #extra_lepton_veto, 
        ##double_lepton_veto,
        tau_veto_from_dy,
        jet_veto_map,
        "PuppiMET.{pt,phi}",
        "Jet.{pt,eta,phi,neEmEF,chEmEF}",
        "Flag.ecalBadCalibFilter",
        make_column_valid,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        trigger_selection,
        ##IF_DATASET_IS_DY(genZ_selection),
        IF_MC(genZ_selection),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        ##jet_cleaning,
        #get_categories, 
        process_ids,
        #extra_lepton_veto, 
        ##double_lepton_veto,
        make_column_valid,
    },
    exposed=False,
)
def select_base(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    if self.dataset_inst.has_tag("is_dy_m50"):
        events, tau_veto_result = self[tau_veto_from_dy](events, **kwargs)
        results += tau_veto_result
    
    results += SelectionResult(steps={"starts_with": np.ones(len(events), dtype=bool)})
        
    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    else:
        results += SelectionResult(steps={"json": np.ones(len(events), dtype=bool)})


    events, veto_result = self[jet_veto_map](events, **kwargs)
    results += veto_result

    # ############################ #
    #     met filter selection     #
    # ############################ #
    events, met_filter_results = self[met_filters](events, **kwargs)
    if self.dataset_inst.is_data:
        # BadCalibrationFilter is meant to reject events with noise which are related to bad crystals in Ecal.
        # The issue was present in some run ranges in data (end of 2022 & early 2023) only. The fraction of bad
        # data in terms of lumi is small. So, one can ignore the recipe in MC.
        # Also, DO NOT USE THIS FLAG AT ALL.
        #BadCalibrationFilter = ak.values_astype(events.Flag.ecalBadCalibFilter, bool)
        BadCalibrationFilter = events.event > 0 # all True
        met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")
        jet = ak.with_name(events.Jet, "PtEtaPhiMLorentzVector")
        BadCalibrationFilter_perjet_mask = (
            (events.PuppiMET.pt > 100)
            & (events.Jet.pt > 50)
            & (events.Jet.eta > -0.5) & (events.Jet.eta < -0.1)
            & (events.Jet.phi > -2.1) & (events.Jet.phi < -1.8)
            & ((events.Jet.neEmEF > 0.9) | (events.Jet.chEmEF > 0.9))
            & ak.all(met.metric_table(jet, metric=lambda a,b: a.delta_phi(b)) > 2.9, axis=1)
        )
        BadCalibrationFilter_mask = ak.values_astype(ak.any(BadCalibrationFilter_perjet_mask, axis=1), bool)
        # Keep the default BadCalibrationFilter if outside the run-range, else use the calibrationFilter_mask
        BadCalibrationFilter = ak.where((events.run >= 362433) & (events.run <= 367144), ~BadCalibrationFilter_mask, BadCalibrationFilter)
        
        BadCalibrationFilter_results = SelectionResult(steps={"met_filter:BadCalibration" : BadCalibrationFilter}) # negate the entire mask
        met_filter_results += BadCalibrationFilter_results
        
    results += met_filter_results    
    
    # trigger selection
    events, trigger_results = self[trigger_selection](events, channel = self.config_inst.x.channel)
    results += trigger_results
    
    # Get genZ collection for Zpt reweighting
    #if self.dataset_inst.has_tag("is_dy") or self.dataset_inst.has_tag("is_w") or self.dataset_inst.has_tag("is_signal"):
    if self.dataset_inst.is_mc:
        events = self[genZ_selection](events, **kwargs)
       
    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ] 
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           #call_force=False,
                                                                                                           **kwargs)
    results += ele_results

    
    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ] 
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           #call_force=False,
                                                                                                           **kwargs)
    results += muon_results


    # tau selection
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ]
    events, tau_results, good_tau_indices = self[tau_selection](events,
                                                                #call_force=True,
                                                                **kwargs)
    results += tau_results


    # jet selection
    events, jet_results, jet_indices = self[jet_selection](events,
                                                           #call_force=True,
                                                           **kwargs)
    results += jet_results

    # -------- Sel : b-veto -------- #
    # jet selection
    # -------------------------------------------- #
    # this is moved here, because now the jets are
    # cleaned against the tau cadidates of hacnd
    # -------------------------------------------- #
    #events, jet_clean_result, ditaujet_jet_indices = self[jet_cleaning](events,
    #                                                                    jet_indices,
    #                                                                    jet_results,
    #                                                                    ditaujet_jet_indices,
    #                                                                    call_force=True, 
    #                                                                    **kwargs)
    #results += jet_clean_result

    events = self[process_ids](events, **kwargs)

    events = self[make_column_valid](events)
    
    #print("Select_base done")
    
    return events, \
        results, \
        good_ele_indices, \
        veto_ele_indices, \
        good_muon_indices, \
        veto_muon_indices, \
        good_tau_indices, \
        jet_indices

