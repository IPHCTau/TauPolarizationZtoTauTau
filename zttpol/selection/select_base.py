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

from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import (
    EMPTY_FLOAT, Route, set_ak_column, IF_DATA, IF_MC
)
from columnflow.hist_util import create_hist_from_variables, fill_hist

from zttpol.selection.physics_objects import (
    muon_selection,
    electron_selection,
    tau_selection,
    jet_selection,
    genZ_selection
)
from zttpol.selection.lepton_veto import tau_veto_from_dy
from zttpol.selection.trigger import trigger_selection

from zttpol.production.columnvalid import make_column_valid

from zttpol.util import (
    IF_DATASET_IS_DY,
    IF_DATASET_IS_W,
    IF_DATASET_IS_SIGNAL,
    IF_RUN2, IF_RUN3,
    IF_DATASET_HAS_LHE_WEIGHTS
)
from columnflow.production.cms.jet import jet_id
from columnflow.production.cms.pileup import pu_weight
#from columnflow.production.cms.pdf import (
#    pdf_weights,
#    pdf_weights_raw
#)
from columnflow.production.cms.scale import (
    murmuf_weights,
    murmuf_weights_raw,
    murmuf_envelope_weights
)
from columnflow.production.cms.parton_shower import ps_weights

from zttpol.production.helper import jet_id_manual
#from zttpol.production.helper import assign_helicity
#from zttpol.production.sample_split_helicity import split_dy
from zttpol.production.weights import scale_mc_weight

from columnflow.types import TYPE_CHECKING

logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
hist = maybe_import("hist")



@selector(
    uses={
        attach_coffea_behavior,
        json_filter, 
        met_filters, 
        process_ids,
        trigger_selection,
        IF_MC(genZ_selection),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        tau_veto_from_dy,
        IF_RUN3(jet_veto_map),
        IF_RUN2(jet_id_manual),
        IF_RUN3(jet_id),
        "PuppiMET.{pt,phi}",
        "Jet.{pt,eta,phi,neEmEF,chEmEF}",
        "Flag.ecalBadCalibFilter",
        make_column_valid,
        IF_MC(pu_weight),
        IF_MC(ps_weights),
        #IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights, pdf_weights_raw,
        #                           murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),
        IF_DATASET_HAS_LHE_WEIGHTS(murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),
        #assign_helicity,
        #split_dy,
        IF_MC(scale_mc_weight),
    },
    produces={
        trigger_selection,
        IF_MC(genZ_selection),
        IF_RUN2(jet_id_manual),
        IF_RUN3(jet_id),
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        process_ids,
        make_column_valid,
        IF_MC(pu_weight),
        IF_MC(ps_weights),
        #IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights, pdf_weights_raw,
        #                           murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),
        IF_DATASET_HAS_LHE_WEIGHTS(murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),        
        #assign_helicity,
        #split_dy,
        IF_MC(scale_mc_weight),
    },
    exposed=False,
)
def select_base(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, ak.Array, SelectionResult]:

    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)

    # scale MC-weight : -1 or +1
    if self.dataset_inst.is_mc:
        events = self[scale_mc_weight](events, **kwargs)
    
    # no_sel
    no_sel_mask = np.ones(len(events), dtype=bool)
    
    # prepare the selection results recommended to be updated at every step
    results = SelectionResult()
    results += SelectionResult(steps={"starts_with": no_sel_mask})

    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    else:
        results += SelectionResult(steps={"json": np.ones(len(events), dtype=bool)})

    # if dy to 2 tau datasets are being used, veto the dy to 2 tau part from dy inclusive samples
    if self.dataset_inst.has_tag("is_dy_lep") and self.dataset_inst.has_tag("drop_tautau_from_dy_incl"):
        events, tau_veto_result = self[tau_veto_from_dy](events, **kwargs)
        results += tau_veto_result

    # set JetID
    if self.config_inst.campaign.x.run == 2:
        events = self[jet_id_manual](events, **kwargs)
    elif self.config_inst.campaign.x.run == 3:
        events = self[jet_id](events, **kwargs)
    else:
        raise RuntimeError(f'wrong run : {self.config_inst.campaign.x.run}')

    #from IPython import embed; embed()
    
    # jet veto map (only for Run 3)
    if self.has_dep(jet_veto_map):
        events, veto_result = self[jet_veto_map](events, **kwargs)
        results += veto_result

    # met filter
    events, met_filter_results = self[met_filters](events, **kwargs)
    # extra met filter for Run 3 (applicable for 2022 & 2023 only)
    if self.dataset_inst.is_data and (self.config_inst.campaign.x.run == 3):
        # BadCalibrationFilter is meant to reject events with noise which are related to bad crystals in Ecal.  #
        # The issue was present in some run ranges in data (end of 2022 & early 2023) only. The fraction of bad #
        # data in terms of lumi is small. So, one can ignore the recipe in MC.                                  #
        # Also, DO NOT USE THIS FLAG AT ALL.                                                                    #
        BadCalibrationFilter = events.event > 0
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
    events, trigger_results = self[trigger_selection](events)
    results += trigger_results
    
    # Get genZ collection for Zpt reweighting
    # todo: USE ONLY FOR PROCESSES INVOLVING W & Z
    #if self.dataset_inst.has_tag("is_dy") or self.dataset_inst.has_tag("is_w"):
    if self.dataset_inst.is_mc:
        events = self[genZ_selection](events, **kwargs)
       
    # electron selection
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           **kwargs)
    results += ele_results

    
    # muon selection
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           **kwargs)
    results += muon_results


    # tau selection
    events, tau_results, good_tau_indices = self[tau_selection](events,
                                                                **kwargs)
    results += tau_results

    # jet selection
    events, jet_results, jet_indices = self[jet_selection](events,
                                                           **kwargs)
    results += jet_results

    
    # build nLepton mask
    has_possible_pairs = events.event >= 0
    if self.config_inst.x.channel == "emu":
        has_possible_pairs = (ak.num(good_muon_indices, axis=1) > 0) & (ak.num(good_ele_indices, axis=1) > 0)
    elif self.config_inst.x.channel == "etau":
        has_possible_pairs = (ak.num(good_ele_indices, axis=1) > 0) & (ak.num(good_tau_indices, axis=1) > 0)
    elif self.config_inst.x.channel == "mutau":
        has_possible_pairs = (ak.num(good_muon_indices, axis=1) > 0) & (ak.num(good_tau_indices, axis=1) > 0)
    elif self.config_inst.x.channel == "tautau":
        has_possible_pairs = ak.num(good_tau_indices, axis=1) >= 2
    elif self.config_inst.x.channel == "ee":
        has_possible_pairs = ak.num(good_ele_indices, axis=1) >= 2
    elif self.config_inst.x.channel == "mumu":
        has_possible_pairs = ak.num(good_muon_indices, axis=1) >= 2  
    else:
        raise RuntimeError(f"Wrong Channel : {self.config_inst.x.channel}")

    results += SelectionResult(steps={"at least 2 leptons": has_possible_pairs})

    # generate processIDs
    events = self[process_ids](events, **kwargs)
        
    #if self.dataset_inst.has_tag("is_dy_tautau"):
    #    events = self[assign_helicity](events)
    #    events = self[split_dy](events,**kwargs)

    
    # take care of the NaN values in some coulmns
    events = self[make_column_valid](events, debug=self.config_inst.x.verbose.selection.main)

    # mc-only functions
    if self.dataset_inst.is_mc:

        # pdf weights
        #for pdf_cls in [pdf_weights, pdf_weights_raw]:
        #    if self.has_dep(pdf_cls):
        #        events = self[pdf_cls](
        #            events,
        #            outlier_log_mode="debug",
        #            invalid_weights_action="ignore" if self.dataset_inst.has_tag("partial_lhe_weights") else "raise",
        #            **kwargs,
        #        )

        # renormalization/factorization scale weights
        for murmuf_cls in [murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights]:
            if self.has_dep(murmuf_cls):
                events = self[murmuf_cls](events, **kwargs)

        # parton shower weights
        events = self[ps_weights](events, invalid_weights_action="ignore_one", **kwargs)
        
        # pileup weights
        events = self[pu_weight](events, **kwargs)


    
    # to be used in the channel specific selectors
    return events, \
        no_sel_mask, \
        results, \
        good_ele_indices, \
        veto_ele_indices, \
        good_muon_indices, \
        veto_muon_indices, \
        good_tau_indices, \
        jet_indices

