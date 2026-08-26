# coding: utf-8

from __future__ import annotations

import law

from typing import Optional
from operator import and_
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats

from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import (
    EMPTY_FLOAT,
    Route,
    set_ak_column,
    IF_DATA,
    IF_MC
)
from columnflow.hist_util import create_hist_from_variables, fill_hist

from zttpol.util import (
    IF_DATASET_IS_DY,
    IF_DATASET_IS_W,
    IF_DATASET_IS_SIGNAL,
    IF_RUN2,
    IF_RUN3,
    IF_DATASET_HAS_LHE_WEIGHTS
)
from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.pdf import (
    pdf_weights,
    pdf_weights_raw
)
from columnflow.production.cms.scale import (
    murmuf_weights,
    murmuf_weights_raw,
    murmuf_envelope_weights
)
from columnflow.production.cms.parton_shower import ps_weights

logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
hist = maybe_import("hist")

@selector(
    uses={
        IF_MC(pu_weight),
        IF_MC(ps_weights),
        #IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights, pdf_weights_raw,
        #                           murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),
        IF_DATASET_HAS_LHE_WEIGHTS(murmuf_weights, murmuf_weights_raw, murmuf_envelope_weights),
        "process_id",
        increment_stats,
    },
    exposed=False,
)
def custom_increment_stats(
        self: Selector,
        events: ak.Array,
        task: law.Task,
        results: SelectionResult,
        stats: defaultdict,
        no_sel: np.ndarray | ak.Array,
        event_sel: np.ndarray | ak.Array,
        event_sel_variations: dict[str, np.ndarray | ak.Array] | None = None, 
        **kwargs
) -> tuple[ak.Array, SelectionResult]:

    if event_sel_variations is None:
        event_sel_variations = {}
    event_sel_variations = {n: s for n, s in event_sel_variations.items() if s is not None}

    # when a shift was requested, skip all other systematic variations
    skip_shifts = task.global_shift_inst != "nominal"

    # start creating a "stats map"
    # - keys: names of histograms to be created
    # - values: (weight array, selection array)
    # note that only a subset of entries end up in the stats dictionary, but all are used for histograms
    stats_map: dict[str, np.ndarray | ak.Array | tuple[np.ndarray | ak.Array, np.ndarray | ak.Array]] = {}
    keys_for_stats = []

    # helper to cast to float64
    f64 = lambda a: ak.values_astype(a, np.float64)

    def add(key, sel, weight=None, for_stats=False, for_hists=True):
        stats_map[key] = sel if weight is None else (weight, sel)
        if for_stats and key not in keys_for_stats:
            keys_for_stats.append(key)

    # basic event counts
    add("num_events", no_sel, for_stats=True)
    add("num_events_selected", event_sel, for_stats=True)
    for var_name, var_sel in event_sel_variations.items():
        add(f"num_events_selected_{var_name}", var_sel, for_stats=True)

    # add mc info
    if self.dataset_inst.is_mc:
        add("sum_mc_weight", no_sel, events.mc_weight, for_stats=True)
        add("sum_mc_weight_selected", event_sel, events.mc_weight, for_stats=True)
        for var_name, var_sel in event_sel_variations.items():
            add(f"sum_mc_weight_selected_{var_name}", var_sel, events.mc_weight, for_stats=True)

        # pu weights with variations
        for route in sorted(self[pu_weight].produced_columns):
            add(f"sum_mc_weight_{route}", no_sel, events.mc_weight * route.apply(events))

        # pdf weights with variations
        #if self.has_dep(pdf_weights):
        #    for v in (("",) if skip_shifts else ("", "_up", "_down")):
        #        add(f"sum_pdf_weight{v}", no_sel, events[f"pdf_weight{v}"])
        #        add(f"sum_pdf_weight{v}_selected", event_sel, events[f"pdf_weight{v}"])

        # mur/muf weights with variations
        if self.has_dep(murmuf_weights):
            for v in (("",) if skip_shifts else ("", "_up", "_down")):
                add(f"sum_murmuf_weight{v}", no_sel, events[f"murmuf_weight{v}"])
                add(f"sum_murmuf_weight{v}_selected", event_sel, events[f"murmuf_weight{v}"])

        # parton shower weights with variations
        if self.has_dep(ps_weights):
            for v in (("",) if skip_shifts else ("", "_up", "_down")):
                add(f"sum_isr_weight{v}", no_sel, events[f"isr_weight{v}"])
                add(f"sum_isr_weight{v}_selected", event_sel, events[f"isr_weight{v}"])
                add(f"sum_fsr_weight{v}", no_sel, events[f"fsr_weight{v}"])
                add(f"sum_fsr_weight{v}_selected", event_sel, events[f"fsr_weight{v}"])

    group_map = {
        "process": {
            "values": events.process_id,
            "mask_fn": (lambda v: events.process_id == v),
        },
    }
        
    events, results = self[increment_stats](events=events,
                                            results=results,
                                            stats=stats,
                                            weight_map=stats_map,
                                            group_map=group_map,
                                            **kwargs)
                
    return events, results
