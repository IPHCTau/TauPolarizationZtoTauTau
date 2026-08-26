# coding: utf-8

"""
Column production methods related to generic event weights.
"""

from __future__ import annotations

import re
import copy
import functools

import law

from columnflow.production import Producer, producer
from columnflow.production.normalization import stitched_normalization_weights

from columnflow.production.cms.electron import electron_weights
from columnflow.production.cms.muon import muon_weights
from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.pdf import pdf_weights
from columnflow.production.cms.scale import murmuf_weights
from columnflow.production.cms.parton_shower import ps_weights

from columnflow.util import maybe_import, safe_div
from columnflow.columnar_util import set_ak_column
from columnflow.types import Any

from zttpol.production.tau import tau_weights

ak = maybe_import("awkward")
np = maybe_import("numpy")


# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

# ------------------------------------- #
#     Scale MC Weight : -1 or +1        #
# ------------------------------------- #
@producer(
    uses={"mc_weight"},
    produces={"mc_weight"},
    mc_only=True,
)
def scale_mc_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    mc_weight = np.sign(events.mc_weight)
    events = set_ak_column(events, "mc_weight", mc_weight, value_type=np.float32)

    return events



# ------------------------------------- #
#            Electron Weight            #
# ------------------------------------- #
electron_id_weights = electron_weights.derive(
    "electron_id_weights",
    cls_dict={
        "get_electron_file": (lambda self, external_files: external_files.electron_sf),
        "get_electron_config": (lambda self: self.config_inst.x.electron_id_sf),
        "use_supercluster_eta": True,
        "weight_name": "electron_id_weight",
    },
)
electron_reco_weights = electron_weights.derive(
    "electron_reco_weights",
    cls_dict={
        "get_electron_file": (lambda self, external_files: external_files.electron_sf),
        "get_electron_config": (lambda self: self.config_inst.x.electron_reco_sf),
        "use_supercluster_eta": True,
        "weight_name": "electron_reco_weight",
    },
)
# ------------------------------------- #
#              Muon Weight              #
# ------------------------------------- #
muon_id_weights = muon_weights.derive(
    "muon_id_weights",
    cls_dict={
        "get_muon_file": (lambda self, external_files: external_files.muon_sf),
        "get_muon_config": (lambda self: self.config_inst.x.muon_id_sf),
        "weight_name": "muon_id_weight",
    },
)
muon_iso_weights = muon_weights.derive(
    "muon_iso_weights",
    cls_dict={
        "get_muon_file": (lambda self, external_files: external_files.muon_sf),
        "get_muon_config": (lambda self: self.config_inst.x.muon_iso_sf),
        "weight_name": "muon_iso_weight",
    },
)
muon_IsoMu24_trigger_weights = muon_weights.derive(
    "muon_single_trigger_weights",
    cls_dict={
        "get_muon_file": (lambda self, external_files: external_files.muon_sf),
        "get_muon_config": (lambda self: self.config_inst.x.muon_single_trigger_sf),
        "weight_name": "muon_single_trigger_weight"
    }
)

# Use muon_IsoMu24_trigger_weights in the func below
# to use the pt and trigger_id condition
@producer(
    uses={
        "Muon.{pt,eta,phi,mass}",
        "single_triggered",
        muon_IsoMu24_trigger_weights,
    },
    produces={
        muon_IsoMu24_trigger_weights,
    },
)
def muon_single_trigger_weights(
    self: Producer,
    events: ak.Array,
    **kwargs,
) -> ak.Array:
    """
    Calculate the single-muon trigger SF weight.
    """

    # Trigger SF is only defined for pT >= 26 GeV.
    # Clamp pT only in the temporary array used for SF evaluation.
    trigger_sf_events = set_ak_column_f32(
        events,
        "Muon.pt",
        ak.where(events.Muon.pt >= 26.0, events.Muon.pt, 26.0),
    )

    trigger_sf_events = self[muon_IsoMu24_trigger_weights](
        trigger_sf_events,
        **kwargs,
    )

    for route in self[muon_IsoMu24_trigger_weights].produced_columns:
        events = set_ak_column_f32(
            events,
            route,
            ak.where(
                events.single_triggered,
                route.apply(trigger_sf_events),
                1.0,
            ),
        )

    del trigger_sf_events

    return events


@producer(
    uses={
        muon_single_trigger_weights,
    },
    produces={
        "trigger_weight",
        "trigger_weight_up",
        "trigger_weight_down",
    },
)
def muon_trigger_weights(
    self: Producer,
    events: ak.Array,
    **kwargs,
) -> ak.Array:

    events = self[muon_single_trigger_weights](events, **kwargs)

    for postfix in ["", "_up", "_down"]:
        events = set_ak_column_f32(
            events,
            f"trigger_weight{postfix}",
            events[f"muon_single_trigger_weight{postfix}"],
        )

    return events

"""


@producer(
    uses={
        "Muon.{pt,eta,phi,mass}",
        "single_triggered",
        muon_IsoMu24_trigger_weights,
    },
    produces={
        *[f"trigger_weight{tag}" for tag in ["", "_up", "_down"]],
    },
)
def muon_trigger_weights(self: Producer,
                         events: ak.Array,
                         **kwargs) -> ak.Array:
    # compute muon trigger SF weights (NOTE: trigger SFs are only defined for muons with
    # pt > 26 GeV, so create a copy of the events array with with all muon pt < 26 GeV set to 26 GeV)
    trigger_sf_events = set_ak_column_f32(events, "Muon.pt", ak.where(events.Muon.pt >= 26., events.Muon.pt, 26.))
    trigger_sf_events = self[muon_IsoMu24_trigger_weights](trigger_sf_events, **kwargs)
    for route in self[muon_IsoMu24_trigger_weights].produced_columns:
        events = set_ak_column_f32(events, route, ak.where(events.single_triggered,
                                                           route.apply(trigger_sf_events),
                                                           1.0))
    # memory cleanup
    del trigger_sf_events

    return events
"""

# ------------------------------------- #
#              TauID Weight             #
# ------------------------------------- #
tau_id_weights = tau_weights.derive(
    "tau_id_weights",
    cls_dict={
        "get_tau_file": (lambda self, external_files: external_files.tau_sf),
        "get_tau_config": (lambda self: self.config_inst.x.tau_id_sf),
        "weight_name": "tau_id_weight",
    },
)




# ------------------------------------- #
#       Normalized PileUp Weight        #
# ------------------------------------- #
@producer(
    uses={
        pu_weight.PRODUCES,
        "process_id"
    },
    mc_only=True,
)
def normalized_pu_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    for weight_name in self.pu_weight_names:
        # create a weight vector starting with ones
        norm_weight_per_pid = np.ones(len(events), dtype=np.float32)

        # fill weights with a new mask per unique process id (mostly just one)
        for pid in self.unique_process_ids:
            pid_mask = events.process_id == pid
            norm_weight_per_pid[pid_mask] = self.ratio_per_pid[weight_name][pid]

        # multiply with actual weight
        norm_weight_per_pid = norm_weight_per_pid * events[weight_name]

        # store it
        norm_weight_per_pid = ak.values_astype(norm_weight_per_pid, np.float32)
        events = set_ak_column_f32(events, f"normalized_{weight_name}", norm_weight_per_pid)

    return events


@normalized_pu_weight.post_init
def normalized_pu_weight_post_init(self: Producer, task: law.Task, **kwargs) -> None:
    super(normalized_pu_weight, self).post_init_func(task=task, **kwargs)

    # remember pu columns to read and produce
    self.pu_weight_names = {
        weight_name
        for weight_name in map(str, self[pu_weight].produced_columns)
        if (
            weight_name.startswith("pu_weight") and
            (task.global_shift_inst.is_nominal or not weight_name.endswith(("_up", "_down")))
        )
    }
    # adjust columns
    self.uses -= {pu_weight.PRODUCES}
    self.uses |= self.pu_weight_names
    self.produces |= {f"normalized_{weight_name}" for weight_name in self.pu_weight_names}


@normalized_pu_weight.requires
def normalized_pu_weight_requires(self: Producer, task: law.Task, reqs: dict, **kwargs) -> None:
    super(normalized_pu_weight, self).requires_func(task=task, reqs=reqs, **kwargs)

    from columnflow.tasks.selection import MergeSelectionStats
    reqs["selection_stats"] = MergeSelectionStats.req_different_branching(
        task,
        branch=-1 if task.is_workflow() else 0,
    )


@normalized_pu_weight.setup
def normalized_pu_weight_setup(self: Producer, task: law.Task, inputs: dict, **kwargs) -> None:
    super(normalized_pu_weight, self).setup_func(task=task, inputs=inputs, **kwargs)

    
    # load the selection stats
    #selection_stats = task.cached_value(
    #    key="selection_stats",
    #    func=lambda: inputs["selection_stats"]["collection"][0]["stats"].load(formatter="json"),
    #)
    selection_stats = task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    )

    # get the unique process ids in that dataset
    key = "sum_mc_weight_pu_weight_per_process"
    self.unique_process_ids = list(map(int, selection_stats[key].keys()))

    # helper to get numerators and denominators
    def numerator_per_pid(pid):
        key = "sum_mc_weight_per_process"
        return selection_stats[key].get(str(pid), 0.0)

    def denominator_per_pid(weight_name, pid):
        key = f"sum_mc_weight_{weight_name}_per_process"
        return selection_stats[key].get(str(pid), 0.0)

    # extract the ratio per weight and pid
    self.ratio_per_pid = {
        weight_name: {
            pid: safe_div(numerator_per_pid(pid), denominator_per_pid(weight_name, pid))
            for pid in self.unique_process_ids
        }
        for weight_name in (str(route) for route in self[pu_weight].produced_columns)
        if weight_name.startswith("pu_weight")
    }


# ------------------------------------- #
#         Normalized Pdf Weight         #
# ------------------------------------- #
@producer(
    uses={
        pdf_weights.PRODUCES
    },
    mc_only=True,
)
def normalized_pdf_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    for weight_name in self.pdf_weight_names:
        # create the normalized weight
        avg = self.average_pdf_weights[weight_name]
        normalized_weight = events[weight_name] / avg

        # store it
        events = set_ak_column_f32(events, f"normalized_{weight_name}", normalized_weight)

    return events


@normalized_pdf_weight.post_init
def normalized_pdf_weight_post_init(self: Producer, task: law.Task, **kwargs) -> None:
    super(normalized_pdf_weight, self).post_init_func(task=task, **kwargs)

    # remember pdf columns to read and produce
    self.pdf_weight_names = {
        weight_name
        for weight_name in map(str, self[pdf_weights].produced_columns)
        if (
            weight_name.startswith("pdf_weight") and
            (task.global_shift_inst.is_nominal or not weight_name.endswith(("_up", "_down")))
        )
    }
    # adjust columns
    self.uses.clear()
    self.uses |= self.pdf_weight_names
    self.produces |= {f"normalized_{weight_name}" for weight_name in self.pdf_weight_names}


@normalized_pdf_weight.requires
def normalized_pdf_weight_requires(self: Producer, task: law.Task, reqs: dict, **kwargs) -> None:
    super(normalized_pdf_weight, self).requires_func(task=task, reqs=reqs, **kwargs)

    from columnflow.tasks.selection import MergeSelectionStats
    reqs["selection_stats"] = MergeSelectionStats.req_different_branching(
        task,
        branch=-1 if task.is_workflow() else 0,
    )


@normalized_pdf_weight.setup
def normalized_pdf_weight_setup(self: Producer, task: law.Task, inputs: dict, **kwargs) -> None:
    super(normalized_pdf_weight, self).setup_func(task=task, inputs=inputs, **kwargs)

    # load the selection stats
    selection_stats = task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    )
    
    # save average weights
    self.average_pdf_weights = {
        weight_name: safe_div(selection_stats[f"sum_{weight_name}"], selection_stats["num_events"])
        for weight_name in self.pdf_weight_names
    }


# variation of the pdf weights producer that does not store up and down shifted weights
# but that stores all available pdf weights for the full treatment based on histograms
all_pdf_weights = pdf_weights.derive("all_pdf_weights", cls_dict={"store_all_weights": True})


# ------------------------------------- #
#    Normalized renorm and fact scale   #
# ------------------------------------- #
@producer(
    uses={murmuf_weights.PRODUCES},
    mc_only=True,
)
def normalized_murmuf_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    for weight_name in self.mu_weight_names:
        # create the normalized weight
        avg = self.average_mu_weights[weight_name]
        normalized_weight = events[weight_name] / avg

        # store it
        events = set_ak_column_f32(events, f"normalized_{weight_name}", normalized_weight)

    return events


@normalized_murmuf_weight.post_init
def normalized_murmuf_weight_post_init(self: Producer, task: law.Task, **kwargs) -> None:
    super(normalized_murmuf_weight, self).post_init_func(task=task, **kwargs)

    # remember mur/muf columns to read and produce
    self.mu_weight_names = {
        weight_name
        for weight_name in map(str, self[murmuf_weights].produced_columns)
        if (
            weight_name.startswith("murmuf_weight") and
            (task.global_shift_inst.is_nominal or not weight_name.endswith(("_up", "_down")))
        )
    }
    # adjust columns
    self.uses.clear()
    self.uses |= self.mu_weight_names
    self.produces |= {f"normalized_{weight_name}" for weight_name in self.mu_weight_names}


@normalized_murmuf_weight.requires
def normalized_murmuf_weight_requires(self: Producer, task: law.Task, reqs: dict, **kwargs) -> None:
    super(normalized_murmuf_weight, self).requires_func(task=task, reqs=reqs, **kwargs)

    from columnflow.tasks.selection import MergeSelectionStats
    reqs["selection_stats"] = MergeSelectionStats.req_different_branching(
        task,
        branch=-1 if task.is_workflow() else 0,
    )


@normalized_murmuf_weight.setup
def normalized_murmuf_weight_setup(self: Producer, task: law.Task, inputs: dict, **kwargs) -> None:
    super(normalized_murmuf_weight, self).setup_func(task=task, inputs=inputs, **kwargs)

    # load the selection stats
    selection_stats = task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    )
    
    # save average weights
    self.average_mu_weights = {
        weight_name: safe_div(selection_stats[f"sum_{weight_name}"], selection_stats["num_events"])
        for weight_name in self.mu_weight_names
    }



# ------------------------------------- #
#            isr and fsr weights        #
# ------------------------------------- #
@producer(
    uses={ps_weights.PRODUCES},
    mc_only=True,
)
def normalized_ps_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    for weight_name in self.ps_weight_names:
        # create the normalized weight
        avg = self.average_ps_weights[weight_name]
        normalized_weight = events[weight_name] / avg

        # store it
        events = set_ak_column_f32(events, f"normalized_{weight_name}", normalized_weight)

    return events


@normalized_ps_weight.post_init
def normalized_ps_weight_post_init(self: Producer, task: law.Task, **kwargs) -> None:
    super(normalized_ps_weight, self).post_init_func(task=task, **kwargs)

    # remember ps weight columns to read and produce
    self.ps_weight_names = {
        weight_name
        for weight_name in map(str, self[ps_weights].produced_columns)
        if (
            "weight" in weight_name and
            (task.global_shift_inst.is_nominal or not weight_name.endswith(("_up", "_down")))
        )
    }
    # adjust columns
    self.uses.clear()
    self.uses |= self.ps_weight_names
    self.produces |= {f"normalized_{weight_name}" for weight_name in self.ps_weight_names}


@normalized_ps_weight.requires
def normalized_ps_weight_requires(self: Producer, task: law.Task, reqs: dict, **kwargs) -> None:
    super(normalized_ps_weight, self).requires_func(task=task, reqs=reqs, **kwargs)

    from columnflow.tasks.selection import MergeSelectionStats
    reqs["selection_stats"] = MergeSelectionStats.req_different_branching(
        task,
        branch=-1 if task.is_workflow() else 0,
    )


@normalized_ps_weight.setup
def normalized_ps_weight_setup(self: Producer, task: law.Task, inputs: dict, **kwargs) -> None:
    super(normalized_ps_weight, self).setup_func(task=task, inputs=inputs, **kwargs)

    # load the selection stats
    selection_stats = task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    )

    # save average weights
    self.average_ps_weights = {
        weight_name: safe_div(selection_stats[f"sum_{weight_name}"], selection_stats["num_events"])
        for weight_name in self.ps_weight_names
    }
