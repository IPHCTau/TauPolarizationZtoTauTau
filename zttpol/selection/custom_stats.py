# coding: utf-8

"""
Custom stats for HLepRare samples
"""

import law

from typing import Optional
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats

from columnflow.production.processes import process_ids
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column


logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")




@selector(uses={"process_id", optional("mc_weight")}) #TODO: Move it to utils
def custom_increment_stats(
    self: Selector,
    events: ak.Array,
    results: SelectionResult,
    stats: dict,
    **kwargs,
) -> ak.Array:
    """
    Unexposed selector that does not actually select objects but instead increments selection
    *stats* in-place based on all input *events* and the final selection *mask*.
    This function saves the nevents/nfiles (from campaign) per file, so that at 
    the time of merging, the total nevents produced can be used for normalization.
    instead of nevents, the sum_mc_weights will be used actually.
    ** the same nevents or sum_wt will be saved for each process id.
    e.g. In an inclusive dataset, if three different processes are there,
    all of those will have the same number. After merging, all processes
    eventually will have same numbers.
    """
    # get event masks
    event_mask = results.event

    # get a list of unique process ids present in the chunk
    unique_process_ids = np.unique(events.process_id)

    # increment plain counts
    n_evt_per_file = self.dataset_inst.aux['n_events']/self.dataset_inst.n_files # new
    sumwt_per_file = self.dataset_inst.n_events/self.dataset_inst.n_files     # new

    stats["num_events"] = int(n_evt_per_file)
    stats.setdefault(f"num_events_per_process", defaultdict(int))
    for p in unique_process_ids:
        # for splitting, each process will have the same number of events
        stats[f"num_events_per_process"][int(p)] = int(n_evt_per_file)
        
        
    if self.dataset_inst.is_mc:
        stats[f"sum_mc_weight"] = sumwt_per_file
        stats.setdefault(f"sum_mc_weight_per_process", defaultdict(float))
        for p in unique_process_ids:
            # for splitting, each process will have the same sumwt 
            stats[f"sum_mc_weight_per_process"][int(p)] = sumwt_per_file
        
    return events, results
