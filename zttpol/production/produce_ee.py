# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer

from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

from columnflow.config_util import get_events_from_categories
from zttpol.production.ReArrangeZcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from zttpol.production.ProduceObservables import ProduceRecoObservables, ProduceGenObservables

from zttpol.production.weights import (
    electron_id_weights,
    electron_reco_weights,
    #electron_single_trigger_weights,
)

from zttpol.util import (
    IF_DATASET_HAS_LHE_WEIGHTS,
    IF_DATASET_IS_DY,
    IF_DATASET_IS_W,
    IF_DATASET_IS_SIGNAL,
    IF_DATASET_IS_TT
)
from zttpol.util import (
    IF_RUN2,
    IF_RUN3,
    IF_ALLOW_STITCHING,
    IF_GENMATCH,
    IF_GENMATCH_ON_FOR_SIGNAL,
    transverse_mass
)

from zttpol.production.produce_base import produce_base


np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)



@producer(
    uses={
        produce_base,
        # -- electron -- #
        electron_id_weights,
        electron_reco_weights,
        #electron_single_trigger_weights,
    },
    produces={
        produce_base,
        # -- electron -- #
        electron_id_weights,
        electron_reco_weights,
        #electron_single_trigger_weights,
    },
)
def produce_ee(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[produce_base](events, **kwargs)

    
    if self.dataset_inst.is_mc:
        events = self[electron_id_weights](events, **kwargs)
        events = self[electron_reco_weights](events, **kwargs)
        #events = self[muon_single_trigger_weights](events, **kwargs)

    #events = self[ff_weight](events, **kwargs)        
    
    return events
