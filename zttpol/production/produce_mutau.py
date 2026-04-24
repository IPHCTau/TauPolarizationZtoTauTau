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
from zttpol.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP
#from zttpol.production.weights import tauspinner_weight
#from zttpol.production.extra_weights import ff_weight, classify_events

from zttpol.production.muon_weights import muon_id_weights, muon_iso_weights
from zttpol.production.muon_weights import muon_trigger_weights, muon_xtrigger_weights
#from zttpol.production.electron_weights import electron_idiso_weights, electron_trigger_weights, electron_xtrigger_weights
from zttpol.production.tau_weights import tau_id_weights

from zttpol.production.sample_split import split_dy


#from zttpol.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi
from zttpol.util import IF_DATASET_HAS_LHE_WEIGHTS, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL, IF_DATASET_IS_TT
from zttpol.util import IF_RUN2, IF_RUN3, IF_ALLOW_STITCHING, IF_GENMATCH, IF_GENMATCH_ON_FOR_SIGNAL, transverse_mass

from zttpol.production.applyFastMTT import apply_fastMTT
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
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- tau -- #
        tau_id_weights,
        #IF_DATASET_IS_SIGNAL(tauspinner_weights),
        #ff_weight,
        #classify_events,
        reArrangeDecayProducts,
        #ProduceDetPhiCP,
        IF_GENMATCH(reArrangeGenDecayProducts),
        #IF_GENMATCH(ProduceGenPhiCP),
        apply_fastMTT,
    },
    produces={
        produce_base,
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- tau -- #
        tau_id_weights,
        #IF_DATASET_IS_SIGNAL(tauspinner_weights),
        #ff_weight,
        #classify_events,
        #ProduceDetPhiCP,
        #IF_GENMATCH(ProduceGenPhiCP),
        apply_fastMTT,
    },
)
def produce_mutau(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[produce_base](events, **kwargs)

    # ################## #
    #     Run FastMTT    #
    # ################## #
    logger.info(" >>>--- FastMTT-Wiktors --->>> [Not as fast as you think]")
    events = self[apply_fastMTT](events, run_fmtt=self.config_inst.x.enable_fastMTT)

    # ########################### #
    # -------- For PhiCP -------- #
    # ########################### #
    events, P4_dict = self[reArrangeDecayProducts](events)
    #events   = self[ProduceDetPhiCP](events, P4_dict)
    
    if self.config_inst.x.extra_tags.genmatch:
        events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
        #events = self[ProduceGenPhiCP](events, P4_gen_dict) 

    
    #logger.info(" >>>--- Evaluate Classifier Models (IC) --->>> [In extra_weights.py and processes.py]")
    #events = self[classify_events](events, **kwargs)
    
    if self.dataset_inst.is_mc:
        events = self[muon_id_weights](events, **kwargs)
        events = self[muon_iso_weights](events, **kwargs)
        events = self[muon_trigger_weights](events, **kwargs)
        events = self[muon_xtrigger_weights](events, **kwargs)
        events = self[tau_id_weights](events, do_syst=True, **kwargs)

        
    #events = self[ff_weight](events, **kwargs)        
    
    return events
