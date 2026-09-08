# coding: utf-8
"""                                                                                                                                                                                                               
Column production methods related to higher-level features.
"""
import functools

import law
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column, Route, EMPTY_FLOAT
from columnflow.util import maybe_import
from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

logger = law.logger.get_logger(__name__)


"""
@producer(
    uses={
        # nano columns
        "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
        "PuppiMET.significance",
    },
    produces={
        "PuppiMET.covXX", "PuppiMET.covXY", "PuppiMET.covYY",
        "PuppiMET.significance",
    },
)
def make_column_valid(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = set_ak_column(events, "PuppiMET.covXX", ak.nan_to_num(events.PuppiMET.covXX, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.covXY", ak.nan_to_num(events.PuppiMET.covXY, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.covYY", ak.nan_to_num(events.PuppiMET.covYY, nan=-9999.9))
    events = set_ak_column(events, "PuppiMET.significance", ak.nan_to_num(events.PuppiMET.significance, nan=-9999.9))
    
    return events
"""



@producer()
def make_column_valid(
        self: Producer,
        events: ak.Array,
        debug=False,
        **kwargs,
) -> ak.Array:

    for route in self.produced_columns:
        column = route.apply(events)

        try:
            finite = ak.fill_none(np.isfinite(column), True)
        except (TypeError, ValueError):
            continue

        n_bad = int(ak.sum(~finite, axis=None))
        
        if n_bad:
            if debug:
                logger.warning(
                    f"{route}: replacing {n_bad} non-finite values"
                )
            
            fixed = ak.where(
                finite,
                column,
                ak.full_like(column, EMPTY_FLOAT),
            )

            events = set_ak_column(
                events,
                route,
                fixed,
            )


    return events
    

@make_column_valid.init
def make_column_valid_init(self: Producer, **kwargs) -> None:

    self.validity_columns = {
        "PuppiMET.{covXX,covXY,covYY,significance}",
        "Tau.{IPx,IPy,IPz,ipLengthSig,rawUParTVSmu,rawUParTVSjet,rawUParTVSe,rawPNetVSmu,rawPNetVSjet,rawPNetVSe,qConfUParT,probDM2UParT,probDM2PNet,probDM1UParT,probDM1PNet,probDM11UParT,probDM11PNet,probDM10UParT,probDM10PNet,probDM0UParT,probDM0PNet}",
    }
        
        
    self.uses |= self.validity_columns
    self.produces |= self.validity_columns
