# coding: utf-8

"""
Exemplary calibration methods.
"""

from __future__ import annotations

import law
import functools

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.tau import tec

from columnflow.util import maybe_import

from columnflow.columnar_util import set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import IF_DATA, IF_MC

from columnflow.production.util import attach_coffea_behavior

from zttpol.calibration.calibrate_base import calibrate_base

from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)



# ----------------------------------------- #
#             tautau Calibrator             #
# ----------------------------------------- #
@calibrator(
    uses={
        calibrate_base,
    },
    produces={
        calibrate_base,
    },
    exposed=True,
)
def calibrate_tautau(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    """
      Calibrations used so far:
        - Base
        - Tau Energy Scale Correction (MC only)
    """

    # base calibrator
    events = self[calibrate_base](events)

    # data/mc specific calibrations
    if self.dataset_inst.is_mc:
        # full tec
        events = self[self.tec_full_cls](events, **kwargs)

    return events



@calibrate_tautau.init
def calibrate_tautau_init(self: Calibrator, **kwargs) -> None:
    # set the name of the met collection to use
    met_name = self.config_inst.x.met_name
    raw_met_name = self.config_inst.x.raw_met_name

    # derive calibrators to add settings once
    flag = f"custom_calibs_registered_{self.cls_name}"
    if not self.config_inst.x(flag, False):
        def add_calib_cls(name, base, cls_dict=None):
            self.config_inst.set_aux(f"calib_{name}_cls", base.derive(name, cls_dict=cls_dict or {}))
        # derive tec calibrators
        add_calib_cls("tec_full", tec, cls_dict={
            "met_name": met_name,
            "propagate_met": False,  # not needed after JET-to-MET propagation
        })

        # change the flag
        self.config_inst.set_aux(flag, True)


    # store references to classes
    self.tec_full_cls = self.config_inst.x.calib_tec_full_cls
    
    # collect derived calibrators and add them to the calibrator uses and produces
    derived_calibrators = {
        self.tec_full_cls,
    }

    self.uses |= derived_calibrators
    self.produces |= derived_calibrators        
