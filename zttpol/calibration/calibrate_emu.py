# coding: utf-8

"""
Exemplary calibration methods.
"""

from __future__ import annotations

import law
import functools

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.jets import jets
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from zttpol.calibration.electron import electron_smearing_scaling
from zttpol.calibration.tau import tau_energy_scale
from zttpol.calibration.calibrate_base import calibrate_base

from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)



# ----------------------------------------- #
#                emu Calibrator             #
# ----------------------------------------- #
@calibrator(
    uses={
        calibrate_base,
        electron_smearing_scaling,
    },
    produces={
        calibrate_base,
        electron_smearing_scaling,
    },
    exposed=True,
)
def calibrate_emu(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    """
      Calibrations used so far:
        - Base
        - Electron Scale & Smearing Correction (Smearing: MC only, Scale: DATA only)
    """

    # base calibrator
    events = events[calibrate_base](events)
    
    # electron scale and smearing correction
    events = self[electron_smearing_scaling](events)
        
    return events
