import law

import os
from typing import Optional
from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional


from zttpol.production.helper import to_cartesian



np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

logger = law.logger.get_logger(__name__)



class piHelper:
    def __init__(self,
                 tau_p4 = None,
                 tau_pi_p4 = None,
                 **kwargs):
        self.LFtauLV    = to_cartesian(tau_p4)
        self.LFtauPiLV  = to_cartesian(tau_pi_p4)

        refFrame = self.getRefFrame()
        self.boostvec = refFrame.pvec

        self.DPFtauLV   = self.LFtauLV.boost(self.boostvec.negative())
        self.DPFtauPiLV = self.LFtauPiLV.boost(self.boostvec.negative())
        

    def getOmegaBar(self):
        efrac = self.LFtauPiLV.energy/self.LFtauLV.energy
        omegabar = 2.0 * ak.where(efrac < 1.0, efrac, 1.0) - 1.0
        return omegabar

    def getRefFrame(self):
        return self.LFtauPiLV
