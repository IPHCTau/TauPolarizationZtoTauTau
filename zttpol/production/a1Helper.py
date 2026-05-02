import law

import os
from typing import Optional
from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional


from zttpol.production.helper import to_cartesian, rotate, getPolarimetricVector



np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

logger = law.logger.get_logger(__name__)



class a1Helper:
    def __init__(self,
                 tau_p4 = None,
                 tau_pi_p4 = None,
                 **kwargs):
        """
        """
        
        self.LFtauLV = to_cartesian(tau_p4)
        self.LFtauA1PiLV = to_cartesian(tau_pi_p4)
        self.LFosPionLV  = self.LFtauA1PiLV[:,0:1]
        self.LFss1pionLV = self.LFtauA1PiLV[:,1:2]
        self.LFss2pionLV = self.LFtauA1PiLV[:,2:3]
        self.LFa1LV = self.LFosPionLV + self.LFss1pionLV + self.LFss2pionLV

        self.RotVector = self.getRotationVec()
        self.RefFrameLV = self.getRefFrame()

        RefFrameRotated = rotate(self.RefFrameLV, self.RotVector)
        

        self.tau_rot = rotate(self.LFtauLV, self.RotVector)
        self.os_rot  = rotate(self.LFosPionLV, self.RotVector)
        self.ss1_rot = rotate(self.LFss1pionLV, self.RotVector)
        self.ss2_rot = rotate(self.LFss2pionLV, self.RotVector)

        self._tauAlongZLabFrame = self.tau_rot
        
        self.boostvec = RefFrameRotated.boostvec

        self.RFtauLV = self.tau_rot.boost(self.boostvec.negative())
        self.RFosPionLV = self.os_rot.boost(self.boostvec.negative())
        self.RFss1PionLV = self.ss1_rot.boost(self.boostvec.negative())
        self.RFss2PionLV = self.ss2_rot.boost(self.boostvec.negative())
        
        self.debug  = kwargs.get('debug', False)

        

    def getRefFrame(self):
        """
        everythings to be boosted in the tau rest frame
        """
        return self.LFtauLV


    def getRotationVec(self):
        return self.LFtauLV.pvec
                                            

    def getOmegaBar(self):
        nMinusTauAlongZLabFrame = -self.nTZLFr()
        pions = ak.concatenate([self.RFosPionLV,
                                self.RFss1PionLV,
                                self.RFss2PionLV], axis=1)
        h = getPolarimetricVector(tauP4  = self.RFtauLV, 
                                  pionP4 = pions,
                                  tauch  = self.LFtauLV.charge,
                                  leg    = 'a1')
        hunit = h.unit
        omegabar = nMinusTauAlongZLabFrame.dot(hunit)

        return omegabar
        
    
    def nTZLFr(self):
        vec = self._tauAlongZLabFrame.pvec
        mag = vec.mag
        #return self._tauAlongZLabFrame.pvec.unit
        return vec/mag
    
    
    def getOmega(self):
        pass
