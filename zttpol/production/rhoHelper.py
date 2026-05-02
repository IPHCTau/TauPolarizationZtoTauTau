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



class rhoHelper:
    def __init__(self,
                 tau_p4 = None,
                 tau_pi_p4 = None,
                 tau_pi0_p4 = None,
                 **kwargs):
        """
        Helper class to produce angular observables for rho in case of 
        tau to rho decay (DM = 1 : tau- -> rho- nu -> pi- pi0 nu)
        
        Arguments:
        1. ak.Array of LorentzVector of one tau decaying to rho 
           --> could be from generator level (true tau), regreseed by SVfit/FastMTT/anyMLbasedAlgo
        2. ak.Array of LorentzVector of charged prong (pion) decaying from tau (i.e. rho)
           --> reco is quite close to gen, so further regression is not needed 
        3. ak.Array of LorentzVector of pizero decaying from tau (i.e. rho)
           --> if possible to improve pi0, could be very useful [perhaps using ML]
        4. in **kwargs:
           debug = True or False (default) to print the frequency of unphysical solutions

        """

        self.LFtauLV       = to_cartesian(tau_p4)
        self.LFtauRhoPiLV  = to_cartesian(tau_pi_p4)
        self.LFtauRhoPi0LV = to_cartesian(tau_pi0_p4)

        self.LFproductLV   = self.getRefFrame()

        self.mpi    = 0.13957018 #GeV
        self.mpi0   = 0.1349766  #GeV
        self.mtau   = 1.77687    #GeV
        self.coscab = 0.975
        self.mrho   = 0.773      # GeV

        self.debug  = kwargs.get('debug', False)

        self.boostvec = self.LFproductLV.boostvec


        """
        self.InvisibleLV = self.TauLV - self.ProductLV
        self.DPF_InvisibleLV = self.InvisibleLV.Boost(self.boostvec.negative())

        # rotation
        RotVector = self.DPF_TauLV.pvec
        self.DPF_TauLV_rot = rotate(self.DPF_TauLV, RotVector)
        self.DPF_TauRhoPiLV_rot = rotate(self.DPF_TauRhoPiLV, RotVector)
        self.DPF_TauRhoPi0LV_rot = rotate(self.DPF_TauRhoPi0LV, RotVector)
        self.ProductLV_rot = rotate(self.ProductLV, RotVector)
        """
        
        
    def getRefFrame(self):
        """
          Sum of P4s : Pi + Pi0 
        """
        return self.LFtauRhoPiLV + self.LFtauRhoPi0LV


    def getOmegaBar(self):
        """
        omegaBar : Most optimized observable, but true tau P4 is very difficult to obtain
                   which is crucial to construct the polarimetric vector
        """

        # rotate lab-frame vectors using tau lab direction
        Rot1 = self.LFtauLV.pvec

        tauLabR1 = rotate(self.LFtauLV, Rot1)
        piLabR1  = rotate(self.LFtauRhoPiLV, Rot1)
        pi0LabR1 = rotate(self.LFtauRhoPi0LV, Rot1)

        # boost into tau rest frame
        #logger.warning('change boostvec extracted from rotated tau to lab tau (to be discussed)')
        #boostvec = self.TauLV.boostvec
        boostvec = tauLabR1.boostvec

        # Decay Product Frame
        DPFtauRotLV       = self.LFtauLV.boost(boostvec.negative())
        DPFtauRhoPiRotLV  = self.LFtauRhoPiLV.boost(boostvec.negative())
        DPFtauRhoPi0RotLV = self.LFtauRhoPi0LV.boost(boostvec.negative())

        
        hvec = getPolarimetricVector(tauP4    = DPFtauRotLV,
                                     pionP4   = DPFtauRhoPiRotLV,
                                     pizeroP4 = DPFtauRhoPi0RotLV,
                                     leg      = 'rho')
        #from IPython import embed; embed(); exit()
        hunit = hvec.unit
        #omegabar = hunit.dot(tauLabR1.pvec.unit) # basically cos(theta) type

        tauLabR1_pvec = tauLabR1.pvec
        tauLabR1_pvec_unit = tauLabR1_pvec/tauLabR1_pvec.mag
        omegabar = hunit.dot(tauLabR1_pvec_unit)
        
        near_boundary = (np.abs(omegabar) > 1) & (np.abs(omegabar) < 1.01)
        omegabar = ak.where(near_boundary & (omegabar > 0),
                            1.0,
                            ak.where(near_boundary & (omegabar < 0),
                                     -1.0,
                                     omegabar))
        
        bad_omegabar = np.abs(omegabar) > 1.0
        if (self.debug == True) and ak.any(bad_omegabar):
            logger.warning(f"no. of unphysical omegabar : {ak.sum(bad_omegabar)}")


        return omegabar
        
    

    def getOmega(self, cosbeta = None, costheta = None, cospsi = None):
        """
        omega : simplified version of omega from pure analytic approach (Not used ?)
        """
        QQ = self.LFproductLV.mass2

        sintheta = self.getSintheta(costheta = costheta)
        sinpsi   = self.getSinpsiLF(cospsi = cospsi)

        
        Be = 0.5 * (3 * cosbeta**2 - 1)
        Ps = 0.5 * (3 * cospsi - 1)

        RR = self.mtau**2 / QQ
        R  = np.sqrt(RR)
        
        num = (-2 + RR + 2 * (1 + RR) * Ps * Be) * costheta  +  3 * R * Be * sintheta * 2 * cospsi * sinpsi
        den = 2 + RR - 2 * (1 - RR) * Ps * Be

        valid = (QQ > 0) & (den != 0)

        if self.debug == True:
            logger.warning(f"no. of unphysical omega : {ak.sum(~valid)}")
            
        omega = ak.where(valid, num / den, -999.0)
        
        return omega

    
        
    def getCosbeta(self):
        """
        cos beta : Angle between the direction of charged pion 
                   and direction of rho in the rest frame of rho 
                   (AN2018_051_v8 : L700)
                   This is visible omega (used in 2016 analysis)
        """
        efrac   = (self.LFtauRhoPiLV.energy - self.LFtauRhoPi0LV.energy) / (self.LFtauRhoPiLV.energy + self.LFtauRhoPi0LV.energy)
        cosbeta = (self.mrho/np.sqrt(self.mrho * self.mrho - 4 * self.mpi * self.mpi))  *   efrac

        # make sure cosbeta should be within -1 and +1
        # rounding off values to -1 or +1 if tolerance < 0.01
        # else give warning
        near_boundary = (np.abs(cosbeta) > 1) & (np.abs(cosbeta) < 1.01)
        cosbeta = ak.where(near_boundary & (cosbeta > 0),
                           1.0,
                           ak.where(near_boundary & (cosbeta < 0),
                                    -1.0,
                                    cosbeta))
        
        bad_cosbeta = np.abs(cosbeta) > 1.0
        if (self.debug == True) and ak.any(bad_cosbeta):
            logger.warning(f"no. of unphysical cosbeta : {ak.sum(bad_cosbeta)}")

        return cosbeta
            
        
    def getUltrarel_cospsiLF(self, costheta = None):
        """
        Lorentz Transformation: maps costheta into cospsi
        """
        QQ = self.LFproductLV.mass2
        ct_rho = costheta
        num = ct_rho * (self.mtau**2 + QQ) + (self.mtau**2 - QQ)
        den = ct_rho * (self.mtau**2 - QQ) + (self.mtau**2 + QQ)
        
        cospsi = num / den

        near = (np.abs(cospsi) > 1) & (np.abs(cospsi) < 1.01)
        cospsi = ak.where(near & (cospsi > 0),
                          1.0,
                          ak.where(near & (cospsi < 0),
                                   -1.0,
                                   cospsi)
                          )
        bad_cospsi = np.abs(cospsi) > 1
        if (self.debug == True) and ak.any(bad_cospsi):
            logger.warning(f"no. of unphysical cospsi : {ak.sum(bad_cospsi)}")
            logger.warning("clipping cospsi within -1.0 and +1.0 (to be discussed)")            

        # clipping unphysical cospsi
        cospsi = ak.where(cospsi < -1.0,
                          -1.0,
                          ak.where(cospsi > 1.0,
                                   1.0,
                                   cospsi))            
        return cospsi
        

    def getCostheta(self):
        """
        cos theta : 
        """
        QQ = self.LFproductLV.mass2
        x  = self.LFproductLV.energy / self.LFtauLV.energy
        s  = 4 * self.LFtauLV.energy**2
        
        root_arg = 1 - 4 * self.mtau**2 / s
        isvalid = root_arg > 0

        #from IPython import embed; embed()
        
        costheta_raw = (2 * x * self.mtau**2 - self.mtau**2 - QQ) / ((self.mtau**2 - QQ) * np.sqrt(root_arg))
        costheta = ak.where(isvalid, costheta_raw, -999.0)

        near = (np.abs(costheta) > 1) & (np.abs(costheta) < 1.01)
        costheta = ak.where(near & (costheta > 0),
                            1.0,
                            ak.where(near & (costheta < 0),
                                     -1.0,
                                     costheta))
        bad_costheta = np.abs(costheta) > 1
        if (self.debug == True) and ak.any(bad_costheta):
            logger.warning(f"no. of unphysical costheta : {ak.sum(bad_costheta)}")
            logger.warning("clipping costheta within -1.0 and +1.0 (to be discussed)")

        # clipping unphysical costheta
        costheta = ak.where(costheta < -1.0,
                            -1.0,
                            ak.where(costheta > 1.0,
                                     1.0,
                                     costheta))
        
        return costheta
        
        
    
    def getSintheta(self, costheta = None):
        sintheta = np.sqrt(1 - costheta**2)
        return sintheta

        
    def getSinpsiLF(self, cospsi = None):
        sinpsi = np.sqrt(1 - cospsi**2)
        return sinpsi
