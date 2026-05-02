import os
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

from zttpol.production.PolarimetricA1 import PolarimetricA1
from zttpol.production.helper import getlistofobservables, get_OMEGA, get_OMEGAVIS, get_OMEGABAR
from zttpol.production.piHelper import piHelper
from zttpol.production.rhoHelper import rhoHelper
from zttpol.production.a1Helper import a1Helper




def get_observables_emu(ztollP4 : dict,
                        leg1    : str,
                        leg2    : str,
                        **kwargs):
    pass



def get_observables_etau(ztollP4 : dict,
                         leg1    : str,
                         leg2    : str,
                         **kwargs):
    pass




# Returning dictionary contains all of the above
def get_observables_mutau(ztollP4 : dict,
                          leg1    : str,
                          leg2    : str,
                          **kwargs):
    regalgo = kwargs.get('regress_algo', '')

    obs_vars = getlistofobservables()
    obs_temp = {var: None for var in obs_vars}

    
    if leg1 != 'mu':
        raise RuntimeError("sanity check : leg1 must be mu")

    if leg2 == 'pi':
        # run svfit here
        piObsCalc = piHelper(tau_p4 = ztollP4['p4z2'],
                             tau_pi_p4 = ztollP4['p4z2pi'],
                             debug = True)
        omegabar_2 = piObsCalc.getOmegaBar()

        obs_temp['omegabar_2'] = omegabar_2
        
    
    elif leg2 == 'rho':
        rhoObsCalc = rhoHelper(tau_p4 = ztollP4['p4z2'],
                               tau_pi_p4 = ztollP4['p4z2pi'],
                               tau_pi0_p4 = ztollP4['p4z2pi0'],
                               debug = True)
        #from IPython import embed; embed()
        omegavis_2 = rhoObsCalc.getCosbeta()
        costheta_2 = rhoObsCalc.getCostheta()
        omega_2    = rhoObsCalc.getOmega(cosbeta = omegavis_2,
                                         costheta = costheta_2,
                                         cospsi = rhoObsCalc.getUltrarel_cospsiLF(costheta = costheta_2))

        # return costheta_2, cosbeta_2, omegavis_2
        obs_temp['costheta_2'] = costheta_2
        obs_temp['omegavis_2'] = omegavis_2
        obs_temp['omega_2']    = omega_2

        # SVfit can be used here
        #  remove the empty lists, but keep the sequence in mind
        #  Run SVfit producer
        #  It should return regressed values
        #  Wrap the arrays with the old sequence
        #  Call rhoHelper
        #  Get any observable

        #print('FastMTT')
        #events, pt_reg, _ = self[apply_fastMTT](events)
        
        omegabar_2 = rhoObsCalc.getOmegaBar()
        obs_temp['omegabar_2'] = omegabar_2

        
    elif leg2 == 'a1':
        a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z2'],
                             tau_pi_p4 = ztollP4['p4z2pi'],
                             debug = True)
        omegabar_2 = a1ObsCalc.getOmegaBar()
        obs_temp['omegabar_2'] = omegabar_2
        
    
    else:
        raise RuntimeError(f"WRONG LEG for mutau : {leg2}")



    return obs_temp




def get_observables_tautau(ztollP4 : dict,
                           leg1    : str,
                           leg2    : str,
                           **kwargs):
    pass
