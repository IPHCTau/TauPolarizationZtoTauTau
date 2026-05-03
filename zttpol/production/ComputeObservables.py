import os
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

from zttpol.production.PolarimetricA1 import PolarimetricA1
from zttpol.production.helper import getlistofobservables, getCombOMEGA
from zttpol.production.piHelper import piHelper
from zttpol.production.rhoHelper import rhoHelper
from zttpol.production.a1Helper import a1Helper




def get_observables_emu(ztollP4 : dict,
                        leg1    : str,
                        leg2    : str,
                        **kwargs):
    obs_vars = getlistofobservables()
    obs_temp = {var: None for var in obs_vars}

    if leg1 not in {'e','mu'}:
        raise RuntimeError("sanity check : leg1 must be e/mu")
    if leg2 not in {'e','mu'}:
        raise RuntimeError("sanity check : leg2 must be e/mu")

    massvis = (ztollP4['p4z1'] + ztollP4['p4z2']).mass
    obs_temp['massvis'] = massvis
    
    return obs_temp
    



def get_observables_etau(ztollP4 : dict,
                         leg1    : str,
                         leg2    : str,
                         **kwargs):

    regalgo = kwargs.get('regress_algo', '')

    obs_vars = getlistofobservables()
    obs_temp = {var: None for var in obs_vars}

    
    if leg1 != 'e':
        raise RuntimeError("sanity check : leg1 must be mu")

    if leg2 == 'pi':
        piObsCalc = piHelper(tau_p4 = ztollP4['p4z2'],
                             tau_pi_p4 = ztollP4['p4z2pi'],
                             debug = True)
        omegabar_2 = piObsCalc.getOmegaBar()
        obs_temp['omegabar_2'] = omegabar_2
        
        massvis = (ztollP4['p4z1'] + ztollP4['p4z2pi']).mass
        obs_temp['massvis'] = massvis
        
    
    elif leg2 == 'rho':
        rhoObsCalc = rhoHelper(tau_p4 = ztollP4['p4z2'],
                               tau_pi_p4 = ztollP4['p4z2pi'],
                               tau_pi0_p4 = ztollP4['p4z2pi0'],
                               debug = True)
        omegavis_2 = rhoObsCalc.getCosbeta()
        costheta_2 = rhoObsCalc.getCostheta()
        omega_2    = rhoObsCalc.getOmega(cosbeta = omegavis_2,
                                         costheta = costheta_2,
                                         cospsi = rhoObsCalc.getUltrarel_cospsiLF(costheta = costheta_2))

        # return costheta_2, cosbeta_2, omegavis_2
        obs_temp['costheta_2'] = costheta_2
        obs_temp['omegavis_2'] = omegavis_2
        obs_temp['omega_2']    = omega_2

        massvis = (ztollP4['p4z1'] + rhoObsCalc.LFproductLV).mass
        obs_temp['massvis'] = massvis

        omegabar_2 = rhoObsCalc.getOmegaBar()
        obs_temp['omegabar_2'] = omegabar_2

        
    elif leg2 == 'a1':
        a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z2'],
                             tau_pi_p4 = ztollP4['p4z2pi'],
                             debug = True)
        omegabar_2 = a1ObsCalc.getOmegaBar()
        obs_temp['omegabar_2'] = omegabar_2

        massvis = (ztollP4['p4z1'] + a1ObsCalc.LFa1LV).mass
        obs_temp['massvis'] = massvis
                
    
    else:
        raise RuntimeError(f"WRONG LEG for mutau : {leg2}")

    return obs_temp





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

        massvis = (ztollP4['p4z1'] + ztollP4['p4z2pi']).mass
        obs_temp['massvis'] = massvis
        
    
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

        massvis = (ztollP4['p4z1'] + rhoObsCalc.LFproductLV).mass
        obs_temp['massvis'] = massvis
        
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

        massvis = (ztollP4['p4z1'] + a1ObsCalc.LFa1LV).mass
        obs_temp['massvis'] = massvis
        
    
    else:
        raise RuntimeError(f"WRONG LEG for mutau : {leg2}")



    return obs_temp




def get_observables_tautau(ztollP4 : dict,
                           leg1    : str,
                           leg2    : str,
                           **kwargs):

    regalgo = kwargs.get('regress_algo', '')

    obs_vars = getlistofobservables()
    obs_temp = {var: None for var in obs_vars}


    VALID_LEGS = {"pi", "rho", "a1"}
    
    if leg1 not in VALID_LEGS:
        raise RuntimeError(f"sanity check : Invalid leg1 : {leg1}. Must be one of : {VALID_LEGS}")
    if leg2 not in VALID_LEGS:
        raise RuntimeError(f"sanity check : Invalid leg2 : {leg2}. Must be one of : {VALID_LEGS}")


    # ------------------------------------- #
    #        pi-pi / pi-rho / pi-a1         #
    # ------------------------------------- #
    if leg1 == "pi":
        # run svfit here

        # tau-1
        piObsCalc = piHelper(tau_p4 = ztollP4['p4z1'],
                             tau_pi_p4 = ztollP4['p4z1pi'],
                             debug = True)
        omegabar_1 = piObsCalc.getOmegaBar()        
        obs_temp['omegabar_1'] = omegabar_1

        if leg2 == "pi":
            # tau-2
            piObsCalc = piHelper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = piObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2
            massvis = (ztollP4['p4z1pi'] + ztollP4['p4z2pi']).mass
            obs_temp['massvis'] = massvis
            
        elif leg2 == "rho":
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
            massvis = (ztollP4['p4z1pi'] + rhoObsCalc.LFproductLV).mass
            obs_temp['massvis'] = massvis
            
            omegabar_2 = rhoObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2


        elif leg2 == "a1":
            #from IPython import embed; embed()
            a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = a1ObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2
            massvis = (ztollP4['p4z1pi'] + a1ObsCalc.LFa1LV).mass
            obs_temp['massvis'] = massvis

        else:
            raise RuntimeWarning(f"wrong leg2 : {leg2}")
            
    # ------------------------------------- #
    #        rho-pi / rho-rho / rho-a1      #
    # ------------------------------------- #
    elif leg1 == "rho":
        rhoObsCalc = rhoHelper(tau_p4 = ztollP4['p4z1'],
                               tau_pi_p4 = ztollP4['p4z1pi'],
                               tau_pi0_p4 = ztollP4['p4z1pi0'],
                               debug = True)
        #from IPython import embed; embed()
        omegavis_1 = rhoObsCalc.getCosbeta()
        costheta_1 = rhoObsCalc.getCostheta()
        omega_1    = rhoObsCalc.getOmega(cosbeta = omegavis_1,
                                         costheta = costheta_1,
                                         cospsi = rhoObsCalc.getUltrarel_cospsiLF(costheta = costheta_1))
        
        # return costheta_2, cosbeta_2, omegavis_2
        obs_temp['costheta_1'] = costheta_1
        obs_temp['omegavis_1'] = omegavis_1
        obs_temp['omega_1']    = omega_1

        rhoLV = rhoObsCalc.LFproductLV
        
        omegabar_1 = rhoObsCalc.getOmegaBar()
        obs_temp['omegabar_1'] = omegabar_1

        if leg2 == "pi":
            # tau-2
            piObsCalc = piHelper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = piObsCalc.getOmegaBar()            
            obs_temp['omegabar_2'] = omegabar_2
            
            massvis = (rhoLV + ztollP4['p4z2pi']).mass
            obs_temp['massvis'] = massvis
            
        elif leg2 == "rho":
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
            massvis = (rhoLV + rhoObsCalc.LFproductLV).mass
            obs_temp['massvis'] = massvis
            
            omegabar_2 = rhoObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2

        elif leg2 == "a1":
            a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = a1ObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2

            massvis = (rhoLV + a1ObsCalc.LFa1LV).mass
            obs_temp['massvis'] = massvis

        else:
            raise RuntimeWarning(f"wrong leg2 : {leg2}")

        
    # ------------------------------------- #
    #         a1-pi / a1-rho / a1-a1        #
    # ------------------------------------- #
    elif leg1 == "a1":
        a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z1'],
                             tau_pi_p4 = ztollP4['p4z1pi'],
                             debug = True)
        omegabar_1 = a1ObsCalc.getOmegaBar()
        obs_temp['omegabar_1'] = omegabar_1

        a1LV = a1ObsCalc.LFa1LV
        
        if leg2 == "pi":
            # tau-2
            piObsCalc = piHelper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = piObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2

            massvis = (a1LV + ztollP4['p4z2pi']).mass
            obs_temp['massvis'] = massvis
            
        elif leg2 == "rho":
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
            massvis = (a1LV + rhoObsCalc.LFproductLV).mass
            obs_temp['massvis'] = massvis
            
            omegabar_2 = rhoObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2

        elif leg2 == "a1":
            a1ObsCalc = a1Helper(tau_p4 = ztollP4['p4z2'],
                                 tau_pi_p4 = ztollP4['p4z2pi'],
                                 debug = True)
            omegabar_2 = a1ObsCalc.getOmegaBar()
            obs_temp['omegabar_2'] = omegabar_2

            massvis = (a1LV + a1ObsCalc.LFa1LV).mass
            obs_temp['massvis'] = massvis

        else:
            raise RuntimeWarning(f"wrong leg2 : {leg2}")

            
    else:
        raise RuntimeError(f"WRONG LEG combination for tautau : {leg1}-{leg2}")


    # construction of OMEGA
    omega_1 = obs_temp['omega_1']
    omega_2 = obs_temp['omega_2']
    if omega_1 is not None and omega_2 is not None:
        OMEGA = getCombOMEGA(omega_1, omega_2, debug=True)
        obs_temp['OMEGA'] = OMEGA

    omegavis_1 = obs_temp['omegavis_1']
    omegavis_2 = obs_temp['omegavis_2']
    if omegavis_1 is not None and omegavis_2 is not None:
        OMEGAVIS = getCombOMEGA(omegavis_1, omegavis_2, debug=True)
        obs_temp['OMEGAVIS'] = OMEGAVIS
        
    omegabar_1 = obs_temp['omegabar_1']
    omegabar_2 = obs_temp['omegabar_2']
    if omegabar_1 is not None and omegabar_2 is not None:
        OMEGABAR = getCombOMEGA(omegabar_1, omegabar_2, debug=True)
        obs_temp['OMEGABAR'] = OMEGABAR


        
    return obs_temp
