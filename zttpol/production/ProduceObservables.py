# coding: utf-8
"""
A producer to create Acoplanarity angles for differennt methods
"""

import os
import copy
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from zttpol.production.PhiCP_Estimator import GetPhiCP
#from httcp.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional

from zttpol.production.helper import getlistofobservables, wrap, clean_rearranged_dict, unwrap
from zttpol.production.ComputeObservables import get_observables_mutau


import law
np = maybe_import("numpy")
ak = maybe_import("awkward")

logger = law.logger.get_logger(__name__)




@producer(
    uses={
        "channel_id",
        "zcand.{pt,eta,phi,mass,decayMode,charge,IPx,IPy,IPz}",
        optional("zcand.pt_fastMTT"), optional("zcand.eta_fastMTT"), optional("zcand.phi_fastMTT"), optional("zcand.mass_fastMTT"),
        "zcandprod.{pt,eta,phi,mass,pdgId}",
        optional("GenTau.pt"), optional("GenTau.eta"), optional("GenTau.phi"), optional("GenTau.mass"),
        optional("GenTau.IPx"), optional("GenTau.IPy"), optional("GenTau.IPz"),
        optional("GenTau.decayMode"), optional("GenTau.charge"), 
        optional("GenTauProd.pt"), optional("GenTauProd.eta"), optional("GenTauProd.phi"), optional("GenTauProd.mass"),
    },
)
def ProduceObservables(
        self: Producer,
        events: ak.Array,
        p4zcandinfo: dict,
        **kwargs,
) -> tuple[ak.Array,ak.Array,ak.Array,ak.Array,ak.Array]:

    channel = self.config_inst.x.channel
    
    PrepareP4 = lambda p4dict, mask : {key: ak.where(mask, val, val[:,:0]) for key, val in p4dict.items()}

    def extract_observables(observables, p4zcand, mask, leg1, leg2, observable_func, keylist=['p4z1']):
        zcandP4 = PrepareP4(p4zcand, mask)
        zcandP4, count = clean_rearranged_dict(zcandP4, keylist=keylist)
        observable_dict = observable_func(zcandP4, leg1, leg2) # get_observables_mutau
        observable_dict = unwrap(observable_dict, count)
        observables = wrap(mask, observable_dict, observables)
        return observables



    
    is_e      = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -1, axis=1), False)
    is_mu     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -2, axis=1), False)
    is_pi     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  0, axis=1), False)
    is_rho    = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  1, axis=1), False)
    is_a1DM2  = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  2, axis=1), False)
    is_a1DM10 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 10, axis=1), False)
    is_a1DM11 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 11, axis=1), False)
    is_a1     = is_a1DM10 # for now

    p4h1      = p4zcandinfo["p4z1"]
    p4h2      = p4zcandinfo["p4z2"]


    dummyPhiCP = ak.from_regular(ak.values_astype(events.event[:,None][:,:0], np.float32))
    dummy = ak.from_regular(ak.values_astype(events.event[:,None][:,:0], np.float32))


    variables = getlistofobservables()
    dummy_observables = {var: dummy for var in variables}


    # -------------------------------------------------------------------- #
    #                                 e-mu                                 #
    # -------------------------------------------------------------------- #
    if channel == "emu":

        mask_e_mu      = is_e(p4z1)   &   is_mu(p4z2)
        mask_mu_e      = is_mu(p4z1)  &   is_e(p4z2)

        
    # -------------------------------------------------------------------- #
    #                                e-tau                                 #
    # -------------------------------------------------------------------- #
    elif channel == "etau":
        
        # etau
        mask_e_pi      = is_pi(p4z2)
        mask_e_rho     = (is_rho(p4z2) | is_a1DM2(p4z2))
        mask_e_a1      = (is_a1DM10(p4z2) | is_a1DM11(p4h2))
        

    # -------------------------------------------------------------------- #
    #                               mu-tau                                 #
    # -------------------------------------------------------------------- #
    elif channel == "mutau":

        observables = copy.deepcopy(dummy_observables)
        
        mask_mu_pi     = is_pi(p4h2)
        mask_mu_rho    = (is_rho(p4h2)  | is_a1DM2(p4h2))
        mask_mu_a1     = (is_a1DM10(p4h2) | is_a1DM11(p4h2))



        # mu-pi
        logger.info('mu-pi')
        #zcandP4 = PrepareP4(p4zcandinfo, mask_mu_pi)
        #observable_dict = get_observables_mutau(zcandP4, 'mu', 'pi')
        #observables = wrap(mask_mu_pi, observable_dict, observables)        
        observables = extract_observables(observables, p4zcandinfo, mask_mu_pi, 'mu', 'pi', get_observables_mutau, keylist=['p4z2'])
        
        # mu-rho
        logger.info('mu-rho')
        observables = extract_observables(observables, p4zcandinfo, mask_mu_rho, 'mu', 'rho', get_observables_mutau, keylist=['p4z2'])
        #zcandP4 = PrepareP4(p4zcandinfo, mask_mu_rho)
        #zcandP4, count = clean_rearranged_dict(zcandP4, keylist=['p4z1'])
        #observable_dict = get_observables_mutau(zcandP4, 'mu', 'rho')
        #from IPython import embed; embed()        
        #observable_dict = unwrap(observable_dict, count)
        #observables = wrap(mask_mu_rho, observable_dict, observables)
        # mu-rho
        logger.info('mu-a1')
        observables = extract_observables(observables, p4zcandinfo, mask_mu_a1, 'mu', 'a1', get_observables_mutau, keylist=['p4z2','p4z2pi'])

        
    # -------------------------------------------------------------------- #
    #                              tau-tau                                 #
    # -------------------------------------------------------------------- #
    elif channel == "tautau":

        pass
        
    
    else:
        raise RuntimeError(f"WRONG {channel}")
    
    return events, observables




# ------------ DETECTOR LEVEL ----------- #
@producer(
    uses={
        ProduceObservables, #ProduceDetCosPsi,
    },
    produces={
        *getlistofobservables(),
    },
)
def ProduceRecoObservables(
        self: Producer,
        events: ak.Array,
        p4zcandinfo: dict,
        **kwargs,
) -> ak.Array:

    logger.info("Reco level")
    
    events, observables = self[ProduceObservables](events, p4zcandinfo)

    vars = getlistofobservables()
    for var in vars:
        if var not in observables:
            raise KeyError(f"Observable '{var}' not found in observables")
        val = observables[var]
            
        events = set_ak_column(events,
                               var,
                               ak.values_astype(val, np.float32))

    return events



# ------------ GENERATOR LEVEL ----------- #
@producer(
    uses={
        ProduceObservables, #ProduceDetCosPsi,
    },
    produces={
        *[f'gen_{var}' for var in getlistofobservables()],
    },
)
def ProduceGenObservables(
        self: Producer,
        events: ak.Array,
        p4zcandinfo: dict,
        **kwargs,
) -> ak.Array:

    logger.info("Gen level")

    events, observables = self[ProduceObservables](events, p4zcandinfo)

    vars = getlistofobservables()
    for var in vars:
        if var not in observables:
            raise KeyError(f"Observable '{var}' not found in observables")
        val = observables[var]
            
        events = set_ak_column(events,
                               f'gen_{var}',
                               ak.values_astype(val, np.float32))

    return events
