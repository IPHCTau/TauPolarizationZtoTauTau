# coding: utf-8

"""
Exemplary calibration methods.
"""

from __future__ import annotations

import law
import functools

from columnflow.calibration import Calibrator, calibrator

from columnflow.calibration.cms.met import met_phi_run2, met_phi
from columnflow.calibration.cms.jets import jec, jer
from columnflow.calibration.cms.tau import tec

from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.cms.seeds import (
    deterministic_event_seeds, deterministic_jet_seeds, deterministic_electron_seeds,
)
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import

from columnflow.columnar_util import set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import IF_DATA, IF_MC

from zttpol.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)



# ----------------------------------------- #
#               Base Calibrator             #
# ----------------------------------------- #
@calibrator(
    uses={
        IF_MC(mc_weight),
        deterministic_event_seeds,
        deterministic_jet_seeds,
        deterministic_electron_seeds,
    } | {optional("GenJet.*")} | {attach_coffea_behavior},
    produces={
        IF_MC(mc_weight),
        deterministic_event_seeds,
        deterministic_jet_seeds,
        deterministic_electron_seeds,
    },
    exposed=False,
)
def calibrate_base(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    #task = kwargs["task"]
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)    

    # add deterministic seeds that could (e.g.) be used for smearings
    # seed producers
    # !! as this is the first step, the object collections should still be pt-sorted,
    # !! so no manual sorting needed here (but necessary if, e.g., jec is applied before)
    events = self[deterministic_event_seeds](events, **kwargs)
    events = self[deterministic_jet_seeds](events, **kwargs)
    events = self[deterministic_electron_seeds](events, **kwargs)

    # data/mc specific calibrations
    if self.dataset_inst.is_data:
        # nominal jec
        events = self[self.jec_nominal_cls](events, **kwargs)
    else:
        # for mc, when the nominal shift is requested, apply calibrations with uncertainties (i.e. full), otherwise
        # invoke calibrators configured not to evaluate and save uncertainties
        #if task.global_shift_inst.is_nominal:

        # full jec and jer
        events = self[self.jec_full_cls](events, **kwargs)
            
    return events



@calibrate_base.init
def calibrate_base_init(self: Calibrator, **kwargs) -> None:
    # set the name of the met collection to use
    met_name = self.config_inst.x.met_name
    raw_met_name = self.config_inst.x.raw_met_name

    # derive calibrators to add settings once
    flag = f"custom_calibs_registered_{self.cls_name}"
    if not self.config_inst.x(flag, False):
        def add_calib_cls(name, base, cls_dict=None):
            self.config_inst.set_aux(f"calib_{name}_cls", base.derive(name, cls_dict=cls_dict or {}))
        # jec calibrators
        add_calib_cls("jec_full", jec, cls_dict={
            "mc_only": True,
            "met_name": met_name,
            "raw_met_name": raw_met_name,
        })
        add_calib_cls("jec_nominal", jec, cls_dict={
            "uncertainty_sources": [],
            "met_name": met_name,
            "raw_met_name": raw_met_name,
        })
        add_calib_cls("met_phi", met_phi_run2 if self.config_inst.campaign.x.run == 2 else met_phi)

        # change the flag
        #self.config_inst.set_aux(flag, True)


    # store references to classes
    self.jec_full_cls = self.config_inst.x.calib_jec_full_cls
    self.jec_nominal_cls = self.config_inst.x.calib_jec_nominal_cls
    self.met_phi_cls = self.config_inst.x.calib_met_phi_cls
    
    # collect derived calibrators and add them to the calibrator uses and produces
    derived_calibrators = {
        self.jec_full_cls,
        self.jec_nominal_cls,
        self.met_phi_cls,
    }

    self.uses |= derived_calibrators
    self.produces |= derived_calibrators        
