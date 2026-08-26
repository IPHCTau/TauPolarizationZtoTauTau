# coding: utf-8

"""
Exemplary calibration methods.
"""

from __future__ import annotations

import law
import functools

from columnflow.calibration import Calibrator, calibrator

from columnflow.calibration.cms.met import met_phi_run2, met_phi
from columnflow.calibration.cms.jets import jec, jer_horn_handling

from columnflow.production.cms.mc_weight import mc_weight
#from zttpol.production.extra_weights import scale_mc_weight
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import

from columnflow.columnar_util import set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import IF_DATA, IF_MC

from zttpol.util import IF_RUN2, IF_RUN3, IF_RUN_3_2024

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)

# https://github.com/uhh-cms/hh2bbtautau/blob/master/hbt/calibration/default.py
# pt clamping for the evaluation of the L2L3Residual jec corrections in 2024
def jec_clamp_2024_data_l2l3residual(calibrator, corrector, variable_map):
    # apply only to L2L3Residual for data in 2024
    if (
        calibrator.config_inst.campaign.x.year == 2024 and
        calibrator.dataset_inst.is_data and
        corrector.level == "L2L3Residual"
    ):
        min_eta, max_eta = 2.0, 2.5
        clamp_pt = np.float32(35.0)
        clamp_mask = (
            ((abs_eta := abs(variable_map["JetEta"])) > min_eta) &
            (abs_eta < max_eta) &
            (variable_map["JetPt"] > clamp_pt)
        )
        variable_map["JetPt"] = ak.where(clamp_mask, clamp_pt, variable_map["JetPt"])
    return variable_map



# ----------------------------------------- #
#               Base Calibrator             #
# ----------------------------------------- #
@calibrator(
    uses={
        IF_MC(mc_weight),
    },
    produces={
        IF_MC(mc_weight),
    },
    exposed=False,
)
def calibrate_base(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    task = kwargs["task"]
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)    
        

    # data/mc specific calibrations
    if self.dataset_inst.is_data:
        # nominal jec
        events = self[self.jec_nominal_cls](events, **kwargs)
    else:
        # for mc, when the nominal shift is requested, apply calibrations with uncertainties (i.e. full), otherwise
        # invoke calibrators configured not to evaluate and save uncertainties
        if task.global_shift_inst.is_nominal:
            # full jec and jer
            events = self[self.jec_full_cls](events, **kwargs)
            events = self[self.jer_jec_full_cls](events, **kwargs)
        else:
            # full jec and jer
            # nominal jec and jer
            events = self[self.jec_nominal_cls](events, **kwargs)
            events = self[self.jer_jec_nominal_cls](events, **kwargs)

    # apply met phi correction
    if self.has_dep(self.met_phi_cls):
        events = self[self.met_phi_cls](events, **kwargs)

            
    return events



@calibrate_base.init
def calibrate_base_init(self: Calibrator, **kwargs) -> None:
    super(calibrate_base, self).init_func(**kwargs)
    
    # set the name of the met collection to use
    met_name = self.config_inst.x.met_name
    raw_met_name = self.config_inst.x.raw_met_name

    # derive calibrators to add settings once
    flag = f"custom_calibs_registered_{self.cls_name}"
    if not self.config_inst.x(flag, False):
        def add_calib_cls(name, base, cls_dict=None):
            self.config_inst.set_aux(f"calib_{name}_cls", base.derive(name, cls_dict=cls_dict or {}))
        # jec calibrators
        update_jec_corrector_variables = (
            jec_clamp_2024_data_l2l3residual
            if self.config_inst.campaign.x.year == 2024
            else None
        )
        add_calib_cls("jec_full", jec, cls_dict={
            "mc_only": True,
            "met_name": met_name,
            "raw_met_name": raw_met_name,
            "update_corrector_variables": update_jec_corrector_variables,
        })
        add_calib_cls("jec_nominal", jec, cls_dict={
            "uncertainty_sources": [],
            "met_name": met_name,
            "raw_met_name": raw_met_name,
            "update_corrector_variables": update_jec_corrector_variables,
        })
        add_calib_cls("jer_jec_full", jer_horn_handling, cls_dict={
            "met_name": met_name,
        })
        add_calib_cls("jer_jec_nominal", jer_horn_handling, cls_dict={
            "met_name": met_name,
            "jec_uncertainty_sources": [],
        })
        add_calib_cls("met_phi", met_phi_run2 if self.config_inst.campaign.x.run == 2 else met_phi)

        # change the flag
        self.config_inst.set_aux(flag, True)


    # store references to classes
    self.jec_full_cls = self.config_inst.x.calib_jec_full_cls
    self.jec_nominal_cls = self.config_inst.x.calib_jec_nominal_cls
    self.jer_jec_full_cls = self.config_inst.x.calib_jer_jec_full_cls
    self.jer_jec_nominal_cls = self.config_inst.x.calib_jer_jec_nominal_cls
    self.met_phi_cls = self.config_inst.x.calib_met_phi_cls
    
    # collect derived calibrators and add them to the calibrator uses and produces
    derived_calibrators = {
        self.jec_full_cls,
        self.jec_nominal_cls,
        self.jer_jec_full_cls,
        self.jer_jec_nominal_cls,
        # TODO: 2024: remove condition when met phi corrections are made available
        ~IF_RUN_3_2024(self.met_phi_cls),
    }

    self.uses |= derived_calibrators
    self.produces |= derived_calibrators        
