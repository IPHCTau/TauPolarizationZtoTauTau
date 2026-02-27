# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.normalization import stitched_normalization_weights

from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.pdf import pdf_weights
from columnflow.production.cms.top_pt_weight import top_pt_weight

from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

from columnflow.config_util import get_events_from_categories
from httcp.production.ReArrangeZcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP
#from httcp.production.weights import tauspinner_weight
#from httcp.production.extra_weights import zpt_reweight,zpt_reweight_v2, ff_weight, met_recoil_corr, top_pt_weight, classify_events

from httcp.production.extra_weights import zpt_reweight,zpt_reweight_v2, met_recoil_corr, top_pt_weight #, classify_events
#from httcp.production.muon_weights import muon_id_weights, muon_iso_weights, muon_trigger_weights, muon_xtrigger_weights
#from httcp.production.electron_weights import electron_idiso_weights, electron_trigger_weights, electron_xtrigger_weights
#from httcp.production.tau_weights import tau_all_weights, tauspinner_weights
#from httcp.production.tau_weights import tauspinner_weights


from httcp.production.sample_split import split_dy
from httcp.production.processes import build_abcd_masks

from httcp.production.columnvalid import make_column_valid

#from httcp.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi

from httcp.util import IF_DATASET_HAS_LHE_WEIGHTS, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL, IF_DATASET_IS_TT
from httcp.util import IF_RUN2, IF_RUN3, IF_ALLOW_STITCHING, IF_GENMATCH, IF_GENMATCH_ON_FOR_SIGNAL, transverse_mass

from httcp.production.applyFastMTT import apply_fastMTT

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


@producer(
    uses={
        # nano columns
        "zcand.*", optional("GenTau.*"), optional("GenTauProd.*"),
        "Jet.pt",
        "PuppiMET.pt", "PuppiMET.phi",
        reArrangeDecayProducts,
        #IF_GENMATCH_ON_FOR_SIGNAL(reArrangeGenDecayProducts),
        #IF_GENMATCH_ON_FOR_SIGNAL(ProduceGenPhiCP), ####ProduceGenCosPsi, 
        IF_GENMATCH(reArrangeGenDecayProducts),
        IF_GENMATCH(ProduceGenPhiCP), ####ProduceGenCosPsi, 
        ProduceDetPhiCP, ####ProduceDetCosPsi,
        apply_fastMTT,
    },
    produces={
        # new columns
        "zcand_invm",
        "zcand_dr",
        "zcand_dphi",
        "n_jet",
        #IF_GENMATCH_ON_FOR_SIGNAL(ProduceGenPhiCP), ####ProduceGenCosPsi,
        IF_GENMATCH(ProduceGenPhiCP), ####ProduceGenCosPsi,
        ProduceDetPhiCP, ####ProduceDetCosPsi,
        "dphi_met_z1", "dphi_met_z2",
        "met_var_qcd_z1", "met_var_qcd_z2",
        "hT",
        "pt_tt", "pt_vis",
        "mt_1", "mt_2", "mt_lep", "mt_tot",
        apply_fastMTT,
    },
)
def zcand_features(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    zcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    zcand1 = zcand_[:,0:1]
    zcand2 = zcand_[:,1:2]

    vis_p4 = zcand1 + zcand2
    
    met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")

    mass = (zcand1 + zcand2).mass
    dr = zcand1.delta_r(zcand2)
    dphi = zcand1.delta_phi(zcand2)

    # deltaPhi between MET and zcand1 & 2
    dphi_met_z1 = met.delta_phi(zcand1)
    dphi_met_z2 = met.delta_phi(zcand2)
    met_var_qcd_z1 = met.pt * np.cos(dphi_met_z1)/zcand1.pt
    met_var_qcd_z2 = met.pt * np.cos(dphi_met_z2)/zcand2.pt

    # scalar sum pt  
    hT = ak.sum(events.Jet.pt, axis=1)
    #hT = ak.where(ak.num(events.Jet.pt, axis=1) > 0, hT, events.Jet.pt[:,:1])
    events = set_ak_column_f32(events, "hT",  hT)

    events = set_ak_column_f32(events, "zcand_invm",  mass)
    events = set_ak_column_f32(events, "zcand_dr",    dr)
    events = set_ak_column_f32(events, "zcand_dphi",  dphi)
    events = set_ak_column_f32(events, "dphi_met_z1", np.abs(dphi_met_z1))
    events = set_ak_column_f32(events, "dphi_met_z2", np.abs(dphi_met_z2))
    events = set_ak_column_f32(events, "met_var_qcd_z1", met_var_qcd_z1)
    events = set_ak_column_f32(events, "met_var_qcd_z2", met_var_qcd_z2)
    
    events = set_ak_column_i32(events, "n_jet", ak.num(events.Jet.pt, axis=1))

    # pt_ll
    trasnverse_momentum = lambda obj1, obj2 : np.sqrt(obj1.pt**2 + obj2.pt**2 + 2 * obj1.pt * obj2.pt * np.cos(obj1.delta_phi(obj2)))

    pt_vis = trasnverse_momentum(zcand1, zcand2)
    events = set_ak_column_f32(events, "pt_vis", pt_vis)
    pt_tt  = trasnverse_momentum(vis_p4, met)
    events = set_ak_column_f32(events, "pt_tt", pt_tt)

    mt_1 = transverse_mass(zcand1, met)
    mt_2 = transverse_mass(zcand2, met)
    mt_lep = transverse_mass(zcand1, zcand2)
    mt_tot = np.sqrt(mt_1**2 + mt_2**2 + mt_lep**2)
    events = set_ak_column_f32(events, "mt_1", mt_1)
    events = set_ak_column_f32(events, "mt_2", mt_2)
    events = set_ak_column_f32(events, "mt_lep", mt_lep)
    events = set_ak_column_f32(events, "mt_tot", mt_tot)
        
    # ################## #
    #     Run FastMTT    #
    # ################## #
    logger.info(" >>>--- FastMTT-Wiktors --->>> [Not as fast as you think]")
    events = self[apply_fastMTT](events, run_fmtt=self.config_inst.x.enable_fastMTT, **kwargs)

    # ########################### #
    # -------- For PhiCP -------- #
    # ########################### #
    events, P4_dict = self[reArrangeDecayProducts](events)
    events   = self[ProduceDetPhiCP](events, P4_dict)
    ###events  = self[ProduceDetCosPsi](events, P4_dict) # for CosPsi only
    
    if self.config_inst.x.extra_tags.genmatch:
        #if "is_signal" in list(self.dataset_inst.aux.keys()):
        events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
        events = self[ProduceGenPhiCP](events, P4_gen_dict) 
        #events = self[ProduceGenCosPsi](events, P4_gen_dict) # for CosPsi only
    # ########################### #
    
    return events




@producer(
    uses={
        make_column_valid,
        attach_coffea_behavior,
        normalization_weights,
        "zcand.decayMode",
        #IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        #muon_id_weights,
        #muon_iso_weights,
        #muon_trigger_weights,
        #muon_xtrigger_weights,
        # -- electron -- #
        #electron_idiso_weights,
        #electron_trigger_weights,
        #electron_xtrigger_weights,
        # -- tau -- #
        #tau_all_weights,
        #IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_DY(zpt_reweight_v2),
        IF_DATASET_IS_TT(top_pt_weight),
        met_recoil_corr,
        zcand_features,
        category_ids,
        build_abcd_masks,
        "channel_id",
        #ff_weight,
        "process_id",
        #classify_events,
    },
    produces={
        make_column_valid,
        attach_coffea_behavior,
        normalization_weights,
        IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        #muon_id_weights,
        #muon_iso_weights,
        #muon_trigger_weights,
        #muon_xtrigger_weights,
        # -- electron -- #
        #electron_idiso_weights,
        #electron_trigger_weights,
        #electron_xtrigger_weights,
        # -- tau -- #
        #tau_all_weights,
        IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_TT(top_pt_weight),
        IF_DATASET_IS_DY(zpt_reweight_v2),
        met_recoil_corr,
        zcand_features,
        category_ids,
        build_abcd_masks,
        #ff_weight,
        "process_id",
        #classify_events,        
    },
)
def produce_base(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[make_column_valid](events, **kwargs)
    
    events = self[attach_coffea_behavior](events, **kwargs)
    
    #if self.dataset_inst.is_mc:
    if self.dataset_inst.has_tag("is_dy") or self.dataset_inst.has_tag("is_w") or self.dataset_inst.has_tag("is_signal"):
        events = self[met_recoil_corr](events, **kwargs)

    events = self[zcand_features](events, **kwargs)

    #logger.info(" >>>--- Evaluate Classifier Models (IC) --->>> [In extra_weights.py and processes.py]")
    #events = self[classify_events](events, **kwargs)
    
    logger.warning("NO b-veto cut for tautau categories : Imperial")
    events = self[build_abcd_masks](events, **kwargs)
    # building category ids
    events, category_ids_debug_dict = self[category_ids](events, debug=False)

    # debugging categories
    if self.config_inst.x.verbose.production.main:
        from zttpol.production.debug import category_flow
        #print(category_ids_debug_dict)
        od_cats_in_config = self.config_inst.categories
        cat = [cat for cat in od_cats_in_config if cat.name == self.config_inst.x.channel][0]
        cat_events = get_events_from_categories(events, cat)
        logger.info(f"Analysis : {self.config_inst.x.channel}")
        category_flow(self.config_inst.x.channel, cat_events)


    if self.dataset_inst.is_mc:
        # allow stitching is applicable only when datasets are DY or wjets, only if the stitching booleans are true in config
        allow_stitching = bool(ak.any([(self.dataset_inst.has_tag("is_dy_m50") and self.config_inst.x.allow_dy_stitching),
                                       (self.dataset_inst.has_tag("is_w") and self.config_inst.x.allow_w_stitching)]))
        #from IPython import embed; embed()
        if allow_stitching:
            events = self[stitched_normalization_weights](events, **kwargs)
        else:
            events = self[normalization_weights](events, **kwargs)
        
        # TODO : pileup weight is constrained to max value 10
        # TODO : check columnflow production/pileup
        events = self[pu_weight](events, **kwargs)
        if self.has_dep(pdf_weights):
            events = self[pdf_weights](events, **kwargs)

        # -- Z-pT reweighting with corrections from Imperial (Danny) and Kansas (Dennis : v2)
        #if self.has_dep(zpt_reweight):
        if self.has_dep(zpt_reweight_v2):
            #events = self[zpt_reweight](events, **kwargs)
            events = self[zpt_reweight_v2](events, **kwargs)

        #from IPython import embed; embed()
        #processes = self.dataset_inst.processes.names()
        
        #if self.dataset_inst.has_tag("is_dy"):
        #    logger.warning(f"splitting Drell-Yan dataset <{self.dataset_inst.name}>")
        #    events = self[split_dy](events,**kwargs)

        # top pt weight
        if self.has_dep(top_pt_weight):
            events = self[top_pt_weight](events, **kwargs)
        else:
            logger.warning(f"No top pt reweihting for <{self.dataset_inst.name}> dataset")

    #events = self[ff_weight](events, **kwargs)        
    
    return events
