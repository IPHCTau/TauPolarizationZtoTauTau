# coding: utf-8

"""
Example histogram producer.
"""
import law

from columnflow.histogramming import HistProducer
from columnflow.histogramming.default import cf_default
from columnflow.util import maybe_import
from columnflow.config_util import get_shifts_from_sources
from columnflow.columnar_util import Route

ak = maybe_import("awkward")
np = maybe_import("numpy")

logger = law.logger.get_logger(__name__)


# extend columnflow's default hist producer
@cf_default.hist_producer(
    mc_only=False,
    keep_weights=None,
    drop_weights=None,
)
def main(self: HistProducer, events: ak.Array, **kwargs) -> ak.Array:
    # build the full event weight
    weight = ak.Array(np.ones(len(events), dtype=np.float32))
    logger.info("The following weights will be applied for this dataset")
    #if self.dataset_inst.is_mc and len(events):
    if len(events):
        for column in self.weight_columns:
            logger.info(column)
            weight = weight * Route(column).apply(events)

        logger.info("and, these are the syst variations :")
        logger.info(self.shifts)
            
    return events, weight


@main.init
def main_init(self: HistProducer) -> None:
    #self.weight_columns = {}
    weights = {}
    self.weight_columns = []
    if self.config_inst.x.channel == "mutau":
        weights = self.config_inst.x.event_mutau_weights
    elif self.config_inst.x.channel == "etau":
        weights = self.config_inst.x.event_etau_weights
    elif self.config_inst.x.channel == "tautau":
        weights = self.config_inst.x.event_tautau_weights
    #elif self.config_inst.x.channel == "emu":
    #    weights = self.config_inst.x.event_emu_weights
    else:
        raise RuntimeError(f'Wrong channel : {self.config_inst.x.channel}')

    # helpers to match to kept or dropped weights
    do_keep = pattern_matcher(self.keep_weights) if self.keep_weights else (lambda _: True)
    do_drop = pattern_matcher(self.drop_weights) if self.drop_weights else (lambda _: False)

    for weight_name in weights:
        if not do_keep(weight_name) or do_drop(weight_name):
            continue

        # manually skip weights for samples that do not have lhe info
        if getattr(self, "dataset_inst", None) is not None:

            #if (weight_name != "ff_weight" or weight_name != "closure_weight") and self.dataset_inst.is_data:    
            #    continue
            if self.dataset_inst.is_data:
                if not weight_name in ["ff_weight","ff_cls_corr_weight","ff_ext_corr_weight"]:
                    continue

            # skip pdf weights for samples that dont have lhe weight
            is_lhe_weight = any(
                shift_inst.has_tag("pdf_weight")
                for shift_inst in weights[weight_name] #self.config_inst.x.event_weights[weight_name]
            )
            if self.dataset_inst.has_tag("no_lhe_weights"):
                if weight_name in ["pdf_weight"] or is_lhe_weight:
                    continue

            # zpt weight only for DY samples
            is_zpt_reweight = any(
                shift_inst.has_tag("zpt_reweight")
                for shift_inst in weights[weight_name] #self.config_inst.x.event_weights[weight_name]
            )
            if not self.dataset_inst.has_tag("is_dy"):
                if weight_name in ["zpt_reweight"] or is_zpt_reweight:
                    continue

            # top pt weight only for ttbar samples
            is_top_pt_reweight = any(
                shift_inst.has_tag("top_pt_weight")
                for shift_inst in weights[weight_name] #self.config_inst.x.event_weights[weight_name]
            )
            if not self.dataset_inst.has_tag("is_tt"):
                if weight_name in ["top_pt_weight"]:
                    continue

            # tau-spinner weights are only for signal samples
            is_tauspinner_weight = any(
                shift_inst.has_tag("tauspinner_weight")
                for shift_inst in weights[weight_name] #self.config_inst.x.event_weights[weight_name] 
            )
            if not (self.dataset_inst.has_tag("is_ggf_signal") or self.dataset_inst.has_tag("is_vh_signal") or self.dataset_inst.has_tag("is_vbf_signal")):
                if weight_name in ["tauspinner_weight"] or is_tauspinner_weight:
                    continue

        self.weight_columns.append(weight_name)
        self.uses.add(weight_name)
        self.shifts |= {
            shift_inst.name
            for shift_inst in weights[weight_name] #self.config_inst.x.event_weights[weight_name]
        }


    """
    if self.dataset_inst.is_data:
        return

    if self.config_inst.x.channel == "mutau":
        self.weight_columns = self.config_inst.x.event_mutau_weights
    

    # store column names referring to weights to multiply
    #self.weight_columns |= {"normalization_weight", "muon_weight"}
    self.uses |= self.weight_columns

    # declare shifts that the produced event weight depends on
    shift_sources = {"mu"}
    self.shifts |= set(get_shifts_from_sources(self.config_inst, *shift_sources))
    """
