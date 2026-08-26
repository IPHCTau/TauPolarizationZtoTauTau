# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools
import itertools

import os
import re
import law
import yaml
from glob import glob
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.columnar_util import ColumnCollection, EMPTY_FLOAT, skip_column
from columnflow.cms_util import CMSDatasetInfo

from columnflow.tasks.external import ExternalFile as Ext

from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category,
    add_shift_aliases, 
    get_shifts_from_sources,
    verify_config_processes,
)

logger = law.logger.get_logger(__name__)


ak = maybe_import("awkward")


#thisdir = os.path.dirname(os.path.abspath(__file__))
basecorrdir  = "/eos/project/i/iphctau/public/common/JSONPOG_Corrections"
corrdir      = os.path.join(basecorrdir, "metadata_26052026")
tooldir      = os.path.join(basecorrdir, "tools_26052026")


#thisdir = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/TauPolarizationZtoTauTau/zttpol/config"
#tooldir = os.path.join(basecorrdir, "data/tools")
#corrdir = "/eos/project/i/iphctau/public/common/JSONPOG_Corrections/metadata_26052026"

def add_config (ana: od.Analysis,
                campaign: od.Campaign,
                config_name           = None,
                config_id             = None,
                limit_dataset_files   = None,
                channel               = None,
) -> od.Config :

    #from IPython import embed; embed()
    logger.info(f"Channel to be analyzed: {channel}")
    if channel not in {'ee','mumu','emu','etau','mutau','tautau'}:
        raise RuntimeError(f"WRONG CHANNEL : {channel}")
    
    # gather campaign data
    run     = campaign.x.run
    year    = campaign.x.year
    year2   = year % 100
    vnano   = campaign.x.version
    postfix = campaign.x.postfix
    
    # some validations
    assert run in {2,3}
    assert year in {2016,2017,2018,2022,2023,2024,2025}

    year_postfix = ""
    if year == 2022:
        year_postfix = "EE" if postfix == "postEE" else ""
    elif year == 2023:
        year_postfix = "BPix" if postfix == "postBPix" else ""

    dummy_proc_id = int(f'{year2}{vnano}3000') # because, all the datasetIDs are like 18150101, after 1815, its zero.. for customIDs its 3 instead of 0
        
    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    #for p in procs: print(p)
    
    # --------------------------------------------------------------------------- #
    #  create a config by passing the campaign, so id and name will be identical  #
    # --------------------------------------------------------------------------- #

    cfg = ana.add_config(campaign,
                         name  = config_name,
                         id    = config_id)

    # add channel
    cfg.x.channel = channel

    cfg.add_channel(name="emu",    id=1)
    cfg.add_channel(name="etau",   id=2)
    cfg.add_channel(name="mutau",  id=3)
    cfg.add_channel(name="tautau", id=4)
    cfg.add_channel(name="ee",     id=5)
    cfg.add_channel(name="mumu",   id=6)


    # combination of processes
    cfg.add_process(
        name="multiboson",
        id=dummy_proc_id + 1,
        label="VV(V)",
        processes=[procs.n.ww_dl, procs.n.ww_sl,
                   procs.n.wz_wlnu_zll, procs.n.wz_wqq_zll, procs.n.wz_wlnu_zqq,
                   procs.n.zz_zqq_zll, procs.n.zz_zll_znunu, procs.n.zz_zll_zll,
                   procs.n.www, procs.n.wwz, procs.n.wzz, procs.n.zzz,
                   procs.n.ewk_wp_lnu_m50toinf, procs.n.ewk_wm_lnu_m50toinf, procs.n.ewk_z_ll_m50toinf],
    )
    cfg.add_process(
        name="top",
        id=dummy_proc_id + 11, 
        label="top",
        processes=[procs.n.tt,
                   procs.n.st],
    )
    cfg.add_process(
        name="higgs",
        id=dummy_proc_id + 21, 
        label="higgs",
        processes=[procs.n.h_ggf_htt,
                   procs.n.zh_htt,
                   procs.n.wph_htt,
                   procs.n.wmh_htt],
    )

    # ---------------------------------- #
    # add processes we are interested in #
    # ---------------------------------- #
    process_names = [
        ## Data
        "data",
        ## Drell-Yan
        "dy_m50toinf",
        "dy_2tau_m50toinf_nj",
        "dy_2tau_m50toinf_nj_LHEspin_minus",
        "dy_2tau_m50toinf_nj_LHEspin_plus",
        ## tt
        "tt",
        ## st
        "st",
        ## Wlnu
        "w_lnu",
        ## VV
        "ww_dl",
        "ww_sl",
        "wz_wlnu_zll",
        "wz_wqq_zll",
        "wz_wlnu_zqq",
        "zz_zqq_zll",
        "zz_zll_znunu",
        "zz_zll_zll",
        ## VVV
        "www",
        "wwz",
        "wzz",
        "zzz",
        ## EWK
        "ewk_wp_lnu_m50toinf",
        "ewk_wm_lnu_m50toinf",
        "ewk_z_ll_m50toinf",
        ## higgs
        "h_ggf_htt",
        "zh_htt",
        "wph_htt",
        "wmh_htt",
    ]

    for process_name in process_names:
        if process_name in procs:
            proc = procs.get(process_name)
        elif process_name == "qcd":
            # qcd is not part of procs since there is no dataset registered for it
            from cmsdb.processes.qcd import qcd
            proc = cfg.add_process(qcd)
        else:
            logger.warning(f"WARNING: {process_name} not in campaign root processes")
            continue

        # add tags to processes
        if re.match(r"^tt(|_.+)$", process_name):
            proc.add_tag({"ttbar", "tt"})


        # add the process
        cfg.add_process(proc)

        
    # configuration of colors, labels, etc. for plotting
    from zttpol.config.styles import stylize_processes
    stylize_processes(cfg)    

    
    # ----------------------------- #
    # add datasets we need to study #
    # ----------------------------- #

    dataset_names = [
        ## Drell-Yan
        "dy_lep_m50_amcatnlo",
        "dy_2tau_m50_amcatnlo",
        ## tt
        "tt_sl_powheg",
        "tt_dl_powheg",
        "tt_fh_powheg",
        ## st
        "st_tchannel_t_powheg",
        "st_tchannel_tbar_powheg",
        "st_tw_t_powheg",
        "st_tw_tbar_powheg",
        ## wlnu
        "wj_incl_madgraph",
        ## vv
        "ww_dl_powheg",
        "ww_sl_powheg",
        "wz_wlnu_zll_powheg",
        "wz_wqq_zll_powheg",
        "wz_wlnu_zqq_powheg",
        "zz_zqq_zll_powheg",
        "zz_zll_zll_powheg",
        ## vvv
        "www_amcatnlo",
        "wwz_amcatnlo",
        "wzz_amcatnlo",
        "zzz_amcatnlo",
        ## ewk
        "ewk_wp_lnu_madgraph",
        "ewk_wm_lnu_madgraph",
        "ewk_z_ll_madgraph",
        ## higgs
        "h_ggf_htt_powheg",
        "zh_htt_powheg",
        "wph_htt_powheg",
        "wmh_htt_powheg"
    ]

    if channel in {'ee','mumu'}:
        names_to_remove = [
            'dy_2tau_m50_amcatnlo',
            'wj_incl_madgraph',
            'ewk_wp_lnu_madgraph',
            'ewk_wm_lnu_madgraph',
            'h_ggf_htt_powheg',
            'zh_htt_powheg',
            'wph_htt_powheg',
            'wmh_htt_powheg',
        ]
        dataset_names = [x for x in dataset_names if x not in names_to_remove]
        

    
    # Not relevant here
    # Todo : Need to use the filter efficiency for DY dataset for tau polarization analysis
    cfg.x.signal_filter_efficiency = DotDict.wrap({
        "h_ggf_tautau_uncorrelatedDecay_SM_Filtered_ProdAndDecay": 0.3847,
        "h_ggf_tautau_uncorrelatedDecay_MM_Filtered_ProdAndDecay": 0.3848,
        "h_ggf_tautau_uncorrelatedDecay_CPodd_Filtered_ProdAndDecay": 0.3848,
        "wph_tautau_uncorrelatedDecay_Filtered": 0.3743,
        "wmh_tautau_uncorrelatedDecay_Filtered": 0.3944,
        "zh_tautau_uncorrelatedDecay_Filtered": 0.3933,
        "h_vbf_tautau_uncorrelatedDecay_Filtered": 0.4091,
    })

    local_file_path = {
        '2018': '/eos/cms/store/group/phys_tau/alebihan/taupola_noskim/Run2_2018',
    }[f'{year}{postfix}']
    
    #datasets_data = []
    #if year == 2018:
    datasets_data = {
        'ee'     : {
            '2018':         ["data_egamma_A","data_egamma_B","data_egamma_C","data_egamma_D"],
            '2022preEE':    ["data_e_C","data_e_D"],
            '2022postEE':   ["data_e_E","data_e_F","data_e_G"],
            '2023preBPix':  ["data_e0_C","data_e1_C"],
            '2023postBPix': ["data_e0_D","data_e1_D"],
        },
        'mumu'   : {
            '2018':         ["data_single_mu_A", "data_single_mu_B", "data_single_mu_C", "data_single_mu_D"],
            '2022preEE':    ["data_single_mu_C","data_mu_C","data_mu_D"],
            '2022postEE':   ["data_mu_E","data_mu_F","data_mu_G"],
            '2023preBPix':  ["data_mu0_C","data_mu1_C"],
            '2023postBPix': ["data_mu_E","data_mu_F","data_mu_G"],
        },
        'emu'    : {
            '2018':         ["data_muoneg_A", "data_muoneg_B", "data_muoneg_C", "data_muoneg_D"],
            '2022preEE':    [],
            '2022postEE':   [],
            '2023preBPix':  [],
            '2023postBPix': [],
        },
        'etau'   : {
            '2018':         ["data_egamma_A","data_egamma_B","data_egamma_C","data_egamma_D"],
            '2022preEE':    ["data_e_C","data_e_D"],
            '2022postEE':   ["data_e_E","data_e_F","data_e_G"],
            '2023preBPix':  ["data_e0_C","data_e1_C"],
            '2023postBPix': ["data_e0_D","data_e1_D"],
        },
        'mutau'  : {
            '2018':         ["data_single_mu_A", "data_single_mu_B", "data_single_mu_C", "data_single_mu_D"],
            '2022preEE':    ["data_single_mu_C","data_mu_C","data_mu_D"],
            '2022postEE':   ["data_mu_E","data_mu_F","data_mu_G"],
            '2023preBPix':  ["data_mu0_C","data_mu1_C"],
            '2023postBPix': ["data_mu_E","data_mu_F","data_mu_G"],
        },
        'tautau' : {
            '2018':         ["data_tau_A","data_tau_B","data_tau_C","data_tau_D"],
            '2022preEE':    ["data_tau_C","data_tau_D"],
            '2022postEE':   ["data_tau_E","data_tau_F","data_tau_G"],
            '2023preBPix':  ["data_tau_C"],
            '2023postBPix': ["data_tau_D"],
        },
    }[channel][f'{year}{postfix}']

    dataset_names = datasets_data + dataset_names

    logger.info(f"DATASET NAMES : {dataset_names}")
    
    for dataset_name in dataset_names:
        # development switch in case datasets are not _yet_ there
        if dataset_name not in campaign.datasets:
            logger.warning(f"WARNING: {dataset_name} not in cmsdb campaign")
            continue
        
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # datasets that are allowed to contain some events with missing lhe infos
        # (known to happen for amcatnlo)
        if dataset.name.endswith("_amcatnlo"):
            dataset.add_tag("partial_lhe_weights")
        # multiboson samples should not have any lhe weights
        if re.match(r"^(ww|wz|zz)(_|$)", dataset.name):
            dataset.add_tag("no_lhe_weights")
            dataset.remove_tag("partial_lhe_weights")
        # tags for DY
        if re.match(r"^dy_.*$", dataset.name):
            dataset.add_tag("is_dy")
        if re.match(r"^dy_lep_.*$", dataset.name):
            dataset.add_tag("is_dy_lep")
            if channel not in {'ee','mumu'}:
                dataset.add_tag("drop_tautau_from_dy_incl")
        if re.match(r"^dy_2tau_.*$", dataset.name):
            dataset.add_tag("is_dy_tautau") 
        # WJets tags
        if re.match(r"^wj.*$", dataset.name):
            dataset.add_tag("is_w")
        # tt
        if dataset.name.startswith("tt_"):
            dataset.add_tag({"has_top", "ttbar", "tt"})
        # st
        if dataset.name.startswith("st_"):
            dataset.add_tag({"has_top", "single_top", "st"})
        
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    
    # --------------------------------------------------------------------------------------------- #
    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    # --------------------------------------------------------------------------------------------- #

    cfg.x.default_calibrator      = f"calibrate_{channel}"
    logger.info(f"Calibrator : {cfg.x.default_calibrator}")
    cfg.x.default_selector        = f"select_{channel}"
    logger.info(f"Selector : {cfg.x.default_selector}")    
    cfg.x.default_selector_steps  = []
    cfg.x.default_reducer         = "cf_default"
    cfg.x.default_producer        = f"produce_{channel}"
    logger.info(f"Producer : {cfg.x.default_producer}")
    cfg.x.default_hist_producer   = "main"
    cfg.x.default_ml_model        = None
    cfg.x.default_inference_model = "main"
    cfg.x.default_categories      = "main"
    cfg.x.default_variables       = ("event","channel_id")
    #cfg.x.default_weight_producer = "main" # "normalization_only"
    
    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "backgrounds": (backgrounds := [
            "w_lnu",
            "tt",
            "dy_m50toinf", "dy_m10to50",
            "st",
            "vv",
        ]),
        "bkg_sig"       : (bkg_sig       := [*backgrounds, "h_ggf_htt"]),
        "data_bkg"      : (data_bkg      := [*backgrounds, "data"]),
        "data_bkg_sig"  : (data_bkg_sig  := [*backgrounds, "data", "h_ggf_htt"]),
    }
    cfg.x.process_settings_groups = {
        "unstack_processes": {proc: {"unstack": True, "scale": 10.0} for proc in ("h_ggf_htt")},
    }
    cfg.x.general_settings_groups = {
        "compare_shapes": {"skip_ratio": False, "shape_norm": False, "yscale": "log"} #"cms_label": "simpw"},
    }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json",
                    "met_filter",
                    "trigger",
                    "b_veto",
                    "has_2_or_more_leps_with_at_least_1_tau",
                    "dilepton_veto",
                    "has_at_least_1_pair_before_trigobj_matching",
                    "has_at_least_1_pair_after_trigobj_matching",
                    "extra_lepton_veto",
                    "One_higgs_cand_per_event",
                    "has_proper_tau_decay_products"],
    }

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    cfg.x.allow_dy_stitching = False
    cfg.x.allow_dy_stitching_for_plotting = False
    """
    # define inclusive datasets for the stitched process identification with corresponding leaf processes
    # drell-yan [NLO]
    cfg.x.dy_stitching = {
        "dy": {
            "inclusive_dataset": cfg.datasets.n.dy_lep_m50_amcatnlo,
            "leaf_processes": [
                # the following processes cover the full njet phasespace
                *(
                    procs.get(f"dy_m50toinf_{nj}j")
                    for nj in [0, 1, 2]
                ),
            ],
        },
    }

    """

    cfg.x.allow_w_stitching = False
    cfg.x.allow_w_stitching_for_plotting = False
    """
    cfg.x.w_stitching = {
        "wj": {
            "inclusive_dataset": cfg.datasets.n.wj_incl_amcatnlo,
            "leaf_processes": [
                # the following processes cover the full njet and pt phasespace
                *(
                    procs.get(f"w_lnu_{nj}j") # njet from NLO samples
                    for nj in [0,1,2]
                ),
            ],
        },
    }
    """

    # plotting overwrites
    from zttpol.config.styles import setup_plot_styles
    setup_plot_styles(cfg)
    
    # --------------------------------------------------------------------------------------------- #
    # Luminosity and Normalization
    # lumi values in inverse fb
    # TODO later: preliminary luminosity using norm tag. Must be corrected, when more data is available
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # --------------------------------------------------------------------------------------------- #

    if year == 2018:
        cfg.x.luminosity = Number(59_830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })
    
    elif year == 2022:
        if postfix == "preEE":
            cfg.x.luminosity = Number(7_980.4541, {
                "lumi_13p6TeV_2022": 0.014j,
                "lumi_13p6TeV_22_23_24": 0.0138j,
            })
        elif postfix == "postEE":
            cfg.x.luminosity = Number(26_671.6097, {
                "lumi_13p6TeV_2022": 0.014j,
                "lumi_13p6TeV_22_23_24": 0.0138j,
            })
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2023:
        if postfix == "preBPix":
            cfg.x.luminosity = Number(18_062.6591, {
                "lumi_13p6TeV_2023": 0.013j,
                "lumi_13p6TeV_22_23_24": 0.0017j,
                "lumi_13p6TeV_23_24": 0.0127j,
            })
        elif postfix == "postBPix":
            cfg.x.luminosity = Number(9_693.1301, {
                "lumi_13p6TeV_2023": 0.013j,
                "lumi_13p6TeV_22_23_24": 0.0017j,
                "lumi_13p6TeV_23_24": 0.0127j,
            })
    elif year == 2024:
        cfg.x.luminosity = Number(109_948.177486, {
            "lumi_13p6TeV_2024": 0.016j,
            "lumi_13p6TeV_22_23_24": 0.0020j,
            "lumi_13p6TeV_23_24": 0.0068j,
            "lumi_13p6TeV_24": 0.0144j,
        })
    else:
        raise RuntimeError(f"Wrong year: {year}")


    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    # TODO later: Error for run three not available yet. Using error from run 2.
    cfg.x.minbias_xs = Number(69.2, 0.046j)


    # --------------------------------------------------------------------------------------------- #
    # Adding triggers
    # --------------------------------------------------------------------------------------------- #

    from zttpol.config import (
        triggers_emu,
        triggers_etau,
        triggers_mutau,
        triggers_tautau,
        triggers_ee,
        triggers_mumu,
    )

    def add_triggers(cfg, year, channel, postfix):
        modules = {
            "emu": triggers_emu,
            "etau": triggers_etau,
            "mutau": triggers_mutau,
            "tautau": triggers_tautau,
            "ee" : triggers_ee,
            "mumu" : triggers_mumu,
        }

        try:
            module = modules[channel]
        except KeyError:
            raise ValueError(f"Unknown channel '{channel}'")
        
        funcs = {
            2018: module.add_triggers_2018,
            2022: module.add_triggers_2022,
            2023: module.add_triggers_2023,
        }
        
        try:
            funcs[year](cfg, postfix)
        except KeyError:
            raise ValueError(f"Unsupported year '{year}'")
        

    add_triggers(cfg, year, channel, postfix)
        


    # ------------------------------------------------------------- #
    # bjet settings                                                 #
    # b-tag working points                                          #
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/      #
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/    #
    # TODO later: complete WP when data becomes available           #
    # ------------------------------------------------------------- #

    btag_key = f"{year}{year_postfix}"
    if run == 2:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP?rev=6
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP?rev=8
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17?rev=15
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18?rev=18
        cfg.x.btag_working_points = DotDict.wrap({
            "deepcsv": {
                "loose": {"2016APV": 0.2027, "2016": 0.1918, "2017": 0.1355, "2018": 0.1208}[btag_key],
                "medium": {"2016APV": 0.6001, "2016": 0.5847, "2017": 0.4506, "2018": 0.4168}[btag_key],
                "tight": {"2016APV": 0.8819, "2016": 0.8767, "2017": 0.7738, "2018": 0.7665}[btag_key],
            },
            # https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005 chapter 4.5 in v12
            "particleNetMD": {
                "hp": {"2016APV": 0.9883, "2016": 0.9883, "2017": 0.9870, "2018": 0.9880}[btag_key],
                "mp": {"2016APV": 0.9737, "2016": 0.9735, "2017": 0.9714, "2018": 0.9734}[btag_key],
                "lp": {"2016APV": 0.9088, "2016": 0.9137, "2017": 0.9105, "2018": 0.9172}[btag_key],
            },
            # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2018NanoAODv15/#ak4-b-tagging
            "upart": {
                "loose"   : {'2018': 0.0308}[btag_key],
                "medium"  : {'2018': 0.1610}[btag_key],
                "tight"   : {'2018': 0.5405}[btag_key],
                "vtight"  : {'2018': 0.6992}[btag_key],
                "vvtight" : {'2018': 0.9655}[btag_key]
            }
        })
        # btag columns and working point values for easy use throughout the code
        cfg.x.btag_upart = DotDict(
            btv_name="UParTAK4",
            jet_column="btagUParTAK4B",
            wp=cfg.x.btag_working_points.upart.tight,
            weight_column="normalized_njet_btag_weight_upart",
        )
        cfg.x.btag_default = cfg.x.btag_upart

    elif run == 3:
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer23
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer23BPix
        cfg.x.btag_working_points = DotDict.wrap({
            "deepjet": {
                "loose": {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0479, "2023BPix": 0.048, "2024": None}[btag_key],
                "medium": {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.2431, "2023BPix": 0.2435, "2024": None}[btag_key],
                "tight": {"2022": 0.7183, "2022EE": 0.73, "2023": 0.6553, "2023BPix": 0.6563, "2024": None}[btag_key],
                "xtight": {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.7667, "2023BPix": 0.7671, "2024": None}[btag_key],
                "xxtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.9459, "2023BPix": 0.9483, "2024": None}[btag_key],  # noqa: E501
            },
            "particleNet": {
                "loose": {"2022": 0.047, "2022EE": 0.0499, "2023": 0.0358, "2023BPix": 0.0359, "2024": None}[btag_key],
                "medium": {"2022": 0.245, "2022EE": 0.2605, "2023": 0.1917, "2023BPix": 0.1919, "2024": None}[btag_key],  # noqa: E501
                "tight": {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.6172, "2023BPix": 0.6133, "2024": None}[btag_key],  # noqa: E501
                "xtight": {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.7515, "2023BPix": 0.7544, "2024": None}[btag_key],  # noqa: E501
                "xxtight": {"2022": 0.961, "2022EE": 0.9664, "2023": 0.9659, "2023BPix": 0.9688, "2024": None}[btag_key],  # noqa: E501
            },
            # 2024 wps can be taken from "UParTAK4_wp_values" correction set in BTV correctionlib file
            "upart": {
                "loose": {"2022": None, "2022EE": None, "2023": None, "2023BPix": None, "2024": 0.0246}[btag_key],
                "medium": {"2022": None, "2022EE": None, "2023": None, "2023BPix": None, "2024": 0.1272}[btag_key],
                "tight": {"2022": None, "2022EE": None, "2023": None, "2023BPix": None, "2024": 0.4648}[btag_key],
                "xtight": {"2022": None, "2022EE": None, "2023": None, "2023BPix": None, "2024": 0.6298}[btag_key],
                "xxtight": {"2022": None, "2022EE": None, "2023": None, "2023BPix": None, "2024": 0.9739}[btag_key],
            },
        })
        cfg.x.btag_pnet = DotDict(
            btv_name="PNet",
            jet_column="btagPNetB",
            wp=cfg.x.btag_working_points.particleNet.medium,
            weight_column="normalized_njet_btag_weight_pnet",
        )
        cfg.x.btag_upart = DotDict(
            btv_name="UParT",
            jet_column="btagUParTAK4B",
            wp=cfg.x.btag_working_points.upart.medium,
            weight_column="btag_weight",  # no need for normalization in wp based method
        )
        cfg.x.btag_default = cfg.x.btag_upart if year == 2024 else cfg.x.btag_pnet
    else:
        assert False

    

    # --------------- #
    #   MET settings  #
    # --------------- #

    if run == 2:
        #cfg.x.met_name = "MET"
        #cfg.x.raw_met_name = "RawMET"
        cfg.x.met_name = "PuppiMET"
        cfg.x.raw_met_name = "RawPuppiMET"
        
        # met phi correction config
        from columnflow.calibration.cms.met import METPhiConfigRun2
        cfg.x.met_phi_correction = METPhiConfigRun2(
            met_name=cfg.x.met_name,
            correction_set_template="{variable}_metphicorr_pfmet_{data_source}",
            keep_uncorrected=True,
        )
    elif run == 3:
        cfg.x.met_name = "PuppiMET"
        cfg.x.raw_met_name = "RawPuppiMET"

        # met phi correction config
        from columnflow.calibration.cms.met import METPhiConfig
        cfg.x.met_phi_correction = METPhiConfig(
            met_name=cfg.x.met_name,
            met_type=cfg.x.met_name,
            correction_set="met_xy_corrections",
            keep_uncorrected=True,
            pt_phi_variations={
                "stat_xdn": "metphi_statx_down",
                "stat_xup": "metphi_statx_up",
                "stat_ydn": "metphi_staty_down",
                "stat_yup": "metphi_staty_up",
            },
            variations={
                "pu_dn": "minbias_xs_down",
                "pu_up": "minbias_xs_up",
            },
        )
    else:
        assert False        


    # ---------------------- #
    #   Adding met filters   #
    # ---------------------- #
    
    from zttpol.config.met_filters import add_met_filters
    add_met_filters(cfg)


    # --------------------------------------------------------------------------------------------- #
    # jet settings                                                                                  #
    # common jec/jer settings configuration                                                         #
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201                                    #
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107                                #
    # --------------------------------------------------------------------------------------------- #

    from columnflow.calibration.cms.jets import JECConfig, JERConfig
    
    # common jec/jer settings configuration
    if run == 2:
        # https://cms-jerc.web.cern.ch/Recommendations/#run-2
        # https://cms-jerc.web.cern.ch/Recommendations/#2018
        # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=204
        # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=109
        jec_campaign = f"Summer20UL{year2}{postfix}NanoV{vnano}"
        #jec_version = {2016: "V7", 2017: "V5", 2018: "V5"}[year]
        jec_version = {
            2016: "V7",
            2017: "V5",
            2018: "V1"
        }[year]
        jer_campaign = f"Summer{'20' if year == 2016 else '19'}UL{year2}{postfix}"
        jer_version = "JR" + {
            2016: "V3",
            2017: "V2",
            2018: "V2"
        }[year]
        #jet_type = "AK4PFchs"
        jet_type = "AK4PFPuppi"
        
        
    elif run == 3:
        # https://cms-jerc.web.cern.ch/Recommendations/#2022
        jerc_postfix = {
            (2022, ""): "_22Sep2023",
            (2022, "EE"): "_22Sep2023",
            (2023, ""): "Prompt23",
            (2023, "BPix"): "Prompt23",
            (2024, ""): "Prompt24",
        }[(year, year_postfix)]
        jec_campaign = f"Summer{year2}{year_postfix}{jerc_postfix}"
        jec_version = {
            (2022, ""): "V4",
            (2022, "EE"): "V4",
            (2023, ""): "V4",
            (2023, "BPix"): "V4",
            (2024, ""): "V5",
        }[(year, year_postfix)]
        jer_campaign = f"Summer{year2}{year_postfix}{jerc_postfix}"

        # special "Run" fragment in 2023 jer campaign
        if year == 2023:
            jer_campaign += f"_Run{'Cv1234' if postfix == 'preBPix' else 'D'}"
        jer_version = "JR" + {
            2022: "V2",
            2023: "V3",
            2024: "V2",
        }[year]
        jet_type = "AK4PFPuppi"
    else:
        assert False


    # full list of jec sources in a fixed order that is used to assign consistent ids across configs
    # (please add new sources at the bottom to preserve the order of existing ones)
    # the boolean flag decides whether to use them in the JEC config and if shifts should be created for them
    # https://cms-jerc.web.cern.ch/Recommendations/#uncertainites-and-correlations
    jec_source_era = f"{year}{year_postfix}"
    all_jec_sources = {
        "AbsoluteFlavMap": False,
        "AbsoluteMPFBias": False,
        "AbsoluteSample": False,
        "AbsoluteScale": False,
        "AbsoluteStat": False,
        "FlavorPhotonJet": False,
        "FlavorPureBottom": False,
        "FlavorPureCharm": False,
        "FlavorPureGluon": False,
        "FlavorPureQuark": False,
        "FlavorQCD": False,
        "FlavorZJet": False,
        "Fragmentation": False,
        "PileUpDataMC": False,
        "PileUpEnvelope": False,
        "PileUpMuZero": False,
        "PileUpPtBB": False,
        "PileUpPtEC1": False,
        "PileUpPtEC2": False,
        "PileUpPtHF": False,
        "PileUpPtRef": False,
        "RelativeBal": False,
        "RelativeFSR": False,
        "RelativeJEREC1": False,
        "RelativeJEREC2": False,
        "RelativeJERHF": False,
        "RelativePtBB": False,
        "RelativePtEC1": False,
        "RelativePtEC2": False,
        "RelativePtHF": False,
        "RelativeSample": False,
        "RelativeStatEC": False,
        "RelativeStatFSR": False,
        "RelativeStatHF": False,
        "SinglePionECAL": False,
        "SinglePionHCAL": False,
        "SubTotalAbsolute": False,
        "SubTotalMC": False,
        "SubTotalPileUp": False,
        "SubTotalPt": False,
        "SubTotalRelative": False,
        "SubTotalScale": False,
        "TimePtEta": False,
        "Total": True,
        "TotalNoFlavor": False,
        "TotalNoFlavorNoTime": False,
        "TotalNoTime": False,
        "CorrelationGroupFlavor": False,
        "CorrelationGroupIntercalibration": False,
        "CorrelationGroupMPFInSitu": False,
        "CorrelationGroupUncorrelated": False,
        "CorrelationGroupbJES": False,
        "Regrouped_Absolute": True,
        f"Regrouped_Absolute_{jec_source_era}": True,
        "Regrouped_BBEC1": True,
        f"Regrouped_BBEC1_{jec_source_era}": True,
        "Regrouped_EC2": True,
        f"Regrouped_EC2_{jec_source_era}": True,
        "Regrouped_FlavorQCD": True,
        "Regrouped_HF": True,
        f"Regrouped_HF_{jec_source_era}": True,
        "Regrouped_RelativeBal": True,
        f"Regrouped_RelativeSample_{jec_source_era}": True,
        "Regrouped_Total": True,
    }

    bjec_config = None
    # JEC
    cfg.x.jec = DotDict.wrap({
        "Jet": JECConfig(
            jet_name="Jet",
            jet_type=jet_type,
            campaign=jec_campaign,
            version=jec_version,
            levels=["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            levels_for_type1_met=["L1FastJet"],
            data_per_era=False,  # no more era-dependence in latest jec campaigns
            uncertainty_sources=[src for src, flag in all_jec_sources.items() if flag],
            bjec_config=bjec_config,
        ),
    })

    # JER
    cfg.x.jer = DotDict.wrap({
        "Jet": JERConfig(
            jet_name="Jet",
            jet_type=jet_type,
            campaign=jer_campaign,
            version=jer_version,
            use_jer_tool=True,
        ),
    })

    """
    cfg.x.jec = DotDict.wrap({
        "Jet": {
            "campaign": jec_campaign,
            "version": jec_version,
            "data_per_era": year == 2022,  # 2022 JEC has the era in the correction set name
            "jet_type": jet_type,
            "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            "levels_for_type1_met": ["L1FastJet"],
            "uncertainty_sources": [src for src, flag in all_jec_sources.items() if flag],
        },
    })

    # JER
    cfg.x.jer = DotDict.wrap({
        "Jet": {
            "campaign": jer_campaign,
            "version": jer_version,
            "jet_type": jet_type,
        },
    })
    """
    # updated jet id
    from columnflow.production.cms.jet import JetIdConfig
    cfg.x.jet_id = JetIdConfig(
        corrections={
            "AK4PUPPI_Tight": 2,
            "AK4PUPPI_TightLeptonVeto": 3
        }
    )
    cfg.x.fatjet_id = JetIdConfig(
        corrections={
            "AK8PUPPI_Tight": 2,
            "AK8PUPPI_TightLeptonVeto": 3
        }
    )

    # trigger sf corrector
    #cfg.x.jet_trigger_corrector = "jetlegSFs"
    # ditau + jet trigger
    cfg.x.jet_trigger_corrector = "jetleg60"




    # ---------------------- #
    # tau settings           #
    # tau-id working points  #
    # ---------------------- #

    #cfg.x.tautag_working_points = DotDict.wrap({
    cfg.x.tauIDWPs = DotDict.wrap({        
        "DeepTau2018v2p5": {
            "vs_e": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4,
                     "Medium": 5,
                     "Tight": 6, "VTight": 7, "VVTight": 8}, # VVVL,VVL,VL,L,M,T,VT,VVT
            "vs_m": {"VLoose": 1, "Loose": 2,
                     "Medium": 3,
                     "Tight": 4}, # VL,L,M,T
            "vs_j": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4,
                     "Medium": 5,
                     "Tight": 6, "VTight": 7, "VVTight": 8}  # VVVL,VVL,VL,L,M,T,VT,VVT
        },
        "PNet": {
            "vs_e": {"VVVLoose": 0.1266, "VVLoose": 0.3547, "VLoose": 0.6997, "Loose": 0.9345,
                     "Medium": 0.9791,
                     "Tight": 0.9897, "VTight": 0.9946, "VVTight": 0.9971}, # VVVL,VVL,VL,L,M,T,VT,VVT
            "vs_m": {"VLoose": 0.2399, "Loose": 0.6037,
                     "Medium": 0.8697,
                     "Tight": 0.9451}, # VL,L,M,T
            "vs_j": {"VVVLoose": 0.0565, "VVLoose": 0.1774, "VLoose": 0.3810, "Loose": 0.6857,
                     "Medium": 0.8347,
                     "Tight": 0.9059, "VTight": 0.9494, "VVTight": 0.9737} # VVVL,VVL,VL,L,M,T,VT,VVT
        }
    })
    
    #cfg.x.deep_tau_wp = DotDict.wrap({
    cfg.x.tauIDWPs_config = DotDict.wrap({
        "DeepTau2018v2p5": {
            "vs_e": {
                "emu"    : "Tight", # need to check
                "etau"   : "Tight",
                "mutau"  : "VVLoose",
                "tautau" : "VVLoose",
                "mumu"   : "Tight",
                "ee"     : "Tight",
            },
            "vs_m": {
                "emu"    : "Tight", # need to check
                "etau"   : "Tight",
                "mutau"  : "Tight",
                "tautau" : "Tight", #"VLoose",  #"Tight", #"VLoose",
                "mumu"   : "Tight",
                "ee"     : "Tight",
            },
            "vs_j": {
                "emu"    : "VTight", # need to check
                "etau"   : "VTight",
                "mutau"  : "VTight", ##"Medium" : OLD,
                "tautau" : "VTight", #"VLoose", #"Medium", #"VTight", ## VTight : Proposed by Imperial, was Medium in Run2
                "mumu"   : "Tight",
                "ee"     : "Tight",
            },
        }
    })
    
    """    
    cfg.x.deep_tau_info = DotDict.wrap({
        "DeepTau2018v2p5": {
            "wp": {
                "vs_e": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
                "vs_m": {"VVVLoose": 1, "VVLoose": 1, "VLoose": 1, "Loose": 2, "Medium": 3, "Tight": 4, "VTight": 4, "VVTight": 4},
                "vs_j": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
            },
            "vs_e": {
                "emu"    : "Tight", # need to check
                "etau"   : "Tight",
                "mutau"  : "VVLoose",
                "tautau" : "VVLoose",
            },
            # All DeepTauVsMu WPs are changed to TIGHT WPs (studied by IC)
            "vs_m": {
                "emu"    : "Tight", # need to check
                "etau"   : "Tight",
                "mutau"  : "Tight",
                "tautau" : "Tight", #"VLoose",  #"Tight", #"VLoose",
            },
            "vs_j": {
                "emu"    : "VTight", # need to check
                "etau"   : "VTight",
                "mutau"  : "VTight", ##"Medium" : OLD,
                "tautau" : "VTight", #"VLoose", #"Medium", #"VTight", ## VTight : Proposed by Imperial, was Medium in Run2
                # W A R N I N G !!! Medium is being used for ML, Change it to VTight for CP analysis
                #"tautau" : "Medium",
            },
        },
    })
    """

    cfg.x.deep_tau_tagger = "DeepTau2018v2p5"
       
    # tec config
    from columnflow.calibration.cms.tau import TECConfig
    corrector_kwargs = {
        ("ee", 2)    : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['ee'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['ee']},        
        ("mumu", 2)  : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['mumu'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['mumu']},        
        ("mutau", 2) : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['mutau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['mutau']},
        ("etau", 2)  : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['etau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['etau']},
        ("tautau", 2): {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['tautau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['tautau']},
        ("etau", 3)  : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['etau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['etau']},
        ("mutau", 3) : {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['mutau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['mutau']},
        ("tautau", 3): {"wp": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_j['tautau'],
                        "wp_VSe": cfg.x.tauIDWPs_config[cfg.x.deep_tau_tagger].vs_e['tautau']},
    }[(channel, run)]
    #print(corrector_kwargs)
    cfg.x.tec = TECConfig(tagger=cfg.x.deep_tau_tagger,
                          #correction_set=f'tau_es_dm_{cfg.x.deep_tau_tagger}_{year}_{postfix}',
                          corrector_kwargs=corrector_kwargs)



    

    # --------------------------------------------------------- #
    #   Adding external files e.g. goodlumi, normtag, SFs etc   #
    #  restructure the postfix names to build appropriate tags  #
    # --------------------------------------------------------- #

    cfg.x.external_files = DotDict()

    # helpers
    def wrap_ext(obj):
        if isinstance(obj, Ext):
            return obj
        if isinstance(obj, tuple):
            if len(obj) != 2:
                raise ValueError(f"cannot wrap tuple '{obj}' into ExternalFile, expected length 2")
            return Ext(location=obj[0], version=obj[1])
        return law.util.map_struct(wrap_ext, obj, map_tuple=False)
    
    def add_external(name, value):
        if isinstance(value, dict):
            value = DotDict.wrap(value)
        cfg.x.external_files[name] = wrap_ext(value)
        return cfg.x.external_files[name]

        
    # adding Goldenlumi JSON and Normtag JSON
    normtagjson = {
        '2018' : ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
        '2022' : ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json", "v1"),
        '2023' : ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json", "v1"),
        '2024' : ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json", "v1")
    }[f'{year}']
    goldenjson = {
        '2018' : ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt","v1"),
        '2022' : ("/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json","v1"),
        '2023' : ("/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json","v1"),
        '2024' : ("/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json""v1")
    }[f'{year}']

    logger.info(f"Normtag JSON : {normtagjson[0]}")
    logger.info(f"Golden JSON  : {goldenjson[0]}")
    add_external(
        "lumi", {
            "golden"  : goldenjson,
            "normtag" : normtagjson
        }
    )

    # pileup weight corrections
    # https://cms-analysis-corrections.docs.cern.ch/corrections/LUM/
    puwtjson = {
        '2018'         : (f"{corrdir}/LUM/Run2-2018-UL-NanoAODv9/latest/puWeights.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/LUM/Run3-22CDSep23-Summer22-NanoAODv12/latest/puWeights.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/LUM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/puWeights.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/LUM/Run3-23CSep23-Summer23-NanoAODv12/latest/puWeights.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/LUM/Run3-23DSep23-Summer23-NanoAODv12/latest/puWeights.json.gz", "v1"),
        '2024'         : (f"{corrdir}/LUM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/puWeights_BCDEFGHI.json.gz", "v1")
    }[f'{year}{postfix}']

    logger.info(f"PileupWt JSON : {puwtjson[0]}")
    add_external("pu_sf", puwtjson)

    # jet energy corrections
    # https://cms-analysis-corrections.docs.cern.ch/corrections/JME/
    jecjson = {
        '2018'         : (f"{corrdir}/JME/Run2-2018-UL-NanoAODv15/latest/jet_jerc.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jet_jerc.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jet_jerc.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/JME/Run3-23DSep23-Summer23-NanoAODv12/latest/jet_jerc.json.gz", "v1"),
        '2024'         : (f"{corrdir}/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jet_jerc.json.gz", "v1")
    }[f'{year}{postfix}']

    logger.info(f"JEC JSON : {jecjson[0]}")
    add_external("jet_jerc", jecjson)

    # jer smearing tool
    add_external("jer_tool", (f"{tooldir}/central_jme_files/jer_smear.json.gz", "v1"))

    logger.warning("No jet_veto_map for Run 2")
    if run == 3:
        # jet veto map
        jetvetomap = {
            #'2018'         : (f"{corrdir}/JME/Run2-2018-UL-NanoAODv15/latest/jetvetomaps.json.gz", "v1"),
            '2022PreEE'    : (f"{corrdir}/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jetvetomaps.json.gz", "v1"),
            '2022PostEE'   : (f"{corrdir}/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jetvetomaps.json.gz", "v1"),
            '2023PreBPix'  : (f"{corrdir}/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jetvetomaps.json.gz", "v1"),
            '2023PostBPix' : (f"{corrdir}/JME/Run3-23DSep23-Summer23-NanoAODv12/latest/jetvetomaps.json.gz", "v1"),
            '2024'         : (f"{corrdir}/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jetvetomaps.json.gz", "v1")        
        }[f'{year}{postfix}']
        
        logger.info(f"JetVetoMap JSON : {jetvetomap[0]}")
        add_external("jet_veto_map", jetvetomap)

    logger.warning("No JetID from JetID JSON for Run 2 ... it will be done manually")
    if run == 3:
        jetidjson = {
            '2022PreEE'    : (f"{corrdir}/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jetid.json.gz", "v1"),
            '2022PostEE'   : (f"{corrdir}/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jetid.json.gz", "v1"),
            '2023PreBPix'  : (f"{corrdir}/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jetid.json.gz", "v1"),
            '2023PostBPix' : (f"{corrdir}/JME/Run3-23DSep23-Summer23-NanoAODv12/latest/jetid.json.gz", "v1"),
            '2024'         : (f"{corrdir}/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jetid.json.gz", "v1")
        }[f'{year}{postfix}']

        logger.info(f"JetID JSON : {jetidjson[0]}")
        add_external("jet_id", jetidjson)
    
    # met phi corrections
    # WARNING : Remember that the met phi correction for 2018 is for v9 instead of v15
    metphijson = {
        '2018'         : (f"{corrdir}/JME/Run2-2018-UL-NanoAODv9/latest/met.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/met_xyCorrections_2022_2022.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/met_xyCorrections_2022_2022EE.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/", "v1"),
        '2023PostBPix' : (f"{corrdir}/JME/Run3-23DSep23-Summer23-NanoAODv12/latest/", "v1")
    }[f'{year}{postfix}']

    logger.info(f"METphi JSON : {metphijson[0]}")
    add_external("met_phi_corr", metphijson)
    
    # muonID SF
    # WARNING : for 2018, muonID SF is from v9, v15 not available yet
    muonIDSF = {
        '2018'         : (f"{corrdir}/MUO/Run2-2018-UL-NanoAODv9/latest/muon_Z.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/MUO/Run3-22CDSep23-Summer22-NanoAODv12/latest/muon_Z.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/MUO/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/muon_Z.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/MUO/Run3-23CSep23-Summer23-NanoAODv12/latest/muon_Z.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/MUO/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/muon_Z.json.gz", "v1"),
        '2024'         : (f"{corrdir}/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_Z.json.gz", "v1")
    }[f'{year}{postfix}']

    logger.info(f"MuonIDSF JSON : {muonIDSF[0]}")
    add_external("muon_sf", muonIDSF)

    logger.warning("No Muon scale and smearing for Run 2, Not Available?")
    if run == 3:
        # muon energy (scale and resolution) corrections and helper tools
        muonscalesmear = {
            '2022PreEE'    : (f"{corrdir}/MUO/Run3-22CDSep23-Summer22-NanoAODv12/latest/muon_scalesmearing.json.gz", "v1"),
            '2022PostEE'   : (f"{corrdir}/MUO/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/muon_scalesmearing.json.gz", "v1"),
            '2023PreBPix'  : (f"{corrdir}/MUO/Run3-23CSep23-Summer23-NanoAODv12/latest/muon_scalesmearing.json.gz", "v1"),
            '2023PostBPix' : (f"{corrdir}/MUO/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/muon_scalesmearing.json.gz", "v1"),
            '2024'	       : (f"{corrdir}/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_scalesmearing.json.gz", "v1")
        }[f'{year}{postfix}']
        
        logger.info(f"MuonSS JSON : {muonscalesmear[0]}")
        add_external("muon_sr", muonscalesmear)
        add_external("muon_sr_tools", Ext(
            f"{tooldir}/muonscarekit-1c7426b5.tar.gz",
            subpaths="muonscarekit-master/scripts/MuonScaRe.py",
            version="v1",
        ))
    
    # electronID SF
    eleIDSF = {
        '2018'         : (f"{corrdir}/EGM/Run2-2018-UL-NanoAODv15/latest/electron.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electron.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electron.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electron.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electron.json.gz", "v1"),
        '2024'         : (f"{corrdir}/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electron.json.gz", "v1")
    }[f'{year}{postfix}']

    logger.info(f"EleIDSF JSON : {eleIDSF[0]}")
    add_external("electron_sf", eleIDSF)
    
    # electron scale and smearing
    elescalesmear = {
        '2018'         : (f"{corrdir}/EGM/Run2-2018-UL-NanoAODv15/latest/electronSS_EtDependent.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
        '2024'         : (f"{corrdir}/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electronSS_EtDependent.json.gz", "v1")
    }[f'{year}{postfix}']

    logger.info(f"EleSS JSON : {elescalesmear[0]}")
    add_external("electron_ss", elescalesmear)
    # EGM tool for simplified ss application
    add_external("egm_tool", (f"{tooldir}/custom_egm_files/egm_tools.json.gz", "v1"))
    
    
    # Not available for Run2?
    if run == 3:
        eletrigsf = {
            '2022PreEE'    : (f"{corrdir}/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electronHlt.json.gz", "v1"),
            '2022PostEE'   : (f"{corrdir}/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
            '2023PreBPix'  : (f"{corrdir}/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
            '2023PostBPix' : (f"{corrdir}/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electronSS_EtDependent.json.gz", "v1"),
            '2024'         : (f"{corrdir}/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electronSS_EtDependent.json.gz", "v1")
        }[f'{year}{postfix}']

        logger.info(f"EleHLT JSON : {eletrigsf[0]}")
        add_external("electron_trig_sf", eletrigsf)

    # tauID SF + TES + TriggerSF
    tauSF = {
        '2018'         : (f"{corrdir}/TAU/Run2-2018-UL-NanoAODv15/latest/tau.json.gz", "v1"),
        '2022PreEE'    : (f"{corrdir}/TAU/Run3-22CDSep23-Summer22-NanoAODv12/latest/tau.json.gz", "v1"),
        '2022PostEE'   : (f"{corrdir}/TAU/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/tau.json.gz", "v1"),
        '2023PreBPix'  : (f"{corrdir}/TAU/Run3-23CSep23-Summer23-NanoAODv12/latest/tau.json.gz", "v1"),
        '2023PostBPix' : (f"{corrdir}/TAU/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/tau.json.gz", "v1"),
        '2024'         : (f"{corrdir}/TAU/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/tau.json.gz", "v1"),
    }[f'{year}{postfix}']

    logger.info(f"TauSF JSON : {tauSF[0]}")
    add_external("tau_sf", tauSF)

    

    """
    cfg.x.external_files = DotDict.wrap({
        "muon_sr_tools"     : (f"{tooldir}/MuonScaReKIT/scripts/MuonScaRe.py", ""),
        # .. https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossMuTauHlt.json?ref_type=heads
        "muon_xtrig_sf"     : (f"{corrdir}/MUO/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/CrossMuTauHlt.json",  "v1"), # Mu xTrig SF
        # .. https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossEleTauHlt.json?ref_type=heads
        "electron_xtrig_sf" : (f"{corrdir}/EGM/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/CrossEleTauHlt.json", "v1"), # Ele xTrig SF
        # https://cms-analysis-corrections.docs.cern.ch/corrections/TAU/
        "tau_sf"            : (f"{corrdir}/TAU/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12/tau.json.gz", "v1"), # TEC and ID SF from POG (VsJet Medium WP only)
        "gen_tau_sf"        : (f"{corrdir}/TAU/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/tau_sf_pt-dm_DeepTau2018v2p5VSjet_{year}_{postfix}.json.gz", "v1"), # Tau ID SF (NEW FROM IC)
        "tes_sf"            : (f"{corrdir}/TAU/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/tau_es_dm_DeepTau2018v2p5_{year}_{postfix}.json.gz", "v1"), # Tau ID SF (NEW FROM IC)
        "tau_trig_sf"       : (f"{corrdir}/TAU/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/tau_trigger_DeepTau2018v2p5_{year}_{postfix}.json.gz", "v1"), # Tau ID SF (NEW FROM IC)
        # .. https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/ditaujet_jetleg_SFs_preEE.json?ref_type=heads
        "ditau_jet_trig_sf" : (f"{corrdir}/TAU/Run{run}_{year}{year_postfix}_Summer{year2}_NanoAODv12_Extra/ditaujet_jetleg_SFs_{postfix}.json",   "v1"),
        # https://gitlab.cern.ch/dwinterb/HiggsDNA/-/tree/master/higgs_dna/systematics/ditau/ROOT/Zpt?ref_type=heads
        "zpt_rewt_v1_sf"    : (f"{corrdir}/{year}{postfix}/Zpt/myZptCorrections.json.gz",   "v1"), # Zpt Rewt
        # https://indico.cern.ch/event/1360909/contributions/6000616/attachments/2875911/5036473/HLepRare_24.06.12.pdf
        # https://indico.cern.ch/event/489921/contributions/2000259/attachments/1248156/1839106/Recoil_20160323.pdf
        # /afs/cern.ch/user/d/dmroy/public/DY_pTll_recoil_corrections.json.gz
        #"zpt_rewt_v2_sf"    : (f"{external_path_parent}/Run3/Zpt/DY_pTll_recoil_corrections.json.gz",   "v1"), # Zpt Rewt
        "zpt_rewt_v2_sf"    : (f"{corrdir}/Run3/Zpt/DY_pTll_weights_v3.json.gz",                        "v1"), # Zpt Rewt        
        #"tautau_ff"         : (f"{external_path}/Fake_tautau/fake_factor_{year}_{postfix}.json",                       "v1"),
        #"tautau_ff"         : (f"{external_path}/Fake_tautau/fake_factor_{year}_{postfix}.json",                       "v1"),
        #"tautau_ff"         : (f"{external_path_parent}/FF_TauTau_combined/fake_factor_20222023.json",                 "v1"),
        #"tautau_ff"         : (f"{external_path_parent}/Run3/fake_factor_2023_HPS_DM_0_to_10.json",                    "v1"),
        "tautau_ff"         : (f"{corrdir}/Run3/FFtautauJSONs/fake_factor_20222023.json",                 "v1"),
        "tautau_ff_closure" : (f"{corrdir}/Run3/FFtautauJSONs/closure_correction_2023_postBPix_merged.json",   "v1"),
        #"tautau_ff0"        : (f"{external_path}/Fake_tautau/fake_factor_{year}_{postfix}_0cat_v2.json",               "v1"),
        #"tautau_ext_corr"   : (f"{external_path}/Fake_tautau/extrapolation_correction_inclusive.json",                 "v1"),
        #"btag_sf_corr": (f"{json_mirror}/POG/BTV/{year}_Summer{year2}{year_postfix}/btagging.json.gz",                "v1"),
        #"met_phi_corr": (f"{json_mirror}/POG/JME/2018_UL/met.json.gz",                                                "v1"), #met phi, unavailable Run3
        #"met_recoil"        : (f"{external_path_parent}/Run3/Recoil_corrections.json.gz",                              "v1"),
        "met_recoil"        : (f"{corrdir}/Run3/MET_Recoil/Recoil_corrections_v3.json.gz",                "v1"),
        "model_et_EVEN"     : (f"{corrdir}/Run3/ClassifierModels/Model_IC_etau/model_EVEN.json",            ""),
        "model_et_ODD"      : (f"{corrdir}/Run3/ClassifierModels/Model_IC_etau/model_ODD.json",             ""),        
        "model_mt_EVEN"     : (f"{corrdir}/Run3/ClassifierModels/Model_IC_mutau/model_EVEN.json",           ""),
        "model_mt_ODD"      : (f"{corrdir}/Run3/ClassifierModels/Model_IC_mutau/model_ODD.json",            ""),        
        "model_tt_EVEN"     : (f"{corrdir}/Run3/ClassifierModels/Model_IC_tautau/model_EVEN.json",          ""),
        "model_tt_ODD"      : (f"{corrdir}/Run3/ClassifierModels/Model_IC_tautau/model_ODD.json",           ""),
    })
    """





    
    
    # ------------------------------------------------------ #
    # electron settings                                      #
    # names of electron correction sets and working points   #
    # (used in the electron_sf producer)                     #
    # ------------------------------------------------------ #

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3

    # names of electron correction sets and working points
    from columnflow.production.cms.electron import ElectronSFConfig
    from columnflow.calibration.cms.egamma import EGammaCorrectionConfig
    if run == 2:
        # SFs
        e_postfix = ""
        if year == 2016:
            e_postfix = {"APV": "preVFP", "": "postVFP"}[year_postfix]
        #cfg.x.electron_sf_names = ElectronSFConfig(
        #    correction="UL-Electron-ID-SF",
        #    campaign=f"{year}{e_postfix}",
        #    working_point="wp80iso",
        #)
        cfg.x.electron_id_sf = ElectronSFConfig(
            correction="UL-Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point="wp80iso",
        )
        cfg.x.electron_reco_sf = ElectronSFConfig(
            correction="UL-Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point={
                "RecoBelow20": (lambda variables: variables["pt"] < 20.0),
                "RecoAbove20": (lambda variables: variables["pt"] >= 20.0),
            },
        )
        # TODO: disabled for now, need to adapt once new run 2 processing is done
        # # eec and eer
        # cfg.x.eec = EGammaCorrectionConfig(
        #     correction_set="Scale",
        #     value_type="total_correction",
        #     uncertainty_type="total_uncertainty",
        # )
        # cfg.x.eer = EGammaCorrectionConfig(
        #     correction_set="Smearing",
        #     compound=False,
        #     value_type="rho",
        #     uncertainty_type="err_rho",
        # )
        # electron scale and smearing (eec and eer)
        logger.warning("For ESS, egm tool is disabled because Electron.seediEtaOriX & Electron.seediPhiOriY are only available in NanoAODv15")
        cfg.x.ess = EGammaCorrectionConfig(
            scale_correction_set="Scale",
            scale_compound=True,
            smear_syst_correction_set="SmearAndSyst",
            systs=["scale_down", "scale_up", "smear_down", "smear_up"],
            use_egm_tool=True,
        )
    elif run == 3:
        # SFs
        if year == 2022:
            e_postfix = {"": "Re-recoBCD", "EE": "Re-recoE+PromptFG"}[campaign.x.postfix]
        elif year == 2023:
            e_postfix = {"": "PromptC", "BPix": "PromptD"}[campaign.x.postfix]
        elif year == 2024:
            e_postfix = "Prompt"
        else:
            assert False
        cfg.x.electron_id_sf = ElectronSFConfig(
            correction="Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point="wp80iso",
        )
        cfg.x.electron_reco_sf = ElectronSFConfig(
            correction="Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point={
                "RecoBelow20": (lambda variables: variables["pt"] < 20.0),
                "Reco20to75": (lambda variables: (variables["pt"] >= 20.0) & (variables["pt"] < 75.0)),
                "RecoAbove75": (lambda variables: variables["pt"] >= 75.0),
            },
        )
        cfg.x.electron_trigger_sf_names = ElectronSFConfig(
            correction="Electron-HLT-SF",
            campaign=f"{year}{e_postfix}",
            hlt_path="HLT_SF_Ele30_TightID",
        )
        cfg.x.single_trigger_electron_data_effs_cfg = ElectronSFConfig(
            correction="Electron-HLT-DataEff",
            campaign=f"{year}{e_postfix}",
            hlt_path="HLT_SF_Ele30_TightID",
        )
        cfg.x.single_trigger_electron_mc_effs_cfg = ElectronSFConfig(
            correction="Electron-HLT-McEff",
            campaign=f"{year}{e_postfix}",
            hlt_path="HLT_SF_Ele30_TightID",
        )
        cfg.x.cross_trigger_electron_data_effs_cfg = ElectronSFConfig(
            correction="Electron-HLT-DataEff",
            campaign=f"{year}{e_postfix}",
            hlt_path="HLT_SF_Ele24_TightID",
        )
        cfg.x.cross_trigger_electron_mc_effs_cfg = ElectronSFConfig(
            correction="Electron-HLT-McEff",
            campaign=f"{year}{e_postfix}",
            hlt_path="HLT_SF_Ele24_TightID",
        )
        # electron scale and smearing (eec and eer)
        cfg.x.ess = EGammaCorrectionConfig(
            scale_correction_set="Scale",
            scale_compound=True,
            smear_syst_correction_set="SmearAndSyst",
            systs=["scale_down", "scale_up", "smear_down", "smear_up"],
            use_egm_tool=True,
        )
    else:
        assert False

    """
    electron_sf_tag = ""
    if year == 2022:
        electron_sf_tag = "2022Re-recoE+PromptFG" if year_postfix else "2022Re-recoBCD"
    elif year == 2023:
        electron_sf_tag = "2023PromptD" if year_postfix else "2023PromptC"
    elif year == 2024:
        raise RuntimeWarning("too early")
    else:
        raise RuntimeError("wrong year")
    
    cfg.x.electron_sf_names = (
        "Electron-ID-SF",
        electron_sf_tag,
        "wp80iso",
    )
    cfg.x.electron_trig_sf_names = (
        "Electron-HLT-SF",
        electron_sf_tag,
        "HLT_SF_Ele30_TightID",
    )
    cfg.x.electron_xtrig_sf_names = (
        "Electron-HLT-SF",
        electron_sf_tag,
        "HLT_SF_Ele24_TightID",
    )
    """

    # -------------------------------------------------- #
    # tau settings                                       #
    # names of channels, tau ID working points, systs    #
    # (used in the zttpol tau producer)                  #
    # -------------------------------------------------- #
    from zttpol.production.tau import TauSFConfig
    if run == 2:
        cfg.x.tau_id_sf = TauSFConfig(
            tagger_vsJet="DeepTau2018v2p5",
            correction_vsJet="",
            syst_vsJet=[(i, f"syst_TES_{year}_dm{i}") for i in [0,1,10,11]],
            tagger_vsEle="DeepTau2018v2p5",
            correction_vsEle="",
            syst_vsEle=[], # default : up/down
            tagger_vsMu="DeepTau2018v2p5",
            correction_vsMu=[],
        )

    
    # -------------------------------------------------- #
    # muon settings                                      #
    # names of muon correction sets and working points   #
    # (used in the muon producer)                        #
    # -------------------------------------------------- #

    # names of muon correction sets and working points
    # (used in the muon producer)
    from columnflow.production.cms.muon import MuonSFConfig
    if run == 2:
        cfg.x.muon_id_sf = MuonSFConfig(correction="NUM_MediumID_DEN_TrackerMuons")
        cfg.x.muon_iso_sf = MuonSFConfig(correction="NUM_TightRelIso_DEN_MediumID")
        cfg.x.muon_single_trigger_sf = MuonSFConfig(correction="NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight")
    elif run == 3:
        # id and iso
        cfg.x.muon_id_sf = MuonSFConfig(correction="NUM_TightID_DEN_TrackerMuons", min_pt=15.0)
        cfg.x.muon_id_sf_lowpt = MuonSFConfig(correction="NUM_TightID_DEN_TrackerMuons")  # producer uses min_pt above
        cfg.x.muon_iso_sf = MuonSFConfig(correction="NUM_TightPFIso_DEN_TightID", min_pt=15.0)

        # trigger
        cfg.x.muon_trigger_sf_names = MuonSFConfig(
            correction="NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        )

        cfg.x.single_trigger_muon_data_effs_cfg = MuonSFConfig(
            correction="NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        )
        cfg.x.single_trigger_muon_mc_effs_cfg = MuonSFConfig(
            correction="NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        )

        cfg.x.cross_trigger_muon_data_effs_cfg = MuonSFConfig(
            correction="NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_DATAeff",
        )
        cfg.x.cross_trigger_muon_mc_effs_cfg = MuonSFConfig(
            correction="NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_MCeff",
        )

        # mec/mer
        from columnflow.calibration.cms.muon import MuonSRConfig
        cfg.x.muon_sr = MuonSRConfig(
            systs=["scale_up", "scale_down", "res_up", "res_down"],
        )
    else:
        assert False

    """
    cfg.x.muon_id_sf_names  = (
        "NUM_MediumID_DEN_TrackerMuons",
        f"{year}_{postfix}"
    )
    cfg.x.muon_iso_sf_names = (
        "NUM_TightPFIso_DEN_MediumID",
        f"{year}_{postfix}"
    )
    cfg.x.muon_IsoMu24_trigger_sf_names = (
        "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium",
        f"{year}_{postfix}"
    )
    cfg.x.muon_xtrig_sf_names = (
        "NUM_IsoMu20_DEN_CutBasedIdMedium_and_PFIsoMedium",
        f"{year}_{postfix}"
    )
    """

    
    ################################################################################################
    # dataset / process specific methods
    ################################################################################################

    """
    # top pt reweighting
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting?rev=31
    from columnflow.production.cms.top_pt_weight import TopPtWeightConfig
    cfg.x.top_pt_weight = TopPtWeightConfig(
        params={
            "a": 0.0615,
            "a_up": 0.0615 * 1.5,
            "a_down": 0.0615 * 0.5,
            "b": -0.0005,
            "b_up": -0.0005 * 1.5,
            "b_down": -0.0005 * 0.5,
        },
        pt_max=500.0,
    )

    cfg.x.top_pt_reweighting_params = {
        "a": 0.0615,
        "a_up": 0.0615 * 1.5,
        "a_down": 0.0615 * 0.5,
        "b": -0.0005,
        "b_up": -0.0005 * 1.5,
        "b_down": -0.0005 * 0.5,
    }
    """
    # https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/master/higgs_dna/systematics/ditau/event_weight_systematics.py#L124-125
    cfg.x.top_pt_reweighting_params = {
        "a": 0.103,
        "a_up": 0.103 * 1.5,
        "a_down": 0.103 * 0.5,
        "b": -0.0018,
        "b_up": -0.0018 * 1.5,
        "b_down": -0.0018 * 0.5,
        "c": -0.000134,
        "c_up": -0.000134 * 1.5,
        "c_down": -0.000134 * 0.5,
        "d": 0.973,
        "d_up": 0.973 * 1.5,
        "d_down": 0.973 * 0.5,
        "e": 0.991,
        "e_up": 0.991 * 1.5,
        "e_down": 0.991 * 0.5,
        "f": 0.000075,
        "f_up": 0.000075 * 1.5,
        "f_down": 0.000075 * 0.5,   
    }


    
    
    # --------------------------------------------------------------------------------------------- #
    # register shifts
    # --------------------------------------------------------------------------------------------- #

    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})

    cfg.add_shift(name="hdamp_up", id=3, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="hdamp_down", id=4, type="shape", tags={"disjoint_from_nominal"})

    cfg.add_shift(name="mtop_up", id=5, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="mtop_down", id=6, type="shape", tags={"disjoint_from_nominal"})

    # --- >>> PU weight <<< --- #
    cfg.add_shift(name="minbias_xs_up", id=7, type="shape")
    cfg.add_shift(name="minbias_xs_down", id=8, type="shape")
    add_shift_aliases(
        cfg,
        "minbias_xs",
        {
            "pu_weight": "pu_weight_{name}",
            "normalized_pu_weight": "normalized_pu_weight_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        },
    )

    cfg.add_shift(name="top_pt_up", id=9, type="shape")
    cfg.add_shift(name="top_pt_down", id=10, type="shape")
    add_shift_aliases(
        cfg,
        "top_pt",
        {
            "top_pt_weight": "top_pt_weight_{direction}",
        },
    )

    for i, (jec_source, flag) in enumerate(all_jec_sources.items()):
        if not flag:
            continue
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5001 + 2 * i,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5002 + 2 * i,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        add_shift_aliases(
            cfg,
            f"jec_{jec_source}",
            {
                "Jet.pt": "Jet.pt_{name}",
                "Jet.mass": "Jet.mass_{name}",
                f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )

    cfg.add_shift(name="jer_up", id=6001, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6002, type="shape", tags={"jer"})
    add_shift_aliases(
        cfg,
        "jer",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        },
    )

    for i, (match, dm) in enumerate(itertools.product(["tau", "e", "mu"], [0, 1, 10, 11])):
        cfg.add_shift(name=f"tec_{match}_dm{dm}_up", id=20 + 2 * i, type="shape", tags={"tec"})
        cfg.add_shift(name=f"tec_{match}_dm{dm}_down", id=21 + 2 * i, type="shape", tags={"tec"})
        add_shift_aliases(
            cfg,
            f"tec_{match}_dm{dm}",
            {
                "Tau.pt": "Tau.pt_{name}",
                "Tau.mass": "Tau.mass_{name}",
                # no MET propagation needed for tec shifts
                # f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                # f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )
    
    # --- >>> tau weight <<< --- #
    if channel not in ['ee','mumu']:
        tau_uncerts = [f'tau_dm{i}' for i in [0,1,10,11]]
        e_uncerts = [f'e_dm{i}' for i in [0,1,10,11]]
        m_uncerts = [f'mu_{region}' for region in ['0p0To0p4', '0p4To0p8', '0p8To1p2', '1p2To1p7', '1p7To2p3']]
        cfg.x.tau_unc_names = {
            'etau': e_uncerts + tau_uncerts,
            'mutau': m_uncerts + tau_uncerts,
            'tautau': tau_uncerts, 
        }[channel]
        for i, unc in enumerate(cfg.x.tau_unc_names):
            cfg.add_shift(name=f"tau_{unc}_up", id=51 + 2 * i, type="shape")
            cfg.add_shift(name=f"tau_{unc}_down", id=52 + 2 * i, type="shape")
            add_shift_aliases(cfg, f"tau_{unc}", {"tau_id_weight": f"tau_id_weight_{unc}_{{direction}}"})
    
    # --- >>> tau weight <<< --- #
    #tau_uncerts = [f'tau_dm{i}' for i in [0,1,10,11]]
    #e_uncerts = [f'e_dm{i}' for i in [0,1,10,11]]
    #m_uncerts = [f'mu_{region}' for region in ['0p0To0p4', '0p4To0p8', '0p8To1p2', '1p2To1p7', '1p7To2p3']]
    #cfg.x.tau_unc_names = tau_uncerts + e_uncerts + m_uncerts
    #for i, unc in enumerate(cfg.x.tau_unc_names):
    #    cfg.add_shift(name=f"tau_{unc}_up", id=51 + 2 * i, type="shape")
    #    cfg.add_shift(name=f"tau_{unc}_down", id=52 + 2 * i, type="shape")
    #    add_shift_aliases(cfg,f"tau_{unc}", {"tau_id_weight": f"tau_id_weight_{unc}_{{direction}}"})

        
    #cfg.add_shift(name=f"tau_trig_up", id=91, type="shape")
    #cfg.add_shift(name=f"tau_trig_down", id=92, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "tau_trig",
    #    {
    #        "tau_trigger_weight": "tau_trigger_weight_{direction}",
    #    },
    #)

    # --- >>> ele weight <<< --- #
    cfg.add_shift(name="e_id_up", id=101, type="shape")
    cfg.add_shift(name="e_id_down", id=102, type="shape")
    add_shift_aliases(cfg,"e_id", {"electron_id_weight": "electron_id_weight_{direction}"})
    
    cfg.add_shift(name="e_reco_up", id=103, type="shape")
    cfg.add_shift(name="e_reco_down", id=104, type="shape")
    add_shift_aliases(cfg, "e_reco",{"electron_reco_weight": "electron_reco_weight_{direction}"})
    
    # electron scale and smearing
    cfg.add_shift(name="eec_up", id=105, type="shape", tags={"eec"})
    cfg.add_shift(name="eec_down", id=106, type="shape", tags={"eec"})
    add_shift_aliases(cfg, "eec", {"Electron.pt": "Electron.pt_scale_{direction}"})

    cfg.add_shift(name="eer_up", id=107, type="shape", tags={"eer"})
    cfg.add_shift(name="eer_down", id=108, type="shape", tags={"eer"})
    add_shift_aliases(cfg, "eer", {"Electron.pt": "Electron.pt_smear_{direction}"})

      
    #cfg.add_shift(name="e_trig_up", id=92, type="shape")
    #cfg.add_shift(name="e_trig_down", id=93, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "e_trig",
    #    {
    #        "electron_Ele30_WPTight_trigger_weight": "electron_Ele30_WPTight_trigger_weight_{direction}",
    #    },
    #)
    #cfg.add_shift(name="e_xtrig_up", id=94, type="shape")
    #cfg.add_shift(name="e_xtrig_down", id=95, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "e_xtrig",
    #    {
    #        "electron_xtrig_weight": "electron_xtrig_weight_{direction}",
    #    },
    #)

    # --- >>> mu weight <<< --- #
    cfg.add_shift(name="mu_id_up", id=121, type="shape")
    cfg.add_shift(name="mu_id_down", id=122, type="shape")
    add_shift_aliases(cfg, "mu_id", {"muon_id_weight": "muon_id_weight_{direction}"})
    
    cfg.add_shift(name="mu_iso_up", id=123, type="shape")
    cfg.add_shift(name="mu_iso_down", id=124, type="shape")
    add_shift_aliases(cfg,"mu_iso", {"muon_iso_weight": "muon_iso_weight_{direction}"})

    # muon scale and resolution
    cfg.add_shift(name="mec_up", id=125, type="shape", tags={"mec"})
    cfg.add_shift(name="mec_down", id=126, type="shape", tags={"mec"})
    add_shift_aliases(cfg, "mec", {"Muon.pt": "Muon.pt_scale_{direction}"})
    
    cfg.add_shift(name="mer_up", id=127, type="shape", tags={"mer"})
    cfg.add_shift(name="mer_down", id=128, type="shape", tags={"mer"})
    add_shift_aliases(cfg, "mer", {"Muon.pt": "Muon.pt_res_{direction}"})
    
    
    #cfg.add_shift(name="mu_trig_up", id=104, type="shape")
    #cfg.add_shift(name="mu_trig_down", id=105, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "mu_trig",
    #    {
    #        "muon_single_trigger_weight": "muon_single_trigger_weight_{direction}",
    #    },
    #)
    #cfg.add_shift(name="mu_xtrig_up", id=106, type="shape")
    #cfg.add_shift(name="mu_xtrig_down", id=107, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "mu_xtrig",
    #    {
    #        "muon_xtrig_weight": "muon_xtrig_weight_{direction}",
    #    },
    #)

    # --- >>> pdf weight <<< --- #
    cfg.add_shift(name="pdf_up", id=141, type="shape")
    cfg.add_shift(name="pdf_down", id=142, type="shape")
    add_shift_aliases(cfg,"pdf", {
        "pdf_weight": "pdf_weight_{direction}",
        "normalized_pdf_weight": "normalized_pdf_weight_{direction}",
    })
    # --- >>> renormalization and factorization weight <<< --- #
    cfg.add_shift(name="murmuf_up", id=151, type="shape", tags={"lhe_weight"})
    cfg.add_shift(name="murmuf_down", id=152, type="shape", tags={"lhe_weight"})
    add_shift_aliases(cfg,"murmuf", {
        "murmuf_weight": "murmuf_weight_{direction}",
        "normalized_murmuf_weight": "normalized_murmuf_weight_{direction}",
    })
    # --- >>> parton-shower (isr) weight <<< --- #
    cfg.add_shift(name="isr_up", id=161, type="shape")
    cfg.add_shift(name="isr_down", id=162, type="shape")
    add_shift_aliases(cfg,"isr", {
        "isr_weight": "isr_weight_{direction}",
        "normalized_isr_weight": "normalized_isr_weight_{direction}",
    })
    # --- >>> parton-shower (fsr) weight <<< --- #    
    cfg.add_shift(name="fsr_up", id=171, type="shape")
    cfg.add_shift(name="fsr_down", id=172, type="shape")
    add_shift_aliases(cfg,"fsr",{
        "fsr_weight": "fsr_weight_{direction}",
        "normalized_fsr_weight": "normalized_fsr_weight_{direction}",
    })

    # --- >>> zpt weight <<< --- #    
    cfg.add_shift(name="zpt_up", id=181, type="shape")
    cfg.add_shift(name="zpt_down", id=182, type="shape")
    add_shift_aliases(cfg,"zpt",{
        "zpt_reweight": "zpt_reweight_{direction}",
    })


    # --- >>> fake factor weight <<< --- # D U M M Y !!!!!!!!   
    cfg.add_shift(name="ff_up", id=191, type="shape")
    cfg.add_shift(name="ff_down", id=192, type="shape")
    add_shift_aliases(cfg,"ff", {"ff_weight": "ff_weight_{direction}"})
    
    cfg.add_shift(name="ff_cls_corr_up", id=193, type="shape")
    cfg.add_shift(name="ff_cls_corr_down", id=194, type="shape")
    add_shift_aliases(cfg,"ff_cls_corr", {"ff_cls_corr_weight": "ff_cls_corr_weight_{direction}"})

    cfg.add_shift(name="ff_ext_corr_up", id=195, type="shape")
    cfg.add_shift(name="ff_ext_corr_down", id=196, type="shape")
    add_shift_aliases(cfg,"ff_ext_corr", {"ff_ext_corr_weight": "ff_ext_corr_weight_{direction}"})



    # --- >>> tau weight <<< --- #
    cfg.x.trigger_legs = {
        'ee'    : ['e'],
        'etau'  : ['e','tau_dm0','tau_dm1','tau_dm10','tau_dm11'],
        'emu'   : ['emu'],
        'mumu'  : ['mu'],
        'mutau' : ['mu','tau_dm0','tau_dm1','tau_dm10','tau_dm11'],
        'tautau': ['tau_dm0','tau_dm1','tau_dm10','tau_dm11','jet'] if run == 3 else ['tau_dm0','tau_dm1','tau_dm10','tau_dm11'],
    }[channel]
    for i, leg in enumerate(cfg.x.trigger_legs):
        cfg.add_shift(name=f"trigger_{leg}_up", id=201 + 2 * i, type="shape")
        cfg.add_shift(name=f"trigger_{leg}_down", id=202 + 2 * i, type="shape")
        add_shift_aliases(cfg, f"trigger_{leg}", {"trigger_weight": f"trigger_weight_{leg}_{{direction}}"})

    

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    

    #---------------------------------------------------------------------------------------------#
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    # get_shifts = functools.partial(get_shifts_from_sources, cfg)
    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    #---------------------------------------------------------------------------------------------#

    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight"           : [],
        "normalized_pu_weight"           : get_shifts("minbias_xs"),
        "trigger_weight"                 : get_shifts(*(f"trigger_{leg}" for leg in cfg.x.trigger_legs)),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })
    if channel in {'mumu','emu','mutau'}:
        # Muon related corrections
        cfg.x.event_weights["muon_id_weight"] = get_shifts("mu_id")
        cfg.x.event_weights["muon_iso_weight"] = get_shifts("mu_iso")
    
    if channel in {'ee','emu','etau'}:
        cfg.x.event_weights["electron_id_weight"] = get_shifts("e_id")
        cfg.x.event_weights["electron_reco_weight"] = get_shifts("e_reco")
        
    if channel in {'etau','mutau','tautau'}:
        cfg.x.event_weights["tau_id_weight"] = get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names))
        #cfg.x.event_weights["tau_trigger_weight"] = get_shifts("tau_trig")

    """
    cfg.x.event_mumu_weights = DotDict({
        "normalization_weight"                  : [],
        "normalized_pu_weight"                  : get_shifts("minbias_xs"),
        "muon_id_weight"                        : get_shifts("mu_id"),
        "muon_iso_weight"                       : get_shifts("mu_iso"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })
    cfg.x.event_emu_weights = DotDict({
        "normalization_weight"                  : [],
        "pu_weight"                             : get_shifts("minbias_xs"),
        "electron_idiso_weight"                 : get_shifts("e"),
        "electron_Ele30_WPTight_trigger_weight" : get_shifts("e_trig"),
        "electron_xtrig_weight"                 : get_shifts("e_xtrig"),
        "muon_id_weight"                        : get_shifts("mu_id"),
        "muon_iso_weight"                       : get_shifts("mu_iso"),
        "muon_IsoMu24_trigger_weight"           : get_shifts("mu_trig"),
        "muon_xtrig_weight"                     : get_shifts("mu_xtrig"),
        #"ff_weight"                             : [],
        #"ff_cls_corr_weight"                    : [],
        ##"ff_ext_corr_weight"                    : [],
        #"tauspinner_weight"                     : get_shifts("tauspinner"),
        #"pdf_weight"                            : [],
        "zpt_reweight"                          : get_shifts("zpt"),
        "top_pt_weight"                         : [],
    })
    cfg.x.event_etau_weights = DotDict({
        "normalization_weight"                  : [],
        "normalized_pu_weight"                  : get_shifts("minbias_xs"),
        "electron_id_weight"                    : get_shifts("e_id"),
        "electron_reco_weight"                  : get_shifts("e_reco"),
        #"electron_Ele30_WPTight_trigger_weight" : get_shifts("e_trig"),
        #"electron_xtrig_weight"                 : get_shifts("e_xtrig"),
        "tau_id_weight"                         : get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })
    cfg.x.event_mutau_weights = DotDict({
        "normalization_weight"                  : [],
        "normalized_pu_weight"                  : get_shifts("minbias_xs"),
        "muon_id_weight"                        : get_shifts("mu_id"),
        "muon_iso_weight"                       : get_shifts("mu_iso"),
        #"muon_IsoMu24_trigger_weight"           : get_shifts("mu_trig"),
        #"muon_xtrig_weight"                     : get_shifts("mu_xtrig"),
        "tau_id_weight"                         : get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })
    cfg.x.event_tautau_weights = DotDict({
        "normalization_weight"                  : [],
        "normalized_pu_weight"                  : get_shifts("minbias_xs"),
        "tau_id_weight"                         : get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })
    """
    #---------------------------------------------------------------------------------------------#
    # No Idea
    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    #---------------------------------------------------------------------------------------------#
    
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev1'

    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }

    
    # ------------------------------------------------------------- #
    #  LFN settings                                                 #
    #  Campaigns must have creator names as desy or IPHC            #
    #  the main path is in the campaign __init__, the basepath      #
    #  is mentioned in the dataset info                             #
    # ------------------------------------------------------------- #

    #redirector_name = "wlcg_fs_eoscms_redirector"
    #redirector_name = "wlcg_fs_imperial_redirector"
    #redirector_name = "wlcg_fs_global_redirector"
    local_redirector_name = "wlcg_fs_eoscms_redirector"
    global_redirector_name = "wlcg_fs_iphc_redirector"
    
    campaign_tag = cfg.campaign.x("custom").get("creator")
    if campaign_tag == "IPHC":
        def get_dataset_lfns(dataset_inst: od.Dataset,
                             shift_inst: od.Shift,
                             dataset_key: str) -> list[str]:

            # destructure dataset_key into parts and create the lfn base directory
            logger.info(f"Creating custom get_dataset_lfns for {config_name}")

            if not dataset_inst.name.startswith("dy_2tau_m50"):

                store_path = CMSDatasetInfo.from_key(dataset_key).store_path.lstrip("/")            
                logger.info(f'path : {store_path}')
            
                dir_cls = law.wlcg.WLCGDirectoryTarget
                lfn_base = dir_cls(store_path, fs=global_redirector_name)
                lfn_num_bases = [lfn_base.child(d, type="d") for d in lfn_base.listdir() if d.isnumeric()]
            
                lfns = sum((
                    [
                        "/" + lfn_num_base.child(basename, type="f").path.lstrip("/")
                        for basename in lfn_num_base.listdir(pattern="*.root")
                    ]
                    for lfn_num_base in lfn_num_bases
                ), [])
                

            else:
                logger.warning(f"For {dataset_inst.name} dataset, use files stored locally in {local_file_path}")
                basepath = local_file_path
                logger.info(f"Location : {basepath}")
                lfn_base = law.wlcg.WLCGDirectoryTarget(
                    f"{basepath}{dataset_key}",
                    fs=local_redirector_name,
                )
                logger.info(f"lfn basedir:{lfn_base}")
                # loop though files and interpret paths as lfns
                lfns = [
                    lfn_base.child(basename, type="f").path
                    for basename in lfn_base.listdir(pattern="*.root")
                ]
                
            return sorted(lfns)

        
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        #cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: redirector_name
        #cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: [global_redirector_name]
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: [local_redirector_name,
                                                                  global_redirector_name]
        
    #---------------------------------------------------------------------------------------------#
    # Add categories described in categorization.py
    #---------------------------------------------------------------------------------------------#

    from importlib import import_module

    module = import_module(f"zttpol.config.categories_{channel}")
    module.add_categories(cfg)

    
    #cfg.x.ff_apply_id_map = DotDict.wrap({
    #    "etau"  : {},
    #    "mutau" : {},
    #    "tautau" : {
    #        "id_for_B"  : [cfg.get_category("tautau").id, cfg.get_category("real_1").id, cfg.get_category("hadB").id],
    #        "id_for_C"  : [cfg.get_category("tautau").id, cfg.get_category("real_1").id, cfg.get_category("hadC").id],  # category_id for AR C 
    #        "id_for_C0" : [cfg.get_category("tautau").id, cfg.get_category("real_1").id, cfg.get_category("hadC0").id], # category_id for AR C0
    #    },
    #})
    
    #---------------------------------------------------------------------------------------------#
    # Add variables described in variables.py
    #---------------------------------------------------------------------------------------------#
    
    from zttpol.config.variables import add_variables
    add_variables(cfg)

    #---------------------------------------------------------------------------------------------#
    # columns to keep after reduce events, MergeSelectionMasks and UniteColumns tasks
    #---------------------------------------------------------------------------------------------#

    
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            ## mandatory
            ColumnCollection.MANDATORY_COFFEA,
            "LHEPart.*",
            ## TauProds
            "TauProd.*",
            ## PV/SV
            "PV.{x,y,z,npvs,npvsGood}",
            "PVBS.*",
            "SV.{x,y,z}",
            ## Gen
            "GenPart.{pt,eta,phi,mass,status,pdgId,statusFlags}",
            "GenZ.{pt,eta,phi,mass}",
            "GenZvis.{pt,eta,phi,mass}",
            ## MET
            f"{cfg.x.met_name}.{{pt,phi,significance,covXX,covXY,covYY}}",
            # variations created during met phi calibration and that are not registered shifts to the selector
            f"{cfg.x.met_name}.{{pt,phi}}_{{unsmeared,metphi_*,minbias_xs_*}}",
            ## Jet
            "Jet.{pt,eta,phi,mass,hadronFlavour,puId,btagDeepFlavB,chHEF,neHEF,chEmEF,neEmEF,muEF,chMultiplicity,neMultiplicity}",
            "bJet.{pt,eta,phi,mass,btagDeepFlavB,hadronFlavour}",
            "trigJet.{pt,eta,phi,mass}",
            "metRecoilJet.{pt,eta,phi,mass}",
            ## Tau
            "RawTau.{pt,eta,phi,mass,decayMode,rawIdx}",
            "PreSelTau.{pt,eta,phi,mass,decayMode,rawIdx}",
            "Tau.*",
            ## Muon
            "RawMuon.{pt,eta,phi,mass,dxy,decayMode,rawIdx}",
            "PreSelMuon.{pt,eta,phi,mass,dxy,decayMode,rawIdx}",
            "Muon.*",
            ## Electron
            "RawElectron.{pt,eta,phi,mass,dxy,decayMode,rawIdx}",
            "PreSelElectron.{pt,eta,phi,mass,dxy,decayMode,rawIdx}",
            "Electron.*",    
            # keep all columns added during selection and reduction, but skip cutflow features
            ColumnCollection.ALL_FROM_SELECTOR,
            skip_column("cutflow.*"),
        },
        "cf.MergeSelectionMasks": {
            "cutflow.*",
        },
        "cf.UniteColumns": {
            # all columns except for shifts
            "all": {
                "*",
                *skip_column("*_{up,down}"),
            },
        },
    })

    # --------------------------------------------------------------------------------------------- #
    # Adding hist hooks
    # --------------------------------------------------------------------------------------------- #

    #cfg.x.regions_to_extrapolate_fake = "AB" # "AB" or "CD" or "C0D0"
    #cfg.x.save_qcd = True
    #from zttpol.config.hist_hooks import add_hist_hooks
    #add_hist_hooks(cfg)

    # fastMTT helper
    cfg.x.enable_fastMTT = False
    cfg.x.enable_fastMTT_for_phiCP = False # PV only + should be False automatically if not cfg.x.enable_fastMTT
    
    #---------------------------------------------------------------------------------------------#
    # Helper switch for debugging
    #---------------------------------------------------------------------------------------------#

    cfg.x.verbose = DotDict.wrap({
        "calibration": {
            "main"                    : False,
            "tau"                     : False,
        },
        "selection": {
            "main"                    : False,
            "trigobject_matching"     : False,
            "extra_lep_veto"          : False,
            "dilep_veto"              : False,
            "higgscand"               : False,
        },
        "production": {
            "main"                    : False,
        },
    })


    cfg.x.extra_tags = DotDict.wrap({
        "genmatch"       : False,
    })


    cfg.x.is_channel_specific = False
    cfg.x.channel_specific_info = DotDict.wrap({
        "emu"    : False,
        "etau"   : False,
        "mutau"  : False,
        "tautau" : True,
        "ee"     : False,
        "mumu"   : False,
    })
