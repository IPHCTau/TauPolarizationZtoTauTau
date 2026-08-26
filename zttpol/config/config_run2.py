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
from columnflow.columnar_util import ColumnCollection, EMPTY_FLOAT
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
thisdir = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/TauPolarizationZtoTauTau/zttpol/config"
tooldir = os.path.join(os.path.dirname(thisdir), "data/tools")
corrdir = os.path.join(os.path.dirname(thisdir), "data/corrections")

def add_config (ana: od.Analysis,
                campaign: od.Campaign,
                config_name           = None,
                config_id             = None,
                limit_dataset_files   = None,
                channel               = None,
) -> od.Config :

    #from IPython import embed; embed()
    logger.info(f"Channel to be analyzed: {channel}")
    if channel not in {'emu','etau','mutau','tautau'}:
        raise RuntimeError(f"WRONG CHANNEL : {channel}")
    
    # gather campaign data
    run = campaign.x.run
    year = campaign.x.year
    postfix = campaign.x.postfix
    
    # some validations
    assert run in {2}
    assert year in {2016,2017,2018}

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)

    
    # --------------------------------------------------------------------------------------------- #
    # create a config by passing the campaign, so id and name will be identical    
    # --------------------------------------------------------------------------------------------- #

    cfg = ana.add_config(campaign,
                         name  = config_name,
                         id    = config_id)

    # add channel
    cfg.x.channel = channel

    
    cfg.add_channel(name="emu",    id=1)
    cfg.add_channel(name="etau",   id=2)
    cfg.add_channel(name="mutau",  id=3)
    cfg.add_channel(name="tautau", id=4)

    
    # combination of processes
    cfg.add_process(
        name="multiboson",
        id=7999, # vv proc is 8000 and vvv is 9000 in cmsdb processes
        label="VV(V)",
        processes=[procs.n.vv, procs.n.vvv],
    )
    cfg.add_process(
        name="top",
        id=1999, 
        label="top",
        processes=[procs.n.tt, procs.n.st],
    )
    
    # --------------------------------------------------------------------------------------------- #
    # add processes we are interested in
    # --------------------------------------------------------------------------------------------- #
    process_names = [
        #"dy",
        "dy_2tau_m50toinf_amcatnloFXFX",
        "dy_2tau_m50toinf_amcatnloFXFX_LHEspin_minus",
        "dy_2tau_m50toinf_amcatnloFXFX_LHEspin_plus",
    ]

    for process_name in process_names:
        # development switch in case datasets are not _yet_ there
        #if process_name not in procs:
        #    logger.warning(f"WARNING: {process_name} not in cmsdb processes")
        #    continue
        # add the process
        #if process_name == "h_ggf_tautau":
        #    procs.get(process_name).is_signal = True
        #proc = cfg.add_process(procs.get(process_name))
        #print(procs.get(process_name))
        #if proc.name == "h_ggf_tautau":
        #    proc.is_signal = True
        if process_name == "qcd":
            # qcd is not part of procs since there is no dataset registered for it
            from cmsdb.processes.qcd import qcd
            proc = cfg.add_process(qcd)
        elif process_name not in procs:
            logger.warning(f"WARNING: {process_name} not in cmsdb processes")
            continue
        else:
            proc = cfg.add_process(procs.get(process_name))

        if re.match(r"^tt(|_.+)$", process_name):
            proc.add_tag({"ttbar", "tt"})

            
    # configuration of colors, labels, etc. can happen here
    from zttpol.config.styles import stylize_processes
    stylize_processes(cfg)    

    
    # --------------------------------------------------------------------------------------------- #
    # add datasets we need to study
    # --------------------------------------------------------------------------------------------- #

    dataset_names = [
        ##Drell-Yan
        "dy_2tau_m50_njets_amcatnlo",
    ]

    cfg.x.signal_filter_efficiency = DotDict.wrap({
        "h_ggf_tautau_uncorrelatedDecay_SM_Filtered_ProdAndDecay": 0.3847,
        "h_ggf_tautau_uncorrelatedDecay_MM_Filtered_ProdAndDecay": 0.3848,
        "h_ggf_tautau_uncorrelatedDecay_CPodd_Filtered_ProdAndDecay": 0.3848,
        "wph_tautau_uncorrelatedDecay_Filtered": 0.3743,
        "wmh_tautau_uncorrelatedDecay_Filtered": 0.3944,
        "zh_tautau_uncorrelatedDecay_Filtered": 0.3933,
        "h_vbf_tautau_uncorrelatedDecay_Filtered": 0.4091,
    })

    
    datasets_data = []
    if year == 2018:
        datasets_data = ["data_e_C",   "data_e_D",
                         "data_single_mu_C",
                         "data_mu_C",  "data_mu_D", 
                         "data_tau_C", "data_tau_D"]
    else:
        raise RuntimeError(f"Check Year in __init__.py in cmsdb campaign: {year}")


    dataset_names = datasets_data + dataset_names
    for dataset_name in dataset_names:
        # development switch in case datasets are not _yet_ there
        if dataset_name not in campaign.datasets:
            logger.warning(f"WARNING: {dataset_name} not in cmsdb campaign")
            continue
        
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if re.match(r"^(ww|wz|zz|www|wwz|wzz|zzz)$", dataset.name):
            dataset.add_tag("no_lhe_weights")
        #elif re.match(r"^dy_lep_m.*$", dataset.name):
        #    dataset.add_tag("is_dy")
        if re.match(r"^dy_.*$", dataset.name):
            dataset.add_tag("is_dy")
        if re.match(r"^dy_lep_m50.*$", dataset.name):
            dataset.add_tag("is_dy_m50")
        if re.match(r"^dy_2tau_m50.*$", dataset.name):
            dataset.add_tag("is_dy_tautau") 
        elif re.match(r"^wj.*$", dataset.name):
            dataset.add_tag("is_w")
        elif re.match(r"^h_ggf_tautau.*$", dataset.name):
            dataset.add_tag("is_signal")
            dataset.add_tag("is_ggf_signal")
        elif re.match(r"^h_vbf_tautau.*$", dataset.name):
            dataset.add_tag("is_signal")
            dataset.add_tag("is_vbf_signal")
        elif re.match(r"^(.*)h_tautau(.*)Filtered$", dataset.name):
            dataset.add_tag("is_signal")
            dataset.add_tag("is_vh_signal")
        elif dataset.name.startswith("tt_"):
            dataset.add_tag("is_tt")
        elif dataset.name.startswith("st_"):
            dataset.add_tag("is_st")
            #dataset.add_tag({"has_top", "ttbar", "tt"})
        #elif dataset.name.startswith("st_"):
        #    dataset.add_tag({"has_top", "single_top", "st"})
        
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

    cfg.x.default_hist_producer   = "cf_default"
    cfg.x.default_ml_model        = None
    cfg.x.default_inference_model = "main"
    cfg.x.default_categories      = "main"
    #cfg.x.default_variables = ("n_jet", "jet1_pt")
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

    
    # --------------------------------------------------------------------------------------------- #
    # Luminosity and Normalization
    # lumi values in inverse fb
    # TODO later: preliminary luminosity using norm tag. Must be corrected, when more data is available
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # --------------------------------------------------------------------------------------------- #

    if year == 2016:
        if postfix == "preVFP":
            cfg.x.luminosity = Number(19_500, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
        elif postfix == "postVFP":
            cfg.x.luminosity = Number(16_800, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
        else:
            raise RuntimeError(f"Wrong postfix: {postfix}")
    elif year == 2017:
        cfg.x.luminosity = Number(41_480, {
            "lumi_13TeV_2017": 0.02j,
            "lumi_13TeV_1718": 0.006j,
            "lumi_13TeV_correlated": 0.009j,
        })
    elif year == 2018:
        cfg.x.luminosity = Number(59_830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })
    else:
        raise RuntimeError(f"Wrong year: {year}")


    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    # TODO later: Error for run three not available yet. Using error from run 2.
    cfg.x.minbias_xs = Number(69.2, 0.046j)

    year_postfix = ""
    if year == 2016:
        year_postfix = "VFP" if postfix == "postVFP" else ""


    ################################################################################################
    # met settings
    ################################################################################################

    # --------------------------------------------------------------------------------------------- #
    # Adding triggers
    # --------------------------------------------------------------------------------------------- #

    if year == 2018:
        if channel == 'emu':
            from zttpol.config.triggers_emu import add_triggers_2018
        elif channel == 'etau':
            from zttpol.config.triggers_etau import add_triggers_2018
        elif channel == 'mutau':
            from zttpol.config.triggers_mutau import add_triggers_2018
        elif channel ==	'tautau':
            from zttpol.config.triggers_tautau import add_triggers_2018

        add_triggers_2018(cfg, postfix)

    
    else:
        raise RuntimeError(f"Wrong year: {year}. Check __init__.py in cmsdb campaign")


    # --------------------------------------------------------------------------------------------- #
    # Adding met filters
    # --------------------------------------------------------------------------------------------- #
    
    from zttpol.config.met_filters import add_met_filters
    add_met_filters(cfg)


    # --------------------------------------------------------------------------------------------- #
    # Adding external files e.g. goodlumi, normtag, SFs etc
    # restructure the postfix names to build appropriate tags
    # --------------------------------------------------------------------------------------------- #

    year2 = year%100 # 16, 17, 18
    #external_path_parent = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/corrections")
    #external_path_parent = os.path.join(corrdir, "corrections")
    #external_path_tail   = f"{year}{postfix}" if postfix else f"{year}"
    #external_path        = os.path.join(external_path_parent, f"{external_path_tail}")
    normtag_path         = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags"

    #print(f"external_path_parent : {external_path_parent}")
    #print(f"external_path        : {external_path}")

    #json_mirror = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/jsonpog-integration")
    #json_mirror = os.path.join(corrdir, "jsonpog-integration")
    #print(f"json_mirror          : {json_mirror}")

    #normtagjson = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json"
    normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
    goldenjson = None
    if year == 2016:
        goldenjson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
    elif year == 2017:
        goldenjson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
    elif year == 2018:
        goldenjson = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
    else:
        raise RuntimeError(f"Check year : {year}")

    #goldenjson  = glob(f"{external_path}/Lumi/*.json")[0]

    #print(f"GoldenJSON           : {goldenjson}")
    #print(f"NormtagJSON          : {normtagjson}")
    #print(f"sdcsdcs sdcascascsd cscsdc : {external_path}/Zpt/myZptCorrections.json.gz")
    
    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden"  : (goldenjson,  "v1"),  # noqa
            "normtag" : (normtagjson, "v1"),
        },
        # https://gitlab.cern.ch/cms-analysis-corrections/LUM
        "pu_sf"             : (f"{corrdir}/LUM/Run{run}-{year}{year_postfix}-UL-NanoAODv9/puWeights.json.gz",   "v1"), # PU
        # https://gitlab.cern.ch/cms-analysis-corrections/JME
        "jet_jerc"          : (f"{corrdir}/JME/Run{run}-{year}{year_postfix}-UL-NanoAODv15/jet_jerc.json.gz",    "v1"), # JEC
        "jet_veto_map"      : (f"{corrdir}/JME/Run{run}-{year}{year_postfix}-UL-NanoAODv15/jetvetomaps.json.gz", "v1"), # JetVeto
        # https://gitlab.cern.ch/cms-analysis-corrections/MUO
        "muon_sf"           : (f"{corrdir}/MUO/Run{run}-{year}{year_postfix}-UL-NanoAODv9/muon_Z.json.gz", "v1"), # Mu POG SF
        # -- https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossMuTauHlt.json?ref_type=heads
        # https://gitlab.cern.ch/cms-analysis-corrections/EGM
        "electron_sf"       : (f"{corrdir}/EGM/Run{run}-{year}{year_postfix}-UL-NanoAODv15/electron.json.gz",    "v1"), # Ele POG SF
        "electron_ss"       : (f"{corrdir}/EGM/Run{run}-{year}{year_postfix}-UL-NanoAODv15/electronSS_EtDependent.json.gz",  "v1"), # Ele Scale Smearing
        # -- https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/blob/main/data/TriggerScaleFactors/2022preEE/CrossEleTauHlt.json?ref_type=heads
        # https://gitlab.cern.ch/cms-analysis-corrections/TAU
        "tau_sf"            : (f"{corrdir}/TAU/Run{run}-{year}{year_postfix}-UL-NanoAODv15/tau.json.gz", "v1"), # TEC and ID SF from POG (VsJet Medium WP only)
        "gen_tau_sf"        : (f"{corrdir}/TAU/Run{run}-{year}{year_postfix}-UL-NanoAODv15/tau.json.gz", "v1"), # TEC and ID SF from POG (VsJet Medium WP only)
        "tes_sf"            : (f"{corrdir}/TAU/Run{run}-{year}{year_postfix}-UL-NanoAODv15/tau.json.gz", "v1"), # Tau ID SF (NEW FROM IC)
    })

    # --------------------------------------------------------------------------------------------- #
    # electron settings
    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    # --------------------------------------------------------------------------------------------- #
    
    electron_sf_tag = ""
    if year == 2018:
        electron_sf_tag = "2018"
    else:
        raise RuntimeError("wrong year")
    
    cfg.x.electron_sf_names = (
        "UL-Electron-ID-SF",
        electron_sf_tag,
        "wp80iso",
    )
    
    # --------------------------------------------------------------------------------------------- #
    # muon settings
    # names of muon correction sets and working points
    # (used in the muon producer)
    # --------------------------------------------------------------------------------------------- #

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
    

    # --------------------------------------------------------------------------------------------- #
    # jet settings
    # TODO: keep a single table somewhere that configures all settings: btag correlation, year
    #       dependence, usage in calibrator, etc
    # common jec/jer settings configuration
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107
    # --------------------------------------------------------------------------------------------- #

    # common jec/jer settings configuration
    if run == 2:
        # https://cms-jerc.web.cern.ch/Recommendations/#run-2
        # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=204
        # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=109
        jec_campaign = f"Summer19UL{year2}{postfix}"
        jec_version = {2016: "V7", 2017: "V5", 2018: "V5"}[year]
        jer_campaign = f"Summer{'20' if year == 2016 else '19'}UL{year2}{postfix}"
        jer_version = "JR" + {2016: "V3", 2017: "V2", 2018: "V2"}[year]
        jet_type = "AK4PFchs"
    
    elif run == 3:
        # https://cms-jerc.web.cern.ch/Recommendations/#2022
        jerc_postfix = {
            (22, ""): "_22Sep2023",
            (22, "EE"): "_22Sep2023",
            (23, ""): "Prompt23",
            (23, "BPix"): "Prompt23",
        }[(year2, year_postfix)]
        jec_campaign = f"Summer{year2}{year_postfix}{jerc_postfix}"
        jec_version = {
            (22, ""): "V3",
            (22, "EE"): "V3",
            (23, ""): "V2",
            (23, "BPix"): "V3",
        }[(year2, year_postfix)]
        jer_campaign = f"Summer{year2}{year_postfix}{jerc_postfix}"

        if year == 2023:
            jer_campaign += f"_Run{'Cv1234' if postfix == 'preBPix' else 'D'}"
        jer_version = "JR" + {2022: "V1", 2023: "V1", 2024: "V1"}[year]
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

    # updated jet id
    from columnflow.production.cms.jet import JetIdConfig
    cfg.x.jet_id = JetIdConfig(corrections={"AK4PUPPI_Tight": 2, "AK4PUPPI_TightLeptonVeto": 3})
    cfg.x.fatjet_id = JetIdConfig(corrections={"AK8PUPPI_Tight": 2, "AK8PUPPI_TightLeptonVeto": 3})

    # trigger sf corrector
    cfg.x.jet_trigger_corrector = "jetlegSFs"
    

    # --------------------------------------------------------------------------------------------- #
    # bjet settings
    # b-tag working points
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
    # TODO later: complete WP when data becomes available
    # --------------------------------------------------------------------------------------------- #

    btag_key = f"{year}{year_postfix}"
    cfg.x.btag_working_points = DotDict.wrap({
        "deepjet": {
            "loose"  : {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0583, "2023BPix": 0.0583, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.3086, "2023BPix": 0.3086, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.7183, "2022EE": 0.73,   "2023": 0.7183, "2023BPix": 0.7183, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.8111, "2023BPix": 0.8111, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.9512, "2023BPix": 0.9512, "2024": 0.0}[btag_key],
        },
        "robustParticleTransformer": {
            "loose"  : {"2022": 0.0849, "2022EE": 0.0897, "2023": 0.0849, "2023BPix": 0.0849, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.4319, "2022EE": 0.451,  "2023": 0.4319, "2023BPix": 0.4319, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.8482, "2022EE": 0.8604, "2023": 0.8482, "2023BPix": 0.8482, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.9151, "2022EE": 0.9234, "2023": 0.9151, "2023BPix": 0.9151, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9874, "2022EE": 0.9893, "2023": 0.9874, "2023BPix": 0.9874, "2024": 0.0}[btag_key],
        },
        "particleNet": {
            "loose"  : {"2022": 0.047,  "2022EE": 0.0499, "2023": 0.0499, "2023BPix": 0.0499, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.245,  "2022EE": 0.2605, "2023": 0.2605, "2023BPix": 0.2605, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.6915, "2023BPix": 0.6915, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.8033, "2023BPix": 0.8033, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.961,  "2022EE": 0.9664, "2023": 0.9664, "2023BPix": 0.9664, "2024": 0.0}[btag_key],
        },
    })
    # 2023 is dummy ... CORRECT IT LATER ===>> TODO


    # --------------------------------------------------------------------------------------------- #
    # tau settings
    # tau-id working points
    # --------------------------------------------------------------------------------------------- #

    cfg.x.deep_tau_tagger = "DeepTau2018v2p5"
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
                "tautau" : "Tight", #"VLoose",
            },
            "vs_j": {
                "emu"    : "VTight", # need to check
                "etau"   : "VTight",
                "mutau"  : "VTight", ##"Medium" : OLD,
                "tautau" : "VTight", ## VTight : Proposed by Imperial, was Medium in Run2
                # W A R N I N G !!! Medium is being used for ML, Change it to VTight for CP analysis
                #"tautau" : "Medium",
            },
        },
    })


    # tec config
    from columnflow.calibration.cms.tau import TECConfig
    corrector_kwargs = {
        ("etau", 3)  : {"wp": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_j['etau'], "wp_VSe": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_e['etau']},
        ("mutau", 3) : {"wp": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_j['mutau'], "wp_VSe": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_e['mutau']},
        ("tautau", 3): {"wp": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_j['tautau'], "wp_VSe": cfg.x.deep_tau_info[cfg.x.deep_tau_tagger].vs_e['tautau']},
    }[(channel, run)]
    #print(corrector_kwargs)
    cfg.x.tec = TECConfig(tagger=cfg.x.deep_tau_tagger,
                          correction_set=f'tau_es_dm_{cfg.x.deep_tau_tagger}_{year}_{postfix}',
                          corrector_kwargs=corrector_kwargs)


    ################################################################################################
    # electron settings
    ################################################################################################

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3

    # names of electron correction sets and working points
    from columnflow.production.cms.electron import ElectronSFConfig
    from columnflow.calibration.cms.egamma import EGammaCorrectionConfig


    # electron scale and smearing (eec and eer)
    cfg.x.ess = EGammaCorrectionConfig(
        scale_correction_set="Scale",
        scale_compound=True,
        smear_syst_correction_set="SmearAndSyst",
        systs=["scale_down", "scale_up", "smear_down", "smear_up"],
    )


    ################################################################################################
    # muon settings
    ################################################################################################

    # names of muon correction sets and working points
    # (used in the muon producer)
    from columnflow.production.cms.muon import MuonSFConfig

    # mec/mer
    from columnflow.calibration.cms.muon import MuonSRConfig
    cfg.x.muon_sr = MuonSRConfig(
        systs=["scale_up", "scale_down", "res_up", "res_down"],
    )
    
    
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
    # met settings
    # --------------------------------------------------------------------------------------------- #

    
    
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
            #"normalized_pu_weight": "normalized_pu_weight_{name}",
        },
    )


    for i, (jec_source, flag) in enumerate(all_jec_sources.items()):
        if not flag:
            continue
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5000 + 2 * i,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5001 + 2 * i,
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

    cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
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
    cfg.add_shift(name=f"tau_up", id=50, type="shape")
    cfg.add_shift(name=f"tau_down", id=51, type="shape")
    add_shift_aliases(
        cfg,
        "tau",
        {
            "tau_weight": "tau_weight_{direction}",
        },
    )
    cfg.add_shift(name=f"tau_trig_up", id=52, type="shape")
    cfg.add_shift(name=f"tau_trig_down", id=53, type="shape")
    add_shift_aliases(
        cfg,
        "tau_trig",
        {
            "tau_trigger_weight": "tau_trigger_weight_{direction}",
        },
    )
    #cfg.add_shift(name=f"tes_up", id=54, type="shape")
    #cfg.add_shift(name=f"tes_down", id=55, type="shape")
    #add_shift_aliases(
    #    cfg,
    #    "tes",
    #    {
    #        "tes_weight": "tes_weight_{direction}",
    #    },
    #)

    # --- >>> ele weight <<< --- #
    cfg.add_shift(name="e_up", id=90, type="shape")
    cfg.add_shift(name="e_down", id=91, type="shape")
    add_shift_aliases(
        cfg,
        "e",
        {
            "electron_idiso_weight": "electron_idiso_weight_{direction}",
        },
    )
    cfg.add_shift(name="e_trig_up", id=92, type="shape")
    cfg.add_shift(name="e_trig_down", id=93, type="shape")
    add_shift_aliases(
        cfg,
        "e_trig",
        {
            "electron_Ele30_WPTight_trigger_weight": "electron_Ele30_WPTight_trigger_weight_{direction}",
        },
    )
    cfg.add_shift(name="e_xtrig_up", id=94, type="shape")
    cfg.add_shift(name="e_xtrig_down", id=95, type="shape")
    add_shift_aliases(
        cfg,
        "e_xtrig",
        {
            "electron_xtrig_weight": "electron_xtrig_weight_{direction}",
        },
    )

    # --- >>> mu weight <<< --- #
    cfg.add_shift(name="mu_id_up", id=100, type="shape")
    cfg.add_shift(name="mu_id_down", id=101, type="shape")
    add_shift_aliases(
        cfg,
        "mu_id",
        {
            "muon_id_weight": "muon_id_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_iso_up", id=102, type="shape")
    cfg.add_shift(name="mu_iso_down", id=103, type="shape")
    add_shift_aliases(
        cfg,
        "mu_iso",
        {
            "muon_iso_weight": "muon_iso_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_trig_up", id=104, type="shape")
    cfg.add_shift(name="mu_trig_down", id=105, type="shape")
    add_shift_aliases(
        cfg,
        "mu_trig",
        {
            "muon_IsoMu24_trigger_weight": "muon_IsoMu24_trigger_weight_{direction}",
        },
    )
    cfg.add_shift(name="mu_xtrig_up", id=106, type="shape")
    cfg.add_shift(name="mu_xtrig_down", id=107, type="shape")
    add_shift_aliases(
        cfg,
        "mu_xtrig",
        {
            "muon_xtrig_weight": "muon_xtrig_weight_{direction}",
        },
    )

    # --- >>> pdf weight <<< --- #    
    cfg.add_shift(name="pdf_up", id=130, type="shape")
    cfg.add_shift(name="pdf_down", id=131, type="shape")
    add_shift_aliases(
        cfg,
        "pdf",
        {
            "pdf_weight": "pdf_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )
    # --- >>> tau spinner weight <<< --- #
    cfg.add_shift(name="tauspinner_up",   id=150, type="shape") # cp-even
    cfg.add_shift(name="tauspinner_down", id=151, type="shape") # cp-odd
    add_shift_aliases(
        cfg,
        "tauspinner",
        {
            "tauspinner_weight": "tauspinner_weight_{direction}",
        },
    )
    # --- >>> zpt weight <<< --- #    
    cfg.add_shift(name="zpt_up", id=160, type="shape")
    cfg.add_shift(name="zpt_down", id=161, type="shape")
    add_shift_aliases(
        cfg,
        "zpt",
        {
            "zpt_reweight": "zpt_reweight_{direction}",
        },
    )


    # --- >>> fake factor weight <<< --- # D U M M Y !!!!!!!!   
    cfg.add_shift(name="ff_up", id=170, type="shape")
    cfg.add_shift(name="ff_down", id=171, type="shape")
    add_shift_aliases(
        cfg,
        "ff",
        {
            "ff_weight": "ff_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )
    cfg.add_shift(name="ff_cls_corr_up", id=172, type="shape")
    cfg.add_shift(name="ff_cls_corr_down", id=173, type="shape")
    add_shift_aliases(
        cfg,
        "ff_cls_corr",
        {
            "ff_cls_corr_weight": "ff_cls_corr_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )

    cfg.add_shift(name="ff_ext_corr_up", id=180, type="shape")
    cfg.add_shift(name="ff_ext_corr_down", id=181, type="shape")
    add_shift_aliases(
        cfg,
        "ff_ext_corr",
        {
            "ff_ext_corr_weight": "ff_ext_corr_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )

    cfg.add_shift(name="top_pt_up", id=190, type="shape")
    cfg.add_shift(name="top_pt_down", id=191, type="shape")
    add_shift_aliases(
        cfg,
        "top_pt",
        {
            "top_pt_weight": "top_pt_weight_{direction}",
        },
    )

    # muon scale and resolution
    cfg.add_shift(name="mec_up", id=201, type="shape", tags={"mec"})
    cfg.add_shift(name="mec_down", id=202, type="shape", tags={"mec"})
    add_shift_aliases(cfg, "mec", {"Muon.pt": "Muon.pt_scale_{direction}"})

    cfg.add_shift(name="mer_up", id=203, type="shape", tags={"mer"})
    cfg.add_shift(name="mer_down", id=204, type="shape", tags={"mer"})
    add_shift_aliases(cfg, "mer", {"Muon.pt": "Muon.pt_res_{direction}"})


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
        "pu_weight"                             : get_shifts("minbias_xs"),
        "electron_idiso_weight"                 : get_shifts("e"),
        "electron_Ele30_WPTight_trigger_weight" : get_shifts("e_trig"),
        "electron_xtrig_weight"                 : get_shifts("e_xtrig"),
        "tau_weight"                            : get_shifts("tau"),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        "zpt_reweight"                          : get_shifts("zpt"),
        "top_pt_weight"                         : [],
    })
    cfg.x.event_mutau_weights = DotDict({
        "normalization_weight"                  : [],
        "pu_weight"                             : get_shifts("minbias_xs"),
        "muon_id_weight"                        : get_shifts("mu_id"),
        "muon_iso_weight"                       : get_shifts("mu_iso"),
        "muon_IsoMu24_trigger_weight"           : get_shifts("mu_trig"),
        "muon_xtrig_weight"                     : get_shifts("mu_xtrig"),
        "tau_weight"                            : get_shifts("tau"),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        "zpt_reweight"                          : get_shifts("zpt"),
        "top_pt_weight"                         : [],
    })
    """
    cfg.x.event_tautau_weights = DotDict({
        "normalization_weight"                  : [],
        "pu_weight"                             : get_shifts("minbias_xs"),
        "tau_weight"                            : get_shifts("tau"),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        "zpt_reweight"                          : get_shifts("zpt"),
        "top_pt_weight"                         : [],
    })
    """
    cfg.x.event_tautau_weights = DotDict({
        "normalization_weight"                  : [],
        "pu_weight"                             : get_shifts("minbias_xs"),
        #"tau_weight"                            : get_shifts("tau"),
        #"tau_trigger_weight"                    : get_shifts("tau_trig"),
        #"ff_weight"                             : [],
        ##"ff_cls_corr_weight"                    : [],
        ###"ff_ext_corr_weight"                    : [],
        ##"tauspinner_weight"                     : get_shifts("tauspinner"),
        ##"pdf_weight"                            : [],
        #"zpt_reweight"                          : get_shifts("zpt"),
        #"top_pt_weight"                         : [],
    })

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

    
    #---------------------------------------------------------------------------------------------#
    # Channels [possible to avoid as category id does the same]
    #---------------------------------------------------------------------------------------------#


    #---------------------------------------------------------------------------------------------#
    # LFN settings
    # Campaigns must have creator names as desy or IPHC
    # the main path is in the campaign __init__, the basepath is mentioned in the dataset info
    #---------------------------------------------------------------------------------------------#

    campaign_tag = cfg.campaign.x("custom").get("creator")
    if campaign_tag == "desy" or campaign_tag == "IPHC":
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            logger.info(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                logger.warning("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eoscms_redirector",
                #fs="wlcg_fs_imperial_redirector",
            )
            logger.info(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eoscms_redirector"
        #cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_imperial_redirector"
        
    #---------------------------------------------------------------------------------------------#
    # Add categories described in categorization.py
    #---------------------------------------------------------------------------------------------#

    from importlib import import_module

    module = import_module(f"zttpol.config.categories_{channel}")
    module.add_categories(cfg)

    """
    if channel == 'emu':
        from zttpol.config.categories_emu import add_categories
    elif channel == 'etau':
        from zttpol.config.categories_etau import add_categories
    elif channel == 'mutau':
        from zttpol.config.categories_mutau import add_categories
    elif channel == 'tautau':
        from zttpol.config.categories_tautau import add_categories

    add_categories(cfg)
    """

    
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
            # TauProds
            "TauProd.*",
            "GenTop_pt",
            # general event info
            "run", "luminosityBlock", "event", "LHEPdfWeight",
            "PV.x","PV.y","PV.z","PV.npvs","PV.npvsGood",
            "mc_weight",
            "PVBS.*",
            "SV.x","SV.y","SV.z",#"SV.pAngle",
            "Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP",
            #"trigger_ids",
            "single_triggered", "cross_triggered",
            #"single_e_triggered", "cross_e_triggered",
            #"single_mu_triggered", "cross_mu_triggered",
            "cross_tau_triggered", "cross_tau_jet_triggered",
            # --- new for categorization --- #
            "is_os", "is_iso_1", "is_iso_2", "is_low_mt", "is_b_veto",
            "is_real_1", "is_real_2", "is_fake_1", "is_fake_2",
            "is_lep_1", "is_pi_1", "is_pi_2", "is_rho_1", "is_rho_2",
            "is_a1_1pr_2pi0_1", "is_a1_1pr_2pi0_2",
            "is_a1_3pr_0pi0_1", "is_a1_3pr_0pi0_2",
            "is_a1_3pr_1pi0_1", "is_a1_3pr_1pi0_2",
            "is_ipsig_0to1_1",
            "has_0jet", "has_1jet", "has_2jet",
            #"LHE.NpLO", "LHE.NpNLO", "LHE.Vpt", "LHE.Njets",
        } | {
            f"GenPart.{var}" for var in [
                "pt", "eta", "phi", "mass",
                "status", "pdgId", "statusFlags",
            ]
        } | {
            f"GenZ.{var}" for var in [
                "pt", "eta", "phi", "mass",
            ]
        } | {
            f"GenZvis.{var}" for var in [
                "pt", "eta", "phi", "mass",
            ]            
        } | {
            f"PuppiMET.{var}" for var in [
                "pt_no_tes", "phi_no_tes",
                "pt_no_corr", "phi_no_corr",
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]
            #} | {
            #f"MET.{var}" for var in [
            #    "pt", "phi", "significance",
            #    "covXX", "covXY", "covYY",
            #]
        } | {
            f"Jet.{var}" for var in [
                "pt_no_corr", "eta_no_corr",
                "phi_no_corr", "mass_no_corr",
                "pt", "eta", "phi", "mass",
                "btagDeepFlavB", "hadronFlavour"
            ]
        } | {
            f"bJet.{var}" for var in [
                "pt", "eta", "phi", "mass",
                "btagDeepFlavB", "hadronFlavour"
            ]
        } | {
            f"trigJet.{var}" for var in [
                "pt", "eta", "phi", "mass",
            ]
        } | {
            f"metRecoilJet.{var}" for var in [
                "pt", "eta", "phi", "mass",
            ]
        } | { # raw tau in events before any selection
            f"RawTau.{var}" for var in [
                "pt","eta","phi","mass",
                "decayMode",
                "rawIdx",
            ]
        } | { # taus after good selection
            f"GoodTau.{var}" for var in [
                "pt","eta","phi","mass",
                "decayMode",
                "rawIdx",
            ]
        } | { # taus in hcand 
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", "IPx", "IPy", "IPz",
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu",
                "decayMode",
                "decayModeHPS",
                "decayModePNet",
                "genPartFlav",
                "rawIdx",
                "pt_no_tes", "mass_no_tes","SVx","SVy","SVz",
            ]
        } | {
            f"TauSpinner.weight_cp_{var}" for var in [
                "0", "0_alt", "0p25", "0p25_alt", "0p375",
                "0p375_alt", "0p5", "0p5_alt", "minus0p25", "minus0p25_alt"
            ]
        } | { # muon before any selection 
            f"RawMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # good muons after selection
            f"GoodMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # to check if extra lep veto work correctly
            f"VetoMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # to check if dilep veto work correctly
            f"DoubleVetoMuon.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx"
            ]
        } | { # muons from hcand
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", "IPx", "IPy", "IPz",
                "decayMode", "pfRelIso04_all","mT", "rawIdx","SVx","SVy","SVz",
            ]
        } | { # eles before any selection
            f"RawElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # eles after good selection
            f"GoodElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # to check if extra lep veto work correctly
            f"VetoElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # to check if dilep veto work correctly
            f"DoubleVetoElectron.{var}" for var in [
                "pt","eta","phi","mass","dxy",
                "decayMode",
                "rawIdx",
            ]
        } | { # electrons from hcand
            f"Electron.{var}" for var in [
                "pt_no_ss", "pt","eta","phi","mass",
                "dxy","dz", "charge", "IPx", "IPy", "IPz",
                "decayMode", "pfRelIso03_all", "mT", "rawIdx",
                "deltaEtaSC","SVx","SVy","SVz",
            ]
        } | {
            f"TrigObj.{var}" for var in [
                "id", "pt", "eta", "phi", "filterBits",
            ]
        } | {
            f"zcand.{var}" for var in [
                "pt","eta","phi","mass", "charge",
                "decayMode", "rawIdx", "isolation",
                "IPx", "IPy", "IPz", "IPsig",
                "SVx","SVy","SVz",
            ]
        } | {
            "GenTau.*", "GenTauProd.*",
        } | {
            f"zcandprod.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
                "pdgId", "tauIdx",
            ]
        } | {ColumnCollection.ALL_FROM_SELECTOR},
        "cf.MergeSelectionMasks": {
            "normalization_weight",
            "cutflow.*", "process_id", "category_ids",
        },
        "cf.UniteColumns": {
            "run","luminosityBlock","event","LHEPdfWeight",
            "PV.x","PV.y","PV.z","PV.npvs","PV.npvsGood",
            "PVBS.x", "PVBS.y", "PVBS.z",
            "mc_weight",
            #"SV.x","SV.y","SV.z",#"SV.pAngle",
            "Pileup.nTrueInt","Pileup.nPU","genWeight",
            "single_triggered", "cross_triggered",
            "cross_tau_triggered", "cross_tau_jet_triggered",
            # --- new for categorization --- #
            "is_os", "is_iso_1", "is_iso_2", "is_low_mt", "is_b_veto",
            "is_real_1", "is_real_2", "is_fake_1", "is_fake_2",
            "is_lep_1", "is_pi_1", "is_pi_2", "is_rho_1", "is_rho_2",
            "is_a1_1pr_2pi0_1", "is_a1_1pr_2pi0_2",
            "is_a1_3pr_0pi0_1", "is_a1_3pr_0pi0_2",
            "is_a1_3pr_1pi0_1", "is_a1_3pr_1pi0_2",
            "is_ipsig_0to1_1",
            "has_0jet", "has_1jet", "has_2jet",
            "PuppiMET.*", "Jet.*",
            #"TauSpinner.*",
            "tauspinner_weight",
            "tauspinner_weight_cpeven",    # 0
            "tauspinner_weight_cpeven_alt",# 0
            "tauspinner_weight_cpmix",     # 0.25
            "tauspinner_weight_cpmix_alt", # 0.25_alt
            "tauspinner_weight_cpmixm",    # -0.25
            "tauspinner_weight_cpmixm_alt",# -0.25_alt
            "tauspinner_weight_cpalpha0p375",     # 0.375
            "tauspinner_weight_cpalpha0p375_alt", # 0.375_alt
            "tauspinner_weight_cpodd",     # 0.5
            "tauspinner_weight_cpodd_alt", # 0.5
            "zcand.*", "zcandprod.*", "TauProd.*",
            "zcand_invm_fastMTT",
            "GenTau.*", "GenTauProd.*",
            "process_id", "category_ids",
            "PhiCP*",
            
        },
    })

    # --------------------------------------------------------------------------------------------- #
    # Adding hist hooks
    # --------------------------------------------------------------------------------------------- #

    cfg.x.regions_to_extrapolate_fake = "AB" # "AB" or "CD" or "C0D0"
    cfg.x.save_qcd = True
    from zttpol.config.hist_hooks import add_hist_hooks
    add_hist_hooks(cfg)

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
            "main"                    : True,
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
        "genmatch"       : True,
    })


    cfg.x.is_channel_specific = False
    cfg.x.channel_specific_info = DotDict.wrap({
        "etau"   : False,
        "mutau"  : False,
        "tautau" : True,
    })
