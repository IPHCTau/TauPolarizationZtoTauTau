import law
logger = law.logger.get_logger(__name__)

from zttpol.config.config_zttpol import add_config

from cmsdb.campaigns.run2_2018_taupol_nano_v15 import campaign_run2_2018_taupol_nano_v15
#from cmsdb.campaigns.run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2
#from cmsdb.campaigns.run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2
#from cmsdb.campaigns.run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2
#from cmsdb.campaigns.run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2 import campaign_run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2

def get_cfg_type_from_cmd(luigi_parser, law_parser):
    parsed_args, _  = law_parser.parse_known_args(luigi_parser.cmdline_args)
    root_task = parsed_args.root_task
    cfg_name_in_cmd = parsed_args.configs if 'PlotVariables' in root_task else parsed_args.config
    islimited = True if cfg_name_in_cmd.endswith('_limited') else False
    isfull = False if islimited else True
    return islimited, isfull



def build_analysis(analysis=None,
                   era=None,
                   postfix=None,
                   channel=None,
                   islimited=True, isfull=False):
    campaign = {
        "2018"          : campaign_run2_2018_taupol_nano_v15,
        #"2022preEE"    : {"campaign" : campaign_run3_2022_preEE_tau_pol_nano_v15, "islimited": True, "isfull": True},
        #"2022preEE"    : {"campaign" : campaign_run3_2022_preEE_nano_cp_tau_v14_nanoprod_2024_v2,    "islimited": True, "isfull": True},
        #"2022postEE"   : {"campaign" : campaign_run3_2022_postEE_nano_cp_tau_v14_nanoprod_2024_v2,   "islimited": True, "isfull": False},
        #"2023preBPix"  : {"campaign" : campaign_run3_2023_preBPix_nano_cp_tau_v14_nanoprod_2024_v2,  "islimited": True, "isfull": False},
        #"2023postBPix" : {"campaign" : campaign_run3_2023_postBPix_nano_cp_tau_v14_nanoprod_2024_v2, "islimited": True, "isfull": False},
        #"2024"          : {"campaign" : campaign_run3_2024_taupol_nano_v15, "islimited": True, "isfull": False},
    }[f"{era}{postfix}"]

    if isfull:
        logger.warning(f"Full Campaign for {era}{postfix} : <{campaign.name}> - <{campaign.id}>")
        add_config(
            analysis,
            campaign.copy(),
            config_name=campaign.name,
            config_id=int(campaign.id),
            channel=channel)
    if islimited:
        logger.warning(f"Limited Campaign for {era}{postfix} : <{campaign.name}> - <{campaign.id+10}> - % only 1 root file per dataset will be considered %")
        add_config(
            analysis,
            campaign.copy(),
            config_name=f"{campaign.name}_limited",
            config_id=int(campaign.id)+10,
            limit_dataset_files=1,
            channel=channel)
    
