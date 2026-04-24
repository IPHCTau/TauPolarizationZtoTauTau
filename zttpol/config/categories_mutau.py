# coding: utf-8

# To handle many categories:
# https://github.com/columnflow/columnflow/commit/3a104d633fa47a8efc789f7aba054ed967017347#diff-01da7ecfbc4b8bb83460821147201da604f9825e32cbd51e365bdfcbc7cb0912R61
# https://github.com/columnflow/columnflow/issues/547
# https://github.com/columnflow/columnflow/issues/559

import law
import order as od

from zttpol.util import call_once_on_config
from columnflow.config_util import (
    add_category, create_category_combinations, CategoryGroup
)

logger = law.logger.get_logger(__name__)




@call_once_on_config()
def add_RealOrFake_categories(config: od.Config) -> None:
    add_category(config, name="real_2", id=1000000, selection="cat_real_2", label="prompt",     tags={"tau2isRealMC"})
    add_category(config, name="fake_2", id=2000000, selection="cat_fake_2", label="non-prompt", tags={"tau2isFakeMC"})

    
@call_once_on_config()
def add_njet_categories(config: od.Config) -> None:
    add_category(config, name="has_0j", id=100000, selection="cat_0j", label=r"$0$ jet",           tags={"has0j"})
    add_category(config, name="has_1j", id=200000, selection="cat_1j", label=r"$1$ jet",           tags={"has1j"})
    add_category(config, name="has_2j", id=300000, selection="cat_2j", label=r"$\geq{2}$ jets",    tags={"has2j"})
    

@call_once_on_config()
def add_ABCD_categories(config: od.Config) -> None:
    """
     just before the final leaf
     keep ids from 1500 with 1500 interval, up to 30000
    """
    # DESY
    add_category(config,name="DRnum",  id=10000,  selection="cat_os_noniso1_iso2_lowmt",      label="dr_num",   tags={"os","noniso1", "iso2", "lowmt" })
    add_category(config,name="DRden",  id=20000,  selection="cat_ss_noniso1_iso2_lowmt",      label="dr_den",   tags={"ss","noniso1", "iso2", "lowmt" })
    add_category(config,name="AR",     id=30000,  selection="cat_ss_iso1_iso2_lowmt",         label="ar",       tags={"ss","iso1",    "iso2", "lowmt" })
    add_category(config,name="SR",     id=40000,  selection="cat_os_iso1_iso2_lowmt",         label="sr",       tags={"os","iso1",    "iso2", "lowmt" })
    

@call_once_on_config()
def add_classifier_categories(config: od.Config) -> None:
    """
     just before the final leaf
     keep ids from 100 with 50 interval, up to 1400
    """
    # hardonic
    add_category(config,name="nodeDY",    id=100,  selection="cat_node_dy",            label="DY_node",    tags={"dy_node"})
    add_category(config,name="nodeFake",  id=200,  selection="cat_node_fake",          label="Fake_node",  tags={"fake_node"})
    add_category(config,name="nodeHiggs", id=300,  selection="cat_node_higgs",         label="Higgs_node", tags={"higgs_node"})


    
@call_once_on_config()
def add_DM_categories(config: od.Config) -> None:
    """
     final chain of categories
     reserve the ids from 1 to 50
    """
    # leptonic
    add_category(config, name="pi_2",          id=1,  selection="cat_pi_2",            label=r"$\tau_{h}\to\pi$",                                    tags={"tau2pi"    }) # h2 -> pi
    add_category(config, name="rho_2",         id=2,  selection="cat_rho_2",           label=r"$\tau_{h}\to\rho$",                                   tags={"tau2rho"   }) # h2 -> rho
    add_category(config, name="a1dm2_2",       id=3,  selection="cat_a1dm2_2",         label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                  tags={"tau2a1DM2" }) # h2 -> a1
    add_category(config, name="a1dm10_2",      id=4,  selection="cat_a1dm10_2",        label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                  tags={"tau2a1DM10"}) # h2 -> a1
    add_category(config, name="a1dm11_2",      id=5,  selection="cat_a1dm11_2",        label=r"$\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                  tags={"tau2a1DM11"}) # h2 -> a1


    
# ################### #
# main categorization #
# ################### #

@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    
    add_category(config,
                 name="mutau",
                 id=10000000,
                 selection="cat_mutau",
                 label=r"$\mu\tau_{h}$",
                 tags={"mutau"})

    #add_njet_categories(config)
    add_RealOrFake_categories(config)
    
    add_ABCD_categories(config)
    add_DM_categories(config)

    #add_classifier_categories(config)

    # ############################################################### #
    # To create combinations of categories                            #
    # the entire root -> leaf categories will be in the category_ids  #
    # other combinatorics will be produced in createHistograms task,  #
    # if mentioned                                                    #
    # ############################################################### #

    def name_fn(categories: dict[str, od.Category]) -> str:
        return "__".join(cat.name for cat in categories.values() if cat)


    def kwargs_fn(categories: dict[str, od.Category]):
        return {
            "id": sum([c.id for c in categories.values()]),
            "label": "+".join([c.label for c in categories.values()]),
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
        }

    main_categories = {
        "channel": CategoryGroup(['mutau'], is_complete=True, has_overlap=False),
        "RorF"   : CategoryGroup(['real_2'], is_complete=False, has_overlap=False),
        "abcd"   : CategoryGroup(['DRnum','DRden','AR','SR'], is_complete=True, has_overlap=False),
        "cp"     : CategoryGroup(['pi_2','rho_2','a1dm2_2','a1dm10_2','a1dm11_2'], is_complete=True, has_overlap=False),
    }

    create_category_combinations(config=config,
                                 categories=main_categories,
                                 name_fn=name_fn,
                                 parent_mode="safe",
                                 kwargs_fn=kwargs_fn,
                                 skip_existing=False)
    

    all_cats = [cat.name for cat, _, _ in config.walk_categories()]
    logger.warning(f"{len(all_cats)} categories created for mutau channel")
    
    
