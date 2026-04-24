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
    add_category(config, name="real_1", id=1000000, selection="cat_real_1", label="prompt",     tags={"tau1isRealMC"})
    add_category(config, name="fake_1", id=2000000, selection="cat_fake_1", label="non-prompt", tags={"tau1isFakeMC"})


    
@call_once_on_config()
def add_njet_categories(config: od.Config) -> None:
    """
    """
    add_category(config, name="has_0j", id=100000, selection="cat_0j", label=r"$0$ jet",           tags={"has0j"})
    add_category(config, name="has_1j", id=200000, selection="cat_1j", label=r"$1$ jet",           tags={"has1j"})
    add_category(config, name="has_2j", id=300000, selection="cat_2j", label=r"$\geq{2}$ jets",    tags={"has2j"})
    

@call_once_on_config()
def add_ABCD_categories(config: od.Config) -> None:
    """
     just before the final leaf
     keep ids from 1500 with 1500 interval, up to 30000
    """
    # hardonic
    add_category(config,name="hadA",  id=10000,  selection="cat_ss_iso1_iso2_bveto",       label="A",   tags={"ss","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadB",  id=15000,  selection="cat_ss_noniso1_iso2_bveto",    label="B",   tags={"ss","noniso1","iso2",   "bveto"})
    add_category(config,name="hadA0", id=20000,  selection="cat_ss_iso1_noniso2_bveto",    label="A0",  tags={"ss","iso1",   "noniso2","bveto"})
    add_category(config,name="hadB0", id=25000,  selection="cat_ss_noniso1_noniso2_bveto", label="B0",  tags={"ss","noniso1","noniso2","bveto"})
    add_category(config,name="hadD0", id=30000,  selection="cat_os_iso1_noniso2_bveto",    label="D0",  tags={"os","iso1",   "noniso2","bveto"})
    add_category(config,name="hadC0", id=35000,  selection="cat_os_noniso1_noniso2_bveto", label="C0",  tags={"os","noniso1","noniso2","bveto"})
    add_category(config,name="hadD",  id=40000, selection="cat_os_iso1_iso2_bveto",       label="Signal Region",  tags={"os","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadC",  id=45000, selection="cat_os_noniso1_iso2_bveto",    label="Application Region",  tags={"os","noniso1","iso2",   "bveto"})
    

@call_once_on_config()
def add_classifier_categories(config: od.Config) -> None:
    """
     just before the final leaf
     keep ids from 100 with 50 interval, up to 1400
    """
    # hardonic
    add_category(config,name="nodeDY",    id=100,  selection="cat_tautau_node_dy",            label="DY_node",    tags={"dy_node"})
    add_category(config,name="nodeFake",  id=200,  selection="cat_tautau_node_fake",          label="Fake_node",  tags={"fake_node"})
    add_category(config,name="nodeHiggs", id=300,  selection="cat_tautau_node_higgs",         label="Higgs_node", tags={"higgs_node"})


    
@call_once_on_config()
def add_DM_categories(config: od.Config) -> None:
    """
     final chain of categories
     reserve the ids from 1 to 50
    """
    # hadronic
    add_category(config, name="pi_1",          id=1,  selection="cat_pi_1",            label=r"$\tau_{h}^{1}\to\pi$",                                                  tags={"tau1pi"       })  # h1 -> pi
    add_category(config, name="rho_1",         id=2,  selection="cat_rho_1",           label=r"$\tau_{h}^{1}\to\rho$",                                                 tags={"tau1rho"      })  # h1 -> rho
    add_category(config, name="a1dm2_1",       id=3,  selection="cat_a1dm2_1",         label=r"$\tau_{h}^{1}\to a_{1}(1\pi-2\pi^{0})$",                                tags={"tau1a1DM2"    })  # h1 -> a1
    add_category(config, name="a1dm10_1",      id=4,  selection="cat_a1dm10_1",        label=r"$\tau_{h}^{1}\to a_{1}(3\pi-0\pi^{0})$",                                tags={"tau1a1DM10"   })  # h1 -> a1
    add_category(config, name="a1dm11_1",      id=5,  selection="cat_a1dm11_1",        label=r"$\tau_{h}^{1}\to a_{1}(3\pi-1\pi^{0})$",                                tags={"tau1a1DM11"   })  # h1 -> a1
    add_category(config, name="pi_pi",         id=6,  selection="cat_pi_pi",           label=r"$\tau_{h}\to\pi-\tau_{h}\to\pi$",                                       tags={"pi_pi"        })  # 
    add_category(config, name="pi_rho",        id=7,  selection="cat_pi_rho",          label=r"$\tau_{h}\to\pi-\tau_{h}\to\rho$",                                      tags={"pi_rho"       })  # 
    add_category(config, name="pi_a1dm2",      id=8,  selection="cat_pi_a1dm2",        label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                     tags={"pi_a1DM2"     })  # 
    add_category(config, name="pi_a1dm10",     id=9,  selection="cat_pi_a1dm10",       label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                     tags={"pi_a1DM10"    })  # 
    add_category(config, name="pi_a1dm11",     id=10, selection="cat_pi_a1dm11",       label=r"$\tau_{h}\to\pi-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                     tags={"pi_a1DM11"    })  # 
    add_category(config, name="rho_rho",       id=11, selection="cat_rho_rho",         label=r"$\tau_{h}\to\rho-\tau_{h}\to\rho$",                                     tags={"rho_rho"      })  # 
    add_category(config, name="rho_a1dm2",     id=12, selection="cat_rho_a1dm2",       label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",                    tags={"rho_a1DM2"    })  # 
    add_category(config, name="rho_a1dm10",    id=13, selection="cat_rho_a1dm10",      label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",                    tags={"rho_a1DM10"   })  # 
    add_category(config, name="rho_a1dm11",    id=14, selection="cat_rho_a1dm11",      label=r"$\tau_{h}\to\rho-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",                    tags={"rho_a1DM11"   })  # 
    add_category(config, name="a1dm2_a1dm2",   id=15, selection="cat_a1dm2_a1dm2",     label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(1\pi-2\pi^{0})$",   tags={"a1DM2_a1DM2"  })  # 
    add_category(config, name="a1dm2_a1dm10",  id=16, selection="cat_a1dm2_a1dm10",    label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",   tags={"a1DM2_a1DM10" })  # 
    add_category(config, name="a1dm2_a1dm11",  id=17, selection="cat_a1dm2_a1dm11",    label=r"$\tau_{h}\to a_{1}(1\pi-2\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM2_a1DM11" })  # 
    add_category(config, name="a1dm10_a1dm10", id=18, selection="cat_a1dm10_a1dm10",   label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})-\tau_{h}\to a_{1}(3\pi-0\pi^{0})$",   tags={"a1DM10_a1DM10"})  # 
    add_category(config, name="a1dm10_a1dm11", id=19, selection="cat_a1dm10_a1dm11",   label=r"$\tau_{h}\to a_{1}(3\pi-0\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM10_a1DM11"})  # 
    add_category(config, name="a1dm11_a1dm11", id=20, selection="cat_a1dm11_a1dm11",   label=r"$\tau_{h}\to a_{1}(3\pi-1\pi^{0})-\tau_{h}\to a_{1}(3\pi-1\pi^{0})$",   tags={"a1DM11_a1DM11"})  # 

    
# ################### #
# main categorization #
# ################### #

@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_category(config,
             name="tautau",
             id=10000000,
             selection="cat_tautau",
             label=r"$\tau_{h} \tau_{h}$",
             tags={"tautau"})

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
        "channel": CategoryGroup(['tautau'], is_complete=True, has_overlap=False),
        "RorF"   : CategoryGroup(['real_1'], is_complete=False, has_overlap=False),
        "abcd"   : CategoryGroup(['hadA','hadB','hadA0','hadB0','hadC0','hadD0','hadC','hadD'], is_complete=True, has_overlap=False),
        "cp"     : CategoryGroup(['pi_pi','pi_rho','pi_a1dm2','pi_a1dm10','pi_a1dm11',
                                  'rho_rho','rho_a1dm2','rho_a1dm10','rho_a1dm11',
                                  'a1dm2_a1dm10','a1dm2_a1dm11',
                                  'a1dm10_a1dm10','a1dm10_a1dm11',
                                  'a1dm11_a1dm11'], is_complete=True, has_overlap=False),
    }

    create_category_combinations(config=config,
                                 categories=main_categories,
                                 name_fn=name_fn,
                                 parent_mode="safe",
                                 kwargs_fn=kwargs_fn,
                                 skip_existing=False)
    

    all_cats = [cat.name for cat, _, _ in config.walk_categories()]
    logger.warning(f"{len(all_cats)} categories created for tautau channel")
    
