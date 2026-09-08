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
    add_category(config, name="real_2", id=4000000, selection="cat_real_2", label="prompt",     tags={"tau2isRealMC"})
    add_category(config, name="fake_2", id=7000000, selection="cat_fake_2", label="non-prompt", tags={"tau2isFakeMC"})

@call_once_on_config()
def add_ABCD_categories(config: od.Config) -> None:
    """
     just before the final leaf
     keep ids from 1500 with 1500 interval, up to 30000
    """
    add_category(config,name="A",   id=17000,   selection="cat_os_noniso1", label="A",  tags={"os","noniso"})
    add_category(config,name="B",   id=18000,   selection="cat_ss_noniso1", label="B",  tags={"ss","noniso"})
    add_category(config,name="C",   id=19000,   selection="cat_ss_iso1",    label="C",  tags={"ss","iso"})
    add_category(config,name="D",   id=20000,   selection="cat_os_iso1",    label="SR", tags={"os","iso"})

# ################### #
# main categorization #
# ################### #

@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_category(config,
                 name="emu",
                 id=10000000,
                 selection="cat_emu",
                 label=r"$e\mu$",
                 tags={"emu"})

    add_RealOrFake_categories(config)
    add_ABCD_categories(config)

    # ############################################################### #
    # To create combinations of categories                            #
    # the entire root -> leaf categories will be in the category_ids  #
    # other combinatorics will be produced in createHistograms task,  #
    # if mentioned                                                    #
    # ############################################################### #
    def name_fn(categories: dict[str, od.Category]) -> str:
        return "__".join(cat.name for cat in categories.values() if cat)


    def kwargs_fn(categories: dict[str, od.Category], add_qcd_group: bool = True):
        # build auxiliary information
        aux = {}
        # DRnum, DRden, AR and SR belonging to the same combination
        # receive an identical qcd_group.
        if add_qcd_group and "abcd" in categories:
            aux["qcd_group"] = name_fn({
                name: cat
                for name, cat in categories.items()
                if name != "abcd"
            })

        # return the desired kwargs
        return {
            "id": sum([c.id for c in categories.values()]),
            "label": "+".join([c.label for c in categories.values()]),
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
            "aux": aux,
        }

    main_categories = {
        "channel": CategoryGroup(['emu'], is_complete=True, has_overlap=False),
        "RorF1"  : CategoryGroup(['real_1'], is_complete=False, has_overlap=False),
        "RorF2"  : CategoryGroup(['real_2'], is_complete=False, has_overlap=False),
        "abcd"   : CategoryGroup(['A','B','C','D'], is_complete=True, has_overlap=False),
    }

    create_category_combinations(config=config,
                                 categories=main_categories,
                                 name_fn=name_fn,
                                 parent_mode="safe",
                                 kwargs_fn=kwargs_fn,
                                 skip_existing=False)


    all_cats = [cat.name for cat, _, _ in config.walk_categories()]
    logger.warning(f"{len(all_cats)} categories created for emu channel")
 
