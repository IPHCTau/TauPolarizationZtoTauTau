# coding: utf-8

"""
Style definitions.
"""

from __future__ import annotations

import re
import random
from collections import defaultdict

import order as od

from columnflow.util import DotDict, try_int
from columnflow.types import Callable


def setup_plot_styles(config: od.Config) -> None:
    """
    Setup plot styles.
    """
    # general settings
    config.x.default_general_settings = {
        "cms_label": "wip", "whitespace_fraction": 0.31,
    }

    # default component configs
    gridspec = {
        "height_ratios": [3, 0.9],
    }
    legend = {
        "borderpad": 0, "borderaxespad": 1.2, "columnspacing": 1.8, "labelspacing": 0.28, "fontsize": 16,
        "cf_line_breaks": True, "cf_short_labels": False,
    }
    ratio = {
        "yloc": "center",
    }
    annotate = {
        "fontsize": 18, "style": "italic", "xycoords": "axes fraction", "xy": (0.035, 0.955),
    }

    # wide legend
    # - 3 columns, backgrounds in first 2 columns
    # - shortened process labels
    # - changed annotation (channel) position to fit right under legend
    wide_legend = {
        **legend,
        "ncols": 3, "loc": "upper left", "cf_entries_per_column": legend_entries_per_column, "cf_short_labels": True,
    }
    annotate_wide = {
        **annotate,
        "xy": (0.035, 0.765),
    }

    # wide extended legend, same as wide legend except
    # - process labels are not shortened
    # - annotation (channel) moved slightly down to fut under (now taller) legend
    wide_ext_legend = {
        **wide_legend,
        "cf_short_labels": False,
    }
    annotate_wide_ext = {
        **annotate_wide,
        "xy": (0.035, 0.750),
    }

    # extended y-range on ratio pad
    extended_ratio = {
        **ratio,
        "ylim": (0.42, 1.58),
    }

    # construct named style configs
    config.x.custom_style_config_groups = {
        "default": (default_cfg := {
            "gridspec_cfg": gridspec,
            "rax_cfg": ratio,
            "legend_cfg": legend,
            "annotate_cfg": annotate,
        }),
        "wide_legend": (wide_legend_cfg := {
            **default_cfg,
            "legend_cfg": wide_legend,
            "annotate_cfg": annotate_wide,
        }),
        "wide_ext_legend": {
            **wide_legend_cfg,
            "legend_cfg": wide_ext_legend,
            "annotate_cfg": annotate_wide_ext,
        },
        "ext_ratio": {
            "rax_cfg": extended_ratio,
        },
    }

    config.x.default_custom_style_config = "wide_legend"  # also via --custom-style-config
    config.x.default_blinding_threshold = 0  # also via --blinding-threshold



def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
    cfg = config

    # color scheme, see https://cms-analysis.docs.cern.ch/guidelines/plotting/colors
    cfg.x.colors = DotDict(
        # cms recommended 10-color scheme
        blue       = "#3f90da",
        yellow     = "#ffa90e",
        red        = "#bd1f01",
        grey       = "#94a4a2",
        purple     = "#832db6",
        brown      = "#a96b59",
        orange     = "#e76300",
        olive      = "#b9ac70",
        dark_grey  = "#717581",
        teal       = "#92dadd",
        # additional 10 colors
        navy       = "#004488",
        cyan       = "#00a6d6",
        magenta    = "#cc79a7",
        dark_green = "#117733",
        pink       = "#e78ac3",
        indigo     = "#5e4fa2",
        green      = "#44aa55",
        wine       = "#8c2d4f",
        sand       = "#d6b56c",
        slate      = "#536878",
        # some more colors (for unstacked plots)
        deep_blue  = "#01579b",
        vermillion = "#c43e00",
        forest_green = "#006644",
        royal_purple = "#6a1b94"
    )

    cfg.x.color_names = [
        "blue", "yellow", "red", "grey", "purple", "brown", "orange", "olive", "dark_grey", "teal",
        "navy", "cyan", "magenta", "dark_green", "pink", "indigo", "green", "wine", "sand", "slate"
    ]
    cfg.x.get_color_from_sequence = lambda i: cfg.x.colors[cfg.x.color_names[i % len(cfg.x.color_names)]]

    cfg.x.line_styles = ["solid", "dotted", "dashed", "dashdot"]
    cfg.x.get_line_style_from_sequence = lambda i: cfg.x.line_styles[i % len(cfg.x.line_styles)]
    
    cfg.x.proc_col_label = DotDict.wrap({
        'dy_m50toinf'                       : [cfg.x.colors.blue,   r"$Z \to e^+e^-/\mu^+\mu^-/\tau^+\tau^-$"],
        'dy_2e_or_2mu_m50toinf_nj'          : [cfg.x.colors.blue,   r"$Z \to e^+e^-/\mu^+\mu^-$"], # because tautau is vetoed from gen-level
        'dy_2tau_m50toinf_nj_LHEspin_minus' : [cfg.x.colors.yellow, r"$Z \to \tau_L^- \tau_R^+$"],
        'dy_2tau_m50toinf_nj_LHEspin_plus'  : [cfg.x.colors.red,    r"$Z \to \tau_R^- \tau_L^+$"],
        'dy_2tau'                           : [cfg.x.colors.yellow, r"$Z \to \tau^+\tau^-$"],
        'tt'                                : [cfg.x.colors.grey,   r"$t\bar{t}$"],
        'st'                                : [cfg.x.colors.purple, r"$st(\bar{t})(W)$"],
        'top'                               : [cfg.x.colors.grey,   r"$t\bar{t} + st(\bar{t})(W)$"],
        'w_lnu'                             : [cfg.x.colors.purple, r"$W \to \ell\nu$"],
        'vv'                                : [cfg.x.colors.brown,  r"$VV (W,Z)$"],
        'vvv'                               : [cfg.x.colors.orange, r"$W \to \ell\nu$"],
        'ewk'                               : [cfg.x.colors.olive,  r"$EWK-W/Z$"],
        'multiboson'                        : [cfg.x.colors.olive,  r"$VV(V)$"],
        'higgs'                             : [cfg.x.colors.teal,   r"$(W/Z)h \to \tau^+\tau^-$"],
        'qcd'                               : [cfg.x.colors.pink,   r"$QCD$"],
    })
    

    for proc_name, (color, label) in cfg.x.proc_col_label.items():
        if (p := config.get_process(proc_name, default=None)):
            p.color1 = color
            p.label = label
    



def legend_entries_per_column(ax, handles: list, labels: list, n_cols: int) -> list[int]:
    """
    Controls number of entries such that backgrounds are in the first n - 1 columns, and everything
    else in the last one.
    """
    # get number of background and remaining entries
    n_backgrounds = sum(1 for handle in handles if handle.__class__.__name__ == "StepPatch")
    n_other = len(handles) - n_backgrounds

    # fill number of entries per column
    entries_per_col = n_cols * [0]
    n_bkg_cols = n_cols
    # set last column if non-backgrounds are present
    if n_other:
        entries_per_col[-1] = n_other
        n_bkg_cols -= 1
    # fill background columns
    for i in range(n_bkg_cols):
        entries_per_col[i] = n_backgrounds // n_bkg_cols + (n_backgrounds % n_bkg_cols > i)

    return entries_per_col


def update_handles_labels_factory(remove_mc_stat_label: bool = False) -> Callable:
    """
    Factory to generate a function that updates legend handles and labels given some conditions passed as arguments.
    """
    def remove_mc_stat_label(handles: list, labels: list) -> None:
        for i, label in enumerate(labels):
            if re.match(r"^MC stat\.? unc.*$", label):
                labels.pop(i)
                handles.pop(i)
                break
            
