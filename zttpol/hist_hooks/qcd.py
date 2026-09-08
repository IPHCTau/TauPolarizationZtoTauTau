# coding: utf-8

"""
Histogram hooks for QCD estimation.
"""

from __future__ import annotations

import collections
import functools

import law
import order as od
import scinum as sn

from columnflow.util import maybe_import, DotDict
from columnflow.types import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    hist = maybe_import("hist")


logger = law.logger.get_logger(__name__)


# Z -> tau tau ABCD category names. The original HBT implementation expects
# the four regions to be represented by adjacent ``os/ss`` and
# ``iso/noniso`` category components. This analysis instead uses one combined
# component per region.
QCD_REGION_NAMES = {
    "os_iso": "SR",
    "ss_iso": "AR",
    "os_noniso": "DRnum",
    "ss_noniso": "DRden",
}


def get_qcd_group(cat_inst: od.Category) -> str:
    """Return and validate the QCD group assigned during category construction."""
    qcd_group = getattr(cat_inst.x, "qcd_group", None)
    if not qcd_group:
        raise ValueError(
            f"category '{cat_inst.name}' has QCD-region tags but no 'qcd_group' auxiliary value",
        )
    return qcd_group


def replace_qcd_region(category_name: str, region_key: str) -> str:
    """Replace SR/AR/DRnum/DRden in a combined category name."""
    if region_key not in QCD_REGION_NAMES:
        raise ValueError(f"unknown QCD region key '{region_key}'")

    parts = category_name.split("__")
    region_names = set(QCD_REGION_NAMES.values())
    region_indices = [index for index, part in enumerate(parts) if part in region_names]

    if len(region_indices) != 1:
        raise ValueError(
            f"cannot replace the QCD region in category '{category_name}': "
            f"expected exactly one of {sorted(region_names)}",
        )

    parts[region_indices[0]] = QCD_REGION_NAMES[region_key]
    return "__".join(parts)


# helper to convert a histogram to a number object containing bin values and uncertainties
# from variances stored in an array of values
def hist_to_num(h: hist.Histogram, unc_name=str(sn.DEFAULT)) -> sn.Number:
    return sn.Number(h.values(), {unc_name: h.variances()**0.5})


# helper to integrate values stored in an array based number object
def integrate_num(num: sn.Number, axis=None) -> sn.Number:
    return sn.Number(
        nominal=num.nominal.sum(axis=axis),
        uncertainties={
            unc_name: (
                (unc_values_up**2).sum(axis=axis)**0.5,
                (unc_values_down**2).sum(axis=axis)**0.5,
            )
            for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
        },
    )


# helper to ensure that a specific category exists on the "category" axis of a histogram
def ensure_category(h: hist.Histogram, category_name: str) -> hist.Histogram:
    cat_axis = h.axes["category"]
    if category_name in cat_axis:
        return h
    dummy_fill = {ax.name: ax[0] for ax in h.axes if ax.name != "category"}
    h.fill(**dummy_fill, category=category_name, weight=0.0)
    return h


def add_hooks(analysis_inst: od.Analysis) -> None:
    """
    Add histogram hooks to a analysis.
    """
    def qcd_estimation_per_config(
        task: law.Task,
        config_inst: od.Config,
        hists: dict[od.Process, Any],
        requested_category: str | None = None,
        empty_bin_value: float = 0.0,
        # which source(s) of variance to keep
        variance_strategy: Literal["mc_data", "mc", "data", "none"] = "mc_data",
        # fill *empty_bin_value* (zero) into values (variances) in case one of the two integrals for the transfer factor
        # calculation is negative
        fill_empty_negative_norms: bool = True,
        # residual filling of *empty_bin_value* (zero) into values (variances) where the bin content is <= 0
        fill_empty_residual: bool = True,
        # whether to do sum leaf categories first for all data and mc histograms and then performing the ABCD method,
        # opposed to doing the ABCD method per leaf category and then summing up the resulting qcd histograms
        sum_leaves_first: bool = True,
        # strategy for shape transfer
        shape_transfer: Literal["from_os_noniso", "from_ss_iso"] = "from_os_noniso",
        **kwargs,
    ) -> dict[od.Process, Any]:
        import numpy as np
        import hist

        if variance_strategy not in (known_variance_strategies := {"mc_data", "mc", "data", "none"}):
            raise ValueError(
                f"unknown variance strategy '{variance_strategy}', known strategies are: "
                f"{', '.join(known_variance_strategies)}",
            )

        if shape_transfer not in (known_shape_transfers := {"from_os_noniso", "from_ss_iso"}):
            raise ValueError(
                f"unknown shape transfer strategy '{shape_transfer}', known strategies are: "
                f"{', '.join(known_shape_transfers)}",
            )

        # get the qcd process
        qcd_proc = config_inst.get_process("qcd", default=None)
        if not qcd_proc:
            return hists

        # extract all unique category names and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_names = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            category_names.update(list(cat_ax))

        def get_region_key(cat_inst: od.Category) -> str | None:
            """Return the ABCD region represented by a category's tags."""
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                return "os_iso"
            if cat_inst.has_tag({"os", "noniso"}, mode=all):
                return "os_noniso"
            if cat_inst.has_tag({"ss", "iso"}, mode=all):
                return "ss_iso"
            if cat_inst.has_tag({"ss", "noniso"}, mode=all):
                return "ss_noniso"
            return None

        # get ABCD categories corresponding to the requested category
        requested_group = None
        requested_region = None
        cat_inst_req = None
        # determine the group corresponding to the requested category
        # (disable this check if the qcd estimation should be done on granular leaf categories again)
        if sum_leaves_first and requested_category:
            def find_cat(cat_name, cat_tags, check_group=None):
                cat_inst = config_inst.get_category(cat_name)
                if not cat_inst.has_tag(cat_tags, mode=all):
                    raise ValueError(
                        f"requested category {cat_name} does not have the required tags {cat_tags} for the "
                        f"ABCD method.",
                    )
                if check_group and get_qcd_group(cat_inst) != check_group:
                    raise ValueError(
                        f"requested category {cat_name} is not part of the same ABCD group {check_group} as the "
                        f"other categories.",
                    )
                return cat_inst

            requested_cat_inst = config_inst.get_category(requested_category)

            if requested_cat_inst.has_tag({"os", "iso"}, mode=all):
                # The requested category is the SR category itself.
                requested_region = "os_iso"
                cat_inst_req = DotDict()
                cat_inst_req.os_iso = find_cat(requested_category, {"os", "iso"})
                requested_group = get_qcd_group(cat_inst_req.os_iso)
                cat_inst_req.ss_iso = find_cat(
                    replace_qcd_region(requested_category, "ss_iso"),
                    {"ss", "iso"},
                    requested_group,
                )
                cat_inst_req.os_noniso = find_cat(
                    replace_qcd_region(requested_category, "os_noniso"),
                    {"os", "noniso"},
                    requested_group,
                )
                cat_inst_req.ss_noniso = find_cat(
                    replace_qcd_region(requested_category, "ss_noniso"),
                    {"ss", "noniso"},
                    requested_group,
                )
            elif requested_cat_inst.has_tag({"ss", "iso"}, mode=all):
                requested_region = "ss_iso"
                requested_group = get_qcd_group(requested_cat_inst)
            elif requested_cat_inst.has_tag({"os", "noniso"}, mode=all):
                requested_region = "os_noniso"
                requested_group = get_qcd_group(requested_cat_inst)
            elif requested_cat_inst.has_tag({"ss", "noniso"}, mode=all):
                requested_region = "ss_noniso"
                requested_group = get_qcd_group(requested_cat_inst)
            else:
                # The requested category is a parent such as mutau__real_2.
                # Its name is identical to the qcd_group assigned to the four
                # generated children SR, AR, DRnum and DRden.
                requested_region = "all"
                requested_group = requested_category
                logger.info(
                    "requested parent category %s; use ABCD group %s",
                    requested_category,
                    requested_group,
                )

            # For parent and control-region requests, find the four generated
            # ABCD category roots through their common qcd_group.  Their
            # descendants are collected from the histogram axes below.
            if cat_inst_req is None:
                grouped_categories = collections.defaultdict(list)
                for cat_inst, _, _ in config_inst.walk_categories():
                    if getattr(cat_inst.x, "qcd_group", None) != requested_group:
                        continue
                    if region_key := get_region_key(cat_inst):
                        grouped_categories[region_key].append(cat_inst)

                missing_regions = set(QCD_REGION_NAMES) - set(grouped_categories)
                if missing_regions:
                    raise ValueError(
                        f"QCD group '{requested_group}' is missing category roots for regions "
                        f"{sorted(missing_regions)}",
                    )

                # Prefer the least granular category when multiple candidates
                # exist.  It owns all finer descendants needed for summation.
                cat_inst_req = DotDict({
                    region_key: min(cats, key=lambda cat: cat.name.count("__"))
                    for region_key, cats in grouped_categories.items()
                })

        # create qcd groups
        qcd_groups: dict[str, dict[str, list[od.Category]]] = collections.defaultdict(DotDict)
        for cat_name in category_names:
            cat_inst = config_inst.get_category(cat_name)
            # store references to the four category objects
            region_key = get_region_key(cat_inst)
            if region_key is None:
                continue

            # Standalone base categories such as SR, AR, DRnum and DRden carry
            # the region tags as well, but they are not generated combinations
            # and therefore have no qcd_group auxiliary value.  Only grouped
            # combination categories participate in the ABCD estimate.
            qcd_group = getattr(cat_inst.x, "qcd_group", None)
            if not qcd_group:
                continue

            qcd_groups[qcd_group].setdefault(region_key, []).append(cat_inst)

            # store the group corresponding to the requested category (if set)
            if requested_group and cat_inst_req is not None:
                if any(c.has_category(cat_name, deep=True) for c in cat_inst_req.values()):
                    grouped_cats = qcd_groups[requested_group].setdefault(region_key, [])
                    if cat_inst not in grouped_cats:
                        grouped_cats.append(cat_inst)

        # get complete qcd groups, potentially only selecting the one corresponding to the requested category
        if requested_group:
            complete_groups = [requested_group] if len(qcd_groups[requested_group]) == 4 else []
        else:
            complete_groups = [
                name for name, cats in qcd_groups.items()
                if len(cats) == 4
            ]

        # nothing to do if there are no complete groups
        if not complete_groups:
            discovered_regions = sorted(qcd_groups.get(requested_group, {}))
            task.logger.warning(
                "no complete ABCD groups found, skipping QCD estimation; "
                f"requested group={requested_group}, discovered regions={discovered_regions}",
            )
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        if not mc_hists:
            task.logger.warning("no MC histograms found, skipping QCD estimation")
            return hists
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not data_hists:
            task.logger.warning("no data histograms found, skipping QCD estimation")
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]

            if not requested_group:
                for key, cats in group.items():
                    if len(cats) > 1:
                        raise ValueError(f"ABCD group {group_name} has multiple categories for region {key}")

            # get the corresponding histograms and convert them to number objects, each one storing an array of values
            # with uncertainties
            # shapes: (SHIFT, VAR)
            def get_hist(h: hist.Histogram, region_name: str) -> hist.Histogram:
                # define intermediate categories to sum over if necessary
                for cat in group[region_name]:
                    h = ensure_category(h, cat.name)
                h = h[{"category": [hist.loc(cat.name) for cat in group[region_name]]}]
                return h[{"category": sum}]

            os_noniso_mc = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

            # data will always have a single shift whereas mc might have multiple,
            # broadcast numbers in-place manually if necessary
            if (n_shifts := mc_hist.axes["shift"].size) > 1:
                def broadcast_data_num(num: sn.Number) -> None:
                    num._nominal = np.repeat(num.nominal, n_shifts, axis=0)
                    for name, (unc_up, unc_down) in num._uncertainties.items():
                        num._uncertainties[name] = (
                            np.repeat(unc_up, n_shifts, axis=0),
                            np.repeat(unc_down, n_shifts, axis=0),
                        )
                broadcast_data_num(os_noniso_data)
                broadcast_data_num(ss_noniso_data)
                broadcast_data_num(ss_iso_data)

            # helper to warn about negative bins
            def warn_negative_bins(neg_mask: np.ndarray, region_name: str) -> None:
                if neg_mask.any():
                    shift_ids = [sid for neg, sid in zip(neg_mask, mc_hist.axes["shift"]) if neg]
                    shifts = list(map(config_inst.get_shift, shift_ids))
                    logger.warning(
                        f"negative QCD integral in {region_name} region for group {group_name} and shifts: "
                        f"{', '.join(shift.name for shift in shifts)}",
                    )

            # re-assign histograms to BCD regions depending on shape transfer strategy
            # (B: shape transfer, C/D: transfer factor)
            if shape_transfer == "from_os_noniso":
                b_data, b_mc = os_noniso_data, os_noniso_mc
                c_data, c_mc = ss_iso_data, ss_iso_mc
                d_data, d_mc = ss_noniso_data, ss_noniso_mc
                c_region_name, d_region_name = "ss_iso", "ss_noniso"
            else:  # from_ss_iso
                b_data, b_mc = ss_iso_data, ss_iso_mc
                c_data, c_mc = os_noniso_data, os_noniso_mc
                d_data, d_mc = ss_noniso_data, ss_noniso_mc
                c_region_name, d_region_name = "os_noniso", "ss_noniso"

            # Data-minus-MC residuals in the three control regions.  These are
            # the QCD contributions that are stored when AR, DRnum or DRden is
            # requested directly.
            residual_nums = {
                "os_noniso": os_noniso_data - os_noniso_mc,
                "ss_iso": ss_iso_data - ss_iso_mc,
                "ss_noniso": ss_noniso_data - ss_noniso_mc,
            }

            # Get the SR QCD shape from B.
            # Shapes: (SHIFT, VAR)
            qcd_sr: sn.Number = b_data - b_mc

            # get integrals to compute the transfer factor
            # shapes: (SHIFT,)
            c_int: sn.Number = integrate_num(c_data, axis=1) - integrate_num(c_mc, axis=1)
            d_int: sn.Number = integrate_num(d_data, axis=1) - integrate_num(d_mc, axis=1)

            # check negative integrals in shift bins
            c_neg_mask = c_int() <= 0
            warn_negative_bins(c_neg_mask, c_region_name)
            d_neg_mask = d_int() <= 0
            warn_negative_bins(d_neg_mask, d_region_name)

            # ABCD method
            # shape: (SHIFT, VAR)
            qcd_sr = qcd_sr * ((c_int / d_int)[:, None])

            def finite(arr, val=0.0):
                arr = arr.copy()
                arr[~np.isfinite(arr)] = val
                return arr

            def number_to_arrays(
                num: sn.Number,
                mc_unc_names: list[str],
                data_unc_names: list[str],
            ) -> tuple[np.ndarray, np.ndarray]:
                """Convert a residual number to value and variance arrays."""
                values = finite(num())
                variances = num(sn.UP, sn.ALL, unc=True)**2

                if variance_strategy == "mc":
                    variances = num(sn.UP, mc_unc_names, unc=True)**2
                elif variance_strategy == "data":
                    variances = num(sn.UP, data_unc_names, unc=True)**2
                elif variance_strategy == "none":
                    variances = np.zeros_like(values)

                return values, finite(variances)

            all_mc_unc_names = [
                "os_noniso_mc",
                "ss_iso_mc",
                "ss_noniso_mc",
            ]
            all_data_unc_names = [
                "os_noniso_data",
                "ss_iso_data",
                "ss_noniso_data",
            ]

            region_results = {}
            region_results["os_iso"] = number_to_arrays(
                qcd_sr,
                all_mc_unc_names,
                all_data_unc_names,
            )
            region_results["os_noniso"] = number_to_arrays(
                residual_nums["os_noniso"],
                ["os_noniso_mc"],
                ["os_noniso_data"],
            )
            region_results["ss_iso"] = number_to_arrays(
                residual_nums["ss_iso"],
                ["ss_iso_mc"],
                ["ss_iso_data"],
            )
            region_results["ss_noniso"] = number_to_arrays(
                residual_nums["ss_noniso"],
                ["ss_noniso_mc"],
                ["ss_noniso_data"],
            )

            # A non-positive C or D normalization invalidates only the SR
            # prediction; the observed control-region residuals remain valid.
            sr_values, sr_variances = region_results["os_iso"]
            neg_int_mask = c_neg_mask | d_neg_mask
            if fill_empty_negative_norms:
                sr_values[neg_int_mask, :] = empty_bin_value
                sr_variances[neg_int_mask, :] = 0.0

            # Apply the configured residual zero filling independently in every
            # region before a possible parent-category sum.
            if fill_empty_residual:
                for region_key, (values, variances) in region_results.items():
                    zero_mask = values <= 0
                    if region_key == "os_iso" and not fill_empty_negative_norms:
                        zero_mask &= ~neg_int_mask[:, None]
                    values[zero_mask] = empty_bin_value
                    variances[zero_mask] = 0.0

            # Select what should be stored in the requested histogram category.
            # Parent categories contain the sum over the four mutually
            # exclusive ABCD regions.
            if requested_region == "all":
                qcd_values = sum(values for values, _ in region_results.values())
                qcd_variances = sum(variances for _, variances in region_results.values())
            else:
                output_region = requested_region or "os_iso"
                qcd_values, qcd_variances = region_results[output_region]

            logger.info(
                "QCD ABCD group %s: C(data-MC)=%s, D(data-MC)=%s, "
                "stored region=%s, stored QCD yield=%s",
                group_name,
                c_int().tolist(),
                d_int().tolist(),
                requested_region or "os_iso",
                qcd_values.sum(axis=1).tolist(),
            )

            # ensure that the requested category exists in the qcd histogram (if set)
            if requested_group:
                qcd_hist = ensure_category(qcd_hist, requested_category)

            # insert values into the qcd histogram
            # For a requested parent category, group.os_iso can contain several
            # descendant categories and its ordering is not deterministic.  In
            # that case, always fill the category that the plotting task
            # actually requested.
            output_category = (
                requested_category
                if requested_group
                else group.os_iso[0].name
            )
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == output_category:
                    qcd_hist.view().value[cat_index, ...] = qcd_values
                    qcd_hist.view().variance[cat_index, ...] = qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} for category "
                    f"{output_category}",
                )

        return hists

    def qcd_estimation(
        task: law.Task,
        hists: dict[od.Config, dict[od.Process, Any]],
        category_name: str,
        variable_name: str,
        **kwargs,
    ) -> dict[od.Config, dict[od.Process, Any]]:
        return {
            config_inst: qcd_estimation_per_config(
                task,
                config_inst,
                hists[config_inst],
                requested_category=category_name,
                **kwargs,
            )
            for config_inst in hists.keys()
        }

    # add different hook variations
    analysis_inst.x.hist_hooks.qcd = qcd_estimation
    analysis_inst.x.hist_hooks.qcd_zerofill = functools.partial(qcd_estimation, empty_bin_value=1e-5)
    analysis_inst.x.hist_hooks.qcd_from_ss_iso = functools.partial(qcd_estimation, shape_transfer="from_ss_iso")
    analysis_inst.x.hist_hooks.qcd_raw = functools.partial(
        qcd_estimation,
        fill_empty_negative_norms=False,
        fill_empty_residual=False,
    )
