from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..utils.corrections.function_alignment import build_norm_with_fallback
from ..utils.curve_fitting import (
    exp_decay,
    exp_growth,
    fit_and_select_best,
    linear_func,
    parabolic_func,
    prediction_std,
)


def _append_to_log(log_file_path: Optional[str], log_message):
    if not log_file_path:
        return
    with open(log_file_path, "a", encoding="utf-8", errors="replace") as log_file:
        print(f" {log_message}", file=log_file)


def _truthy_mask(series: pd.Series) -> pd.Series:
    normalized = series.astype(str).str.strip().str.lower()
    return normalized.isin(["true", "1", "1.0", "yes", "y", "t"])


def _require_cols(df: pd.DataFrame, name: str, cols):
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required column(s): {', '.join(missing)}")


def linearity_confirm(log_file_path: Optional[str] = None):
    while True:
        corr_sel = input(
            "\nChoose which linearity corrections to apply:\n"
            "1) Do not apply any correction.\n"
            "2) Apply correction for Nitrogen.\n"
            "3) Apply correction for Carbon.\n"
            "4) Apply correction for both Nitrogen and Carbon.\n"
            "Enter 1, 2, 3 or 4:\n"
        ).strip()
        if corr_sel == "1":
            _append_to_log(log_file_path, "User choice: no EA linearity correction.")
            return False, False
        if corr_sel == "2":
            _append_to_log(log_file_path, "User choice: apply EA linearity correction for Nitrogen.")
            return True, False
        if corr_sel == "3":
            _append_to_log(log_file_path, "User choice: apply EA linearity correction for Carbon.")
            return False, True
        if corr_sel == "4":
            _append_to_log(log_file_path, "User choice: apply EA linearity correction for Nitrogen and Carbon.")
            return True, True
        print("Wrong input. Try again.")


def _normalize_name(value) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _get_area_column(df: pd.DataFrame) -> str:
    for col in ("Area All", "area", "Peak Area", "Peak area"):
        if col in df.columns:
            return col
    raise ValueError("No peak-area column found. Expected one of: Area All, area, Peak Area, Peak area.")


def _get_standards_root() -> Path:
    package_root = Path(__file__).resolve().parent.parent
    preferred = package_root / "standards_manager" / "EA" / "data"
    fallback = package_root / "standards_manager" / "data"
    return preferred if preferred.exists() else fallback


def get_ea_linearity_standards_path(tag: str) -> Path:
    root = _get_standards_root()
    if "N" in tag.upper():
        filename = "ea_RS_dN.csv"
    elif "C" in tag.upper():
        filename = "ea_RS_dC.csv"
    else:
        raise ValueError(f"Unknown isotope tag: {tag}")

    path = root / filename
    if not path.exists():
        raise FileNotFoundError(f"EA linearity standards file not found: {path}")
    return path


def load_ea_linearity_metadata(tag: str, log_file_path: Optional[str] = None) -> pd.DataFrame:
    path = get_ea_linearity_standards_path(tag)
    df = pd.read_csv(path)

    if "Linearity Column" in df.columns:
        mask = _truthy_mask(df["Linearity Column"])
    elif "Linearity Standard" in df.columns:
        mask = _truthy_mask(df["Linearity Standard"])
        _append_to_log(log_file_path, f"- Using 'Linearity Standard' as the EA linearity selector for {path.name}.")
    elif "Amount (linearity)" in df.columns:
        mask = df["Amount (linearity)"].notna() & (df["Amount (linearity)"].astype(str).str.strip() != "")
        _append_to_log(log_file_path, f"- Falling back to non-empty 'Amount (linearity)' rows as EA linearity selectors for {path.name}.")
    else:
        raise ValueError(
            f"{path.name} does not contain 'Linearity Column', 'Linearity Standard', or 'Amount (linearity)'."
        )

    selected = df.loc[mask].copy()
    if selected.empty:
        raise ValueError(f"No standards marked for EA linearity correction in {path.name}.")
    name_col = _get_metadata_name_column(selected)
    selected_names = selected[name_col].dropna().astype(str).str.strip()
    selected_names = selected_names[selected_names != ""].tolist()
    _append_to_log(
        log_file_path,
        f"- EA linearity standards selected from {path.name}: {', '.join(selected_names)}.",
    )
    return selected


def _get_metadata_name_column(df: pd.DataFrame) -> str:
    for col in ("EA-IRMS Standards", "ID", "Identifier 1", "Name"):
        if col in df.columns:
            return col
    return df.columns[0]


def prepare_ea_linearity_standards(
    df: pd.DataFrame,
    standards_meta: pd.DataFrame,
    *,
    identifier_column: str,
    target_column: str,
    area_column: str,
    element_type: str,
    component: str,
    log_file_path: Optional[str] = None,
) -> pd.DataFrame:
    _require_cols(df, "ea_df", [identifier_column, target_column, area_column, "Component", "Element Type"])
    name_col = _get_metadata_name_column(standards_meta)

    component_df = df[
        (df["Element Type"] == element_type) &
        (df["Component"] == component) &
        df[target_column].notna() &
        df[area_column].notna()
    ].copy()
    if component_df.empty:
        raise ValueError(f"No valid {element_type} observations with both isotope values and peak areas were found.")

    component_df["_norm_name"] = component_df[identifier_column].map(_normalize_name)
    standards_meta = standards_meta.copy()
    standards_meta["_norm_name"] = standards_meta[name_col].map(_normalize_name)

    matched_parts = []
    missing = []
    for _, row in standards_meta.iterrows():
        norm_name = row["_norm_name"]
        sub = component_df[component_df["_norm_name"] == norm_name].copy()
        if sub.empty:
            sub = component_df[
                component_df["_norm_name"].apply(
                    lambda value: norm_name in value or value in norm_name
                )
            ].copy()
        if sub.empty:
            missing.append(str(row[name_col]))
            continue
        sub["linearity_standard_name"] = row[name_col]
        matched_parts.append(sub)

    if missing:
        raise ValueError(
            f"No input-file observations matched the selected {element_type} linearity standards: {', '.join(missing)}"
        )

    lin_std = pd.concat(matched_parts, ignore_index=False)
    if lin_std[area_column].nunique(dropna=True) < 2:
        raise ValueError(f"The selected {element_type} linearity standards do not span usable peak areas.")
    if len(lin_std) < 3:
        raise ValueError(f"Too few valid {element_type} linearity-standard observations to fit a model (n={len(lin_std)}).")

    _append_to_log(
        log_file_path,
        f"- Matched {len(lin_std)} {element_type} linearity-standard observations across {lin_std['linearity_standard_name'].nunique()} standards.",
    )
    return lin_std


def determine_linearity_model(norm: pd.DataFrame, area_column: str, target_column: str, include_parabolic: bool = False):
    x = norm[area_column].to_numpy(float)
    y = norm[target_column].to_numpy(float)
    best_model, popt, sse, pcov = fit_and_select_best(x, y, include_parabolic=include_parabolic)
    return best_model, popt, pcov, sse, x, y


def generate_model_curve(x, best_model, popt):
    if best_model == "linear":
        return linear_func(x, *popt)
    if best_model == "decay":
        return exp_decay(x, *popt)
    if best_model == "growth":
        return exp_growth(x, *popt)
    if best_model == "parabolic":
        return parabolic_func(x, *popt)
    raise ValueError(f"Unknown linearity model: {best_model}")


def build_linearity_model(
    standards_df: pd.DataFrame,
    *,
    target_column: str,
    area_column: str,
    group_column: str,
    log_file_path: Optional[str] = None,
    include_parabolic: bool = False,
) -> Dict[str, object]:
    _require_cols(standards_df, "standards_df", [target_column, area_column, group_column])
    norm, norm_meta = build_norm_with_fallback(
        standards_df,
        y_col=target_column,
        area_col=area_column,
        group_col=group_column,
        name_col="Identifier 1",
        log_file_path=log_file_path,
    )
    best_model, popt, pcov, sse, x, y = determine_linearity_model(norm, area_column, target_column, include_parabolic)
    y_fit = generate_model_curve(x, best_model, popt)
    tss = np.sum((y - y.mean()) ** 2)
    r_squared = 1.0 if tss == 0 else 1 - (sse / tss)

    lin_top_sort = norm.sort_values(by=area_column, ascending=False)
    top_count = max(int(len(lin_top_sort) * 0.2), 1)
    lin_reference = lin_top_sort.head(top_count)[target_column].median()

    _append_to_log(log_file_path, f"- Best fit EA linearity model type: {best_model}")
    _append_to_log(log_file_path, f"- EA linearity model parameters: {popt}")
    _append_to_log(log_file_path, f"- EA linearity model stats: R² = {r_squared:.3f} | SSE = {sse:.3f}")

    return {
        "best_model": best_model,
        "popt": popt,
        "pcov": pcov,
        "sse": float(sse),
        "r_squared": float(r_squared),
        "lin_reference": float(lin_reference),
        "norm": norm,
        "norm_meta": norm_meta,
        "target_column": target_column,
        "area_column": area_column,
        "group_column": group_column,
    }


def apply_linearity_model(
    df: pd.DataFrame,
    model_info: Dict[str, object],
    *,
    input_column: str,
    area_column: str,
    output_column: str,
    error_column: str,
    row_mask=None,
    log_file_path: Optional[str] = None,
    label: str = "data",
) -> pd.DataFrame:
    _require_cols(df, label, [input_column, area_column])
    out = df.copy()
    mask = pd.Series(True, index=out.index) if row_mask is None else pd.Series(row_mask, index=out.index).fillna(False)
    if not mask.any():
        out[output_column] = out.get(output_column, np.nan)
        out[error_column] = out.get(error_column, np.nan)
        _append_to_log(log_file_path, f"- Linearity correction applied to {label}: n=0")
        return out

    x = out.loc[mask, area_column].to_numpy(float)
    yhat = generate_model_curve(x, model_info["best_model"], model_info["popt"])
    offset = float(model_info["lin_reference"]) - yhat

    out.loc[mask, output_column] = out.loc[mask, input_column].to_numpy(float) + offset
    out.loc[mask, error_column] = prediction_std(model_info["best_model"], x, model_info["popt"], model_info["pcov"], nsigma=2)
    _append_to_log(log_file_path, f"- Linearity correction applied to {label}: n={int(mask.sum())}")
    return out


def plot_linearity_diagnostics(
    model_info: Dict[str, object],
    *,
    fig_path: Optional[str],
    figure_name: str,
    y_label: str,
    title: str,
):
    if not fig_path:
        return

    norm = model_info["norm"].copy()
    area_col = model_info["area_column"]
    target_col = model_info["target_column"]
    x = norm[area_col].to_numpy(float)
    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = generate_model_curve(x_fit, model_info["best_model"], model_info["popt"])

    fig, ax = plt.subplots(figsize=(6, 4))
    for group, sub in norm.groupby(model_info["group_column"]):
        ax.scatter(sub[area_col], sub[target_col], alpha=0.6, ec="k", label=str(group))

    ax.plot(x_fit, y_fit, "r-", label=model_info["best_model"])
    ax.set_xlabel("Peak Area")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    plt.savefig(os.path.join(fig_path, figure_name), bbox_inches="tight")
    plt.show()
