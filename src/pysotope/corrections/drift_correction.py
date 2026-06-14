from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


EA_STANDARDS_DIR_CANDIDATES = [
    Path(__file__).resolve().parent.parent / "standards_manager" / "EA" / "data",
    Path(__file__).resolve().parent.parent / "standards_manager" / "data",
]
EA_DRIFT_STANDARD_FILES = {
    "N": "ea_RS_dN.csv",
    "C": "ea_RS_dC.csv",
}


def _append_to_log(log_file_path: Optional[str], log_message):
    if not log_file_path:
        return
    with open(log_file_path, "a", encoding="utf-8", errors="replace") as log_file:
        print(f" {log_message}", file=log_file)


def _pos_response(response: str) -> bool:
    return response.lower() in {"yes", "y", "true", "t", ""}


def _neg_response(response: str) -> bool:
    return response.lower() in {"no", "n", "false", "f"}


def _require_cols(df: pd.DataFrame, name: str, cols):
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required column(s): {', '.join(missing)}")


def q_drift(log_file_path: Optional[str] = None) -> str:
    valid = {"yes", "y", "true", "t", "no", "n", "false", "f"}
    while True:
        choice = input("\nApply drift correction? (Y/N)\n").strip().lower()
        if choice in valid:
            _append_to_log(log_file_path, f"- Drift correction applied: {choice}")
            return choice[0]
        print("\n[ERROR] Invalid response. Please type Y or N.\n")


def drift_confirm(log_file_path: Optional[str] = None):
    while True:
        corr_sel = input(
            "\nChoose which corrections to apply:\n"
            "1) Do not apply any correction.\n"
            "2) Apply correction for Nitrogen.\n"
            "3) Apply correction for Carbon.\n"
            "4) Apply correction for both Nitrogen and Carbon.\n"
            "Enter 1, 2, 3 or 4:\n"
        ).strip()
        if corr_sel == "1":
            _append_to_log(log_file_path, "User choice: no EA drift correction.")
            return False, False
        if corr_sel == "2":
            _append_to_log(log_file_path, "User choice: apply EA drift correction for Nitrogen.")
            return True, False
        if corr_sel == "3":
            _append_to_log(log_file_path, "User choice: apply EA drift correction for Carbon.")
            return False, True
        if corr_sel == "4":
            _append_to_log(log_file_path, "User choice: apply EA drift correction for Nitrogen and Carbon.")
            return True, True
        print("Wrong input. Try again.")


def get_isotope(tag):
    if "N" in tag.upper():
        return "d 15N/14N", r"$\delta^{15}\mathrm{N}$", "Nitrogen", "N2"
    if "C" in tag.upper():
        return "d 13C/12C", r"$\delta^{13}\mathrm{C}$", "Carbon", "CO2"
    raise ValueError(f"Unknown isotope tag: {tag}")


def get_sorghum(df: pd.DataFrame, log_file_path: Optional[str] = None):
    sorghum = df[df["Identifier 1"].astype(str).str.lower() == "sorghum"].copy()
    if sorghum.empty:
        _append_to_log(log_file_path, "Error: No 'SORGHUM' standards found in the data.")
        raise ValueError("No 'SORGHUM' standards found.")

    sorghum_n = sorghum[sorghum["Element Type"] == "Nitrogen"].copy()
    sorghum_c = sorghum[sorghum["Element Type"] == "Carbon"].copy()
    _append_to_log(log_file_path, "'SORGHUM' standards found in the data.")
    return sorghum_n, sorghum_c


def _truthy_mask(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "1.0", "yes", "y", "t"})
    )


def _get_ea_standards_dir() -> Path:
    for path in EA_STANDARDS_DIR_CANDIDATES:
        if path.exists():
            return path
    return EA_STANDARDS_DIR_CANDIDATES[-1]


def _load_ea_drift_metadata(tag: str) -> pd.DataFrame:
    isotope_tag = tag.upper()
    if isotope_tag not in EA_DRIFT_STANDARD_FILES:
        raise ValueError(f"Unknown EA drift isotope tag: {tag}")

    path = _get_ea_standards_dir() / EA_DRIFT_STANDARD_FILES[isotope_tag]
    if not path.exists():
        raise FileNotFoundError(f"EA drift standards file not found: {path}")

    metadata = pd.read_csv(path, encoding="unicode_escape")
    if "ID" not in metadata.columns:
        raise ValueError(f"{path.name} is missing required column: ID")
    if "drift" not in metadata.columns:
        raise ValueError(
            f"{path.name} is missing required column: drift. "
            "Add a 'drift' checkbox column with pysotope.standard_editor()."
        )
    return metadata[_truthy_mask(metadata["drift"])].copy()


def get_ea_drift_standards(df: pd.DataFrame, tag: str, log_file_path: Optional[str] = None):
    metadata = _load_ea_drift_metadata(tag)
    _, _, element, _ = get_isotope(tag)

    if metadata.empty:
        raise ValueError(f"No {element} EA standards are selected for drift correction.")

    ids = metadata["ID"].dropna().astype(str).str.strip()
    ids = ids[ids != ""].tolist()
    if not ids:
        raise ValueError(f"No valid {element} EA drift standard IDs found.")

    id_values = df["Identifier 1"].astype(str).str.strip().str.lower()
    selected_ids = {std_id.lower() for std_id in ids}
    drift_standards = df[
        id_values.isin(selected_ids) & (df["Element Type"] == element)
    ].copy()

    if drift_standards.empty:
        msg = f"No selected {element} EA drift standards found in the data: {', '.join(ids)}."
        _append_to_log(log_file_path, f"Error: {msg}")
        raise ValueError(msg)

    _append_to_log(log_file_path, f"Selected EA {element} drift standard(s): {', '.join(ids)}.")
    return drift_standards


def _standard_name_column(df: pd.DataFrame) -> str:
    for col in ("EA-IRMS Standards", "ID", "Identifier 1", "Name"):
        if col in df.columns:
            return col
    raise ValueError("EA standards table does not contain a recognizable standard-name column.")


def _normalized_names(values) -> dict[str, str]:
    names = {}
    for value in values:
        if pd.isna(value):
            continue
        name = str(value).strip()
        if name:
            names[name.lower()] = name
    return names


def _parse_indices(idx_str: str) -> list[int]:
    return [int(idx.strip()) for idx in idx_str.split(",") if idx.strip()]


def review_ea_reference_standards(
    df: pd.DataFrame,
    carbon_standards: pd.DataFrame,
    nitrogen_standards: pd.DataFrame,
    log_file_path: Optional[str] = None,
):
    while True:
        choice = input("Check and/or remove any standard measurements? [Y/N]\n").strip().lower()
        if choice in {"yes", "y", "true", "t"}:
            break
        if choice in {"no", "n", "false", "f", ""}:
            _append_to_log(log_file_path, "No reference standard measurements removed")
            return
        print("Wrong input. Please enter Y or N.")

    removed_rows = []
    for tag, standards in (("N", nitrogen_standards), ("C", carbon_standards)):
        col, _, element, _ = get_isotope(tag)
        name_col = _standard_name_column(standards)
        standard_names = _normalized_names(standards[name_col])

        present_keys = (
            df.loc[df["Element Type"] == element, "Identifier 1"]
            .astype(str)
            .str.strip()
            .str.lower()
            .unique()
        )
        present = [standard_names[key] for key in standard_names if key in present_keys]

        for standard_name in present:
            mask = (
                (df["Element Type"] == element) &
                (df["Identifier 1"].astype(str).str.strip().str.lower() == standard_name.lower())
            )
            sub = df.loc[mask, ["Identifier 1", "Time", col]].copy()
            if sub.empty:
                continue

            print("")
            print(f"{element} reference standard: {standard_name}")
            print(sub)
            idx_str = input(
                "Enter indices to remove in a comma-separated manner. "
                "(Enter None if no values are to be removed)\n"
            ).strip()
            if not idx_str or idx_str.lower() == "none":
                continue

            try:
                idxs = _parse_indices(idx_str)
            except ValueError:
                print("Invalid index list. Please use comma-separated integer index values.")
                continue

            valid_idxs = [idx for idx in idxs if idx in df.index and idx in sub.index]
            invalid_idxs = sorted(set(idxs).difference(valid_idxs))
            if invalid_idxs:
                print(f"Ignoring invalid index value(s): {invalid_idxs}")

            if valid_idxs:
                removed = df.loc[valid_idxs, ["Identifier 1", "Time", col]].copy()
                removed.insert(0, "Element Type", element)
                removed_rows.append(removed)
                df.drop(valid_idxs, inplace=True, errors="ignore")

    if removed_rows:
        removed_df = pd.concat(removed_rows).sort_index()
        _append_to_log(log_file_path, "Removed reference standard runs")
        _append_to_log(log_file_path, removed_df)
    else:
        _append_to_log(log_file_path, "Removed reference standard runs: None")


def drop_standards(std_df: pd.DataFrame, df: pd.DataFrame, tag: str, log_file_path: Optional[str] = None):
    col, _, el, _ = get_isotope(tag)
    sub = std_df[["Identifier 1", "Time", col]].copy()
    print("")
    print(sub)
    idx_str = input(
        "Enter indices that you wish to remove in a comma-separated manner. "
        "(Enter None if no values are to be removed)\n"
    ).strip()
    if idx_str and idx_str.lower() != "none":
        idxs = [int(idx.strip()) for idx in idx_str.split(",") if idx.strip()]
        std_df.drop(idxs, inplace=True, errors="ignore")
        df.drop(idxs, inplace=True, errors="ignore")
        msg = f"{el} values at {idx_str} have been removed."
        print(msg)
        _append_to_log(log_file_path, msg)


def _normalize_standards(
    standards_df: pd.DataFrame,
    target_column: str,
    index_column: str,
    weight_column: Optional[str] = None,
    group_column: Optional[str] = None,
) -> pd.DataFrame:
    cols = [target_column, index_column]
    if weight_column:
        cols.append(weight_column)
    if group_column:
        cols.append(group_column)

    norm = standards_df[cols].copy()
    if group_column:
        norm[target_column] = norm.groupby(group_column)[target_column].transform(lambda s: s - s.mean())
    else:
        norm[target_column] = norm[target_column] - norm[target_column].mean()
    return norm


def build_drift_model(
    standards_df: pd.DataFrame,
    *,
    target_column: str,
    index_column: str,
    weight_column: Optional[str] = None,
    group_column: Optional[str] = None,
    log_file_path: Optional[str] = None,
) -> Dict[str, object]:
    required = [target_column, index_column]
    if weight_column:
        required.append(weight_column)
    if group_column:
        required.append(group_column)
    _require_cols(standards_df, "standards_df", required)

    norm = _normalize_standards(
        standards_df,
        target_column=target_column,
        index_column=index_column,
        weight_column=weight_column,
        group_column=group_column,
    )
    index_mean = norm[index_column].mean()
    norm["t_ctr"] = norm[index_column] - index_mean

    exog = sm.add_constant(norm["t_ctr"], has_constant="add")
    y = norm[target_column].astype(float).to_numpy()

    if weight_column:
        weights = norm[weight_column].astype(float).to_numpy()
        max_weight = np.nanmax(weights)
        if not np.isfinite(max_weight) or max_weight == 0:
            weights = np.ones_like(weights, dtype=float)
        else:
            weights = weights / max_weight
        model = sm.WLS(y, exog, weights=weights).fit()
    else:
        model = sm.OLS(y, exog).fit()

    intercept, slope = model.params
    p_value = float(model.pvalues.iloc[1]) if len(model.pvalues) > 1 else np.nan
    std_error = float(model.bse.iloc[1]) if len(model.bse) > 1 else np.nan
    r_squared = float(getattr(model, "rsquared_adj", getattr(model, "rsquared", np.nan)))

    norm["corrected_norm"] = norm[target_column] - slope * norm["t_ctr"]

    _append_to_log(
        log_file_path,
        f"- Drift model fit completed for {target_column}: b0={intercept:.4f}, b1={slope:.4f}, R²={r_squared:.3f}",
    )

    return {
        "model": model,
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": r_squared,
        "p_value": p_value,
        "std_error": std_error,
        "index_mean": float(index_mean),
        "target_column": target_column,
        "index_column": index_column,
        "weight_column": weight_column,
        "group_column": group_column,
        "normalized_standards": norm,
    }


def apply_drift_model(
    df: pd.DataFrame,
    model_info: Dict[str, object],
    *,
    target_column: str,
    index_column: str,
    output_column: str,
    error_column: str,
    row_mask=None,
    log_file_path: Optional[str] = None,
    label: str = "data",
    error_nsigma: float = 2.0,
) -> pd.DataFrame:
    _require_cols(df, label, [target_column, index_column])
    out = df.copy()
    mask = pd.Series(True, index=out.index) if row_mask is None else pd.Series(row_mask, index=out.index).fillna(False)
    if not mask.any():
        out[output_column] = out.get(output_column, np.nan)
        out[error_column] = out.get(error_column, np.nan)
        _append_to_log(log_file_path, f"- Drift correction applied to {label}: n=0")
        return out

    t_ctr = out.loc[mask, index_column].astype(float) - float(model_info["index_mean"])
    exog = sm.add_constant(t_ctr, has_constant="add")
    prediction = model_info["model"].get_prediction(exog)
    pred_error = prediction.se_obs * error_nsigma

    out.loc[mask, output_column] = out.loc[mask, target_column].astype(float) - float(model_info["slope"]) * t_ctr
    out.loc[mask, error_column] = pred_error
    _append_to_log(log_file_path, f"- Drift correction applied to {label}: n={int(mask.sum())}")
    return out


def plot_drift_diagnostics(
    model_info: Dict[str, object],
    *,
    fig_path: Optional[str],
    figure_name: str,
    y_label: str,
    x_label: str,
):
    if not fig_path:
        return

    norm = model_info["normalized_standards"].copy()
    fig, ax = plt.subplots(figsize=(6, 4))

    group_column = model_info.get("group_column")
    if group_column and group_column in norm.columns:
        markers = ["o", "s", "^", "D", "v", "P", "X"]
        for idx, group in enumerate(norm[group_column].dropna().unique()):
            marker = markers[idx % len(markers)]
            temp = norm[norm[group_column] == group]
            ax.scatter(temp[model_info["index_column"]], temp[model_info["target_column"]], marker=marker, c="k", ec="k", alpha=0.5, label=f"{group} Raw")
            ax.scatter(temp[model_info["index_column"]], temp["corrected_norm"], marker=marker, c="red", ec="k", alpha=0.7, label=f"{group} Corrected")
    else:
        ax.scatter(norm[model_info["index_column"]], norm[model_info["target_column"]], c="k", ec="k", alpha=0.5, label="Raw")
        ax.scatter(norm[model_info["index_column"]], norm["corrected_norm"], c="red", ec="k", alpha=0.7, label="Corrected")

    t_line = np.linspace(norm[model_info["index_column"]].min(), norm[model_info["index_column"]].max(), 50)
    y_line = model_info["slope"] * (t_line - model_info["index_mean"]) + model_info["intercept"]
    ax.plot(t_line, y_line, "--k", lw=1, label="Drift fit")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()
    fig.tight_layout()
    plt.savefig(os.path.join(fig_path, figure_name), bbox_inches="tight")
    plt.show()
