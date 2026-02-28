import numpy as np
import os
import pandas as pd
from .figures import total_dD_correction_plot
from .config import CorrectionConfig
from .base_functions import append_to_log

def combine_errors(df, config, log_file_path, isotope):

    error_terms = []

    ref_err_col = f"RS error"
    if ref_err_col in df.columns:
        error_terms.append(df[ref_err_col])

    if config.linearity_applied:
        col = f"Linearity {isotope} error"
        if col in df.columns:
            error_terms.append(df[col])

    if config.drift_applied:
        col = f"Drift {isotope} error"
        if col in df.columns:
            error_terms.append(df[col])

    if f"Methanol error {isotope}" in df.columns:
        error_terms.append(df[f"Methanol error {isotope}"])

    if error_terms:
        df["Total Error"] = np.sqrt(sum(err.astype(float)**2 for err in error_terms))
    else:
        df["Total Error"] = np.nan

    if f"Corrected {isotope}" in df.columns:
        summary = df.groupby("Component").agg(
            mean_corrected=(f"Corrected {isotope}", "mean"),
            mean_total_error=("Total Error", "mean")
        ).reset_index()

        append_to_log(log_file_path, "Summary of standards")
        append_to_log(log_file_path, summary)

def mean_values_with_uncertainty(data, cfg, iso,
                                 sample_name_header="Identifier 1",
                                 chain_header="chain"):

    def pick_first_existing(cols, candidates):
        for c in candidates:
            if c in cols:
                return c
        return None

    cols = set(data.columns)

    if iso not in cols:
        raise KeyError(f"Expected raw isotope column '{iso}' not found.")

    iso_col = iso

    drift_corr_col = f"drift_corrected_{iso}"
    lin_corr_col   = f"linearity_corrected_{iso}"
    rs_col         = f"RS_{iso}"
    meth_col       = f"methanol_{iso}"
    pame_col       = f"PAME_methanol_{iso}"

    grouped = data.groupby([sample_name_header, chain_header], dropna=False)

    agg_dict = {}

    if "area" in cols:
        agg_dict["area"] = ["mean"]

    agg_dict[iso_col] = ["mean", "std", "count"]

    if cfg.drift_applied and drift_corr_col in cols:
        agg_dict[drift_corr_col] = ["mean"]

    if cfg.linearity_applied and lin_corr_col in cols:
        agg_dict[lin_corr_col] = ["mean"]

    if "drift_error" in cols:
        agg_dict["drift_error"] = ["mean"]

    if "linearity_error" in cols:
        agg_dict["linearity_error"] = ["mean"]

    if rs_col in cols:
        agg_dict[rs_col] = ["mean", "std"]

    if "RS_error" in cols:
        agg_dict["RS_error"] = ["mean"]

    if meth_col in cols:
        agg_dict[meth_col] = ["mean", "std"]

    if "methanol_error" in cols:
        agg_dict["methanol_error"] = ["mean"]

    if pame_col in cols:
        agg_dict[pame_col] = ["mean", "std"]

    stats = grouped.agg(agg_dict).reset_index()
    stats.columns = ["_".join(filter(None, c)).strip("_")
                     for c in stats.columns.values]

    # SEM
    n_col = f"{iso}_count"
    std_col = pick_first_existing(
        set(stats.columns),
        [f"{rs_col}_std", f"{meth_col}_std", f"{iso}_std"]
    )

    sem_name = f"replicate_{iso}_sem"

    if std_col and n_col in stats.columns:
        stats[sem_name] = stats[std_col] / np.sqrt(
            stats[n_col].replace(0, np.nan)
        )
    else:
        stats[sem_name] = np.nan

    # Total uncertainty
    error_cols = []

    for col in ["drift_error_mean",
                "linearity_error_mean",
                "RS_error_mean",
                "methanol_error_mean"]:
        if col in stats.columns:
            error_cols.append(col)

    error_cols.append(sem_name)

    if error_cols:
        stats["total_uncertainty"] = np.sqrt(
            sum(stats[c].astype(float)**2 for c in error_cols if c in stats.columns)
        )
    else:
        stats["total_uncertainty"] = np.nan

    return stats

def output_results(raw_unknown, unknown, sd, unknown_pame,
                   folder_path, fig_path, res_path,
                   isotope, pame, log_file_path, cfg):

    rs_label = "RS"

    # Dynamic renaming map
    column_name_mapping = {
        'Identifier 1': 'Sample Name',
        'chain': 'Chain Length',
        'area_mean': 'Mean Area',

        f'{isotope}_mean': f'Average raw {isotope}',
        f'{isotope}_std': f'{isotope} Std Dev',
        f'{isotope}_count': 'Number of Replicates',

        f'drift_corrected_{isotope}_mean': f'Drift Corrected {isotope}',
        'drift_error_mean': 'Drift Error',

        f'linearity_corrected_{isotope}_mean': f'Linearity Corrected {isotope}',
        'linearity_error_mean': 'Linearity Error',

        f'RS_{isotope}_mean': f'{rs_label} Corrected {isotope}',
        'RS_error_mean': f'{rs_label} Error',

        f'methanol_{isotope}_mean': f'Methanol Corrected {isotope}',
        'methanol_error_mean': f'Methanol Error',

        f'replicate_{isotope}_sem': f'Mean replicate std {isotope}',
        'total_uncertainty': 'Total Uncertainty'
    }

    unknown_renamed = unknown.rename(columns=column_name_mapping)

    if pame and unknown_pame is not None:
        unknown_pame_renamed = unknown_pame.rename(columns=column_name_mapping)

    # Determine final corrected column
    priority_cols = [
        f'Methanol Corrected {isotope}',
        f'{rs_label} Corrected {isotope}',
        f'Linearity Corrected {isotope}',
        f'Drift Corrected {isotope}'
    ]

    final_corrected = None
    for col in priority_cols:
        if col in unknown_renamed.columns:
            final_corrected = unknown_renamed[col]
            break

    if final_corrected is None:
        raise KeyError(
            f"No corrected column found for {isotope}. "
            f"Available columns: {list(unknown_renamed.columns)}"
        )

    unknown_renamed.insert(
        len(unknown_renamed.columns),
        f'Corrected {isotope}',
        final_corrected
    )

    unknown_renamed.to_csv(
        os.path.join(res_path, 'Results - sample mean.csv'),
        index=False
    )

    if pame and unknown_pame is not None:
        unknown_pame_renamed.to_csv(
            os.path.join(res_path, 'PAME Results - sample mean.csv'),
            index=False
        )

    # Standards renaming
    standards_map = {
        f'{isotope}': f'Raw {isotope}',
        f'drift_corrected_{isotope}': f'Drift corrected {isotope}',
        'drift_error': f'Drift {isotope} error',
        f'linearity_corrected_{isotope}': f'Linearity corrected {isotope}',
        'linearity_error': f'Linearity {isotope} error',
        f'RS_{isotope}': f'{rs_label} corrected {isotope}',
        'RS_error': f'{rs_label} error',
        f'methanol_{isotope}': f'Methanol corrected {isotope}',
        'total_uncertainty': 'Total uncertainty'
    }

    standards_categorized = sd.rename(columns=standards_map)

    standards_categorized.to_csv(
        os.path.join(res_path, 'Results - standards.csv'),
        index=False
    )

    combine_errors(standards_categorized, cfg, log_file_path, isotope)

    print("\nCorrections complete :)")