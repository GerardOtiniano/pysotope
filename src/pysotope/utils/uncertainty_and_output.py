import numpy as np
import os
import pandas as pd
from .figures import total_dD_correction_plot
from .config import CorrectionConfig
from .base_functions import append_to_log

# def output_results(raw_unknown, unknown, sd, unknown_pame, folder_path, fig_path, res_path, isotope, pame, log_file_path, cfg):
#     # Reference scale label for OUTPUT headers only
#     ref_label = "VSMOW" if isotope == "dD" else "PDB"

#     # Define column name mappings for output CSV files
#     column_name_mapping = {
#         'Identifier 1': 'Sample Name',
#         'chain': 'Chain Length',
#         'area_mean': 'Mean Area',

#         # raw replicate stats (still keyed as dD_* in your current stats output)
#         'dD_mean': f'Average raw {isotope}',
#         'dD_std': f'{isotope} Std Dev',
#         'dD_count': 'Number of Replicates',

#         'drift_corrected_dD_mean': f'Drift Corrected {isotope}',
#         'drift_error_mean': 'Drift Error',
#         'linearity_corrected_dD_mean': f'Linearity Corrected {isotope}',
#         'linearity_error_mean': 'Linearity Error',

#         # *** OUTPUT LABEL changes here ***
#         'VSMOW_dD_mean': f'{ref_label} Corrected {isotope}',
#         'VSMOW_error_mean': f'{ref_label} Error',

#         'methanol_dD_mean': f'Methanol Corrected {isotope}',
#         'methanol_error_mean': f'Methanol Error',

#         'PAME_methanol_dD_mean': f'PAME Methanol Calculated {isotope}',
#         'PAME_methanol_dD_std': f'PAME Methanol Calculated Error {isotope}',

#         'replicate_dD_sem': f'Mean replicate std {isotope}',
#         'total_uncertainty': 'Total Uncertainty'
#     }

#     # Rename columns in unknown dataframe
#     unknown_renamed = unknown.rename(columns=column_name_mapping)
#     if pame:
#         unknown_pame_renamed = unknown_pame.rename(columns=column_name_mapping)

#     if isotope == 'dD':
    
#         meth_col = 'Methanol Corrected ' + str(isotope)
    
#         if meth_col in unknown_renamed.columns:
#             # Methanol correction was applied
#             corrected_values = unknown_renamed[meth_col]
#         else:
#             # Methanol correction was skipped
#             # Use the pre-methanol corrected isotope value instead
#             corrected_values = unknown_renamed[str(isotope)]
    
#         unknown_renamed.insert(
#             len(unknown_renamed.columns) - 1,
#             'Corrected ' + str(isotope),
#             corrected_values)

#     # Save unknown dataframe to CSV
#     unknown_renamed.reset_index(drop=True, inplace=True)
#     unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample mean.csv'), index=False)
#     if pame:
#         unknown_pame_renamed.reset_index(drop=True, inplace=True)
#         unknown_pame_renamed.to_csv(os.path.join(res_path, 'PAME Results - sample mean.csv'), index=False)

#     # Select columns for standards dataframe (internal keys unchanged)
#     columns_to_select = [
#         "Date", "Time", "Identifier 1", "chain", "Rt", "area", "Area 2", "Area 3",
#         "Ampl  2", "Ampl  3", "BGD 2", 'BGD 3', "time_rel", "dD", "drift_corrected_dD",
#         "drift_error", "linearity_corrected_dD", "linearity_error", "VSMOW_dD", "VSMOW_error",
#         "total_uncertainty"
#     ]

#     existing_columns = [col for col in columns_to_select if col in sd.columns]
#     standards_selected = sd[existing_columns].copy()

#     # Filter and categorize standards
#     lin_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C20|C28")].copy()
#     lin_std_temp['Standard Type'] = "Linearity"

#     drift_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C218|C24")].copy()
#     drift_std_temp['Standard Type'] = "Drift"

#     standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

#     # Rename columns in standards dataframe (OUTPUT LABEL changes here too)
#     column_rename_map = {
#         "Identifier 1": "Standard ID",
#         "chain": "Component",
#         "Rt": "Retention time",
#         "area": "Peak area",
#         "Area 2": "Peak area 2",
#         "Area 3": "Peak area 3",
#         "Ampl  2": "Amplitude 2",
#         "Ampl  3": "Amplitude 3",
#         "BGD 2": "Background 2",
#         "BGD 3": "Background 3",
#         "time_rel": "Time relative",

#         "dD": f"Raw {isotope}",
#         "drift_corrected_dD": f"Drift corrected {isotope}",
#         "drift_error": f"Drift {isotope} error",
#         "linearity_corrected_dD": f"Linearity {isotope}",
#         "linearity_error": f"Linearity {isotope} error",

#         # *** OUTPUT LABEL changes here ***
#         "VSMOW_dD": f"{ref_label} corrected {isotope}",
#         "VSMOW_error": f"{ref_label} error",

#         "methanol_dD": f"Methanol corrected {isotope}",
#         "methanol_error": f"Methanol error {isotope}",
#         "total_uncertainty": "Total uncertainty"
#     }

#     standards_categorized = standards_categorized.rename(columns=column_rename_map)

#     # Insert corrected column using ref_label
#     standards_categorized.insert(
#         len(standards_categorized.columns) - 1,
#         'Corrected ' + str(isotope),
#         standards_categorized[f"{ref_label} corrected {isotope}"]
#     )

#     standards_categorized.to_csv(os.path.join(res_path, 'Results - standards.csv'), index=False)

#     # Rename columns in raw_unknown dataframe
#     raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

#     if 'total_error' in raw_unknown.columns:
#         raw_unknown_renamed = raw_unknown_renamed.drop('total_error', axis=1)

#     raw_unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample replicates.csv'), index=False)

#     total_dD_correction_plot(raw_unknown_renamed, unknown_renamed, folder_path, fig_path, isotope)
#     combine_errors(standards_categorized, cfg, log_file_path, isotope)
#     print("\nCorrections complete :)")
    
def combine_errors(df, config, log_file_path, isotope):

    ref_label = "VSMOW" if isotope == "dD" else "PDB"

    error_terms = []

    ref_err_col = f"{ref_label} error"
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
        df["Total Error"] = np.sqrt(sum(err**2 for err in error_terms))
    else:
        df["Total Error"] = np.nan

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

    drift_corr_col = pick_first_existing(cols,
        [f"drift_corrected_{iso}", "drift_corrected_dD", "drift_corrected_dC"])

    lin_corr_col = pick_first_existing(cols,
        [f"linearity_corrected_{iso}", "linearity_corrected_dD", "linearity_corrected_dC"])

    drift_err_col = pick_first_existing(cols, ["drift_error"])
    lin_err_col = pick_first_existing(cols, ["linearity_error"])

    vsmow_col = pick_first_existing(cols,
        [f"VSMOW_{iso}", "VSMOW_dD", "VSMOW_dC"])

    vsmow_err_col = pick_first_existing(cols, ["VSMOW_error"])

    # FIX: methanol now isotope-aware
    meth_col = pick_first_existing(cols,
        [f"methanol_{iso}", "methanol_dD", "methanol_dC"])

    meth_err_col = pick_first_existing(cols, ["methanol_error"])

    pame_col = pick_first_existing(cols,
        [f"PAME_methanol_{iso}", "PAME_methanol_dD"])

    grouped = data.groupby([sample_name_header, chain_header], dropna=False)

    agg_dict = {}

    if "area" in cols:
        agg_dict["area"] = ["mean"]

    agg_dict[iso_col] = ["mean", "std", "count"]

    if getattr(cfg, "drift_applied", False) and drift_corr_col:
        agg_dict[drift_corr_col] = ["mean"]
    if getattr(cfg, "linearity_applied", False) and lin_corr_col:
        agg_dict[lin_corr_col] = ["mean"]
    if drift_err_col:
        agg_dict[drift_err_col] = ["mean"]
    if lin_err_col:
        agg_dict[lin_err_col] = ["mean"]

    if vsmow_col:
        agg_dict[vsmow_col] = ["mean", "std"]
    if vsmow_err_col:
        agg_dict[vsmow_err_col] = ["mean"]

    if meth_col:
        agg_dict[meth_col] = ["mean", "std"]
    if meth_err_col:
        agg_dict[meth_err_col] = ["mean"]
    if pame_col:
        agg_dict[pame_col] = ["mean", "std"]

    stats = grouped.agg(agg_dict).reset_index()
    stats.columns = ["_".join(filter(None, c)).strip("_")
                     for c in stats.columns.values]

    # SEM logic
    n_col = f"{iso_col}_count"

    std_candidates = []
    if meth_col:
        std_candidates.append(f"{meth_col}_std")
    if pame_col:
        std_candidates.append(f"{pame_col}_std")
    if vsmow_col:
        std_candidates.append(f"{vsmow_col}_std")
    std_candidates.append(f"{iso_col}_std")

    sem_std_col = pick_first_existing(set(stats.columns), std_candidates)

    sem_name = f"replicate_{iso}_sem"

    if sem_std_col and n_col in stats.columns:
        stats[sem_name] = stats[sem_std_col] / np.sqrt(
            stats[n_col].replace(0, np.nan)
        )
    else:
        stats[sem_name] = np.nan

    # Total uncertainty
    error_cols = []

    for err in ["drift_error_mean",
                "linearity_error_mean",
                "VSMOW_error_mean",
                "methanol_error_mean"]:
        if err in stats.columns:
            error_cols.append(err)

    error_cols.append(sem_name)

    if error_cols:
        sq = sum(stats[c].astype(float)**2 for c in error_cols if c in stats.columns)
        stats["total_uncertainty"] = np.sqrt(sq)
    else:
        stats["total_uncertainty"] = np.nan

    return stats

def output_results(raw_unknown, unknown, sd, unknown_pame,
                   folder_path, fig_path, res_path,
                   isotope, pame, log_file_path, cfg):

    # Reference scale label for OUTPUT headers only
    ref_label = "VSMOW" if isotope == "dD" else "PDB"

    # ---------------------------------------------
    # Column name mapping (output labels only)
    # ---------------------------------------------
    column_name_mapping = {
        'Identifier 1': 'Sample Name',
        'chain': 'Chain Length',
        'area_mean': 'Mean Area',

        'dD_mean': f'Average raw {isotope}',
        'dD_std': f'{isotope} Std Dev',
        'dD_count': 'Number of Replicates',

        'drift_corrected_dD_mean': f'Drift Corrected {isotope}',
        'drift_error_mean': 'Drift Error',
        'linearity_corrected_dD_mean': f'Linearity Corrected {isotope}',
        'linearity_error_mean': 'Linearity Error',

        'VSMOW_dD_mean': f'{ref_label} Corrected {isotope}',
        'VSMOW_error_mean': f'{ref_label} Error',

        'methanol_dD_mean': f'Methanol Corrected {isotope}',
        'methanol_error_mean': f'Methanol Error',

        'PAME_methanol_dD_mean': f'PAME Methanol Calculated {isotope}',
        'PAME_methanol_dD_std': f'PAME Methanol Calculated Error {isotope}',

        'replicate_dD_sem': f'Mean replicate std {isotope}',
        'total_uncertainty': 'Total Uncertainty'
    }

    unknown_renamed = unknown.rename(columns=column_name_mapping)

    if pame:
        unknown_pame_renamed = unknown_pame.rename(columns=column_name_mapping)

    # -------------------------------------------------
    # Determine which column is the FINAL corrected one
    # -------------------------------------------------

    meth_col = f'Methanol Corrected {isotope}'
    ref_col = f'{ref_label} Corrected {isotope}'

    if meth_col in unknown_renamed.columns:
        final_corrected = unknown_renamed[meth_col]
    elif ref_col in unknown_renamed.columns:
        final_corrected = unknown_renamed[ref_col]
    else:
        raise KeyError(
            f"No corrected column found for {isotope}. "
            f"Expected '{meth_col}' or '{ref_col}'. "
            f"Available columns: {list(unknown_renamed.columns)}"
        )

    unknown_renamed.insert(
        len(unknown_renamed.columns) - 1,
        f'Corrected {isotope}',
        final_corrected
    )

    # -------------------------------------------------
    # Save unknown sample means
    # -------------------------------------------------

    unknown_renamed.reset_index(drop=True, inplace=True)
    unknown_renamed.to_csv(
        os.path.join(res_path, 'Results - sample mean.csv'),
        index=False
    )

    if pame:
        unknown_pame_renamed.reset_index(drop=True, inplace=True)
        unknown_pame_renamed.to_csv(
            os.path.join(res_path, 'PAME Results - sample mean.csv'),
            index=False
        )

    # -------------------------------------------------
    # Standards processing
    # -------------------------------------------------

    columns_to_select = [
        "Date", "Time", "Identifier 1", "chain", "Rt", "area",
        "Area 2", "Area 3", "Ampl  2", "Ampl  3",
        "BGD 2", 'BGD 3', "time_rel", "dD",
        "drift_corrected_dD", "drift_error",
        "linearity_corrected_dD", "linearity_error",
        "VSMOW_dD", "VSMOW_error",
        "methanol_dD",
        "total_uncertainty"
    ]

    existing_columns = [c for c in columns_to_select if c in sd.columns]
    standards_selected = sd[existing_columns].copy()

    lin_std_temp = standards_selected[
        standards_selected["Identifier 1"].str.contains("C20|C28", na=False)
    ].copy()
    lin_std_temp['Standard Type'] = "Linearity"

    drift_std_temp = standards_selected[
        standards_selected["Identifier 1"].str.contains("C218|C24", na=False)
    ].copy()
    drift_std_temp['Standard Type'] = "Drift"

    standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

    # Rename standards columns
    column_rename_map = {
        "Identifier 1": "Standard ID",
        "chain": "Component",
        "Rt": "Retention time",
        "area": "Peak area",

        "dD": f"Raw {isotope}",
        "drift_corrected_dD": f"Drift corrected {isotope}",
        "drift_error": f"Drift {isotope} error",
        "linearity_corrected_dD": f"Linearity {isotope}",
        "linearity_error": f"Linearity {isotope} error",

        "VSMOW_dD": f"{ref_label} corrected {isotope}",
        "VSMOW_error": f"{ref_label} error",

        "methanol_dD": f"Methanol corrected {isotope}",
        "total_uncertainty": "Total uncertainty"
    }

    standards_categorized = standards_categorized.rename(columns=column_rename_map)

    # Determine corrected standard column
    std_meth_col = f"Methanol corrected {isotope}"
    std_ref_col = f"{ref_label} corrected {isotope}"

    if std_meth_col in standards_categorized.columns:
        std_corrected = standards_categorized[std_meth_col]
    elif std_ref_col in standards_categorized.columns:
        std_corrected = standards_categorized[std_ref_col]
    else:
        raise KeyError(
            f"No corrected standard column found for {isotope}"
        )

    standards_categorized.insert(
        len(standards_categorized.columns) - 1,
        f'Corrected {isotope}',
        std_corrected
    )

    standards_categorized.to_csv(
        os.path.join(res_path, 'Results - standards.csv'),
        index=False
    )

    # -------------------------------------------------
    # Raw replicate export
    # -------------------------------------------------

    raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

    if 'total_error' in raw_unknown_renamed.columns:
        raw_unknown_renamed = raw_unknown_renamed.drop('total_error', axis=1)

    raw_unknown_renamed.to_csv(
        os.path.join(res_path, 'Results - sample replicates.csv'),
        index=False
    )

    total_dD_correction_plot(
        raw_unknown_renamed,
        unknown_renamed,
        folder_path,
        fig_path,
        isotope
    )

    combine_errors(standards_categorized, cfg, log_file_path, isotope)

    print("\nCorrections complete :)")