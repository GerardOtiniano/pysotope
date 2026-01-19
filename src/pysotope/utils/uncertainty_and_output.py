import numpy as np
import os
import pandas as pd
from .figures import total_dD_correction_plot
from .config import CorrectionConfig
from .base_functions import append_to_log

# def output_results(raw_unknown, unknown, sd, unknown_pame, folder_path, fig_path, res_path, isotope, pame, log_file_path, cfg):
#     # Define column name mappings for output CSV files
#     column_name_mapping = {
#         'Identifier 1': 'Sample Name',
#         'chain': 'Chain Length',
#         'area_mean': 'Mean Area',
#         'dD_mean': f'Average raw {isotope}',
#         'dD_std': f'{isotope} Std Dev',
#         'dD_count': 'Number of Replicates',
#         'drift_corrected_dD_mean': f'Drift Corrected {isotope}',
#         'drift_error_mean': 'Drift Error',
#         'linearity_corrected_dD_mean': f'Linearity Corrected {isotope}',
#         'linearity_error_mean': 'Linearity Error',
#         'VSMOW_dD_mean': f'VSMOW Corrected {isotope}',
#         'VSMOW_error_mean': 'VSMOW Error',
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
#     if isotope=='dD':
#         unknown_renamed.insert(len(unknown_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_renamed['Methanol Corrected '+str(isotope)])
#         if pame:
#             unknown_pame_renamed.insert(len(unknown_pame_renamed.columns) - 1, 'Corrected Methanol '+str(isotope), unknown_pame_renamed['PAME Methanol Calculated '+str(isotope)])
#     else:
#         unknown_renamed.reset_index(drop=True, inplace=True)
#         unknown_renamed.insert(len(unknown_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_renamed['VSMOW Corrected dC'])
#         if pame:
#             unknown_pame_renamed.reset_index(drop=True, inplace=True)
#             unknown_pame_renamed.insert(len(unknown_pame_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_pame_renamed['VSMOW Corrected dC'])
#     # Save unknown dataframe to CSV
#     unknown_renamed.reset_index(drop=True, inplace=True)
#     unknown_renamed.to_csv(os.path.join(res_path,'Results - sample mean.csv'), index=False)
#     if pame:
#         unknown_pame_renamed.reset_index(drop=True, inplace=True)
#         unknown_pame_renamed.to_csv(os.path.join(res_path,'PAME Results - sample mean.csv'), index=False)
#     # Select columns for standards dataframe
#     columns_to_select = [
#         "Date", "Time", "Identifier 1", "chain", "Rt", "area", "Area 2", "Area 3",
#         "Ampl  2", "Ampl  3", "BGD 2", 'BGD 3', "time_rel", "dD", "drift_corrected_dD",
#         "drift_error", "linearity_corrected_dD", "linearity_error", "VSMOW_dD", "VSMOW_error",
#         "total_uncertainty"]

#     # Select only existing columns from sd dataframe
#     existing_columns = [col for col in columns_to_select if col in sd.columns]
#     standards_selected = sd[existing_columns].copy()

#     # Filter and categorize standards
#     lin_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C20|C28")].copy()
#     lin_std_temp['Standard Type'] = "Linearity"
    
#     drift_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C218|C24")].copy()
#     drift_std_temp['Standard Type'] = "Drift"
    
#     standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

#     # Rename columns in standards dataframe
#     column_rename_map = {
#         "Identifier 1": "Standard ID", "chain": "Component", "Rt": "Retention time",
#         "area": "Peak area", "Area 2": "Peak area 2", "Area 3": "Peak area 3",
#         "Ampl  2": "Amplitude 2", "Ampl  3": "Amplitude 3", "BGD 2": "Background 2",
#         "BGD 3": "Background 3", "time_rel": "Time relative", "dD": f"Raw {isotope}",
#         "drift_corrected_dD": f"Drift corrected {isotope}", "drift_error": f"Drift {isotope} error",
#         "linearity_corrected_dD": f"Linearity {isotope}", "linearity_error": f"Linearity {isotope} error",
#         "VSMOW_dD": f"VSMOW corrected {isotope}", "vsmow_error": "VSMOW error",
#         "methanol_dD": f"Methanol corrected {isotope}", "methanol_error": f"Methanol error {isotope}",
#         "total_uncertainty": "Total uncertainty"}
#     standards_categorized = standards_categorized.rename(columns=column_rename_map)
#     standards_categorized.insert(len(standards_categorized.columns) - 1, 'Corrected '+str(isotope), standards_categorized[f"VSMOW corrected {isotope}"])
#     # Save standards dataframe to CSV
#     standards_categorized.to_csv(os.path.join(res_path, 'Results - standards.csv'), index=False)

#     # Rename columns in raw_unknown dataframe
#     raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

#     # Drop 'total_error' column if exists
#     if 'total_error' in raw_unknown.columns:
#         raw_unknown_renamed = raw_unknown_renamed.drop('total_error', axis=1)

#     # Save raw_unknown dataframe to CSV
#     raw_unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample replicates.csv'), index=False)
#     total_dD_correction_plot(raw_unknown_renamed, unknown_renamed, folder_path, fig_path, isotope)
#     combine_errors(standards_categorized, cfg, log_file_path, isotope)
#     print("\nCorrections complete :)")

def output_results(raw_unknown, unknown, sd, unknown_pame, folder_path, fig_path, res_path, isotope, pame, log_file_path, cfg):
    # Reference scale label for OUTPUT headers only
    ref_label = "VSMOW" if isotope == "dD" else "PDB"

    # Define column name mappings for output CSV files
    column_name_mapping = {
        'Identifier 1': 'Sample Name',
        'chain': 'Chain Length',
        'area_mean': 'Mean Area',

        # raw replicate stats (still keyed as dD_* in your current stats output)
        'dD_mean': f'Average raw {isotope}',
        'dD_std': f'{isotope} Std Dev',
        'dD_count': 'Number of Replicates',

        'drift_corrected_dD_mean': f'Drift Corrected {isotope}',
        'drift_error_mean': 'Drift Error',
        'linearity_corrected_dD_mean': f'Linearity Corrected {isotope}',
        'linearity_error_mean': 'Linearity Error',

        # *** OUTPUT LABEL changes here ***
        'VSMOW_dD_mean': f'{ref_label} Corrected {isotope}',
        'VSMOW_error_mean': f'{ref_label} Error',

        'methanol_dD_mean': f'Methanol Corrected {isotope}',
        'methanol_error_mean': f'Methanol Error',

        'PAME_methanol_dD_mean': f'PAME Methanol Calculated {isotope}',
        'PAME_methanol_dD_std': f'PAME Methanol Calculated Error {isotope}',

        'replicate_dD_sem': f'Mean replicate std {isotope}',
        'total_uncertainty': 'Total Uncertainty'
    }

    # Rename columns in unknown dataframe
    unknown_renamed = unknown.rename(columns=column_name_mapping)
    if pame:
        unknown_pame_renamed = unknown_pame.rename(columns=column_name_mapping)

    if isotope == 'dD':
        unknown_renamed.insert(
            len(unknown_renamed.columns) - 1,
            'Corrected ' + str(isotope),
            unknown_renamed['Methanol Corrected ' + str(isotope)]
        )
        if pame:
            unknown_pame_renamed.insert(
                len(unknown_pame_renamed.columns) - 1,
                'Corrected Methanol ' + str(isotope),
                unknown_pame_renamed['PAME Methanol Calculated ' + str(isotope)]
            )
    else:
        unknown_renamed.reset_index(drop=True, inplace=True)

        # *** was: 'VSMOW Corrected dC' ***
        unknown_renamed.insert(
            len(unknown_renamed.columns) - 1,
            'Corrected ' + str(isotope),
            unknown_renamed[f'{ref_label} Corrected {isotope}']
        )
        if pame:
            unknown_pame_renamed.reset_index(drop=True, inplace=True)
            unknown_pame_renamed.insert(
                len(unknown_pame_renamed.columns) - 1,
                'Corrected ' + str(isotope),
                unknown_pame_renamed[f'{ref_label} Corrected {isotope}']
            )

    # Save unknown dataframe to CSV
    unknown_renamed.reset_index(drop=True, inplace=True)
    unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample mean.csv'), index=False)
    if pame:
        unknown_pame_renamed.reset_index(drop=True, inplace=True)
        unknown_pame_renamed.to_csv(os.path.join(res_path, 'PAME Results - sample mean.csv'), index=False)

    # Select columns for standards dataframe (internal keys unchanged)
    columns_to_select = [
        "Date", "Time", "Identifier 1", "chain", "Rt", "area", "Area 2", "Area 3",
        "Ampl  2", "Ampl  3", "BGD 2", 'BGD 3', "time_rel", "dD", "drift_corrected_dD",
        "drift_error", "linearity_corrected_dD", "linearity_error", "VSMOW_dD", "VSMOW_error",
        "total_uncertainty"
    ]

    existing_columns = [col for col in columns_to_select if col in sd.columns]
    standards_selected = sd[existing_columns].copy()

    # Filter and categorize standards
    lin_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C20|C28")].copy()
    lin_std_temp['Standard Type'] = "Linearity"

    drift_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C218|C24")].copy()
    drift_std_temp['Standard Type'] = "Drift"

    standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

    # Rename columns in standards dataframe (OUTPUT LABEL changes here too)
    column_rename_map = {
        "Identifier 1": "Standard ID",
        "chain": "Component",
        "Rt": "Retention time",
        "area": "Peak area",
        "Area 2": "Peak area 2",
        "Area 3": "Peak area 3",
        "Ampl  2": "Amplitude 2",
        "Ampl  3": "Amplitude 3",
        "BGD 2": "Background 2",
        "BGD 3": "Background 3",
        "time_rel": "Time relative",

        "dD": f"Raw {isotope}",
        "drift_corrected_dD": f"Drift corrected {isotope}",
        "drift_error": f"Drift {isotope} error",
        "linearity_corrected_dD": f"Linearity {isotope}",
        "linearity_error": f"Linearity {isotope} error",

        # *** OUTPUT LABEL changes here ***
        "VSMOW_dD": f"{ref_label} corrected {isotope}",
        "VSMOW_error": f"{ref_label} error",

        "methanol_dD": f"Methanol corrected {isotope}",
        "methanol_error": f"Methanol error {isotope}",
        "total_uncertainty": "Total uncertainty"
    }

    standards_categorized = standards_categorized.rename(columns=column_rename_map)

    # Insert corrected column using ref_label
    standards_categorized.insert(
        len(standards_categorized.columns) - 1,
        'Corrected ' + str(isotope),
        standards_categorized[f"{ref_label} corrected {isotope}"]
    )

    standards_categorized.to_csv(os.path.join(res_path, 'Results - standards.csv'), index=False)

    # Rename columns in raw_unknown dataframe
    raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

    if 'total_error' in raw_unknown.columns:
        raw_unknown_renamed = raw_unknown_renamed.drop('total_error', axis=1)

    raw_unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample replicates.csv'), index=False)

    total_dD_correction_plot(raw_unknown_renamed, unknown_renamed, folder_path, fig_path, isotope)
    combine_errors(standards_categorized, cfg, log_file_path, isotope)
    print("\nCorrections complete :)")


    
    
    

def combine_errors(df, config, log_file_path, isotope):
    """
    Combine error terms using sum of squares based on the applied corrections,
    always including VSMOW_error.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing error columns.
        config (CorrectionConfig): Configuration object specifying which corrections were applied.

    Returns:
        pd.Series: A Series with the combined error for each row.
    """
    if isotope == "dD":
        error_terms = [df['VSMOW error']]  # Always include VSMOW_error
    else:
        error_terms = [df["PDB error"]]
    if config.linearity_applied:
        error_terms.append(df[f'Linearity {isotope} error'])

    if config.drift_applied:
        error_terms.append(df[f'Drift {isotope} error'])
    # Combine using root sum of squares
    squared = [err**2 for err in error_terms]
    df['Total Error'] = np.sqrt(sum(squared))
    summary = df.groupby('Component').agg(
        mean_corrected_dD=(f'Corrected {isotope}', 'mean'),
        sum_squared_error=('Total Error', lambda x: np.mean(x))).reset_index()
    append_to_log(log_file_path, "Summary of standards")
    append_to_log(log_file_path, summary)
    

# def mean_values_with_uncertainty(data, cfg, iso, sample_name_header="Identifier 1", chain_header="chain"):
#     """
#     Group by sample & chain and compute means + uncertainties,
#     only aggregating drift/linearity columns if those corrections
#     were applied (per cfg).
#     """
#     grouped = data.groupby([sample_name_header, chain_header])

#     # 1) build the basic aggregation dict
#     agg_dict = {"area": ["mean"]}

#     # 2) drift (if applied)
#     if cfg.drift_applied:
#         agg_dict.update({
#             "drift_corrected_dD": "mean",
#             "drift_error":        "mean"
#         })

#     # 3) linearity (if applied)
#     if cfg.linearity_applied:
#         agg_dict.update({
#             "linearity_corrected_dD": "mean",
#             "linearity_error":        "mean"
#         })

#     # 4) always include VSMOW for whichever iso
#     if iso == "dC":
#         agg_dict.update({
#             "VSMOW_dC":     ["mean", "std"],
#             "VSMOW_error":  "mean",
#             "dC":           ["mean", "std", "count"],
#         })
#     else:  # dD
#         agg_dict.update({
#             "VSMOW_dD":    ["mean", "std"],
#             "VSMOW_error": "mean",
#             "dD":          ["mean", "std", "count"],
#         })
#         # methanol / PAME if present
#         if "methanol_dD" in data.columns:
#             agg_dict.update({
#                 "methanol_dD":   ["mean", "std"],
#                 "methanol_error": "mean"
#             })
#         if "PAME_methanol_dD" in data.columns:
#             agg_dict.update({
#                 "PAME_methanol_dD": ["mean", "std"]
#             })

#     # 5) perform the aggregation
#     stats = grouped.agg(agg_dict).reset_index()

#     # 6) flatten MultiIndex
#     stats.columns = [
#         "_".join(filter(None, col)).strip("_")
#         for col in stats.columns.values
#     ]

#     # 7) compute replicate SEM
#     if iso == "dD":
#         # choose whichever std column exists
#         if "methanol_dD_std" in stats.columns and "dD_count" in stats.columns:
#             stats["replicate_dD_sem"] = stats["methanol_dD_std"] / np.sqrt(stats["dD_count"])
#         elif "PAME_methanol_dD_std" in stats.columns and "dD_count" in stats.columns:
#             stats["replicate_dD_sem"] = stats["PAME_methanol_dD_std"] / np.sqrt(stats["dD_count"])
#         else:
#             stats["replicate_dD_sem"] = np.nan
#     else:  # dC
#         if "dC_count" in stats.columns:
#             stats["replicate_dC_sem"] = stats["VSMOW_error_mean"] / np.sqrt(stats["dC_count"])
#         else:
#             stats["replicate_dC_sem"] = np.nan

#     # 8) build error list for total uncertainty
#     error_cols = []
#     if cfg.drift_applied:
#         error_cols.append("drift_error_mean")
#     if cfg.linearity_applied:
#         error_cols.append("linearity_error_mean")

#     error_cols.append("VSMOW_error_mean")

#     if iso == "dD":
#         if "methanol_error_mean" in stats.columns:
#             error_cols.append("methanol_error_mean")
#         # replicate SEM already in replicate_dD_sem
#         error_cols.append("replicate_dD_sem")
#     else:
#         error_cols.append("replicate_dC_sem")

#     # 9) compute total_uncertainty = sqrt(sum(errors²))
#     stats["total_uncertainty"] = np.sqrt(
#         np.sum([stats[col]**2 for col in error_cols], axis=0)
#     )

#     return stats


def mean_values_with_uncertainty(data, cfg, iso, sample_name_header="Identifier 1", chain_header="chain"):
    """
    Group by sample & chain and compute means + uncertainties.

    Compatible with both dD and dC even if upstream code still produces
    dD-suffixed correction columns during dC processing (e.g., drift_corrected_dD).

    Outputs (same *type* of info as before):
      - area_mean
      - raw iso mean/std/count (e.g., dC_mean, dC_std, dC_count)
      - drift/linearity corrected means + mean errors (if applied)
      - VSMOW mean/std + VSMOW_error_mean (if available)
      - methanol / PAME for dD if present
      - replicate_{iso}_sem
      - total_uncertainty
    """

    def pick_first_existing(cols, candidates):
        for c in candidates:
            if c in cols:
                return c
        return None

    # -----------------------------
    # 0) decide which columns to use
    # -----------------------------
    cols = set(data.columns)

    iso_col = iso  # must exist: "dD" or "dC"
    if iso_col not in cols:
        raise KeyError(f"Expected raw isotope column '{iso_col}' not found. Columns: {list(data.columns)}")

    # Correction column candidates: prefer iso-specific, fall back to legacy dD names
    drift_corr_col = pick_first_existing(cols, [f"drift_corrected_{iso}", "drift_corrected_dD", "drift_corrected_dC"])
    lin_corr_col   = pick_first_existing(cols, [f"linearity_corrected_{iso}", "linearity_corrected_dD", "linearity_corrected_dC"])

    drift_err_col  = pick_first_existing(cols, ["drift_error"])
    lin_err_col    = pick_first_existing(cols, ["linearity_error"])

    # VSMOW column: prefer iso-specific, fall back to whatever exists (your dC run currently has VSMOW_dD)
    vsmow_col = pick_first_existing(cols, [f"VSMOW_{iso}", "VSMOW_dD", "VSMOW_dC"])
    vsmow_err_col = pick_first_existing(cols, ["VSMOW_error"])

    # Methanol/PAME are dD-specific in your codebase
    meth_col = pick_first_existing(cols, ["methanol_dD"])
    meth_err_col = pick_first_existing(cols, ["methanol_error"])
    pame_col = pick_first_existing(cols, ["PAME_methanol_dD"])

    # -----------------------------
    # 1) group and aggregate
    # -----------------------------
    grouped = data.groupby([sample_name_header, chain_header], dropna=False)

    agg_dict = {}

    # area
    if "area" in cols:
        agg_dict["area"] = ["mean"]

    # drift (if applied + columns exist)
    if getattr(cfg, "drift_applied", False):
        if drift_corr_col is not None:
            agg_dict[drift_corr_col] = ["mean"]
        if drift_err_col is not None:
            agg_dict[drift_err_col] = ["mean"]

    # linearity (if applied + columns exist)
    if getattr(cfg, "linearity_applied", False):
        if lin_corr_col is not None:
            agg_dict[lin_corr_col] = ["mean"]
        if lin_err_col is not None:
            agg_dict[lin_err_col] = ["mean"]

    # VSMOW (if present)
    if vsmow_col is not None:
        agg_dict[vsmow_col] = ["mean", "std"]
    if vsmow_err_col is not None:
        agg_dict[vsmow_err_col] = ["mean"]

    # raw isotope always
    agg_dict[iso_col] = ["mean", "std", "count"]

    # methanol / PAME if present (same behavior as before: only if columns exist)
    if meth_col is not None:
        agg_dict[meth_col] = ["mean", "std"]
    if meth_err_col is not None:
        agg_dict[meth_err_col] = ["mean"]
    if pame_col is not None:
        agg_dict[pame_col] = ["mean", "std"]

    # Run aggregation
    stats = grouped.agg(agg_dict).reset_index()

    # Flatten MultiIndex columns
    stats.columns = ["_".join(filter(None, col)).strip("_") for col in stats.columns.values]

    # -----------------------------
    # 2) compute replicate SEM
    # -----------------------------
    # Decide which *std* column represents the measurement scatter best.
    # For dD: prefer methanol std (if present), else PAME std, else VSMOW std, else raw dD std.
    # For dC: prefer VSMOW std (if present), else raw dC std.
    n_col = f"{iso_col}_count"

    if iso == "dD":
        std_candidates = []
        if meth_col is not None:
            std_candidates.append(f"{meth_col}_std")
        if pame_col is not None:
            std_candidates.append(f"{pame_col}_std")
        if vsmow_col is not None:
            std_candidates.append(f"{vsmow_col}_std")
        std_candidates.append(f"{iso_col}_std")  # fallback

        sem_std_col = pick_first_existing(set(stats.columns), std_candidates)
        sem_name = "replicate_dD_sem"
    else:
        std_candidates = []
        if vsmow_col is not None:
            std_candidates.append(f"{vsmow_col}_std")
        std_candidates.append(f"{iso_col}_std")
        sem_std_col = pick_first_existing(set(stats.columns), std_candidates)
        sem_name = "replicate_dC_sem"

    if sem_std_col is not None and n_col in stats.columns:
        stats[sem_name] = stats[sem_std_col] / np.sqrt(stats[n_col].replace(0, np.nan))
    else:
        stats[sem_name] = np.nan

    # -----------------------------
    # 3) build error list for total uncertainty
    # -----------------------------
    error_cols = []

    # drift/linearity errors are means after flattening: "<col>_mean"
    if getattr(cfg, "drift_applied", False) and drift_err_col is not None:
        c = f"{drift_err_col}_mean"
        if c in stats.columns:
            error_cols.append(c)

    if getattr(cfg, "linearity_applied", False) and lin_err_col is not None:
        c = f"{lin_err_col}_mean"
        if c in stats.columns:
            error_cols.append(c)

    # VSMOW error (if present)
    if vsmow_err_col is not None:
        c = f"{vsmow_err_col}_mean"
        if c in stats.columns:
            error_cols.append(c)

    # methanol error only relevant for dD if present
    if iso == "dD" and meth_err_col is not None:
        c = f"{meth_err_col}_mean"
        if c in stats.columns:
            error_cols.append(c)

    # replicate SEM always included
    if sem_name in stats.columns:
        error_cols.append(sem_name)

    # -----------------------------
    # 4) compute total uncertainty
    # -----------------------------
    if len(error_cols) == 0:
        stats["total_uncertainty"] = np.nan
    else:
        # robust: treat missing values as nan; sum squares across columns row-wise
        sq = np.zeros(len(stats), dtype=float)
        sq[:] = 0.0
        for c in error_cols:
            sq += (stats[c].astype(float) ** 2).to_numpy()
        stats["total_uncertainty"] = np.sqrt(sq)

    return stats