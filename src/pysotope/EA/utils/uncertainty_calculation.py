import numpy as np
import pandas as pd

from . .base_functions import append_to_log
from .VPDB_correction import get_isotope as VPDB_get_isotope
from .ea_drift_correction import get_isotope as drift_get_isotope

def uncertainty_calculation(df, cfg, log_file):
    """
    Calculate total uncertainty using combine_errors(), append the total uncertainty for
    each sample into the DataFrame and output the final processed DataFrame as a .csv file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing drift corrected (if applied) and VPDB calibrated EA-IRMS data.
    cfg : CorrectionConfig
        A configuration object holding boolean values for which corrections were applied.
    log_file : str
         Path to the log file where processing details will be logged.

     Returns
    -------
    pandas.DataFrame
        DataFrame containing drift corrected (if applied) and VPDB calibrated EA-IRMS data
        with calculated combined uncertainty.
    """

    N_df = df[df['Element Type'] == 'Nitrogen'].copy()
    C_df = df[df['Element Type'] == 'Carbon'].copy()
    combine_errors(N_df,cfg,log_file,'N')
    combine_errors(C_df, cfg, log_file, 'C')
    final_df = pd.concat([N_df, C_df], ignore_index=False)
    final_df.sort_values(by="Seconds Since Start", inplace=True)

    return final_df

def combine_errors(df, cfg, log_file, tag):
    """
    Combine standard error for drift correction (if applied) and VPDB calibration using sum
    of squares, for both Nitrogen and Carbon.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing drift corrected (if applied) and VPDB calibrated EA-IRMS data.
    cfg : CorrectionConfig
        A configuration object holding boolean values for which corrections were applied.
    log_file : str
         Path to the log file where processing details will be logged.
    tag : str
        Tag for the isotope being processed for uncertainty calculation.
        - 'N' → Nitrogen
        - 'C' → Carbon

    Notes
    -----
    - Logs the statistical summary of uncertainty calculation and corrections to the log file.
    """

    vpdb_input_col, vpdb_col, actual_col, iso_label, el = VPDB_get_isotope(tag, cfg)
    drift_input_col, iso_label, el, comp = drift_get_isotope(tag)
    error_terms = [df[f"{vpdb_col}_se"]]  # Always include VSMOW_error

    if vpdb_input_col.endswith("_corr"):
        error_terms.append(df[f"{drift_input_col}_se"])
    # Combine using root sum of squares
    squared = [err ** 2 for err in error_terms]
    df['sum_squared_error'] = np.sqrt(sum(squared))

    summary_df = pd.DataFrame({f'mean_{vpdb_col}': [df[vpdb_col].mean()]})

    if vpdb_input_col.endswith("_corr"):
        summary_df[f'mean_{vpdb_input_col}'] = [df[vpdb_input_col].mean()]

    summary_df["total_error"] = [df['sum_squared_error'].mean()]

    append_to_log(log_file, f"\nSummary of {el} samples")
    append_to_log(log_file, summary_df)


# def mean_values_with_uncertainty(data, cfg, tag, sample_name_header="Identifier 1"):
#     grouped = data.groupby([sample_name_header])
#     vpdb_input_col, vpdb_col, actual_col, iso_label, el = VPDB_get_isotope(tag, cfg)
#     drift_input_col, iso_label, el, comp = drift_get_isotope(tag)
#
#
#     # 1) build the basic aggregation dict
#     agg_dict = {vpdb_col: ["mean","count"], f"{vpdb_col}_se" : ["mean"]}
#
#     # 2) drift (if applied)
#     if vpdb_input_col.endswith("_corr"):
#         agg_dict.update({
#             vpdb_input_col: ["mean"],
#             f"{drift_input_col}_se": ["mean"]
#         })
#
#     # 5) perform the aggregation
#     stats = grouped.agg(agg_dict).reset_index()
#     stats.columns = ['_'.join(col).strip('_') for col in stats.columns.values]
#
#     # 7) compute replicate SEM
#     if vpdb_input_col.endswith("_corr"):
#         stats[f"replicate_{drift_input_col}_sem"] = stats[f"{vpdb_input_col}_mean"] / np.sqrt(stats[f"{vpdb_col}_count"])
#         stats[f"replicate_{drift_input_col}_error_sem"] = stats[f"{drift_input_col}_se_mean"] / np.sqrt(stats[f"{vpdb_col}_count"])
#
#     stats[f"replicate_{vpdb_col}_sem"] = stats[f"{vpdb_col}_mean"] / np.sqrt(stats[f"{vpdb_col}_count"])
#     stats[f"replicate_{vpdb_col}_error_sem"] = stats[f"{vpdb_col}_se_mean"] / np.sqrt(stats[f"{vpdb_col}_count"])
#
#     # 8) build error list for total uncertainty
#     error_cols = []
#     if vpdb_input_col.endswith("_corr"):
#         error_cols.append(f"replicate_{drift_input_col}_error_sem")
#
#     error_cols.append(f"replicate_{vpdb_col}_error_sem")
#
#     # 9) compute total_uncertainty = sqrt(sum(errors²))
#     stats[f"total_uncertainty_{tag}"] = np.sqrt(np.sum([stats[col] ** 2 for col in error_cols], axis=0))
#
#     return stats