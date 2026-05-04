# src/pyosotope/EA
import os
import numpy as np

from .base_functions import create_folder, append_to_log
from .utils.config import CorrectionConfig
from .utils.VPDB_correction import VPDB_correction
from .utils.uncertainty_calculation import uncertainty_calculation
from .utils.import_data import load_ea_standards, import_EA_data
from ..corrections.drift_correction import (
    apply_drift_model,
    build_drift_model,
    drift_confirm,
    drop_standards,
    get_isotope,
    get_sorghum,
    plot_drift_diagnostics,
)
from ..corrections.linearity_correction import (
    _get_area_column,
    apply_linearity_model,
    build_linearity_model,
    linearity_confirm,
    load_ea_linearity_metadata,
    plot_linearity_diagnostics,
    prepare_ea_linearity_standards,
)

def ea_process():
    """
    Process Elemental Analyzer (EA-IRMS) isotope data through
    drift correction, VPDB calibration, and uncertainty propagation.

    This function executes the full EA isotope processing workflow.
    It imports raw EA-IRMS data, applies instrumental drift correction,
    performs VPDB normalization using reference standards, calculates
    analytical uncertainties, and writes intermediate and final results
    to disk.

    Processing Steps
    ----------------
    1. Create project output directory structure.
    2. Load EA reference standards.
    3. Import raw EA data file.
    4. Apply time-based drift correction.
    5. Apply VPDB scale calibration.
    6. Compute propagated analytical uncertainties.
    7. Export intermediate and final processed datasets.

    Workflow Details
    ----------------
    Drift Correction
        Corrects for temporal instrument drift using standard measurements
        throughout the analytical run. The correction parameters are saved
        for use in subsequent calibration steps.

    VPDB Calibration
        Converts drift-corrected isotope values to the VPDB scale using
        regression-based normalization against certified reference standards.

    Uncertainty Propagation
        Calculates final analytical uncertainty by combining:
        - Measurement precision
        - Drift model uncertainty
        - Calibration regression uncertainty

    Outputs
    -------
    Drift_Results.csv
        Dataset after drift correction.

    VPDB_Results.csv
        Dataset after VPDB normalization.

    EA_processed_<project_name>.csv
        Final dataset including calibrated isotope values and
        propagated analytical uncertainties.

    Returns
    -------
    None
        Results are written directly to the output directory.
        The function does not return a DataFrame.

    Notes
    -----
    - This function assumes EA reference standards are properly
      defined in the EA standards configuration file.
    - All processing steps append metadata and results to a log file
      for traceability and reproducibility.
    - The processing pipeline is deterministic and does not include
      interactive user input.
    - Errors are propagated analytically using first-order
      variance propagation methods.
    """
    # Setup Output Folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder()
    # standards = load_ea_standards()
    c_std, n_std = load_ea_standards()
    # append_to_log(log_file_path, standards, True)
    append_to_log(log_file_path, c_std, True)
    append_to_log(log_file_path, "")
    append_to_log(log_file_path, n_std, True)

    # Import Data
    df = import_EA_data(loc)

    area_column = _get_area_column(df)

    # Drift Correction
    append_to_log(log_file_path, "\n\n\nDrift Correction:")
    sorghum_n, sorghum_c = get_sorghum(df, log_file_path=log_file_path)
    drop_standards(sorghum_n, df, "N", log_file_path=log_file_path)
    drop_standards(sorghum_c, df, "C", log_file_path=log_file_path)

    n_model = build_drift_model(
        sorghum_n,
        target_column="d 15N/14N",
        index_column="Seconds Since Start",
        log_file_path=log_file_path,
    )
    c_model = build_drift_model(
        sorghum_c,
        target_column="d 13C/12C",
        index_column="Seconds Since Start",
        log_file_path=log_file_path,
    )
    plot_drift_diagnostics(
        n_model,
        fig_path=fig_path,
        figure_name="Drift_Nitrogen.png",
        y_label=r"$\delta^{15}\mathrm{N}$",
        x_label="Seconds Since Start",
    )
    plot_drift_diagnostics(
        c_model,
        fig_path=fig_path,
        figure_name="Drift_Carbon.png",
        y_label=r"$\delta^{13}\mathrm{C}$",
        x_label="Seconds Since Start",
    )

    apply_n, apply_c = drift_confirm(log_file_path=log_file_path)
    cfg = CorrectionConfig(drift_N_applied=apply_n, drift_C_applied=apply_c)
    df["d 15N/14N_corr"] = np.nan
    df["d 15N/14N_se"] = np.nan
    df["d 13C/12C_corr"] = np.nan
    df["d 13C/12C_se"] = np.nan

    if apply_n:
        _, _, el, component = get_isotope("N")
        df = apply_drift_model(
            df,
            n_model,
            target_column="d 15N/14N",
            index_column="Seconds Since Start",
            output_column="d 15N/14N_corr",
            error_column="d 15N/14N_se",
            row_mask=(df["Component"] == component) & df["d 15N/14N"].notna(),
            log_file_path=log_file_path,
            label=el,
            error_nsigma=2.0,
        )
    if apply_c:
        _, _, el, component = get_isotope("C")
        df = apply_drift_model(
            df,
            c_model,
            target_column="d 13C/12C",
            index_column="Seconds Since Start",
            output_column="d 13C/12C_corr",
            error_column="d 13C/12C_se",
            row_mask=(df["Component"] == component) & df["d 13C/12C"].notna(),
            log_file_path=log_file_path,
            label=el,
            error_nsigma=2.0,
        )
    if not (apply_n or apply_c):
        append_to_log(log_file_path, "No drift correction applied.")

    res_name = "Drift_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

#     # Linearity Correction
#     append_to_log(log_file_path, "\n\n\nLinearity Correction:")
#     apply_lin_n, apply_lin_c = linearity_confirm(log_file_path=log_file_path)
#     df["d 15N/14N_lin_corr"] = np.nan
#     df["d 15N/14N_lin_se"] = np.nan
#     df["d 13C/12C_lin_corr"] = np.nan
#     df["d 13C/12C_lin_se"] = np.nan
#
#     if apply_lin_n:
#         try:
#             n_meta = load_ea_linearity_metadata("N", log_file_path=log_file_path)
#             n_stds = prepare_ea_linearity_standards(
#                 df,
#                 n_meta,
#                 identifier_column="Identifier 1",
#                 target_column=cfg.dN_col,
#                 area_column=area_column,
#                 element_type="Nitrogen",
#                 component="N2",
#                 log_file_path=log_file_path,
#             )
#             n_lin_model = build_linearity_model(
#                 n_stds,
#                 target_column=cfg.dN_col,
#                 area_column=area_column,
#                 group_column="linearity_standard_name",
#                 log_file_path=log_file_path,
#             )
#             plot_linearity_diagnostics(
#                 n_lin_model,
#                 fig_path=fig_path,
#                 figure_name="Linearity_Nitrogen.png",
#                 y_label=r"$\delta^{15}\mathrm{N}$",
#                 title="Nitrogen Linearity Standards",
#             )
#             df = apply_linearity_model(
#                 df,
#                 n_lin_model,
#                 input_column=cfg.dN_col,
#                 area_column=area_column,
#                 output_column="d 15N/14N_lin_corr",
#                 error_column="d 15N/14N_lin_se",
#                 row_mask=(df["Component"] == "N2") & df[cfg.dN_col].notna() & df[area_column].notna(),
#                 log_file_path=log_file_path,
#                 label="Nitrogen",
#             )
#             cfg.linearity_N_applied = True
#         except Exception as exc:
#             cfg.linearity_N_applied = False
#             warning = f"Skipping Nitrogen linearity correction: {exc}"
#             print(warning)
#             append_to_log(log_file_path, warning)
#     if apply_lin_c:
#         try:
#             c_meta = load_ea_linearity_metadata("C", log_file_path=log_file_path)
#             c_stds = prepare_ea_linearity_standards(
#                 df,
#                 c_meta,
#                 identifier_column="Identifier 1",
#                 target_column=cfg.dC_col,
#                 area_column=area_column,
#                 element_type="Carbon",
#                 component="CO2",
#                 log_file_path=log_file_path,
#             )
#             c_lin_model = build_linearity_model(
#                 c_stds,
#                 target_column=cfg.dC_col,
#                 area_column=area_column,
#                 group_column="linearity_standard_name",
#                 log_file_path=log_file_path,
#             )
#             plot_linearity_diagnostics(
#                 c_lin_model,
#                 fig_path=fig_path,
#                 figure_name="Linearity_Carbon.png",
#                 y_label=r"$\delta^{13}\mathrm{C}$",
#                 title="Carbon Linearity Standards",
#             )
#             df = apply_linearity_model(
#                 df,
#                 c_lin_model,
#                 input_column=cfg.dC_col,
#                 area_column=area_column,
#                 output_column="d 13C/12C_lin_corr",
#                 error_column="d 13C/12C_lin_se",
#                 row_mask=(df["Component"] == "CO2") & df[cfg.dC_col].notna() & df[area_column].notna(),
#                 log_file_path=log_file_path,
#                 label="Carbon",
#             )
#             cfg.linearity_C_applied = True
#         except Exception as exc:
#             cfg.linearity_C_applied = False
#             warning = f"Skipping Carbon linearity correction: {exc}"
#             print(warning)
#             append_to_log(log_file_path, warning)
#     if not (apply_lin_n or apply_lin_c):
#         append_to_log(log_file_path, "No linearity correction applied.")

    # res_name = "Linearity_Results.csv"
    # res = os.path.join(results_path, res_name)
    # df.to_csv(res)

    # VPDB Calibration
    # df = VPDB_correction(df, standards, cfg, log_file_path, fig_path)
    df = VPDB_correction(df, c_std, n_std, cfg, log_file_path, fig_path)
    res_name = "VPDB_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # Uncertainty Calculation
    final_df = uncertainty_calculation(df,cfg,log_file_path)
    project_name = str(os.path.basename(loc))
    res = os.path.join(results_path, f"EA_processed_{project_name}")
    final_df.to_csv(res)
