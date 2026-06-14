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
    get_ea_drift_standards,
    get_isotope,
    plot_drift_diagnostics,
    review_ea_reference_standards,
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


def _format_standard_names(std_df):
    names = std_df["Identifier 1"].dropna().astype(str).str.strip()
    names = names[names != ""].drop_duplicates()
    return ", ".join(names)


def ea_process():
    """
    Process Elemental Analyzer (EA-IRMS) isotope data through
    drift correction, reference standard calibration, and uncertainty propagation.

    This function executes the full EA isotope processing workflow.
    It imports raw EA-IRMS data, applies instrumental drift correction,
    performs reference standard normalization, calculates
    analytical uncertainties, and writes intermediate and final results
    to disk.

    Processing Steps
    ----------------
    1. Create project output directory structure.
    2. Load EA reference standards.
    3. Import raw EA data file.
    4. Apply time-based drift correction.
    5. Apply reference standard calibration.
    6. Compute propagated analytical uncertainties.
    7. Export intermediate and final processed datasets.

    Workflow Details
    ----------------
    Drift Correction
        Corrects for temporal instrument drift using standard measurements
        throughout the analytical run. The correction parameters are saved
        for use in subsequent calibration steps.

    Reference Standard Calibration
        Converts drift-corrected isotope values to the appropriate reference
        scale using regression-based normalization against certified standards.

    Uncertainty Propagation
        Calculates final analytical uncertainty by combining:
        - Measurement precision
        - Drift model uncertainty
        - Calibration regression uncertainty

    Outputs
    -------
    Drift_Results.csv
        Dataset after drift correction.

    Reference Standard Results.csv
        Dataset after reference standard normalization.

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
    review_ea_reference_standards(df, c_std, n_std, log_file_path=log_file_path)

    # Drift Correction
    append_to_log(log_file_path, "\n\n\nDrift Correction:")

    def prepare_drift_diagnostic(tag, figure_name):
        col, y_label, el, component = get_isotope(tag)
        drift_standards = get_ea_drift_standards(df, tag, log_file_path=log_file_path)
        drift_model = build_drift_model(
            drift_standards,
            target_column=col,
            index_column="Seconds Since Start",
            group_column="Identifier 1",
            log_file_path=log_file_path,
        )
        plot_drift_diagnostics(
            drift_model,
            fig_path=fig_path,
            figure_name=figure_name,
            y_label=y_label,
            x_label="Seconds Since Start",
        )
        return {
            "tag": tag,
            "column": col,
            "element": el,
            "component": component,
            "standards": drift_standards,
            "model": drift_model,
        }

    drift_n = prepare_drift_diagnostic("N", "Drift_Nitrogen.png")
    drift_c = prepare_drift_diagnostic("C", "Drift_Carbon.png")

    apply_n, apply_c = drift_confirm(log_file_path=log_file_path)
    cfg = CorrectionConfig(drift_N_applied=apply_n, drift_C_applied=apply_c)
    df["d 15N/14N_corr"] = np.nan
    df["d 15N/14N_se"] = np.nan
    df["d 13C/12C_corr"] = np.nan
    df["d 13C/12C_se"] = np.nan

    def apply_selected_drift(drift_info, output_column, error_column):
        col = drift_info["column"]
        el = drift_info["element"]
        component = drift_info["component"]
        append_to_log(
            log_file_path,
            f"{el} drift correction reference standard(s) used: {_format_standard_names(drift_info['standards'])}.",
        )
        return apply_drift_model(
            df,
            drift_info["model"],
            target_column=col,
            index_column="Seconds Since Start",
            output_column=output_column,
            error_column=error_column,
            row_mask=(df["Component"] == component) & df[col].notna(),
            log_file_path=log_file_path,
            label=el,
            error_nsigma=2.0,
        )

    if apply_n:
        df = apply_selected_drift(drift_n, "d 15N/14N_corr", "d 15N/14N_se")
    else:
        append_to_log(log_file_path, "Nitrogen drift model was previewed but not applied.")
    if apply_c:
        df = apply_selected_drift(drift_c, "d 13C/12C_corr", "d 13C/12C_se")
    else:
        append_to_log(log_file_path, "Carbon drift model was previewed but not applied.")
    if not (apply_n or apply_c):
        append_to_log(log_file_path, "No drift correction applied.")

    res_name = "Drift_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # Linearity Correction
    append_to_log(log_file_path, "\n\n\nLinearity Correction:")
    apply_lin_n, apply_lin_c = linearity_confirm(log_file_path=log_file_path)
    df["d 15N/14N_lin_corr"] = np.nan
    df["d 15N/14N_lin_se"] = np.nan
    df["d 13C/12C_lin_corr"] = np.nan
    df["d 13C/12C_lin_se"] = np.nan

    def confirm_linearity_application(element):
        while True:
            response = input(f"Apply {element} linearity correction? (Y/N)\n").strip().lower()
            if response in {"yes", "y", "true", "t"}:
                append_to_log(log_file_path, f"User accepted EA linearity correction for {element}.")
                return True
            if response in {"no", "n", "false", "f", ""}:
                append_to_log(log_file_path, f"User rejected EA linearity correction for {element}.")
                return False
            print("Wrong input. Please enter Y or N.")

    def process_selected_linearity(tag, output_column, error_column, figure_name):
        input_col, input_desc = cfg.dN_col if tag == "N" else cfg.dC_col
        _, y_label, element, component = get_isotope(tag)
        append_to_log(
            log_file_path,
            f"{element} linearity correction will use {input_desc.lower()} isotope values from column '{input_col}'.",
        )
        metadata = load_ea_linearity_metadata(tag, log_file_path=log_file_path)
        standards = prepare_ea_linearity_standards(
            df,
            metadata,
            identifier_column="Identifier 1",
            target_column=input_col,
            area_column=area_column,
            element_type=element,
            component=component,
            log_file_path=log_file_path,
        )
        model = build_linearity_model(
            standards,
            target_column=input_col,
            area_column=area_column,
            group_column="linearity_standard_name",
            log_file_path=log_file_path,
        )
        plot_linearity_diagnostics(
            model,
            fig_path=fig_path,
            figure_name=figure_name,
            y_label=y_label,
            title=f"{element} Linearity Standards",
        )
        print(
            f"\n{element} linearity model: {model['best_model']}",
            f"R² = {model['r_squared']:.3f} | SSE = {model['sse']:.3f}",
            sep="\n",
            end="\n\n",
        )
        if not confirm_linearity_application(element):
            return False
        nonlocal_df = apply_linearity_model(
            df,
            model,
            input_column=input_col,
            area_column=area_column,
            output_column=output_column,
            error_column=error_column,
            row_mask=(df["Component"] == component) & df[input_col].notna() & df[area_column].notna(),
            log_file_path=log_file_path,
            label=element,
        )
        return nonlocal_df

    if apply_lin_n:
        try:
            lin_result = process_selected_linearity(
                "N",
                "d 15N/14N_lin_corr",
                "d 15N/14N_lin_se",
                "Linearity_Nitrogen.png",
            )
            if lin_result is not False:
                df = lin_result
                cfg.linearity_N_applied = True
        except Exception as exc:
            cfg.linearity_N_applied = False
            warning = f"Skipping Nitrogen linearity correction: {exc}"
            print(warning)
            append_to_log(log_file_path, warning)
    else:
        append_to_log(log_file_path, "Nitrogen linearity correction was not selected.")

    if apply_lin_c:
        try:
            lin_result = process_selected_linearity(
                "C",
                "d 13C/12C_lin_corr",
                "d 13C/12C_lin_se",
                "Linearity_Carbon.png",
            )
            if lin_result is not False:
                df = lin_result
                cfg.linearity_C_applied = True
        except Exception as exc:
            cfg.linearity_C_applied = False
            warning = f"Skipping Carbon linearity correction: {exc}"
            print(warning)
            append_to_log(log_file_path, warning)
    else:
        append_to_log(log_file_path, "Carbon linearity correction was not selected.")

    if not (cfg.linearity_N_applied or cfg.linearity_C_applied):
        append_to_log(log_file_path, "No linearity correction applied.")

    res_name = "Linearity_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # Reference standard calibration
    # df = VPDB_correction(df, standards, cfg, log_file_path, fig_path)
    df = VPDB_correction(df, c_std, n_std, cfg, log_file_path, fig_path)
    res_name = "Reference Standard Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # Uncertainty Calculation
    final_df = uncertainty_calculation(df,cfg,log_file_path)
    project_name = str(os.path.basename(loc))
    res = os.path.join(results_path, f"EA_processed_{project_name}")
    final_df.to_csv(res)
