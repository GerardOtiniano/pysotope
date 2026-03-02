# src/OSIBL_correction
import os

from .utils.corrections.drift import *
from .utils.corrections.linearity import *
from .utils.corrections.methanol import *
from .utils.corrections.vsmow import *
from .utils.outliers.outliers import *
from .utils.queries import *
from .utils.regression import *
from .utils.uncertainty_and_output import *
from .utils.figures import *
from .utils.base_functions import *
from .utils.config import CorrectionConfig
from .utils.corrections.pame import *

def iso_process(user_linearity_conditions = False, min_area_threshold = None, include_parabolic=False):
    """
    Main processing pipeline for GC-IRMS isotope data.

    This function executes the complete correction workflow for
    compound-specific isotope measurements, including:

        • Instrument drift correction
        • Linearity (area-dependent) correction
        • Reference standard normalization (e.g., VSMOW or VPDB)
        • Optional methanol/PAME correction
        • Outlier removal
        • Replicate averaging with propagated uncertainty
        • Figure generation and export
        • Final results export to disk

    The function is file-based and does not return data objects.
    All intermediate and final results are written to structured
    output folders along with diagnostic figures and a processing log.

    Parameters
    ----------
    user_linearity_conditions : dict or bool, default=False
        Optional user-defined configuration for linearity correction.
        If provided, overrides automatically determined regression
        conditions (e.g., selected standards or model form).

    min_area_threshold : float or None, default=None
        Minimum peak area required for inclusion in processing.
        Peaks below this threshold are removed prior to correction.
        If None, no area-based filtering is applied.

    include_parabolic : bool, default=False
        If True, includes a second-order (parabolic) term in the
        linearity correction model. If False, a linear model is used.

    Outputs
    -------
    This function writes the following to disk:

    • Corrected (drift, linearity, reference standard, methylation) sample results
    • Mean of corrected sample results
    • Final processed dataset with propagated uncertainties
    • Diagnostic plots (drift, linearity, standards, outliers)
    • Processing log file documenting all correction steps
    • Correction Figures

    No objects are returned.

    Notes
    -----
    Processing sequence:

        1. Isotope system selection
        2. Output directory creation
        3. Standards loading
        4. Data import and classification
        5. Area-based filtering (optional)
        6. Drift correction (time-based regression)
        7. Linearity correction (area-based regression)
        8. Reference standard normalization
        9. Optional methanol/PAME correction
        10. Optional hydrogen or carbon isotope calculation of methanol from PAME
        11. Outlier removal
        12. Replicate averaging with uncertainty propagation
        13. Final export and figure generation

    Uncertainty Propagation
    -----------------------
    Analytical uncertainty is propagated analytically through:

        • Drift regression parameter covariance
        • Linearity regression parameter covariance
        • Reference standard scaling uncertainty
        • Methanol correction uncertainty (if applicable)
        • Replicate variability

    All variance components are combined using first-order
    Taylor series propagation and reported in the final dataset.
    """
     
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import os
    from IPython.display import clear_output
    from scipy.stats import linregress
    from matplotlib.dates import date2num
    from datetime import datetime, timedelta
    import time
    from scipy.stats import zscore
    from sklearn.linear_model import HuberRegressor
    import scipy.stats as stats

    # Query isotope system
    isotope = isotope_type()
    cfg = CorrectionConfig(isotope)

    # Setup output folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder(isotope)

    # Set standards
    standards_df = load_standards(isotope)#query_stds(alt_stds, isotope)
    append_to_log(log_file_path, standards_df)

    # Import data
    lin_std, drift_std, samples, correction_log, pame = import_data(loc, folder_path, log_file_path, isotope, standards_df)
    if min_area_threshold is not None:
        lin_std = lin_std[lin_std["area"]>min_area_threshold]
        drift_std = drift_std[drift_std["area"]>min_area_threshold]
        samples = samples[samples["area"]>min_area_threshold]
    uncorrected_samples = samples.copy()

    # Run standard plots for area
    std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path,isotope=isotope, dD=isotope)

    # Drift Correction
    samples, lin_std, drift_std, dD_temp, correction_log = process_drift_correction(cfg, samples, lin_std, drift_std, correction_log, log_file_path=log_file_path, fig_path=fig_path,isotope=isotope)

    # # Show plots again
    # std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path, dD=dD_temp,isotope=isotope)

    # Linearity (area) correction
    drift_std, correction_log, lin_std, samples = process_linearity_correction(cfg, samples, drift_std, lin_std, dD_temp, correction_log,
                                                                               folder_path, fig_path, isotope, user_linearity_conditions,
                                                                               log_file_path=log_file_path, include_parabolic=include_parabolic)

    # Reference standard correction
    samples, standards = reference_standard_correction(cfg, samples, lin_std, drift_std, correction_log, folder_path, fig_path, log_file_path, isotope, standards_df)

    # PAME
    if pame:
        samples, pame_unknown, pame_derv_meth, pame_derv_meth_unc  = calculate_methanol_dD(cfg, samples, isotope, log_file_path)
        samples, standards = q_methylation(samples, standards, isotope, log_file_path, 
                                           meth_val = pame_derv_meth, meth_std = pame_derv_meth_unc); # methanol values from PAME
    else: # methylation correction without PAME
        samples, standards = q_methylation(samples, standards, isotope, log_file_path);

    # Remove outliers
    samples, excluded_samples = outlier_removal(samples, fig_path, log_file_path, isotope)
    raw_samples = samples
    
    print("standards", list(standards))
    print("samples", list(samples))
    # Calculate mean values of replicate analyses
    samples = mean_values_with_uncertainty(samples, cfg, sample_name_header="Identifier 1", chain_header="chain", iso=isotope)
    if pame:
        pame_unknown = mean_values_with_uncertainty(pame_unknown, cfg = cfg, sample_name_header="Identifier 1", chain_header="chain", iso=isotope)
    else:
        pame_unknown = None
    # Final Data Correction and Plot
    output_results(raw_samples, samples, standards, pame_unknown, folder_path, fig_path, results_path, isotope, pame, log_file_path, cfg)


