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

def iso_process(pame=False, user_linearity_conditions = False, min_area_threshold = None, include_parabolic=False):
    """
     Main processing pipeline for isotope data measured on the OSIBL GC-IRMS.
    
     This function orchestrates the full correction workflow, including:
     chain assignment, drift correction, linearity correction, VSMOW scaling,
     optional methanol and PAME corrections, uncertainty propagation,
     figure generation, and final result export.
    
     Parameters
     ----------
     pame : bool, default=False
         If True, applies PAME-specific handling and enables methanol correction
         for methylated n-alkanoic acids. When False, standard processing
         for non-methylated compounds is performed.
    
     user_linearity_conditions : dict or None, default=None
         Optional user-specified configuration for linearity correction.
         If provided, overrides automatically detected linearity parameters.
         Expected to contain regression settings and/or standard selections
         used during linearity calibration.
    
     min_area_threshold : bool or float, default=False
         If False, no minimum peak area filtering is applied.
         If True, a default minimum area threshold is used.
         If a float is provided, it is interpreted as the minimum
         integrated peak area required for inclusion in processing.
    
     include_parabolic : bool, default=False
         If True, includes a parabolic (second-order) term in the linearity
         correction model. If False, a linear model is used.
    
     Returns
     -------
     raw_samples : pandas.DataFrame
         The original imported dataset with chain assignments and
         intermediate annotations added, prior to full correction.
    
     samples : pandas.DataFrame
         Fully corrected sample data including:
         - Drift-corrected isotope values
         - Linearity-corrected values
         - VSMOW-scaled isotope values
         - Optional methanol/PAME corrections
         - Propagated analytical uncertainties
         - Standard-normalized results
    
     standards : pandas.DataFrame
         Processed reference standards used in drift, linearity,
         and VSMOW normalization, including fitted parameters
         and residual diagnostics.
    
     pame_unknown : pandas.DataFrame or None
         If `pame=True`, contains methanol-corrected isotope values
         for methylated compounds. Otherwise, None.
    
     Notes
     -----
     This function controls the full correction sequence in the following order:
    
     1. Configuration loading
     2. Data import and chain assignment
     3. Outlier handling (if applicable)
     4. Drift correction
     5. Linearity correction
     6. VSMOW normalization
     7. Optional methanol/PAME correction
     8. Uncertainty propagation
     9. Figure generation and output export
    
     All intermediate correction steps are logged internally to ensure
     traceability and reproducibility.
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


