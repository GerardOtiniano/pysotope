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




def iso_process(pame=False, user_linearity_conditions = False):
    cfg = CorrectionConfig()
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

    """
    to do:
        - remove isotope argument throughout
    """
    # Setup output folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder()

    # Set standards - Gerard
    standards_df = load_standards(isotope)
    append_to_log(log_file_path, standards_df)

    # Import data - Gerard
    std, samples, correction_log = import_data(loc, folder_path, log_file_path, standards_df)
    uncorrected_samples = samples.copy() # Function will return standard dataframe, sample dataframe, and correction_log

    # Run standard plots for area
    std_plot(std, folder_path=folder_path, fig_path=fig_path) # Modify to work with standards from import data

    # Drift Correction
    samples, std, dD_temp, correction_log = process_drift_correction(cfg, samples, std, correction_log, log_file_path=log_file_path, fig_path=fig_path) # Linear regression

    # VPD correction
    samples, standards = vsmow_correction(cfg, samples, std, correction_log, folder_path, fig_path, log_file_path)

    # Remove outliers
    samples, excluded_samples = outlier_removal(samples, fig_path, log_file_path)
    raw_samples = samples

    # Calculate mean values of replicate analyses
    samples = mean_values_with_uncertainty(samples, cfg, sample_name_header="Identifier 1", chain_header="chain")

    # Final Data Correction and Plot
    output_results(raw_samples, samples, standards, pame_unknown, folder_path, fig_path, results_path, log_file_path, cfg)





