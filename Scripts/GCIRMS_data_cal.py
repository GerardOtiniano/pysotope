import os
from correction.drift import *
from correction.linearity import *
from correction.methanol import *
from correction.vsmow import *
from outliers.raw_standards import *
from outliers.sample import *
from base_functions import *
from queries import *
from regression import *
from uncertainty_and_output import *
from figures import *


def GCIRMS_data_cal(directory=os.getcwd(), pame=False):
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
    project_name = query_project_name()
    folder_path, fig_path, results_path, loc, log_file_path = create_folder(project_name, isotope, directory)

    lin_std, drift_std, samples, correction_log = import_data(loc, folder_path, log_file_path, isotope=isotope)
    uncorrected_samples = samples.copy()

    # Run standard plots for area
    std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path, isotope=isotope)

    # Drift Correction
    samples, lin_std, drift_std, dD_temp, correction_log = process_drift_correction(samples, lin_std, drift_std, correction_log, log_file_path=log_file_path, fig_path=fig_path, isotope=isotope)

    # Show plots again
    std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path, dD=dD_temp, isotope=isotope)

    # Linearity (area) correction
    drift_std, correction_log, lin_std, samples = process_linearity_correction(samples, drift_std, lin_std, dD_temp, correction_log, folder_path, fig_path, isotope, log_file_path=log_file_path)

    # VSMOW correction
    samples, standards = vsmow_correction(samples, lin_std, drift_std, correction_log, folder_path, fig_path, log_file_path, isotope)

    # Methylation Correction
    if isotope == "dD":
        samples, standards = q_methylation(samples, standards, log_file_path)

    # Remove outliers
    samples, excluded_samples = outlier_removal(samples, fig_path, log_file_path)
    raw_samples = samples

    # Calculate mean values of replicate analyses
    samples = mean_values_with_uncertainty(samples, sample_name_header="Identifier 1", chain_header="chain")

    # Final Data Correction and Plot
    output_results(raw_samples, samples, standards, folder_path, fig_path, results_path, isotope)
