# src/pyosotope/EA
import os
from base_functions import create_folder, append_to_log
from utils.import_data import load_ea_standards, import_EA_data
from utils.ea_drift_correction import drift_correction
from utils.VPD_correction import plot_measured_vs_actual

def ea_process():
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from IPython.display import clear_output
    from scipy.stats import linregress
    from matplotlib.dates import date2num
    from datetime import datetime, timedelta
    import time
    from scipy.stats import zscore
    from sklearn.linear_model import HuberRegressor
    import scipy.stats as stats

    # Setup output folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder()
    standards = load_ea_standards()
    append_to_log(log_file_path, standards)

    # Import data
    df = import_EA_data(loc)
    
    # Drift correction
    df = drift_correction(df) # Need to apply correction to Carbon
    df.to_csv("test.csv")
    
    # VPD correction
    # plot_measured_vs_actual(df, standards, 'Nitrogen')

# /Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv 