# src/pyosotope/EA
import os

from base_functions import create_folder, append_to_log
from src.pysotope.EA.utils.VPDB_correction import VPDB_correction
from src.pysotope.EA.utils.uncertainty_calculation import uncertainty_calculation
from utils.import_data import load_ea_standards, import_EA_data
from utils.ea_drift_correction import drift_correction

def ea_process():
    # Setup Output Folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder()
    standards = load_ea_standards()
    append_to_log(log_file_path, standards, True)

    # Import Data
    df = import_EA_data(loc)

    # Drift Correction
    df, cfg = drift_correction(df, log_file_path, fig_path)
    res_name = "Drift_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # VPDB Calibration
    df = VPDB_correction(df, standards, cfg, log_file_path, fig_path)
    res_name = "VPDB_Results.csv"
    res = os.path.join(results_path, res_name)
    df.to_csv(res)

    # Uncertainty Calculation
    final_df = uncertainty_calculation(df,cfg,log_file_path)
    project_name = str(os.path.basename(loc))
    res = os.path.join(results_path, f"EA_processed_{project_name}")
    final_df.to_csv(res)

#/Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv