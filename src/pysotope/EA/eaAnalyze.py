# src/pyosotope/EA
import os

from .base_functions import create_folder, append_to_log
from .utils.VPDB_correction import VPDB_correction
from .utils.uncertainty_calculation import uncertainty_calculation
from .utils.import_data import load_ea_standards, import_EA_data
from .utils.ea_drift_correction import drift_correction

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