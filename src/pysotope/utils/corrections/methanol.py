import pandas as pd
import numpy as np
import datetime
import os

def append_to_log(log_file_path, log_message):
    """
    Add entry to log file.
    """
    with open(log_file_path, "a") as log_file:
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        initial_message = f"Log file created at {current_datetime}\n"
        log_file.write(log_message + "; " + str(current_datetime) + "\n")

def methyl_correction(unknown, stds, isotope, log_file_path, mdD = -72.5, mdD_err = 3.1):
    """
    Correct FAMES for δD of methyl groups introduced during methylation.
    Equations for correction based on Pilecky et al. (20221). Generalized as
    d2h = [(nH_total * d2H_deriv) - (nH_deriv * d2H_reagent)] / (nH_native)
    where,
        nH_tot      = total number of hydrogens
        nH_native   = number of hydrogens in underivatized compound
        nH_deriv    = hydrogens introduced by derivatizing reagent (methanol)
        d2H_deriv   = measured isotope value of derivatized product
        d2H_reagent = isotope composition of methnol

    ~GAO~ 12/4/2023
    """
    # Check if samples are empty
    mask = unknown['chain'] != "PAME"
    if unknown.empty or not mask.any():
        append_to_log(
            log_file_path,
            f"Methanol correction skipped: no derivatized compounds present.")
        # Ensure downstream compatibility
        unknown[f"methanol_{isotope}"] = np.nan
        unknown["methanol_error"] = np.nan
        return unknown
    # Extract the number of carbons from the 'chain' column
    c_n = unknown.loc[unknown['chain']!="PAME", 'chain'].str.extract(r'C(\d+)').astype(int).squeeze()
    # Apply the correction formula
    if isotope == "dD":
        unknown.loc[unknown['chain']!="PAME", f'methanol_{isotope}'] = ((unknown[f'RS_{isotope}'] * (2*c_n + 2)) - (mdD * 3)) / (2 * c_n)
        unknown.loc[unknown['chain']!="PAME",'methanol_error'] = mdD_err
        append_to_log(log_file_path, f"Methanol correction equation for {isotope}: [(RS_{isotope})(2 x #carbon + 2) - (methanol_{isotope}x3)] / [(2 x #carbon)+2]")

    elif isotope == "dC":
        unknown.loc[unknown['chain']!="PAME", f'methanol_{isotope}'] = (((c_n+1))*(unknown[f'RS_{isotope}'])-(mdD))/(c_n)
        unknown.loc[unknown['chain']!="PAME",'methanol_error'] = mdD_err
        append_to_log(log_file_path, f"Methanol correction equation for {isotope}: [(RS_{isotope})(#carbon + 1) - (methanol_{isotope})] / [#carbon]")
    return unknown
