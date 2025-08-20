#src/pysotope/EA/utils/
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import os

from . .base_functions import append_to_log
from .config import CorrectionConfig

log_file_path = fig = None

def drift_correction(df, log_file, fig_dir):
    """
    Perform isotopic drift correction on isotopic ratios for Nitrogen and/or Carbon.

    This function:
    1. Filters 'SORGHUM' standards and builds drift correction models for Nitrogen
       and Carbon using those standards.
    2. Prompts the user to choose whether to apply correction for Nitrogen,
       Carbon, both, or neither.
    3. Applies the selected corrections to the corresponding isotope ratio columns
       in the DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements.
    log_file : str
        Path to the log file where processing details will be logged.
    fig_dir : str
        Path to the directory for storing any figures that are generated.

    Returns
    -------
    tuple
        df : pandas.DataFrame
            The EA-IRMS DataFrame with drift-corrected isotope ratios.
            ('d 15N/14N_corr', 'd 15N/14N_se', 'd 13C/12C_corr', 'd 13C/12C_se').
        cfg : CorrectionConfig
            A configuration object holding boolean values for which corrections were applied.

    Notes
    -----
    - Makes log_file and fig_dir global for access by all functions.
    - If neither correction is applied, the DataFrame is returned unchanged
      (except for initialized correction columns).
    - Necessary details are logged into the log file.
    """

    global log_file_path
    global fig
    fig = fig_dir
    log_file_path = log_file
    append_to_log(log_file_path, "\n\n\nDrift Correction:")

    # Build model
    N, C = get_sorghum(df)
    drop_standards(N, df, 'N')
    N_model, N_relative = drift_model(N, "N")
    drop_standards(C, df, 'C')
    C_model, C_relative = drift_model(C, "C")

    # Confirm application with user
    N_con, C_con = drift_confirm()
    cfg = CorrectionConfig(drift_N_applied = N_con, drift_C_applied = C_con)
    df["d 15N/14N_corr"] = df["d 15N/14N_se"] = df["d 13C/12C_corr"] = df["d 13C/12C_se"] = np.nan
    if not (N_con or C_con):
        append_to_log(log_file_path, "No drift correction applied.")
        return df, cfg
    else:
        if N_con:
            df = apply_drift_model(df, N_model, N_relative, 'N')
        if C_con:
            df = apply_drift_model(df, C_model, C_relative, 'C')
        return df, cfg

def drift_confirm():
    """
    Prompt the user to choose whether to apply correction for Nitrogen, Carbon, both, or
    neither.

    Returns
    ----------
    tuple:
        N_corr_applied : bool
            True if Nitrogen correction is selected, False otherwise.
        C_corr_applied : bool
            True if Carbon correction is selected, False otherwise.

    Notes
    -----
    - The function will repeatedly prompt until a valid option (1–4) is entered.
    - The choice made by the user is logged into the log file.
    """

    input_user = False
    while not(input_user):
        corr_sel = input("\nChoose which corrections to apply:\n1) Do not apply any correction.\n2) Apply correction for Nitrogen.\n3) Apply correction for Carbon.\n4) Apply correction for both Nitrogen and Carbon.\n Enter 1, 2, 3 or 4:")
        input_user = True
        if corr_sel == '1':
            append_to_log(log_file_path, "\nUser choice: 1) Do not apply any correction.")
            N_corr_applied = False
            C_corr_applied = False
        elif corr_sel == '2':
            append_to_log(log_file_path, "\nUser choice: 2) Apply correction for Nitrogen.")
            N_corr_applied = True
            C_corr_applied = False
        elif corr_sel == '3':
            append_to_log(log_file_path, "\nUser choice: 3) Apply correction for Carbon.")
            N_corr_applied = False
            C_corr_applied = True
        elif corr_sel == '4':
            append_to_log(log_file_path, "\nUser choice: 4) Apply correction for both Nitrogen and Carbon.")
            N_corr_applied = True
            C_corr_applied = True
        else:
            input_user = False
            print("Wrong input. Try again.")
    return N_corr_applied, C_corr_applied

def get_sorghum(df):
    """
    Extracts 'SORGHUM' standards for both Nitrogen and Carbon from the EA-IRMS dataset.

    Parameters
    ----------
    df : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements.

    Returns
    -------
    tuple
        sorghum_N : pandas.DataFrame
            EA-IRMS data filtered with 'SORGHUM' standard measurements only for Nitrogen.
        sorghum_C : pandas.DataFrame
            EA-IRMS data filtered with 'SORGHUM' standard measurements only for Carbon.

    Notes
    -----
    - Raises ValueError if no 'SORGHUM' standards are found in the EA-IRMS dataset.
    - Whether 'SORGHUM' standards are found or not is logged into the log file.
    """

    sorghum = df[df['Identifier 1'].str.lower() == 'sorghum'].copy()
    if sorghum.empty:
        append_to_log(log_file_path, "\nError: No 'SORGHUM' standards found in the data.")
        raise ValueError("No 'SORGHUM' standards found.")

    sorghum_N = sorghum[sorghum['Element Type'] == 'Nitrogen'] # Only consider nitrogen for sorghum
    sorghum_C = sorghum[sorghum['Element Type'] == 'Carbon']
    append_to_log(log_file_path, "\n'SORGHUM' standards found in the data.")
    return sorghum_N, sorghum_C

def drop_standards(std_df, df, tag):
    """
    Display the relevant 'SORGHUM' standard measurements (Nitrogen or Carbon) and
    prompts the user to enter the indices of rows to remove. The specified
    rows are dropped from both the standards DataFrame and the main dataset.

    Parameters
    ----------
    std_df : pandas.DataFrame
        DataFrame containing only the 'SORGHUM' standard measurements for the selected
        isotope (Nitrogen or Carbon).
    df : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements.
    tag : str
        Indicator for which isotope is being processed:
        - 'N' → Nitrogen
        - 'C' → Carbon

    Notes
    -----
    - If the user enters "None", no rows are removed.
    - Indices must be entered as integers, separated by commas (e.g., "210,416").
    - Removed indices are logged into the log file.
    """

    if tag == 'N':
        col = "d 15N/14N"
        el = 'Nitrogen'
    elif tag == 'C':
        col = "d 13C/12C"
        el = 'Carbon'
    sub = std_df[["Identifier 1", "Time", col]].copy()
    print('')
    print(sub)
    idx_str = input("Enter indices that you wish to remove in a comma-separated manner. (Enter None if no values are to be removed)\n ")
    if idx_str:
        idxs = idx_str.split(",")
        idxs = [int(idx) for idx in idxs]
        std_df.drop(idxs, inplace=True)
        df.drop(idxs, inplace=True)
        msg = f"{el} values at {idx_str} have been removed."
        print(msg)
        append_to_log(log_file_path, msg)

def get_isotope(tag):
    """
    Map an isotope tag to its data column name, plot label, and element name.

    Parameters
    ----------
    tag : str
        Tag for the isotope for which information needs to be retrieved.
        - 'N' → Nitrogen
        - 'C' → Carbon

    Returns
    -------
    tuple
        col_name : str
            Name of the column in the dataset containing corresponding isotope ratio
            (e.g., "d 15N/14N", "d 13C/12C").
        plot_label : str
            Formatted plot label for corresponding isotope ratio.
        element : str
            Name of the chemical element associated with the inputted isotope tag
            ("Nitrogen" or "Carbon").
        component : str
            Component associated with the element. (N2 for 'N', CO2 for 'C')

    Notes
    -----
    - Isotope tag matching is case-insensitive.
    - Raises ValueError if the tag entered isn't 'N' or 'C'.
    - If incorrect tag is passed as an argument, logs it into the log file.
    """

    if "N" in tag.upper():
        return "d 15N/14N", r"$\delta^{15}\mathrm{N}$", "Nitrogen", "N2"
    elif "C" in tag.upper():
        return "d 13C/12C", r"$\delta^{13}\mathrm{C}$", "Carbon", "CO2"
    else:
        append_to_log(log_file_path, f"Unknown isotope tag: {tag}")
        raise ValueError(f"Unknown isotope tag: {tag}")
    
def drift_model(df, tag):
    """
    Build and evaluate a linear drift correction model for the specified isotope.

    Builds an ordinary least squares (OLS) linear regression model of the isotope
    ratio versus time (seconds since start) using 'SORGHUM' standard values. Plots the
    drift correction applied to 'SORGHUM' to display the correction to the user.

    Parameters
    ----------
    df : pandas.DataFrame
        EA-IRMS data filtered with 'SORGHUM' standard measurements for the specified isotope.
    tag : str
        Isotope tag to specify which isotope to model. Used to retrieve the corresponding
        isotope columns and labels.
        - 'N' → Nitrogen
        - 'C' → Carbon

    Returns
    -------
    tuple
        model : statsmodels.regression.linear_model.RegressionResultsWrapper
            The fitted OLS regression model object.
        mean_value : float
            The mean of the original isotope ratio values before correction.

    Notes
    ------------
    - Prints the linear regression equation and model statistics to the console
    - Logs the linear regression equation and model statistics to the log file.
    - Calls 'plot_drift_correction' to generate for the drift correction.
    """

    iso_col, iso_label, el, comp = get_isotope(tag)
    append_to_log(log_file_path, f"\nDrift Correction Model for {el}: ")

    # Model derivation
    secs = df["Seconds Since Start"].to_numpy()
    X    = sm.add_constant(secs.reshape(-1, 1))
    y    = df[iso_col].to_numpy()

    model   = sm.OLS(y, X).fit()
    y_hat   = model.predict(X)
    y_corr  = y - (y_hat - y.mean())
    
    # Linear equation information
    intercept = model.params[0]
    slope     = model.params[1]
    
    # stats
    adj_r2 = model.rsquared_adj
    rmse   = np.sqrt(model.mse_resid)
    fit_line = f"y = {slope:.2f}·x + {intercept:.2f}"
    r_sq = f"{adj_r2:.3f}"
    rmse_msg = f"{rmse:.3f}"
    print(f"\n{iso_col}:  {fit_line}", f"Adj. R$^{{2}}$ = {r_sq}", f"RMSE = {rmse_msg}", sep='\n', end='\n')
    append_to_log(log_file_path, f"Equation of fitted line: {fit_line}")
    append_to_log(log_file_path, f"R-squared: {r_sq}")
    append_to_log(log_file_path, f"Root Mean Squared Error (RMSE): {rmse_msg}")

    plot_drift_correction(secs, y, y_corr, iso_label, el)
    return model, y.mean()
    
def plot_drift_correction(X, y, y_corr, iso_label, el):
    """
    Plot uncorrected and corrected isotope ratio values for 'SORGHUM' standard over time
    to visualize drift correction.

    Creates a scatter plot showing uncorrected isotope measurements and drift-corrected
    values against elapsed time ("Seconds Since Start"). Saves the figure to a file
    and displays it.

    Parameters
    ----------
    X : array-like
        Time points in seconds since the start at which measurements were taken.
    y : array-like
        Uncorrected isotope ratio values.
    y_corr : array-like
        Drift-corrected isotope ratio values.
    iso_label : str
        Formatted plot label for axis with isotope ratio values.
    el : str
        Name of the chemical element whose isotope ratio values are being processed
        and plotted ("Nitrogen" or "Carbon").

    Notes
    -----
    - Saves the plot as `"Drift_<Element>.png"` in the directory specified by the global
      `fig` variable.
    """

    plt.figure(figsize=(6, 4))

    plt.scatter(X, y,      c="k",   s=200, alpha=0.5,
                label=f"Uncorrected {iso_label}")
    plt.scatter(X, y_corr, c="red", s=200, alpha=0.5, edgecolor="k",
                label=f"Corrected {iso_label}")

    plt.xlabel("Seconds Since Start")
    plt.ylabel(iso_label)
    plt.legend()

    plt.tight_layout()

    fig_name = f"Drift_{el}.png"
    fig_path = os.path.join(fig, fig_name)
    plt.savefig(fig_path)

    plt.show()
    
def apply_drift_model(df_corr, model, rel, tag):
    """
    Apply isotopic drift correction to the EA-IRMS DataFrame using a fitted linear
    regression model.

    Uses the isotope tag to determine which isotope column and component to process.
    Predicts drift based on elapsed time, corrects isotope ratio values, and
    records prediction uncertainty. Updates the EA-IRMS DataFrame in place.

    Parameters
    ----------
    df_corr : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements.
    model : statsmodels.regression.linear_model.RegressionResultsWrapper
        Fitted linear regression model used to predict drift values.
    rel : float
        Relative mean isotope value used as a baseline for drift correction.
    tag : str
        Isotope tag used to identify the isotope for correction.
        - 'N' → Nitrogen
        - 'C' → Carbon

    Returns
    -------
    pandas.DataFrame
        Updated EA-IRMS DataFrame with new columns added for the isotope:
        - "<isotope_column>_corr": corrected isotope ratios,
        - "<isotope_column>_se": standard error of predictions.

    Notes
    -----
    - Prints the number of corrected rows to the console.
    - Logs the correction summary to the log file.
    """

    col, iso_label, el, comp = get_isotope(tag)
    secs = df_corr["Seconds Since Start"].to_numpy().reshape(-1, 1)
    mask = (df_corr["Component"] == comp) & df_corr[col].notna()
    if mask.any():
        X      = sm.add_constant(secs[mask])
        pred   = model.get_prediction(X)
        drift  = pred.predicted_mean - rel
        se     = pred.se_mean
        df_corr.loc[mask, f"{col}_corr"]  = (df_corr.loc[mask, col].to_numpy() - drift) # Corrected values
        df_corr.loc[mask, f"{col}_se"]    = se #Prediction uncertainty

    msg = f"\nDrift correction applied to {el}:\n  {tag} rows corrected: {mask.sum()}"
    print(msg)
    append_to_log(log_file_path, msg)

    return df_corr