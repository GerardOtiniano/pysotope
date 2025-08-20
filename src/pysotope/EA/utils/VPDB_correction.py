import os
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd

from . .base_functions import append_to_log

EA_stds = ["Sorghum", "Acetanilide", "Wheat Flour_NML", "BL Sediment"]
identifiers = ['sorghum', 'acetanilide', 'wheat flour', "bl sediment", "sample"]
log_file_path = fig = None

def VPDB_correction(df, standards, cfg, log_file, fig_dir):
    """
    Perform VPDB calibration on isotopic ratios for Nitrogen and Carbon.

    This function:
    1. Splits the dataset by element type (Nitrogen and Carbon).
    2. Retrieves a DataFrame with the values for the three standards "SORGHUM",
       "ACETANILIDE" and "WHEAT FLOUR", along with their actual values for both
       elements.
    3. Retrieves a DataFrame with the values for the "BL SEDIMENT" along with
       the actual values for both elements.
    4. Builds VPDB calibration models for Nitrogen and Carbon based on the
       "SORGHUM", "ACETANILIDE", "WHEAT FLOUR" standards and tests it
       against the "BL SEDIMENT" samples.
    5. Applies the calibration model to the respective datasets.
    6. Combines back the calibrated Nitrogen and Carbon datasets and sorts
       it by measurement time.

    Parameters
    ----------
    df : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements (Post drift correction).
    standards : pandas.DataFrame
        EA-IRMS standards DataFrame containing isotopic information about the standards.
    cfg : CorrectionConfig
        A configuration object holding boolean values for which corrections were applied.
    log_file : str
         Path to the log file where processing details will be logged.
    fig_dir : str
        Path to the directory for storing any figures that are generated.

    Returns
    -------
    pandas.DataFrame
        Combined calibrated EA-IRMS dataset for Nitrogen and Carbon isotopes,
        sorted by "Seconds Since Start".
    """

    global log_file_path
    global fig
    fig = fig_dir
    log_file_path = log_file

    msg = "\n\nVPDB Calibration:\n"
    append_to_log(log_file_path, msg)
    print(msg)

    N_df = df[df['Element Type'] == 'Nitrogen'].copy() #EA-IRMS DataFrame for Nitrogen
    C_df = df[df['Element Type'] == 'Carbon'].copy() #EA-IRMS DataFrame for Carbon

    N_std_df = standardize(N_df, standards, 'N', 'std') #Standards DataFrame for Nitrogen
    C_std_df = standardize(C_df, standards, 'C', 'std') #Standards DataFrame for Carbon

    N_bl_df = standardize(N_df, standards, 'N', 'BL') #Testing sample DataFrame for Nitrogen
    C_bl_df = standardize(C_df, standards, 'C', 'BL') #Testing sample DataFrame for Carbon

    N_model, N_rel, N_slope, N_int = VPDB_model(N_std_df, cfg, 'N', N_bl_df)
    C_model, C_rel, C_slope, C_int = VPDB_model(C_std_df, cfg, 'C', C_bl_df)

    # Calibrated EA-IRMS DataFrame for Nitrogen
    cal_N_df = apply_vpdb_calibration(N_df, N_model, N_rel, cfg, 'N', N_slope, N_int)

    # Calibrated EA-IRMS DataFrame for Carbon
    cal_C_df = apply_vpdb_calibration(C_df, C_model, C_rel, cfg, 'C', C_slope, C_int)

    # Calibrated EA-IRMS DataFrame (Nitrogen and Carbon combined)
    final_df = pd.concat([cal_N_df, cal_C_df], ignore_index=False)
    final_df.sort_values(by="Seconds Since Start", inplace=True)

    return final_df

def standardize(df, standards, tag, option):
    """
    Create standards or testing sample DataFrame with actual values and add those actual values
    to the main EA-IRMS DataFrame.

    Depending on the 'tag' ('N' or 'C') and the 'option' ('std' or 'BL'), this function creates standard
    or testing sample DataFrame with corresponding values based on the element being processed. Also adds
    the actual values retrieved from the EA-IRMS standards DataFrame into the main EA-IRMS DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        The EA-IRMS dataset containing isotopic measurements (Post drift correction).
    standards : pandas.DataFrame
        EA-IRMS standards DataFrame containing isotopic information about the standards.
    tag : str
        Indicator for which isotope is being processed:
        - 'N' → Nitrogen
        - 'C' → Carbon
    option : str
        Indicator for which DataFrame needs to be created.
        - 'std': Standards DataFrame ("SORGHUM", "ACETANILIDE", "WHEAT FLOUR")
        - 'BL': Testing sample DataFrame ("BL SEDIMENT")

    Returns
    -------
    pandas.DataFrame
        Subset of the input DataFrame, either standards or testing sample created based on
        option type with actual values retrived from EA-IRMS standards DataFrame.

    Notes
    ------------
    - Modifies the input DataFrame in place by adding actual isotopic ratio values for standards
      or testing sample.
    - Raises ValueError if any standard or the testing sample is not found in the EA-IRMS dataset.
    - Whether any standard or the testing sample is not found is logged into the log file.
    """

    if tag == 'N':
        col = 'd15N(AIR) value'
        el = 'Nitrogen'
    elif tag == 'C':
        col = 'd13C(VPDB) value'
        el = 'Carbon'

    if option == "std":
        std_str = ""
        for i in range(0,3):
            std = df[df['Identifier 1'].str.lower() == identifiers[i]].copy()
            if std.empty:
                append_to_log(log_file_path, f"Error: No '{identifiers[i].upper()}' standards found for {el} in the data.")
                raise ValueError(f"No '{identifiers[i].upper()}' standards found for {el}.")

            std_val = standards.loc[standards['EA-IRMS Standards'] == EA_stds[i], col].values[0]
            df.loc[df['Identifier 1'].str.lower() == identifiers[i], col] = std_val
            std_str += f"'{identifiers[i].upper()}',"
        df_out = df[df['Identifier 1'].str.lower().isin(identifiers[0:3])].copy()
        append_to_log(log_file_path, f"{std_str.rstrip(',')} standards found for {el} in the data.")

    elif option == "BL":
        bl = df[df['Identifier 1'].str.lower() == identifiers[3]].copy()
        if bl.empty:
            append_to_log(log_file_path, f"\nError: No '{identifiers[3].upper()}' samples found for {el} in the data.")
            raise ValueError(f"No '{identifiers[3].upper()}' samples found for {el}.")

        std_val = standards.loc[standards['EA-IRMS Standards'] == EA_stds[3], col].values[0]
        df.loc[df['Identifier 1'].str.lower() == "bl sediment", col] = std_val
        df_out = df[(df['Identifier 1'].str.lower() == identifiers[3])].copy()
        append_to_log(log_file_path, f"'{identifiers[3].upper()}' samples found for {el} in the data.")

    return df_out

def get_isotope(tag, cfg):
    """
    Map an isotope tag to its data column name from the configuration, output column name,
    actual column name, plot label, and element name.

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
            Name of the column in the dataset containing corresponding isotope ratio.
            Extracted from 'cfg' depending on drift correction.
        output_col : str
            Name of the columns where the VPDB calibrated values will be outputted into
            (e.g., "VPDB_d15N/14N", "VPDB_d13C/12C").
        actual_col : str
            Name of the column containing actual standard values
            (e.g., 'd15N(AIR) value', 'd13C(VPDB) value').
        plot_label : str
            Formatted plot label for corresponding isotope ratio.
        element : str
            Name of the chemical element associated with the inputted isotope tag
            ("Nitrogen" or "Carbon").

    Notes
    -----
    - Isotope tag matching is case-insensitive.
    - Raises ValueError if the tag entered isn't 'N' or 'C'.
    - If incorrect tag is passed as an argument, logs it into the log file.
    """

    if "N" in tag.upper():
        return cfg.dN_col, "VPDB_d15N/14N", 'd15N(AIR) value', r"$\delta^{15}\mathrm{N}$", "Nitrogen"
    elif "C" in tag.upper():
        return cfg.dC_col, "VPDB_d13C/12C", 'd13C(VPDB) value', r"$\delta^{13}\mathrm{C}$", "Carbon"
    else:
        append_to_log(log_file_path, f"Unknown isotope tag: {tag}")
        raise ValueError(f"Unknown isotope tag: {tag}")

def VPDB_model(df, cfg, tag, bl_df):
    """
    Build and evaluate a VPDB calibration model for the specified isotope.

    Fits an ordinary least squares (OLS) linear regression model of actual isotopic ratio
    versus measured isotopic ratio using "SORGHUM", "ACETANILIDE", "WHEAT FLOUR" standards.
    Applies the model on the "BL SEDIMENT" testing sample DataFrame to test the model. Plots
    the calibration applied to the standards and testing sample to display the calibration
    to the user.

    Parameters
    ----------
    df : pandas.DataFrame
        Standards DataFrame containing isotope measurements for the three standards along
        with actual values for the specified isotope.
    cfg : CorrectionConfig
        Configuration object that provides isotope column names based on drift correction.
    tag : str
        Isotope tag to specify which isotope to model. Used to retrieve the corresponding
        isotope columns and labels.
        - 'N' → Nitrogen
        - 'C' → Carbon
    bl_df : pandas.DataFrame
        Testing sample DataFrame containing isotope measurements for the "BL SEDIMENT"
        along with actual values for the specified isotope.

    Returns
    -------
    tuple
        model : statsmodels.regression.linear_model.RegressionResultsWrapper
            The fitted OLS regression model object.
        mean : float
            Combined mean of the actual isotope values from standards and testing sample.

    Side Effects
    ------------
    - Prints model statistics and baseline testing results to the console.
    - Logs detailed calibration information to the global log file.
    - Calls `plot_VPDB_calibration` to generate calibration plots.

    Notes
    -----
    - Bias and RMSE are calculated on testing sample to assess model accuracy.
    - Prints the linear regression equation and model statistics to the console
    - Logs the linear regression equation and model statistics to the log file.
    - Calls 'plot_drift_correction' to generate calibration plots.
    """

    input_col, output_col, actual_col, iso_label, el = get_isotope(tag, cfg)
    append_to_log(log_file_path, f"\nVPDB Calibration Model for {el}: ")

    #Model derivation
    vals = df[input_col].to_numpy()
    X = sm.add_constant(vals.reshape(-1, 1))
    y = df[actual_col].to_numpy()

    model = sm.OLS(y, X).fit()

    y_pred = model.predict(X)

    df[output_col] = y_pred

    # Linear equation information
    intercept = model.params[0]
    slope = model.params[1]

    # stats
    adj_r2 = model.rsquared_adj
    rmse = np.sqrt(model.mse_resid)

    fit_line = f"y = {slope:.2f}·x + {intercept:.2f}"
    r_sq = f"{adj_r2:.3f}"
    rmse_msg = f"{rmse:.3f}"

    print(f"{output_col}:  {fit_line}", f"Adj. R$^{{2}}$ = {r_sq}", f"RMSE = {rmse_msg}", sep='\n', end='\n\n')
    append_to_log(log_file_path, f"Equation of fitted line: {fit_line}")
    append_to_log(log_file_path, f"R-squared: {r_sq}")
    append_to_log(log_file_path, f"Root Mean Squared Error (RMSE): {rmse_msg}\n")

    append_to_log(log_file_path, f"Testing VPDB Calibration on '{identifiers[3].upper()}' samples:")
    vals_bl = bl_df[input_col].to_numpy()
    y_bl_actual = bl_df[actual_col].to_numpy()
    X_bl = sm.add_constant(vals_bl.reshape(-1, 1))
    y_bl_pred = model.predict(X_bl)

    bl_df[output_col] = y_bl_pred

    # stats
    bl_bias = np.mean(y_bl_pred - y_bl_actual)
    bl_rmse = np.sqrt(np.mean((y_bl_pred - y_bl_actual) ** 2))

    mean_pred = f"{np.mean(y_bl_pred):.3f}"
    bias = f"{bl_bias:.3f}"
    rmse_msg = f"{bl_rmse:.3f}"

    print(f"BL Sediment {output_col}", f"Mean Prediction: {mean_pred}", f"Bias: {bias}", f"RMSE: {rmse_msg}", sep='\n', end='\n\n')
    append_to_log(log_file_path, f"Mean Prediction: {mean_pred}")
    append_to_log(log_file_path, f"Bias: {bias}")
    append_to_log(log_file_path, f"Root Mean Squared Error (RMSE): {rmse_msg}")

    mean = (y.sum() + y_bl_actual.sum()) / (y.size + y_bl_actual.size)

    combined_df = pd.concat([df, bl_df], ignore_index=False)

    plot_VPDB_calibration(combined_df, identifiers[0:4], el, actual_col, output_col, iso_label, slope, intercept)

    return model, mean, slope, intercept

def plot_VPDB_calibration(df, identifiers, el, actual_col, output_col, iso_label, slope = None, intercept = None, input_col = None):
    """
    Plot predicted versus measured/actual isotope ratio values to visualize VPDB calibration.

    Creates a scatter plot comparing predicted isotope ratio values from the calibration
    model against measured/actual values for EA_IRMS data. Includes a 1:1 reference line
    and uses different markers and colors for standards, testing samples, and target samples.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing isotope data including predicted and measured/actual values.
    identifiers : list of str
        List containing identifiers for different standards, testing sample and target sample category.
    el : str
        Name of the chemical element associated with the inputted isotope tag
        ("Nitrogen" or "Carbon").
    actual_col : str
        Name of the column containing actual standard values
        (e.g., 'd15N(AIR) value', 'd13C(VPDB) value').
    output_col : str
        Name of the columns where the VPDB calibrated values will be outputted into
        (e.g., "VPDB_d15N/14N", "VPDB_d13C/12C").
    iso_label : str
        Formatted plot label for corresponding isotope ratio.
    input_col : str, optional
        Name of the column in the dataset containing corresponding isotope ratio. Extracted from
        'cfg' depending on drift correction. Only provided, if plotting is being done for samples.

    Notes
    -----
    - The legend differentiates standard samples, testing sample, and target samples by marker style
      and color.
    - Saves the plot as "<id>_<el>.png" in the directory specified by the global 'fig' variable.
    """

    plt.figure(figsize=(6, 4))

    x = df[output_col].to_numpy()

    # model removed from parameters
    # x_line = np.linspace(x.min(), x.max(), 100)
    # X_line = sm.add_constant(x_line.reshape(-1, 1))
    # y_line = model.predict(X_line)
    # plt.plot(x_line, y_line, linestyle='--', color='k')

    # plt.plot([x.min()-1, x.max()+1], [x.min()-1, x.max()+1], linestyle='--', color='k')
    if slope:
        y1 = x.min()*slope + intercept
        y2 = x.max()*slope + intercept
        plt.plot([x.min(),x.max()], [y1,y2], linestyle='--', color='k')

    colors = ["orange", "blue", "green", "red", "k"]
    for i, identifier in enumerate(identifiers):
        sub_df = df[df['Identifier 1'].str.lower() == identifier].copy()
        x = sub_df[output_col].to_numpy()
        y = sub_df[actual_col].to_numpy()
        label = identifier.upper()
        ec = "k"
        id = "VPDB_Model"
        if identifier == 'bl sediment':
            c = colors[3]
            m = 's'
            label = f"TEST: {identifier.upper()}"
        elif identifier == 'sample':
            sub_df = df[~df['Identifier 1'].str.lower().isin(identifiers[0:4])].copy()
            x = sub_df[output_col].to_numpy()
            y = sub_df[input_col].to_numpy()
            c = colors[4]
            m = 'x'
            ec = None
            id = "VPDB"
        else:
            c = colors[i]
            m = 'o'

        plt.scatter(x, y, label=label, marker=m, color=c, alpha=0.6, s=100, edgecolor=ec)

    plt.xlabel(f"Predicted {iso_label} (‰, VPDB)")
    plt.ylabel(f"Measured/Actual {iso_label} (‰, VPDB)")
    plt.legend()
    plt.tight_layout()

    fig_name = f"{id}_{el}.png"
    fig_path = os.path.join(fig, fig_name)
    plt.savefig(fig_path)

    plt.show()

def apply_vpdb_calibration(df, model, rel, cfg, tag, slope, intercept):
    """
    Apply VPDB calibration to the EA-IRMS DataFrame using a fitted linear regression model.

    Predicts calibrated isotope values and their uncertainties for the target samples based on
    specified isotope tag and records prediction uncertainty, updates the DataFrame with these
    values, and generates a calibration plot.

    Parameters
    ----------
    df : pandas.DataFrame
        The EA-IRMS DataFrame (either Nitrogen or Carbon) containing isotopic measurements.
    model : statsmodels.regression.linear_model.RegressionResultsWrapper
        Fitted linear regression model used for VPDB calibration.
    rel : float
        Relative mean isotope value used for calibration difference calculation.
    cfg : CorrectionConfig
        Configuration object that provides isotope column names based on drift correction.
    tag : str
        Isotope tag used to identify the isotope for calibration.
        - 'N' → Nitrogen
        - 'C' → Carbon

    Returns
    -------
    pandas.DataFrame
        Updated EA-IRMS DataFrame (either Nitrogen or Carbon) with new columns added:
        - "<output_col>": calibrated isotope values,
        - "<output_col>_diff": calibration difference from reference mean,
        - "<output_col>_se": standard error of the calibration predictions.

    Notes
    ------------
    - Prints the number of calibrated rows to the console.
    - Logs the calibration summary to the log file.
    - Calls `plot_VPDB_calibration` to generate a plot of calibration results.
    """

    append_to_log(log_file_path, "")
    input_col, output_col, actual_col, iso_label, el = get_isotope(tag, cfg)
    vals = df[input_col].to_numpy()
    X = sm.add_constant(vals.reshape(-1, 1))
    y = model.get_prediction(X)
    y_pred = y.predicted_mean
    cal_diff = y_pred - rel
    se = y.se_mean
    df[output_col] = y_pred
    df[f"{output_col}_diff"] = cal_diff
    df[f"{output_col}_se"] = se

    plot_VPDB_calibration(df, identifiers, el, actual_col, output_col, iso_label, slope, intercept, input_col)

    rows = df[output_col].notna().sum()
    msg = f"VPDB calibration applied to {el}:\n  {tag} rows calibrated: {rows}"
    print(msg)
    append_to_log(log_file_path, msg)

    if tag == 'N':
        print("")

    return df