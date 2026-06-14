import os
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pandas as pd

from . .base_functions import append_to_log

log_file_path = fig = None


def _scale_name(tag):
    return "AIR" if "N" in tag.upper() else "VPDB"


def VPDB_correction(df, C_std, N_std, cfg, log_file, fig_dir):
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

    msg = "\n\nReference Standard Calibration:\n"
    append_to_log(log_file_path, msg)
    print(msg)

    #EA-IRMS DataFrame for Carbon
    C_df = df[df['Element Type'] == 'Carbon'].copy()
    C_vpdb_ids = _detect_vpdb_standards(C_df, C_std)
    C_std_df = standardize(C_df, C_std, 'C', 'std', C_vpdb_ids)
    C_check_df = standardize(C_df, C_std, 'C', 'check', C_vpdb_ids)
    C_model, C_rel, C_slope, C_int = VPDB_model(C_std_df, cfg, 'C', C_check_df, C_vpdb_ids)
    cal_C_df = apply_vpdb_calibration(C_df, C_model, C_rel, cfg, 'C', C_slope, C_int, C_vpdb_ids)

    #EA-IRMS DataFrame for Nitrogen
    N_df = df[df['Element Type'] == 'Nitrogen'].copy()
    N_vpdb_ids = _detect_vpdb_standards(N_df, N_std)
    N_std_df = standardize(N_df, N_std, 'N', 'std', N_vpdb_ids)
    N_check_df = standardize(N_df, N_std, 'N', 'check', N_vpdb_ids)
    N_model, N_rel, N_slope, N_int = VPDB_model(N_std_df, cfg, 'N', N_check_df, N_vpdb_ids)
    cal_N_df = apply_vpdb_calibration(N_df, N_model, N_rel, cfg, 'N', N_slope, N_int, N_vpdb_ids)

    # Calibrated EA-IRMS DataFrame (Nitrogen and Carbon combined)
    final_df = pd.concat([cal_N_df, cal_C_df], ignore_index=False)
    final_df.sort_values(by="Seconds Since Start", inplace=True)

    return final_df

def standardize(df, standards, tag, option, vpdb_ids):
    """
    Create the VPDB model-standard or RS-accuracy-check dataframe and attach
    true isotope values from the EA standards table.
    """

    if tag == 'N':
        col = 'd15N(AIR) value'
        el = 'Nitrogen'
    elif tag == 'C':
        col = 'd13C(VPDB) value'
        el = 'Carbon'
    else:
        raise ValueError(f"Unknown isotope tag: {tag}")

    if col not in standards.columns:
        raise ValueError(f"The standards dataframe is missing required column: '{col}'")

    if option == "std":
        selected_ids = vpdb_ids["model"]
        label = f"{_scale_name(tag)} calibration"
    elif option == "check":
        selected_ids = vpdb_ids["check"]
        label = "RS accuracy check"
    else:
        raise ValueError("option must be either 'std' or 'check'")

    df = df.copy()
    df["_id_norm"] = df["Identifier 1"].map(_norm_name)

    standards = standards.copy()
    standards["_std_norm"] = standards["EA-IRMS Standards"].map(_norm_name)

    if len(selected_ids) == 0:
        append_to_log(
            log_file_path,
            f"No {label} standards marked True for {el}; skipping."
        )
        return df.iloc[0:0].drop(columns=["_id_norm"])

    found = []
    skipped = []
    valid_ids = []

    for std_id in selected_ids:
        std_name = vpdb_ids["name_lookup"].get(std_id, std_id)
        mask = df["_id_norm"] == std_id

        if not mask.any():
            skipped.append(std_name)
            append_to_log(
                log_file_path,
                f"{std_name} not present in {el}; skipping as {label} standard."
            )
            continue

        std_match = standards[standards["_std_norm"] == std_id]

        if std_match.empty:
            skipped.append(std_name)
            append_to_log(
                log_file_path,
                f"{std_name} not found in {el} standards file; skipping as {label} standard."
            )
            continue

        std_val = pd.to_numeric(std_match[col].iloc[0], errors="coerce")

        if pd.isna(std_val):
            skipped.append(std_name)
            append_to_log(
                log_file_path,
                f"{std_name} has no valid {col}; skipping as {label} standard for {el}."
            )
            continue

        df.loc[mask, col] = std_val
        found.append(std_name)
        valid_ids.append(std_id)

    if len(valid_ids) == 0:
        raise ValueError(
            f"No valid {label} standards were found for {el}. "
            f"Checked standards: {', '.join(selected_ids)}"
        )

    df_out = df[df["_id_norm"].isin(valid_ids)].copy()
    df_out.drop(columns=["_id_norm"], inplace=True)

    append_to_log(
        log_file_path,
        f"{label} standard(s) used for {el}: {', '.join(found)}"
    )

    if skipped:
        append_to_log(
            log_file_path,
            f"{label} standard(s) skipped for {el}: {', '.join(skipped)}"
        )

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
            (e.g., "AIR_corr_d15N/14N", "VPDB_corr_d13C/12C").
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
        return cfg.dN_col[0], cfg.dN_col[1], "AIR_corr_d15N/14N", 'd15N(AIR) value', r"$\delta^{15}\mathrm{N}$", "Nitrogen"
    elif "C" in tag.upper():
        return cfg.dC_col[0], cfg.dC_col[1], "VPDB_corr_d13C/12C", 'd13C(VPDB) value', r"$\delta^{13}\mathrm{C}$", "Carbon"
    else:
        append_to_log(log_file_path, f"Unknown isotope tag: {tag}")
        raise ValueError(f"Unknown isotope tag: {tag}")

# def VPDB_model(df, cfg, tag, bl_df):
def VPDB_model(df, cfg, tag, check_df, vpdb_ids):
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
    input_col, input_description, output_col, actual_col, iso_label, el = get_isotope(tag, cfg)
    scale = _scale_name(tag)
    append_to_log(log_file_path, f"\n{scale} Calibration Model for {el}: ")

    # WLS model
    meas = df[input_col].to_numpy()
    act  = df[actual_col].to_numpy()

    X_full = sm.add_constant(meas.reshape(-1, 1))  # keep a full matrix for later prediction

    # drop infinite
    fit_mask = np.isfinite(meas) & np.isfinite(act)
    X = X_full[fit_mask]
    y = act[fit_mask]
    df_fit = df.loc[fit_mask].copy()

    # Weights - based on standard deviation from replicate median
    g = df_fit['Identifier 1'].astype('string').str.lower()
    x = df_fit[input_col].to_numpy()

    # per-group median and MAD of measured values
    med = df_fit.groupby(g, dropna=False)[input_col].transform('median').to_numpy()
    mad = (df_fit.groupby(g, dropna=False)[input_col]
              .transform(lambda s: np.median(np.abs(s - np.median(s))))).to_numpy() # mean absolute deviation - less senstivie to outlier than standard deviation

    eps = 1e-9 # avoid divide by zero
    z = np.abs(x - med) / (1.4826 * mad + eps) # ((measured_value) - (median_value))/(scaled_mad) # Scaled mad assumes gaussian distribution
    # wDown weight if far from median
    w = 1.0 / (1.0 + z**2)
    w = np.clip(w, 0.05, 1.0) # safety
    w[~np.isfinite(w)] = 0.05

    model = sm.WLS(y, X, weights=w).fit()

    y_pred = model.predict(X_full)
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

    if check_df.empty:
        append_to_log(log_file_path, f"No RS accuracy-check standards available for {scale} model testing.")

        mean = y.mean()
        combined_df = df.copy()
    else:
        check_names = [
            vpdb_ids["name_lookup"].get(x, x).upper()
            for x in vpdb_ids["check"]
        ]

        append_to_log(
            log_file_path,
            f"Testing {scale} Calibration on RS accuracy-check standard(s): {', '.join(check_names)}"
        )

        vals_check = check_df[input_col].to_numpy()
        y_check_actual = check_df[actual_col].to_numpy()
        X_check = sm.add_constant(vals_check.reshape(-1, 1))
        y_check_pred = model.predict(X_check)

        check_df = check_df.copy()
        check_df[output_col] = y_check_pred

        check_bias = np.mean(y_check_pred - y_check_actual)
        check_rmse = np.sqrt(np.mean((y_check_pred - y_check_actual) ** 2))

        mean_pred = f"{np.mean(y_check_pred):.3f}"
        bias = f"{check_bias:.3f}"
        rmse_msg = f"{check_rmse:.3f}"

        print(
            f"RS accuracy check {output_col}",
            f"Mean Prediction: {mean_pred}",
            f"Bias: {bias}",
            f"RMSE: {rmse_msg}",
            sep='\n',
            end='\n\n'
        )

        append_to_log(log_file_path, f"Mean Prediction: {mean_pred}")
        append_to_log(log_file_path, f"Bias: {bias}")
        append_to_log(log_file_path, f"Root Mean Squared Error (RMSE): {rmse_msg}")

        mean = (y.sum() + y_check_actual.sum()) / (y.size + y_check_actual.size)

    plot_ids = vpdb_ids["model"] + vpdb_ids["check"]
    plot_VPDB_calibration(
        # combined_df,
        df,
        check_df,
        plot_ids,
        el,
        actual_col,
        output_col,
        iso_label,
        scale,
        slope,
        intercept,
    )

    return model, mean, slope, intercept

def plot_VPDB_calibration(
    df,
    check_df,
    identifiers,
    el,
    actual_col,
    output_col,
    iso_label,
    scale,
    slope=None,
    intercept=None,
    input_col=None,
):
    """
    Plot actual standard values versus VPDB model-predicted values.

    This plot is only for VPDB calibration/model diagnostics.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing model/calibration standards only.
    check_df : pandas.DataFrame
        DataFrame containing RS accuracy-check standards only.
    identifiers : list of str
        Retained for compatibility, but not used for deciding what to plot.
    el : str
        Element name, e.g. "Nitrogen" or "Carbon".
    actual_col : str
        Column containing true standard values.
    output_col : str
        Column containing VPDB model-predicted values.
    iso_label : str
        Isotope label for plot axes.
    slope, intercept : float, optional
        Calibration model slope/intercept. Included for compatibility.
    input_col : str, optional
        Not used here. Included for compatibility.
    """

    plt.figure(figsize=(6, 4))

    colors = plt.cm.tab10.colors
    color_i = 0

    id = "VPDB_Model"

    # Plot model standards from df
    for standard_id in df["Identifier 1"].map(_norm_name).dropna().unique():
        sub_df = df[df["Identifier 1"].map(_norm_name) == standard_id].copy()

        if sub_df.empty:
            continue

        if actual_col not in sub_df.columns or output_col not in sub_df.columns:
            continue

        x = sub_df[actual_col].to_numpy()
        y = sub_df[output_col].to_numpy()

        label = str(sub_df["Identifier 1"].iloc[0])

        plt.scatter(
            x,
            y,
            label=label,
            marker="o",
            color=colors[color_i % len(colors)],
            alpha=0.6,
            s=100,
            edgecolor="k",
        )
        color_i += 1

    # Plot check standards from check_df
    if check_df is not None and not check_df.empty:
        for check_standard_id in check_df["Identifier 1"].map(_norm_name).dropna().unique():
            sub_df = check_df[
                check_df["Identifier 1"].map(_norm_name) == check_standard_id
            ].copy()

            if sub_df.empty:
                continue

            if actual_col not in sub_df.columns or output_col not in sub_df.columns:
                continue

            x = sub_df[actual_col].to_numpy()
            y = sub_df[output_col].to_numpy()

            check_label = str(sub_df["Identifier 1"].iloc[0])

            plt.scatter(
                x,
                y,
                label=f"Check standard: {check_label}",
                marker="s",
                color=colors[color_i % len(colors)],
                alpha=0.8,
                s=120,
                edgecolor="k",
            )

            color_i += 1

    # 1:1 line
    all_x = []
    all_y = []

    if actual_col in df.columns and output_col in df.columns:
        all_x.extend(df[actual_col].dropna().to_numpy())
        all_y.extend(df[output_col].dropna().to_numpy())

    if check_df is not None and not check_df.empty:
        if actual_col in check_df.columns and output_col in check_df.columns:
            all_x.extend(check_df[actual_col].dropna().to_numpy())
            all_y.extend(check_df[output_col].dropna().to_numpy())

    if len(all_x) > 0 and len(all_y) > 0:
        lower = min(np.min(all_x), np.min(all_y))
        upper = max(np.max(all_x), np.max(all_y))
        pad = 0.05 * (upper - lower) if upper > lower else 1.0

        plt.plot(
            [lower - pad, upper + pad],
            [lower - pad, upper + pad],
            linestyle="--",
            color="black",
            linewidth=1,
            alpha=0.7,
            label="1:1",
        )

    plt.xlabel(f"Actual {iso_label} (‰, {scale})")
    plt.ylabel(f"Predicted {iso_label} (‰, {scale})")
    plt.legend()
    plt.tight_layout()

    fig_name = f"{id}_{el}.png"
    fig_path = os.path.join(fig, fig_name)
    plt.savefig(fig_path)

    plt.show()

def plot_VPDB_calibration_applied(df, identifiers, el, actual_col, output_col, iso_label, input_col, input_description):
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
    colors = plt.cm.tab10.colors
    scale = "AIR" if output_col.startswith("AIR_") else "VPDB"
    plot_id = scale
    standard_ids = [x for x in identifiers if x != "sample"]

    for i, identifier in enumerate(identifiers):
        if identifier == "sample":
            sub_df = df[ ~df["Identifier 1"].map(_norm_name).isin(standard_ids)].copy()
            label = "Samples"
            marker = "x"
            edgecolor = None
        else:
            sub_df = df[df["Identifier 1"].map(_norm_name) == identifier].copy()
            label = str(sub_df["Identifier 1"].iloc[0]) if not sub_df.empty else identifier
            marker = "o"
            edgecolor = "k"
        if sub_df.empty:
            continue
        x = sub_df[input_col].to_numpy()
        y = sub_df[output_col].to_numpy()
        plt.scatter(x, y, label=label, marker=marker, color=colors[i % len(colors)], alpha=0.75, s=100, ec=edgecolor)

    plt.xlabel(f"Measured ({input_description}) {iso_label} (‰)")
    plt.ylabel(f"{scale}-calibrated {iso_label} (‰)")
    plt.legend()
    plt.tight_layout()
    fig_name = f"{plot_id}_{el}.png"
    fig_path = os.path.join(fig, fig_name)
    plt.savefig(fig_path)
    plt.show()

# def apply_vpdb_calibration(df, model, rel, cfg, tag, slope, intercept):
def apply_vpdb_calibration(df, model, rel, cfg, tag, slope, intercept, vpdb_ids):
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
    input_col, input_description, output_col, actual_col, iso_label, el = get_isotope(tag, cfg)
    vals = df[input_col].to_numpy()
    X = sm.add_constant(vals.reshape(-1, 1))
    y = model.get_prediction(X)
    y_pred = y.predicted_mean
    cal_diff = y_pred - rel
    se = y.se_mean
    df[output_col] = y_pred
    df[f"{output_col}_diff"] = cal_diff
    df[f"{output_col}_se"] = se

    # plot_VPDB_calibration_applied(df, identifiers, el, actual_col, output_col, iso_label, input_col)
    plot_ids = vpdb_ids["model"] + vpdb_ids["check"] + ["sample"]
    plot_VPDB_calibration_applied(df, plot_ids, el,
        actual_col,
        output_col,
        iso_label,
        input_col,
        input_description)
    rows = df[output_col].notna().sum()
    msg = f"{_scale_name(tag)} calibration applied to {el}:\n  {tag} rows calibrated: {rows}"
    print(msg)
    append_to_log(log_file_path, msg)

    if tag == 'N':
        print("")

    return df

def _norm_name(x):
     """
     Normalize standard/sample names for matching between the standards file
     and the measurement file.
     """
     if pd.isna(x):
         return ""
     return str(x).strip().lower()


def _coerce_bool_series(s):
     """
     Coerce common CSV boolean representations to True/False.
     """
     if s.dtype == bool:
         return s

     return (
         s.astype(str)
          .str.strip()
          .str.lower()
          .map({
              "true": True,
              "t": True,
              "yes": True,
              "y": True,
              "1": True,
              "1.0": True,
              "false": False,
              "f": False,
              "no": False,
              "n": False,
              "0": False,
              "0.0": False,
              "": False,
              "nan": False,
          })
          .fillna(False)
          .astype(bool)
     )


def _detect_vpdb_standards(df, standards):
     """
     Identify which standards from the standards table are present in the
     measurement dataframe.

     Matching is currently case-insensitive and whitespace-trimmed.

     Returns
     -------
     dict
         {
             "model": list of normalized identifiers used to fit the VPDB model,
             "check": list of normalized identifiers used as RS accuracy checks,
             "all": list of normalized identifiers for all matched standards,
             "name_lookup": dict mapping normalized identifier -> original standard name
         }
     """
     required_cols = {"EA-IRMS Standards"}
     missing = required_cols.difference(standards.columns)
     if missing:
         raise ValueError(
             f"The standards dataframe is missing required column(s): {sorted(missing)}"
         )

     if "Identifier 1" not in df.columns:
         raise ValueError("The measurement dataframe is missing required column: 'Identifier 1'")

     standards = standards.copy()
     standards["_std_norm"] = standards["EA-IRMS Standards"].map(_norm_name)

     measured_ids = set(df["Identifier 1"].map(_norm_name).dropna())
     matched = standards[standards["_std_norm"].isin(measured_ids)].copy()

     if matched.empty:
         raise ValueError(
             "No standards from standards['EA-IRMS Standards'] were found in "
             "df['Identifier 1'] after lowercase matching."
         )

#      if "Use as Standard" in matched.columns:
#          matched["Use as Standard"] = _coerce_bool_series(matched["Use as Standard"])
#          model_rows = matched[matched["Use as Standard"]].copy()
#      else:
#          model_rows = matched.copy()
#
     # if "RS accuracy check" in matched.columns:
#          matched["RS accuracy check"] = _coerce_bool_series(matched["RS accuracy check"])
#          check_rows = matched[matched["RS accuracy check"]].copy()
#      else:
#          check_rows = matched.iloc[0:0].copy()
#
#      if model_rows.empty:
#          raise ValueError(
#              "Standards were found in the measurement data, but none were marked "
#              "'Use as Standard' == True."
#          )
     if "RS accuracy check" in matched.columns:
        matched["RS accuracy check"] = _coerce_bool_series(matched["RS accuracy check"])
     else:
        matched["RS accuracy check"] = False
     check_rows = matched[matched["RS accuracy check"]].copy()

     if "Use as Standard" in matched.columns:
        matched["Use as Standard"] = _coerce_bool_series(matched["Use as Standard"])
     else:
        matched["Use as Standard"] = True
     model_rows = matched[
        matched["Use as Standard"] & ~matched["RS accuracy check"]].copy()

     name_lookup = dict(zip(matched["_std_norm"], matched["EA-IRMS Standards"]))

     return {
         "model": model_rows["_std_norm"].tolist(),
         "check": check_rows["_std_norm"].tolist(),
         "all": matched["_std_norm"].tolist(),
         "name_lookup": name_lookup,
     }
