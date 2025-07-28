#src/pysotope/EA/utils/
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error

def drift_correction(df):
    """
    Filter standards with name 'SORGHUM' and plot tiime vs isotope values
    for Nitrogen.

    Parameters:
        df (pd.DataFrame): Full EA dataset with labeled element type.
    
    Currently working:
        The code derives the regression models for the carbon and nitrogen 
        isotopes suing sorghum, displaying the uncorrected and corrected 
        values for only sorghum.
    To do:
        • Upon seeing the corrected sorghum values, the user should select whether 
            to apply to the drift corrections to Nitrogen, Carbon, both, or 
            neither
        • Add standard and correction information to log file (based on iso_process)
        • Add regresion information to log (linear equation, fit values)
        • Save figures to output folder
    """
    # Build model
    N, C = get_sorghum(df)
    N_model, N_relative = drift_model(N, "N")
    C_model, C_relative = drift_model(C, "C")
    
    # Confirm application with user
    
    # Apply correction to samples
    df = apply_drift_model(df, N_model, N_relative, C_model, C_relative)
    return df
def get_sorghum(df):
    sorghum = df[df['Identifier 1'].str.lower() == 'sorghum'].copy()
    if sorghum.empty:
        print("No 'SORGHUM' standards found.")
        return
    sorghum_N = sorghum[sorghum['Element Type'] == 'Nitrogen'] # Only consider nitrogen for sorghum
    sorghum_C = sorghum[sorghum['Element Type'] == 'Carbon']
    return sorghum_N, sorghum_C
    
def drift_model(df, isotope_tag):
    iso_col, iso_label = get_isotope(isotope_tag)
    
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
    print(f"{iso_col}:  y = {slope:.2f}·x + {intercept:.2f}")
    print(f"Adj. R$^{{2}}$ = {adj_r2:.3f}")
    print(f"RMSE = {rmse:.3f}")

    plot_drift_correction(secs, y, y_corr, iso_label)
    return model, y.mean()
    
def get_isotope(name):
    if "C" in name.upper():
        # raw string avoids \d warning
        return "d 13C/12C", r"$\delta^{13}\mathrm{C}$"
    elif "N" in name.upper():
        return "d 15N/14N", r"$\delta^{15}\mathrm{N}$"
    else:
        raise ValueError(f"Unknown isotope tag: {name}")
    
def plot_drift_correction(X, y, y_corr, iso_label):
    plt.figure(figsize=(6, 4))

    plt.scatter(X, y,      c="k",   s=200, alpha=0.5,
                label=f"Uncorrected {iso_label}")
    plt.scatter(X, y_corr, c="red", s=200, alpha=0.5, edgecolor="k",
                label=f"Corrected {iso_label}")

    plt.xlabel("Seconds Since Start")
    plt.ylabel(iso_label)
    plt.legend()

    plt.tight_layout()
    plt.show()
    
def apply_drift_model(df_corr, N_model, N_rel, C_model, C_rel):
    secs = df_corr["Seconds Since Start"].to_numpy().reshape(-1, 1)
    mask_N = (df_corr["Component"] == "N2") & df_corr["d 15N/14N"].notna()
    if mask_N.any():
        X_N       = sm.add_constant(secs[mask_N])
        pred_N    = N_model.get_prediction(X_N)
        drift_N   = pred_N.predicted_mean - N_rel
        se_N      = pred_N.se_mean                               
        df_corr.loc[mask_N, "d 15N/14N_corr"]  = (
            df_corr.loc[mask_N, "d 15N/14N"].to_numpy() - drift_N) # Corrected values
        df_corr.loc[mask_N, "d 15N/14N_se"]    = se_N # Prediction uncertainty
    mask_C = (df_corr["Component"] == "CO2") & df_corr["d 13C/12C"].notna()
    if mask_C.any():
        X_C       = sm.add_constant(secs[mask_C])
        pred_C    = C_model.get_prediction(X_C)
        drift_C   = pred_C.predicted_mean - C_rel
        se_C      = pred_C.se_mean
        df_corr.loc[mask_C, "d 13C/12C_corr"]  = (
            df_corr.loc[mask_C, "d 13C/12C"].to_numpy() - drift_C)
        df_corr.loc[mask_C, "d 13C/12C_se"]    = se_C

    print("Drift correction applied:"
          f"\n  N rows corrected: {mask_N.sum()}"
          f"\n  C rows corrected: {mask_C.sum()}")

    return df_corr