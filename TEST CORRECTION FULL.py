import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm


# def add_seconds_since_start(df):
#     """
#     Combines 'Date' and 'Time' columns into a datetime object and computes
#     seconds since the first timestamp.

#     Parameters:
#         df (pd.DataFrame): Input dataframe with 'Date' and 'Time' columns.

#     Returns:
#         pd.DataFrame: Modified dataframe with 'Seconds Since Start' column.
#     """
#     df = df.copy()

#     # Combine 'Date' and 'Time' into a single datetime column
#     df['Datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], errors='coerce')

#     # Drop rows where datetime couldn't be parsed (optional)
#     df = df[df['Datetime'].notna()]

#     # Calculate seconds since the first timestamp
#     t0 = df['Datetime'].iloc[0]
#     df['Seconds Since Start'] = (df['Datetime'] - t0).dt.total_seconds()

#     return df
# df = pd.read_csv('/Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv')
# # ---------------------------------------------------------------
# # 1.  predictor list must match the calibration model
# # ---------------------------------------------------------------
# predictors = ["Seconds Since Start", "Ampl. 44", "Ampl. 45", "Ampl. 46"]

# # ---------------------------------------------------------------
# # 2.  design matrix for ALL rows (standards + samples)
# # ---------------------------------------------------------------
# df = add_seconds_since_start(df)             # your helper
# df = df[df["Component"] == "CO2"]            # keep δ13C rows
# X_all = sm.add_constant(df[predictors])


# # ---------------------------------------------------------------
# # 3.  predicted error for every row
# # ---------------------------------------------------------------
# df["pred_err"] = wls_model.predict(X_all)
# print(wls_model.summary())
# # corrected / predicted true value
# df["δ13C_pred_true"] = df["d 13C/12C"] - df["pred_err"]

# # ---------------------------------------------------------------
# # 4.  build a comparison table
# # ---------------------------------------------------------------
# true_vals = {
#     "SORGHUM":      -13.780,
#     "UREA":         -43.260,
#     "ACETANILIDE":  -28.842,
#     "WHEAT FLOUR": -26.871
# }
# df["δ13C_actual"] = df["Identifier 1"].map(true_vals)
# # BL SEDIMENT target:
# BL_TARGET   = -28.012
# BL_SIGMA    =  0.084

# # standards dataframe
# stds_plot = df[df["Identifier 1"].isin(true_vals.keys())].copy()
# # BL SEDIMENT dataframe
# bl_plot   = df[df["Identifier 1"]=="BL SEDIMENT"].copy()
# bl_plot["δ13C_actual"] = BL_TARGET   # same target for all replicates

# # ---------------------------------------------------------------
# # 5.  scatter-plot
# # ---------------------------------------------------------------
# fig, ax = plt.subplots(figsize=(5,5))

# # 1-to-1 reference line
# lims = [min(df["δ13C_actual"].min(), df["δ13C_pred_true"].min()) - 1,
#         max(df["δ13C_actual"].max(), df["δ13C_pred_true"].max()) + 1]
# ax.plot(lims, lims, ls="--", c="grey")
# ax.set_xlim(lims); ax.set_ylim(lims)

# # standards (black ×)
# for name in stds_plot['Identifier 1'].unique():
#     temp = stds_plot[stds_plot['Identifier 1']==name]
#     ax.scatter(temp["δ13C_actual"],
#                temp["δ13C_pred_true"],
#                marker="o", ec='k', s=150, alpha = 0.5, label=name)

# # BL SEDIMENT (red •)
# ax.scatter(bl_plot["δ13C_actual"],
#            bl_plot["δ13C_pred_true"],
#            marker="o", color="red", ec='k', s=200, alpha = 0.5, label="TEST: BL SEDIMENT")

# # ±1 σ band for BL target
# ax.errorbar(BL_TARGET, bl_plot["δ13C_pred_true"].mean(),
#             xerr=BL_SIGMA, fmt="none", ecolor="red", capsize=4)

# ax.set_xlabel("Actual δ$^{13}$C (‰, VPDB)")
# ax.set_ylabel("Predicted δ$^{13}$C (‰, VPDB)")
# ax.legend()
# ax.set_xlim(-30,-11)
# ax.set_ylim(-30,-11)
# plt.tight_layout()
# plt.show()

# # ---------------------------------------------------------------
# # 6.  Numeric diagnostics
# # ---------------------------------------------------------------
# BL_TARGET = -28.012
# BL_SIGMA  =  0.084        # 1-σ certified uncertainty

# # ----------  A.  BL-SEDIMENT performance -----------------------
# bl_resid = bl_plot["δ13C_pred_true"] - BL_TARGET
# n_bl     = len(bl_resid)

# mean_pred_bl = bl_plot["δ13C_pred_true"].mean()
# bias_bl      = mean_pred_bl - BL_TARGET               # systematic offset
# rmse_bl      = np.sqrt((bl_resid**2).mean())          # root-mean-square error
# sd_bl        = bl_resid.std(ddof=1) if n_bl > 1 else 0
# chi2_red_bl  = ((bl_resid / BL_SIGMA)**2).sum() / (n_bl - 1) if n_bl > 1 else np.nan

# print(f"\nBL SEDIMENT (n = {n_bl})")
# print(f"  Mean predicted δ13C  : {mean_pred_bl:+8.3f} ‰")
# print(f"  Bias (mean-target)   : {bias_bl:+8.3f} ‰")
# print(f"  RMSE                 : {rmse_bl:8.3f} ‰")
# print(f"  1-σ scatter          : {sd_bl:8.3f} ‰")
# print(f"  Reduced χ² vs σ=0.084: {chi2_red_bl:8.2f}")

# # ----------  B.  Overall standard-set performance --------------
# std_resid   = stds_plot["δ13C_pred_true"] - stds_plot["δ13C_actual"]
# rmse_stand  = np.sqrt((std_resid**2).mean())
# mean_abs_std= std_resid.abs().mean()

# print("\nCalibration standards (cross-check)")
# print(f"  RMSE (all CRMs)      : {rmse_stand:8.3f} ‰")
# print(f"  Mean |residual|      : {mean_abs_std:8.3f} ‰")
import pandas as pd, numpy as np, statsmodels.api as sm
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# helper to add Seconds-Since-Start
# ------------------------------------------------------------------
def add_seconds_since_start(df):
    df = df.copy()
    df["Datetime"] = pd.to_datetime(df["Date"] + " " + df["Time"], errors="coerce")
    df = df[df["Datetime"].notna()]
    t0 = df["Datetime"].iloc[0]
    df["Seconds Since Start"] = (df["Datetime"] - t0).dt.total_seconds()
    return df

# ------------------------------------------------------------------
# 0.  load & restrict to δ13C rows
# ------------------------------------------------------------------
df = pd.read_csv("/Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv")
df = add_seconds_since_start(df)
df = df[df["Component"] == "CO2"]           # δ13C only

# ------------------------------------------------------------------
# 1.  calibration standards and their certified data
# ------------------------------------------------------------------
true_vals  = {  # ‰ VPDB
    "SORGHUM"     : -13.780,
    "UREA"        : -43.260,
    "ACETANILIDE" : -28.842,
    "WHEAT FLOUR" : -26.871,
}
sigma_vals = {  # 1-σ (‰)
    "SORGHUM"     : 0.085,
    "UREA"        : 0.020,
    "ACETANILIDE" : 0.015,
    "WHEAT FLOUR" : 0.050,
}

crm = df[df["Identifier 1"].isin(true_vals)].copy()
crm["true_13C"]  = crm["Identifier 1"].map(true_vals)
crm["sigma_13C"] = crm["Identifier 1"].map(sigma_vals)
crm["delta_err"] = crm["d 13C/12C"] - crm["true_13C"]   # instrument error (meas−true)

# ------------------------------------------------------------------
# 2.  fit weighted least-squares model   error ~ Seconds + Ampl44 + Ampl45 + Ampl46
# ------------------------------------------------------------------
predictors = ["Seconds Since Start", "Ampl. 44"]
X_crm = sm.add_constant(crm[predictors])
w_crm = 1.0 / (crm["sigma_13C"] ** 2)                    # weights = 1/σ²

wls_model = sm.WLS(crm["delta_err"], X_crm, weights=w_crm).fit()
print("\nWeighted calibration model\n")
print(wls_model.summary())

# ------------------------------------------------------------------
# 3.  apply correction to every row (standards + samples)
# ------------------------------------------------------------------
X_all = sm.add_constant(df[predictors])
df["pred_err"]       = wls_model.predict(X_all)
df["δ13C_pred_true"] = df["d 13C/12C"] - df["pred_err"]   # corrected (on VPDB scale)
df["δ13C_actual"]    = df["Identifier 1"].map(true_vals)  # NaN for unknowns

# ------------------------------------------------------------------
# 4.  set up BL SEDIMENT test
# ------------------------------------------------------------------
BL_TARGET = -28.012
BL_SIGMA  =  0.084

bl_plot   = df[df["Identifier 1"] == "BL SEDIMENT"].copy()
bl_plot["δ13C_actual"] = BL_TARGET

stds_plot = df[df["Identifier 1"].isin(true_vals)].copy()  # CRMs for plotting

# ------------------------------------------------------------------
# 5.  1-to-1 plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5,5))

lims = [min(df["δ13C_actual"].min(), df["δ13C_pred_true"].min()) - 1,
        max(df["δ13C_actual"].max(), df["δ13C_pred_true"].max()) + 1]
ax.plot(lims, lims, "--", color="grey", zorder=1)
ax.set_xlim(lims); ax.set_ylim(lims)

# CRMs (black ◌)
for name, g in stds_plot.groupby("Identifier 1"):
    ax.scatter(g["δ13C_actual"], g["δ13C_pred_true"],
               marker="o", ec="k", fc="none", s=120, alpha=0.6, label=name)

# BL SEDIMENT (red ●)
ax.scatter(bl_plot["δ13C_actual"], bl_plot["δ13C_pred_true"],
           marker="o", color="red", ec="k", s=140, alpha=0.6, label="BL SEDIMENT")

# ±1 σ band for BL target
ax.errorbar(BL_TARGET, bl_plot["δ13C_pred_true"].mean(),
            xerr=BL_SIGMA, fmt="none", ecolor="red", capsize=4)

ax.set_xlabel("Actual δ$^{13}$C (‰ VPDB)")
ax.set_ylabel("Predicted δ$^{13}$C (‰ VPDB)")
ax.legend()
plt.tight_layout(); plt.show()

# ------------------------------------------------------------------
# 6.  numeric diagnostics
# ------------------------------------------------------------------
bl_resid = bl_plot["δ13C_pred_true"] - BL_TARGET
print("\nBL SEDIMENT diagnostics")
print(f"  n                    : {len(bl_resid)}")
print(f"  mean bias            : {bl_resid.mean():+7.3f} ‰")
print(f"  RMSE                 : {np.sqrt((bl_resid**2).mean()):7.3f} ‰")
print(f"  1-σ scatter          : {bl_resid.std(ddof=1):7.3f} ‰")

std_resid = stds_plot["δ13C_pred_true"] - stds_plot["δ13C_actual"]
print("\nCRM cross-check")
print(f"  RMSE (all CRMs)      : {np.sqrt((std_resid**2).mean()):7.3f} ‰")
print(f"  mean |residual|      : {std_resid.abs().mean():7.3f} ‰")

# %%

import pandas as pd, numpy as np, statsmodels.api as sm
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 0.  housekeeping -------------------------------------------------
# ------------------------------------------------------------------
df = pd.read_csv("/Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv")
df = add_seconds_since_start(df)
df = df[df["Component"] == "CO2"]     # work only with δ13C rows

# ------------------------------------------------------------------
# 1.  drift model from SPHAGNUM standards --------------------------
# ------------------------------------------------------------------
sph = df[df["Identifier 1"] == "SORGHUM"].copy()

t0          = sph["Seconds Since Start"].min()           # earliest run
base_value  = sph.loc[sph["Seconds Since Start"] == t0, "d 13C/12C"].iloc[0]

sph["t_rel"]   = sph["Seconds Since Start"] - t0         # drift clock starts at 0
sph["delta13C"] = sph["d 13C/12C"] - base_value          # baseline set to 0 ‰

drift_mod = sm.OLS(sph["delta13C"], sm.add_constant(sph["t_rel"])).fit()
print("\nDrift model:\n", drift_mod.summary())

# predicted drift offset for every row
df["t_rel"]   = df["Seconds Since Start"] - t0
df["drift_off"]= drift_mod.predict(sm.add_constant(df["t_rel"]))
df["d13C_corr"]= df["d 13C/12C"] - df["drift_off"]       # drift-corrected

# ------------------------------------------------------------------
# 2.  scale calibration on three CRMs ------------------------------
# ------------------------------------------------------------------
true_vals = {
    "SORGHUM":     -13.780,
    "ACETANILIDE": -28.842,
    "WHEAT FLOUR": -26.871
}
df["true_13C"] = df["Identifier 1"].map(true_vals) 
crm = df[df["Identifier 1"].isin(true_vals.keys())].copy()
crm["true_13C"] = crm["Identifier 1"].map(true_vals)

# linear calibration: true = a + b·(drift-corrected)
cal_mod = sm.OLS(crm["true_13C"], sm.add_constant(crm["d13C_corr"])).fit()
print("\nScale-transfer model:\n", cal_mod.summary())

# apply calibration to every row
df["d13C_pred_true"] = cal_mod.predict(sm.add_constant(df["d13C_corr"]))
crm["d13C_pred_true"] = df.loc[crm.index, "d13C_pred_true"]
# ------------------------------------------------------------------
# 3.  evaluate on BL SEDIMENT --------------------------------------
# ------------------------------------------------------------------
BL_TARGET = -28.012;  BL_SIGMA = 0.084
bl = df[df["Identifier 1"] == "BL SEDIMENT"].copy()
bl["true_13C"] = BL_TARGET
bl_err = bl["d13C_pred_true"] - BL_TARGET

print(f"\nBL SEDIMENT n={len(bl)}  bias={bl_err.mean():+.3f} ‰  "
      f"RMSE={np.sqrt((bl_err**2).mean()):.3f} ‰  σ={bl_err.std(ddof=1):.3f} ‰")

# ------------------------------------------------------------------
# 4.  visual comparison --------------------------------------------
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5,5.5))

# 1:1 line
lims = [df[["true_13C","d13C_pred_true"]].min().min()-1,
        df[["true_13C","d13C_pred_true"]].max().max()+1]
ax.plot(lims, lims, '--', color='grey', zorder=1)
ax.set_xlim(lims); ax.set_ylim(lims)

# plot CRMs
for name, g in crm.groupby("Identifier 1"):
    ax.scatter(g["true_13C"], g["d13C_pred_true"],
               marker="x", s=100, label=name, zorder=2, edgecolor='k')

# plot BL SEDIMENT
ax.scatter(bl["true_13C"], bl["d13C_pred_true"],
           marker="o", s=120, color='red', edgecolor='k',
           label="BL SEDIMENT", zorder=3)
ax.errorbar(BL_TARGET, bl["d13C_pred_true"].mean(),
            xerr=BL_SIGMA, fmt="none", ecolor="red", capsize=4)
ax.set_xlim(-30,-10)
ax.set_ylim(-30,-10)
ax.set_xlabel("Actual δ$^{13}$C (‰ VPDB)")
ax.set_ylabel("Predicted δ$^{13}$C (‰ VPDB)")
ax.legend()
plt.tight_layout(); plt.show()