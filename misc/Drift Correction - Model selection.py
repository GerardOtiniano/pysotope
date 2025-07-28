#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 10:48:23 2025

@author: gerard
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.linear_model import RidgeCV

# ---------------------------------------------------------------------
# 1.  build design matrix X and response y
# ---------------------------------------------------------------------
predictors = ['Ampl. 44', 'Ampl. 45', 'Ampl. 46', 'Seconds Since Start']#["Ampl. 28", "Ampl. 29"]#, "Area 28",  "Area All","Ampl. 28", "Ampl. 29", "Area 29"]
response   = 'd 13C/12C'#"d 15N/14N"
df = pd.read_csv('/Users/gerard/Documents/GitHub/pysotope/src/pysotope/EA/example_raw.csv')
df = add_seconds_since_start(df)
samples = df.loc[(df['Identifier 2']=="SAMPLE")&(df['Identifier 1']!="SORGHUM")]
df = df[df['Identifier 1']=='SORGHUM']
df = df[df['Identifier 2']=="STANDARD"]
df = df[df['Component']=="CO2"]#"N2"]

df_clean = df.dropna(subset=predictors + [response])      # drop rows with NaNs
X = sm.add_constant(df_clean[predictors])                 # add intercept column
y = df_clean[response]

def vif_table(X):
    """Return a DataFrame with VIF values for each column in X (no constant)."""
    return pd.DataFrame({
        "variable": X.columns,
        "VIF": [variance_inflation_factor(X.values, i)
                for i in range(X.shape[1])]
    }).sort_values("VIF", ascending=False)

# build X without constant for VIF check
X0 = df_clean[predictors]
print(vif_table(X0))

X_sc = df_clean[predictors].values        
y_sc = df_clean[response].values

ridge = RidgeCV(alphas=np.logspace(-3, 3, 50)).fit(X_sc, y_sc)
print("Best α:", ridge.alpha_)
print("R²:", ridge.score(X_sc, y_sc))

# ------------------------------------------------------------------
# 2.  Fit ordinary least-squares model
# ------------------------------------------------------------------
model = sm.OLS(y, X).fit()

# ------------------------------------------------------------------
# 3.  Print the classic regression summary
# ------------------------------------------------------------------
print(model.summary())

# ------------------------------------------------------------------
# 4.  Apply the fitted correction
#     (subtract fitted effect but keep the grand mean)
# ------------------------------------------------------------------
y_hat       = model.fittedvalues
mean_y      = y.mean()
df_clean["d 15N/14N_corrected"] = y - (y_hat - mean_y)

# ------------------------------------------------------------------
# 5.  Visual check: original vs. corrected
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4))

ax.scatter(df_clean["Seconds Since Start"], df_clean[response],
           label="Original", color='k',  s=200, alpha=0.6, edgecolors="k")
ax.scatter(df_clean["Seconds Since Start"], df_clean["d 15N/14N_corrected"],
           label="Corrected", color='red', s=200, alpha=0.6, edgecolors="k")

ax.set_xlabel("Time Since Start (s)")
ax.set_ylabel("δ13C/12C")#"δ¹⁵N/¹⁴N")
ax.legend()
plt.tight_layout()
# plt.savefig(f'/Users/gerard/Documents/GitHub/pysotope/Figures/Sorghum Nitrogen Standard/Carbon Correction Amplitudes and Time.png')
plt.show()

# ------------------------------------------------------------------
# 6.  Predict δ13C/12C for unknown samples
# ------------------------------------------------------------------
# --- make sure 'samples' has the same preprocessing ---------------
samples = add_seconds_since_start(samples)          # already done above
samples_clean = samples.dropna(subset=predictors)   # keep rows w/ all X vars

# ----------—  OLS prediction --------------------------------------
X_new_ols = sm.add_constant(samples_clean[predictors])
samples_clean["d 13C/12C_pred_OLS"] = model.predict(X_new_ols)

# optional: corrected = measured – (fitted – overall mean from *standards*)
overall_mean = y.mean()
samples_clean["d 13C/12C_corr_OLS"] = (
    samples_clean[response]                       # measured value in samples
    - (samples_clean["d 13C/12C_pred_OLS"] - overall_mean)
)

# ----------—  Ridge prediction ------------------------------------
X_new_ridge = samples_clean[predictors].values
samples_clean["d 13C/12C_pred_RIDGE"] = ridge.predict(X_new_ridge)

# ------------------------------------------------------------------
# 7.  Inspect or save
# ------------------------------------------------------------------
print(samples_clean[
    ["Identifier 1", "Identifier 2", "Seconds Since Start",
     response, "d 13C/12C_pred_OLS", "d 13C/12C_pred_RIDGE"]
].head())
# %%

fig = plt.figure()
plt.scatter(samples_clean['d 13C/12C'], samples_clean['d 13C/12C_corr_OLS'])
plt.show()
# %%

def add_seconds_since_start(df):
    """
    Combines 'Date' and 'Time' columns into a datetime object and computes
    seconds since the first timestamp.

    Parameters:
        df (pd.DataFrame): Input dataframe with 'Date' and 'Time' columns.

    Returns:
        pd.DataFrame: Modified dataframe with 'Seconds Since Start' column.
    """
    df = df.copy()

    # Combine 'Date' and 'Time' into a single datetime column
    df['Datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], errors='coerce')

    # Drop rows where datetime couldn't be parsed (optional)
    df = df[df['Datetime'].notna()]

    # Calculate seconds since the first timestamp
    t0 = df['Datetime'].iloc[0]
    df['Seconds Since Start'] = (df['Datetime'] - t0).dt.total_seconds()

    return df

# %%

import statsmodels.api as sm
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# stepwise selection ---------------------------------------------------
# ---------------------------------------------------------------------
def stepwise_selection(df, response, candidates,
                       criterion="AIC", verbose=True,
                       direction="both", max_iter=None):
    """
    Forward-backward stepwise regression based on AIC or BIC.

    Parameters
    ----------
    df : DataFrame
        Data containing response and candidate predictors.
    response : str
        Column name of response variable.
    candidates : list[str]
        List of column names to consider as predictors.
    criterion : {"AIC", "BIC"}, default "AIC"
        Information criterion to minimise.
    verbose : bool
        If True, prints the chosen step each iteration.
    direction : {"forward", "backward", "both"}
        Search strategy.
    max_iter : int or None
        Optional cap on number of steps.

    Returns
    -------
    model : fitted statsmodels OLS object with selected predictors.
    selected : list[str] names of predictors in the final model.
    """
    def fit_sm(cols):
        X = sm.add_constant(df[cols])
        y = df[response]
        model = sm.OLS(y, X).fit()
        return model

    selected = []                    # predictors in current model
    remaining = candidates.copy()
    current_score = np.inf
    best_new_score = np.inf
    step = 0

    while True:
        changed = False
        step += 1
        if max_iter and step > max_iter:
            break

        # --- forward step -------------------------------------------------
        if direction in ("forward", "both"):
            scores_with_candidates = []
            for cand in remaining:
                model = fit_sm(selected + [cand])
                score = model.aic if criterion.upper() == "AIC" else model.bic
                scores_with_candidates.append((score, cand, model))
            if scores_with_candidates:
                best_new_score, best_cand, best_model = min(scores_with_candidates, key=lambda x: x[0])
                if best_new_score < current_score:
                    selected.append(best_cand)
                    remaining.remove(best_cand)
                    current_score = best_new_score
                    changed = True
                    if verbose:
                        print(f"Step {step:2d}: + {best_cand:<20} {criterion}={current_score:.3f}")

        # --- backward step -------------------------------------------------
        if direction in ("backward", "both") and selected:
            scores_with_candidates = []
            for cand in selected:
                cols = selected.copy()
                cols.remove(cand)
                if cols:
                    model = fit_sm(cols)
                    score = model.aic if criterion.upper() == "AIC" else model.bic
                else:
                    score = np.inf
                scores_with_candidates.append((score, cand))
            worst_new_score, worst_cand = min(scores_with_candidates, key=lambda x: x[0])
            if worst_new_score < current_score:
                selected.remove(worst_cand)
                remaining.append(worst_cand)
                current_score = worst_new_score
                changed = True
                if verbose:
                    print(f"Step {step:2d}: - {worst_cand:<20} {criterion}={current_score:.3f}")

        if not changed:
            break

    final_model = fit_sm(selected) if selected else fit_sm([])
    return final_model, selected

# predictors to screen
cands = ['Seconds Since Start', 'Ampl. 28', 'Ampl. 29',
         'Amt%', 'Area 28', 'Area 29', 'Area All']

# run stepwise (both directions, minimise AIC)
model, chosen = stepwise_selection(df_clean, 'd 15N/14N', cands,
                                   criterion="AIC", direction="both",
                                   verbose=True)

print("\nSelected predictors:", chosen)
print(model.summary())

# %%
for x in ['Ampl. 29', 'Ampl. 28']:
    fig = plt.figure()
    plt.scatter(df_clean[x], df_clean['d 15N/14N'], c='k')
    plt.xlabel(x)
    plt.ylabel("δ15N/14N")
    plt.savefig(f'/Users/gerard/Documents/GitHub/pysotope/Figures/Sorghum Nitrogen Standard/{x} vs dN.png')
    plt.show()