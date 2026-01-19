# -*- coding: utf-8 -*-
import pandas as pd
from . .queries import *
from . .queries import neg_response
from . .figures import *
from . .base_functions import create_subfolder
from IPython.display import clear_output
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from . .curve_fitting import *
from .function_alignment import build_norm_with_fallback

def lin_response(log_file_path):
    valid_responses = ['yes', 'y', 'true', 't', 'no', 'n', 'false', 'f']
    while True:
        response = input("\nAssign a linearity correction? (Y/N)\n").lower()
        if response in valid_responses:
            append_to_log(log_file_path, "- Linearity application application: "+str(response))
            return response
        else:
            print("\nInvalid response. Try again.\n")

def process_linearity_correction(cfg, samp, drift, lin_std, user_choice, correction_log, folder_path, fig_path, isotope, 
                                 user_linearity_conditions, log_file_path, include_parabolic):
    append_to_log(log_file_path, "Linearity correction")
    ex = pd.DataFrame()
    dD_id = cfg.dD_col
    norm, norm_meta = build_norm_with_fallback(lin_std, y_col = dD_id, area_col="area", chain_col="chain", log_file_path = log_file_path)
    norm[dD_id] = norm[dD_id]-norm[dD_id].min()+1
    response = lin_response(log_file_path)
    if neg_response(response):
        print("\nSkipping linearity correction.\n")
        return drift, correction_log, lin_std, samp
    else:
        if user_linearity_conditions:
            while True:
                chain_corr_val = input("\nEnter peak area cutoff value for C_{20} and C_{28} linearity stds\n")
                try:
                    chain_corr_vals = [float(num) for num in chain_corr_val.split()]  # Split the input string into a list
                    if len(chain_corr_vals) == 1:
                        correction_log.loc["Linearity"] = chain_corr_vals[0]  # Save the single value
                        break  # Break the loop since the condition is met
                    else:
                        print("Please enter one number.\n")
                except ValueError:
                    print("Invalid input. Please enter a numerical value.\n")
        else:
            chain_corr_val = 0
            correction_log.loc["Linearity"] = chain_corr_val
    chain_corr_val = float(chain_corr_val)
    verify_lin_plot(norm, samp, fig_path, dD_id,log_file_path, cutoff_line=chain_corr_val, isotope=isotope, include_parabolic=include_parabolic)

    user_input = input("\nDoes this look correct? (Y/N)\n").lower()
    if pos_response(user_input):
        cfg.linearity_applied = True
        lin_std, drift, samp, excluded_drift, excluded_lin_std, excluded_samp  = linearity_correction(drift, samp, lin_std, norm, chain_corr_val,
                                                                                                      dD_id, folder_path, fig_path, log_file_path=log_file_path,
                                                                                                      include_parabolic=include_parabolic)
        append_to_log(log_file_path, "- Minimum peak area to derive linearity correction: "+str(chain_corr_val))
        append_to_log(log_file_path, "- Number of drift standards excluded because of below threshold area: "+str(len(excluded_drift)))
        append_to_log(log_file_path, "- Number of linearity standards excluded because of below threshold area: "+str(len(excluded_lin_std)))
        append_to_log(log_file_path, "- Number of samples excluded because of below threshold area: "+str(len(excluded_samp)))
        # Save excluded samples if any
        if not excluded_drift.empty or not excluded_lin_std.empty or not excluded_samp.empty:
            subfolder = create_subfolder(folder_path, "Excluded_Data")
            if not excluded_drift.empty:
                excluded_drift.to_csv(os.path.join(subfolder, 'Drift_standards_excluded_peak_area.csv'), index=False)
            if not excluded_lin_std.empty:
                excluded_lin_std.to_csv(os.path.join(subfolder, 'Linearity_standards_excluded_peak_area.csv'), index=False)
            if not excluded_samp.empty:
                excluded_samp.to_csv(os.path.join(subfolder, 'Samples_excluded_by_peak_area.csv'), index=False)
            #log_excluded_samples(subfolder, excluded_drift, excluded_lin_std, excluded_samp)

        # clear_output(wait=True)
        return drift, correction_log, lin_std, samp

    elif neg_response(user_input):
        print("\nSkipping linearity correction.\n")
        time.sleep(0)  # Wait for 1 second
        clear_output(wait=True)
        return drift, correction_log, lin_std, samp

    else:
        time.sleep(0)
        clear_output(wait=True)
        print("\nInvalid response. Try again.\n")

def linearity_correction(
    drift, samp, lin_std, lin_norm,
    area_cutoff, dD_id, folder_path, fig_path, log_file_path, fig=False, include_parabolic=False):
    area_cutoff = float(area_cutoff)
    filtered_lin_norm = lin_norm.loc[lin_norm['area'] >= area_cutoff].copy()
    filtered_lin_std  = lin_std.loc[lin_std['area']   >= area_cutoff].copy()
    filtered_drift    = drift.loc[drift['area']       >= area_cutoff].copy()
    filtered_samp     = samp.loc[samp['area']         >= area_cutoff].copy()

    excluded_drift    = drift.loc[drift['area']       < area_cutoff].copy()
    excluded_lin_std  = lin_std.loc[lin_std['area']   < area_cutoff].copy()
    excluded_samp     = samp.loc[samp['area']         < area_cutoff].copy()
    xdata = filtered_lin_norm['area'].to_numpy(float)
    ydata = filtered_lin_norm[dD_id].to_numpy(float)

    used_eiv = False
    try:
        # best_model, popt, pcov, odr_details = fit_and_select_best_eiv(xdata, ydata)
        trip = summarize_triplicates(lin_norm, xcol="area", ycol=dD_id, groupcols=("chain","Identifier 1"), use_se=True)
        x = trip["xbar"].to_numpy(float)
        y = trip["ybar"].to_numpy(float)
        sx = trip["sx"].to_numpy(float)   # ← this is your σx from triplicates
        sy = trip["sy"].to_numpy(float)   # ← this is your σy from triplicates

        best_model, popt, pcov, details = fit_and_select_best_eiv(x, y, sx=sx, sy=sy, include_parabolic=include_parabolic)
        if best_model == "linear":
            y_fit = linear_func(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} x + {popt[1]:.6g}"
        elif best_model == "decay":
            y_fit = exp_decay(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} · exp(-{popt[1]:.6g} x) + {popt[2]:.6g}"
        elif best_model == "growth":
            y_fit = exp_growth(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} · (1 - exp(-{popt[1]:.6g} x)) + {popt[2]:.6g}"
        elif best_model == "parabolic":
            y_fit = parabolic(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} x² + {popt[1]:.6g} x + {popt[2]:.6g}"

        sse_like = details['sse_ortho']
        tss = np.sum((ydata - ydata.mean()) ** 2)
        r_squared = 1.0 if tss == 0 else 1 - (np.sum((ydata - y_fit)**2) / tss)
        append_to_log(log_file_path, f"- Best fit model type (EIV/ODR): {best_model}")
        append_to_log(log_file_path, f"- {best_model} model: {parameter_text}")
        append_to_log(log_file_path, f"- {best_model} model stats (ODR): pseudo-R²: {r_squared:.3f} | Ortho SSE: {sse_like:.3f}")
        used_eiv = True

    except Exception:
        best_model, popt, sse, pcov = fit_and_select_best(xdata, ydata, include_parabolic=include_parabolic)
        if best_model == "linear":
            y_fit = linear_func(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} x + {popt[1]:.6g}"
        elif best_model == "decay":
            y_fit = exp_decay(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} · exp(-{popt[1]:.6g} x) + {popt[2]:.6g}"
        elif best_model == "growth":  # growth
            y_fit = exp_growth(xdata, *popt)
            parameter_text = f"y = {popt[0]:.6g} · (1 - exp(-{popt[1]:.6g} x)) + {popt[2]:.6g}"
        elif best_model == "parabolic":  # growth
            y_fit = parabolic_func(xdata, *popt)
            parameter_text =  f"y = {popt[0]:.6g} · x^2 + {popt[1]:.6g} · x + {popt[2]:.6g}"

        tss = np.sum((ydata - ydata.mean()) ** 2)
        r_squared = 1.0 if tss == 0 else 1 - (sse / tss)

        append_to_log(log_file_path, f"- Best fit model type (OLS): {best_model}")
        append_to_log(log_file_path, f"- {best_model} model: {parameter_text}")
        append_to_log(log_file_path, f"- {best_model} model stats (OLS): R²: {r_squared:.3f} | SSE: {sse:.3f}")

    lin_top_sort   = filtered_lin_norm.sort_values(by='area', ascending=False)
    top_count      = max(int(len(lin_top_sort) * 0.2), 1)
    lin_top_qt     = lin_top_sort.head(top_count)
    lin_reference  = lin_top_qt[dD_id].median()  # robust anchor

    # Optional: quick asymptote check for decay/growth
    if best_model in ("decay", "growth"):
        a, b, c = popt
        asymptote = (a + c) if best_model == "growth" else c
        append_to_log(log_file_path, f"- Asymptote check: {asymptote:.3f} vs ref {lin_reference:.3f}")

    filtered_lin_std = apply_corr(filtered_lin_std, best_model, popt, pcov, lin_reference, dD_id,
                                  used_eiv=used_eiv)
    filtered_drift   = apply_corr(filtered_drift,   best_model, popt, pcov, lin_reference, dD_id,
                                  used_eiv=used_eiv)
    filtered_samp    = apply_corr(filtered_samp,    best_model, popt, pcov, lin_reference, dD_id,
                                  used_eiv=used_eiv)

    return filtered_lin_std, filtered_drift, filtered_samp, excluded_drift, excluded_lin_std, excluded_samp


def apply_corr(df, model_name, popt, pcov, lin_reference, dD_id, used_eiv=False):
    """
    If used_eiv=True, compute prediction uncertainty with x-error contribution
    (requires optional column 'sigma_area' for per-row sx; defaults to 0).
    """
    x = df['area'].to_numpy(float)

    # predicted model value at each x
    if model_name == "linear":
        yhat = linear_func(x, *popt)
    elif model_name == "decay":
        yhat = exp_decay(x, *popt)
    elif model_name == "growth":
        yhat = exp_growth(x, *popt)
    elif model_name == "parabolic":
        yhat = parabolic_func(x, *popt)

    offset = lin_reference - yhat

    out = df.copy()
    out['linearity_corrected_dD'] = out[dD_id].to_numpy(float) + offset

    if used_eiv:
        sx = out['sigma_area'].to_numpy(float) if 'sigma_area' in out.columns else None
        out['linearity_error'] = prediction_std_eiv(model_name, x, popt, pcov, sx=sx, nsigma=2)
    else:
        out['linearity_error'] = prediction_std(model_name, x, popt, pcov, nsigma=2)

    return out

def summarize_triplicates(df, xcol="area", ycol="dD", groupcols=("chain","Identifier 1"), use_se=True,
                          rel_floor=0.002, abs_floor_x=1e-9, abs_floor_y=1e-9):
    """
    Collapse triplicates (or n-licates) to means and uncertainties.
    If use_se=True: returns SE for sx, sy (good for means). Otherwise returns SD.
    Applies small floors so weights don't explode.
    """
    g = df.groupby(list(groupcols), dropna=False)
    out = g.agg(
        xbar=(xcol, "mean"),
        ybar=(ycol, "mean"),
        sx  =(xcol, "std"),
        sy  =(ycol, "std"),
        n   =(ycol, "size"),
    ).reset_index()

    # replace NaN std for n<2 with 0
    out["sx"] = out["sx"].fillna(0.0)
    out["sy"] = out["sy"].fillna(0.0)

    if use_se:
        out["sx"] = out.apply(lambda r: (r["sx"]/np.sqrt(r["n"])) if r["n"]>=2 else 0.0, axis=1)
        out["sy"] = out.apply(lambda r: (r["sy"]/np.sqrt(r["n"])) if r["n"]>=2 else 0.0, axis=1)

    # apply small relative/absolute floors
    out["sx"] = np.maximum(out["sx"], np.maximum(rel_floor*np.abs(out["xbar"]), abs_floor_x))
    out["sy"] = np.maximum(out["sy"], np.maximum(rel_floor*np.abs(out["ybar"]), abs_floor_y))

    return out