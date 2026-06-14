# # -*- coding: utf-8 -*-
import pandas as pd
from . .queries import *
from . .queries import neg_response
from . .figures import *
from . .base_functions import create_subfolder
from IPython.display import clear_output
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from . .curve_fitting import *
from .function_alignment import build_norm_with_fallback

def apply_corr(df, model_name, popt, pcov, lin_reference, dD_id, isotope, used_eiv=False):
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
    out[f'linearity_corrected_{isotope}'] = out[dD_id].to_numpy(float) + offset
    if used_eiv:
        sx = out['sigma_area'].to_numpy(float) if 'sigma_area' in out.columns else None
        out['linearity_error'] = prediction_std_eiv(model_name, x, popt, pcov, sx=sx, nsigma=2)
    else:
        out['linearity_error'] = prediction_std(model_name, x, popt, pcov, nsigma=2)
    return out


def lin_response(log_file_path):
    valid_responses = ['yes', 'y', 'true', 't', 'no', 'n', 'false', 'f']
    while True:
        response = input("\nAssign a linearity correction? (Y/N)\n").lower()
        if response in valid_responses:
            append_to_log(log_file_path, "- Linearity application: " + str(response))
            return response
        else:
            print("\nInvalid response. Try again.\n")

def determine_linearity_model(norm, area_cutoff, dD_id, include_parabolic, force_linearity_model=None):
    """
    Determine the best linearity model using OLS fitting.
    """
    filtered = norm.loc[norm["area"] >= area_cutoff].copy()
    x = filtered["area"].to_numpy(float)
    y = filtered[dD_id].to_numpy(float)
    best_model, popt, sse, pcov = fit_and_select_best(
        x, y, include_parabolic=include_parabolic, force_linearity_model=force_linearity_model)
    return best_model, popt, pcov, sse, x, y

def generate_model_curve(x, best_model, popt):
    if best_model == "linear":
        return linear_func(x, *popt)
    elif best_model == "decay":
        return exp_decay(x, *popt)
    elif best_model == "growth":
        return exp_growth(x, *popt)
    elif best_model == "parabolic":
        return parabolic_func(x, *popt)
    else:
        raise ValueError("Unknown linearity model.")

def process_linearity_correction(
    cfg,
    samp,
    drift,
    lin_std,
    user_choice,
    correction_log,
    folder_path,
    fig_path,
    isotope,
    user_linearity_conditions,
    log_file_path,
    include_parabolic,
    force_linearity_model=None):
    append_to_log(log_file_path, "Linearity correction")
    append_to_log(log_file_path, f"- Linearity applied to {user_choice} isotope values")
    dD_id = user_choice#isotope
    norm, norm_meta = build_norm_with_fallback(
        lin_std,
        y_col=dD_id,
        # y_col=user_choice,
        area_col="area",
        chain_col="chain",
        log_file_path=log_file_path)
    # shift values positive for fitting stability
    # norm[dD_id] = norm[dD_id] - norm[dD_id].min() + 1
    response = lin_response(log_file_path)
    if neg_response(response):
        print("\nSkipping linearity correction.\n")
        return drift, correction_log, lin_std, samp
    if user_linearity_conditions:
        while True:
            val = input(
                "\nEnter peak area cutoff value for C20 and C28 linearity standards\n")
            try:
                area_cutoff = float(val)
                correction_log.loc["Linearity"] = area_cutoff
                break
            except ValueError:
                print("Invalid input. Enter a numerical value.\n")
    else:
        area_cutoff = 0
        correction_log.loc["Linearity"] = area_cutoff
    best_model, popt, pcov, sse, x, y = determine_linearity_model(
        norm,
        area_cutoff,
        dD_id,
        include_parabolic,
        force_linearity_model = force_linearity_model)
    y_fit = generate_model_curve(x, best_model, popt)
    verify_lin_plot(
        norm,
        samp,
        fig_path,
        dD_id,
        log_file_path,
        cutoff_line=area_cutoff,
        isotope=isotope,
        best_model=best_model,
        popt=popt,
        sse=sse)
    # User confirm
    user_input = input("\nDoes this look correct? (Y/N)\n").lower()
    if neg_response(user_input):
        print("\nSkipping linearity correction.\n")
        clear_output(wait=True)
        return drift, correction_log, lin_std, samp
    else:
        cfg.linearity_applied = True

    lin_std, drift, samp, excluded_drift, excluded_lin_std, excluded_samp = (
        linearity_correction(
            drift,
            samp,
            lin_std,
            norm,
            area_cutoff,
            dD_id,
            folder_path,
            fig_path,
            log_file_path,
            isotope,
            best_model,
            popt,
            pcov))
    append_to_log(log_file_path, "- Linearity model applied: " + best_model)
    return drift, correction_log, lin_std, samp

def linearity_correction(
    drift,
    samp,
    lin_std,
    lin_norm,
    area_cutoff,
    dD_id,
    folder_path,
    fig_path,
    log_file_path,
    isotope,
    best_model,
    popt,
    pcov):
    area_cutoff = float(area_cutoff)

    filtered_lin_norm = lin_norm.loc[lin_norm["area"] >= area_cutoff].copy()
    filtered_lin_std = lin_std.loc[lin_std["area"] >= area_cutoff].copy()
    filtered_drift = drift.loc[drift["area"] >= area_cutoff].copy()
    filtered_samp = samp.loc[samp["area"] >= area_cutoff].copy()

    excluded_drift = drift.loc[drift["area"] < area_cutoff].copy()
    excluded_lin_std = lin_std.loc[lin_std["area"] < area_cutoff].copy()
    excluded_samp = samp.loc[samp["area"] < area_cutoff].copy()

    xdata = filtered_lin_norm["area"].to_numpy(float)
    ydata = filtered_lin_norm[dD_id].to_numpy(float)

    y_fit = generate_model_curve(xdata, best_model, popt)

    tss = np.sum((ydata - ydata.mean()) ** 2)
    sse = np.sum((ydata - y_fit) ** 2)

    r_squared = 1.0 if tss == 0 else 1 - (sse / tss)

    append_to_log(log_file_path, f"- Best fit model type: {best_model}")
    append_to_log(log_file_path, f"- Model parameters: {popt}")
    append_to_log(log_file_path, f"- Model stats: R² = {r_squared:.3f} | SSE = {sse:.3f}")

    # Reference value
    lin_top_sort = filtered_lin_norm.sort_values(by="area", ascending=False)
    top_count = max(int(len(lin_top_sort) * 0.2), 1)
    lin_reference = lin_top_sort.head(top_count)[dD_id].median()

    # Corrections
    filtered_lin_std = apply_corr(
        filtered_lin_std,
        best_model,
        popt,
        pcov,
        lin_reference,
        dD_id,
        isotope)
    filtered_drift = apply_corr(
        filtered_drift,
        best_model,
        popt,
        pcov,
        lin_reference,
        dD_id,
        isotope)
    filtered_samp = apply_corr(
        filtered_samp,
        best_model,
        popt,
        pcov,
        lin_reference,
        dD_id,
        isotope)

    # print(list(filtered_lin_std))
    # fig = plt.figure()
    # plt.scatter(filtered_lin_std['area'], filtered_lin_std['linearity_corrected_dD'])
    # plt.show()
    # for i in filtered_lin_std['chain'].unique():
    #     fig = plt.figure()
    #     temp = filtered_lin_std[filtered_lin_std['chain']==i]
    #     plt.scatter(temp['area'], temp['dD'], c='k')
    #     plt.scatter(temp['area'], temp['linearity_corrected_dD'], c='red')

    #     plt.show()

    return (
        filtered_lin_std,
        filtered_drift,
        filtered_samp,
        excluded_drift,
        excluded_lin_std,
        excluded_samp)