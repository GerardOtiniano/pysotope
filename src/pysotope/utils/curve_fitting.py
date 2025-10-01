from scipy.optimize import curve_fit
from scipy import odr
import numpy as np
import pandas as pd


def linear_func(x, m, b):
    return m*x + b

def exp_decay(x, a, b, c):
    # a·exp(−b·x) + c  → decays from (a+c) to c as x↑
    return a * np.exp(-b * x) + c

def exp_growth(x, a, b, c):
    # a·(1–exp(−b·x)) + c → grows from c to (a+c) as x↑
    return a * (1 - np.exp(-b * x)) + c

def _odr_linear(beta, x):
    m, b = beta
    return linear_func(x, m, b)

def _odr_decay(beta, x):
    a, b, c = beta
    return exp_decay(x, a, b, c)

def _odr_growth(beta, x):
    a, b, c = beta
    return exp_growth(x, a, b, c)

def guess_linear_params(x, y):
    m, b = np.polyfit(x, y, 1)
    return (m, b)

def guess_decay_params(x, y):
    # assume f(0)=a+c≈y[0], f(∞)=c≈y[-1]
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float) 
    c0 = y[-1]
    a0 = y[0] - c0
    # b0 ≈ 1 / span
    b0 = 1.0 / (x.max() - x.min() + 1e-6)
    return (a0 if a0>0 else 1.0, b0, c0)

def guess_growth_params(x, y):
    # assume f(0)=c≈y[0], f(∞)=a+c≈y[-1]
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float) 
    c0 = y[0]
    a0 = y[-1] - c0
    b0 = 1.0 / (x.max() - x.min() + 1e-6)
    return (a0 if a0>0 else 1.0, b0, c0)


# def fit_and_select_best(x, y):
#     # 1) Linear
#     p0_lin = guess_linear_params(x, y)
#     popt_lin, pcov_lin = curve_fit(linear_func, x, y, p0=p0_lin, maxfev=2_000_000)
#     resid_lin = y - linear_func(x, *popt_lin)
#     sse_lin   = np.sum(resid_lin**2)

#     # 2) Exponential decay
#     p0_dec = guess_decay_params(x, y)
#     popt_dec, pcov_dec = curve_fit(
#         exp_decay, x, y, p0=p0_dec,
#         bounds=([0, 0, -np.inf], [np.inf, np.inf, np.inf]),
#         maxfev=2_000_000
#     )
#     resid_dec = y - exp_decay(x, *popt_dec)
#     sse_dec   = np.sum(resid_dec**2)

#     # 3) Exponential growth
#     p0_gro = guess_growth_params(x, y)
#     popt_gro, pcov_gro = curve_fit(
#         exp_growth, x, y, p0=p0_gro,
#         bounds=([0, 0, -np.inf], [np.inf, np.inf, np.inf]),
#         maxfev=2_000_000
#     )
#     resid_gro = y - exp_growth(x, *popt_gro)
#     sse_gro   = np.sum(resid_gro**2)

#     # Compare SSEs
#     sse_list   = [sse_lin, sse_dec, sse_gro]
#     model_list = ["linear", "decay", "growth"]
#     popt_list  = [popt_lin, popt_dec, popt_gro]
#     pcov_list  = [pcov_lin, pcov_dec, pcov_gro]

#     idx = int(np.argmin(sse_list))
#     return model_list[idx], popt_list[idx], sse_list[idx], pcov_list[idx]

def fit_and_select_best(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    order = np.argsort(x)          # sort by area
    x = x[order]; y = y[order]

    p0_lin = guess_linear_params(x, y)
    popt_lin, pcov_lin = curve_fit(linear_func, x, y, p0=p0_lin, maxfev=2_000_000)
    sse_lin = np.sum((y - linear_func(x, *popt_lin))**2)

    p0_dec = guess_decay_params(x, y)
    popt_dec, pcov_dec = curve_fit(
        exp_decay, x, y, p0=p0_dec,
        bounds=([0, 0, -np.inf], [np.inf, np.inf, np.inf]),
        maxfev=2_000_000)
    sse_dec = np.sum((y - exp_decay(x, *popt_dec))**2)

    p0_gro = guess_growth_params(x, y)
    popt_gro, pcov_gro = curve_fit(
        exp_growth, x, y, p0=p0_gro,
        bounds=([0, 0, -np.inf], [np.inf, np.inf, np.inf]),
        maxfev=2_000_000)
    sse_gro = np.sum((y - exp_growth(x, *popt_gro))**2)

    sse_list   = [sse_lin, sse_dec, sse_gro]
    model_list = ["linear", "decay", "growth"]
    popt_list  = [popt_lin, popt_dec, popt_gro]
    pcov_list  = [pcov_lin, pcov_dec, pcov_gro]
    idx = int(np.argmin(sse_list))
    return model_list[idx], popt_list[idx], sse_list[idx], pcov_list[idx]

def fit_linear_model(x,y):
    p0_lin = guess_linear_params(x, y)
    popt_lin, pcov_lin = curve_fit(linear_func, x, y, p0=p0_lin, maxfev=2000000, method = 'dogbox')
    residuals_lin = y - linear_func(x, *popt_lin)
    sse_lin = np.sum(residuals_lin**2)
    return popt_lin, sse_lin, pcov_lin


def prediction_std(model_name, x, popt, pcov, nsigma=1):
    """
    Return 1‑sigma prediction uncertainty for each x, given the fitted parameters
    and their covariance matrix from ``curve_fit``.

    Parameters
    ----------
    model_name : str  ('linear' | 'decay' | 'growth')
    x          : array‑like
    popt       : fitted parameter vector
    pcov       : parameter‑covariance matrix (k×k)
    """
    x = np.asarray(x, dtype=float)

    if model_name == "linear":
        # y = m x + b
        m, b = popt
        J = np.column_stack([x, np.ones_like(x)])       # ∂y/∂m, ∂y/∂b

    elif model_name == "decay":
        # y = a · exp(−b x) + c
        a, b, c = popt
        e = np.exp(-b * x)
        J = np.column_stack([
            e,                # ∂y/∂a
            -a * x * e,       # ∂y/∂b
            np.ones_like(x)   # ∂y/∂c
        ])

    elif model_name == "growth":
        # y = a · (1 − exp(−b x)) + c
        a, b, c = popt
        e = np.exp(-b * x)
        J = np.column_stack([
            1 - e,            # ∂y/∂a
            +a * x * e,       # ∂y/∂b
            np.ones_like(x)   # ∂y/∂c
        ])

    else:
        raise ValueError(f"Unknown model '{model_name}'")

    # Var(ŷ) = J Σ Jᵀ  → one value per row
    var_pred = np.einsum('ij,jk,ik->i', J, pcov, J)
    sigma = np.sqrt(var_pred)          # 1‑sigma
    return nsigma * sigma

# ODR method for linearity (incorporate uncertainty)
def fit_with_odr(x, y, sx=None, sy=None, model="linear", beta0=None, maxit=2000):
    """
    Errors-in-variables fit using ODR (handles σx and σy).
    x, y : 1D arrays
    sx, sy : 1-sigma arrays (or None)
    model : 'linear' | 'decay' | 'growth'
    returns dict with beta, cov_beta, sse_ortho, etc.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    if sx is None:
        sx = np.ones_like(x)
    if sy is None:
        sy = np.ones_like(y)

    if model == "linear":
        f = odr.Model(_odr_linear)
        if beta0 is None:
            m, b = np.polyfit(x, y, 1)
            beta0 = [m, b]
        k = 2
    elif model == "decay":
        f = odr.Model(_odr_decay)
        if beta0 is None:
            # similar to your guesser
            c0 = y[-1]
            a0 = max(y[0] - c0, 1.0)
            b0 = 1.0 / (x.max() - x.min() + 1e-6)
            beta0 = [a0, b0, c0]
        k = 3
    elif model == "growth":
        f = odr.Model(_odr_growth)
        if beta0 is None:
            c0 = y[0]
            a0 = max(y[-1] - c0, 1.0)
            b0 = 1.0 / (x.max() - x.min() + 1e-6)
            beta0 = [a0, b0, c0]
        k = 3
    else:
        raise ValueError("model must be 'linear'|'decay'|'growth'")

    data = odr.RealData(x, y, sx=sx, sy=sy)
    odr_inst = odr.ODR(data, f, beta0=beta0, maxit=maxit)
    out = odr_inst.run()

    return {
        "beta": out.beta,
        "cov_beta": out.cov_beta,
        "k": k,
        "sse_ortho": out.sum_square,
        "res_var": out.res_var,
        "info": out.info,
        "stopreason": out.stopreason,
        "success": out.info in (1, 2, 3),}

def fit_and_select_best_eiv(x, y, sx=None, sy=None):
    """
    ODR-based model selection. Returns (best_model, popt, pcov, details)
    Uses AIC-like score based on orthogonal SSE.
    """
    fits = {}
    for name in ("linear", "decay", "growth"):
        try:
            res = fit_with_odr(x, y, sx=sx, sy=sy, model=name)
            n = len(x)
            k = res["k"]
            sse = res["sse_ortho"]
            # AIC-like score (constants cancel in comparisons)
            score = n * np.log(sse / max(n, 1) + 1e-12) + 2 * k
            fits[name] = (res, score)
        except Exception:
            continue

    if not fits:
        raise RuntimeError("ODR fits failed for all models.")

    best = min(fits.items(), key=lambda kv: kv[1][1])[0]
    res = fits[best][0]
    return best, np.asarray(res["beta"]), np.asarray(res["cov_beta"]), res

def prediction_std_eiv(model_name, x, popt, pcov, sx=None, nsigma=1):
    """
    Combine parameter uncertainty (J * pcov * J^T) with x-uncertainty via slope^2 * sx^2.
    """
    x = np.asarray(x, float)
    if sx is None:
        sx = np.zeros_like(x)

    # Jacobian wrt parameters (same structure as your prediction_std)
    if model_name == "linear":
        m, b = popt
        J = np.column_stack([x, np.ones_like(x)])
        slope = np.full_like(x, m)
    elif model_name == "decay":
        a, b, c = popt
        e = np.exp(-b * x)
        J = np.column_stack([e, -a * x * e, np.ones_like(x)])
        slope = -a * b * e
    elif model_name == "growth":
        a, b, c = popt
        e = np.exp(-b * x)
        J = np.column_stack([1 - e, +a * x * e, np.ones_like(x)])
        slope = a * b * e
    else:
        raise ValueError(f"Unknown model '{model_name}'")

    var_param = np.einsum("ij,jk,ik->i", J, pcov, J)
    var_x = (slope ** 2) * (np.asarray(sx, float) ** 2)
    sigma = np.sqrt(var_param + var_x)
    return nsigma * sigma