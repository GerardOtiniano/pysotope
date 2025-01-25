from scipy.optimize import curve_fit
import numpy as np
import pandas as pd

def log_func(x, a, b, c):
    """
    a * exp(-b * x) + c
    a: amplitude (scales the height of the curve)
    b: decay rate
    c: offset (vertical shift)
    """
    return a * np.exp(b * x) + c

def exp_func(x, a, b, c):
    """
    Returns a * exp(b * x) + c
    """
    return a * np.exp(b * x) + c

def guess_exponential_params(x, y):
    # x and y are 1D numpy arrays
    x = np.array(x)
    y = np.array(y)
    c0 = y.min()  # offset guess
    y_adjusted = np.array(y - c0)
    # To avoid log of zero or negative, ensure the first and last points are above c0
    # If they are not, you may need a more robust method
    if (y_adjusted <= 0).any():
        # Fallback: set offset to 0 if min(y) is negative or zero
        c0 = 0
        y_adjusted = np.array(y)
    a0 = y_adjusted[0]  # amplitude guess
    x = list(x)
    if len(x) > 1 and (y_adjusted[0] > 0) and (y_adjusted[-1] > 0):
        b0 = np.log(y_adjusted[-1] / y_adjusted[0]) / (x[-1] - x[0])
    else:
        b0 = 0.0  # fallback if we have too few points or invalid ratio
    
    return (a0, b0, c0)

def guess_log_params(x, y):
    x = np.array(x)
    y = np.array(y)
    # x must be > 0 for log
    b0 = 1.0
    # Handle case where x[0], x[-1] might be the same or cause log(0):
    if len(x) > 1 and x[0] > 0 and x[-1] > 0 and x[0] != x[-1]:
        a0 = (y[-1] - y[0]) / (np.log(x[-1]) - np.log(x[0]))
        c0 = y[0] - a0 * np.log(b0 * x[0])
    else:
        a0 = 1.0
        c0 = 0.0
    
    return (a0, b0, c0)


def fit_and_select_best(x, y):
    # Exponential fit
    p0_exp = guess_exponential_params(x, y)
    popt_exp, pcov_exp = curve_fit(exp_func, x, y, p0=p0_exp, maxfev=20000)
    # Evaluate SSE for exponential
    residuals_exp = y - exp_func(x, *popt_exp)
    sse_exp = np.sum(residuals_exp**2)

    # Log fit
    p0_log = guess_log_params(x, y)
    popt_log, pcov_log = curve_fit(log_func, x, y, p0=p0_log, maxfev=20000)
    # Evaluate SSE for log
    residuals_log = y - log_func(x, *popt_log)
    sse_log = np.sum(residuals_log**2)

    # Compare SSE (lower is better)
    if sse_exp < sse_log:
        return "exponential", popt_exp, sse_exp, pcov_exp
    else:
        return "log", popt_log, sse_log, pcov_log

def prediction_std(best_model, x, popt, pcov):
    """
    Compute the standard error (1-sigma) of the predictions from exp_func
    given parameters popt and covariance pcov.
    """
    if best_model == "exponential":
        a, b, c = popt
    
        # Jacobian partial derivatives at each x:
        d_da = np.exp(b*x)       # ∂f/∂a
        d_db = a * x * np.exp(b*x)  # ∂f/∂b
        d_dc = np.ones_like(x)   # ∂f/∂c
        
        # Stack into shape (n_points, 3)
        # Each row: [d_da, d_db, d_dc]
        J = np.vstack([d_da, d_db, d_dc]).T
        
        # For each x_i, Var(f_i) = J_i * pcov * J_i^T
        # We'll compute row by row:
        var_pred = []
        for row in J:
            # row shape is (3,)
            # reshape to (1, 3) for matrix multiplication
            jac_row = row.reshape(1, -1)
            var_i = jac_row @ pcov @ jac_row.T  # scalar
            var_pred.append(var_i[0, 0])
        
        return np.sqrt(var_pred)  # 1-sigma errors
    else:
    
        """
        Computes the 1-sigma standard error of predictions for the model:
            f(x) = a * ln(b*x) + c
    
        Parameters
        ----------
        x : array-like
            The x-values at which to compute prediction errors (must be > 0).
        popt : array-like
            Best-fit parameters [a, b, c] from curve_fit.
        pcov : 2D array
            Covariance matrix of the parameters from curve_fit.
    
        Returns
        -------
        errors : np.ndarray
            1D array of predicted 1-sigma standard errors for each x.
    
        Notes
        -----
        The partial derivatives are:
            ∂f/∂a = ln(b*x)
            ∂f/∂b = a / b
            ∂f/∂c = 1
        so the Jacobian row for each x_i is [ln(b*x_i),  a/b,  1].
        """
        a, b, c = popt
        x = np.array(x, dtype=float)
    
        # Jacobian partial derivatives for each x
        d_da = np.log(b * x)         # ∂f/∂a = ln(b*x)
        d_db = np.full_like(x, a/b)  # ∂f/∂b = a / b (constant w.r.t x)
        d_dc = np.ones_like(x)       # ∂f/∂c = 1
    
        # Construct Jacobian: shape (n_points, 3)
        J = np.vstack([d_da, d_db, d_dc]).T
    
        var_pred = []
        for row in J:
            # row is shape (3,) => reshape to (1, 3) for matrix multiplication
            jac_row = row.reshape(1, -1)
            var_i = jac_row @ pcov @ jac_row.T  # scalar (1x1 matrix)
            var_pred.append(var_i[0, 0])
    
        return np.sqrt(var_pred)  # 1-sigma standard error