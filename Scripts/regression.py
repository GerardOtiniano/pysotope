from base_functions import append_to_log
from scipy.stats import linregress
import numpy as np
from sklearn.linear_model import HuberRegressor
import scipy.stats as stats


def linear_regression(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return slope, intercept, r_value, p_value, std_err


def robust_linear_regression(x, y):
    x_reshaped = np.array(x).reshape(-1, 1)
    huber = HuberRegressor()
    huber.fit(x_reshaped, y)

    slope = huber.coef_[0]
    intercept = huber.intercept_
    r_squared = huber.score(x_reshaped, y)

    # Estimating the standard error of the regression coefficients
    # Predictions
    y_pred = huber.predict(x_reshaped)

    # Residuals
    residuals = y - y_pred

    # Standard deviation of the errors
    sigma = np.sqrt(np.sum(residuals**2) / (len(y) - 2))

    # Standard error of coefficients
    s_xx = np.sum((x - np.mean(x)) ** 2)
    std_err = sigma / np.sqrt(s_xx)

    # T-statistic for the slope
    t_statistic = slope / std_err

    # P-value for the slope
    p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), df=len(x) - 2))

    return slope, intercept, r_squared, p_value, std_err


# def gls_regression(x, y, log_file_path, omega_inv=None):
#     """
#     Generalized Least Squares Regression Model Function.

#     Parameters:
#     - x: Independent variable(s), should be a 1D or 2D array-like structure.
#     - y: Dependent variable, should be a 1D array-like structure.
#     - omega_inv: Inverse of the covariance matrix of the errors. If None, OLS is performed.

#     Returns:
#     - slope: Slope coefficient(s) of the model.
#     - intercept: Intercept of the model.
#     - r_squared: Coefficient of determination.
#     - p_value: p-value associated with the slope(s).
#     - std_err: Standard error of the slope coefficient(s).
#     """
#     append_to_log(log_file_path, "Applied generalized least squares regression function")
#     # Add constant to the independent variable(s) for the intercept
#     x_with_const = sm.add_constant(x)

#     # If no omega_inv is provided, use OLS as a special case of GLS
#     if omega_inv is None:
#         model = sm.OLS(y, x_with_const)
#     else:
#         model = sm.GLS(y, x_with_const, sigma=omega_inv)

#     results = model.fit()

#     # Extracting model parameters
#     intercept, slope = results.params
#     r_squared = results.rsquared
#     p_value = results.pvalues.iloc[1]  # Assuming the slope is the second parameter
#     std_err = results.bse.iloc[1]  # Standard error for the slope


#     return slope, intercept, r_squared, p_value, std_err, results


def gls_regression(x, y, log_file_path, weights=None):  # should be wls_regression - changing it here so that wls can be used without renaming ever instance
    """
    Weighted Least Squares Regression Model Function.

    Parameters:
    - x: Independent variable(s), should be a 1D or 2D array-like structure.
    - y: Dependent variable, should be a 1D array-like structure.
    - weights: An array of weights, should be the same length as y. If None, performs OLS.

    Returns:
    - slope: Slope coefficient(s) of the model.
    - intercept: Intercept of the model.
    - r_squared: Coefficient of determination.
    - p_value: p-value associated with the slope(s).
    - std_err: Standard error of the slope coefficient(s).
    - results: The results object from the regression model.
    """

    append_to_log(log_file_path, "Applied weighted least squares regression function")
    x_with_const = sm.add_constant(x)
    weights = 1 / (np.abs(x) + 1)  # Adding 1 to avoid division by zero, assuming x is your independent variable array
    model = sm.WLS(y, x_with_const, weights=weights)
    results = model.fit()

    intercept, slope = results.params
    r_squared = results.rsquared
    p_value = results.pvalues.iloc[1]  # Assuming the slope is the second parameter
    std_err = results.bse.iloc[1]  # Standard error for the slope

    return slope, intercept, r_squared, p_value, std_err, results

    # df = pd.DataFrame([x,y])
    # df['residuals'] = model.resid
    # df['squared_residuals'] =  df['residuals']**2
    # df['x_group']  = pd.cut(df[x],bins=5)
    # grouped_variance = df.groupby('x_group')['squared_residuals'].mean()
    # df['variance_estimate'] = df['x_group'].map(grouped_variance)
