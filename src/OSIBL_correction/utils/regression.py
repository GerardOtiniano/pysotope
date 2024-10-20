from utils.base_functions import append_to_log
import statsmodels.api as sm
import numpy as np

def wls_regression(x, y, log_file_path, weights=None): #should be wls_regression - changing it here so that wls can be used without renaming ever instance
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