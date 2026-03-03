import pandas as pd
from . .queries import *
from . .base_functions import append_to_log
# from .utils.base_functions import *
def get_total_analytical_error(df, cfg):
    # base error (injection / raw measurement)
    sigma2 = 0

    if cfg.drift_applied:
        sigma2 += df[f"drift_error"]**2

    if cfg.linearity_applied:
        sigma2 += df[f"linearity_error"]**2

    # reference standard normalization always applied
    sigma2 += df[f"RS_error"]**2

    return np.sqrt(sigma2)

def calculate_methanol_dD(cfg, unknown, isotope, log_file_path):
    """
    Calculation of dD and dC of methanol from PAME analysis
    PAME correction equation from Lee et al., (2017). This equatoin solves
    for the hydrogen (carbon) isotope compostion of the methyl group,
    not the bulk methanol. 
    
    PAME error calucalted using inverse-variance weighting of actual PAME 
    uncertainty and error from corrections.
    """
    pame_unknown =  unknown[unknown['chain']=="PAME"]
    samples =  unknown[unknown['chain']!="PAME"]
    pame_uk = pame_unknown.copy()
    pame_uk.loc[:,f'PAME_methanol_{isotope}'] = pd.NA
    phthalic_original, pame_std = q_original_phthalic_value(isotope)
    if isotope == "dD":
        X = pame_uk.loc[pame_uk['chain']=="PAME", f'RS_{isotope}']
        Y = phthalic_original
        pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}'] = ((10 * X - 4 * Y) / 6) # Methanol value
        sigma_X = get_total_analytical_error(pame_uk, cfg) # Propogated error from corrections
        sigma_Y = pame_std # Known phthalic uncertainty
        # Propagated methanol correction uncertainty
        sigma_methanol = np.sqrt((10/6)**2 * sigma_X**2 + (4/6)**2 * sigma_Y**2)
        pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}_error'] = sigma_methanol
    elif isotope == "dC":
        X = pame_uk.loc[pame_uk['chain']=="PAME", f'RS_{isotope}']
        Y = phthalic_original
        pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}'] = ((10 * X - 8 * Y) / 2) # Methanol value
        sigma_X = get_total_analytical_error(pame_uk, cfg) # Propogated error from corrections
        sigma_Y = pame_std # Known phthalic uncertainty
        # Propagated methanol correction uncertainty
        sigma_methanol = np.sqrt((10/2)**2 * sigma_X**2 + (8/2)**2 * sigma_Y**2)
        pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}_error'] = sigma_methanol
    
    # Calcualte mean and error of calucalted methanol isotope values
    values = pame_uk[f'PAME_methanol_{isotope}'].astype(float)
    errors = pame_uk[f'PAME_methanol_{isotope}_error'].astype(float)
    mask = (~values.isna()) & (~errors.isna()) # remove NaN
    values = values[mask]
    errors = errors[mask]
    weights = 1 / errors**2
    weighted_mean = (weights * values).sum() / weights.sum() # Weighted mean
    sigma_analytical = np.sqrt(1 / weights.sum()) # Analytical uncertainty of mean
    sigma_rep = values.std(ddof=1) # Between-replicate scatter
    sigma_final = np.sqrt(sigma_analytical**2 + sigma_rep**2) # Final combined uncertainty
    print(f" Mean (std) calculated methanol value from PAME: {values.mean().round(2)} ({sigma_final.round(2)}), (n = {len(values)})")
    append_to_log(log_file_path, f" - Mean (std) calculated methanol {isotope} value from PAME: {values.mean().round(2)} ({sigma_final.round(2)}), (n = {len(values)})")
    # pame_uk = pame_uk.dropna(axis=1, how='all')
    return samples, pame_uk, values.mean(), sigma_final