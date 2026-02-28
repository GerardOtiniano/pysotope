import pandas as pd

def methyl_correction(unknown, stds, isotope, mdD = -72.5, mdD_err = 3.1):
    """
    Correct FAMES for δD of methyl groups introduced 
    during methylation.
    ~GAO~ 12/4/2023
    """
    # Extract the number of carbons from the 'chain' column
    c_n = unknown.loc[unknown['chain']!="PAME", 'chain'].str.extract(r'C(\d+)').astype(int).squeeze()
    # Apply the correction formula
    unknown.loc[unknown['chain']!="PAME", f'methanol_{isotope}'] = ((unknown[f'RS_{isotope}'] * (2 * c_n + 2)) - (mdD * 3)) / (2 * c_n)
    unknown.loc[unknown['chain']!="PAME",'methanol_error'] = mdD_err
    return unknown
