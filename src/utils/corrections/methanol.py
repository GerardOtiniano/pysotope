import pandas as pd

def methyl_correction(unknown, stds, mdD = -72.5, mdD_err = 3.1):
    """
    Correct FAMES for δD of methyl groups introduced 
    during methylation.
    ~GAO~ 12/4/2023
    """
    # Extract the number of carbons from the 'chain' column
    #c_n = unknown.loc[unknown['chain']!="Phthalic acid", 'chain'].str.extract('C(\d+)').astype(int).squeeze()
    c_n = unknown.loc[unknown['chain']!="Phthalic acid", 'chain'].str.extract(r'C(\d+)').astype(int).squeeze()
    # Apply the correction formula
    unknown.loc[unknown['chain']!="Phthalic acid", 'methanol_dD'] = ((unknown['VSMOW_dD'] * (2 * c_n + 2)) - (mdD * 3)) / (2 * c_n)
    unknown.loc[unknown['chain']!="Phthalic acid",'methanol_error'] = mdD_err
    return unknown
