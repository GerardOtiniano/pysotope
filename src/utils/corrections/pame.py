import pandas as pd
from utils.base_functions import *

def calculate_methanol_dD(unknown, isotope, log_file_path):
    """
    Calculation of dD and dC of methanol from PAME analysis
    """
    # Split samples and pames into seperate dataframes
    pame_unknown =  unknown[unknown['chain']=="Phthalic acid"]
    samples =  unknown[unknown['chain']!="Phthalic acid"]
    pame_uk = pame_unknown.copy()
    pame_uk.loc[:,'PAME_methanol_dD'] = pd.NA
    phthalic_original = float(q_original_phthalic_value())
    pame_uk.loc[pame_uk['chain']=="Phthalic acid", 'PAME_methanol_dD'] = (10 * pame_uk.loc[pame_uk['chain']=="Phthalic acid", 'VSMOW_dD'] - 4 * phthalic_original) / 6
    append_to_log(log_file_path, "Calculated methanol values from PAMEs: "+str(pame_uk.loc[pame_uk['chain']=="Phthalic acid", 'PAME_methanol_dD']))
    pame_uk = pame_uk.dropna(axis=1, how='all')
    return samples, pame_uk