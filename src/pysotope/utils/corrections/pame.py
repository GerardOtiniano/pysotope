import pandas as pd
from . .queries import *
from . .base_functions import append_to_log
# from .utils.base_functions import *

def calculate_methanol_dD(unknown, isotope, log_file_path):
    """
    Calculation of dD and dC of methanol from PAME analysis
    """
    # Split samples and pames into seperate dataframes
    pame_unknown =  unknown[unknown['chain']=="PAME"]
    samples =  unknown[unknown['chain']!="PAME"]
    pame_uk = pame_unknown.copy()
    pame_uk.loc[:,f'PAME_methanol_{isotope}'] = pd.NA
    phthalic_original = float(q_original_phthalic_value())
    pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}'] = (10 * pame_uk.loc[pame_uk['chain']=="PAME", f'RS_{isotope}'] - 4 * phthalic_original) / 6
    append_to_log(log_file_path, "Calculated methanol values from PAMEs: "+str(pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}']))
    pames = pame_uk.loc[pame_uk['chain']=="PAME", f'PAME_methanol_{isotope}']
    print(f" Mean (std) calcualted methanol value from PAME: {pames.mean().round(2)} ({pames.std().round(2)})")
    pame_uk = pame_uk.dropna(axis=1, how='all')
    return samples, pame_uk