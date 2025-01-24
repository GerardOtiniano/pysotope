import pandas as pd
import time
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit
from IPython.display import clear_output
from scipy.stats import zscore
from . .queries import *
from . .base_functions import *

# Raw Standards
def remove_standards_PA(l,d):
    print("Identifying standards with isotope values outside the 95% confidence interval.")
    time.sleep(0)
    ci95_remove = pd.DataFrame()
    for chains in ["C18","C20","C24","C28"]:
        temp = pd.concat([l,d])
        mask = temp.chain==chains
        temp.loc[mask,'zscore'] = zscore(temp.loc[mask,'dD'])
        temp = temp.loc[mask]
        plt.scatter(temp[np.abs(temp.zscore)<2].area,temp[np.abs(temp.zscore)<2].zscore, c = 'grey', ec='k')
        ci95_remove = pd.concat([ci95_remove,temp[np.abs(temp.zscore)>2]])
        plt.scatter(temp[np.abs(temp.zscore)>2].area,temp[np.abs(temp.zscore)>2].zscore, c = 'red', ec='k',
                   label = "Standards outside 95% CI")
    
    plt.axhline(2,c='k',alpha=0.2,linestyle='--', label = "95% Confidence Interval")
    plt.xlabel("Peak Area")
    plt.ylabel("Z score")
    plt.axhline(-2,c='k',alpha=0.2,linestyle='--')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=1)
    plt.show()

    print(ci95_remove)
    zscore_stds = pd.concat([l,d])
    zscore_stds['zs'] = zscore(zscore_stds.dD)
    mask = np.abs(zscore_stds.zs)<2
    print(zscore_stds.loc[mask, ["Identifier 1", "chain", "dD", 'area']])
    if pos_response(input("Remove standards outside the 95% confidence interval of measured isotope values? (Y/N)\n")):
        stds  = pd.concat([l,d])
        stds = stds[~mask]
        return stds
    else:
        return pd.concat([l,d])
    
# Samples

def outlier_detect(unknown, log_file_path):
    while True:
        try:
            unique_outlier = float(input("Input the standard deviation acceptable above which samples will be considered highly variable (best to use 0):\n"))
            append_to_log(log_file_path, f"Z-score limit used to detect sample outliers: {unique_outlier}")
            break
        except ValueError:
            print("Please enter a valid number.")
    # Calculate Z-scores of the standard deviations
    unknown['Std Dev Z-Score'] = zscore(unknown['VSMOW_dD'])

    # Identify samples with high standard deviation
    high_std_dev = unknown[np.abs(unknown['Std Dev Z-Score']) > np.float64(unique_outlier)]
    cond1 = high_std_dev['Identifier 1'] + high_std_dev['chain']
    cond2 = unknown["Identifier 1"]+unknown["chain"]
    high_dD = unknown[cond2.isin(cond1)]
    return high_dD

def outlier_removal(unknown, fig_path, log_file_path):
    #unknown = total_uncertainty(unknown)
    excluded = pd.DataFrame()
    id_outliers = input("Check for outliers? (Y/N):\n")
    unknown_final = unknown.copy() # moved outside forloop ~GAO~ 26/01/24
    if pos_response(id_outliers):
        outlier_path = create_subfolder(fig_path, 'Outliers')
        high_unknown = outlier_detect(unknown, log_file_path).copy()
        unique_groups = [(x, y) for x in high_unknown["Identifier 1"].unique() for y in high_unknown[high_unknown["Identifier 1"] == x]["chain"].unique()]
        for x, y in unique_groups:
            temp2 = high_unknown[(high_unknown["Identifier 1"] == x) & (high_unknown["chain"] == y)]
            print(temp2[["Identifier 1", "chain", "VSMOW_dD", 'area']]) # print the replicate analyses of interest
            temp = unknown[unknown.chain == y]
            plt.scatter(temp.area, temp.VSMOW_dD, marker='x', c='k', label=str(y) + " samples", alpha = 0.5)
            plt.scatter(temp2.area, temp2.VSMOW_dD, c='red', ec='k', label=str(x) + ", " + str(y), alpha = 0.75)
            plt.legend()
            plt.savefig(os.path.join(outlier_path, str(x)+' '+str(y)+'.png'), dpi=300, bbox_inches='tight')
            plt.show()
            del_vals = input("Input comma seperated list of index values for row(s) that should be removed. Type 'None' if no samples should be removed:\n")

            if del_vals.lower() != 'none' or del_vals =="":
                try:
                    del_indices = [int(val.strip()) for val in del_vals.split(',')]
                    if not all(idx in unknown.index for idx in del_indices):
                        raise ValueError("One or more index values are invalid. Please enter valid index values.")

                    rows_to_exclude = unknown.loc[del_indices]
                    append_to_log(log_file_path, "Outlier sample removed: "+str(rows_to_exclude["Identifier 1"]+" "+str(rows_to_exclude["chain"])))
                    excluded = pd.concat([excluded, rows_to_exclude])
                    unknown_final = unknown_final.drop(del_indices, errors='ignore')
                    high_unknown = high_unknown.drop(del_indices, errors='ignore')
                except ValueError as ve:
                    print(ve)
                except Exception as e:
                    print(f"An error occurred: {e}")
            clear_output(wait=True)  # Clear the output
        print("All outliers have been processed.")
        return unknown_final, excluded

    else:
        print("Outlier check skipped.")
        return unknown_final, pd.DataFrame()

    time.sleep(0)  # Wait for 1 second