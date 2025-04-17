
import pandas as pd
from datetime import datetime, timedelta
import time
import os
import numpy as np
from matplotlib.dates import date2num
from .queries import *
from .queries import query_file_location, query_stds
import matplotlib.pyplot as plt

def make_correction_df():
    correction_log_data = {
        "type": ["Drift", "Linearity", "VSMOW", "Methylation"],
        "sample": [0, 0, 0, 0],  # Default values
    }
    correction_log = pd.DataFrame(correction_log_data)
    correction_log = correction_log.set_index(['type'])
    return correction_log

def try_parse_date(date_str):
    # List of date formats to try
    formats = ["%m/%d/%Y %H:%M:%S", "%m/%d/%y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d/%m/%Y %H:%M:%S"]

    for fmt in formats:
        try:
            return datetime.datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None  # Return None if all formats fail

def create_log_file(folder_path):
    """
    Create log file.
    """
    import platform
    import pandas as pd
    import matplotlib
    import numpy as np
    import scipy
    import statsmodels
    import sklearn
    import IPython
    # Ensure the folder exists
    os.makedirs(folder_path, exist_ok=True)
    # Create the full path for the log file
    log_file_path = os.path.join(folder_path, 'Log file.txt')
    # Create the log file and write the initial message
    with open(log_file_path, 'w') as log_file:
        current_datetime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        initial_message = "Log file created at "+str(current_datetime)+"\n"
        log_file.write(initial_message)
        log_file.write(f"Python version: {platform.python_version()}\n")
        log_file.write(f"pandas version: {pd.__version__}\n")
        log_file.write(f"matplotlib version: {matplotlib.__version__}\n")
        log_file.write(f"numpy version: {np.__version__}\n")
        log_file.write(f"scipy version: {scipy.__version__}\n")
        log_file.write(f"statsmodels version: {statsmodels.__version__}\n")
        log_file.write(f"sklearn (scikit-learn) version: {sklearn.__version__}\n")
        log_file.write(f"IPython version: {IPython.__version__}\n\n\n")  
    return log_file_path

# def append_to_log(log_file_path, log_message):
#     """
#     Add entry to log file.
#     """
#     with open(log_file_path, 'a') as log_file:
#         current_datetime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
#         initial_message = f"\nLog file created at {current_datetime}\n"
#         log_file.write("\n"+str(log_message))# + "; "+str(current_datetime)+"\n")
        
def append_to_log(log_file_path, log_message):
    # timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file_path, 'a') as log_file:
        print(f" {log_message}", file=log_file)
        
def chain_subsetrer(std_df, std_meta, std_type):
    chains = list(std_meta[std_meta['type']==std_type]['chain length']) 
    df = std_df[std_df['Identifier 1'].str.contains(chains[0]) & std_df['Identifier 1'].str.contains(chains[1])]
    df = df[df.chain.isin(chains)]
    return df, chains
        
def import_data(data_location, folder_path, log_file_path, isotope, standards_df):
    """
    Import .csv file from GCIRMs - default .csv file from GCIRMS creates issues with header. The function assigns new header names,
    creates a date-time format for linear regression, identifieds standards, and isolates standards and samples.
    Outputs: 
        df             - original dataframe
        linearirty_std - dataframe with linearity standards
        drif_std       - dataframe with drift standards
        unknown        - dataframe with sample data
        pame           - boolian identifying presence of pames in dataset
    ~GAO~ 12/4/2023
    """
    # Create log file
    
    df = pd.read_csv(data_location)
    new_name = [str(isotope),'area','chain']; x = 0; 
    if isotope == "dD": iso_rat = "d 2H/1H"
    elif isotope == "dC": iso_rat = "d 13C/12C"
    else: raise ValueError("Unsupported isotope system.")
    
    column_found = False  # Flag to check if expected column is found and renamed
    for name in [str(iso_rat),'Area All','Component']:
        if name in df.columns:
            df = df.rename(columns={df.columns[df.columns.str.contains(name)][0]: new_name[x]})
            x += 1
            column_found = True
        else:
            df[new_name[x]] = np.nan
            x += 1

    if not column_found:
        raise ValueError("The expected header for the isotope system was not found in the csv file. Please verify the isotope system of interest.")
    
    #df['date-time_true'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], format='%m/%d/%y %H:%M:%S')
    df['date-time_true'] = df.apply(lambda row: try_parse_date(row['Date'] + ' ' + row['Time']), axis=1)
    df['date-time'] = date2num(df['date-time_true'])
    # shift time relative to maximum (or minimum)
    #df['time_rel']=df['date-time']-df['date-time'].max()
    df['time_rel']=df['date-time']-df['date-time'].min()+1 # add 1 to avoid zero value for logarithmic regression correction

    # Check if PAME is in dataset
    # if df['chain'].astype(str).str.contains("phthalic", case=False, na=False).any():
    if df['chain'].astype(str).str.contains("phthalic", case=False, na=False).any():
        pame = True
        append_to_log(log_file_path, 'PAME detected in analysis')
        print("PAMEs detected.\nThe calculated methanol value from the PAMEs will be displayed to the user and stored in the log file.\n")
    else: pame = False
    
    # Seperate samples, H3+, drift, and linearity standards
    linearity_std, linearity_chain_lengths = chain_subsetrer(df, standards_df, "linearity")
    append_to_log(log_file_path, f"Number of linearity standards analyzed: {len(linearity_std[linearity_std.chain == linearity_chain_lengths[1]])}")
    drift_std, drift_chain_lengths = chain_subsetrer(df, standards_df, "drift")
    append_to_log(log_file_path, f"Number of Drift standards analyzed: {len(drift_std[drift_std.chain == drift_chain_lengths[1]])}")

    # Remove first two drift runs
    drift_std = drift_std.sort_values('date-time_true')
    unique_time_signatures = drift_std["date-time"].unique() # identify unique drift runs
    time_signatures_to_remove = unique_time_signatures[:2] # Modified Jan 7, 2024 - line above is original method, but didnt work?
    drift_std = drift_std[~drift_std["date-time"].isin(time_signatures_to_remove)] # Remove first two runs - OSIBL ignores for variance
    append_to_log(log_file_path, "First two drift standards ignored.")
    
    mask    = (df['Identifier 1'].str.contains(linearity_chain_lengths[0]) & df['Identifier 1'].str.contains(linearity_chain_lengths[1]))
    unknown = df[~mask]
    mask    = (unknown['Identifier 1'].str.contains(drift_chain_lengths[0]) & unknown['Identifier 1'].str.contains(drift_chain_lengths[1]))
    unknown = unknown[~mask]
    unknown = unknown[~unknown['Identifier 1'].str.contains('H3+')]
    rt_dict = ask_user_for_rt(log_file_path)
    if rt_dict:
        unknown   = process_dataframe(unknown, rt_dict, folder_path, log_file_path)
        unknown   = unknown[unknown.chain!="None"]
        linearity_std   = process_dataframe(linearity_std, rt_dict, folder_path, log_file_path)
        drift_std = process_dataframe(drift_std, rt_dict, folder_path, log_file_path)
    else: unknown = unknown[unknown.chain.isin(['Phthalic acid','C16',"C18","C20","C22","C24","C24","C26","C28","C30","C32"])]
    for i in [unknown, drift_std, linearity_std]:
        i = i[~i.chain.isna()]
    correction_log = make_correction_df()
    return linearity_std, drift_std, unknown, correction_log, pame

def create_folder(isotope):
    input_file = query_file_location() # Location of input file
    project_name = "Output "+str(os.path.basename(input_file))
    directory = os.path.dirname(input_file)
    folder_path = os.path.join(directory, project_name)
    log_file_path = create_log_file(folder_path)
    if isotope == "dD": iso_name = "δD"
    else: iso_name = "δC"
    append_to_log(log_file_path, "Isotope type: "+str(iso_name))
    os.makedirs(folder_path, exist_ok=True)
    
    # Make output folders
    fig_path = os.path.join(folder_path, 'Figures')
    os.makedirs(fig_path, exist_ok=True)
    
    results_path = os.path.join(folder_path, 'Results')
    os.makedirs(results_path, exist_ok=True)
    return folder_path, fig_path, results_path, input_file, log_file_path

def create_subfolder(folder_path, name):
    subf_path = os.path.join(folder_path, name)
    os.makedirs(subf_path, exist_ok=True)
    return subf_path


def closest_rt(df, time_val, target_rt, threshold=0.05):
    """
    Find the closest retention time(s) to the target.
    If two values are almost equally close (within a threshold), return both.
    """
    sample_df = df[df['Time'] == time_val]
    differences = (sample_df['Rt'] - target_rt).abs()
    min_diff = differences.min()
    closest_rows = sample_df[differences <= min_diff * (1 + threshold)]
    return closest_rows

def ask_user_for_rt(log_file_path):
    chain_lengths = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']
    while True:
        response = input("Do you want to detect components in this dataset by retention time? (Y/N):\n").strip().lower()
        if pos_response(response):
            append_to_log(log_file_path, "User opted to identify chains.")
            rt_values = input("Enter retention times for " + ", ".join(chain_lengths) + " separated by commas (type 'none' for any you don't want to use):\n")
            rt_values = rt_values.split(',')
            if len(rt_values) == len(chain_lengths):
                rt_dict = {chain: (None if rt.strip().lower() == 'none' else float(rt.strip())) for chain, rt in zip(chain_lengths, rt_values)}
                return rt_dict
            else:
                print("Invalid input. Please provide the correct number of values.\n")
        elif neg_response(response):
            append_to_log(log_file_path, "User opted not to identify chains.")
            print("Component detection not selected.\n")
            return None
        else:
            print("Invalid response. Please answer 'yes' or 'no'.\n")
            
def process_dataframe(df, rt_dict, folder_path, log_file_path):
    if rt_dict is None:
        return df
    rt_path = create_subfolder(folder_path, 'Retention time figures')
    df['chain'] = None
    unique_times = df['Time'].unique()
    for time_val in unique_times:
        sample_id = df.loc[df['Time'] == time_val, 'Identifier 1'].iloc[0]

        # For 'standard' types, use chains mentioned in Identifier 1 or all chains if none are mentioned
        filtered_rt_dict = rt_dict
        for chain, rt in filtered_rt_dict.items():
            if rt is not None:
                closest_rows = closest_rt(df, time_val, rt)
                if len(closest_rows) == 1:
                    # Only one clear closest match
                    correct_rt = closest_rows.iloc[0]['Rt']
                    df.loc[(df['Time'] == time_val) & (df['Rt'] == correct_rt), 'chain'] = chain
                elif len(closest_rows) > 1:
                    # Two closely matched peaks, prompt the user
                    # clear_output(wait=True)
                    plt.figure()
                    sample_df = df[df['Time'] == time_val]
                    plt.scatter(sample_df['Rt'], sample_df['Area All'], label=sample_id, color='red', ec='k')
                    plt.plot(sample_df['Rt'], sample_df['Area All'], label=sample_id, linestyle='--', c='k')
                    x=0; lim = -999
                    for index, (_, row) in enumerate(closest_rows.iterrows(), start=1):
                        if x == 0: lim=row['Rt']
                        plt.axvline(x=row['Rt'], color='red', linestyle='--', alpha = 0.5)
                        plt.text(row['Rt'], sample_df['Area All'].mean()+x, str(index), color='k', fontsize=12, verticalalignment='bottom')
                        x=x+5
                    plt.xlabel('Retention Time')
                    plt.ylabel('Area')
                    plt.title(f'Close Matches for {sample_id} ({time_val}) - {chain}')
                    if lim != -999:
                        if lim > row['Rt']:
                            x_min = row['Rt']-50
                            x_max = lim+50
                        else:
                            x_min = lim-50
                            x_max = row['Rt']+50
                    else:
                        x_min = 450
                        x_max = row['Rt']+50
                    plt.xlim(x_min, x_max)
                    plt.legend()
                    plt.savefig(os.path.join(rt_path, 'Sample '+str(sample_id)+'Chain '+str(chain)+' rt '+str(rt)+'.png'), dpi=300, bbox_inches='tight')
                    plt.show()

                    choice = input(f"Enter the number associated with the correct retention time for {chain} in sample {sample_id} ({time_val}), or type 'none' to skip:\n").strip().lower()
                    if choice == 'none':
                        continue
                    choice = int(choice)
                    correct_rt = closest_rows.iloc[choice - 1]['Rt']
                    df.loc[(df['Time'] == time_val) & (df['Rt'] == correct_rt), 'chain'] = chain
                else:
                    df.loc[closest_rows.index, 'chain'] = chain
    export_df = df[df['chain'].isin(['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32'])]
    append_to_log(log_file_path, f"Chain lengths identified by user: {export_df.chain.unique()}")
    return export_df