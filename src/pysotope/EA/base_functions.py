
import pandas as pd
from datetime import datetime, timedelta
import time
import os
import numpy as np
from matplotlib.dates import date2num
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

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
        # current_datetime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        current_datetime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
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


def append_to_log(log_file_path, log_message):
    # timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file_path, 'a', encoding='utf-8', errors='replace') as log_file:
        print(f" {log_message}", file=log_file)


def query_file_location():
    while True:
        loc = input("\nProvide the full path of the raw EA-IRMS datafile (as .csv).\n")
        
        # Remove single quotes if present at both the start and end of the input
        if loc.startswith("'") and loc.endswith("'"):
            loc = loc[1:-1]
        
        if os.path.isfile(loc) and loc.endswith(".csv"):
            return loc
        else:
            print("\nFile does not exist or is not a .csv file. Try again.\n")
            
def create_folder():
    # Remove isotope argument
    input_file = query_file_location() # Location of input file
    project_name = "Output "+str(os.path.basename(input_file))
    directory = os.path.dirname(input_file)
    folder_path = os.path.join(directory, project_name)
    log_file_path = create_log_file(folder_path)
    append_to_log(log_file_path, "Compound type: Nitrogen and Carbon") # Place holder
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