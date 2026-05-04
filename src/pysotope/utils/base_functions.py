import os
from importlib import import_module
from pathlib import Path
from datetime import datetime as dt
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import append_to_log
from .queries import *
from .queries import query_file_location

_gc_importer = import_module('..import.gc_importer', __package__)


def make_correction_df():
    return _gc_importer.make_correction_df()


def chain_subsetrer(std_df, std_meta, std_type):
    return _gc_importer.chain_subsetrer(std_df, std_meta, std_type)


def import_data(data_location, folder_path, log_file_path, isotope, standards_df):
    return _gc_importer.import_data(data_location, folder_path, log_file_path, isotope, standards_df)


def closest_rt(df, time_val, target_rt, threshold=0.05):
    return _gc_importer.closest_rt(df, time_val, target_rt, threshold=threshold)


def process_dataframe(df, rt_dict, folder_path, log_file_path):
    return _gc_importer.process_dataframe(df, rt_dict, folder_path, log_file_path)


def create_log_file(folder_path):
    """Create log file."""
    import platform
    import matplotlib
    import scipy
    import statsmodels
    import sklearn
    import IPython

    os.makedirs(folder_path, exist_ok=True)
    log_file_path = os.path.join(folder_path, 'Log file.txt')
    with open(log_file_path, 'w') as log_file:
        current_datetime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        initial_message = 'Log file created at ' + str(current_datetime) + '\n'
        log_file.write(initial_message)
        log_file.write(f'Python version: {platform.python_version()}\n')
        log_file.write(f'pandas version: {pd.__version__}\n')
        log_file.write(f'matplotlib version: {matplotlib.__version__}\n')
        log_file.write(f'numpy version: {np.__version__}\n')
        log_file.write(f'scipy version: {scipy.__version__}\n')
        log_file.write(f'statsmodels version: {statsmodels.__version__}\n')
        log_file.write(f'sklearn (scikit-learn) version: {sklearn.__version__}\n')
        log_file.write(f'IPython version: {IPython.__version__}\n\n\n')
    return log_file_path


def create_folder(isotope):
    input_file = query_file_location()
    project_name = 'Output ' + str(os.path.basename(input_file))
    directory = os.path.dirname(input_file)
    folder_path = os.path.join(directory, project_name)
    log_file_path = create_log_file(folder_path)
    iso_name = 'dD' if isotope == 'dD' else 'dC'
    append_to_log(log_file_path, 'Isotope type: ' + str(iso_name))
    os.makedirs(folder_path, exist_ok=True)

    fig_path = os.path.join(folder_path, 'Figures')
    os.makedirs(fig_path, exist_ok=True)

    results_path = os.path.join(folder_path, 'Results')
    os.makedirs(results_path, exist_ok=True)
    return folder_path, fig_path, results_path, input_file, log_file_path


def create_subfolder(folder_path, name):
    subf_path = os.path.join(folder_path, name)
    os.makedirs(subf_path, exist_ok=True)
    return subf_path


def load_standards(isotope: str = 'dD') -> pd.DataFrame:
    from ..standards_manager.editor import get_standard_path, standard_editor

    path = get_standard_path(f'RS_{isotope}.csv')
    if not path.exists():
        print("Didn't find standards")
        return standard_editor()

    df = pd.read_csv(path, dtype={'type': str, 'chain length': str})
    df = df[df['Use as Standard'] == True]
    df = df[df['Use as Standard'] != False]
    df['RS accuracy check'] = df['RS accuracy check'].astype(str).str.lower() == 'true'
    df['Use as Standard'] = df['Use as Standard'].astype(str).str.lower() == 'true'
    return df
