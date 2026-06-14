import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import date2num

from ..utils.chains.chains import get_selected_chains
from ..utils.common import append_to_log, id_mask, try_parse_date
from ..utils.queries import ask_user_for_rt


GC_CHAIN_EXPORT = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']


def make_correction_df() -> pd.DataFrame:
    correction_log_data = {
        'type': ['Drift', 'Linearity', 'VSMOW', 'Methylation'],
        'sample': [0, 0, 0, 0],
    }
    correction_log = pd.DataFrame(correction_log_data)
    return correction_log.set_index(['type'])


def chain_subsetrer(std_df: pd.DataFrame, std_meta: pd.DataFrame, std_type: str):
    chains = list(std_meta[std_meta['type'] == std_type]['chain length'])
    ids = list(std_meta[std_meta['type'] == std_type]['ID'])

    mask1 = id_mask(std_df, ids, col='Identifier 1', mode='any')
    df = std_df[mask1]

    if 'Identifier 2' in df.columns and df['Identifier 2'].notna().any():
        mask2 = df['Identifier 2'].astype(str).str.lower() == 'standard'
        df = df[mask2]

    if len(chains) > 0:
        df = df[df.chain.isin(chains)]

    return df, chains, ids


def closest_rt(df: pd.DataFrame, time_val, target_rt, threshold=0.05) -> pd.DataFrame:
    sample_df = df[df['Time'] == time_val]
    differences = (sample_df['Rt'] - target_rt).abs()
    min_diff = differences.min()
    return sample_df[differences <= min_diff * (1 + threshold)]


def _create_rt_subfolder(folder_path):
    rt_path = os.path.join(folder_path, 'Retention time figures')
    os.makedirs(rt_path, exist_ok=True)
    return rt_path


def process_dataframe(df: pd.DataFrame, rt_dict, folder_path, log_file_path) -> pd.DataFrame:
    if rt_dict is None:
        return df

    rt_path = _create_rt_subfolder(folder_path)
    out = df.copy()
    out['chain'] = None
    unique_times = out['Time'].unique()

    for time_val in unique_times:
        sample_id = out.loc[out['Time'] == time_val, 'Identifier 1'].iloc[0]

        for chain, rt in rt_dict.items():
            if rt is None:
                continue

            closest_rows = closest_rt(out, time_val, rt)
            if len(closest_rows) == 1:
                correct_rt = closest_rows.iloc[0]['Rt']
                out.loc[(out['Time'] == time_val) & (out['Rt'] == correct_rt), 'chain'] = chain
            elif len(closest_rows) > 1:
                plt.figure()
                sample_df = out[out['Time'] == time_val]
                plt.scatter(sample_df['Rt'], sample_df['Area All'], label=sample_id, color='red', ec='k')
                plt.plot(sample_df['Rt'], sample_df['Area All'], label=sample_id, linestyle='--', c='k')
                x_offset = 0
                limit = -999
                last_row = None
                for index, (_, row) in enumerate(closest_rows.iterrows(), start=1):
                    last_row = row
                    if x_offset == 0:
                        limit = row['Rt']
                    plt.axvline(x=row['Rt'], color='red', linestyle='--', alpha=0.5)
                    plt.text(
                        row['Rt'],
                        sample_df['Area All'].mean() + x_offset,
                        str(index),
                        color='k',
                        fontsize=12,
                        verticalalignment='bottom',
                    )
                    x_offset += 5
                plt.xlabel('Retention Time')
                plt.ylabel('Area')
                plt.title(f'Close Matches for {sample_id} ({time_val}) - {chain}')
                if limit != -999 and last_row is not None:
                    if limit > last_row['Rt']:
                        x_min = last_row['Rt'] - 50
                        x_max = limit + 50
                    else:
                        x_min = limit - 50
                        x_max = last_row['Rt'] + 50
                else:
                    x_min = 450
                    x_max = last_row['Rt'] + 50 if last_row is not None else 500
                plt.xlim(x_min, x_max)
                plt.legend()
                plt.savefig(
                    os.path.join(rt_path, f'Sample {sample_id}Chain {chain} rt {rt}.png'),
                    dpi=300,
                    bbox_inches='tight',
                )
                plt.show()

                choice = input(
                    f"Enter the number associated with the correct retention time for {chain} "
                    f"in sample {sample_id} ({time_val}), or type 'none' to skip:\n"
                ).strip().lower()
                if choice == 'none':
                    continue
                choice = int(choice)
                correct_rt = closest_rows.iloc[choice - 1]['Rt']
                out.loc[(out['Time'] == time_val) & (out['Rt'] == correct_rt), 'chain'] = chain
            else:
                out.loc[closest_rows.index, 'chain'] = chain

    export_df = out[out['chain'].isin(GC_CHAIN_EXPORT)]
    append_to_log(log_file_path, f"Chain lengths identified by user: {export_df.chain.unique()}")
    return export_df


def import_data(data_location, folder_path, log_file_path, isotope, standards_df):
    """Read and classify raw GC-IRMS data into standards and unknowns."""
    df = pd.read_csv(data_location)
    new_names = [str(isotope), 'area', 'chain']

    if isotope == 'dD':
        iso_rat = 'd 2H/1H'
    elif isotope == 'dC':
        iso_rat = 'd 13C/12C'
    else:
        raise ValueError('Unsupported isotope system.')

    required = ["Date", "Time", "Identifier 1", iso_rat, "Area All", "Component"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            "The raw GC-IRMS file is missing required column(s): "
            f"{', '.join(missing)}"
        )

    for idx, name in enumerate([str(iso_rat), 'Area All', 'Component']):
        df = df.rename(columns={name: new_names[idx]})

    df['date-time_true'] = df.apply(lambda row: try_parse_date(row['Date'] + ' ' + row['Time']), axis=1)
    df['date-time'] = date2num(df['date-time_true'])
    df['time_rel'] = df['date-time'] - df['date-time'].min() + 1

    if df['chain'].astype(str).str.contains(r'phthalic|PAME', case=False, na=False).any():
        pame = True
        append_to_log(log_file_path, 'PAME detected in analysis')
        print(
            'PAMEs detected.\n'
            'The calculated methanol value from the PAMEs will be displayed to the user and stored in the log file.\n'
        )
    else:
        pame = False

    linearity_std, _, linearity_ids = chain_subsetrer(df, standards_df, 'linearity')
    append_to_log(log_file_path, f'Number of linearity standards analyzed: {len(linearity_std)}')
    drift_std, _, drift_ids = chain_subsetrer(df, standards_df, 'drift')
    append_to_log(log_file_path, f'Number of Drift standards analyzed: {len(drift_std)}')

    drift_std = drift_std.sort_values('date-time_true')
    unique_time_signatures = drift_std['date-time'].unique()
    time_signatures_to_remove = unique_time_signatures[:2]
    drift_std = drift_std[~drift_std['date-time'].isin(time_signatures_to_remove)]
    append_to_log(log_file_path, 'First two drift standards ignored.')

    mask = id_mask(df, linearity_ids, col='Identifier 1', mode='all')
    unknown = df[~mask]
    mask = id_mask(unknown, drift_ids, col='Identifier 1', mode='all')
    unknown = unknown[~mask]
    unknown = unknown[~unknown['Identifier 1'].str.contains('H3+', regex=False, na=False)]

    rt_dict = ask_user_for_rt(log_file_path, df, isotope)
    if rt_dict:
        unknown = process_dataframe(unknown, rt_dict, folder_path, log_file_path)
        unknown = unknown[unknown.chain != 'None']
        linearity_std = process_dataframe(linearity_std, rt_dict, folder_path, log_file_path)
        drift_std = process_dataframe(drift_std, rt_dict, folder_path, log_file_path)
    else:
        selected_chains = get_selected_chains()
        unknown = unknown[unknown.chain.isin(['Phthalic acid'] + selected_chains)]

    unknown = unknown[~unknown.chain.isna()]
    drift_std = drift_std[~drift_std.chain.isna()]
    linearity_std = linearity_std[~linearity_std.chain.isna()]
    correction_log = make_correction_df()
    return linearity_std, drift_std, unknown, correction_log, pame
