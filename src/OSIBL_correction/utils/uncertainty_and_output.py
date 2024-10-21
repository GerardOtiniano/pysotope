import numpy as np
import os
import pandas as pd
from .figures import total_dD_correction_plot

def output_results(raw_unknown, unknown, sd, unknown_pame, folder_path, fig_path, res_path, isotope, pame):
    # Define column name mappings for output CSV files
    column_name_mapping = {
        'Identifier 1': 'Sample Name',
        'chain': 'Chain Length',
        'area_mean': 'Mean Area',
        'dD_mean': f'Average raw {isotope}',
        'dD_std': f'{isotope} Std Dev',
        'dD_count': 'Number of Replicates',
        'drift_corrected_dD_mean': f'Drift Corrected {isotope}',
        'drift_error_mean': 'Drift Error',
        'linearity_corrected_dD_mean': f'Linearity Corrected {isotope}',
        'linearity_error_mean': 'Linearity Error',
        'VSMOW_dD_mean': f'VSMOW Corrected {isotope}',
        'VSMOW_error_mean': 'VSMOW Error',
        'methanol_dD_mean': f'Methanol Corrected {isotope}',
        'methanol_error_mean': f'Methanol Error',
        'PAME_methanol_dD_mean': f'PAME Methanol Calculated {isotope}',
        'PAME_methanol_dD_std': f'PAME Methanol Calculated Error {isotope}',
        'replicate_dD_sem': f'Mean replicate std {isotope}',
        'total_uncertainty': 'Total Uncertainty'
    }
    
    # Rename columns in unknown dataframe
    unknown_renamed = unknown.rename(columns=column_name_mapping)
    if pame:
        unknown_pame_renamed = unknown_pame.rename(columns=column_name_mapping)
    if isotope=='dD':
        unknown_renamed.insert(len(unknown_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_renamed['Methanol Corrected '+str(isotope)])
        if pame:
            unknown_pame_renamed.insert(len(unknown_pame_renamed.columns) - 1, 'Corrected Methanol '+str(isotope), unknown_pame_renamed['PAME Methanol Calculated '+str(isotope)])
    else:
        unknown_renamed.reset_index(drop=True, inplace=True)
        unknown_renamed.insert(len(unknown_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_renamed['VSMOW Corrected dC'])
        if pame:
            unknown_pame_renamed.reset_index(drop=True, inplace=True)
            unknown_pame_renamed.insert(len(unknown_pame_renamed.columns) - 1, 'Corrected '+str(isotope), unknown_pame_renamed['VSMOW Corrected dC'])
    # Save unknown dataframe to CSV
    unknown_renamed.reset_index(drop=True, inplace=True)
    unknown_renamed.to_csv(os.path.join(res_path,'Results - sample mean.csv'), index=False)
    if pame:
        unknown_pame_renamed.reset_index(drop=True, inplace=True)
        unknown_pame_renamed.to_csv(os.path.join(res_path,'PAME Results - sample mean.csv'), index=False)
    # Select columns for standards dataframe
    columns_to_select = [
        "Date", "Time", "Identifier 1", "chain", "Rt", "area", "Area 2", "Area 3",
        "Ampl  2", "Ampl  3", "BGD 2", 'BGD 3', "time_rel", "dD", "drift_corrected_dD",
        "drift_error", "linearity_corrected_dD", "linearity_error", "VSMOW_dD", "vsmow_error",
        "total_uncertainty"
    ]

    # Select only existing columns from sd dataframe
    existing_columns = [col for col in columns_to_select if col in sd.columns]
    standards_selected = sd[existing_columns].copy()

    # Filter and categorize standards
    lin_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C20|C28")].copy()
    lin_std_temp['Standard Type'] = "Linearity"
    
    drift_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C218|C24")].copy()
    drift_std_temp['Standard Type'] = "Drift"
    
    standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

    # Rename columns in standards dataframe
    column_rename_map = {
        "Identifier 1": "Standard ID", "chain": "Component", "Rt": "Retention time",
        "area": "Peak area", "Area 2": "Peak area 2", "Area 3": "Peak area 3",
        "Ampl  2": "Amplitude 2", "Ampl  3": "Amplitude 3", "BGD 2": "Background 2",
        "BGD 3": "Background 3", "time_rel": "Time relative", "dD": f"Raw {isotope}",
        "drift_corrected_dD": f"Drift corrected {isotope}", "drift_error": f"Drift {isotope} error",
        "linearity_corrected_dD": f"Linearity {isotope}", "linearity_error": f"Linearity {isotope} error",
        "VSMOW_dD": f"VSMOW corrected {isotope}", "vsmow_error": "VSMOW error",
        "methanol_dD": f"Methanol corrected {isotope}", "methanol_error": f"Methanol error {isotope}",
        "total_uncertainty": "Total uncertainty"
    }
    standards_categorized = standards_categorized.rename(columns=column_rename_map)
    standards_categorized.insert(len(standards_categorized.columns) - 1, 'Corrected '+str(isotope), standards_categorized[f"VSMOW corrected {isotope}"])
    # Save standards dataframe to CSV
    standards_categorized.to_csv(os.path.join(res_path, 'Results - standards.csv'), index=False)

    # Rename columns in raw_unknown dataframe
    raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

    # Drop 'total_error' column if exists
    if 'total_error' in raw_unknown.columns:
        raw_unknown_renamed = raw_unknown_renamed.drop('total_error', axis=1)

    # Save raw_unknown dataframe to CSV
    raw_unknown_renamed.to_csv(os.path.join(res_path, 'Results - sample replicates.csv'), index=False)
    total_dD_correction_plot(raw_unknown_renamed, unknown_renamed, folder_path, fig_path, isotope)
    print("\nCorrections complete :)")

def mean_values_with_uncertainty(data, iso, sample_name_header="Identifier 1", chain_header="chain"):
    # Group by sample name and chain length
    grouped = data.groupby([sample_name_header, chain_header])

    # Start building the aggregation dictionary based on iso
    agg_dict = {
        'area': ['mean'],
        'drift_corrected_dD': 'mean',
        'drift_error': 'mean',
        'linearity_corrected_dD': 'mean',
        'linearity_error': 'mean',
        'VSMOW_dD': ['mean', 'std'],
        'VSMOW_error': 'mean'
    }

    if iso == "dC":
        agg_dict.update({
            'dC': ['mean', 'std', 'count'],
        })
    else:
        agg_dict.update({
            'dD': ['mean', 'std', 'count'],
        })

        # Check if 'methanol_dD' exists in the dataframe
        if 'methanol_dD' in data.columns:
            agg_dict.update({
                'methanol_dD': ['mean', 'std'],
                'methanol_error': 'mean'
            })

        # Check if 'pame_methanol_dD' exists in the dataframe
        if 'PAME_methanol_dD' in data.columns:
            agg_dict.update({
                'PAME_methanol_dD': ['mean', 'std']
            })

    # Apply the aggregation
    stats = grouped.agg(agg_dict).reset_index()

    # Flatten MultiIndex columns
    stats.columns = ['_'.join(col).strip('_') for col in stats.columns.values]
    
    # Calculate SEM for methanol dD (std / sqrt(count)) if iso is 'dD' and methanol_dD exists
    if iso == 'dD':
        if 'methanol_dD_std' in stats.columns and 'dD_count' in stats.columns:
            stats['replicate_dD_sem'] = stats['methanol_dD_std'] / np.sqrt(stats['dD_count'])
        elif 'PAME_methanol_dD_std' in stats.columns and 'dD_count' in stats.columns:
            stats['replicate_dD_sem'] = stats['PAME_methanol_dD_std'] / np.sqrt(stats['dD_count'])
    else:
        stats['replicate_dC_sem'] = stats['VSMOW_error_mean'] / np.sqrt(stats['dC_count'])

    # Calculate total uncertainty
    if iso == 'dD':
        if 'methanol_error_mean' in stats.columns:
            error_columns = ['drift_error_mean', 'linearity_error_mean', 'VSMOW_error_mean', 'methanol_error_mean', 'replicate_dD_sem']
        elif 'PAME_methanol_dD_std' in stats.columns:
            error_columns = ['drift_error_mean', 'linearity_error_mean', 'VSMOW_error_mean', 'replicate_dD_sem']
    else:
        error_columns = ['drift_error_mean', 'linearity_error_mean', 'VSMOW_error_mean', 'replicate_dC_sem']

    stats['total_uncertainty'] = np.sqrt(np.sum([stats[col] ** 2 for col in error_columns], axis=0))

    return stats