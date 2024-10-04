def output_results(raw_unknown, unknown, sd, folder_path, fig_path, res_path, isotope):
    # Define column name mappings for output CSV files
    column_name_mapping = {
        "Identifier 1": "Sample Name",
        "chain": "Chain Length",
        "area_mean": "Mean Area",
        "dD_mean": f"Average raw {isotope}",
        "dD_std": f"{isotope} Std Dev",
        "dD_count": "Number of Replicates",
        "drift_corrected_dD_mean": f"Drift Corrected {isotope}",
        "drift_error_mean": "Drift Error",
        "linearity_corrected_dD_mean": f"Linearity Corrected {isotope}",
        "linearity_error_mean": "Linearity Error",
        "VSMOW_dD_mean": f"VSMOW Corrected {isotope}",
        "VSMOW_error_mean": "VSMOW Error",
        "methanol_dD_mean": f"Final - Methanol Corrected {isotope}",
        "methanol_error_mean": f"Methanol Error",
        "total_uncertainty": "Total Uncertainty",
    }

    # Rename columns in unknown dataframe
    unknown_renamed = unknown.rename(columns=column_name_mapping)

    # Save unknown dataframe to CSV
    unknown_renamed.to_csv(os.path.join(res_path, "Results - sample mean.csv"), index=False)

    # Select columns for standards dataframe
    columns_to_select = ["Date", "Time", "Identifier 1", "chain", "Rt", "area", "Area 2", "Area 3", "Ampl  2", "Ampl  3", "BGD 2", "BGD 3", "time_rel", "dD", "drift_corrected_dD", "drift_error", "linearity_corrected_dD", "linearity_error", "VSMOW_dD", "vsmow_error", "total_uncertainty"]

    # Select only existing columns from sd dataframe
    existing_columns = [col for col in columns_to_select if col in sd.columns]
    standards_selected = sd[existing_columns].copy()

    # Filter and categorize standards
    lin_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C20|C28")].copy()
    lin_std_temp["Standard Type"] = "Linearity"

    drift_std_temp = standards_selected[standards_selected["Identifier 1"].str.contains("C218|C24")].copy()
    drift_std_temp["Standard Type"] = "Drift"

    standards_categorized = pd.concat([lin_std_temp, drift_std_temp])

    # Rename columns in standards dataframe
    column_rename_map = {
        "Identifier 1": "Standard ID",
        "chain": "Component",
        "Rt": "Retention time",
        "area": "Peak area",
        "Area 2": "Peak area 2",
        "Area 3": "Peak area 3",
        "Ampl  2": "Amplitude 2",
        "Ampl  3": "Amplitude 3",
        "BGD 2": "Background 2",
        "BGD 3": "Background 3",
        "time_rel": "Time relative",
        "dD": f"Raw {isotope}",
        "drift_corrected_dD": f"Drift corrected {isotope}",
        "drift_error": f"Drift {isotope} error",
        "linearity_corrected_dD": f"Linearity {isotope}",
        "linearity_error": f"Linearity {isotope} error",
        "VSMOW_dD": f"VSMOW corrected {isotope}",
        "vsmow_error": "VSMOW error",
        "methanol_dD": f"Final - Methanol corrected {isotope}",
        "methanol_error": f"Methanol error {isotope}",
        "total_uncertainty": "Total uncertainty",
    }
    standards_categorized = standards_categorized.rename(columns=column_rename_map)

    # Save standards dataframe to CSV
    standards_categorized.to_csv(os.path.join(res_path, "Results - standards.csv"), index=False)

    # Rename columns in raw_unknown dataframe
    raw_unknown_renamed = raw_unknown.rename(columns=column_rename_map)

    # Drop 'total_error' column if exists
    if "total_error" in raw_unknown.columns:
        raw_unknown_renamed = raw_unknown_renamed.drop("total_error", axis=1)

    # Save raw_unknown dataframe to CSV
    raw_unknown_renamed.to_csv(os.path.join(res_path, "Results - sample replicates.csv"), index=False)
    Total_dD_correction_plot(raw_unknown_renamed, unknown_renamed, folder_path, fig_path, isotope)
    print("\nCorrections complete :)")


def mean_values_with_uncertainty(data, sample_name_header="Identifier 1", chain_header="chain"):
    # Group by sample name and chain length
    grouped = data.groupby([sample_name_header, chain_header])

    # Calculate the mean values and count of replicates
    stats = grouped.agg(
        {"area": ["mean"], "dD": ["mean", "std", "count"], "drift_corrected_dD": "mean", "drift_error": "mean", "linearity_corrected_dD": "mean", "linearity_error": "mean", "VSMOW_dD": "mean", "VSMOW_error": "mean", "methanol_dD": ["mean", "std"], "methanol_error": "mean"}
    ).reset_index()

    # Flatten MultiIndex columns
    stats.columns = ["_".join(col).strip("_") for col in stats.columns.values]

    # Calculate SEM for methanol dD (std / sqrt(count))
    stats["replicate_dD_sem"] = stats["methanol_dD_std"] / np.sqrt(stats["dD_count"])

    # Calculate total uncertainty
    error_columns = ["drift_error_mean", "linearity_error_mean", "VSMOW_error_mean", "methanol_error_mean", "replicate_dD_sem"]
    stats["total_uncertainty"] = np.sqrt(sum(stats[col] ** 2 for col in error_columns))

    return stats
