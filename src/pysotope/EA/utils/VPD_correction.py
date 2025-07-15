import matplotlib.pyplot as plt
import pandas as pd
from difflib import get_close_matches

def plot_measured_vs_actual(df, standards, element):
    """
    Compare measured vs actual isotopic values for a given element.

    Parameters:
        df (pd.DataFrame): Measured values with columns like 'Identifier 1' and isotope columns.
        standards (pd.DataFrame): Standards with known actual values.
        element (str): Either 'Carbon' or 'Nitrogen'.
    
    Returns:
        matplotlib.figure.Figure: The generated scatter plot figure.
    """
    # Define targets (names we care about) and column mappings
    target_names = ['Sorghum', 'Urea', 'Acetanilide', 'Wheat Flour']
    actual_col = {'Carbon': 'd13C(VPDB) value', 'Nitrogen': 'd15N(AIR) value'}
    measured_col = {'Carbon': 'd 13C/12C', 'Nitrogen': 'd 15N/14N'}
    
    if element not in actual_col:
        raise ValueError("Element must be either 'Carbon' or 'Nitrogen'")

    # Build a name map by fuzzy matching
    name_map = {}
    for target in target_names:
        matches = get_close_matches(target.lower(), standards['EA-IRMS Standards'].str.lower(), n=1, cutoff=0.5)
        if matches:
            actual_match = standards[standards['EA-IRMS Standards'].str.lower() == matches[0]]
            if not actual_match.empty:
                name_map[target.upper()] = actual_match.iloc[0]

    # Prepare matched data
    plot_data = []
    for key, row in name_map.items():
        # Look for matching measured entries in df (case-insensitive substring match)
        matches = df[df['Identifier 1'].str.upper().str.contains(key)]
        if matches.empty:
            continue

        for _, m_row in matches.iterrows():
            plot_data.append({
                "standard_name": key.title(),
                "measured": m_row[measured_col[element]],
                "actual": row[actual_col[element]]
            })

    if not plot_data:
        print("No matches found between measured data and standards.")
        return None

    # Convert to DataFrame
    plot_df = pd.DataFrame(plot_data)

    # Plot
    fig, ax = plt.subplots()
    ax.scatter(plot_df['actual'], plot_df['measured'], c='blue', label='Measured vs Actual')
    
    # Identity line
    min_val = min(plot_df['actual'].min(), plot_df['measured'].min())
    max_val = max(plot_df['actual'].max(), plot_df['measured'].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='1:1 Line')

    # Labels
    for _, row in plot_df.iterrows():
        ax.text(row['actual'], row['measured'], row['standard_name'], fontsize=8, ha='right')

    ax.set_xlabel(f"Actual {element} Value")
    ax.set_ylabel(f"Measured {element} Value")
    ax.legend()

    # return fig