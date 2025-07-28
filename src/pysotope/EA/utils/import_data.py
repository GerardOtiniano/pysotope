import pandas as pd
import numpy as np
from pathlib import Path

def label_element_type(df):
    """
    Add a new column 'Element Type' based on which isotope column is non-NaN.

    Parameters:
        df (pd.DataFrame): Input dataframe.

    Returns:
        pd.DataFrame: Modified dataframe with new 'Element Type' column.
    """
    df = df.copy()
    df['Element Type'] = pd.Series([None] * len(df), dtype=object) 
    df.loc[df['d 15N/14N'].notna(), 'Element Type'] = 'Nitrogen'
    df.loc[df['d 13C/12C'].notna(), 'Element Type'] = 'Carbon'
    return df

def load_ea_standards() -> pd.DataFrame:
    """
    Load the EA isotope standards from the package directory.

    Returns:
        pd.DataFrame: DataFrame with columns ['Name', 'd13C', 'd15N', '%C', '%N', 'C/N']
    """
    # Find the path to the EA_standards.csv relative to this file
    HERE = Path(__file__).resolve().parent
    standards_path = HERE / "EA_standards.csv"

    if not standards_path.exists():
        raise FileNotFoundError(f"Standards file not found at: {standards_path}")

    df = pd.read_csv(standards_path, encoding = 'unicode_escape')
    return df

def add_seconds_since_start(df):
    """
    Combines 'Date' and 'Time' columns into a datetime object and computes
    seconds since the first timestamp.

    Parameters:
        df (pd.DataFrame): Input dataframe with 'Date' and 'Time' columns.

    Returns:
        pd.DataFrame: Modified dataframe with 'Seconds Since Start' column.
    """
    df = df.copy()

    # Combine 'Date' and 'Time' into a single datetime column
    df['Datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], errors='coerce')

    # Drop rows where datetime couldn't be parsed (optional)
    df = df[df['Datetime'].notna()]

    # Calculate seconds since the first timestamp
    t0 = df['Datetime'].iloc[0]
    df['Seconds Since Start'] = (df['Datetime'] - t0).dt.total_seconds()

    return df


def import_EA_data(file_path):
    df = pd.read_csv(file_path, encoding='unicode escape')
    df = label_element_type(df)
    df = add_seconds_since_start(df)
    df = df[df['Component'].isin(['N2', 'CO2'])] # Only consider valid peaks
    return df