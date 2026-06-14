from pathlib import Path

import pandas as pd


STANDARDS_DIR_CANDIDATES = [
    Path(__file__).resolve().parent.parent / "standards_manager" / "EA" / "data",
    Path(__file__).resolve().parent.parent / "standards_manager" / "data",
]

EA_CARBON_STANDARDS = "ea_RS_dC.csv"
EA_NITROGEN_STANDARDS = "ea_RS_dN.csv"


def label_element_type(df: pd.DataFrame) -> pd.DataFrame:
    """Label each EA peak by component, with isotope-channel fallback."""
    out = df.copy()
    _require_columns(out, ["Component", "d 15N/14N", "d 13C/12C"], "EA input data")
    component = out["Component"].astype(str).str.strip().str.upper()
    out['Element Type'] = pd.Series([None] * len(out), dtype=object)
    out.loc[component == "N2", 'Element Type'] = 'Nitrogen'
    out.loc[component == "CO2", 'Element Type'] = 'Carbon'

    unlabeled = out['Element Type'].isna()
    out.loc[unlabeled & out['d 15N/14N'].notna(), 'Element Type'] = 'Nitrogen'
    out.loc[unlabeled & out['d 13C/12C'].notna(), 'Element Type'] = 'Carbon'
    return out


def _require_columns(df: pd.DataFrame, cols, context: str):
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise ValueError(f"{context} is missing required column(s): {', '.join(missing)}")


def _get_standards_dir() -> Path:
    for path in STANDARDS_DIR_CANDIDATES:
        if path.exists():
            return path
    return STANDARDS_DIR_CANDIDATES[-1]


def _load_ea_standard_file(filename: str) -> pd.DataFrame:
    path = _get_standards_dir() / filename
    if not path.exists():
        raise FileNotFoundError(f"EA standards file not found: {path}")
    return pd.read_csv(path, encoding='unicode_escape')


def _standard_name_column(df: pd.DataFrame) -> str:
    for col in ("EA-IRMS Standards", "ID", "Identifier 1", "Name"):
        if col in df.columns:
            return col
    raise ValueError("EA standards table does not contain a recognizable standard-name column.")


def load_ea_standards() -> pd.DataFrame:
    """Load a merged EA standards table containing both nitrogen and carbon metadata."""
    c_df = _load_ea_standard_file(EA_CARBON_STANDARDS).copy()
    n_df = _load_ea_standard_file(EA_NITROGEN_STANDARDS).copy()

    c_name_col = _standard_name_column(c_df)
    n_name_col = _standard_name_column(n_df)

    c_df = c_df.rename(columns={c_name_col: "EA-IRMS Standards"})
    n_df = n_df.rename(columns={n_name_col: "EA-IRMS Standards"})

    merged = pd.merge(
        c_df,
        n_df,
        how="outer",
        on="EA-IRMS Standards",
        suffixes=("_C", "_N"),
    )
    # return merged
    return c_df, n_df

def add_seconds_since_start(df: pd.DataFrame) -> pd.DataFrame:
    """Add acquisition datetime and elapsed seconds columns for EA runs."""
    out = df.copy()
    _require_columns(out, ["Date", "Time"], "EA input data")
    datetime_strings = out['Date'].astype(str).str.strip() + ' ' + out['Time'].astype(str).str.strip()
    out['Datetime'] = pd.to_datetime(datetime_strings, format='mixed', errors='coerce')
    out = out[out['Datetime'].notna()]
    if out.empty:
        raise ValueError("No valid EA analysis datetimes could be parsed from the Date and Time columns.")
    t0 = out['Datetime'].iloc[0]
    out['Seconds Since Start'] = (out['Datetime'] - t0).dt.total_seconds()
    return out


def import_EA_data(file_path) -> pd.DataFrame:
    """Read and minimally normalize a raw EA-IRMS export file."""
    df = pd.read_csv(file_path, encoding='unicode escape')
    _require_columns(
        df,
        ["Date", "Time", "Identifier 1", "Component", "d 15N/14N", "d 13C/12C"],
        "EA input data",
    )
    df = label_element_type(df)
    df = add_seconds_since_start(df)
    df = df[df['Component'].isin(['N2', 'CO2'])]
    if df.empty:
        raise ValueError("No EA rows with Component values N2 or CO2 were found.")
    return df
