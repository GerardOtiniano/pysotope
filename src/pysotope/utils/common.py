import datetime

import pandas as pd


def try_parse_date(date_str):
    """Parse several common GC/EA export datetime formats."""
    formats = [
        "%m/%d/%Y %H:%M:%S",
        "%m/%d/%y %H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%d/%m/%Y %H:%M:%S",
    ]

    for fmt in formats:
        try:
            return datetime.datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None


def id_mask(df: pd.DataFrame, ids, col="Identifier 1", mode="all"):
    """Build a boolean mask matching one or more identifiers in a dataframe column."""
    ids = [str(x) for x in ids if pd.notna(x) and str(x).strip() != ""]
    if len(ids) == 0:
        return pd.Series(False, index=df.index)

    series = df[col].astype(str)

    if mode == "all":
        mask = pd.Series(True, index=df.index)
        for value in ids:
            mask &= series.str.contains(value, regex=False, na=False)
        return mask

    if mode == "any":
        mask = pd.Series(False, index=df.index)
        for value in ids:
            mask |= series.str.contains(value, regex=False, na=False)
        return mask

    raise ValueError("mode must be 'all' or 'any'")


def append_to_log(log_file_path, log_message):
    with open(log_file_path, 'a', encoding='utf-8', errors='replace') as log_file:
        print(f" {log_message}", file=log_file)
