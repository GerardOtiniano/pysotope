import pandas as pd
from pathlib import Path
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output
from .utils.queries import isotope_type

HERE    = Path(__file__).resolve().parent
CSV_DIR = HERE / "utils/reference_standards"
CSV_DIR.mkdir(exist_ok=True, parents=True)

DEFAULT_VSMOW = {
    "dD": {
        "type": ["drift","linearity","linearity","drift"],
        "chain length": ["C18","C20","C28","C24"],
        "isotope value": [-206.2, -166.7, -89.28, -179.3],
        "std": [1.7, 0.3, 1.0627, 1.7],
        "n":   [5,    3,    924,    5],
        "RS accuracy check": ["False","False","False","True"]
    },
    "dC": {
        "type": ["drift","linearity","drift"],
        "chain length": ["C18","C20","C24"],
        "isotope value": [-23.24, -30.68, -26.57],
        "std": [0.01,0.02,0.02],
        "n":   [5,3,5],
        "RS accuracy check": ["False","False","True"]
    }
}

def _csv_path(isotope: str) -> Path:
    return CSV_DIR / f"RS_{isotope}.csv"

def standard_editor() -> pd.DataFrame:
    """
    Launch an interactive editor for VSMOW reference standards.

    This function opens a Jupyter-based graphical interface (ipywidgets)
    that allows users to view, modify, add, or remove reference standard
    entries used in isotope normalization. The editor operates directly
    on the CSV files stored within ``utils/reference_standards`` and
    updates them in-place upon user confirmation.

    Parameters
    ----------

    Functionality
    -------------
    The editor provides:

    - Display of current reference standards as an editable table.
    - Modification of:
        - Standard type (e.g., drift, linearity, VSMOW)
        - Chain length designation
        - Accepted isotope value
        - Associated analytical uncertainty (std)
        - Boolean inclusion flag ("Use as Standard")
    - Addition of new rows populated with NaN values.
    - Deletion of existing rows.
    - Saving updates to the corresponding CSV file.
    - Automatic creation of default standards if none exist.

    Data Structure
    --------------
    The standards file contains the following required columns:

    - ``type`` : str
        Role of the standard (e.g., drift, linearity).
    - ``chain length`` : str
        Chain designation (e.g., C16, C18, PAME).
    - ``isotope value`` : float
        Accepted isotope value for calibration.
    - ``std`` : float
        Analytical uncertainty of the accepted value.
    - ``Use as Standard`` : bool
        Flag indicating whether the entry is used in regression fitting.

    Returns
    -------
    None
        This function is interactive and does not return a value.
        All modifications are written directly to disk.

    Notes
    -----
    - This function is intended for use in Jupyter environments.
    - Changes take effect in subsequent runs of ``iso_process``.
    - Incorrect modification of standards may affect drift,
      linearity, and VSMOW normalization accuracy.
    - The editor ensures required columns are preserved but does not
      enforce statistical validation of entered values.

    See Also
    --------
    iso_process : Main isotope processing pipeline.
    utils.reference_standards : Directory containing standards CSV files.
    """
    isotope = isotope_type()
    path = _csv_path(isotope)
    if not path.exists():
        pd.DataFrame(DEFAULT_VSMOW[isotope]).to_csv(path, index=False)

    df = pd.read_csv(path, dtype={"type": str, "chain length": str})

    # Keep these in outer scope so callbacks can modify them
    cell_widgets = []          # 2D list: rows × cols (widget objects)
    row_boxes = []             # list of HBox rows
    table_box = widgets.VBox() # container we can redraw

    def _make_widget_for_value(val, colname=None):
        """Create an appropriate widget for a single cell value."""
        # Treat NaN as "blank" in the UI
        if pd.isna(val):
            # If you prefer blanks everywhere, you can default to Text("")
            # but this keeps numeric columns numeric by default.
            return widgets.Text(value="", layout=widgets.Layout(width="120px"))

        if isinstance(val, (bool, np.bool_)):
            return widgets.Checkbox(value=bool(val), layout=widgets.Layout(width="120px"))
        if isinstance(val, (int, np.integer)):
            return widgets.IntText(value=int(val), layout=widgets.Layout(width="120px"))
        if isinstance(val, (float, np.floating)):
            return widgets.FloatText(value=float(val), layout=widgets.Layout(width="120px"))
        return widgets.Text(value=str(val), layout=widgets.Layout(width="120px"))

    def _make_blank_widget_for_column(col):
        """
        New-row widgets: populate with NaN -> represented as blank-ish widgets.
        We'll use:
          - Checkbox False for the boolean-ish column
          - IntText/FloatText nan isn't allowed; use Text("") or sensible defaults
        """
        if col in bool_cols:
            return widgets.Checkbox(value=False, layout=widgets.Layout(width="120px"))

        # Heuristic: if existing non-null values look numeric, use numeric widget default.
        ser = df[col]
        ser_nonnull = ser.dropna()
        if len(ser_nonnull) > 0 and pd.api.types.is_numeric_dtype(ser_nonnull):
            if pd.api.types.is_integer_dtype(ser_nonnull):
                return widgets.IntText(value=0, layout=widgets.Layout(width="120px"))
            return widgets.FloatText(value=float("nan"), layout=widgets.Layout(width="120px"))

        # Otherwise text blank
        return widgets.Text(value="", layout=widgets.Layout(width="120px"))

    def _build_table():
        """(Re)build the displayed table from current widgets."""
        headers = [widgets.Label(f"{c}", layout=widgets.Layout(width="120px")) for c in df.columns]
        header_box = widgets.HBox(headers)
        table_box.children = [header_box] + row_boxes
    
    bool_cols = {"Use as Standard", "RS accuracy check"} 
    # Initialize widgets from df
    for _, row in df.iterrows():
        row_widgets = []
        for col in df.columns:
            val = row[col]
            # Special-case the boolean-ish column if it was stored as string
            if col in bool_cols:
                # Accept bools or strings like "True"/"False"
                if isinstance(val, str):
                    v = val.strip().lower() == "true"
                else:
                    v = bool(val) if not pd.isna(val) else False
                w = widgets.Checkbox(value=v, layout=widgets.Layout(width="120px"))
            else:
                w = _make_widget_for_value(val, colname=col)
            row_widgets.append(w)
        cell_widgets.append(row_widgets)
        row_boxes.append(widgets.HBox(row_widgets))

    _build_table()

    add_btn  = widgets.Button(description="Add row", button_style="")
    save_btn = widgets.Button(description="Save", button_style="success")
    out      = widgets.Output()

    def _on_add_row(_):
        # Append new widgets (blank/NaN row)
        new_row_widgets = []
        for col in df.columns:
            new_row_widgets.append(_make_blank_widget_for_column(col))
        cell_widgets.append(new_row_widgets)
        row_boxes.append(widgets.HBox(new_row_widgets))
        _build_table()

    def _on_save(_):
        # Harvest values from widgets
        data = {col: [] for col in df.columns}
        for rw in cell_widgets:
            for col, w in zip(df.columns, rw):
                v = w.value
                # Convert empty strings to NaN (so "blank" really is NaN in CSV)
                if isinstance(v, str) and v.strip() == "":
                    v = np.nan
                data[col].append(v)

        new_df = pd.DataFrame(data)

        # Normalize boolean column robustly
        for bcol in (set(new_df.columns) & bool_cols):
            new_df[bcol] = (
                new_df[bcol]
                .astype(str)
                .str.strip()
                .str.lower()
                .isin(["true", "1", "yes", "y", "t"])
            )

        new_df.to_csv(path, index=False)

        # Optional: show what was saved
        with out:
            clear_output()
            display(new_df)

    add_btn.on_click(_on_add_row)
    save_btn.on_click(_on_save)

    controls = widgets.HBox([add_btn, save_btn])
    display(widgets.VBox([table_box, controls, out]))
    return df