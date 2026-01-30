# import pandas as pd
# from pathlib import Path
# import ipywidgets as widgets
# from IPython.display import display, clear_output
# from .utils.queries import isotope_type

# HERE    = Path(__file__).resolve().parent
# CSV_DIR = HERE / "utils/vsmow_standards"
# CSV_DIR.mkdir(exist_ok=True, parents=True)

# DEFAULT_VSMOW = {
#     "dD": {
#         "type": ["drift","linearity","linearity","drift"],
#         "chain length": ["C18","C20","C28","C24"],
#         "isotope value": [-206.2, -166.7, -89.28, -179.3],
#         "std": [1.7, 0.3, 1.0627, 1.7],
#         "n":   [5,    3,    924,    5],
#         "VSMOW accuracy check": ["False","False","False","True"]
#     },
#     "dC": {
#         "type": ["drift","linearity","drift"],
#         "chain length": ["C18","C20","C24"],
#         "isotope value": [-23.24, -30.68, -26.57],
#         "std": [0.01,0.02,0.02],
#         "n":   [5,3,5],
#         "VSMOW accuracy check": ["False","False","True"]
#     }
# }

# def _csv_path(isotope: str) -> Path:
#     return CSV_DIR / f"vsmow_{isotope}.csv"

# def standard_editor() -> pd.DataFrame:
#     """
#     Show an editable grid of the standards CSV using core ipywidgets.
#     On Save, writes back to CSV and displays the updated DataFrame.
#     """
#     isotope = isotope_type()
#     path = _csv_path(isotope)
#     if not path.exists():
#         pd.DataFrame(DEFAULT_VSMOW[isotope]).to_csv(path, index=False)
#     df = pd.read_csv(path, dtype={"type":str, "chain length":str})

#     # Build header row
#     headers = [widgets.Label(f"{c}", layout=widgets.Layout(width="120px"))
#                for c in df.columns]
#     header_box = widgets.HBox(headers)

#     # Build per‐cell widgets
#     cell_widgets = []  # 2D list: rows × cols
#     for _, row in df.iterrows():
#         row_widgets = []
#         for col in df.columns:
#             val = row[col]
#             if isinstance(val, bool):
#                 w = widgets.Checkbox(value=val, layout=widgets.Layout(width="120px"))
#             elif pd.api.types.is_integer_dtype(type(val)) or isinstance(val, int):
#                 w = widgets.IntText(value=int(val), layout=widgets.Layout(width="120px"))
#             elif pd.api.types.is_float_dtype(type(val)) or isinstance(val, float):
#                 w = widgets.FloatText(value=float(val), layout=widgets.Layout(width="120px"))
#             else:
#                 w = widgets.Text(value=str(val), layout=widgets.Layout(width="120px"))
#             row_widgets.append(w)
#         cell_widgets.append(row_widgets)

#     # Pack rows into a VBox
#     row_boxes = [widgets.HBox(r) for r in cell_widgets]
#     table = widgets.VBox([header_box] + row_boxes)

#     save_btn = widgets.Button(description="Save", button_style="success")
#     out      = widgets.Output()

#     def _on_save(_):
#         clear_output(wait=True)
#         data = {col: [] for col in df.columns}
#         for rw in cell_widgets:
#             for col, w in zip(df.columns, rw):
#                 data[col].append(w.value)
#         new_df = pd.DataFrame(data)
#         new_df["VSMOW accuracy check"] = new_df["VSMOW accuracy check"].astype(str).str.lower()=="true"
#         new_df.to_csv(path, index=False)
#         clear_output()
#         # with out:
#         #     clear_output()
#         #     display(new_df)

#     save_btn.on_click(_on_save)
#     display(widgets.VBox([table, save_btn]))
#     return df
import pandas as pd
from pathlib import Path
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output
from .utils.queries import isotope_type

HERE    = Path(__file__).resolve().parent
CSV_DIR = HERE / "utils/vsmow_standards"
CSV_DIR.mkdir(exist_ok=True, parents=True)

DEFAULT_VSMOW = {
    "dD": {
        "type": ["drift","linearity","linearity","drift"],
        "chain length": ["C18","C20","C28","C24"],
        "isotope value": [-206.2, -166.7, -89.28, -179.3],
        "std": [1.7, 0.3, 1.0627, 1.7],
        "n":   [5,    3,    924,    5],
        "VSMOW accuracy check": ["False","False","False","True"]
    },
    "dC": {
        "type": ["drift","linearity","drift"],
        "chain length": ["C18","C20","C24"],
        "isotope value": [-23.24, -30.68, -26.57],
        "std": [0.01,0.02,0.02],
        "n":   [5,3,5],
        "VSMOW accuracy check": ["False","False","True"]
    }
}

def _csv_path(isotope: str) -> Path:
    return CSV_DIR / f"vsmow_{isotope}.csv"

def standard_editor() -> pd.DataFrame:
    """
    Show an editable grid of the standards CSV using core ipywidgets.
    Buttons:
      - Add row: appends a blank (NaN) row to the editor.
      - Save: writes back to CSV.
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
    
    bool_cols = {"Use as Standard", "VSMOW accuracy check"} 
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