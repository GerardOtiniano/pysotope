import pandas as pd
from pathlib import Path
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
    On Save, writes back to CSV and displays the updated DataFrame.
    """
    isotope = isotope_type()
    path = _csv_path(isotope)
    if not path.exists():
        pd.DataFrame(DEFAULT_VSMOW[isotope]).to_csv(path, index=False)
    df = pd.read_csv(path, dtype={"type":str, "chain length":str})

    # Build header row
    headers = [widgets.Label(f"{c}", layout=widgets.Layout(width="120px"))
               for c in df.columns]
    header_box = widgets.HBox(headers)

    # Build per‐cell widgets
    cell_widgets = []  # 2D list: rows × cols
    for _, row in df.iterrows():
        row_widgets = []
        for col in df.columns:
            val = row[col]
            if isinstance(val, bool):
                w = widgets.Checkbox(value=val, layout=widgets.Layout(width="120px"))
            elif pd.api.types.is_integer_dtype(type(val)) or isinstance(val, int):
                w = widgets.IntText(value=int(val), layout=widgets.Layout(width="120px"))
            elif pd.api.types.is_float_dtype(type(val)) or isinstance(val, float):
                w = widgets.FloatText(value=float(val), layout=widgets.Layout(width="120px"))
            else:
                w = widgets.Text(value=str(val), layout=widgets.Layout(width="120px"))
            row_widgets.append(w)
        cell_widgets.append(row_widgets)

    # Pack rows into a VBox
    row_boxes = [widgets.HBox(r) for r in cell_widgets]
    table = widgets.VBox([header_box] + row_boxes)

    save_btn = widgets.Button(description="Save", button_style="success")
    out      = widgets.Output()

    def _on_save(_):
        clear_output(wait=True)
        data = {col: [] for col in df.columns}
        for rw in cell_widgets:
            for col, w in zip(df.columns, rw):
                data[col].append(w.value)
        new_df = pd.DataFrame(data)
        new_df["VSMOW accuracy check"] = new_df["VSMOW accuracy check"].astype(str).str.lower()=="true"
        new_df.to_csv(path, index=False)
        clear_output()
        # with out:
        #     clear_output()
        #     display(new_df)

    save_btn.on_click(_on_save)
    display(widgets.VBox([table, save_btn]))
    return df

# def load_standards(isotope: str="dD") -> pd.DataFrame:
#     """
#     Read the edited CSV (or dump+read defaults first run).
#     """
#     path = _csv_path(isotope)
#     if not path.exists():
#         return standard_editor(isotope)
#     df = pd.read_csv(path, dtype={"type":str, "chain length":str})
#     df["VSMOW accuracy check"] = df["VSMOW accuracy check"].astype(str).str.lower()=="true"
#     return df