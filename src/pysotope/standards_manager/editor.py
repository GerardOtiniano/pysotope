from pathlib import Path

import numpy as np
import pandas as pd
import ipywidgets as widgets
from IPython.display import clear_output, display


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
DATA_DIR.mkdir(exist_ok=True, parents=True)

ISO_FILENAMES = {"rs_dc.csv", "rs_dd.csv"}

DEFAULT_VSMOW = {
    "dD": {
        "ID": ["C18", "C20", "C28", "C24"],
        "type": ["drift", "linearity", "linearity", "drift"],
        "chain length": ["C18", "C20", "C28", "C24"],
        "isotope value": [-206.2, -166.7, -89.28, -179.3],
        "std": [1.7, 0.3, 1.0627, 1.7],
        "n": [5, 3, 924, 5],
        "RS accuracy check": ["False", "False", "False", "True"],
        "Use as Standard": ["True", "True", "True", "True"],
    },
    "dC": {
        "ID": ["C18", "C20", "C24", "C28"],
        "type": ["drift", "linearity", "drift", "linearity"],
        "chain length": ["C18", "C20", "C24", "C28"],
        "isotope value": [-23.24, -30.68, -26.57, -27.35793034],
        "std": [0.01, 0.02, 0.02, 0.26315684],
        "n": [5, 3, 5, 89],
        "RS accuracy check": ["False", "False", "False", "True"],
        "Use as Standard": ["True", "True", "True", "True"],
    },
}

BOOL_COLS = {
    "RS accuracy check",
    "Use as Standard",
    "Linearity Standard",
    "drift",
}

CELL_LAYOUT = widgets.Layout(width="140px")
TABLE_LAYOUT = widgets.Layout(overflow_x="auto")

def _get_bool_cols_for_df(df: pd.DataFrame) -> set[str]:
    return BOOL_COLS.intersection(df.columns)

def list_standard_files():
    return sorted(path.name for path in DATA_DIR.glob("*.csv"))

def get_standard_path(filename: str) -> Path:
    return DATA_DIR / filename

def _is_iso_file(path: Path) -> bool:
    return path.name.lower() in ISO_FILENAMES

def _ensure_default_iso_file(path: Path):
    key = "dD" if path.name.lower() == "rs_dd.csv" else "dC"
    if not path.exists():
        pd.DataFrame(DEFAULT_VSMOW[key]).to_csv(path, index=False)

def _coerce_bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "1.0", "yes", "y", "t"])
    )

# def _make_widget_for_value(val, is_bool=False):
#     if is_bool:
#         if isinstance(val, str):
#             checked = val.strip().lower() == "true"
#         else:
#             checked = bool(val) if not pd.isna(val) else False
#         return widgets.Checkbox(value=checked, layout=CELL_LAYOUT)
#
#     if pd.isna(val):
#         return widgets.Text(value="", layout=CELL_LAYOUT)
#     if isinstance(val, (int, np.integer)):
#         return widgets.IntText(value=int(val), layout=CELL_LAYOUT)
#     if isinstance(val, (float, np.floating)):
#         return widgets.FloatText(value=float(val), layout=CELL_LAYOUT)
#     return widgets.Text(value=str(val), layout=CELL_LAYOUT)

def _make_widget_for_value(val, is_bool=False):
    if is_bool:
        checked = bool(val) if not pd.isna(val) else False
        return widgets.Checkbox(value=checked, layout=CELL_LAYOUT)
    if pd.isna(val):
        return widgets.Text(value="", layout=CELL_LAYOUT)
    if isinstance(val, (int, np.integer)):
        return widgets.IntText(value=int(val), layout=CELL_LAYOUT)
    if isinstance(val, (float, np.floating)):
        return widgets.FloatText(value=float(val), layout=CELL_LAYOUT)
    return widgets.Text(value=str(val), layout=CELL_LAYOUT)

def _make_blank_widget_for_column(df: pd.DataFrame, col: str, bool_cols):
    if col in bool_cols:
        return widgets.Checkbox(value=False, layout=CELL_LAYOUT)

    ser_nonnull = df[col].dropna()
    if len(ser_nonnull) > 0 and pd.api.types.is_numeric_dtype(ser_nonnull):
        if pd.api.types.is_integer_dtype(ser_nonnull):
            return widgets.IntText(value=0, layout=CELL_LAYOUT)
        return widgets.FloatText(value=float("nan"), layout=CELL_LAYOUT)
    return widgets.Text(value="", layout=CELL_LAYOUT)

# def _build_editor(path: Path, df: pd.DataFrame, bool_cols=None):
#     bool_cols = set(bool_cols or [])
#     cell_widgets = []
#     row_boxes = []
#     table_box = widgets.VBox(layout=TABLE_LAYOUT)
#     message = widgets.Output()
#
#     def rebuild_table():
#         headers = [widgets.Label(str(c), layout=CELL_LAYOUT) for c in df.columns]
#         table_box.children = [widgets.HBox(headers)] + row_boxes
#
#     def add_row_widgets(row_widgets):
#         cell_widgets.append(row_widgets)
#         row_boxes.append(widgets.HBox(row_widgets))
#
#     for _, row in df.iterrows():
#         add_row_widgets([
#             _make_widget_for_value(row[col], is_bool=col in bool_cols)
#             for col in df.columns
#         ])
#
#     rebuild_table()
#
#     add_btn = widgets.Button(description="Add row")
#     save_btn = widgets.Button(description="Save", button_style="success")
#
#     def on_add_row(_):
#         add_row_widgets([
#             _make_blank_widget_for_column(df, col, bool_cols)
#             for col in df.columns
#         ])
#         rebuild_table()
#
#     def on_save(_):
#         data = {col: [] for col in df.columns}
#         for row_widget in cell_widgets:
#             for col, widget in zip(df.columns, row_widget):
#                 value = widget.value
#                 if isinstance(value, str) and value.strip() == "":
#                     value = np.nan
#                 data[col].append(value)
#
#         new_df = pd.DataFrame(data, columns=df.columns)
#         for col in bool_cols.intersection(new_df.columns):
#             new_df[col] = _coerce_bool_series(new_df[col])
#
#         new_df.to_csv(path, index=False)
#         with message:
#             clear_output()
#             print(f"Saved {path.name}")
#
#     add_btn.on_click(on_add_row)
#     save_btn.on_click(on_save)
#
#     return widgets.VBox([
#         widgets.HBox([add_btn, save_btn]),
#         table_box,
#         message,
#     ])

def _build_editor(path: Path, df: pd.DataFrame, bool_cols=None):
    bool_cols = set(bool_cols or [])
    cell_widgets = []
    row_boxes = []
    table_box = widgets.VBox(layout=TABLE_LAYOUT)
    message = widgets.Output()
    DELETE_BUTTON_LAYOUT = widgets.Layout(width="32px")
    def rebuild_table():
        headers = [widgets.Label(str(c), layout=CELL_LAYOUT) for c in df.columns]
        spacer = widgets.Label("", layout=DELETE_BUTTON_LAYOUT)
        table_box.children = [widgets.HBox(headers + [spacer])] + row_boxes
    def delete_row(row_widgets, row_box):
        if row_widgets in cell_widgets:
            cell_widgets.remove(row_widgets)
        if row_box in row_boxes:
            row_boxes.remove(row_box)
        rebuild_table()
    def add_row_widgets(row_widgets):
        delete_btn = widgets.Button(
            description="×",
            tooltip="Delete row",
            button_style="danger",
            layout=DELETE_BUTTON_LAYOUT,
        )
        row_box = widgets.HBox(row_widgets + [delete_btn])
        delete_btn.on_click(lambda _, rw=row_widgets, rb=row_box: delete_row(rw, rb))
        cell_widgets.append(row_widgets)
        row_boxes.append(row_box)
    for _, row in df.iterrows():
        add_row_widgets([
            _make_widget_for_value(row[col], is_bool=col in bool_cols)
            for col in df.columns])
    rebuild_table()
    add_btn = widgets.Button(description="Add row")
    save_btn = widgets.Button(description="Save", button_style="success")
    def on_add_row(_):
        add_row_widgets([
            _make_blank_widget_for_column(df, col, bool_cols)
            for col in df.columns])
        rebuild_table()
    def on_save(_):
        data = {col: [] for col in df.columns}
        for row_widget in cell_widgets:
            for col, widget in zip(df.columns, row_widget):
                value = widget.value
                if isinstance(value, str) and value.strip() == "":
                    value = np.nan
                data[col].append(value)
        new_df = pd.DataFrame(data, columns=df.columns)
        for col in bool_cols.intersection(new_df.columns):
            new_df[col] = _coerce_bool_series(new_df[col])
        new_df.to_csv(path, index=False)
        with message:
            clear_output()
            print(f"Saved {path.name}")
    add_btn.on_click(on_add_row)
    save_btn.on_click(on_save)
    return widgets.VBox([
        widgets.HBox([add_btn, save_btn]),
        table_box,
        message,])

def _render_iso_editor(path: Path) -> widgets.Widget:
    _ensure_default_iso_file(path)
    df = pd.read_csv(path, dtype={"type": str, "chain length": str, "ID": str})
    # bool_cols = {"Use as Standard", "RS accuracy check", "Linearity Standard"}
    bool_cols = _get_bool_cols_for_df(df)
    for col in bool_cols:
        df[col] = _coerce_bool_series(df[col])
    return _build_editor(path, df, bool_cols=bool_cols)

def _render_general_editor(path: Path) -> widgets.Widget:
    df = pd.read_csv(path)
    bool_cols = _get_bool_cols_for_df(df)
    for col in bool_cols:
        df[col] = _coerce_bool_series(df[col])
    return _build_editor(path, df, bool_cols=bool_cols)

# def standard_editor() -> pd.DataFrame:
#     files = list_standard_files()
#     if not files:
#         raise FileNotFoundError(f"No standards CSV files found in {DATA_DIR}")
#
#     selector = widgets.Dropdown(options=files, description="File:", layout=widgets.Layout(width="360px"))
#     output = widgets.Output()
#
#     def render_selected(filename):
#         path = get_standard_path(filename)
#         editor = _render_iso_editor(path) if _is_iso_file(path) else _render_general_editor(path)
#         with output:
#             clear_output()
#             display(editor)
#
#     selector.observe(
#         lambda change: render_selected(change["new"]) if change["name"] == "value" else None,
#         names="value",
#     )
#     display(widgets.VBox([selector, output]))
#     render_selected(selector.value)
#     return pd.read_csv(get_standard_path(selector.value))
def standard_editor() -> None:
    files = list_standard_files()
    if not files:
        raise FileNotFoundError(f"No standards CSV files found in {DATA_DIR}")
    selector = widgets.Dropdown(
        options=files,
        description="File:",
        layout=widgets.Layout(width="360px"),
    )
    output = widgets.Output()
    def render_selected(filename):
        path = get_standard_path(filename)
        editor = _render_iso_editor(path) if _is_iso_file(path) else _render_general_editor(path)
        with output:
            clear_output()
            display(editor)
    selector.observe(
        lambda change: render_selected(change["new"]) if change["name"] == "value" else None,
        names="value",
    )
    display(widgets.VBox([selector, output]))
    render_selected(selector.value)
