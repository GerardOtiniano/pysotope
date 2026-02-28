import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.widgets import TextBox
from IPython import get_ipython

def load_standards(isotope: str) -> pd.DataFrame:
    HERE    = Path(__file__).resolve().parent
    CSV_DIR = HERE / "reference_standards"
    CSV_DIR.mkdir(exist_ok=True, parents=True)

    path = CSV_DIR / f"RS_{isotope}.csv"
    df   = pd.read_csv(path, dtype={"type":str, "chain length":str})
    df["RS accuracy check"] = df["RS accuracy check"].astype(str).str.lower()=="true"
    return df

def pick_chain_retention(
    isotope: str,
    df: pd.DataFrame,
    time_col: str       = "Rt",
    amp_col: str        = "area",
    kind_col: str       = "kind",
    chain_order: list   = None,
    picker_tol: float   = 5,
    textbox_rect: tuple = (0.25, 0.90, 0.15, 0.05),
    figsize: tuple      = (8, 4),
) -> dict:
    """
    Interactive picker/TextBox in JupyterLab:
    let the user click points and type chain names.
    """
    # ── 0) Ensure we’re using the widget backend ───────────────────────────
    try:
        ip = get_ipython()
        ip.run_line_magic('matplotlib', 'widget')
        plt.ion()
    except Exception:
        # fallback if not in IPython
        pass

    # ── 1) load standards & split into drift/linearity/sample ─────────────
    standards = load_standards(isotope)
    types     = standards["type"].unique()

    measure = {}
    for t in types:
        chains  = standards.loc[standards["type"]==t, "chain length"].unique().tolist()
        pat     = "|".join(map(re.escape, chains))
        mask    = df["Identifier 1"].astype(str).str.contains(pat, regex=True)
        sub     = df.loc[mask].copy()
        sub[kind_col] = t
        measure[t]    = sub

    all_chains  = standards["chain length"].unique().tolist()
    all_pat     = "|".join(map(re.escape, all_chains))
    mask_all    = df["Identifier 1"].astype(str).str.contains(all_pat, regex=True)
    others      = df.loc[~mask_all].copy()
    others[kind_col] = "sample"
    measure["sample"] = others

    if chain_order is None:
        def keyfn(x):
            m = re.search(r"\d+", x)
            return int(m.group()) if m else x
        chain_order = sorted(all_chains, key=keyfn)

    # ── 2) build plot_df and scatter ─────────────────────────────────────
    plot_df = pd.concat([measure.get("drift",pd.DataFrame()),
                         measure.get("linearity",pd.DataFrame()),
                         measure["sample"]],
                        ignore_index=True)

    fig, ax = plt.subplots(figsize=figsize)
    kind_map = {"drift":"blue","linearity":"red","sample":"black"}
    for k, grp in plot_df.groupby(kind_col):
        ax.scatter(
            grp[time_col], grp[amp_col],
            c      = kind_map.get(k,"gray"),
            label  = k,
            picker = picker_tol,
            alpha  = .7,
            ec     = 'k'
        )
    ax.set_xlabel(time_col)
    ax.set_ylabel(amp_col)
    ax.legend()

    # ── 3) picker callback + TextBox ───────────────────────────────────────
    picked  = {}
    vline   = ax.axvline(np.nan, color="gray", ls="--")
    text_ax = None

    def is_order_ok(chain, rt):
        i   = chain_order.index(chain)
        prev = [picked[c] for c in chain_order[:i]   if c in picked]
        nxt  = [picked[c] for c in chain_order[i+1:] if c in picked]
        if prev and rt <= max(prev): return False
        if nxt  and rt >= min(nxt):  return False
        return True

    def on_pick(event):
        nonlocal text_ax
        if text_ax:
            text_ax.remove()
            text_ax = None

        ind = event.ind[0]
        rt  = plot_df.iloc[ind][time_col]
        vline.set_xdata(rt)
        fig.canvas.draw_idle()

        # show TextBox
        text_ax = fig.add_axes(textbox_rect)
        tb      = TextBox(text_ax, "Chain: ")
        def submit(txt):
            nonlocal text_ax
            text_ax.remove()
            text_ax = None
            fig.canvas.draw_idle()

            chain = txt.strip().upper()
            if chain not in chain_order:
                print(f"⚠ '{chain}' not in {chain_order}")
            elif chain in picked:
                print(f"⚠ '{chain}' already at {picked[chain]:.2f}")
            elif not is_order_ok(chain, rt):
                print(f"⚠ order violation: {chain} @ {rt:.2f}")
            else:
                picked[chain] = rt
                print(f"✅ {chain} → {rt:.2f}")
                if len(picked) == len(chain_order):
                    print("\n🎉 All chains labelled:", picked)

        tb.on_submit(submit)

    fig.canvas.mpl_connect("pick_event", on_pick)
    plt.show()
    return picked