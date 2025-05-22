import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
import os

def assign_chain_length():
    csv_path = input("Provide file location: ").strip()
    
    while csv_path.startswith(("'", '"')) and csv_path.endswith(("'", '"')):
        csv_path = csv_path[1:-1]
        csv_path = csv_path.strip()
    
    if not os.path.isfile(csv_path):
        print(f"File not found: {csv_path}")
        return
    
    df = pd.read_csv(csv_path)
    df = df.loc[(df['Rt'] > 400) & (df['Rt'] < 1400)].copy()

    chain_pattern = re.compile(r'C(\d{1,2})')

    def extract_chains(identifier):
        chains = chain_pattern.findall(str(identifier))
        if len(chains) >= 2:
            return tuple(sorted(map(int, chains[:2])))
        else:
            return None

    df['Chains'] = df['Identifier 1'].apply(extract_chains)

    from matplotlib import colormaps
    cmap = colormaps['tab10']

    standard_chains = df['Chains'].dropna().unique()
    chain_to_color = {tuple(ch): cmap(i) for i, ch in enumerate(standard_chains)}

    def assign_color(chains):
        if chains in chain_to_color:
            return chain_to_color[chains]
        else:
            return 'gray'

    df['Color'] = df['Chains'].apply(assign_color)

    legend_labels = {}
    for ch, color in chain_to_color.items():
        label = f'Standard C{ch[0]} & C{ch[1]}'
        legend_labels[label] = color
    legend_labels['Unknowns'] = 'gray'

    if 'Component' not in df.columns:
        df['Component'] = ''

    chain_lengths = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']
    chain_index = 0

    fig, ax = plt.subplots(figsize=(10, 6))

    for label, color in legend_labels.items():
        if label == 'Unknowns':
            subset = df[df['Color'] == color]
        else:
            parts = label.split()
            chain_nums = (int(parts[1][1:]), int(parts[3][1:]))
            subset = df[df['Chains'] == chain_nums]

        if not subset.empty:
            ax.scatter(subset['Rt'], subset['Area All'], color=color, label=label,
                       alpha=0.75, edgecolor='k', s=50)

    ax.set_xlabel('Retention Time (Rt)')
    ax.set_ylabel('Area All')
    ax.set_title('Rt vs. Area All - Click to Assign Chain Lengths')
    ax.legend(loc='upper right')

    out = widgets.Output()
    display(out)

    vertical_lines = []
    text_labels = []
    assigned_masks = []  # Keep track of masks for undo

    def on_click(event):
        nonlocal chain_index
        if event.inaxes != ax:
            return
        if chain_index >= len(chain_lengths):
            with out:
                out.clear_output()
                print("All chain lengths assigned. You may close the plot.")
            return

        rt_clicked = event.xdata
        current_chain = chain_lengths[chain_index]

        vline = ax.axvline(rt_clicked, color='red', linestyle='--')
        vertical_lines.append(vline)
        txt = ax.text(rt_clicked, ax.get_ylim()[1]*0.95, current_chain,
                      color='red', fontsize=12, rotation=90,
                      verticalalignment='top', horizontalalignment='right')
        text_labels.append(txt)
        fig.canvas.draw_idle()

        mask = (df['Rt'] >= rt_clicked - 5) & (df['Rt'] <= rt_clicked + 5)
        df.loc[mask, 'Component'] = current_chain
        assigned_masks.append(mask)

        with out:
            out.clear_output()
            print(f"Assigned {current_chain} to {mask.sum()} rows near Rt={rt_clicked:.2f} s")
            if chain_index + 1 < len(chain_lengths):
                print(f"Please click the location for {chain_lengths[chain_index + 1]}")

        chain_index += 1

        if chain_index == len(chain_lengths):
            base, ext = os.path.splitext(csv_path)
            new_path = f"{base}_chainID{ext}"
            df.drop(columns=['Color'], inplace=True)
            df.to_csv(new_path, index=False)
            with out:
                print(f"\nAll chain lengths assigned.\nSaved updated dataframe to:\n{new_path}")

    def on_key(event):
        nonlocal chain_index
        if event.key in ['backspace', 'delete']:
            if chain_index == 0:
                with out:
                    out.clear_output()
                    print("No assignments to undo.")
                return
            # Undo last assignment
            chain_index -= 1

            # Remove vertical line and label
            vertical_lines[-1].remove()
            vertical_lines.pop()
            text_labels[-1].remove()
            text_labels.pop()

            # Remove last mask assignment
            last_mask = assigned_masks.pop()
            df.loc[last_mask, 'Component'] = ''

            fig.canvas.draw_idle()

            with out:
                out.clear_output()
                print(f"Undo last assignment. Now please click the location for {chain_lengths[chain_index]}")

    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)

    with out:
        print(f"Please click the location for {chain_lengths[0]}")

    plt.show()