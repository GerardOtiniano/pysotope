import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
import os

def assign_chain_length(output_location=None, rt_min=0, rt_max=2500, chain_lengths = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']):
    csv_path = input("Provide file location: ").strip()

    while csv_path.startswith(("'", '"')) and csv_path.endswith(("'", '"')):
        csv_path = csv_path[1:-1]
        csv_path = csv_path.strip()

    if not os.path.isfile(csv_path):
        print(f"File not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    df = df.loc[(df['Rt'] > rt_min) & (df['Rt'] < rt_max)].copy()
    chain_pattern = re.compile(r'[CN]\d+[A-Za-z0-9]*')

    def extract_chains(identifier):
        matches = chain_pattern.findall(str(identifier))
        if len(matches) >= 2:
            return tuple(sorted(matches[:2]))  # Sort lexicographically (strings)
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
    # for ch, color in chain_to_color.items():
    #     label = f'Standard C{ch[0]} & C{ch[1]}'
    #     legend_labels[label] = color
    for ch, color in chain_to_color.items():
        # ch is now tuple of strings like ('C20N2', 'C28')
        label = f'Standard {ch[0]} & {ch[1]}'
        legend_labels[label] = color
    legend_labels['Unknowns'] = 'gray'

    if 'Component' not in df.columns:
        df['Component'] = ''

    # chain_lengths = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']
    chain_index = 0

    fig, ax = plt.subplots(figsize=(10, 6))

    for label, color in legend_labels.items():
        if label == 'Unknowns':
            subset = df[df['Color'] == color]
            zoro = 0
            alphy = 0.4
        else:
            parts = label.split()
            # parts example: ['Standard', 'C20N2', '&', 'C28']
            chain_1 = parts[1]
            chain_2 = parts[3]
            subset = df[df['Chains'].apply(lambda x: x == (chain_1, chain_2))]
            zoro = 1
            alphy = 0.65

        if not subset.empty:
            ax.scatter(subset['Rt'], subset['Area All'], color=color, label=label,
                       alpha=alphy, edgecolor='k', s=50, zorder=zoro)
    ax.set_xlabel('Retention Time (Rt)')
    ax.set_ylabel('Area All')
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
            # Determine output path
            base_name = os.path.basename(csv_path)
            base, ext = os.path.splitext(base_name)
            if output_location is None:
                new_path = f"{os.path.splitext(csv_path)[0]}_chainID{ext}"
            else:
                # Make sure output_location exists
                os.makedirs(output_location, exist_ok=True)
                new_path = os.path.join(output_location, f"{base}_chainID{ext}")

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