import ipywidgets as widgets
from IPython.display import display
from .chains import load_chain_config, save_chain_config


def edit_chains():

    cfg = load_chain_config()
    sets = cfg["sets"]
    selected = cfg["selected"]

    output = widgets.Output()

    checkboxes = {}
    rows = []

    # ----------------------------------------
    # Enforce single selection behavior
    # ----------------------------------------
    def make_checkbox(key):
        cb = widgets.Checkbox(value=(key == selected), indent=False)

        def on_change(change):
            if change["new"]:
                for k, other_cb in checkboxes.items():
                    if k != key:
                        other_cb.value = False

        cb.observe(on_change, names="value")
        checkboxes[key] = cb
        return cb

    # ----------------------------------------
    # Existing groups
    # ----------------------------------------
    for key, chains in sets.items():

        cb = make_checkbox(key)

        name_label = widgets.Label(
            key,
            layout=widgets.Layout(width="150px")
        )

        chain_label = widgets.Label(
            ", ".join(chains),
            layout=widgets.Layout(width="300px")
        )

        if key.startswith("custom"):
            delete_btn = widgets.Button(
                description="✕",
                layout=widgets.Layout(width="40px")
            )

            def make_delete_callback(k):
                def delete_callback(b):
                    if checkboxes[k].value:
                        output.clear_output()
                        with output:
                            print("Cannot delete active group.")
                        return

                    del sets[k]
                    cfg["sets"] = sets
                    save_chain_config(cfg)

                    output.clear_output()
                    with output:
                        print(f"Deleted group '{k}'. Restart editor to refresh view.")
                return delete_callback

            delete_btn.on_click(make_delete_callback(key))
        else:
            delete_btn = widgets.Label("", layout=widgets.Layout(width="40px"))

        rows.extend([cb, name_label, chain_label, delete_btn])

    # ----------------------------------------
    # New group row
    # ----------------------------------------
    new_key = "__new__"
    cb_new = make_checkbox(new_key)

    new_name = widgets.Text(
        placeholder="Group name",
        layout=widgets.Layout(width="150px")
    )

    new_chains = widgets.Text(
        placeholder="Comma-separated chains",
        layout=widgets.Layout(width="300px")
    )

    rows.extend([
        cb_new,
        new_name,
        new_chains,
        widgets.Label("", layout=widgets.Layout(width="40px"))
    ])

    # ----------------------------------------
    # Grid layout (4 columns)
    # ----------------------------------------
    grid = widgets.GridBox(
        rows,
        layout=widgets.Layout(
            grid_template_columns="40px 150px 300px 40px",
            grid_row_gap="6px",
            align_items="center"
        )
    )

    # ----------------------------------------
    # Save button
    # ----------------------------------------
    save_button = widgets.Button(
        description="Save Selection",
        button_style="success"
    )

    def on_save_clicked(b):

        output.clear_output()

        selected_key = None

        for k, cb in checkboxes.items():
            if cb.value:
                selected_key = k
                break

        if selected_key is None:
            with output:
                print("You must select a group.")
            return

        # Existing group selected
        if selected_key in sets:
            cfg["selected"] = selected_key
            save_chain_config(cfg)
            with output:
                print("User selection updated.")
            return

        # New group selected
        if selected_key == "__new__":
            name = new_name.value.strip()
            chains = [c.strip() for c in new_chains.value.split(",") if c.strip()]

            if not name:
                with output:
                    print("New group must have a name.")
                return

            if not chains:
                with output:
                    print("Chain list cannot be empty.")
                return

            if name in sets:
                with output:
                    print("Group name already exists.")
                return

            sets[name] = chains
            cfg["sets"] = sets
            cfg["selected"] = name
            save_chain_config(cfg)

            with output:
                print("User selection updated.")

    save_button.on_click(on_save_clicked)

    # ----------------------------------------
    # Display
    # ----------------------------------------
    display(grid)
    display(save_button)
    display(output)