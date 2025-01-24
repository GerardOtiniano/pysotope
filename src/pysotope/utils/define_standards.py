import tkinter as tk
import pandas as pd

# Default values
default_stds = {
    "drift_stds": ["C20", "C28"],
    "linearity_stds": ["C18", "C24"]    
}

vsmow_dD = {
    "type":           ["drift", "linearity", "linearity", "drift"],
    "id":             ["C18", "C20", "C28", "C24"],
    "isotope value":  [-206.2, -166.7, -89.28190994, -179.3],
    "std":            [1.7, 0.3, 1.062744893, 1.7],
    "n":              [5, 3, 924, 5],
    "VSMOW accuracy check": ["False", "False", "False", "True"]}

vsmow_dC = {
        "type":           ["drift", "linearity", "drift"],
        "id":             ["C18", "C20", "C24"],
        "isotope value" : [-23.24, -30.68, -26.57],
        "std":            [0.01, 0.02, 0.02],
        "n":              [5,3,5],
        "VSMOW accuracy check": ["False", "False", "True"]}
        
class Editor:
    def __init__(self, alt_stds, isotope):
        self.alt_stds = alt_stds
        self.isotope = isotope
        self.standards_df = pd.DataFrame()  # Initialize the DataFrame
        self.root = tk.Tk()
        self.root.title("Edit Known Values")
        self.create_widgets()
    
    def create_widgets(self):
        if not self.alt_stds:
            # If alt_stds is False, no need to create widgets
            return
        # Select the appropriate data dictionary based on isotope type
        if self.isotope == "dD":
            self.data_dict = vsmow_dD
            isotope_header = "δD VSMOW"
        elif self.isotope == "dC":
            self.data_dict = vsmow_dC
            isotope_header = "δ13C VSMOW"
        else:
            print(f"Unsupported isotope type: {self.isotope}")
            self.root.destroy()
            return
        
        # Define columns
        columns = ["type", "chain length", isotope_header, "standard deviation", "n", "accuracy check"]
        
        # Create column headers
        for j, column in enumerate(columns):
            header = tk.Label(self.root, text=column)
            header.grid(row=0, column=j, padx=5, pady=5)  # Place headers in the first row
    
        self.entries = []  # List to store entry widgets
    
        num_rows = len(self.data_dict["id"])  # Determine number of rows based on data dictionary
    
        for i in range(num_rows):  # Create rows based on the data dictionary
            row_entries = []
            for j in range(len(columns)):
                entry = tk.Entry(self.root)
                entry.grid(row=i+1, column=j, padx=5, pady=5)
    
                # Pre-fill values based on row and column
                if j == 0:  # Type
                    entry.insert(0, self.data_dict["type"][i])
                elif j == 1:  # Chain length
                    entry.insert(0, self.data_dict["id"][i])
                elif j == 2:  # δD or δ13C
                    entry.insert(0, self.data_dict["isotope value"][i])
                elif j == 3:  # std
                    entry.insert(0, self.data_dict["std"][i])
                elif j == 4:  # n
                    entry.insert(0, self.data_dict["n"][i])
                elif j == 5:  # accuracy check
                    entry.insert(0, self.data_dict["accuracy check"][i])
            
                row_entries.append(entry)
    
            self.entries.append(row_entries)  # Store row entries for access later

        # Create a button to save the values
        save_button = tk.Button(self.root, text="Set Values", command=self.save_and_close)
        save_button.grid(row=5, columnspan=len(columns), pady=10)  # Place save button below the table

    def save_and_close(self):
        """Save values and close the editor window."""
        data = {
            "type": [],
            "chain length": [],
            "isotope value": [],  # Will be renamed based on isotope type
            "std": [],
            "n": [],
            "accuracy check": []
        }
    
        # Determine the correct column name based on isotope type
        if self.isotope == "dD":
            isotope_column = "δD"
        elif self.isotope == "dC":
            isotope_column = "δ13C"
        else:
            print("Error: ", f"Unsupported isotope type: {self.isotope}")
            return
    
        try:
            for i in range(len(self.data_dict["id"])):  # Loop through the rows
                data["type"].append(self.entries[i][0].get())
                data["chain length"].append(self.entries[i][1].get())
                data["isotope value"].append(float(self.entries[i][2].get()))
                data["std"].append(float(self.entries[i][3].get()))
                data["n"].append(int(self.entries[i][4].get()))
                
                # Handle "accuracy check" as a boolean
                accuracy_value = self.entries[i][5].get()
                if accuracy_value.lower() == "true":
                    data["accuracy check"].append(True)
                elif accuracy_value.lower() == "false":
                    data["accuracy check"].append(False)
                else:
                    raise ValueError(f"Invalid accuracy check value: {accuracy_value}")
        except ValueError as e:
            print("Input Error: ", f"Invalid input: {e}")
            return
    
        # Rename the isotope column in the DataFrame
        self.standards_df = pd.DataFrame(data).rename(columns={isotope_column: self.isotope})
        self.root.quit()  # Exit the main loop
        self.root.destroy()  # Destroy the window
    
    def run(self):
        self.root.mainloop()  # Start the Tkinter event loop
        return self.standards_df  # Return the DataFrame after the loop exits

def open_editor(alt_stds, isotope):
    if alt_stds:
        editor = Editor(alt_stds, isotope)
        df = editor.run()
        return df
    else:
        # Build the DataFrame directly from the appropriate dictionary
        if isotope == "dD":
            data = {
                "type": vsmow_dD["type"],
                "chain length": vsmow_dD["id"],
                "δD": vsmow_dD["isotope value"],
                "std": vsmow_dD["std"],
                "n": vsmow_dD["n"],
                "accuracy check": [val.lower() == "true" for val in vsmow_dD["accuracy check"]]
            }
            return pd.DataFrame(data)
        elif isotope == "dC":
            data = {
                "type": vsmow_dC["type"],
                "chain length": vsmow_dC["id"],
                "δ13C": vsmow_dC["isotope value"],
                "std": vsmow_dC["std"],
                "n": vsmow_dC["n"],
                "accuracy check": [val.lower() == "true" for val in vsmow_dC["accuracy check"]]
            }
            return pd.DataFrame(data)
        else:
            print(f"Unsupported isotope type: {isotope}")
            return pd.DataFrame()
