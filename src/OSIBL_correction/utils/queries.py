import os
import datetime 
from . .definitions.standards import *

def append_to_log(log_file_path, log_message):
    """
    Add entry to log file.
    """
    with open(log_file_path, "a") as log_file:
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        initial_message = f"Log file created at {current_datetime}\n"
        log_file.write(log_message + "; " + str(current_datetime) + "\n")


def pos_response(response):
    return response.lower() in {"yes", "y", "true", "t", ""}


def neg_response(response):
    return response.lower() in {"no", "n", "false", "f"}


def query_project_name():
    project_name = input("\nProvide the project name.\n")
    return project_name 


def query_file_location():
    while True:
        loc = input("\nProvide the full path of the GC-IRMS datafile (as .csv).\n")
        
        # Remove single quotes if present at both the start and end of the input
        if loc.startswith("'") and loc.endswith("'"):
            loc = loc[1:-1]
        
        if os.path.isfile(loc) and loc.endswith(".csv"):
            return loc
        else:
            print("\nFile does not exist or is not a .csv file. Try again.\n")

def isotope_type():
    """
    Ask the user for their choice isotope.
    """
    while True:  # Start an infinite loop
        choice = input("\nCorrecting carbon or hydrogen? (c/h)\n")

        if choice in ["c", "h"]:
            if choice == "c":
                isotope = "dC"
            else:
                isotope = "dD"
            return isotope
        else:
            print("\nInvalid response\n")

def query_stds(alt_stds):
    
    "Ask user for chain lengths in standards"
    if alt_stds:
        while True:
            usr_lin_std = input("Which chain lengths are included in the linearity standards? (comma-separated, e.g., C20, C28): ")
            # Split the input and check the length
            linearity_chain_lengths = [length.strip() for length in usr_lin_std.split(',')]
            if len(linearity_chain_lengths) == 2:
                break  # Exit loop if exactly two chain lengths are provided
            else:
                print("Error: You must provide exactly two chain lengths. Please try again.")
        while True:
            usr_drift_std = input("Which chain lengths are included in the drift standards? (comma-separated, e.g., C18, C24): ")
            # Split the input and check the length
            drift_chain_lengths = [length.strip() for length in usr_drift_std.split(',')]
            if len(drift_chain_lengths) == 2:
                break  # Exit loop if exactly two chain lengths are provided
            else:
                print("Error: You must provide exactly two chain lengths. Please try again.")
    else:
        linearity_chain_lengths = default_stds["linearity_stds"]
        drift_chain_lengths = default_stds["drift_stds"]
    return linearity_chain_lengths, drift_chain_lengths

def lin_response(log_file_path):
    valid_responses = ["yes", "y", "true", "t", "no", "n", "false", "f"]
    while True:
        response = input("\nAssign a linearity correction? (Y/N)\n").lower()
        if response in valid_responses:
            append_to_log(log_file_path, "Linearity application application: " + str(response))
            return response
        else:
            print("\nInvalid response. Try again.\n")


def q_methylation(unknown, stds, log_file_path):  # , user_choice, response):
    from .corrections.methanol import methyl_correction
    while True:
        response = input("\nMethanol δD is -72.5 ± 3.1 ‰. Is this correct? (Y/N)\n").lower()
        if pos_response(response):
            meth_dD = -72.5
            meth_std = 3.1
            unknown = methyl_correction(unknown, stds)
            break
        elif neg_response(response):
            meth_dD = input("\nMeasured δD value of the methanol used in FAME methylation?\n")
            meth_std = input("\nUncertainty of the methanol δD value?\n")
            unknown, stds = methyl_correction(unknown, stds, mdD=meth_dD, mdD_err=meth_std)
            break
        else:
            print("\nInvalid response. Try again.\n")
    append_to_log(log_file_path, f"Methanol δD: {meth_dD} ± {meth_std} ‰")
    return unknown, stds

def q_original_phthalic_value():
    o_ph_dD = input("Enter isotopic value of phthalic acid: ")
    return o_ph_dD

def q_output():
    o_fp = input("Provide a folder path for the output data:\n")
    return o_fp