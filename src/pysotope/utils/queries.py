import os
import datetime 
from .pick_chain_rt import *

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
        response = input("\nMethanol dD is -72.5 ± 3.1 ‰. Is this correct? (Y/N)\n").lower()
        if pos_response(response):
            meth_dD = -72.5
            meth_std = 3.1
            unknown = methyl_correction(unknown, stds)
            break
        elif neg_response(response):
            meth_dD = input("\nMeasured dD value of the methanol used in FAME methylation?\n")
            meth_std = input("\nUncertainty of the methanol δD value?\n")
            unknown = methyl_correction(unknown, stds, mdD=float(meth_dD), mdD_err=float(meth_std))
            break
        else:
            print("\nInvalid response. Try again.\n")
    append_to_log(log_file_path, f"Methanol dD: {meth_dD} ± {meth_std} ‰")
    return unknown, stds

def q_original_phthalic_value():
    o_ph_dD = input("Enter isotopic value of phthalic acid: ")
    return o_ph_dD

def q_output():
    o_fp = input("Provide a folder path for the output data:\n")
    return o_fp

def ask_user_for_rt(log_file_path, df, isotope):
    chain_lengths = ['C16', 'C18', 'C20', 'C22', 'C24', 'C26', 'C28', 'C30', 'C32']
    while True:
        response = input("Do you want to detect components in this dataset by retention time? (Y/N):\n").strip().lower()
        if pos_response(response):
            # picked = pick_chain_retention(isotope, df)
            append_to_log(log_file_path, "User opted to identify chains.")
            rt_values = input("Enter retention times for " + ", ".join(chain_lengths) + " separated by commas (type 'none' for any you don't want to use):\n")
            rt_values = rt_values.split(',')
            if len(rt_values) == len(chain_lengths):
                rt_dict = {chain: (None if rt.strip().lower() == 'none' else float(rt.strip())) for chain, rt in zip(chain_lengths, rt_values)}
                return rt_dict
            else:
                print("Invalid input. Please provide the correct number of values.\n")
        elif neg_response(response):
            append_to_log(log_file_path, "User opted not to identify chains.")
            print("Component detection not selected.\n")
            return None
        else:
            print("Invalid response. Please answer 'yes' or 'no'.\n")