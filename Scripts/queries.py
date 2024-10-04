from correction.drift import *
from correction.linearity import *
from correction.methanol import *
from correction.vsmow import *
import os
from base_functions import append_to_log
from queries import *


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


def q_drift(log_file_path):
    """
    Ask the user for their choice of correction and keep asking until a valid response is given.
    """
    while True:  # Start an infinite loop
        choice = input("\nApply drift to correction? (Y/N)\n")

        if choice in ["yes", "y", "true", "t", "no", "n", "false", "f"]:
            return choice
            append_to_log(log_file_path, f"Drift application: {choice}")
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
