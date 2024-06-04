import os
from pathlib import Path
import requests
import pandas as pd

import janitor

from . import log

logger = log.get_logger(__file__)


def set_data_dir():
    # Path of the current script
    current_script_path = Path(__file__).resolve()
    # Path to the top-level data directory
    res = current_script_path.parent.parent / "data"
    logger.debug(f"data dir set to {res}")
    return res


data_dir = set_data_dir()


# print(data_dir)

FILE_CHOCIES = [
    "Phosphorylation_site_dataset",
    "Regulatory_sites",
]


def check_files_exist(file_names, directory=data_dir):
    """Check if all required files exist in the specified directory."""
    missing_files = []
    for file_name in file_names:
        if not os.path.exists(os.path.join(directory, file_name)):
            missing_files.append(file_name)
    if missing_files:
        print(f"Missing files: {missing_files}")
        return False
    return True


def load_data(file_name, directory="data"):
    """Load data from a specified file within the given directory."""
    file_path = os.path.join(directory, file_name)
    if os.path.exists(file_path):
        # Example of loading a CSV file, adjust the loading method as needed
        return pd.read_csv(file_path)
    else:
        print(f"File not found: {file_name}")
        return None


def get_psiteplus_file(file, **kwargs):

    skiprows = 3
    if "skiprows" in kwargs:
        skiprows = kwargs.pip("skiprows")

    fullfile = data_dir / file
    df = pd.read_table(fullfile, skiprows=skiprows, **kwargs)
    df = janitor.clean_names(df)
    return df


# # Example Usage
# required_files = ['phosphosite_data.csv', 'another_data_file.csv']
# if check_files_exist(required_files):
#     for file_name in required_files:
#         data = load_data(file_name)
#         # Process data as needed
# else:
#     print("Please ensure all required files are in the data directory.")
