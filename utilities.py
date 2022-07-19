import csv
import os

import pandas as pd


def dict_to_csv(dict, file_name):
    df = pd.DataFrame(dict)
    df.to_csv(file_name)
    return df


def df_to_csv(df, file_name):
    df.to_csv(file_name)
    return df


def create_dir(dir_name):
    results_dir = dir_name
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return dir_name

def write_rows_to_csv(file_path, header, data_rows):
    with open(file_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data_rows)