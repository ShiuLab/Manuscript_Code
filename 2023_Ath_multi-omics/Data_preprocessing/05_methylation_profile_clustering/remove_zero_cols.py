"""
Remove columns with only 0's from clustering output. Makes copies rather
than overwriting the original files.

Quick and dirty add-on to cluster_methylation_features.py

Author: Serena G. Lotreck
"""
import argparse
from os import listdir
from os.path import abspath

import datatable as dt
import pandas as pd


def delete_zero_cols(f, target_dir, out_dir):
    """
    Saves a new file without any columns containing only zeroes.

    parameters:
        f, str: bas name for target file
        target_dir, str: full path to directory containing the file
        out_dir, str: path to save the file

    returns: None
    """
    # Read in file
    df = dt.fread(f'{target_dir}/{f}').to_pandas().set_index('C0')

    # Drop zero columns
    df = df[df.columns[df.sum()!=0]]

    # Save as new file with "no_zero_cols_" prefix
    df.to_csv(f'{out_dir}/no_zero_cols_{f}')


def main(target_dir, output_dir):

    files_to_modify = [f for f in listdir(target_dir) if '.csv' in f]

    for f in files_to_modify:

        print(f'\nModifying file {f}...')
        delete_zero_cols(f, target_dir, output_dir)

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Remove cols with only 0s')

    parser.add_argument('target_dir', type=str,
            help='Directory containing files to edit')
    parser.add_argument('out_dir', type=str,
            help='Directory to save new files')

    args = parser.parse_args()

    args.target_dir = abspath(args.target_dir)
    args.output_dir = abspath(args.target_dir)

    main(args.target_dir, args.output_dir)
