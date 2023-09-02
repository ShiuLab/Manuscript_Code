"""
Calculate bin statistics for methylation data and save as a pickled object.

Author: Serena G. Lotreck
"""
from os.path import abspath
import argparse

from gene import Gene

import datetime
import pandas as pd
import datatable as dt
import pickle


def make_new_cols(df):
    """
    Make mapper for renaming bins with abbreviated names.

    parameters:
        df, pandas df: any one of the 6 methylation stat attribute df's

    returns:
        new_cols, dict: keys are old col names and values are corresponding new
            col names
    """
    new_cols = {}

    for col in df.columns.values.tolist():

        bin_num = col.split('_')[-1]

        if 'pre' in col:
            new_col = f'R_{bin_num}'
        elif 'gene' in col:
            new_col = f'B_{bin_num}'
        elif 'post' in col:
            new_col = f'P_{bin_num}'

        new_cols[col] = new_col

    return new_cols


def bin_and_make_feature_matrices(genes, methylation_datasets, gbb, ppb):
    """
    Make feature matrices for all 6 datasets.

    parameters:
        genes, generator of Gene obj: genes to use to make feature matrices
        methylation_datasets, dict: keys are methylation attribute names,
            values are methylation df's
        gbb, int: number of bins for gene body
        ppb, int: number of bins for pre and post sequences

    returns:
        data, dict: keys are attribute names, values are full feature matrices
    """
    # Initialize variables for feature tables
    print('Initializing loop variables...')
    data = {
        'CG_pres_abs_mean': None,
        'CG_prop_median': None,
        'CHG_pres_abs_mean': None,
        'CHG_prop_median': None,
        'CHH_pres_abs_mean': None,
        'CHH_prop_median': None
    }

    # Loop over genes and add to tables
    print('Looping over genes to build feature matrices...')
    for i, gene in enumerate(genes):

        print(f'Working on gene {i} of {len(genes)}')

        # Do binning 
        for dset_name, dset in methylation_datasets.items():
            gene.set_methylation_attr(dset, dset_name)
            gene.set_bin_stat(dset_name, gbb, ppb)

        # Check and reverse strands if necessary
        gene.reverse_datasets()

        # Make a dict for renaming columns
        new_cols = make_new_cols(gene.CG_pres_abs_mean)

        # If it's the first one, make it the df
        if i == 0:

            for attr_name, value in data.items():

                # Get the df
                df = getattr(gene, attr_name)

                # Rename index 
                df.index.names = ['accession'] ## TODO fix the fact that this
                ## results in a "Unnamed: 0 column" instead of one named
                ## 'accession'

                # Rename cols with shorthand names
                df = df.rename(columns=new_cols)

                # Add prefix
                df = df.add_prefix(f'{gene.name}_{attr_name}_')

                # Assign back to dict
                data[attr_name] = df

        # Otherwise, append
        else:

            for attr_name, value in data.items():
                # Get the new addition
                df = getattr(gene, attr_name)

                # Rename cols with shorthand names
                df = df.rename(columns=new_cols)

                # Add prefix
                df = df.add_prefix(f'{gene.name}_{attr_name}_')

                # Concat onto existing df
                df = pd.concat((data[attr_name], df), axis=1)

                # Put back in dict
                data[attr_name] = df

    return data


def read_methylation_datasets(methylation_datasets):
    """
    Read in methylation datasets from csv's.

    parameters:
        methylation_datasets, dict: keys are methylation attribute
            names, values are paths to the datasets

    returns:
        meethylation_datasets, dict: keys are methylation attribute
            names, values are dataframes with data
    """
    for dset_name, dset in methylation_datasets.items():

        print(f'Reading in {dset_name}...')

        # Read in dataset with datatable and make into a pandas df
        dset_df = dt.fread(dset).to_pandas().set_index(
            'C0').T
        print(dset_df.head())

        # Make multiindex 
        print('Making multiindex...')
        chr_idx = [idx[0] for idx in dset_df.index.str.split('_')]
        bp_idx = [int(idx[1]) for idx in dset_df.index.str.split('_')]
        strand_idx = [idx[3] for idx in dset_df.index.str.split('_')]
        dset_df.index = pd.MultiIndex.from_arrays(
            [chr_idx, bp_idx, strand_idx], names=('seqid', 'bp', 'strand'))

        # Assign back to dict
        methylation_datasets[dset_name] = dset_df

    return methylation_datasets


def get_start_row(gff):
    """
    Determine at what row the data starts in a gff file.

    parameters:
        gff, str: path to gff file with genomic features

    returns:
        start_row, int: the index of the row where the data starts,
            which is equivalent to the number of rows to skip
    """
    start_row = 0
    with open(gff) as f:
        line = f.readline()
        line_num = 0  # Keep track of what line we're on
        while start_row == 0:  # Only read lines until we find the start row
            if line[0] != '#':  # If it doesn't start with a comment, it's the first data row
                if line_num == 0:  # Account for the case where the first line is the first data line
                    start_row = line_num
                    break
                else:
                    start_row = line_num
            else:
                start_row += 1

    return start_row


def get_genes(gff):
    """
    Create Gene class instances for genomic features with the type "gene".

    parameters:
        gff, str: path to gff file with genomic features.

    returns:
        genes, list: list of gene objects
    """
    # Determine at what row the data starts
    print('Getting start row...')
    start_row = get_start_row(gff)

    # Read in the data
    print('Reading in genome data...')
    gff_df = pd.read_csv(gff,
                         sep='\t',
                         skiprows=start_row,
                         names=[
                             'seqid', 'source', 'type', 'start', 'end',
                             'score', 'strand', 'phase', 'attributes'
                         ])

    # Filter df to keep only 'gene' rows
    print('Getting genes...')
    genes_df = gff_df[gff_df['type'] == 'gene']

    # Make gene objects
    print('Making gene objects...')
    genes = []
    for index, row in genes_df.iterrows():
        gene_name = row.at['attributes'][row.at['attributes'].index('Name=') +
                                         5:]
        gene = Gene(row.at['seqid'], int(row.at['start']), int(row.at['end']),
                    row.at['strand'], gene_name)
        genes.append(gene)

    return genes


def main(gff, gbb, ppb, CG_pres_abs, CHG_pres_abs, CHH_pres_abs, CG_prop,
         CHG_prop, CHH_prop, out_loc, file_prefix):

    begin_time = datetime.datetime.now()

    # Make gene objects
    print('\nRetrieving  genes from genome...')
    genes = get_genes(gff)

    # Read in methylation datasets
    print('\nReading in methylation datasets...')
    methylation_datasets = {
        'CG_pres_abs': CG_pres_abs,
        'CHG_pres_abs': CHG_pres_abs,
        'CHH_pres_abs': CHH_pres_abs,
        'CG_prop': CG_prop,
        'CHG_prop': CHG_prop,
        'CHH_prop': CHH_prop
    }
    methylation_datasets = read_methylation_datasets(methylation_datasets)

    # Assign methylation attributes, bin and make feature tables
    print('\nBinning data and making feature tables...')
    feature_matrices = bin_and_make_feature_matrices(genes,
            methylation_datasets, gbb, ppb)

    # Save feature matrices
    print('\nSaving feature matrices...')
    for name, df in feature_matrices.items():
        df.to_csv(f'{out_loc}/{file_prefix}_{name}_feature_matrix.csv')

    end_time = datetime.datetime.now()

    print(f'\nDone columnwise in {end_time - begin_time} seconds')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bin methylated sites')

    parser.add_argument(
        '-gff',
        type=str,
        help='Path to gff3-formatted file with genomic features.')
    parser.add_argument(
        '-gbb',
        type=int,
        default=20,
        help='Number of bins to use for the gene body. Default is 20.')
    parser.add_argument(
        '-ppb',
        type=int,
        default=5,
        help='Number of bins to use for pre and post seqs. Default is 5.')
    parser.add_argument(
        '-CG_pres_abs',
        type=str,
        help='Path to file with presence/absence data for CG methylation.')
    parser.add_argument(
        '-CHG_pres_abs',
        type=str,
        help='Path to file with presence/absence data for CHG methylation.')
    parser.add_argument(
        '-CHH_pres_abs',
        type=str,
        help='Path to file with presence/absence data for CHH methylation.')
    parser.add_argument(
        '-CG_prop',
        type=str,
        help='Path to file with proportion data for CG methylation.')
    parser.add_argument(
        '-CHG_prop',
        type=str,
        help='Path to file with proportion data for CHG methylation.')
    parser.add_argument(
        '-CHH_prop',
        type=str,
        help='Path to file with proportion data for CHH methylation.')
    parser.add_argument('-out_loc',
                        type=str,
                        help='Path to directory to save the output.')
    parser.add_argument('-file_prefix',
                        type=str,
                        help='String to prepend to output files.')

    args = parser.parse_args()

    args.gff = abspath(args.gff)
    args.CG_pres_abs = abspath(args.CG_pres_abs)
    args.CHG_pres_abs = abspath(args.CHG_pres_abs)
    args.CHH_pres_abs = abspath(args.CHH_pres_abs)
    args.CG_prop = abspath(args.CG_prop)
    args.CHG_prop = abspath(args.CHG_prop)
    args.CHH_prop = abspath(args.CHH_prop)
    args.out_loc = abspath(args.out_loc)

    main(**vars(args))
