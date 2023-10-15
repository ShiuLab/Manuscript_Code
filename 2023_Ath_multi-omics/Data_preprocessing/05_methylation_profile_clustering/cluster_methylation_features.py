"""
Make a feature table based on clustering of methylation profiles.

Author: Serena G. Lotreck
"""
import argparse
from os.path import abspath, basename, splitext
from collections import defaultdict
import itertools

from tqdm import tqdm
import datatable as dt
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder


def make_feature_table(encoded, row_ids, new_feat_names, data_type):
    """
    Return the clustered data to the original feature table format, but with
    clusters as features now instead of bins. Final features are of the form
        {gene name}_{data type}_cluster_{cluster number}

    parameters:
        encoded, 2D arr: output of one-hot encoding of the clustered instances
        row_ids, pandas Index: original row labels of the tidy-form data
        new_feat_names, 1D arr: names of the one-hot-encoded clusters
        data_type, str: identifier for this table's data type

    returns:
        feature_table, pandas df: Final feature table
    """
    # Make into pandas df with proper labels
    df = pd.DataFrame(encoded, index=row_ids, columns=new_feat_names)

    # Get rows corresponding to genes
    print('Matching genes to row names...')
    row_names = df.index.values.tolist()
    gene_names = lambda x: x.split('_')[-1]
    row_lists = defaultdict(list)
    for gene_name, row_names in itertools.groupby(row_names, key=gene_names):
        row_lists[gene_name] += list(row_names)

    # Make df
    print('Reorganizing data...')
    for gene_name, row_names in tqdm(row_lists.items()):

        # Make a mapper for the rows
        accesion_names = lambda x: x.split('_')[0]
        row_mapper = {name:accesion_names(name) for name in row_names}

        # Check if this is the first one
        try:
            feature_table
            # If it's not the first, one, subset, rename idxs and concat
            add_to_feat = df.loc[row_names]
            add_to_feat = add_to_feat.rename(index=row_mapper)
            add_to_feat = add_to_feat.add_prefix(f'{gene_name}_{data_type}_')
            feature_table = pd.concat([feature_table, add_to_feat], axis=1)
        except UnboundLocalError:
            # If it is the first, subset and make that into tidy_df, rename
            feature_table = df.loc[row_names]
            feature_table = feature_table.rename(index=row_mapper)
            feature_table = feature_table.add_prefix(f'{gene_name}_{data_type}_')

    return feature_table


def make_tidy_data(df):
    """
    Make a matrix of shape (n_accessions*n_genes, n_bins), where each row
    is {instance (accession) ID}_{gene #}, and the columns are the bins that
    methylation data has been placed into.

    The columns of the input dataset must have a very specific format:
        {gene name}_{methylation & data type}_{R/B/P}_{bin number}
    for example:
        AT1G01010_CG_pres_abs_mean_R_0
    There can be underscores in the methylation and data type, but there
    can NOT be underscores in the gene name, and the bin ID (R/B/P) and number
    must be the last two items when the string is split on underscores.

    parameters:
        df, pandas df: rows are instances and columns are the 30 bins for each
            of 27,000 genes

    returns:
        tidy_df, pandas df: the tidy version of the data
    """
    # Get the lists of column names that correpsond to each gene
    print('Matching columns names to their genes...')
    cols = df.columns.values.tolist()
    gene_names = lambda x: x.split('_')[0] # Lambda func for groupby
    col_lists = defaultdict(list)
    for gene_name, col_names in itertools.groupby(cols, key=gene_names):
        col_lists[gene_name] += list(col_names)

    # Make df
    print('Reorganizing data...')
    for gene_name, col_names in tqdm(col_lists.items()):

        # Make a mapper for the columns
        bin_names = lambda x: '_'.join(x.split('_')[-2:])
        column_mapper = {name:bin_names(name) for name in col_names}

        # Check if this is the first one
        try:
            tidy_df
            # If it's not the first, one, subset, rename idxs and concat
            add_to_tidy = df[col_names]
            add_to_tidy = add_to_tidy.rename(columns=column_mapper)
            add_to_tidy = add_to_tidy.set_index(
                    add_to_tidy.index.astype(str) + f'_{gene_name}')
            tidy_df = pd.concat([tidy_df, add_to_tidy])
        except UnboundLocalError:
            # If it is the first, subset and make that into tidy_df, rename
            tidy_df = df[col_names]
            tidy_df = tidy_df.rename(columns=column_mapper)
            tidy_df = tidy_df.set_index(tidy_df.index.astype(str) + f'_{gene_name}')

    return tidy_df


def main(feature_table_path, data_type, num_clusters, out_loc):

    # Read in data
    print('\nReading in feature table...\n')
    df = dt.fread(feature_table_path).to_pandas().set_index('C0')
    df.index.name = None
    print(df.head())

    # Preprocess data
    print('\nMaking tidy data...')
    df = make_tidy_data(df)
    row_ids = df.index
    X = df.to_numpy()

    # Cluster
    print('\nClustering...')
    clusters = KMeans(n_clusters=num_clusters, random_state=42).fit_predict(X)

    # One hot encode
    print('\nOne hot encoding...')
    encoder = OneHotEncoder()
    encoded = encoder.fit_transform(clusters.reshape((len(row_ids), 1))).toarray()

    # Make feature table
    print('\nMaking feature table...')
    new_feat_names = encoder.get_feature_names_out(input_features=['cluster'])
    print(f'Feature names: '
        '{new_feat_names[:5] if len(new_feat_names) > 5 else new_feat_names}')
    feature_table = make_feature_table(encoded, row_ids, new_feat_names,
            data_type)
    print(f'Snapshot of feature table:\n{feature_table.head()}')

    # Write out
    print('\nWriting out file...')
    file_name = splitext(basename(feature_table_path))[0]
    print(f'{out_loc}/{file_name}_{num_clusters}_clusters.csv')
    feature_table.to_csv(f'{out_loc}/{file_name}_{num_clusters}_clusters.csv')

    print('\nDone!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Cluster methylation data')

    parser.add_argument('feature_table', type=str,
            help='Path to feature table to cluster')
    parser.add_argument('data_type', type=str,
            help='Descriptor for the data type in the table, e.g. '
            'CG_pres_abs_mean')
    parser.add_argument('num_clusters', type=int,
            help='Number of clusters to use for kmeans')
    parser.add_argument('out_loc', type=str,
            help='Path to save files')

    args = parser.parse_args()

    args.feature_table = abspath(args.feature_table)
    args.out_loc = abspath(args.out_loc)

    main(args.feature_table, args.data_type, args.num_clusters, args.out_loc)
