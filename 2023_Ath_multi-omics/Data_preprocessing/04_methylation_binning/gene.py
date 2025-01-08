"""
Gene class definition to use for binning methylation.

Author: Serena G. Lotreck
"""
import pandas as pd
import numpy as np
import operator
import re


class Gene:

    def __init__(self, seqid, start, end, strand, name):
        """
        Instantiate an instance of the Gene class.

        parameters:
            seqid, str: ID from the first column of the gff file
            start, int: start bp
            end, int: end bp
            strand, str: '+' or '-'
            name, str: name of gene

        returns: None
        """
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name

        self.CG_pres_abs = None
        self.CHG_pres_abs = None
        self.CHH_pres_abs = None
        self.CG_prop = None
        self.CHG_prop = None
        self.CHH_prop = None
        self.CG_pres_abs_mean = None
        self.CHG_pres_abs_mean = None
        self.CHH_pres_abs_mean = None
        self.CG_prop_median = None
        self.CHG_prop_median = None
        self.CHH_prop_median = None


    def set_methylation_attr(self, methylation_df, attr_name):
        """
        Filter a methylation dataset by what sites are in this
        gene, or within 500bp up- or downstream, and add the resulting df to the
        attribute corresponding to attr_name.

        parameters:
            methylation_ddf, pandas df: dataframe with methylation data.
            attr_name, str: name of attribute to set.

        returns: None
        """
        # Filter by where the methylation sites are
        in_near_gene_df = methylation_df.loc[
                (self.start-500 < methylation_df.index.get_level_values("bp")) &
                (methylation_df.index.get_level_values("bp") < self.end+500) &
                (methylation_df.index.get_level_values("seqid") == self.seqid)]

        # Set attribute
        setattr(self, attr_name, in_near_gene_df)


    @staticmethod
    def bin_data_calculate_stat(df, start, end, methylation_type, name_part,
            num_bins):
        """
        Put data into base pair bins for a section of the gene (upstream,
        gene body, downstream).

        parameters:
            df, pandas df: the data to bin, accessions are columns, index is
                a 3-level multiindex of seqid, bp, strand
            start, int: what base pair is the start of the region to bin
            end, int: what abse pair is the end of the region to bin
            methylation_type, str: name of methylation type this is for. Name
                should correspond to one of the methylation dataset attributes.
            name_part, str: "pre", "gene_body" or "post".
            num_bins, int: number of bins to split the bp segment into

        returns:
            stat_df, df: accessions are rowns, bins are columns
        """
        # Get the bin boundaries using qcut on all bases in gene
        bin_intervals, bins = pd.qcut(pd.Series(np.arange(start, end, 1)),
                                    num_bins, retbins=True)

        # Get bin names to pass to cut
        bin_names = [f'{name_part}_bin_{i}'
                    for i in range(len(set(bin_intervals)))]

        # Then use bin boundaries with cut on the methylated bases
        bin_idx_values = pd.cut(df.index.get_level_values('bp'),
            bins=bins, labels=bin_names, include_lowest=True)
        bin_df = df.copy(deep=True)
        bin_df.index = bin_idx_values

        # Groupby to get stats
        if 'pres_abs' in methylation_type:
            stat_df = bin_df.groupby(by=bin_df.index).mean().dropna().T
            # T puts bins on columns, acs on rows
        elif 'prop' in methylation_type:
            stat_df = bin_df.groupby(by=bin_df.index).median().dropna().T

        # Add columns of zeros for bins with no C's in them
        num_accs = stat_df.shape[0]
        stat_df = pd.DataFrame({col_name:(np.zeros(num_accs, dtype=float) \
                        if col_name not in stat_df.columns.values.tolist() \
                        else stat_df[col_name].tolist()) for col_name in bin_names},
                        index=stat_df.index)

        # Sort columns
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)',s)]
        sorted_cols = sorted(stat_df.columns.values.tolist(), key=natsort)
        stat_df = stat_df.reindex(sorted_cols, axis=1)

        return stat_df


    def set_bin_stat(self, methylation_type, gb_bins, pp_bins):
        """
        Separate sites into bins and calculate statistic across bins and
        accessions for a given type of methylation, and assign the result
        as an attribute.

        NOTE: If there are no C's in a bin, the values will be set to 0.

        parameters:
            methylation_type, str: name of methylation type to do this for.
                Name should correspond to one of the methylation dataset
                attributes.
            gb_bins, int: number of bins to use for gene body
            pp_bins, int: number of bins to use for pre and post sequences

        returns: None
        """
        # Get the dataset
        df = getattr(self, methylation_type)

        # Separate into in-gene and before/after-gene sites
        gene_body_df = df.loc[(self.start < df.index.get_level_values('bp'))
                            & (df.index.get_level_values('bp') < self.end)]
        upstream_df = df.loc[(self.start > df.index.get_level_values('bp'))
                           & (self.start-500 < df.index.get_level_values('bp'))]
        downstream_df = df.loc[(df.index.get_level_values('bp') > self.end)
                           & (self.end+500 > df.index.get_level_values('bp'))]

        # Bin and calculate within gene bins
        gene_body_stat_df = Gene.bin_data_calculate_stat(gene_body_df, self.start,
                self.end, methylation_type, "gene_body", gb_bins)
        up_stat_df = Gene.bin_data_calculate_stat(upstream_df, self.start-500,
                self.start, methylation_type, "pre", pp_bins)
        down_stat_df = Gene.bin_data_calculate_stat(downstream_df, self.end,
                self.end+500, methylation_type, "post", pp_bins)

        # Combine all df's
        overall_stat_df = pd.concat([up_stat_df, gene_body_stat_df, down_stat_df], axis=1)

        # Set as attribute
        if 'pres_abs' in methylation_type:
            setattr(self, f'{methylation_type}_mean', overall_stat_df)
        elif 'prop' in methylation_type:
            setattr(self, f'{methylation_type}_median', overall_stat_df)


    def reverse_datasets(self):
        """
        Reverses the order of the bins in the summary df's if the gene is on
        the - strand.

        parameters: self

        returns: None
        """
        # Define the attributes to change
        stat_dfs = ['CG_pres_abs_mean', 'CG_prop_median', 'CHG_pres_abs_mean',
                    'CHG_prop_median', 'CHH_pres_abs_mean', 'CHH_prop_median']

        # Check if the gene needs to be reversed
        if self.strand == '-':
            for attr in stat_dfs:

                # Get the attribute df
                df = getattr(self, attr)

                # Check if the attribute is missing
                # Allows test to pass without having to assign all 6, warning
                # is to avoid potential silent errors in the future
                if df is None:
                    print(f'Warning! Attribute {attr} is missing')
                else:
                    # Rename the columns in reverse order
                    df_new_cols = df.columns.values.tolist()
                    df_new_cols.reverse()
                    df = df.set_axis(df_new_cols, axis=1)
                    # Reverse so they're in natsorted order again
                    df = df[df.columns[::-1]]
                    # Reassign attribute
                    setattr(self, attr, df)

