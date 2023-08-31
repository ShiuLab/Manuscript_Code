"""
Spot checks for the Gene class.

Author: Serena G. Lotreck
"""
import unittest
import sys
import os
import shutil
import io

sys.path.append('../data_preprocessing/methylation_binning/')
from gene import Gene

import pandas as pd
from pandas.testing import assert_frame_equal


class TestGene(unittest.TestCase):

    def setUp(self):

        # Create test directory
        self.tmpdir = "tmp"
        os.makedirs(self.tmpdir, exist_ok=True)

        # Create test files
        CG_pres_abs = io.StringIO(',Chr1_23_CG_+,Chr1_553_CG_+,Chr1_872_CG_-,Chr1_1230_CG_+,'
                'Chr2_720_CG_+\nac1,1,0,1,1,0\nac2,'
                '0,0,1,0,1')
        self.CG_pres_abs = pd.read_csv(CG_pres_abs, index_col=0).T
        chr_idx = [idx[0] for idx in self.CG_pres_abs.index.str.split('_')]
        bp_idx = [int(idx[1]) for idx in self.CG_pres_abs.index.str.split('_')]
        strand_idx = [idx[3] for idx in self.CG_pres_abs.index.str.split('_')]
        self.CG_pres_abs.index = pd.MultiIndex.from_arrays([chr_idx, bp_idx, strand_idx],
                                                             names=('seqid','bp','strand'))

        CG_prop = io.StringIO(',Chr1_23_CG_+,Chr1_553_CG_+,Chr1_872_CG_-,Chr1_1230_CG_+,'
                'Chr2_720_CG_+\nac1,0.23,0.45,0.21,0.97,0.11\nac2,'
                '0.65,0.23,0.67,0.19,0.87')
        self.CG_prop = pd.read_csv(CG_prop, index_col=0).T
        chr_idx = [idx[0] for idx in self.CG_prop.index.str.split('_')]
        bp_idx = [int(idx[1]) for idx in self.CG_prop.index.str.split('_')]
        strand_idx = [idx[3] for idx in self.CG_prop.index.str.split('_')]
        self.CG_prop.index = pd.MultiIndex.from_arrays([chr_idx, bp_idx, strand_idx],
                                                             names=('seqid','bp','strand'))

        # Create the Gene instance
        self.gene = Gene('Chr1', 560, 890, '+', 'AT102947')
        self.gene_negative = Gene('Chr1', 560, 890, '-', 'AT345987')


    def tearDown(self):

        shutil.rmtree(self.tmpdir)


    def test_set_methylation_attr(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')

        right_answer = pd.DataFrame({'ac1':[0,1,1], 'ac2':[0,1,0]},
                index=pd.MultiIndex.from_arrays([['Chr1','Chr1','Chr1'],
                                                 [553, 872, 1230],
                                                 ['+', '-', '+']],
                                                names=('seqid', 'bp', 'strand')))

        assert_frame_equal(self.gene.CG_pres_abs, right_answer)


    def test_bin_data_calculate_stat_pre_mean(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')

        pre_df = self.gene.CG_pres_abs.loc[(self.gene.start >
            self.gene.CG_pres_abs.index.get_level_values('bp')) &
            (self.gene.start-500 <
            self.gene.CG_pres_abs.index.get_level_values('bp'))]

        pre_stat_df_mean = self.gene.bin_data_calculate_stat(pre_df,
                self.gene.start-500, self.gene.start, 'CG_pres_abs', 'pre', 5)

        right_answer = pd.DataFrame({f'pre_bin_{i}':[0.0, 0.0] for i in range(5)},
                index=['ac1','ac2'])

        # Don't check dtype because if mean or median functions are used on only
        # one item, if that item is an int, the result is an int
        assert_frame_equal(pre_stat_df_mean, right_answer, check_dtype=False)


    def test_bin_data_calculate_stat_gene_mean(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')

        gene_body_df = self.gene.CG_pres_abs.loc[(self.gene.start <
            self.gene.CG_pres_abs.index.get_level_values('bp')) &
            (self.gene.CG_pres_abs.index.get_level_values('bp') <
            self.gene.end)]

        gene_stat_df_mean = self.gene.bin_data_calculate_stat(gene_body_df,
                self.gene.start, self.gene.end, 'CG_pres_abs', 'gene_body', 20)


        right_ans_dicts = {**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                        **{'gene_body_bin_18':[1.0,1.0]},
                        **{'gene_body_bin_19':[0.0,0.0]}}
        right_answer = pd.DataFrame(right_ans_dicts, index=['ac1','ac2'])

        assert_frame_equal(gene_stat_df_mean, right_answer, check_dtype=False)


    def test_bin_data_calculate_stat_post_median(self):

        self.gene.set_methylation_attr(self.CG_prop, 'CG_prop')

        post_df = self.gene.CG_prop.loc[(self.gene.end <
            self.gene.CG_prop.index.get_level_values('bp'))
            & (self.gene.end+500 > self.gene.CG_prop.index.get_level_values('bp'))]

        post_stat_df_median = self.gene.bin_data_calculate_stat(post_df,
                self.gene.end, self.gene.end+500, 'CG_prop', 'post', 5)


        right_ans_dict = {**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                          **{'post_bin_3':[0.97, 0.19]},
                          **{'post_bin_4':[0.0,0.0]}}
        right_answer = pd.DataFrame(right_ans_dict, index=['ac1','ac2'])

        assert_frame_equal(post_stat_df_median, right_answer, check_dtype=False)


    def test_set_bin_stat_mean(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')

        self.gene.set_bin_stat('CG_pres_abs', 20, 5)

        right_ans_dict = {**{f'pre_bin_{i}':[0.0, 0.0] for i in range(5)},
                          **{**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                             **{'gene_body_bin_18':[1.0,1.0]},
                             **{'gene_body_bin_19':[0.0,0.0]}},
                          **{**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                             **{'post_bin_3':[1.0, 0.0]},
                             **{'post_bin_4':[0.0,0.0]}}}
        right_answer = pd.DataFrame(right_ans_dict, index=['ac1', 'ac2'])

        assert_frame_equal(self.gene.CG_pres_abs_mean, right_answer)


    def test_bin_stat_median(self):

        self.gene.set_methylation_attr(self.CG_prop, 'CG_prop')

        self.gene.set_bin_stat('CG_prop', 20, 5)

        right_ans_dict = {**{**{f'pre_bin_{i}':[0.0, 0.0] for i in range(4)},
                             **{'pre_bin_4':[0.45, 0.23]}},
                          **{**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                             **{'gene_body_bin_18':[0.21,0.67]},
                             **{'gene_body_bin_19':[0.0,0.0]}},
                          **{**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                             **{'post_bin_3':[0.97, 0.19]},
                             **{'post_bin_4':[0.0,0.0]}}}
        right_answer = pd.DataFrame(right_ans_dict, index=['ac1', 'ac2'])

        assert_frame_equal(self.gene.CG_prop_median, right_answer)

    def test_reverse_datasets_pos(self):

        self.gene.set_methylation_attr(self.CG_pres_abs, 'CG_pres_abs')
        self.gene.set_methylation_attr(self.CG_prop, 'CG_prop')

        self.gene.set_bin_stat('CG_pres_abs', 20, 5)
        self.gene.set_bin_stat('CG_prop', 20, 5)

        self.gene.reverse_datasets()

        pa_ans_dict = {**{f'pre_bin_{i}':[0.0, 0.0] for i in range(5)},
                          **{**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                             **{'gene_body_bin_18':[1.0,1.0]},
                             **{'gene_body_bin_19':[0.0,0.0]}},
                          **{**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                             **{'post_bin_3':[1.0, 0.0]},
                             **{'post_bin_4':[0.0,0.0]}}}
        pa_right_answer = pd.DataFrame(pa_ans_dict, index=['ac1', 'ac2'])

        pr_ans_dict = {**{**{f'pre_bin_{i}':[0.0, 0.0] for i in range(4)},
                             **{'pre_bin_4':[0.45, 0.23]}},
                          **{**{f'gene_body_bin_{i}':[0.0,0.0] for i in range(18)},
                             **{'gene_body_bin_18':[0.21,0.67]},
                             **{'gene_body_bin_19':[0.0,0.0]}},
                          **{**{f'post_bin_{i}':[0.0,0.0] for i in range(3)},
                             **{'post_bin_3':[0.97, 0.19]},
                             **{'post_bin_4':[0.0,0.0]}}}
        pr_right_answer = pd.DataFrame(pr_ans_dict, index=['ac1', 'ac2'])

        assert_frame_equal(pa_right_answer, self.gene.CG_pres_abs_mean)
        assert_frame_equal(pr_right_answer, self.gene.CG_prop_median)


    def test_reverse_dataset_neg(self):

        self.gene_negative.set_methylation_attr(self.CG_pres_abs,
                'CG_pres_abs')
        self.gene_negative.set_methylation_attr(self.CG_prop, 'CG_prop')

        self.gene_negative.set_bin_stat('CG_pres_abs', 20, 5)
        self.gene_negative.set_bin_stat('CG_prop', 20, 5)

        self.gene_negative.reverse_datasets()

        pa_ans_dict = {**{**{'pre_bin_0':[0.0, 0.0]},
                          **{'pre_bin_1':[1.0, 0.0]},
                          **{f'pre_bin_{i}':[0.0, 0.0] for i in range(2,5,1)}},
                       **{**{'gene_body_bin_0':[0.0,0.0]},
                          **{'gene_body_bin_1':[1.0, 1.0]},
                          **{f'gene_body_bin_{i}':[0.0,0.0] for i in
                              range(2,20,1)},
                       **{f'post_bin_{i}':[0.0,0.0] for i in range(5)}}}

        pa_right_answer = pd.DataFrame(pa_ans_dict, index=['ac1', 'ac2'])

        pr_ans_dict = {**{**{'pre_bin_0':[0.0, 0.0]},
                          **{'pre_bin_1':[0.97, 0.19]},
                          **{f'pre_bin_{i}':[0.0, 0.0] for i in range(2,5,1)}},
                          **{**{'gene_body_bin_0':[0.0,0.0]},
                             **{'gene_body_bin_1':[0.21,0.67]},
                             **{f'gene_body_bin_{i}':[0.0,0.0]
                                 for i in range(2,20,1)}},
                          **{**{'post_bin_0':[0.45, 0.23]},
                              **{f'post_bin_{i}':[0.0, 0.0] for i in
                                  range(1,5,1)}}}
        pr_right_answer = pd.DataFrame(pr_ans_dict, index=['ac1', 'ac2'])

        assert_frame_equal(pa_right_answer, self.gene_negative.CG_pres_abs_mean)
        assert_frame_equal(pr_right_answer, self.gene_negative.CG_prop_median)


if __name__ == '__main__':
    unittest.main()
