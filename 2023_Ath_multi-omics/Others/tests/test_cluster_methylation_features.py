"""
Spot checks for cluster_methylation_features.py.

Author: Serena G. Lotreck
"""
import unittest
import os
import shutil
import sys

sys.path.append('../data_preprocessing/methylation_profile_clustering')

from cluster_methylation_features import make_tidy_data, make_feature_table
import cluster_methylation_features as cmf
import pandas as pd
from pandas.testing import assert_frame_equal


class TestMakeTidyData(unittest.TestCase):


    def setUp(self):

        self.df = pd.DataFrame({'g1_my_data_type_R_0':[0.5, 1],
                                'g1_my_data_type_R_1':[1, 0.5],
                                'g2_my_data_type_R_0':[0, 0],
                                'g2_my_data_type_R_1':[1, 1]},
                            index=[101, 102])

        self.right_answer = pd.DataFrame({'R_0':[0.5, 1, 0, 0],
                                          'R_1':[1, 0.5, 1, 1]},
                                        index=['101_g1',
                                               '102_g1',
                                               '101_g2',
                                               '102_g2'])


    def test_make_tidy_data(self):

        answer = make_tidy_data(self.df)

        assert_frame_equal(answer, self.right_answer)


class TestMakeFeatureTable(unittest.TestCase):


    def setUp(self):

        self.clusters = [[0, 1],
                         [1, 1],
                         [0, 0],
                         [1, 0]]
        self.row_ids = ['101_g1',
                        '102_g1',
                        '101_g2',
                        '102_g2']
        self.new_feat_names = ['cluster_0', 'cluster_1']
        self.data_type = 'my_data_type'

        self.right_answer = pd.DataFrame({'g1_my_data_type_cluster_0':[0, 1],
                                          'g1_my_data_type_cluster_1':[1, 1],
                                          'g2_my_data_type_cluster_0':[0, 1],
                                          'g2_my_data_type_cluster_1':[0, 0]},
                                        index=['101', '102'])


    def test_make_feature_table(self):

        answer = make_feature_table(self.clusters, self.row_ids,
                self.new_feat_names, self.data_type)

        assert_frame_equal(answer, self.right_answer)


class TestMain(unittest.TestCase):

    def setUp(self):

        self.tmpdir = 'tmp'
        os.makedirs(self.tmpdir, exist_ok=True)

        df = pd.DataFrame({'g1_my_data_type_R_0':[0.5, 1],
                            'g1_my_data_type_R_1':[1, 0.5],
                            'g2_my_data_type_R_0':[0, 0],
                            'g2_my_data_type_R_1':[1, 1]},
                            index=[101, 102])
        self.df_path = f'{self.tmpdir}/df.csv'
        df.to_csv(self.df_path)

        self.data_type = 'my_data_type'
        self.num_clusters = 2

        self.output_path = f'{self.tmpdir}/df_2_clusters.csv'

        self.right_answer = pd.DataFrame({'g1_my_data_type_cluster_0':[1, 0],
                                          'g1_my_data_type_cluster_1':[0, 1],
                                          'g2_my_data_type_cluster_0':[1, 1],
                                          'g2_my_data_type_cluster_1':[0, 0]},
                                        index=[101, 102])


    def tearDown(self):

        shutil.rmtree(self.tmpdir)


    def test_main(self):

        cmf.main(self.df_path, self.data_type, self.num_clusters, self.tmpdir)

        answer = pd.read_csv(self.output_path, index_col=0)

        assert_frame_equal(answer, self.right_answer, check_dtype=False)


if __name__ == "__main__":
    unittest.main()
