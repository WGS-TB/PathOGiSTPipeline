#!/usr/bin/env python3
import sys
sys.path.append('../..')
import unittest
import numpy
import pandas
import pandas.testing
import pathogist
import pathogist.io

class IntegrationTest(unittest.TestCase):
    def setUp(self):
        true_clustering_path = 'tests/integration_tests/test3_data/tb_5_true_clustering.tsv'
        self.true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        self.true_clustering.sort_index(axis=1)
        true_clustering_2_path = 'tests/integration_tests/test3_data/tb_5_true_clustering_2.tsv'
        self.true_clustering_2 = pathogist.io.open_clustering_file(true_clustering_2_path)
        self.true_clustering_2.sort_index(axis=1)
        clustering_path = 'tests/integration_tests/test3_data/tb_5_final_clustering.tsv'
        self.clustering = pathogist.io.open_clustering_file(clustering_path)
        self.clustering.sort_index(axis=1)
       
    def test_integration(self):
        try:
            pandas.testing.assert_frame_equal(self.clustering.sort_index(axis=1),self.true_clustering.sort_index(axis=1))
        except:
            pandas.testing.assert_frame_equal(self.clustering.sort_index(axis=1),self.true_clustering_2.sort_index(axis=1))
