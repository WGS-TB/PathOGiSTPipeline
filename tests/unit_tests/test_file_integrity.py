import numpy
import pandas as pd
import pandas.testing as pt
import sys
import unittest 
sys.path.append('../..')
import pathogist
import pathogist.io
import pathogist.cluster
import yaml
import os

class FileIntegrityTest(unittest.TestCase):

    def setUp(self):
        self.config_path = 'tests/unit_tests/test_data/file_integrity/config.yaml'
        self.forward_path =  'tests/integration_tests/test2_data/tb_2_forward.txt'
        self.reverse_path =  'tests/integration_tests/test2_data/tb_2_reverse.txt'

    def test_config(self):
        with open(self.config_path,'r') as config_stream:
            try:
                config = yaml.load(config_stream) 
            except yaml.YAMLError:
                print(yaml.YAMLError)
                sys.exit(1)
        assert pathogist.io.assert_config(config) == 0

    def test_fastq_input(self):
        forward_reads_paths = {}
        reverse_reads_paths = {}

        with open(self.forward_path,'r') as forwards_file:
            for line in forwards_file:
                path = line.rstrip() 
                # basename of the FASTQ file
                base = os.path.basename(path)
                # remove '_1.fastq'
                accession = os.path.splitext(base)[0].split('_')[0]
                forward_reads_paths[accession] = path

        with open(self.reverse_path,'r') as reverse_file:
            for line in reverse_file:
                path = line.rstrip() 
                # basename of the FASTQ file
                base = os.path.basename(path)
                # remove '_2.fastq' from basename
                accession = os.path.splitext(base)[0].split('_')[0]
                reverse_reads_paths[accession] = path

        assert pathogist.io.check_fastq_input(forward_reads_paths, reverse_reads_paths) == 0
    