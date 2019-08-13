#!/bin/bash

set -e

wget -q https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh -O anaconda.sh
chmod +x anaconda.sh
./anaconda.sh -b
export PATH=/home/travis/anaconda3/bin:$PATH
#ls /home/travis/
conda update --yes conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#conda install kwip --yes
#free -m
#kwip -d tests/integration_tests/temp_dir/kwip_dist.tsv tests/integration_tests/temp_dir/SRR6152717.ct.gz tests/integration_tests/temp_dir/SRR6153024.ct.gz

#conda install conda-build --yes
#conda install --file ../../conda/environment.yaml --yes 

conda create -q --name pathogist --yes 
source activate pathogist 
conda install pathogist --yes -q


PATHOGIST run tests/integration_tests/test_data/pathogist-run_all-test.yaml 
python -m unittest tests/integration_tests/test_integration.py

# Run Genotyped File to Clustering
PATHOGIST run tests/integration_tests/test3_data/pathogist-run_all-test3.yaml 
python -m unittest tests/integration_tests/test3_integration.py

# Run Distance matrices to Clustering
PATHOGIST run tests/integration_tests/test4_data/pathogist-run_all-test4.yaml 
python -m unittest tests/integration_tests/test4_integration.py

# Run Genotyping Software
./PATHOGIST run tests/integration_tests/test2_data/pathogist-run_all-test2.yaml 
python -m unittest tests/integration_tests/test2_integration.py
