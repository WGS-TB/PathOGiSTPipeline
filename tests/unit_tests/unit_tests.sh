#!/bin/bash

python -m unittest tests/unit_tests/test_cluster.py
python -m unittest tests/unit_tests/test_distance.py
python -m unittest tests/unit_tests/test_file_integrity.py