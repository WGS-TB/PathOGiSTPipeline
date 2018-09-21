#!/usr/bin/env python
import os
import sys
import resource
import argparse
import logging
import numpy
import itertools
import re
import collections
import pkg_resources
import shutil
import yaml

import pathogist
import pathogist.cluster
import pathogist.io
import pathogist.distance
import pathogist.visualize

logger = logging.getLogger()

def run_all(param):
    '''
    Run the entire PathOGiST pipeline from distance matrix creation to consensus clustering, or create
    a new configuration file.
    '''
    if param.new_config:
        # Copy the default configuration file to whereever the user has specified
        src_path = pkg_resources.resource_filename(__name__,'pathogist/resources/blank_config.yaml') 
        shutil.copyfile(src_path,param.config) 
        print("New configuration file written at %s" % param.config)
    else:
        with open(param.config,'r') as config_stream:
            try:
                config = yaml.load(config_stream) 
            except yaml.YAMLError:
                print(yaml.YAMLError)
                sys.exit(1)

        # Make sure the configuration file is formatted correctly 
        distance_keys_set = set(config['distances'].keys())
        genotyping_keys_set = set(config['genotyping'].keys())
        threshold_keys_set = set(config['thresholds'].keys())
        fine_clusterings_set = set(config['fine_clusterings'])
        assert( (distance_keys_set & genotyping_keys_set) == set() ),\
            "'distances' and 'genotyping' have a key in common."
        assert( threshold_keys_set == (distance_keys_set | genotyping_keys_set) ),\
            "Set of keys in 'thresholds' not equal to the set of keys in 'genotyping' and 'distances'."
        assert( fine_clusterings_set <= (distance_keys_set | genotyping_keys_set) ),\
            "A value in 'fine_clusterings' does not appear in 'genotyping' or 'distances'."

        # Get genotyping calls
        logger.info('Reading genotyping calls ...')
        read_genotyping_calls = {'SNP': pathogist.io.read_snp_calls,
                                 'MLST': pathogist.io.read_mlst_calls,
                                 'CNV': pathogist.io.read_cnv_calls} 
        calls = {}
        for genotype in config['genotyping'].keys():
            calls[genotype] = read_genotyping_calls[genotype](config['genotyping'][genotype]) 
        # Create distance matrices from calls
        logger.info('Creating distance matrices ...')
        create_genotype_distance = {'SNP': pathogist.distance.create_snp_distance_matrix,
                                    'MLST': pathogist.distance.create_mlst_distance_matrix,
                                    'CNV': pathogist.distance.create_cnv_distance_matrix} 
        distances = {}
        for genotype in calls:
            distances[genotype] = create_genotype_distance[genotype](calls[genotype]) 

        # Read pre-constructed distance matrices
        logger.info('Reading distance matrices ...')
        for genotype in config['distances'].keys():
            distances[genotype] = pathogist.io.open_distance_file(config['distances'][genotype]) 

        # Match the distance matrices if need be
        distance_matrix_samples = [set(distances[key].columns.values) for key in distances]

        if (len(set(distance_matrix_samples)) > 1):
            logger.info('Warning: samples differ across the distance matrices.')
            logger.info('Matching distance matrices ...')
            distances = pathogist.distance.match_distance_matrices(distances)
            
        # dummy variables to make life easier
        genotypes = distances.keys()
        thresholds = config['thresholds']
        all_constraints = config['all_constraints'] 
        output_path = config['output']
        fine_clusterings = config['fine_clusterings']

        clusterings = {}
        for genotype in genotypes:
            logger.info('Clustering %s ...' % genotype)
            clusterings[genotype] = pathogist.cluster.correlation(distances[genotype],
                                                                  thresholds[genotype],
                                                                  all_constraints)     
        logger.info('Performing consensus clustering ...')
        consensus_clustering = pathogist.cluster.consensus(distances,clusterings,fine_clusterings)
        summary_clustering = pathogist.cluster.summarize_clusterings(consensus_clustering,clusterings)
        logger.info('Writing clusterings to file ...')
        pathogist.io.output_clustering(summary_clustering,output_path)
        
         
def correlation(param):
    logger.debug("Opening distance matrix...")
    distance_matrix = pathogist.io.open_distance_file(param.distance_matrix)
    logger.debug("Creating and solving correlation clustering problem ... ")
    clustering = pathogist.cluster.correlation(distance_matrix,param.threshold,param.all_constraints)
    logger.debug("Outputting clustering...")
    pathogist.io.output_clustering(clustering,param.output_path)

def consensus(param):
    logger.info("Getting distance matrices ...")
    distances = collections.OrderedDict()
    with open(param.distance_matrices,'r') as file:
        for line in file:
            name,path = line.rstrip().split('=')
            distances[name] = pathogist.io.open_distance_file(path)

    for cluster1,cluster2 in itertools.combinations(distances.keys(),2):
        columns1 = sorted(list(distances[cluster1].columns.values))
        columns2 = sorted(list(distances[cluster2].columns.values))
        assert( len(columns1) == len(columns2) )
        assert( columns1 == columns2 )
        rows1 = sorted(list(distances[cluster1].index.values))
        rows2 = sorted(list(distances[cluster1].index.values))
        assert( len(rows1) == len(rows2) )
        assert( rows1 == rows2 )

    logger.info("Getting clusterings ...")
    clustering_vectors = collections.OrderedDict()
    clusterings = collections.OrderedDict()
    with open(param.clusterings,'r') as file:
        for line in file:
            cluster,path = line.rstrip().split('=')
            clusterings[cluster] = pathogist.io.open_clustering_file(path)

    for cluster1,cluster2 in itertools.combinations(clusterings.keys(),2):
        columns1 = sorted(list(clusterings[cluster1].columns.values))
        columns2 = sorted(list(clusterings[cluster2].columns.values))
        assert( len(columns1) == len(columns2) )
        assert( columns1 == columns2 )
        rows1 = sorted(list(clusterings[cluster1].index.values))
        rows2 = sorted(list(clusterings[cluster1].index.values))
        assert( len(rows1) == len(rows2) )
        assert( rows1 == rows2 )

    logger.info("Getting other metadata ...")
    fine_clusterings = []
    with open(param.fine_clusterings,'r') as file:
        for line in file:
            fine_clusterings.append( line.rstrip() )
    logger.info("Creating and solving consensus clustering problem ...")
    consensus_clustering = pathogist.cluster.consensus(distances,clusterings,fine_clusterings)
    summary_clustering = pathogist.cluster.summarize_clusterings(consensus_clustering,
                                                                     clusterings)
    logger.info("Writing clusterings to file ...")
    pathogist.io.output_clustering(summary_clustering,param.output_path)

def distance(param):
    logger.info("Creating distance matrix ...")
    distance_matrix = None
    read_genotyping_calls = {'SNP': pathogist.io.read_snp_calls,
                             'MLST': pathogist.io.read_mlst_calls,
                             'CNV': pathogist.io.read_cnv_calls} 
    create_genotype_distance = {'SNP': pathogist.distance.create_snp_distance_matrix,
                                'MLST': pathogist.distance.create_mlst_distance_matrix,
                                'CNV': pathogist.distance.create_cnv_distance_matrix} 
    calls = read_genotyping_calls[param.data_type](param.calls_path)
    distance_matrix = create_genotype_distance[param.data_type](calls)
    if distance_matrix is not None:
        logger.info("Writing distance matrix ...")
        pathogist.io.write_distance_matrix(distance_matrix,param.output_path)
        logger.info("Distance matrix creation complete!")

def visualize(param): 
    logger.info("Visualizing distance matrix ...")
    distance_matrix = pathogist.io.open_distance_file(param.distance_matrix)
    if param.save_pdf:
        pathogist.visualize.visualize(distance_matrix,param.sample_name,pdf = None)
    else:
        pathogist.visualize.visualize(distance_matrix,param.sample_name)

def main():
    MAJOR_VERSION = 1
    MINOR_VERSION = 0

    parser = argparse.ArgumentParser(description=('PathOGiST Version %d.%d\n' +
                    'Copyright (C) 2018 Leonid Chindelevitch, Cedric Chauve, William Hsiao')
                    % (MAJOR_VERSION, MINOR_VERSION), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-ll', '--loglevel', type=str, default="INFO",
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        help='Set the logging level')
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    # command line arguments to run entire pipeline
    all_parser = subparsers.add_parser(name='all', help='run entire PathOGiST pipeline')
    all_parser.add_argument("config", metavar="CONFIG", type=str, 
                     help='path to input configuration file, or path to write a new configuration file')
    all_parser.add_argument("-n","--new_config", action="store_true", default=False,
                            help="write a blank configuration file at path given by CONFIG")

    # Correlation clustering command line arguments
    corr_parser = subparsers.add_parser(name='correlation', help="perform correlation clustering")
    corr_parser.add_argument("distance_matrix", type=str,
                             help="path to the distance matrix file")
    corr_parser.add_argument("threshold", type=float,help="threshold value for correlation")
    corr_parser.add_argument("output_path", type=str, help="path to write cluster output tsv file")
    corr_parser.add_argument("-a", "--all_constraints", action="store_true", default=False,
                            help = "add all constraints to the optimization problem, "
                                 + "not just those with mixed signs.")

    # Consensus clustering command line arguments
    cons_parser = subparsers.add_parser(name='consensus',
                                        help='perform consensus clustering on multiple clusterings')
    cons_parser.add_argument("distance_matrices", type=str,
                             help = "path to file containing paths to distance matrices for different"
                                  + " clusterings")
    cons_parser.add_argument("clusterings", type=str,
                             help = "path to file containing paths to clusterings, represented as"
                                  + " either matrices or lists of clustering assignments")
    cons_parser.add_argument("fine_clusterings", type=str,
                             help = "path to file containing the names of the clusterings which are "
                                  + "the finest")
    cons_parser.add_argument("output_path", type=str, help="path to output tsv file")
    cons_parser.add_argument("-a", "--all_constraints", action="store_true", default=False,
                            help = "add all constraints to the optimization problem, "
                                 + " not just those with mixed signs.")

    # Distance command line arguments
    distance_parser = subparsers.add_parser(name='distance', help = "construct distance matrix from "
                                                                  + "genotyping data")
    distance_parser.add_argument("calls_path", type=str,
                             help = "path to file containing paths to signal calls "
                                  + "(e.g. MLST calls, CNV calls, etc)")
    distance_parser.add_argument("data_type", type=str, choices=['MLST','CNV','SNP'],
                             help = "genotyping data")
    distance_parser.add_argument("output_path", type=str, help="path to output tsv file")

    # Visualization command line arguments
    vis_parser = subparsers.add_parser(name='visualize',help="visualize distance matrix")
    vis_parser.add_argument("distance_matrix",type=str,help="path to distance matrix in tsv format")
    vis_parser.add_argument("-p","--save_pdf",type=str,metavar='DIR',help="save PDF in given directory")
    vis_parser.add_argument("-n","--sample_name",type=str,default="sample",help="name of the sample")

    param = parser.parse_args()

    logging.basicConfig(level=param.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s:%(message)s',
                        datefmt='%I:%M:%S %p')

    if param.subcommand == 'all':
        run_all(param)
    elif param.subcommand == 'correlation':
        correlation(param)
    elif param.subcommand == 'consensus':
        consensus(param)
    elif param.subcommand == 'distance':
        distance(param)
    elif param.subcommand == 'visualize':
        visualize(param)

if __name__ == "__main__":
    main()