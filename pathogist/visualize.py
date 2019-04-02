import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cmx
import logging
import pandas 
import numpy 
from scipy.cluster import hierarchy as hcl
import sklearn
import sklearn.manifold
import networkx
import random
import itertools
import math

logger = logging.getLogger(__name__)

def visualize(distance, name, metadata = None, columns = None, save_path = None):
    '''
    :param distance: distance matrix
    :param name: name of the sample (and name of the pdf)
    :param metadata: location of metadata file
    :param columns: list of columns in metadata file for labels on hierarchical clustering dendogram
    :param create_pdf: variable for the creation of an output pdf
    '''
    distance_histogram(distance, name, save_path)

def distance_histogram(distance, name='SAMPLE', save_path=None):
    '''
    Creates a distance histogram from the distances in the
    lower triangle of the distance matrix

    Used to find the threshold for correlation clustering
    '''

    logger.debug("Creating distance histogram")

    # Obtain distances from lower triangle
    values = distance.values
    distances = values[numpy.tril_indices(len(distance.index), -1)]

    plt.style.use('ggplot')
    plt.figure(0)
    n,bins,patches = plt.hist(distances, bins=500)

    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title('Distance Histogram for %s' % name)
    if save_path is not None:
        plt.savefig(save_path)
    plt.show()

def visualize_clusterings(summary_clusterings,figure_size=(12,9),
                                              output_prefix=None,
                                              mode="spring",
                                              matrix=None):
    '''
    Visualize the consensus clustering and correlation clusterings.
    @param summary_clusterings: A Pandas Dataframe where the index is the sample names, columns are
                   clusterings, and entries are cluster assignments.
    @param mode: method to visualize the clusterings.
                 Options:
                 - spring: Spring layout. Does not require `matrix` parameter. `offset` is not used.
                 - tree: Represent clusters as trees and use spring layout to determine vertex 
                         position. Requires `matrix` parameter. 
    @param matrix: A Pandas Dataframe of the similarity matrix from the consensus clustering step of 
                   the problem
    @param output_prefix: The prefix, including directories, of the save path PNG of the 
                          visualization.
    '''
    plt.figure(figsize=figure_size)
    assert('Consensus' in summary_clusterings.columns.values),\
           'Please ensure your input clusterings contain a consensus clustering.'
    samples = summary_clusterings.index.values.tolist()
    if matrix is not None:
        assert(sorted(matrix.columns.values.tolist()) == sorted(samples))

    # Compute the positions of the samples in the graph based on either a spring layout, or
    # MDS of the consensus similarity matrix
    graph = networkx.Graph() 
    for sample in samples:
        graph.add_node(samples.index(sample))

    if mode == 'spring':
        consensus_clustering = summary_clusterings['Consensus']
        num_clusters = numpy.amax(consensus_clustering.values)
        ## Add an edge between samples that are in the same cluster
        for cluster in range(1,num_clusters+1):
            cluster_samples = consensus_clustering.index[consensus_clustering == cluster].tolist()
            for sample1,sample2 in itertools.combinations(cluster_samples,2):
                graph.add_edge(samples.index(sample1),samples.index(sample2))
        ## Compute the spring layout positions for each sample
        ### networkx default distance between nodes is 1/sqrt(N), where N is the number of vertices
        node_distances = 4/math.sqrt(len(samples)) 
        sample_positions = networkx.spring_layout(graph,dim=2,scale=1,k=node_distances)

        # Visualize all of the clusterings based on the spring layout
        clustering_names = summary_clusterings.columns.values.tolist()
        num_subplots = len(clustering_names)
        num_rows = int(math.ceil(num_subplots/2))
        num_cols = 2
        for clustering_index in range(0,num_subplots):
            subplot_index = clustering_index + 1
            ax = plt.subplot(num_rows,num_cols,subplot_index)
            clustering_name = clustering_names[clustering_index] 
            graph = networkx.Graph()
            for sample in samples:
                graph.add_node(samples.index(sample))
            clustering = summary_clusterings[clustering_name]
            num_clusters = numpy.amax(clustering.values)
            for cluster in range(1,num_clusters+1):
                cluster_samples = clustering.index[clustering == cluster].tolist()
                for sample1,sample2 in itertools.combinations(cluster_samples,2):
                    graph.add_edge(samples.index(sample1),samples.index(sample2))
            ax.set_title("%s Clustering" % clustering_name)
            networkx.draw_networkx(graph,
                                   pos=sample_positions,
                                   node_size=5,
                                   width=0.25,
                                   with_labels=False)
    elif mode == 'tree':
        assert(matrix is not None),'Error: please specify an input matrix.'
        consensus_clustering = summary_clusterings['Consensus']
        for sample1,sample2 in itertools.combinations(samples,2):
            graph.add_edge(samples.index(sample1),
                           samples.index(sample2),
                           weight=matrix[sample1][sample2])
        ## Compute the spring layout positions for each sample
        ### networkx default distance between nodes is 1/sqrt(N), where N is the number of vertices
        node_dists = 1/math.sqrt(len(samples)) 
        graph = networkx.minimum_spanning_tree(graph)
        sample_positions = networkx.spring_layout(graph,dim=2,scale=1,k=node_dists,iterations=500)
        node_colors = [consensus_clustering[sample] for sample in samples]
        networkx.draw(graph,
                      node_color=node_colors,
                      pos=sample_positions,
                      node_size=30,
                      width=0.25,
                      with_labels=False,
                      vmin=0,
                      vmix=max(node_colors),
                      cmap=cmx.tab20)
    else:
        print("Unregnized option.")
        sys.exit(1)

    if output_prefix is None:
        plt.show()
    else:
        output_path = '%s.png' % output_prefix
        plt.savefig(output_path)
        
        
def hierarchical_clustering(distance, name, metadata = None, columns = None, create_pdf = False, pdf = None):
    '''
    Runs hierarchical clustering and creates a dendrogram and MDS plot
    '''

    meta = pandas.read_csv(metadata, header=0, index_col=0, sep="\t")

    # Select a subset of samples common to the distance matrix and the metadata file
    isect = [rowname for rowname in list(distance.index) if rowname in list(meta.index)]
    meta = meta.loc[isect]
    distance = distance.loc[isect, isect]

    # Run hierarchical clustering
    logger.debug("Running hierarchical clustering")
    linkage = hcl.linkage(distance, method='complete')

    # Plotting dendrogram for original labels
    logger.debug("Generating dendrograms")
    plt.figure(0)
    dendro = hcl.dendrogram(linkage, labels=distance.index)

    # Plotting dendrograms for metadata labels
    for index,col in enumerate(columns):
        plt.figure(index+1)
        dendro = hcl.dendrogram(linkage, labels=meta[col])

def plot_nawc_values(nawc_values,title):
    title = species + " " + genotype + " " + algorithm
    x_values = nawc_values.index.values.astype(numpy.float)
    y_values = nawc_values.values
    plt.title(title)
    plt.plot(x_values,y_values)
    plt.ylabel('nAWC')
    plt.xlabel('threshold')
    plt.show() 
