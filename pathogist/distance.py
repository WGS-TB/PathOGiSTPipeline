import logging
import sys
import itertools
import numpy
import pandas
''' These imports are for the suffix array hamming distance function. Not used.
import pathogist.suffix_array_distance
import pathogist.suffix_array_distance.suffix_array
import pathogist.suffix_array_distance.suffix_array_AlgorithmicAlley
import pysais
'''
import time

logger = logging.getLogger(__name__)

def hamming_distance(calls1,calls2):
    '''
    Given two numpy 1-dimensional arrays calls1, calls2 of the same shape, returns the hamming distance.
    '''
    return numpy.count_nonzero(calls1 != calls2)
	
def spoligo_ham_distance(calls1,calls2):
    assert len(calls1) == len(calls2)
    return sum(c1 != c2 for c1, c2 in zip(calls1, calls2))

def l1_norm(calls1,calls2):
    '''
    Given two numpy 1-dimensional arrays calls1, calls2 of the same shape, returns the L1 norm.
    '''
    return numpy.linalg.norm(calls1-calls2,ord=1)

def create_mlst_distance_matrix(calls):
    '''
    Given a dictionary of MLST calls (where sample names are keys to a vector), creates an MLST
    distance matrix represented as a Pandas Dataframe object.
    Distance: hamming distance
    '''
    samples = calls.keys()
    num_samples = len(samples)
    logger.debug("Got " + str(num_samples) + " samples...")
    distance_matrix = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=int),
                                       index=samples,columns=samples,dtype=int)
    for sample1,sample2 in itertools.combinations(samples,2):
        distance_matrix[sample1][sample2] = hamming_distance(calls[sample1],calls[sample2])
        distance_matrix[sample2][sample1] = distance_matrix[sample1][sample2]
    return distance_matrix

def create_cnv_distance_matrix(calls):
    '''
    Given a dictionary of CNV calls (where sample names are keys to a vector), creates a CNV
    distance matrix represented as a Pandas Dataframe object.
    Distance: L1 norm
    '''
    samples = calls.keys()
    num_samples = len(samples)
    distance_matrix = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=int),
                                       index=samples,columns=samples,dtype=float)
    for sample1,sample2 in itertools.combinations(samples,2):
        distance_matrix[sample1][sample2] = l1_norm(calls[sample1],calls[sample2])
        distance_matrix[sample2][sample1] = distance_matrix[sample1][sample2]
    return distance_matrix

def create_snp_distance_matrix(calls):
    '''
    Given a dictionary of SNP calls (where sample names are keys to a vector), creates an SNP
    distance matrix represented as a Pandas Dataframe object.
    Distance: hamming distance
    '''
    # Assume that the SNP distance matrix is created in the same was as the MLST distance matrix.
    return create_mlst_distance_matrix(calls)


def create_spotype_distance_matrix(calls):
    '''
    Given a dictionary of SpoTyping calls (where sample names are keys to a vector), creates a Spoligotyping
    distance matrix represented as a Pandas Dataframe object.
    Distance: hamming distance
    '''
    samples = calls.keys()
    num_samples = len(samples)
    logger.debug("Got " + str(num_samples) + " samples...")
    distance_matrix = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=int),
                                       index=samples,columns=samples,dtype=int)

    for sample1,sample2 in itertools.combinations(samples,2):
        distance_matrix[sample1][sample2] = hamming_distance(calls[sample1],calls[sample2])
        distance_matrix[sample2][sample1] = distance_matrix[sample1][sample2]
    return distance_matrix

def create_spoligo_distance_matrix(calls):
    '''
    Given a dictionary of spoligotyping calls (where samples names are keys to vectors), creates a
    spoligotype distance matrix represented as a pandas dataframe object.
    Distance: hamming distance
    '''
    # Assume that the spoligo distance matrix is created in the same was as the MLST distance matrix.
    return create_mlst_distance_matrix(calls)


def match_distance_matrices(distances):
    '''
    Modify a set of distance_matrices so that they all share the same set of samples
    '''
    samples_intersected = list(set.intersection(*[set(distances[key].columns.values)
                                               for key in distances]))
    matched_distances = {key: distances[key].loc[samples_intersected,samples_intersected]
                         for key in distances.keys()}
    return matched_distances

def aligned(s_i, L, m):
    """
    Let l = si mod m, i.e., the starting position of the suffix si within a profile.
    Then this subroutine returns l/L if l is multiple of L, and −1 otherwise.
    """
    l = s_i % m
    return int(l / L) if l % L == 0 else -1

def HD(p_i, p_j, s, m):
    """
    Given two profiles pi and pj which share a substring of length L,
    this subroutine computes the Hamming distance between pi and pj.
    """
    return hamming_distance(s[p_i * m: p_i * m + m], s[p_j * m: p_j * m + m])

def fast_ham_distance(calls, k):
    """
    An average‐case linear‐time algorithm to compute
    pairwise Hamming distances among a set of taxa,
    under a given Hamming distance threshold.

    Input: A dictionary P of d profiles of length m each; an integer threshold 0 < k < m.
    Output: The matrix X of hamming distance between pairs of profiles that are at Hamming distance at most k
    """
    elapsed_1 = 0
    elapsed_12 = 0
    elapsed_2 = 0
    elapsed_3 = 0

    # Initialization:
    samples = calls.keys()
    P = numpy.array(list(calls.values()))
    d = len(P)
    m = len(P[0])
    s = numpy.concatenate(P)
    n = m * d
    L = int(m / (k+1))
    s_t1 = time.time()
    s_int = numpy.array(pathogist.suffix_array.to_int_keys_best(s), dtype=numpy.int32)
    S = pysais.sais_int(s_int, max(s_int) + 1)
    # S, inv_S = pathogist.suffix_array.suffix_array_best(s)
    # S = pathogist.suffix_array_AlgorithmicAlley.suffix_array_ManberMyers(s)
    inv_S = pathogist.suffix_array.inverse_array(S)
    e_t1 = time.time()
    LCP = pathogist.suffix_array.lcp(s, S, inv_S)
    # LCP, _, _ = pysais.lcp_int(s_int, S)
    e_t12 = time.time()

    # RMQ_lcp = RangeMinimum(LCP)
    HT = dict()
    elapsed_1 = (e_t1 - s_t1)
    elapsed_12 = (e_t12 - e_t1)

    # Candidate pairs enumeration:
    X = numpy.ones((d, d)) * k
    numpy.fill_diagonal(X, 0)
    l_p = -1
    C = [set() for i in range(0, int(m / L) + 1)]
    for i in range(1, n):
        l = LCP[i]
        if l >= L:
            p_i = int(S[i] / m)
            x = aligned(S[i], L, m)
            if x != -1:
                C[x].add(p_i)
            if l_p == -1:
                p_i_1 = int(S[i-1] / m)
                x = aligned(S[i-1], L, m)
                if x != 1:
                    C[x].add(p_i_1)
            l_p = l
        elif l_p != -1:
            # Pairs enumeration:
            s_t3 = time.time()
            for C_t in C:
                for p, q in itertools.combinations(C_t, 2):
                    if (p, q) not in HT:
                        HT[p, q] = 1
                        s_t2 = time.time()
                        delta = HD(p, q, s, m)
                        e_t2 = time.time()
                        elapsed_2 += (e_t2 - s_t2)
                        if delta <= k:
                            X[p][q] = delta
                            X[q][p] = delta
                C_t.clear()
            l_p = -1
            e_t3 = time.time()
            elapsed_3 += (e_t3 - s_t3)

    distance_matrix = pandas.DataFrame(X, index=samples, columns=samples, dtype=int)

    print('Suffix array construction:', elapsed_1)
    print('LCP array construction:', elapsed_12)
    print('HD:', elapsed_2)
    print('Pairs enumeration:', elapsed_3)
    return distance_matrix
