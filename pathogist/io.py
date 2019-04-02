import pandas
import numpy
import re
import os
from collections import defaultdict

def get_sample_name(forward, reverse):
    fastq_list=[re.sub(".*/","",forward), re.sub(".*/","",reverse)]
    prefix = os.path.commonprefix(fastq_list)
    return re.sub("_", "", prefix)


def open_clustering_file(path):
    '''
    Reads a clustering file represented in CSV format.
    Assumes the file has a header, so skips first line.
    '''
    clustering = pandas.read_csv(path,header=0,index_col=0,sep='\t') 
    return clustering

def open_distance_file(filename):
    '''
    Reads distance matrix file represented in CSV format.
    Returns distance matrix as a pandas DataFrame matrix.
    '''
    distance = pandas.read_csv(filename,header=0,index_col=0,sep='\t')
    assert( distance.values.shape[0] == distance.values.shape[1] ),\
        "Distance matrix isn't square."
    return distance

def read_mlst_calls(calls_path):
    '''
    Read MentaLiST MLST calls.
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip().split('=')[0]
            with open(call_path,'r') as call_file:
                # Skip the header
                call_file.readline()
                for line in call_file:
                    columns = line.rstrip().split('\t')
                    sample = columns[0]
                    calls[sample] = numpy.array(columns[1:],dtype=numpy.string_)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of MLST calls."
    return calls

def read_snp_calls(calls_path, bed_path = ''):
    '''
    Read SNP calls from a text file with list of Snippy tsv or snippy-core output.
    '''
    calls = {}
    bed_filter = {}
    if bed_path != '':
        with open(bed_path,'r') as bed_file:
            columns=bed_file.readline().rstrip().split('\t')
            if columns[0].lower() != "chrom":
                for i in range(int(columns[1]),int(columns[2])):
                    filter_pos = columns[0]+"_"+str(i)
                    bed_filter[filter_pos] = 1        
            for line in bed_file:
                columns = line.rstrip().split('\t')
                for i in range(int(columns[1]),int(columns[2])):
                    filter_pos = columns[0]+"_"+str(i)
                    bed_filter[filter_pos] = 1        
            
    if calls_path.endswith("txt"):
        ref = {}
        sample = {}
        pos_count = {}
        #initialize dictionaries
        pos_count = defaultdict(int)
        sample = defaultdict(dict)
        with open(calls_path,'r') as calls_file:
            for line in calls_file:
                # We assume that each line of the calls_file is in the form
                # call_path=sample_name
                call_path = line.rstrip().split('=')[0]
                with open(call_path,'r') as call_file:
                    sample_name=call_file.readline().rstrip()
                    sample[sample_name][-1] = -1
                    for line in call_file:
                        columns = line.rstrip().split('\t')
                        #ensure there is one entry in a sample with no snps calls 
                        sample[sample_name][-1] = -1
                        if columns[2] == "snp":
                            chrom_pos = columns[0]+"_"+columns[1]
                            pos_count[chrom_pos] += 1
                            sample[sample_name][chrom_pos] = columns[4]
                            ref[chrom_pos] = columns[3]
        pos_count = {k:v for (k,v) in pos_count.items() if (v > 1 and v != len(sample.keys()) and k not in bed_filter ) }
        for sample_name in sample.keys():
            snps = []
            sample_name_keys = sample[sample_name].keys()
            for pos in pos_count.keys():
                snps.append(sample[sample_name][pos]) if pos in sample_name_keys else snps.append(ref[pos])
            calls[sample_name] = numpy.array(snps ,dtype="S1")    
    # Using snippy-core output as input
    if calls_path.endswith("tab"):
        snps_union = pandas.read_csv(calls_path, sep='\t')
        snps_union = snps_union.drop(['CHR','POS','Reference','LOCUS_TAG' ,'GENE', 'PRODUCT', 'EFFECT'], axis=1)
        for column in snps_union:
            if column == "CHROM" or column == "POS" or column == "TYPE" or column == "REF":
                continue
            calls[column] = numpy.array(snps_union[column], dtype="S20")
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of SNP calls."
    return calls


def read_cnv_calls(calls_path):
    '''
    Read PRINCE CNV calls.
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip().split('=')[0]
            with open(call_path,'r') as call_file:
                # Skip the header
                call_file.readline()
                for line in call_file:
                    columns = line.rstrip().split('\t')
                    sample = columns[0]
                    calls[sample] = numpy.array(columns[1:],dtype=float)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of CNV calls."
    return calls

def read_spotype_calls(calls_path):
    '''
    Read SpoTyping calls.
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.strip()
            with open(call_path, 'r') as call_file:
                call = call_file.read()
                values = call.split("\t")
                seq = values[0].split("&")
                forward = seq[0]
                reverse = seq[1]
                sample = get_sample_name(forward, reverse)
                spoligotype = []
                if(len(values[1]) == 43):
                    for char in str(values[1]):
                        spoligotype.append(int(char))
                    calls[sample] = numpy.array(spoligotype)
        assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
            "Samples do not have the same number of Spoligotype calls."
    return calls
    


def output_clustering(clustering,output_path):
    '''
    Writes a clustering to file in TSV format.
    '''
    clustering.to_csv(output_path,index=True,sep='\t')

def write_distance_matrix(distance_matrix,output_path):
    '''
    Writes a distance matrix to file in TSV format.
    '''
    distance_matrix.to_csv(output_path,sep='\t')
