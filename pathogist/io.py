import pandas
import numpy
import re
import os
from collections import defaultdict
import subprocess
#import sys

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

def read_mlst_calls(calls_paths):
    '''
    Read MentaLiST MLST calls.
    '''
    # If calls_paths is a string, we assume that it's a path to a file containing the calls paths,
    # and we read the paths from the file.
    if isinstance(calls_paths,str):
        calls_paths_path = calls_paths
        calls_paths = []
        with open(calls_paths_path,'r') as calls_paths_file:
            for line in calls_paths_file:
                calls_paths.append(line.rstrip())
    calls = {}
    #print(calls_paths)
    # Go through the calls file for each sample, and record the MLST calls
    for calls_path in calls_paths:
        calls_path = calls_path.rstrip().split('=')[0]
        with open(calls_path,'r') as calls_file:
            # Skip the header
            calls_file.readline()
            for line in calls_file:
                columns = line.rstrip().split('\t')
                sample = columns[0]
                calls[sample] = numpy.array(columns[1:],dtype=numpy.string_)

    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of MLST calls."

    return calls

def read_snp_calls(calls_paths, bed_path = ''):
    
    ###Read SNP calls from a text file with list of Snippy tsv or snippy-core output.
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
               
    calls = {}
    # If calls_paths is a list of paths or a text file containing the paths, do this. 
    if isinstance(calls_paths,list) or calls_paths.endswith("txt"):
        # collect the paths to the calls into a list
        if isinstance(calls_paths,str):
            calls_paths_path = calls_paths
            calls_paths = []
            with open(calls_paths_path,'r') as calls_paths_file:
                for line in calls_paths_file:
                    calls_path = line.rstrip().split('=')[0]
                    calls_paths.append(calls_path)
        # now get the SNP calls for each of the samples
        ref = {}
        sample = {}
        pos_count = {}
        #initialize dictionaries
        pos_count = defaultdict(int)
        sample = defaultdict(dict)
        for calls_path in calls_paths:
            with open(calls_path,'r') as call_file:
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
                if pos in sample_name_keys:
                    snps.append(sample[sample_name][pos])
                else:
                    snps.append(ref[pos])
            calls[sample_name] = numpy.array(snps ,dtype="S1")    

    # Otherwise, it must be snippy-core output 
    elif calls_paths.endswith("tab"):
        snps_union = pandas.read_csv(calls_path, sep='\t')
        snps_union = snps_union.drop(['CHR','POS','Reference','LOCUS_TAG' ,'GENE', 'PRODUCT', 'EFFECT'], axis=1)
        for column in snps_union:
            if column == "CHROM" or column == "POS" or column == "TYPE" or column == "REF":
                continue
            calls[column] = numpy.array(snps_union[column], dtype="S20")
    return calls

def read_cnv_calls(calls_path):
    '''
    Read PRINCE CNV calls.
    @param calls_path: a string to the PRINCE calls file
    @rvalue calls: a dictionary whose keys are sample names, and values are CNV calls represented by 
                   a numpy array
    '''

    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip().split('=')[0]
            with open(call_path,'r') as call_file:
                # Skip the header
                call_file.readline()
                for line in call_file:
                    if "\t" in line:
                        columns = line.rstrip().split('\t')
                    if "," in line:
                        columns = line.rstrip().split(',')
                    sample = re.sub("_.*", "", columns[0])
                    calls[sample] = numpy.array(columns[1:],dtype=float)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of CNV calls."
    return calls

def read_spotype_calls(calls_paths):
    '''
    Read SpoTyping calls.
    @param calls_paths: a list containing paths to SpoTyping calls files, OR a path to a file
                        containing the paths to the calls 
    @rvalue calls: a dictionary whose keys are sample names, and values are spoligotypes represented
                   by a numpy array 
    '''
    # if calls_paths is a string, we assume its a file containing the paths to the SpoTyping calls
    if isinstance(calls_paths,str):
        calls_paths_path = calls_paths
        calls_paths = []
        with open(calls_paths_path,'r') as calls_paths_file:
            for line in calls_paths_file:
                calls_path = line.rstrip().split('=')[0]
                calls_paths.append(calls_path)
    calls = {}
    for call_path in calls_paths:
        with open(call_path,'r') as call_file:
            # Skip the header
            #call_file.readline()
            for line in call_file:
                values = line.split("\t")
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

def assert_config(config):
    # Makes sure the configuration file is formatted correctly
    #temp directory value assertion
    assert os.path.isdir(config['temp']),\
        "Temp value in the config file is not a real directory"
    #threads value assertion
    assert isinstance(config['threads'], int) and config['threads'] > 0 , \
        "Threads value in the config file must be an integer and greater than 0"
    #run section assertion
    run_genotyping_tools = False
    for tool in config['run']:
        if config['run'][tool] == 1:
           run_genotyping_tools = True
        assert config['run'][tool] == 1 or config['run'][tool] == 0,\
            "The value for the %s under the run section must be 0 or 1 to indicate to not run and run respectively." % tool
    #genotyping section assertion
    for key in config['genotyping']:
        if key == "input_reads": # add assertions
            for type in config['genotyping'][key]:
                reads_path = config['genotyping'][key][type]
                if reads_path == None:
                    if run_genotyping_tools == False:
                        continue
                    else:
                        assert reads_path != None,\
                            "File location values in the %s section under genotyping section must exist when you are running genotyping tools." % type
                else:
                    if run_genotyping_tools == False:
                        continue
                    else:
                        assert os.path.isfile(reads_path),\
                            "File location values in the %s section under genotyping section must exist." % type
        # mentalist section assertion                        
        if key == "mentalist" and config['run'][key] == 1:
            #db_loc section assertion
            db_loc = False
            db_loc_count = 0
            for type in config['genotyping'][key]['db_loc']:
                assert config['genotyping'][key]['db_loc'][type] == 1 or config['genotyping'][key]['db_loc'][type] == 0,\
                    "The value for the %s under the db_loc section must be 0 or 1 to indicate where to locate mentalist db." % type
                if config['genotyping'][key]['db_loc'][type] == 1:
                    db_loc = True
                    db_loc_count += 1
            assert db_loc_count == 1,\
                "Choose only 1 of the options under db_loc for mentalist to obtain mlst database"
            assert db_loc == True, \
                "Choose 1 of the options under db_loc for mentalist to obtain mlst database by inputting 1"
            #local file section assertion
            if config['genotyping'][key]['db_loc']['local_file'] == 1:
                assert config['genotyping']['mentalist']['local_file']['database'] != None,\
                    "Database file value under mentalist section cannot be empty"
                assert os.path.isfile(config['genotyping']['mentalist']['local_file']['database']),\
                    "Database file value must mentalist section exist"
            # build_db assertion
            if config['genotyping'][key]['db_loc']['build_db'] == 1:
                assert isinstance(config['genotyping'][key]['build_db']['options']['k'], int) and config['genotyping'][key]['build_db']['options']['k'] > 0,\
                    "k value in build_db under mentalist must be an integer and greater than 0"
                assert os.path.isfile(config['genotyping'][key]['build_db']['options']['fasta_files']),\
                    "fasta_file value in build_db under mentalist must be a file path that exists"
                assert os.path.isfile(config['genotyping'][key]['build_db']['options']['profile']),\
                    "profile value in build_db under mentalist must be a file path that exists"
            # download_pubmlst assertion
            if config['genotyping'][key]['db_loc']['download_pubmlst'] == 1:
                assert isinstance(config['genotyping'][key]['download_pubmlst']['options']['k'], int) and config['genotyping'][key]['download_pubmlst']['options']['k'] > 0,\
                    "k value in download_pubmlst under mentalist must be an integer and greater than 0"
                assert config['genotyping'][key]['download_pubmlst']['options']['scheme']!= None, \
                    "Scheme cannot be none under download_mlst section"
            # download_cgmlst assertion
            if config['genotyping'][key]['db_loc']['download_cgmlst'] == 1:
                assert isinstance(config['genotyping'][key]['download_cgmlst']['options']['k'], int) and config['genotyping'][key]['download_cgmlst']['options']['k'] > 0,\
                    "k value in download_cgmlst under mentalist must be an integer and greater than 0"
                assert config['genotyping'][key]['download_cgmlst']['options']['scheme']!= None, \
                    "Scheme cannot be none under download_cgmlst section"
            # download_enterobase assertion
            if config['genotyping'][key]['db_loc']['download_enterobase'] == 1:            
                assert isinstance(config['genotyping'][key]['download_enterobase']['options']['k'], int) and config['genotyping'][key]['download_enterobase']['options']['k'] > 0,\
                    "k value in download_enterobase under mentalist must be an integer and greater than 0"
                assert config['genotyping'][key]['download_enterobase']['options']['scheme']!= None, \
                    "Scheme cannot be none under download_enterobase section"
                assert config['genotyping'][key]['download_enterobase']['options']['type'] in ['cg', 'wg'],\
                    "cg and wgs are the only values allowed in the type section of download_enterobase under mentalist"
            # call assertions             
            assert isinstance(config['genotyping'][key]['call']['options']['mutation_threshold'], int) and config['genotyping'][key]['call']['options']['mutation_threshold'] >= 0,\
                "mutation_threshold value in call under mentalist must be an integer and greater than or equal to 0"
            assert isinstance(config['genotyping'][key]['call']['options']['kt'], int) and config['genotyping'][key]['call']['options']['kt'] > 1,\
                "kt value in call under mentalist must be an integer and greater than 0"
            if config['genotyping'][key]['call']['flags'] != None: 
                for flag in config['genotyping'][key]['call']['flags']:
                    assert flag in ['output_votes', 'output_special'], "output_votes and output_special are the only values allowed in the flags of call under mentalist"
        # kwip assertions                    
        if key == "kwip" and config['run'][key] == 1:
            if config['genotyping'][key]['kwip_options'] != None:
                if config['genotyping'][key]['kwip_options']['weights'] != None:
                    assert os.path.isfile(config['genotyping'][key]['kwip_options']['weights']),\
                        "Weights in the kwip section under genotyping section must exist."
            if config['genotyping'][key]['kwip_flags'] != None:
                for flag in config['genotyping'][key]['kwip_flags']:
                    assert flag in ['unweighted', 'calc_weights'],\
                        "unweighted and calc_weights are the only flags allowed in the kwip_flags section"
        # prince assertions                    
        if key == "prince" and config['run'][key] == 1:
            if config['genotyping'][key]['options'] != None:
                if config['genotyping'][key]['options']['templates'] != None:
                    assert os.path.isfile(config['genotyping'][key]['options']['templates']),\
                        "Templates in the prince section under genotyping section must exist."
        # snippy assertions
        if key == "snippy" and config['run'][key] == 1:
            if config['genotyping'][key]['flags'] != None:
                for flag in config['genotyping'][key]['flags']:
                    assert flag in ['unmapped'],\
                        "unmapped is the only flag allowed in the snippy section"
            if config['genotyping'][key]['options'] != None:
                for option in config['genotyping'][key]['options']:
                    if option == "reference":
                        assert os.path.isfile(config['genotyping'][key]['options'][option]),\
                            "Reference value in the snippy section under genotyping section must exist."
                    if option == "mapqual":
                        assert (isinstance(config['genotyping'][key]['options'][option], int) or isinstance(config['genotyping'][key]['options'][option], float))  and  config['genotyping'][key]['options'][option] >= 0,\
                            "mapqual values under snippy section must be integer or float and be greater or equal to 0"
                    if option == "basequal":
                        assert (isinstance(config['genotyping'][key]['options'][option], int) or isinstance(config['genotyping'][key]['options'][option], float))  and  config['genotyping'][key]['options'][option] >= 0,\
                            "basequal values under snippy section must be integer or float and be greater or equal to 0"
                    if option == "mincov":
                        assert (isinstance(config['genotyping'][key]['options'][option], int) or isinstance(config['genotyping'][key]['options'][option], float))  and  config['genotyping'][key]['options'][option] >= 0,\
                            "mincov values under snippy section must be integer or float and be greater or equal to 0"
                    if option == "minfrac":
                        assert (isinstance(config['genotyping'][key]['options'][option], int) or isinstance(config['genotyping'][key]['options'][option], float))  and  config['genotyping'][key]['options'][option] >= 0 and config['genotyping'][key]['options'][option] <= 1,\
                            "minfrac values under snippy section must be integer or float and be between 0 and 1"
        # spotyping assertions
        if key == "spotyping" and config['run'][key] == 1:
            if config['genotyping'][key]['flags'] != None:
                for flag in config['genotyping'][key]['flags']:
                    assert flag in ['seq','noQuery','filter','sorted'],\
                        "seq, noQuery, filter, and sorted are the only flag allowed in the snippy section. Please look at the original config file for formatting"
            if config['genotyping'][key]['options'] != None:
                for option in config['genotyping'][key]['options']:
                    if option == "swift":
                        assert config['genotyping'][key]['options'][option] == "on" or config['genotyping'][key]['options'][option] == "off",\
                            "Swift value in the spotyping section must be on or off." 
                    if option == "min":
                        assert  isinstance(config['genotyping'][key]['options'][option], int) and config['genotyping'][key]['options'][option] >= 0,\
                            "min value in the spotyping section must an integer and greater than 0."         
                    if option == "rmin":
                        assert  isinstance(config['genotyping'][key]['options'][option], int) and config['genotyping'][key]['options'][option] >= 0,\
                            "rmin value in the spotyping section must an integer and greater than 0." 
                    if option == "outdir":
                        assert os.path.isdir(config['genotyping'][key]['options'][option]),\
                            "outdir value under the spotyping section in the config file is not a real directory. If error persists, try using the full directory path"
                    if option == "output":
                        assert "/" not in config['genotyping'][key]['options'][option],\
                            "Output value under the spotyping section cannot contain a forward slash."                     
    # clustering section assertion
    for key in config['clustering']:
        if key == "output_prefix":
            assert config['clustering'][key] != None,\
                "Output_prefix value in the clustering section cannot be None."
        if key == "genotyping":
            if config['clustering'][key] != None:
                for type in config['clustering'][key]:
                    calls_path = config['clustering'][key][type]
                    if calls_path == None:
                       continue
                    if type == 'SNP' and config['run']['snippy'] == 1:
                        assert False,\
                            "You've selected to run snippy in the genotyping section and supplied a list of SNP calls. Please choose only one"
                    if type == 'MLST' and config['run']['mentalist'] == 1:
                        assert False,\
                            "You've selected to run mentalist in the genotyping section and supplied a list of MLST calls. Please choose only one"
                    if type == 'CNV' and config['run']['prince'] == 1:
                        assert False,\
                            "You've selected to run prince in the genotyping section and supplied a list of CNV calls. Please choose only one"
                    if type == 'spoligotyping' and config['run']['spotyping'] == 1:
                        assert False,\
                            "You've selected to run spotyping in the genotyping section and supplied a list of spoligotyping calls. Please choose only one"
                    assert os.path.isfile(calls_path),\
                        "File location values in the genotyping section under clustering section must exist or be empty."
        if key == "genotyping_options":
            if config['clustering'][key]['bed_filter'] != None:
                os.path.isfile(config['clustering'][key]['bed_filter']),\
                    "File location values in the bed_filter section under clustering section must exist or be empty."
        if key == "distances":
            for type in config['clustering'][key]:
                dist_path = config['clustering'][key][type]
                if dist_path == None:
                   continue
                if type == 'SNP' and config['run']['snippy'] == 1:
                    assert False,\
                        "You've selected to run snippy in the genotyping section and supplied a SNP distance matrix. Please choose only one"
                if type == 'MLST' and config['run']['mentalist'] == 1:
                    assert False,\
                        "You've selected to run mentalist in the genotyping section and supplied a MLST distance matrix. Please choose only one"
                if type == 'CNV' and config['run']['prince'] == 1:
                    assert False,\
                        "You've selected to run prince in the genotyping section and supplied a CNV distance matrix. Please choose only one"
                if type == 'spoligotyping' and config['run']['spotyping'] == 1:
                    assert False,\
                        "You've selected to run spotyping in the genotyping section and supplied a spoligotyping distance matrix. Please choose only one"
                if type == 'kWIP' and config['run']['kwip'] == 1:
                    assert False,\
                        "You've selected to run kwip in the genotyping section and supplied a kwip distance matrix. Please choose only one"
                assert os.path.isfile(dist_path),\
                    "File location values in the distances section under clustering section must exist or be empty."
        if key == "fine_clusterings":
            assert config['clustering'][key] != None,\
                "fine_clusterings values must be a combination of at least 1 SNP, kWIP, MLST, CNV, and spoligotyping"
            for type in config['clustering'][key]:
                assert type in ['SNP','kWIP','MLST','CNV','spoligotyping'],\
                    "fine_clusterings values must be a combination of at least 1 SNP, kWIP, MLST, CNV, and spoligotyping"
        if key == "thresholds":
            for type in config['clustering'][key]:
                assert(isinstance(config['clustering'][key][type], int) or isinstance(config['clustering'][key][type], float))  and  config['clustering'][key][type] >= 0,\
                    "Threshold values under clustering section must be integer or float and be greater or equal to 0"
        if key == "all_constraints":
            assert config['clustering'][key] == True or config['clustering'][key] == False,\
                "all_constraints values must be either True or False"
        if key == "method":
            assert config['clustering'][key] == "ILP" or config['clustering'][key] == "C4",\
                "Method values must be either ILP or c4"
        if key == "visualize":
            assert config['clustering'][key] == True or config['clustering'][key] == False,\
                "visualize values must be either True or False"
    return 0
                    
def get_bases_and_reads_number(fastq_path):
    if fastq_path.rstrip().split('.')[-1] == 'gz':
        cmd="zcat " + fastq_path + "|wc -l"
    else:
        cmd="cat " + fastq_path + "|wc -l"
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    lines = output.decode("utf-8").rstrip()
    return lines
    
def check_fastq_input(forward_reads_paths,reverse_reads_paths):
    if len(forward_reads_paths) != len(reverse_reads_paths):
        sys.exit('Error! Check your fastq input files. They do not contain the same number of forward and reverse reads')
    shared_samples = set(forward_reads_paths.keys()).union(reverse_reads_paths.keys())
    if len(forward_reads_paths) != len(shared_samples) or len(reverse_reads_paths) != len(shared_samples):
        sys.exit('Error! Check your fastq input files. The forward and reverse reads contain different samples')
    bad_samples = []
    for sample in shared_samples:
        lines_forward = get_bases_and_reads_number(forward_reads_paths[sample])
        lines_reverse = get_bases_and_reads_number(reverse_reads_paths[sample])
        if abs(int(lines_forward) - int(lines_reverse)) > 1 or int(lines_forward)%4 != 0 or int(lines_reverse)%4 != 0:
            bad_sample.append(sample)
    if bad_samples != [] :
        print(bad_samples)
        sys.exit('Error! Check your fastq input files. The sample contains different number of lines or lines not divisible by four')
    return 0
            
