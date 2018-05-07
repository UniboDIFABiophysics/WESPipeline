import argparse
import pandas as pd
import subprocess as sp
import yaml
import os
import re
import math
import sys

# get home path
homepath = os.path.expanduser('~') + '/'

# get current path
currentpath = os.getcwd() + '/'



#### Parse command line arguments and load CONFIGFILE ####

parser = argparse.ArgumentParser()
parser.add_argument('-c', required=True, dest='configfile', action='store', help='path to config file')

# get arguments
args = parser.parse_args()

### get CONFIG file absolute path
configfile = args.configfile
# in case a relative path was given as input
if not os.path.commonprefix([homepath, configfile]):
    # append to current path
    configfile = currentpath + configfile

## load config file
with open(configfile) as infile:
    config = yaml.load(infile)




### build and create PROCESSPATH
processpath = homepath + 'ANALYSES_WES/' + config['outfolder'] + '/'
if not os.path.exists(processpath):
    os.makedirs(processpath)


# load dataset file and set sample names as indeces
data = homepath + config['dataset']
dataset = pd.read_table(data, sep = '\t', header=0, dtype='str')
dataset = dataset.set_index("sample_id")




##########################
### load sample input list

# move sample list file to processpath
input_list = homepath + config['input_list_launch_snake']
sp.run('cp %s %s' % (input_list, processpath), shell=True)

# load table
input_list = pd.read_table(input_list, sep= '\t', header=0, dtype='str')

# KIT CHECKPOINT
for i, row in input_list.iterrows():

    # get WES library kit for tumor and normal
    kit_tumor = dataset.loc[row.tumor, 'kit_wes']
    kit_normal = dataset.loc[row.normal, 'kit_wes']
    
    if kit_tumor != kit_normal:
        sys.exit('ERROR!! In pair %s - TUMOR and NORMAL have different WES kits!!!' %pair)


# create list of every pair
pairs = ['_'.join([row.tumor, row.normal]) for i, row in input_list.iterrows()]

# create list of every sample
samples = sorted(set(input_list.tumor.tolist() + input_list.normal.tolist()))




### GET SIZE OF EVERY FASTQ PRESENT IN THE DIRECTORIES WERE THE SAMPLES ARE STORED

# get list of every needed fastq path
every_path = list(set(dataset.loc[samples, 'fastq_path_wes']))
# add homepath
every_path = [homepath + path for path in every_path]

# create empty list
ls = []

# loop over fastq paths
for fastq_path in every_path:
    
    # list names and details of every fastq contained in fastq_path
    cmd = 'ls -l %s*fastq.gz' %fastq_path
    ls_temp = sp.run(cmd, shell=True, stdout=sp.PIPE).stdout.decode().splitlines()
    
    # append to growing list
    ls = ls + ls_temp

# remove all empty fields from 'ls -l' output (needed for following step)
ls = [[e2 for e2 in e.split(' ') if e2 != ''] for e in ls]

# extract all sizes
sizes = [float(e[4]) for e in ls]

# extract all fastq names
fastqs = [e[8].split('/')[-1] for e in ls]

# create a SERIES with sizes as data and fastq names as indeces
fastq_size = pd.Series(data=sizes, index=fastqs)


########################################################################





### ORDER PAIRS FOR FASTQ SIZE
### (the size of a pair will be the size of the largest sample in the pair)
### (the size of a sample is the total size of all its fastqs, both R1 and R2)

# create new list for total sizes of fastqs of each sample
tot_sizes = []

# loop over samples
for s in samples:
    
    # get common substring (CS) of every fastq file of current sample
    fastq_cs = dataset.loc[s, 'fastq_common_substring_wes']
    
    # extract the fastq file names that match the CS
    match = [f for f in fastqs if re.search(fastq_cs, f)]
    
    # sum their total size and store in growing list
    tot_sizes.append(sum(fastq_size[match]))

# create a SERIES with total sample sizes as data and sample names as indeces
sample_size = pd.Series(data=tot_sizes, index=samples)

# create empty list
largest_size = []

# loop over pairs
for p in pairs:
    
    # get names of tumor and normal
    names = p.split('_')
    
    # get size of samples, sort them and append the largest (the 2nd one) to growing list
    largest_size.append(sample_size[names].sort_values()[1])

# create SERIES with the largest sample in each pair (with pairs as indeces), sort and get indeces
ordered_pairs = pd.Series(data=largest_size, index=pairs).sort_values().index


########################################################################


# set --keep-going argument
if config['keepgoing']:
    keepgoing = '--keep-going'
else:
    keepgoing = ''



### SPLIT THE ORDERED LIST OF PAIRS INTO BATCHES OF A GIVEN NUMBER OF PAIRS
### LAUNCH SNAKEMAKE ON EACH BATCH, ONE BY ONE

# get total number of pairs
P = len(ordered_pairs)

# set size of batch
n = config['batch_size']

# change current directory to pipeline dir
os.chdir(homepath + 'WESPipeline/')

# loop over ordered pairs, in batches of n pairs
for b in range(0, P, n):
    
    # get pairs for current batch
    pairs_batch = ordered_pairs[b : b + n]
    
    # write to file, which will be the input list for SNAKEMAKE
    with open(processpath + 'pairs_batch', 'w') as outfile:
        outfile.write('\n'.join(pairs_batch))

    # build snakemake command to analyze current batch
    cmd = ' '.join(['snakemake',
                    '-s', config['snakefile'],
                    '--configfile', configfile,
                    '--use-conda',
                    '--cores', config['cores'],
                    '--resources disk=' + config['disk'], 'mem=' + config['resources_mem'],
                    keepgoing
                    ])

    # run snakemake
    sp.run(cmd, shell=True)
    
    # delete last batch input file
    sp.run('rm %spairs_batch' %processpath, shell=True)
    
    if config['custom_storepath']:
        storepath = config['custom_storepath']
        
        print('\nCopying results of batch %s/%s to %s\n' %(b//n+1, math.ceil(P/n), storepath))
        
        # copy results to alternative storepath
        sp.run('rsync -a --ignore-existing %s %s' %(processpath[:-1], storepath), shell=True)



    if config['delete_bam']:
        
        print('\nDeleting bqsr BAMs of batch %s/%s\n' %(b//n+1, math.ceil(P/n)))
        
        # delete .BAM, .BAI in processpath
        sp.run('rm %s03_alignment_genome/02_bqsr/*' %processpath, shell=True)
        sp.run('rm %s05_alignment_MT/02_bqsr/*' %processpath, shell=True)


    if config['delete_trim_logs']:
    
        print('\nDeleting compressed trim logs of batch %s/%s\n' %(b//n+1, math.ceil(P/n)))
        
        # delete gzip compressed trim logs in processpath
        sp.run('rm %s02_fastq_trimmed/logs/*.log.gz' %processpath, shell=True)



    if config['log_runs']:
        
        log_runs_path = homepath + 'ANALYSES_WES/log_runs.tsv'
        columns = ['run', 'tumor', 'normal', 'analysis']
        
        # load table or create DF if it doesn't exist
        try:
            log_runs = pd.read_table(log_runs_path, sep='\t', header=0, dtype='str')
        except:
            log_runs = pd.DataFrame(columns=columns)

        ### create log for current batch ###

        run = [config['outfolder']] * n
        tumor = [p.split('_')[0] for p in pairs_batch]
        normal = [p.split('_')[1] for p in pairs_batch]
        analysis = ['completed'] * n

        if os.path.exists(processpath + 'failed/'):
            failedpath = processpath + 'failed/'
            analysis = ['failed' if p in os.listdir(failedpath) else 'completed' for p in pairs_batch]
        
        df = pd.DataFrame({'run': run,
                          'tumor': tumor,
                          'normal': normal,
                          'analysis': analysis})
            
        # merge with full log
        log_runs = pd.concat([log_runs, df])
        
        # re-order columns
        log_runs = log_runs[columns]
        
        # drop duplicate rows (in case some analysis was repeated)
        log_runs.drop_duplicates(inplace=True)
        
        # write to file
        log_runs.to_csv(log_runs_path, sep='\t', index=False)


print('EVERY BATCH WERE FULLY ANALYZED - CHECK <failed> FOLDER FOR ANY FAILED SAMPLE PAIRS')

#%TODO: print batch in snakemake messages??
