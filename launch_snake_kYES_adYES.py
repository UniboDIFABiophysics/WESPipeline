import argparse
import pandas as pd
import subprocess as sp
import yaml
import os
import re


# get home path
homepath = os.path.expanduser('~') + '/'

# get current path
currentpath = os.getcwd() + '/'



#### Parse command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('-o', required=True, dest='outfolder', action='store', help='name of folder to store current analyses')
parser.add_argument('-i', required=True, dest='input_list', action='store', help='path to input samples list (no need to provide normal samples)')
parser.add_argument('-c', required=True, dest='configfile', action='store', help='path to config file')
parser.add_argument('-p', dest='custom_pair', action='store', help='path to custom TUMOR/NORMAL pairing file')
parser.add_argument('-s', dest='custom_storepath', action='store', help='alternative path to store results and BAM')

# get arguments
args = parser.parse_args()

# get OUTFOLDER
outfolder = args.outfolder

# build and create PROCESSPATH
processpath = homepath + 'ANALYSES_WES/' + outfolder + '/'
if not os.path.exists(processpath):
     os.makedirs(processpath)


### get SAMPLE LIST file absolute path
input_list = args.input_list
# in case a relative path was given as input
if not os.path.commonprefix([homepath, input_list]):
    # append to current path
    input_list = currentpath + input_list


### get CONFIG file absolute path
configfile = args.configfile
# in case a relative path was given as input
if not os.path.commonprefix([homepath, configfile]):
    # append to current path
    configfile = currentpath + configfile


### get absolute path to CUSTOM TUMOR/NORMAL pairING file
custom_pair = args.custom_pair
if custom_pair:
    # in case a relative path was given as input
    if not os.path.commonprefix([homepath, custom_pair]):
        # append to current path
        custom_pair = currentpath + custom_pair


# get alternative store path (if indicated in command line)
if args.custom_storepath:
    storepath = args.custom_storepath + '/'






#### Load data

## Load input sample list
with open(input_list) as infile:
    samples = infile.read().splitlines()
# move sample list file to processpath
sp.run(' '.join(['cp', input_list, processpath]), shell=True)


## load config file for WES
with open(configfile) as infile:
    config = yaml.load(infile)


# load dataset file and set sample names as indeces
data = homepath + config['dataset']
dataset = pd.read_table(data, sep = '\t', header=0, dtype='str')
dataset = dataset.set_index("sample_id")




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
    cmd = ' '.join(['ls', '-l', fastq_path + '*fastq.gz'])
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



##### BUILD TUMOR/NORMAL PAIRING LIST

# create empty list to store tumor_normal name pairs
pair_list = []
checked_pair_list = []

# if -p option was selected
if custom_pair:

    # load custom pairing file
    custom_pair = pd.read_table(custom_pair, sep='\t', header=0, dtype='str')

    # move custom pairing file to processpath
    sp.run(' '.join(['mv', custom_pair, processpath]), shell=True)


    # get list of tumors that are both in sample list and custom list
    tumor = [s for s in samples if s in custom_pair.index]

    # loop over them
    for t in tumor:

        # get normal samples that are paired to t in custom list
        normal = custom_pair.loc[t, 'normal']

        # if there is only one normal, concatenate to tumor and add to tumor_normal list
        if type(normal) == str:
            pair_list.append('_'.join([t, normal]))

        # if there are many normal, concatenate each of them to tumor and add to tumor_normal list
        else:
            for n in normal:
                pair_list.append('_'.join([t, n]))


# if no custom pairing is needed
else:

    # get list of tumor contained in sample list
    tumor = [s for s in samples if dataset.loc[s, 'pairing'] == 'tumor']

    # loop over them
    for t in tumor:

        # get patient ID
        patient = dataset.loc[t, 'patient_id']

        # get paired normal
        normal = dataset[(dataset.patient_id == patient) & (dataset.pairing == 'normal')].index[0]

        # concatenate it to tumor and add to tumor_normal list
        pair_list.append('_'.join([t, normal]))

# LOOP OVER pairs
for pair in pair_list:

        # extract names of current paired samples
        names = pair.split('_')

        # get WES library kit for tumor and normal
        kit_tumor = dataset.loc[names[0], 'kit_wes']
        kit_normal = dataset.loc[names[1], 'kit_wes']

        # KIT CHECKPOINT
        if kit_tumor != kit_normal:
            print('PAIR %s - TUMOR and NORMAL have different WES kits!!!' %pair)
            print('ABORTING')
            break

        # append to growing list
        checked_pair_list.append(pair)

# do unique, to get final list of pairs
pairs = list(set(checked_pair_list))

# split every pair and get list of every sample (both tumor and normal)
samples = [s for p in pairs for s in p.split('_')]




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



### SPLIT THE ORDERED LIST OF PAIRS INTO BATCHES OF 6 PAIRS
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
                    '--config', 'outfolder=' + outfolder, 'adapters=True',
                    '--use-conda',
                    '--cores', config['cores'],
                    '--keep-going'
                   ])

    # run snakemake
    sp.run(cmd, shell=True)

    # delete last batch input file
    sp.run('rm %spairs_batch' %processpath, shell=True)

    # copy results to alternative storepath
    sp.run('rsync -a --ignore-existing %s %s' %(processpath[:-1], storepath), shell=True)
#%TODO: print message

    # delete BAM/BAI in processpath
    sp.run('rm %s03_alignment_genome/02_bqsr/*' %processpath, shell=True)
    sp.run('rm %s05_alignment_MT/02_bqsr/*' %processpath, shell=True)
