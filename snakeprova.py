
shell.executable("/bin/bash")

#configfile: False

import pandas as pd
import os
import yaml
import subprocess as sp


# get main path (home)
homepath = os.path.expanduser('~') + '/'

# # set environments path
# pipelinepath = homepath + 'WESPipeline/'
#
# # get current path
# #currentpath = os.getcwd()
#
outfolder = config['outfolder']

# custom_storepath = config['custom_storepath']
#
# scripts = config['folders']['scripts']
#
# build and create PROCESSPATH & STOREPATH
processpath = homepath + 'ANALYSES_WES/' + outfolder + '/'

if not os.path.exists(processpath):
     os.makedirs(processpath)

# storepath = currentpath
# storepath_tree = storepath
#
# # save tree structure as dictionary
# dirs = {'01_fastq': ['01_preprocess', '02_postprocess', '03_fastqc', 'logs'],
#           '02_fastq_trimmed': ['01_fastqc', 'logs'],
#           '03_alignment_genome': ['01_intermediate', '02_bqsr', 'logs'],
#           '04_alignment_exome': ['01_intermediate', '02_unmapped_fastq', 'logs'],
#           '05_alignment_MT': ['01_intermediate', '02_bqsr', 'logs'],
#           '06_mutect_genome': ['logs'],
#           '07_mutect_MT': ['logs'],
#           '08_varscan_genome': ['01_mpileup', 'logs'],
#           '09_varscan_MT': ['01_mpileup', 'logs'],
#           '10_variants_filter': ['logs'],
#          }
#
# # loop over dictionary items and define each folder and subfolder
# processpath_tree = [['/'.join([processpath, d, d2,'']) for d2 in dirs[d]] for d in dirs]

# # initialize storepath if necessary
# if custom_storepath:
#     storepath = custom_storepath + outfolder
#     # storepath_tree = [['/'.join([storepath, d, d2,'']) for d2 in dirs[d]] for d in dirs]

#/-----------------------------------------------------------------------------------------#/


# # References
# hg = homepath + config['fasta']['genome'] # Human Genome Reference
# MT = homepath + config['fasta']['MT']
# nextera = homepath + config['fasta']['nextera']
# nextera_expanded = homepath + config['fasta']['nextera_expanded']
# truseq = homepath + config['fasta']['truseq']
# truseq_rapid = homepath + config['fasta']['truseq_rapid']
#
# hg_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"] # Indexes for hg
# MT_indexes = [MT+".bwt", MT+".pac", MT+".amb", MT+".ann", MT+".sa"]
# nextera_indexes = [nextera+".bwt", nextera+".pac", nextera+".amb", nextera+".ann", nextera+".sa"]
# nextexp_indexes = [nextera_expanded+".bwt", nextera_expanded+".pac", nextera_expanded+".amb", nextera_expanded+".ann", nextera_expanded+".sa"]
# truseq_indexes = [truseq+".bwt", truseq+".pac", truseq+".amb", truseq+".ann", truseq+".sa"]
# truseq_rapid_indexes = [truseq_rapid+".bwt", truseq_rapid+".pac", truseq_rapid+".amb", truseq_rapid+".ann", truseq_rapid+".sa"]
#
# all_fastq = homepath + config['all_fastq']
# cadaver = config['cadaver']
#
# if not cadaver:
#     all_fastq = "/"
#
# indels_ref = homepath+ config['ref-files']['indels_ref'] # Set of known indels
# dbsnp = homepath+ config['ref-files']['dbsnp'] # SNP database
# cosmic = homepath+ config['ref-files']['cosmic'] # Catalog of somatic mutation in cancer
# humandb = homepath+ config['ref-files']['humandb'] # Annovar human databases folder
# build_ver = config['ref-files']['build_ver'] # Set build version
# dbsnp_ver = config['ref-files']['dbsnp_ver'] # Set SNP database version
# mitochondrial_ver = config['ref-files']['mitochondrial_ver'] # Set parameter command for annotating mitochondria variants
#
# # Softwares
# gatk = homepath+ config['softwares']['gatk']
# muTect = homepath+ config['softwares']['muTect']
# annovar = homepath + config['folders']['annovar']
# annotate = homepath+ config['softwares']['annotate']
# tableannovar = homepath+ config['softwares']['tableannovar']
# adapter_removal = homepath + config['softwares']['adapter_removal']
# softwares = homepath + config['folders']['softwares']
# picard = homepath + config['softwares']['picard']
# varscan = homepath + config['softwares']['varscan']
#
# # Sample details
# platform = config['sample-details']['platform'] # Set platform for alignment
#
# # bed files
# nextexp_bed = homepath + config['bed']['nextera_expanded'] # Set target intervals for exome analysis
# nextera_bed = homepath + config['bed']['nextera']
# MT_bed = homepath + config['bed']['MT']
# truseq_bed = homepath + config['bed']['truseq']
# truseq_rapid_bed = homepath + config['bed']['truseq_rapid']
#
# # Annovar databases
# annovar_dbs = [homepath + config['annovar_dbs']['hg19_refGene'],
#                homepath + config['annovar_dbs']['hg19_cytoBand'],
#                homepath + config['annovar_dbs']['hg19_gSd'],
#                homepath + config['annovar_dbs']['hg19_esp'],
#                homepath + config['annovar_dbs']['hg19_snp138'],
#                homepath + config['annovar_dbs']['hg19_1000g2014oct'],
#                homepath + config['annovar_dbs']['hg19_exac03nontcga'],
#                homepath + config['annovar_dbs']['hg19_ljb26_all'],
#                ]
#
#
# # script rsIDfilter
#
# rsIDfilter = pipelinepath + 'scripts/' + config['rsIDfilter']
#
#
# # Label's parameters
# n_sim = config['n_sim']
# cpu_type = config['cpu_type']
# thrs = config['threads']
# n_cpu = config['n_cpu']
#
# # Max threads per each multithreading rule
# map_thrs = config['map_thrs']
# RT_thrs = config['RT_thrs']
# BaseRecal_thrs = config['BaseRecal_thrs']
# PrintReads_thrs = config['PrintReads_thrs']
#
#
# # choose adapters for trimming
# adapters = config['adapters']
#
# if adapters:
#     pcr1 = ' --pcr1 ' + config['nextera_transposase_2_rc'],
#     pcr2 = ' --pcr2 ' + config['nextera_transposase_1'],
# else:
#     pcr1 = ''
#     pcr2 = ''
#

#### Load data

## Load input sample list
input_list = config['input_list']

with open(input_list) as infile:
    sample_names = infile.read().splitlines()

# # move sample list file to processpath
# sp.call(' '.join(['cp', input_list, processpath]), shell=True)


data = homepath + config['dataset']
dataset = pd.read_table(data, sep = '\t',header=0, dtype='str')
dataset = dataset.set_index("sample_id")

##### BUILD TUMOR/NORMAL MATCHING LIST

# create empty list to store tumor_normal name matches
#match_list = {}
#checked_match_list = {}
match_list = []
checked_match_list = []

custom_match = config['custom_match']

# if -m option was selected
if custom_match:

    # load custom matching file
    custom_match = pd.read_table(custom_match, sep='\t', header=0, dtype='str')

    # move custom matching file to processpath
    sp.call(' '.join(['mv', custom_match, processpath]), shell=True)


    # get list of tumors that are both in sample list and custom list
    tumor = [s for s in sample_names if s in custom_match.index]

    # loop over them
    for t in tumor:

        # get normal samples that are matched to t in custom list
        normal = custom_match.loc[t, 'normal']

        # if there is only one normal, concatenate to tumor and add to tumor_normal list
        if type(normal) == str:
            match_list.append('_'.join([t, normal]))

        # if there are many normal, concatenate each of them to tumor and add to tumor_normal list
        else:
            for n in normal:
                match_list.append('_'.join([t, n]))


# if no custom matching is needed
else:

    # get list of tumor contained in sample list
    tumor = [s for s in sample_names if dataset.loc[s, 'matching'] == 'tumor']

    # loop over them
    for t in tumor:

        # get patient ID
        patient = dataset.loc[t, 'patient_id']
        # get matched normal
        normal = dataset[(dataset.patient_id == patient) & (dataset.matching == 'normal')].index[0]

        # concatenate it to tumor and add to tumor_normal list
        match_list.append('_'.join([t, normal]))
#        match_list[patient] = '_'.join([t, normal])

### LOOP OVER MATCHES
patients_dict={}
for patient in match_list:
#        match = match_list[patient]
        # extract names of current paired samples
        names = patient.split('_')

        # get WES library kit for tumor and normal
        kit_tumor = dataset.loc[names[0], 'kit_wes']
        kit_normal = dataset.loc[names[1], 'kit_wes']

        # KIT CHECKPOINT
        if kit_tumor != kit_normal:
            continue

#        checked_match_list[patient] = match_list[patient]
        checked_match_list.append(patient)
#        patients_dict[patient] = match_list[patient]
        patients_dict[patient] = {}
        patients_dict[patient]['T'] = names[0]
        patients_dict[patient]['N'] = names[1]
        patients_dict[patient]['kit'] = kit_tumor

patients  = [p for p in checked_match_list]
patients = list(set(patients))
#samples = [s for p in patients for s in checked_match_list[p].split('_')]
samples = [s for p in checked_match_list for s in p.split('_')]
samples = list(set(samples))

# print(checked_match_list)
# print(patients_dict)
# print(patients)
# print(samples)

# def get_kit(wildcards):
#     return dataset.loc[wildcards, 'kit_wes']
#
# def get_bed_patient(wildcards):
#     kit = patients_dict[wildcards]['kit']
#     return (homepath + config['bed'][kit] + "_fixed.bed")
#
# def get_fastq_path(wildcards):
#     path = dataset.loc[wildcards,'fastq_path_wes']
#     if not cadaver:
#          path = homepath + '/' + path
#     if path == './':
#         path = ''
#     return path
#
# def get_fastq_id(wildcards):
#     return dataset.loc[wildcards,'fastq_id_wes']

# # get pcr primers
# def get_adapter(wildcards,i):
#     kit = get_kit(wildcards)
#     if i==1:
#         return config['adapter'][kit]['adapter1']
#     else:
#         return config['adapter'][kit]['adapter2']

# def get_bed(wildcards):
#     kit = get_kit(wildcards)
#     return homepath + config['bed'][kit] + "_fixed.bed"
#
# def get_ref(wildcards):
#     kit = get_kit(wildcards)
#     return homepath + config['fasta'][kit]
#
# def get_ref_indexes(wildcards):
#     kit = get_kit(wildcards)
#     r = homepath + config['fasta'][kit]
#     r_indexes = [r+".bwt", r+".pac", r+".amb", r+".ann", r+".sa"]
#     return r_indexes
#
# def get_ref_fai(wildcards):
#     kit = get_kit(wildcards)
#     r = homepath + config['fasta'][kit]
#     return r+'.fai'
#
# def get_ref_dict(wildcards):
#     kit = get_kit(wildcards)
#     r = homepath + config['fasta'][kit]
#     return r.replace('fasta','dict')
#
def get_bam(wildcards, sample_type, Dir):
    return (Dir + patients_dict[wildcards][sample_type] + ".bam")

# def get_bai(wildcards, sample_type, Dir):
#     return (Dir + patients_dict[wildcards][sample_type] + ".bai")
#
# def get_mpileup(wildcards, sample_type, Dir):
#     return (Dir + patients_dict[wildcards][sample_type] + ".mpileup")


# # define directories
# tmp_dir = processpath + 'tmp/'
# if not os.path.exists(tmp_dir):
#      os.makedirs(tmp_dir)
# preprocesspath = processpath + '01_fastq/01_preprocess/'
# postprocesspath = processpath + '01_fastq/02_postprocess/'
# trimmedpath = processpath + '02_fastq_trimmed/'
# fastqcpath = processpath + '01_fastq/03_fastqc/'
# fastqcpath_trim = processpath + '02_fastq_trimmed/01_fastqc/'
#
# alignpath_genome = processpath + '03_alignment_genome/'
# alignpath_exome = processpath + '04_alignment_exome/'
# alignpath_MT = processpath + '05_alignment_MT/'
#
# alignpath_genomebqsr = alignpath_genome + "02_bqsr/"
# alignpath_MTbqsr = alignpath_MT + "02_bqsr/"
# alignpath_exomeunmapped = alignpath_exome + '02_unmapped_fastq/'
# fastq_logs = processpath + '01_fastq/logs/'
# trimmed_lods = processpath + '02_fastq_trimmed/logs/'
# genome_int = alignpath_genome + '01_intermediate/'
# genome_logs = alignpath_genome + 'logs/'
# exome_int = alignpath_exome + '01_intermediate/'
# exome_logs = alignpath_exome + 'logs/'
# MT_int = alignpath_MT + '01_intermediate/'
# MT_logs = alignpath_MT + 'logs/'
#
#
# mutectpath_genome = processpath + '06_mutect_genome/'
# mutectpath_genome_logs = mutectpath_genome + 'logs/'
# mutectpath_MT = processpath + '07_mutect_MT/'
# mutectpath_MT_logs = mutectpath_MT + 'logs/'
#
# varscanpath_genome = processpath + '08_varscan_genome/'
# varscanpath_MT = processpath + '09_varscan_MT/'
# mpileup_varscan_genome = varscanpath_genome + '01_mpileup/'
# mpileup_varscan_MT = varscanpath_MT + '01_mpileup/'
# varscan_genome_logs = varscanpath_genome + 'logs/'
# varscan_MT_logs = varscanpath_MT + 'logs/'
#
# variants_filter = processpath + '10_variants_filter/'
# variants_filter_logs = variants_filter + 'logs/'


# print("("+"|".join(samples)+")")
# print("("+"|".join(patients)+")")

# samples = ['12', '34', '56', '78']
# patients = ['a', 'b', 'c']

# Wildcard costrains necessary for search only certain names
wildcard_constraints:
    sample = "("+"|".join(samples)+")",
    patient = "("+"|".join(patients)+")",

#os.makedirs(processpath + 'campioni/', exist_ok=True)
os.makedirs(processpath + 'persone/', exist_ok=True)

rule all:
    input:
        expand(processpath + '{sample}.bam', sample=samples),
        expand(processpath + 'persone/{patient}.tsv', patient=patients),



rule one:
    output:
        processpath + '{sample}.bam'
    run:
        for sa in samples:
            shell('echo %s > %s%s.bam' %(sa, processpath, sa))



rule two:
    input:
        tumor = lambda wildcards: get_bam(wildcards.patient,'T', processpath),
        normal = lambda wildcards: get_bam(wildcards.patient,'N', processpath)
    output:
        processpath + 'persone/{patient}.tsv'
    shell:
        'cat {input.tumor} {input.normal} > {output}'




# rule all:
#   """
#     PIPELINE ENDING
#   """
#   input:
#     expand(fastqcpath +"{sample}" + "_R1_fastqc.html",sample=samples),
#     expand(fastqcpath +"{sample}" + "_R2_fastqc.html",sample=samples),
#     expand(fastqcpath_trim +"{sample}" + "_R1_fastqc.html",sample=samples),
#     expand(fastqcpath_trim +"{sample}" + "_R2_fastqc.html",sample=samples),
#     expand(fastqcpath +"{sample}" + "_R1_fastqc.zip",sample=samples),
#     expand(fastqcpath +"{sample}" + "_R2_fastqc.zip",sample=samples),
#     expand(fastqcpath_trim +"{sample}" + "_R1_fastqc.zip",sample=samples),
#     expand(fastqcpath_trim +"{sample}" + "_R2_fastqc.zip",sample=samples),
#     expand(genome_logs + "{sample}"+"_recalibrationPlots.pdf",sample=samples),
#     expand(MT_logs + "{sample}"+"_recalibrationPlots.pdf",sample=samples),
#     expand(variants_filter + "{patient}_all_somatic_annotated_filtered_final.tsv",patient=patients),
#   message: " THE END "
#   run:
#     pass
