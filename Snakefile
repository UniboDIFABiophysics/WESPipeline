
shell.executable("/bin/bash")

configfile: "wes_config.yaml"

import pandas as pd
import os
import yaml
import subprocess as sp


# get main path (home)
homepath = os.path.expanduser('~')

# get current path
currentpath = os.getcwd()

folder= config['folders']['resultdir']
scripts =  config['folders']['scripts']

custom_storepath=config['folders']['custom_resultdir']


# build and create PROCESSPATH & STOREPATH
processpath = currentpath + '/wes_analyses' + folder

if not os.path.exists(processpath):
     os.makedirs(processpath)

storepath = currentpath
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

# initialize storepath if necessary
if custom_storepath:
    storepath = currentpath + custom_storepath + folder
    processpath = storepath
    # storepath_tree = [['/'.join([storepath, d, d2,'']) for d2 in dirs[d]] for d in dirs]

#/-----------------------------------------------------------------------------------------#/


# References
hg = homepath + config['fasta']['genome'] # Human Genome Reference
MT = homepath + config['fasta']['MT']
nextera = homepath + config['fasta']['nextera']
nextera_expanded = homepath + config['fasta']['nextera_expanded']
truseq = homepath + config['fasta']['truseq']

hg_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"] # Indexes for hg
MT_indexes = [MT+".bwt", MT+".pac", MT+".amb", MT+".ann", MT+".sa"]
nextera_indexes = [nextera+".bwt", nextera+".pac", nextera+".amb", nextera+".ann", nextera+".sa"]
nextexp_indexes = [nextera_expanded+".bwt", nextera_expanded+".pac", nextera_expanded+".amb", nextera_expanded+".ann", nextera_expanded+".sa"]
truseq_indexes = [truseq+".bwt", truseq+".pac", truseq+".amb", truseq+".ann", truseq+".sa"]

all_fastq = homepath + config['all_fastq']
cadaver = config['cadaver']

if not cadaver:
    all_fastq = "/"

indels_ref = homepath+ config['ref-files']['indels_ref'] # Set of known indels
dbsnp = homepath+ config['ref-files']['dbsnp'] # SNP database
cosmic = homepath+ config['ref-files']['cosmic'] # Catalog of somatic mutation in cancer
humandb = homepath+ config['ref-files']['humandb'] # Annovar human databases folder
build_ver = config['ref-files']['build_ver'] # Set build version
dbsnp_ver = config['ref-files']['dbsnp_ver'] # Set SNP database version
kg_ver = config['ref-files']['kg_ver'] # Set  version # Set ethnicity groups based on time
mitochondrial_ver = config['ref-files']['mitochondrial_ver'] # Set parameter command for annotating mitochondria variants

# Softwares
gatk = homepath+ config['softwares']['gatk']
muTect = homepath+ config['softwares']['muTect']
annovar = homepath+ config['softwares']['annovar']
convert2annovar = homepath+ config['softwares']['convert2annovar']
tableannovar = homepath+ config['softwares']['tableannovar']
adapter_removal = homepath + config['softwares']['adapter_removal']

# Sample details
platform = config['sample-details']['platform'] # Set platform for alignment
library = config['sample-details']['library'] # Set library for alignment

# bed files
nextexp_bed = homepath + config['bed']['nextera_expanded'] # Set target intervals for exome analysis
nextera_bed = homepath + config['bed']['nextera']
MT_bed = homepath + config['bed']['MT']
truseq_bed = homepath + config['bed']['truseq']

# Annovar databases
annovar_dbs = [homepath + config['annovar_dbs']['hg19_db'] , homepath + config['annovar_dbs']['snp138_db']]

# script rsIDfilter

rsIDfilter = scripts + config['rsIDfilter']


# Label's parameters
n_sim = config['n_sim']
cpu_type = config['cpu_type']
thrs = config['threads']
n_cpu = config['n_cpu']

# Max threads per each multithreading rule
max_thrs = int(thrs)
map_thrs = max_thrs
RT_thrs = max_thrs
BQSR1_thrs = max_thrs
BQSR2_thrs = max_thrs
BQSR4_thrs = 4
HC_thrs = 4
HF1_thrs = 1
HF2_thrs = 1
HFC_thrs = 1
M_thrs = 1
LodnT_thrs = 1
Lodnv_thrs = 1





#### Load data

## Load input sample list
input_list = config['input_list']

with open(input_list) as infile:
    sample_names = infile.read().splitlines()

# move sample list file to processpath
sp.call(' '.join(['cp', input_list, processpath]), shell=True)


data = currentpath + config['dataset']
dataset = pd.read_table(data, sep = '\t',header=0, index_col='sample_id', dtype='str')

##### BUILD TUMOR/NORMAL MATCHING LIST

# create empty list to store tumor_normal name matches
match_list = {}
checked_match_list = {}

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
        match_list[patient] = ('_'.join([t, normal]))

### LOOP OVER MATCHES
patients_dict={}
for patient in match_list:
        match = match_list[patient]
        # extract names of current paired samples
        names = match.split('_')

        # get WES library kit for tumor and normal
        kit_tumor = dataset.loc[names[0], 'kit_wes']
        kit_normal = dataset.loc[names[1], 'kit_wes']

        # KIT CHECKPOINT
        if kit_tumor != kit_normal:
            continue

        checked_match_list[patient] = match_list[patient]
        patients_dict[patient] = match_list[patient]
        patients_dict[patient] = {}
        patients_dict[patient]['T'] = names[0]
        patients_dict[patient]['N'] = names[1]
        patients_dict[patient]['kit'] = kit_tumor

patients  = [p for p in checked_match_list]
samples = [s for p in patients for s in checked_match_list[p].split('_')]

def get_kit(wildcards):
    return dataset.loc[wildcards, 'kit_wes']

def get_bed_patient(wildcards):
    kit = patients_dict[wildcards]['kit']
    return (homepath + config['bed'][kit] + "_fixed.bed")

def get_fastq_path(wildcards):
    path = dataset.loc[wildcards,'fastq_path_wes']
    if path == './':
        path = ''
    return path

def get_fastq_id(wildcards):
    return dataset.loc[wildcards,'fastq_id_wes']

# get pcr primers
def get_pcr(wildcards,i):
    kit = get_kit(wildcards)
    if i==1:
        return config['adapter'][kit]['adapter1']
    else:
        return config['adapter'][kit]['adapter2']

def get_bed(wildcards):
    kit = get_kit(wildcards)
    return homepath + config['bed'][kit] + "_fixed.bed"

def get_ref(wildcards):
    kit = get_kit(wildcards)
    return homepath + config['fasta'][kit]

def get_ref_indexes(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    r_indexes = [r+".bwt", r+".pac", r+".amb", r+".ann", r+".sa"]
    return r_indexes

def get_ref_fai(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    return r+'.fai'

def get_ref_dict(wildcards):
    kit = get_kit(wildcards)
    r = homepath + config['fasta'][kit]
    return r.replace('fasta','dict')

def get_bam(wildcards, sample_type, dir):
    return (dir + patients_dict[wildcards][sample_type] + "_recal.bam")

def get_bai(wildcards, sample_type, dir):
    return (dir + patients_dict[wildcards][sample_type] + "_recal.bai")

def get_mpileup(wildcards, sample_type, dir):
    return (dir + patients_dict[wildcards][sample_type] + ".mpileup")


# define align paths

alignpath_genome = processpath + '/03_alignment_genome/'
alignpath_exome = processpath + '/04_alignment_exome/'
alignpath_MT = processpath + '/05_alignment_MT/'


# define directories
tmp_dir_genome = currentpath + '/wes_analyses/tmp/genome/'
tmp_dir_exome = currentpath + '/wes_analyses/tmp/exome/'
tmp_dir_MT = currentpath + '/wes_analyses/tmp/MT/'
tmp_dir_genome_m = currentpath + '/wes_analyses/tmp/genome/m/'
tmp_dir_MT_m = currentpath + '/wes_analyses/tmp/MT/m/'
postprocesspath = processpath + '/01_fastq/02_postprocess/'
trimmedpath = processpath + '/02_fastq_trimmed/'
fastqcpath = processpath + '/01_fastq/03_fastqc/'
fastqcpath_trim = processpath + '/02_fastq_trimmed/01_fastqc/'
alignpath_genomebqsr = alignpath_genome + "02_bqsr/"
alignpath_MTbqsr = alignpath_MT + "02_bqsr/"
alignpath_exomeunmapped = alignpath_exome + '02_unmapped_fastq/'
fastq_logs = processpath + '/01_fastq/logs/'
trimmed_lods = processpath + '/02_fastq_trimmed/logs/'
genome_int = alignpath_genome + '01_intermediate/'
genome_logs = alignpath_genome + 'logs/'
exome_int = alignpath_exome + '01_intermediate/'
exome_logs = alignpath_exome + 'logs/'
MT_int = alignpath_MT + '01_intermediate/'
MT_logs = alignpath_MT + 'logs/'

mutectpath_genome = processpath + '/06_mutect_genome/'
mutectpath_genome_logs = mutectpath_genome + 'logs/'
mutectpath_MT = processpath + '/07_mutect_MT/'
mutectpath_MT_logs = mutectpath_MT + 'logs/'

varscanpath_genome = processpath + '/08_varscan_genome/'
varscanpath_MT = processpath + '/09_varscan_MT/'
mpileup_varscan_genome = varscanpath_genome + '01_mpileup/'
mpileup_varscan_MT = varscanpath_MT + '01_mpileup/'
varscan_genome_logs = varscanpath_genome + 'logs/'
varscan_MT_logs = varscanpath_MT + 'logs/'

variants_filter = processpath + '/10_variants_filter/'
variants_filter_logs = variants_filter + 'logs/'

# Wildcard costrains necessary for search only certain names
wildcard_constraints:
    sample = "("+"|".join(samples)+")",
    patient = "("+"|".join(patients)+")",



##########################################
#           PIPELINE BEGINNING           #
##########################################


rule all:
  """
    PIPELINE ENDING
  """
  input:
    expand(fastqcpath +"{sample}" + "_R1_fastqc.html",sample=samples),
    expand(fastqcpath +"{sample}" + "_R2_fastqc.html",sample=samples),
    expand(fastqcpath_trim +"{sample}" + "_R1_fastqc.html",sample=samples),
    expand(fastqcpath_trim +"{sample}" + "_R2_fastqc.html",sample=samples),
    expand(fastqcpath +"{sample}" + "_R1_fastqc.zip",sample=samples),
    expand(fastqcpath +"{sample}" + "_R2_fastqc.zip",sample=samples),
    expand(fastqcpath_trim +"{sample}" + "_R1_fastqc.zip",sample=samples),
    expand(fastqcpath_trim +"{sample}" + "_R2_fastqc.zip",sample=samples),
    expand(variants_filter + "{patient}_all_somatic_annotated_filtered.tsv",patient=patients),
  run:
    pass

###########################################################################################

rule create_tmpdir:
  output:
    temp(tmp_dir_genome),
    temp(tmp_dir_exome),
    temp(tmp_dir_MT),
    temp(tmp_dir_genome_m),
    temp(tmp_dir_MT_m),
  threads: 1
  run:
    pass

#% TODO: add temp for outpath[1,2] after checking

rule downloadDecompressMergeFastq:
    input:
        all_fastq = all_fastq,
    output:
        outpath1 = postprocesspath+ "{sample}" + "_R1_unchecked.fastq",
        outpath2 = postprocesspath+ "{sample}" + "_R2_unchecked.fastq",
    params:
        fastq_path = lambda wildcards: get_fastq_path(wildcards.sample),
        fastq_id = lambda wildcards: get_fastq_id(wildcards.sample),
        scripts = scripts,
        name = "{sample}",
        homepath = homepath,
        processpath = processpath,
        cadaver = cadaver,
        logpath = fastq_logs,
    threads: 1
    log:
         fastq_logs + '{sample}_cadaver.log',
    script:
        "{params.scripts}"+"downloadDecompressMergeFastq.py"

rule fastq_checkpoint:
    input:
        unchecked1 = postprocesspath+"{sample}" + "_R1_unchecked.fastq",
        unchecked2 = postprocesspath+"{sample}" + "_R2_unchecked.fastq",
    output:
        checked1 = postprocesspath+"{sample}" + "_R1.fastq",
        checked2 = postprocesspath+"{sample}" + "_R2.fastq",
    params:
        name="{sample}",
    run:
        import os
        import subprocess as sp

        size1 = os.stat(input.unchecked1).st_size
        size2 = os.stat(input.unchecked2).st_size
        if size1 != size2:
            sp.call("echo 'Error. "+ params.name +" fastq files have different size!!'",shell=True)
            raise ErrorValue
        else:
            sp.call("mv "+ input.unchecked1 + " " + output.checked1 ,shell=True)
            sp.call("mv "+ input.unchecked2 + " " + output.checked2 ,shell=True)

rule fastqc_R1:
  input:
    fastq1 = postprocesspath+"{sample}" + "_R1.fastq",
  output:
    fastq1_fastqc_zip = fastqcpath +"{sample}" + "_R1_fastqc.zip",
    fastq1_fastqc_html = fastqcpath +"{sample}" + "_R1_fastqc.html",
  params:
    outpath = fastqcpath,
  log:
    fastq_logs + "{sample}" + '_R1_fastqc.log'
  conda:
    "envs/wes_config_conda.yaml"
  #benchmark:
  message: 'fastqc1'
  shell:
    "fastqc -o {params.outpath} {input.fastq1} > {log} 2>&1"

rule fastqc_R2:
  input:
    fastq2 = postprocesspath+"{sample}" + "_R2.fastq",
  output:
    fastq2_fastqc_zip = fastqcpath +"{sample}" + "_R2_fastqc.zip",
    fastq2_fastqc_html = fastqcpath +"{sample}" + "_R2_fastqc.html",
  params:
    outpath = fastqcpath,
  log:
    fastq_logs + "{sample}" + '_R2_fastqc.log'
  conda:
    "envs/wes_config_conda.yaml"
  #benchmark:
  message: 'fastqc2'
  shell:
    "fastqc -o {params.outpath} {input.fastq2} > {log} 2>&1"

#% TODO: Adapter removal issue

rule trim_02:
    input:
        fastq1 = postprocesspath+"{sample}" + "_R1.fastq",
        fastq2 = postprocesspath+"{sample}" + "_R2.fastq",
    output:
        fastq1_trimmed = trimmedpath+"{sample}" + "_R1.fastq",
        fastq2_trimmed = trimmedpath+"{sample}" + "_R2.fastq",
    params:
        pcr1 = lambda wildcards: get_pcr(wildcards.sample,1),
        pcr2 = lambda wildcards: get_pcr(wildcards.sample,2),
        adapter_removal=adapter_removal,
    log:
        log_disc = trimmed_lods + '{sample}_discarded.log',
        log_stats = trimmed_lods + '{sample}_stats.log',
        log_sgtn = trimmed_lods + '{sample}_singleton.log',
    #conda:
    #benchmarks:
    message: 'trimming'
    shell:
        "{params.adapter_removal} --file1 {input.fastq1} --file2 {input.fastq1} --pcr1 {params.pcr1} --pcr2 {params.pcr2} --stats --trimns --trimqualities --minquality 20 --minlength 80 --output1 {output.fastq1_trimmed} --output2 {output.fastq2_trimmed} --discarded {log.log_disc} --outputstats {log.log_stats} --singleton {log.log_sgtn}"


rule fastqc_trimmed_R1:
  input:
    fastq1_trimmed = trimmedpath+"{sample}" + "_R1.fastq",
  output:
    fastq1tr_fastqc_zip = fastqcpath_trim+"{sample}" + "_R1_fastqc.zip",
    fastq1tr_fastqc_html = fastqcpath_trim+"{sample}" + "_R1_fastqc.html",
  params:
    outpath = fastqcpath_trim,
  log:
    trimmed_lods + "{sample}" + '_R1_fastqc.log'
  conda:
    "envs/wes_config_conda.yaml"
  #benchmark:
  message: 'fastqc on trimmed fastq1'
  shell:
    "fastqc -o {params.outpath} {input} > {log} 2>&1"

rule fastqc_trimmed_R2:
  input:
    fastq2_trimmed = trimmedpath+"{sample}" + "_R2.fastq",
  output:
    fastq2tr_fastqc_zip = fastqcpath_trim+"{sample}" + "_R2_fastqc.zip",
    fastq2tr_fastqc_html = fastqcpath_trim+"{sample}" + "_R2_fastqc.html",
  params:
    outpath = fastqcpath_trim,
  log:
    trimmed_lods + "{sample}" + '_R2_fastqc.log'
  conda:
    "envs/wes_config_conda.yaml"
  #benchmark:
  message: 'fastqc on trimmed fastq1'
  shell:
    "fastqc -o {params.outpath} {input} > {log} 2>&1"



##########################
#### GENOME ALIGNMENT ####
##########################

rule map_to_genome:
    """
    This tool maps the samples to the human reference genome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        hg_indexes,
        fastq1_trimmed = trimmedpath+"{sample}" + "_R1.fastq",
        fastq2_trimmed = trimmedpath+"{sample}" + "_R2.fastq",
    output:
        outfile =  genome_int + "{sample}.sam",
    params:
        reference = hg,
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = "{sample}",
    log:
        genome_logs + '{sample}_alignment.log',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_mapping_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    threads: map_thrs
    #resources: mem=6
    #version: 0.1
    message: ">> {wildcards.sample} - Aligning to reference NUCLEAR GENOME"
    shell:
        "bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {params.reference} {input.fastq1_trimmed} {input.fastq2_trimmed} > {output.outfile} 2> {log}"

rule sorting_genome:
    """
    This tool sorts the input SAM or BAM file by coordinate.
    """
    input:
        sam = genome_int + "{sample}.sam",
        tmp = tmp_dir_genome,
    output:
        outdir = genome_int + "{sample}_sorted.bam",
    log:
        alignpath_genome + '{sample}_sorting.log'
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
     #   "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard SortSam INPUT={input.sam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={input.tmp} 2> {log}"

rule marking_genome:
    """
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
    """
    input:
        sorted_bam = genome_int + "{sample}_sorted.bam",
        tmp = tmp_dir_genome_m,
    output:
        out = genome_int+"{sample}"+"_marked.bam",
    log:
        mx = alignpath_genome + '{sample}_metrix.log',
        mark = alignpath_genome + '{sample}_marking.log',
    conda:
        "envs/wes_config_conda.yaml"
   # benchmark:
    #    "benchmarks/benchmark_mark_duplicates_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard MarkDuplicates"
        " INPUT={input.sorted_bam} OUTPUT={output.out} METRICS_FILE={log.mx} "
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR={input.tmp} 2> {log.mark}"

rule indexing_genome:
    """
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    """
    input:
        marked_bam = genome_int+"{sample}"+"_marked.bam",
    output:
        marked_bai = genome_int+"{sample}"+"_marked.bai",
    log:
        alignpath_genome + '{sample}_indexing.log',
    conda:
        "envs/wes_config_conda.yaml"
   # benchmark:
    #    "benchmarks/benchmark_build_bam_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard BuildBamIndex INPUT={input.marked_bam} OUTPUT={output} 2> {log}"


rule RTC_genome:
    """
    This tool defines intervals to target for local realignment.
    """
    input:
        hg.replace('fasta', 'dict'),
        marked_bai = genome_int+"{sample}"+"_marked.bai",
        ref=hg+'.fai',
        indels_ref=indels_ref,
        gatk = gatk,
        marked_bam = genome_int+"{sample}"+"_marked.bam",
        bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        genome_logs+"{sample}"+".intervals",
    params:
        ref=hg,
    log:
        genome_logs + "{sample}" + "_RTC.log"
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
     #   "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{RT_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, RT_thrs=RT_thrs, n_cpu=n_cpu)
    threads: RT_thrs
    shell:
        "java -jar {input.gatk} -T RealignerTargetCreator -R {params.ref} -I {input.marked_bam} -L {input.bed} -ip 50 -known {input.indels_ref} -nt {threads} -o {output} 2> {log}"


rule IndelRealigner_genome:
    """
    This tool performs local realignment of reads around indels.
    """
    input:
        intvs = genome_logs+"{sample}"+".intervals",
        bam = genome_int+"{sample}"+"_marked.bam",
        idx = genome_int+"{sample}"+"_marked.bai",
    output:
        r_bam = genome_int + "{sample}"+"_realigned.bam",
        r_idx = genome_int + "{sample}"+"_realigned.bai",
    params:
        gatk = gatk,
        ref = hg,
        indels_ref=indels_ref,
    log:
        genome_logs + "{sample}" + "_IndelRealigner.log"
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
     #   "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.intvs} -known {params.indels_ref} -ip 50 -o {output.r_bam} 2> {log}"


rule BaseRecal_genome:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    """
    input:
        r_bam = genome_int +"{sample}"+"_realigned.bam",
        r_idx = genome_int +"{sample}"+"_realigned.bai",
        dbsnp = dbsnp,
        bed = lambda wildcards: get_bed(wildcards.sample),
    output:
        outtable = genome_logs + "{sample}"+"_recal_data.table",
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    log:
        genome_logs + "{sample}" + '_recalibrating_01.log'
    conda:
        "envs/wes_config_conda.yaml"
   # benchmark:
    #    "benchmarks/benchmark_BQSR_BaseRecal_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR1_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR1_thrs=BQSR1_thrs,n_cpu=n_cpu)
    threads: BQSR1_thrs
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable} 2> {log}"

rule PrintReads_genome:
    """
    This tool writes out sequence read data.
    """
    input:
        r_bam = genome_int+"{sample}"+"_realigned.bam",
        r_idx = genome_int+"{sample}"+"_realigned.bai",
        outtable = genome_logs + "{sample}"+"_recal_data.table",
    output:
        recal_bam = alignpath_genomebqsr + "{sample}"+"_recal.bam",
        recal_bai = alignpath_genomebqsr + "{sample}"+"_recal.bai",
    params:
        gatk = gatk,
        ref=hg,
    log:
        genome_logs + "{sample}" + '_recalibrating_02.log'
    conda:
        "envs/wes_config_conda.yaml"
   # benchmark:
    #    "benchmarks/benchmark_BQSR_PrintReads_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR4_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR4_thrs=BQSR4_thrs, n_cpu=n_cpu)
    threads: BQSR4_thrs
    shell:
        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable} -nct {threads} -o {output.recal_bam} 2> {log}"






##########################
####  EXOME ALIGNMENT ####
##########################

rule map_to_exome:
    """
    This tool maps the samples to the reference exome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        fai = lambda wildcards: get_ref_fai(wildcards.sample),
        fasta_dict = lambda wildcards: get_ref_dict(wildcards.sample),
        indexes = lambda wildcards: get_ref_indexes(wildcards.sample),
        reference = lambda wildcards: get_ref(wildcards.sample),
        fastq1_trimmed = trimmedpath+"{sample}" + "_R1.fastq",
        fastq2_trimmed = trimmedpath+"{sample}" + "_R2.fastq",
    output:
        outfile =  exome_int + "{sample}.sam",
    params:
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = "{sample}",
    log:
        exome_logs + '{sample}_alignment.log',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_mapping_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    threads: map_thrs
    #resources: mem=6
    #version: 0.1
    message: ">> {wildcards.sample} - - Aligning to reference EXOME"
    shell:
        "bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {input.reference} {input.fastq1_trimmed} {input.fastq2_trimmed} > {output.outfile} 2> {log}"

rule sorting_exome:
    """
    This tool sorts the input SAM or BAM file by coordinate.
    """
    input:
        sam = exome_int + "{sample}.sam",
        tmp = tmp_dir_exome,
    output:
        outdir = exome_int + "{sample}_sorted.bam",
    log:
        exome_logs + '{sample}_sorting.log'
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard SortSam INPUT={input.sam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={input.tmp} 2> {log}"


rule extractUnmapped_extract:
    input:
        exome_int + "{sample}_sorted.bam",
    output:
        exome_int + "{sample}_unmapped.bam",
    log:
        exome_logs + '{sample}_unmapped.log'
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "samtools view -b -f 4 {input} > {output} 2> {log}"

rule extractUnmapped_convert:
    input:
        exome_int + "{sample}_unmapped.bam",
    output:
        unmapped1 = alignpath_exomeunmapped + "{sample}_R1.unmapped.fastq",
        unmapped2 = alignpath_exomeunmapped + "{sample}_R2.unmapped.fastq",
        unpair = alignpath_exomeunmapped + "{sample}_unpaired.unmapped.fastq",
    log:
        exome_logs + '{sample}_unmapped_fastq.log'
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "picard SamToFastq I={input} F={output.unmapped1} F2={output.unmapped2} FU={output.unpair} 2> {log}"



##########################
####   MT ALIGNMENT   ####
##########################

rule map_to_MT:
    """
    This tool maps the samples to the mitochondria genome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        MT_indexes,
        reference = MT,
        unmapped1 = alignpath_exomeunmapped + "{sample}_R1.unmapped.fastq",
        unmapped2 = alignpath_exomeunmapped + "{sample}_R2.unmapped.fastq",
        unpair = alignpath_exomeunmapped + "{sample}_unpaired.unmapped.fastq",
    output:
        outfile =  MT_int + "{sample}.sam",
    params:
        library = lambda wildcards: get_kit(wildcards.sample),
        platform = platform,
        name = "{sample}",
    log:
        MT_logs + '{sample}_alignment.log',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_mapping_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    threads: map_thrs
    #resources: mem=6
    #version: 0.1
    message: ">> {wildcards.sample} - Aligning to reference MITOCHONDRIAL GENOME"
    shell:
        "bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {input.reference} {input.unmapped1} {input.unmapped2} > {output.outfile} 2> {log}"

rule sorting_MT:
    """
    This tool sorts the input SAM or BAM file by coordinate.
    """
    input:
        sam = MT_int + "{sample}.sam",
        tmp = tmp_dir_MT,
    output:
        outdir = MT_int + "{sample}_sorted.bam",
    log:
        alignpath_MT + '{sample}_sorting.log'
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard SortSam INPUT={input.sam} OUTPUT={output.outdir} SORT_ORDER=coordinate TMP_DIR={input.tmp} 2> {log}"

rule marking_MT:
    """
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
    """
    input:
        sorted_bam = MT_int + "{sample}_sorted.bam",
        tmp = tmp_dir_MT_m,
    output:
        out = MT_int+"{sample}"+"_marked.bam",
    log:
        mx = alignpath_MT + '{sample}_metrix.log',
        mark = alignpath_MT + '{sample}_marking.log',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_mark_duplicates_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard MarkDuplicates"
        " INPUT={input.sorted_bam} OUTPUT={output.out} METRICS_FILE={log.mx} "
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR={input.tmp} 2> {log.mark}"

rule indexing_MT:
    """
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    """
    input:
        marked_bam = MT_int+"{sample}"+"_marked.bam",
    output:
        marked_bai = MT_int+"{sample}"+"_marked.bai",
    log:
        alignpath_MT + '{sample}_indexing.log',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_build_bam_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard BuildBamIndex INPUT={input.marked_bam} OUTPUT={output} 2> {log}"


rule RTC_MT:
    """
    This tool defines intervals to target for local realignment.
    """
    input:
        MT.replace('fasta', 'dict'),
        marked_bai = MT_int+"{sample}"+"_marked.bai",
        ref = MT+'.fai',
        indels_ref=indels_ref,
        gatk = gatk,
        marked_bam = MT_int+"{sample}"+"_marked.bam",
        bed = MT_bed + ".bed",
    output:
        MT_logs+"{sample}"+".intervals",
    params:
        ref=MT,
    log:
        MT_logs + "{sample}" + "_RTC.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{RT_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, RT_thrs=RT_thrs, n_cpu=n_cpu)
    threads: RT_thrs
    shell:
        "java -jar {input.gatk} -T RealignerTargetCreator -R {params.ref} -I {input.marked_bam} -L {input.bed} -ip 50 -known {input.indels_ref} -nt {threads} -o {output} 2> {log}"


rule IndelRealigner_MT:
    """
    This tool performs local realignment of reads around indels.
    """
    input:
        intvs = MT_logs+"{sample}"+".intervals",
        bam = MT_int+"{sample}"+"_marked.bam",
        idx = MT_int+"{sample}"+"_marked.bai",
    output:
        r_bam = MT_int+"{sample}"+"_realigned.bam",
        r_idx = MT_int+"{sample}"+"_realigned.bai",
    params:
        gatk = gatk,
        ref = MT,
        indels_ref=indels_ref,
    log:
        MT_logs + "{sample}" + "_IndelRealigner.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.intvs} -known {params.indels_ref} -ip 50 -o {output.r_bam} 2> {log}"


rule BaseRecal_MT:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    """
    input:
        r_bam = MT_int+"{sample}"+"_realigned.bam",
        r_idx = MT_int+"{sample}"+"_realigned.bai",
        dbsnp = dbsnp,
        bed = MT_bed + ".bed",
    output:
        outtable = MT_logs + "{sample}"+"_recal_data.table",
    params:
        gatk = gatk,
        ref=MT,
        indels_ref=indels_ref
    log:
        MT_logs + "{sample}" + '_recalibrating_01.log'
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_BQSR_BaseRecal_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR1_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR1_thrs=BQSR1_thrs,n_cpu=n_cpu)
    threads: BQSR1_thrs
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -L {input.bed} -ip 50 -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable} 2> {log}"

rule PrintReads_MT:
    """
    This tool writes out sequence read data.
    """
    input:
        r_bam = MT_int+"{sample}"+"_realigned.bam",
        r_idx = MT_int+"{sample}"+"_realigned.bai",
        outtable = MT_logs + "{sample}"+"_recal_data.table",
    output:
        recal_bam = alignpath_MTbqsr + "{sample}"+"_recal.bam",
        recal_bai = alignpath_MTbqsr + "{sample}"+"_recal.bai",
    params:
        gatk = gatk,
        ref=MT,
    log:
        MT_logs + "{sample}" + '_recalibrating_02.log'
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_BQSR_PrintReads_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR4_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR4_thrs=BQSR4_thrs, n_cpu=n_cpu)
    threads: BQSR4_thrs
    shell:
        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable} -nct {threads} -o {output.recal_bam} 2> {log}"



#######################
#     Variant Call    #
#######################


##############
###  MUTECT ##
##############


rule muTect_genome:
    """
    This step uses muTect 1.1.4 to call variants.
    """
    input:
        muTect = muTect,
        cosmic = cosmic,
        target = lambda wildcards: get_bed_patient(wildcards.patient),
        normal_bam = lambda wildcards: get_bam(wildcards.patient,'N', alignpath_genomebqsr),
        tumour_bam = lambda wildcards: get_bam(wildcards.patient,'T', alignpath_genomebqsr),
        normal_bai = lambda wildcards: get_bai(wildcards.patient,'N', alignpath_genomebqsr),
        tumour_bai = lambda wildcards: get_bai(wildcards.patient,'T', alignpath_genomebqsr),
    output:
        out = mutectpath_genome + "{patient}.tsv",
        cov = mutectpath_genome + "{patient}_cov_wig.log"
    params:
        dbsnp = dbsnp,
        ref = hg,
    log:
        m = mutectpath_genome_logs + "{patient}" + "_mutect.log",
        err = mutectpath_genome_logs + "{patient}" + "_mutect_err.log"
    conda:
        "envs/wes_config_conda_muTect.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    message: 'MuTect on genomic alignment'
    shell:
        "java -Xmx2g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --dbsnp {params.dbsnp} --intervals {input.target} --input_file:normal {input.normal_bam} --input_file:tumor {input.tumour_bam} --out {output.out} --coverage_file {output.cov} > {log.m} 2> {log.err}"

rule muTect_MT:
    """
    This step uses muTect 1.1.4 to call variants.
    """
    input:
        muTect = muTect,
        cosmic = cosmic,
        target = MT_bed + ".bed",
        normal_bam = lambda wildcards: get_bam(wildcards.patient,'N', alignpath_MTbqsr),
        tumour_bam = lambda wildcards: get_bam(wildcards.patient,'T', alignpath_MTbqsr),
        normal_bai = lambda wildcards: get_bai(wildcards.patient,'N', alignpath_MTbqsr),
        tumour_bai = lambda wildcards: get_bai(wildcards.patient,'T', alignpath_MTbqsr),
    output:
        out = mutectpath_MT + "{patient}.tsv",
        cov = mutectpath_MT + "{patient}_cov_wig.log"
    params:
        dbsnp = dbsnp,
        ref = MT,
    log:
        m = mutectpath_MT_logs + "{patient}" + "_mutect.log",
        err = mutectpath_MT_logs + "{patient}" + "_mutect_err.log"
    conda:
        "envs/wes_config_conda_muTect.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    message: 'MuTect on MT alignment'
    shell:
        "java -Xmx2g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --dbsnp {params.dbsnp} --intervals {input.target} --input_file:normal {input.normal_bam} --input_file:tumor {input.tumour_bam} --out {output.out} --coverage_file {output.cov} > {log.m} 2> {log.err}"






###############
###  VARSCAN ##
###############

rule varscan_genome_sample:
    """
    This step uses Varscan to call variants.
    """
    input:
        target = lambda wildcards: get_bed(wildcards.sample),
        recal_bam = alignpath_genomebqsr + "{sample}"+"_recal.bam",
        recal_bai = alignpath_genomebqsr + "{sample}"+"_recal.bai",
    output:
        out = mpileup_varscan_genome + "{sample}" + ".mpileup",
    params:
        ref = hg,
    log:
        varscan_genome_logs + "{sample}" + "_mpileup.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    shell:
        "samtools mpileup -B -q 1 -f {params.ref} --positions {input.target} {input.recal_bam}"

rule varscan_genome_patient:
    """
    This step uses Varscan to call variants.
    """
    input:
        normal = lambda wildcards: get_mpileup(wildcards.patient,'N', mpileup_varscan_genome),
        tumour = lambda wildcards: get_mpileup(wildcards.patient,'T', mpileup_varscan_genome),
    output:
        snv = varscanpath_genome + "{patient}_snv.tsv",
        indel = varscanpath_genome + "{patient}_indel.tsv",
    params:
        ref = hg,
    log:
        varscan_genome_logs + "{patient}" + "_varscan.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    message: 'VarScan on genomic alignment'
    shell:
        "varscan somatic {input.normal} {input.tumour} --output-snp {output.snv} --output-indel {output.indel} --min-avg-qual 15 --strand_filter 1 --min-var-freq 0.05 --somatic-p-value 0.05 2> {log}"


rule varscan_MT_sample:
    """
    This step uses Varscan to call variants.
    """
    input:
        target = MT_bed + ".bed",
        recal_bam = alignpath_MTbqsr + "{sample}"+"_recal.bam",
        recal_bai = alignpath_MTbqsr + "{sample}"+"_recal.bai",
    output:
        out = mpileup_varscan_MT + "{sample}" + ".mpileup",
    params:
        ref = MT,
    log:
        varscan_MT_logs + "{sample}" + "_mpileup.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    shell:
        "samtools mpileup -B -q 1 -f {params.ref} --positions {input.target} {input.recal_bam}"

rule varscan_MT_patient:
    """
    This step uses Varscan to call variants.
    """
    input:
        normal = lambda wildcards: get_mpileup(wildcards.patient,'N', mpileup_varscan_MT),
        tumour = lambda wildcards: get_mpileup(wildcards.patient,'T', mpileup_varscan_MT),
    output:
        snv = varscanpath_MT + "{patient}_snv.tsv",
        indel = varscanpath_MT + "{patient}_indel.tsv",
    params:
        ref = MT,
    log:
        varscan_MT_logs + "{patient}" + "_varscan.log"
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_muTect_ref_{patient}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    message: 'VarScan on MT alignment'
    shell:
        "varscan somatic {input.normal} {input.tumour} --output-snp {output.snv} --output-indel {output.indel} --min-avg-qual 15 --strand_filter 1 --min-var-freq 0.05 --somatic-p-value 0.05 2> {log}"


######################
###     FILTERING  ###
######################

########## FILTER SOMATIC

rule MergeWriteVariants:
    input:
        mutect_genome = mutectpath_genome + "{patient}.tsv",
        mutect_MT = mutectpath_MT + "{patient}.tsv",
        snv_genome = varscanpath_genome + "{patient}_snv.tsv",
        indel_genome = varscanpath_genome + "{patient}_indel.tsv",
        snv_MT = varscanpath_MT + "{patient}_snv.tsv",
        indel_MT = varscanpath_MT + "{patient}_indel.tsv",
    output:
        m = variants_filter + "{patient}_mutect.tsv",
        snp = variants_filter + "{patient}_varscan_snv_temp1.tsv",
        indel = variants_filter + "{patient}_varscan_indel_temp1.tsv",
        v_all = variants_filter + "{patient}_varscan.tsv",
    params:
        scripts=scripts,
        name = "{patient}",
        outdir = variants_filter,
    # benchmark:
    script:
        "{params.scripts}"+"MergeWriteVariants.py"

rule filterSomatic_SNV:
    input:
        snp = variants_filter + "{patient}_varscan_snv_temp1.tsv",
        indel = variants_filter + "{patient}_varscan_indel_temp1.tsv",
    output:
        variants_filter + "{patient}_varscan_snv_temp2.tsv",
    log:
        variants_filter_logs + "{patient}_somaticFilter_snv_err.log"
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "varscan somaticFilter {input.snp} --indel-file {input.indel} --min-coverage 1 --min-reads2 2 --min-var-freq 0.1 --output-file {output} 2> {log}"

rule filterSomatic_Indel:
    input:
        variants_filter + "{patient}_varscan_indel_temp1.tsv",
    output:
        variants_filter + "{patient}_varscan_indel_temp2.tsv",
    log:
        variants_filter_logs + "{patient}_somaticFilter_indel_err.log"
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "varscan somaticFilter {input} --min-coverage 1 --min-reads2 2 --min-var-freq 0.1 --output-file {output} 2> {log}"


rule KeepOnlyVariants:
    input:
        m = variants_filter + "{patient}_mutect.tsv",
        snv = variants_filter + "{patient}_varscan_snv_temp2.tsv",
        indel = variants_filter + "{patient}_varscan_indel_temp2.tsv",
    output:
        m_tsv = variants_filter +"{patient}_mutect_somatic.tsv",
        vsn_tsv = variants_filter +"{patient}_varscan_somatic.tsv",
    params:
        scripts=scripts,
        name = "{patient}",
        outdir = variants_filter,
    # benchmark:
    script:
        "{params.scripts}"+"KeepOnlyVariants.py"


############# ANNOVAR

############## MUTECT

rule MutectSplitVariants:
    input:
        m_tsv = variants_filter +"{patient}_mutect_somatic.tsv",
    output:
        ann_genome = variants_filter +"{patient}_mutect_annovar_genome.tsv",
        ann_MT = variants_filter +"{patient}_mutect_annovar_mt.tsv",
    params:
        scripts = scripts,
        workdir = variants_filter,
        name = "{patient}"
    script:
        "{params.scripts}"+"MutectSplitVariants.py"


rule annovar_genome_mutect:
    input:
        ann_genome = variants_filter +"{patient}_mutect_annovar_genome.tsv",
        annovar_dbs = annovar_dbs,
    output:
        variants_filter_logs + "{patient}_mutect_genome"
    params:
        protocol = config['protocols']['genome'],
        operation = config['operations']['genome'],
        tableannovar = tableannovar,
        humandb = humandb,
        build_ver = build_ver,
    log:
        variants_filter_logs + "{patient}_mutect_genome_annovar_err.log"
    shell:
        "perl {params.tableannovar} {input.ann_genome} {params.humandb} -buildver {params.build_ver} -remove -protocol {params.protocol} -operation {params.operation} 2> {log}"

rule annovar_MT_mutect:
    input:
        ann_MT = variants_filter +"{patient}_mutect_annovar_mt.tsv",
        annovar_dbs = annovar_dbs,
    output:
        variants_filter_logs + "{patient}_mutect_mt_g"
    params:
        protocol = config['protocols']['MT'],
        operation = config['operations']['MT'],
        tableannovar = tableannovar,
        humandb = humandb,
        mitochondrial_ver = mitochondrial_ver,
    log:
        variants_filter_logs + "{patient}_mutect_mt_g_annovar_err.log"
    shell:
        "perl {params.tableannovar} {input.ann_genome} {params.humandb} -buildver {params.mitochondrial_ver} -remove -protocol {params.protocol} -operation {params.operation} 2> {log}"

######## VARSCAN

rule VarscanSplitVariants:
    input:
        vsn_tsv = variants_filter +"{patient}_varscan_somatic.tsv",
    output:
        ann_genome = variants_filter +"{patient}_varscan_annovar_genome.tsv",
        ann_MT = variants_filter +"{patient}_varscan_annovar_mt.tsv",
    params:
        scripts = scripts,
        workdir = variants_filter,
        name = "{patient}"
    script:
        "{params.scripts}"+"VarscanSplitVariants.py"

rule annovar_genome_varscan:
    input:
        ann_genome = variants_filter +"{patient}_varscan_annovar_genome.tsv",
        annovar_dbs = annovar_dbs,
    output:
        variants_filter_logs + "{patient}_varscan_genome"
    params:
        protocol = config['protocols']['genome'],
        operation = config['operations']['genome'],
        tableannovar = tableannovar,
        humandb = humandb,
        build_ver = build_ver,
    log:
        variants_filter_logs + "{patient}_varscan_genome_annovar_err.log"
    shell:
        "perl {params.tableannovar} {input.ann_genome} {params.humandb} -buildver {params.build_ver} -remove -protocol {params.protocol} -operation {params.operation} 2> {log}"

rule annovar_MT_varscan:
    input:
        ann_MT = variants_filter +"{patient}_varscan_annovar_mt.tsv",
        annovar_dbs = annovar_dbs,
    output:
        variants_filter_logs + "{patient}_varscan_mt_g"
    params:
        protocol = config['protocols']['MT'],
        operation = config['operations']['MT'],
        tableannovar = tableannovar,
        humandb = humandb,
        mitochondrial_ver = mitochondrial_ver,
    log:
        variants_filter_logs + "{patient}_varscan_mt_g_annovar_err.log"
    shell:
        "perl {params.tableannovar} {input.ann_genome} {params.humandb} -buildver {params.mitochondrial_ver} -remove -protocol {params.protocol} -operation {params.operation} 2> {log}"


rule MergeMutectVarscan:
    input:
        variants_filter_logs + "{patient}_mutect_genome",
        variants_filter_logs + "{patient}_mutect_mt_g",
        variants_filter_logs + "{patient}_varscan_genome",
        variants_filter_logs + "{patient}_varscan_mt_g",
    output:
        variants_filter + "{patient}_all_somatic_annotated.tsv",
    params:
        scripts = scripts,
        workdir = variants_filter_logs,
        name = "{patient}",
    script:
        "{params.scripts}"+"MergeMtVar.py"

############ filter on exonic, non-synonymous, polymorphisms


rule create_rsID_table:
    input:
        not_filtered = variants_filter + "{patient}_all_somatic_annotated.tsv",
        maf = "/",
    output:
        variants_filter + "{patient}_rsID.tsv",
    params:
        scripts = scripts,
    script:
        "{params.scripts}" + "fEP_script.py"

rule filterExonicPolymorphic:
    input:
        variants_filter + "{patient}_rsID.tsv",
    output:
        "{patient}_rsID_maf.tsv"
    params:
        rsIDfilter = rsIDfilter,
        workdir = variants_filter,
    log:
        info = variants_filter_logs + "{patient}_rsID_session_info.log",
        err = variants_filter_logs + "{patient}_rsID_err.log",
    conda:
        "envs/wes_config_conda.yaml",
    shell:
        "Rscript {params.rsIDfilter} {params.workdir} > {log.info} 2> {log.err} && "
        "cat {log.info} {log.err} > {log.err}"

rule checkMAFs:
    input:
        not_filtered = variants_filter + "{patient}_all_somatic_annotated.tsv",
        maf = "{patient}_rsID_maf.tsv",
    output:
        variants_filter + "{patient}_all_somatic_annotated_filtered.tsv",
    params:
        scripts = scripts,
    script:
        "{params.scripts}"+"fEP_script.py"

###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

#################################################
#                   GET FASTA                   #
#################################################

rule download_reference:
    """download the hg19 human reference genome from 1000genome"""
    output:
        zipped = temp(hg+'.gz'),
    version: 0.1
    benchmark:
        "benchmarks/benchmark_downloadreference_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && "
        "mv human_g1k_v37.fasta.gz {output.zipped} "

rule gunzip_reference:
    input:
        zipped = hg+'.gz'
    output:
        hg
    benchmark:
        "benchmarks/benchmark_gunzip_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip {input.zipped} || true"

rule get_nextera_fasta:
    input:
        hg = hg,
        bed_fixed = nextera_bed + "_fixed.bed"
    output:
        nextera,
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}"

rule get_nextexp_fasta:
    input:
        hg = hg,
        bed_fixed = nextexp_bed + "_fixed.bed"
    output:
        nextera_expanded,
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}"

rule get_truseq_fasta:
    input:
        hg = hg,
        bed_fixed = truseq_bed + "_fixed.bed"
    output:
        truseq,
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "bedtools getfasta -fi {input.hg} -bed {input.bed_fixed} -fo {output}"

rule get_MT_fasta:
    input:
        hg,
    output:
        MT,
    conda:
        "envs/wes_config_conda.yaml"
    shell:
        "samtools faidx {input} MT > {output}"


########################
#    GET SEVERAL FILES #
########################
rule download_indels_ref:
    """download the indel reference from 1000genome """
    output:
        indel_zipped= temp(indels_ref+'.gz'),
    benchmark:
        "benchmarks/benchmark_downloadindels_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz && "
        "mv Mills_and_1000G_gold_standard.indels.b37.vcf.gz {output.indel_zipped}"

rule gunzip_indelref:
    input:
        indel_zipped = indels_ref+'.gz'
    output:
        indels_ref
    benchmark:
        "benchmarks/benchmark_gunzip_indelref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip {input.indel_zipped} || true"


rule download_dbsnp:
    """download dbsnp_138.b37 from 1000genome """
    output:
        dbsnp_zipped= temp(dbsnp+'.gz'),
    benchmark:
        "benchmarks/benchmark_downloaddbsnp_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz && "
        "mv dbsnp_138.b37.vcf.gz {output.dbsnp_zipped}"

rule gunzip_dbsnp:
    input:
        dbsnp_zipped= dbsnp+'.gz'
    output:
        dbsnp
    benchmark:
        "benchmarks/benchmark_gunzip_dbsnp_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip {input.dbsnp_zipped} || true"


# https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2

rule download_cosmic:
    """download  cosmic from broadinstitute"""
    output:
        cosmic=cosmic,
    shell:
        "wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf && "
        "mv b37_cosmic_v54_120711.vcf {output.cosmic}"

rule download_nextexp_bed:
    """download target from illumina"""
    output:
        nextexp_bed = temp(nextexp_bed + ".bed"),
    shell:
        "wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed && "
        "mv nexterarapidcapture_expandedexome_targetedregions.bed {output.nextexp_bed}"

rule fix_nextexp_bed:
    """fix target: remove prefix chr in first column"""
    input:
        nextexp_bed + ".bed",
    output:
        nextexp_bed + "_fixed.bed",
    params:
        scripts = scripts,
        line_to_skip = None,
    script:
        "{params.scripts}" + "editBEDChromField.py"

rule download_nextera_bed:
    """download target from illumina"""
    output:
        nextera_bed = temp(nextera_bed + ".bed"),
    shell:
        "wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions.bed && "
        "mv nexterarapidcapture_exome_targetedregions.bed {output}"

rule fix_nextera_bed:
    """fix target: remove prefix chr in first column and last column(name)"""
    input:
        nextera_bed + ".bed",
    output:
        nextera_bed + "_fixed.bed",
    params:
        scripts = scripts,
        line_to_skip = None,
    script:
        "{params.scripts}" + "editBEDChromField.py"

rule download_truseq_bed:
    """download target from illumina"""
    output:
        truseq_bed = temp(truseq_bed + ".bed"),
    shell:
        "wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/truseq/truseq-exome-targeted-regions-manifest-v1-2-bed.zip && "
        "gunzip truseq-exome-targeted-regions-manifest-v1-2-bed.zip && "
        "mv truseq-exome-targeted-regions-manifest-v1-2.bed {output}"

rule fix_truseq_bed:
    """fix target: remove prefix chr in first column and last column(name)"""
    input:
        truseq_bed + ".bed",
    output:
        truseq_bed + "_fixed.bed",
    params:
        scripts = scripts,
        line_to_skip = None,
    script:
        "{params.scripts}" + "editBEDChromField.py"

rule create_MTbed:
    output:
        MT_bed = MT_bed + ".bed",
    shell:
        "echo 'MT 1 16569' > {output}"


rule download_annovar_databases:
    input:
        annovar = annovar,
    output:
        annovar_dbs,
    params:
        humandb = humandb,
    shell:
        "{input.annovar} -downdb -buildver hg19 -webfrom annovar snp138 {params.humandb} && "
        "{input.annovar} -downdb 1000g2012apr {params.humandb} -buildver hg19"


#########################################################
#                   REFERENCES INDEXING                  #
#########################################################

###############
##     HG    ##
###############
rule hg_index_bwa:
    """
    Generate the index of the reference genome for the bwa program.
    """
    input:
        hg,
    output:
        hg_indexes,
    conda:
        "envs/wes_config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_hg_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {hg}"

rule index_picard_hg:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        hg=hg,
    output:
        hg.replace('fasta', 'dict'),
    conda:
        "envs/wes_config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input.hg} O={output}"

rule index_samtools_hg:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        hg=hg,
    output:
        hg+'.fai',
    conda:
        "envs/wes_config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input.hg} "

###############
##  Nextera  ##
###############

rule nextera_index_bwa:
    input:
        nextera,
    output:
        nextera_indexes,
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
    #    "benchmarks/benchmark_hg_hg_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {input}"

rule index_picard_nextera:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        nextera,
    output:
        nextera.replace('fasta', 'dict'),
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

rule index_samtools_nextera:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        nextera,
    output:
        nextera+'.fai',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input} "


########################
##  Nextera_expanded  ##
########################

rule nextexp_index_bwa:
    input:
        nextera_expanded,
    output:
        nextexp_indexes,
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
    #    "benchmarks/benchmark_hg_hg_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {input}"

rule index_picard_nextexp:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        nextera_expanded,
    output:
        nextera_expanded.replace('fasta', 'dict'),
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

rule index_samtools_nextexp:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        nextera_expanded,
    output:
        nextera_expanded+'.fai',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input} "

##############
##  Truseq  ##
##############
rule truseq_index_bwa:
    input:
        truseq,
    output:
        truseq_indexes,
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
    #    "benchmarks/benchmark_hg_hg_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {input}"

rule index_picard_truseq:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        truseq,
    output:
        truseq.replace('fasta', 'dict'),
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

rule index_samtools_truseq:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        truseq,
    output:
        truseq+'.fai',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input} "

#############
#     MT    #
#############
rule MT_index_bwa:
    input:
        MT,
    output:
        MT_indexes,
    conda:
        "envs/wes_config_conda.yaml"
    #benchmark:
    #    "benchmarks/benchmark_hg_hg_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {input}"

rule index_picard_MT:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        MT,
    output:
        MT.replace('fasta', 'dict'),
    conda:
        "envs/wes_config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

rule index_samtools_MT:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        MT,
    output:
        MT+'.fai',
    conda:
        "envs/wes_config_conda.yaml"
    # benchmark:
    #     "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input} "


#########################################################
#                   CHECKING REQUIREMENTS               #
#########################################################

rule check_GATK:
    """
        This step check if GATK is present in the directory set in config file.
    """
    output:
        gatk,
    priority: 5
    shell:
        "echo 'Error. Genome Analysis ToolKit not found in softwares directory.' && "
        "exit 1"

rule check_muTect:
    """
        This step check if muTect is present in the directory set in config file.
    """
    output:
        muTect,
    priority: 4
    shell:
        "echo 'Error. muTect not found in softwares directory.' && "
        "exit 1"

rule check_Annovar:
    """
        This step check if annovar is present in the directory set in config file.
    """
    output:
        annovar,
    priority: 3
    shell:
        "echo 'Error. Annovar not found in softwares directory.' && "
        "exit 1"

rule all_fastq:
    output:
        all_fastq,
    priority: 2
    shell:
        "echo 'Error. all_fastq.log not found.' && "
        "exit 1"

#########################################################
#                        THE END                        #
#########################################################
