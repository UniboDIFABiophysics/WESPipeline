#snakemake --dag | dot -Tpdf > dag.pdf
#snakemake --detailed-summary > provenance.tsv
#snakemake --cores 6 --resources mem=12

shell.executable("/bin/bash") # The pipeline works only on /bin/bash

# Open and read configuration file for samples informations
import yaml
with open("samples_config.yaml", "r") as f:
    samples_cfg = yaml.load(f)

# Set configuration file

configfile: "config.yaml"

# Extract the main directory
import os
home = os.path.expanduser("~")

#####################################
#  Set parameters from config file  #
#####################################

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

# References
hg = home + config['ref-files']['hg'] # Human Genome Reference
bwa_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"] # Indexes for hg
indels_ref = home + config['ref-files']['indels_ref'] # Set of known indels
dbsnp = home + config['ref-files']['dbsnp'] # SNP database
cosmic = home + config['ref-files']['cosmic'] # Catalog of somatic mutation in cancer
humandb = home + config['ref-files']['humandb'] # Annovar human databases folder
build_ver = config['ref-files']['build_ver'] # Set build version
dbsnp_ver = config['ref-files']['dbsnp_ver'] # Set SNP database version
kg_ver = config['ref-files']['kg_ver'] # Set  version # Set ethnicity groups based on time
mitochondrial_ver = config['ref-files']['mitochondrial_ver'] # Set parameter command for annotating mitochondria variants

# Folders
scripts = config['folders']['scripts']
resultdir = config['folders']['resultdir']

# Softwares
gatk = home + config['softwares']['gatk']
muTect = home + config['softwares']['muTect']
annovar = home + config['softwares']['annovar']
convert2annovar = home + config['softwares']['convert2annovar']

# Sample details
platform = config['sample-details']['platform'] # Set platform for mapping
library = config['sample-details']['library'] # Set library for mapping
target = home + config['sample-details']['target'] # Set target intervals for exome analysis

# HaplotypeCaller parameters
genotyping_mode = config['hap-caller']['genotyping_mode'] # Set genotyping mode
stand_emit_conf = config['hap-caller']['stand_emit_conf'] # stand_emit_conf 10.0 means that it won’t report any potential SNPs with a quality below 10.0;
stand_call_conf = config['hap-caller']['stand_call_conf'] # but unless they meet the quality threshold set by -stand_call_conf (30.0, in this case), they will be listed as failing the quality filter

# Hard Filter parameters
filter_exp_snps = config['hard-filter']['snps'] # Set command to define thresholds for snps
filter_exp_indels = config['hard-filter']['indels'] # Set command to define thresholds for indels

# Annovar databases
annovar_dbs = [home + config['annovar_dbs']['hg19_db'] , home + config['annovar_dbs']['snp138_db']]

#LODn parameters
min_n_cov = config['LODn']['min_n_cov'] #Set the minimum read depth for normal and tumor samples
min_t_cov = config['LODn']['min_t_cov']

#######################
#   Get fastq files   #
#######################

# Define variables
patients = [p for p in samples_cfg['patients']]
PATHS = [ps for p in samples_cfg['patients'] for ps in samples_cfg['patients'][p]]
SAMPLES = [name.split('/')[-1] for name in PATHS]
sicks = [s for s in samples_cfg['patients'] if len(samples_cfg['patients'][s])==2]

# Wildcard costrains necessary for search only certain names
wildcard_constraints:
    path = "("+"|".join(PATHS)+")",
    sample = "("+"|".join(SAMPLES)+")",
    sick = "("+"|".join(sicks)+")",

# Container to get patient name from sample name
sample_to_patient = {}

for patient in patients:
    for s in samples_cfg['patients'][patient]:
        s = s.split('/')[-1]
        sample_to_patient[s] = patient

# Container to get normal/tumour sample from patient's name
sicks_list = {}
for sick in sicks:
    sicks_list[sick] = {}
    sicks_list[sick]['N'] = samples_cfg['patients'][sick][0].split('/')[-1]
    sicks_list[sick]['T'] = samples_cfg['patients'][sick][1].split('/')[-1]

# Get sample from path
def path_to_sample(wildcards, index = '1'):
    return (PATHS[SAMPLES.index(wildcards)] + "_{index}.fastq".format(index=index))

# Get normal/tumour bam file from patient's name
def get_bam(wildcards, sample_type='N'):
    return (resultdir+sicks_list[wildcards][sample_type]+"_recal.bam")

# Get normal/tumour bai file from patient's name
def get_bai(wildcards, sample_type='N'):
    return (resultdir+sicks_list[wildcards][sample_type]+"_recal.bai")

# Get normal sample name from patient's name
def get_code(wildcards):
    return (sicks_list[wildcards]['N'])

# Get table or vcf tumour sample from patient's name
def get_lodn_infile(wildcards,format,idx=""):
    if format == 'table':
        return (resultdir+sicks_list[wildcards]['T']+".tsv")
    elif format == 'vcf':
        return (resultdir+sicks_list[wildcards]['T']+"_filtered_variants.vcf"+idx)

#######################
# Workflow final rule #
#######################

rule all:
    """
    THE END.
    """
    input:
        expand(resultdir+"{sample}"+".tsv", sample=SAMPLES),
        expand(resultdir+ 'mutect_ann/' + "{sick}"+".tsv", sick=sicks),
        expand(resultdir+"{sick}"+"_lodn.vcf",sick=sicks),
        expand(resultdir+"{sick}"+"_lodn_table.tsv",sick=sicks),
    benchmark:
        "benchmarks/benchmark_rule_all_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        pass

#########################################################
#                      GATK-LODn                        #
#########################################################

#############################
#        gatk_pipe          #
#############################

rule Unzip_sample1:
    """
    Unzip first sample run file.
    """
    input:
        sample1_zipped = "{path}"+"_1.fastq.gz",
    output:
        sample1_unzipped = temp("{path}"+"_1.fastq"),
    shell:
        "gunzip -c {input.sample1_zipped} > {output.sample1_unzipped}"

rule Unzip_sample2:
    """
    Unzip second sample run file.
    """
    input:
        sample2_zipped = "{path}"+"_2.fastq.gz",
    output:
        sample2_unzipped = temp("{path}"+"_2.fastq"),
    shell:
        "gunzip -c {input.sample2_zipped} > {output.sample2_unzipped}"

rule mapping:
    """
    This tool maps the samples to the reference genome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        bwa_indexes,
        sample1 = lambda wildcards: path_to_sample(wildcards.sample,1),
        sample2 = lambda wildcards: path_to_sample(wildcards.sample,2),
    output:
        outfile =  temp(resultdir+"{sample}"+".sam"),
    params:
        reference = hg,
        library = library,
        platform = platform,
        name = "{sample}",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_mapping_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{map_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, map_thrs=map_thrs, n_cpu=n_cpu)
    threads: map_thrs
    resources: mem=6
    version: 0.1
    message: "bwa is aligning the sample '{params.name}' with the reference genome"
    shell:
        "bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {params.reference} {input.sample1} {input.sample2} > {output.outfile}"


rule sort_picard:
    """
    This tool sorts the input SAM or BAM file by coordinate.
    """
    input:
        r = resultdir+"{sample}"+".sam",
    output:
        temp(resultdir+"{sample}"+"_sorted.bam"),
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard SortSam INPUT={input.r} OUTPUT={output}.tmp SORT_ORDER=coordinate"
        " && [ -s {output}.tmp ] && mv {output}.tmp {output}"


rule mark_duplicates:
    """
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
    """
    input:
        r=resultdir+"{sample}"+"_sorted.bam",
    output:
        temp(resultdir+"{sample}"+"_dedup.bam"),
    params:
        metricsfile=resultdir+"{sample}"+"_dedup.metrics.txt",
        outdir=resultdir,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_mark_duplicates_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard MarkDuplicates"
        " INPUT={input.r} OUTPUT={output}.tmp METRICS_FILE={params.metricsfile}"
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000"
        " && [ -s {output}.tmp ] && mv {output}.tmp {output}"


rule build_bam_index:
    """
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    """
    input:
        r=resultdir+"{sample}"+"_dedup.bam",
    output:
        temp(resultdir+"{sample}"+"_dedup.bai"),
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_build_bam_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard BuildBamIndex INPUT={input.r} OUTPUT={output}"


rule realigner_target_creator:
    """
    This tool defines intervals to target for local realignment.
    """
    input:
        hg.replace('fasta', 'dict'),
        idx=resultdir+"{sample}"+"_dedup.bai",
        ref=hg+'.fai',
        indels_ref=indels_ref,
        gatk = gatk,
        seq=resultdir+"{sample}"+"_dedup.bam",
    output:
        temp(resultdir+"{sample}"+".intervals"),
    params:
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{RT_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, RT_thrs=RT_thrs, n_cpu=n_cpu)
    threads: RT_thrs
    shell:
        "java -jar {input.gatk} -T RealignerTargetCreator -R {params.ref} -I {input.seq} -known {input.indels_ref} -nt {threads} -o {output}"


rule IndelRealigner:
    """
    This tool performs local realignment of reads around indels.
    """
    input:
        target=resultdir+"{sample}"+".intervals",
        bam=resultdir+"{sample}"+"_dedup.bam",
        idx=resultdir+"{sample}"+"_dedup.bai",
    output:
        r_bam = temp(resultdir+"{sample}"+"_realigned.bam"),
        r_idx = temp(resultdir+"{sample}"+"_realigned.bai"),
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.target} -known {params.indels_ref} -o {output.r_bam}"

rule BQSR_step_1:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    """
    input:
        r_bam=resultdir+"{sample}"+"_realigned.bam",
        r_idx=resultdir+"{sample}"+"_realigned.bai",
        dbsnp = dbsnp,
    output:
        outtable1 = temp(resultdir+"{sample}"+"_recal_data.table"),
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR1_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR1_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR1_thrs=BQSR1_thrs,n_cpu=n_cpu)
    threads: BQSR1_thrs
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable1}"

rule BQSR_step_2:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a post recalibrated data table.
    """
    input:
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
        r_bam=resultdir+"{sample}"+"_realigned.bam",
        r_idx=resultdir+"{sample}"+"_realigned.bai",
    output:
        outtable2 = temp(resultdir+"{sample}"+"_post_recal_data.table"),
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref,
        dbsnp = dbsnp,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR2_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR2_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR2_thrs=BQSR2_thrs, n_cpu=n_cpu)
    threads: BQSR2_thrs
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {params.dbsnp} -knownSites {params.indels_ref} -BQSR {input.outtable1} -nct {threads} -o {output.outtable2}"

rule BQSR_step_3:
    """
    This tool creates plots to visualize base recalibration results.
    """
    input:
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
        outtable2 = resultdir+"{sample}"+"_post_recal_data.table",
    output:
        plots = resultdir+"{sample}"+"_recalibrationPlots.pdf",
    params:
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR3_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {input.outtable1} -after {input.outtable2} -plots {output.plots}"

rule BQSR_step_4:
    """
    This tool writes out sequence read data.
    """
    input:
        r_bam = resultdir+"{sample}"+"_realigned.bam",
        r_idx=resultdir+"{sample}"+"_realigned.bai",
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
        plots = resultdir+"{sample}"+"_recalibrationPlots.pdf",
    output:
        recal_bam = resultdir+"{sample}"+"_recal.bam",
        recal_bai = resultdir+"{sample}"+"_recal.bai",
    params:
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR4_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{BQSR4_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, BQSR4_thrs=BQSR4_thrs, n_cpu=n_cpu)
    threads: BQSR4_thrs
    shell:
        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable1} -nct {threads} -o {output.recal_bam}"

rule HaplotypeCaller:
    """
    GATK Hard Filtering Variants
    This tool calls germline SNPs and indels via local re-assembly of haplotypes.
    """
    input:
        recal_bam = resultdir+"{sample}"+"_recal.bam",
        recal_bai = resultdir+"{sample}"+"_recal.bai",
    output:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf",
        vcf_idx = resultdir+"{sample}"+"_raw.vcf.idx"
    params:
        gatk = gatk,
        ref=hg,
        genotyping_mode = genotyping_mode,
        stand_emit_conf = stand_emit_conf,
        stand_call_conf = stand_call_conf,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HaplotypeCaller_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HC_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HC_thrs=HC_thrs, n_cpu=n_cpu)
    threads: HC_thrs
    shell:
        "java -jar {params.gatk} -T HaplotypeCaller -R {params.ref} -I {input.recal_bam} --genotyping_mode {params.genotyping_mode} -stand_emit_conf {params.stand_emit_conf} -stand_call_conf {params.stand_call_conf} -nct {threads} -o {output.vcf_raw}"

# Since GATK 3.7 -stand_emit_conf {params.stand_emit_conf} has been deprecated

rule HardFilter_1SNPs:
    """
    GATK Hard Filtering Variants
    This tool extracts the snps from the call set.
    """
    input:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf",
        vcf_idx = resultdir+"{sample}"+"_raw.vcf.idx"
    output:
        raw_snps = temp(resultdir+"{sample}"+"_snps.vcf"),
        raw_snps_idx = temp(resultdir+"{sample}"+"_snps.vcf.idx"),
    params:
        gatk = gatk,
        ref=hg,
        selectType = '-selectType SNP'
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter1SNPs_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HF1_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HF1_thrs=HF1_thrs, n_cpu=n_cpu)
    threads: HF1_thrs
    shell:
        "java -jar {params.gatk} -T SelectVariants -R {params.ref} -V {input.vcf_raw} {params.selectType} -o {output.raw_snps}"

rule HardFilter_2SNPs:
    """
    GATK Hard Filtering Variants
    This tool applies the filter to the SNP call set.
    """
    input:
        raw_snps = resultdir+"{sample}"+"_snps.vcf",
        raw_snps_idx = temp(resultdir+"{sample}"+"_snps.vcf.idx"),
    output:
        filt_snps = temp(resultdir+"{sample}"+"_hard_filtered.snps.vcf"),
        filt_snps_idx = temp(resultdir+"{sample}"+"_hard_filtered.snps.vcf.idx")
    params:
        gatk = gatk,
        ref=hg,
        filter_exp_snps = filter_exp_snps,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter2SNPs_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HF2_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HF2_thrs=HF2_thrs, n_cpu=n_cpu)
    threads: HF2_thrs
    shell:
        "java -jar {params.gatk} -T VariantFiltration -R {params.ref} -V {input.raw_snps} {params.filter_exp_snps} -o {output.filt_snps}"

rule HardFilter_1Indels:
    """
    GATK Hard Filtering Variants
    This tool extracts the Indels from the call set.
    """
    input:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf",
        vcf_idx = resultdir+"{sample}"+"_raw.vcf.idx",
    output:
        raw_indels = temp(resultdir+"{sample}"+"_indels.vcf"),
        raw_indels_idx = temp(resultdir+"{sample}"+"_indels.vcf.idx"),
    params:
        gatk = gatk,
        ref=hg,
        selectType = '-selectType INDEL',
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter1Indels_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HF1_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HF1_thrs=HF1_thrs, n_cpu=n_cpu)
    threads: HF1_thrs
    shell:
        "java -jar {params.gatk} -T SelectVariants -R {params.ref} -V {input.vcf_raw} {params.selectType} -o {output.raw_indels}"

rule HardFilter_2Indels:
    """
    GATK Hard Filtering Variants
    This tool applies the filter to the Indel call set.
    """
    input:
        raw_indels = resultdir+"{sample}"+"_indels.vcf",
        raw_indels_idx = temp(resultdir+"{sample}"+"_indels.vcf.idx"),
    output:
        filt_indels = temp(resultdir+"{sample}"+"_hard_filtered.indels.vcf"),
        filt_indels_idx = temp(resultdir+"{sample}"+"_hard_filtered.indels.vcf.idx"),
    params:
        gatk = gatk,
        ref=hg,
        filter_exp_indels = filter_exp_indels,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter2Indels_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HF2_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HF2_thrs=HF2_thrs, n_cpu=n_cpu)
    threads: HF2_thrs
    shell:
        "java -jar {params.gatk} -T VariantFiltration -R {params.ref} -V {input.raw_indels} {params.filter_exp_indels} -o {output.filt_indels}"

rule HardFilter_Combine:
    """
    GATK Hard Filtering Variants
    This tool combines variants.
    """
    input:
        filt_snps = resultdir+"{sample}"+"_hard_filtered.snps.vcf",
        filt_snps_idx = temp(resultdir+"{sample}"+"_hard_filtered.snps.vcf.idx"),
        filt_indels = resultdir+"{sample}"+"_hard_filtered.indels.vcf",
        filt_indels_idx = temp(resultdir+"{sample}"+"_hard_filtered.indels.vcf.idx"),
    output:
        vcf_filt = resultdir+"{sample}"+"_filtered_variants.vcf",
        vcf_filt_idx = resultdir+"{sample}"+"_filtered_variants.vcf.idx"
    params:
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilterCombine_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{HFC_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, HFC_thrs=HFC_thrs, n_cpu=n_cpu)
    threads: HFC_thrs
    shell:
        "java -jar {params.gatk} -T CombineVariants -R {params.ref} --variant:snps {input.filt_snps} --variant:indels {input.filt_indels} -o {output.vcf_filt} -genotypeMergeOptions PRIORITIZE -priority snps,indels"

rule Annotation:
    """
    This tool uses annovar to annotate.
    """
    input:
        vcf = resultdir+"{sample}"+"_filtered_variants.vcf",
        annovar_dbs = annovar_dbs,
    output:
        k_f = temp(resultdir+"{sample}"+"_rmdup.known.exonic_variant_function"),
        n_f = temp(resultdir+"{sample}"+"_rmdup.novel.exonic_variant_function"),
        mit_rmdup = temp(resultdir+"{sample}"+".mit"),
    params:
        scripts = scripts,
        convert2annovar = convert2annovar,
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        kg_ver = kg_ver,
        mitochondrial_ver = mitochondrial_ver,
        fmt = 'vcf4',
        outdir = resultdir,
        name = "{sample}",
        err_log = True,
        maf = 0.05,
        mutect = False,
    benchmark:
        "benchmarks/benchmark_Annotation_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}"+"Annotation.py"

rule MakeFinalFile:
    """
    This step makes the final file.
    """
    input:
        k_f = resultdir+"{sample}"+"_rmdup.known.exonic_variant_function",
        n_f = resultdir+"{sample}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+"{sample}"+".mit",
    output:
        out = resultdir+"{sample}"+".tsv",
    params:
        scripts = scripts,
        mutect = False,
        sample_order = ['n','t'],
        dbsnp_freq = True,
        dbsnpFreq = None,
        dbsnpAllele = None,
        code = None, # Not necessary except for muTect
        vcf = None,  # Not necessary except for muTect
    benchmark:
        "benchmarks/benchmark_MakeFinalFile_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "MakeFinalFile.py"




##############################
#        mutect_pipe         #
##############################

rule muTect:
    """
    This step uses muTect 1.1.4 to call variants.
    """
    input:
        muTect = muTect,
        cosmic = cosmic,
        fixed_target = target + "_fixed.bed",
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
        tumour_bam = lambda wildcards: get_bam(wildcards.sick,'T'),
        normal_bai = lambda wildcards: get_bai(wildcards.sick,'N'),
        tumour_bai = lambda wildcards: get_bai(wildcards.sick,'T'),
    output:
        vcf = resultdir+"{sick}"+"_mutect.vcf",
        vcf_idx = resultdir+"{sick}"+"_mutect.vcf.idx",
        coverage_out = resultdir+"{sick}"+"_coverage.wig",
    params:
        dbsnp = dbsnp,
        ref = hg,
    conda:
        "envs/config_conda_muTect.yaml"
    benchmark:
        "benchmarks/benchmark_muTect_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{M_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, M_thrs=M_thrs, n_cpu=n_cpu)
    shell:
        "java -Xmx50g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --intervals {input.fixed_target} --input_file:normal {input.normal_bam} --input_file:tumor {input.tumour_bam} --vcf {output.vcf} --coverage_file {output.coverage_out}"

rule Annotation_muTect:
    """
    This tool uses annovar to annotate.
    """
    input:
        vcf = resultdir+"{sick}"+"_mutect.vcf",
        vcf_idx = resultdir+"{sick}"+"_mutect.vcf.idx",
        annovar_dbs = annovar_dbs,
    output:
        k_f = temp(resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known.exonic_variant_function"),
        n_f = temp(resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel.exonic_variant_function"),
        mit_rmdup = temp(resultdir+ 'mutect_ann/' + "{sick}"+".mit"),
    params:
        scripts = scripts,
        convert2annovar = convert2annovar,
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        kg_ver = kg_ver,
        mitochondrial_ver = mitochondrial_ver,
        fmt = 'vcf4old',
        outdir = resultdir + '/mutect_ann',
        name = "{sick}",
        err_log = True,
        maf = 0.05,
        mutect = True,
    threads: M_thrs
    benchmark:
        "benchmarks/benchmark_AnnotationM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}"+"Annotation.py"

rule MakeFinalFile_mutect:
    """
    This step makes the final file.
    """
    input:
        k_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known.exonic_variant_function",
        n_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+".mit",
    output:
        out = resultdir+ 'mutect_ann/' + "{sick}"+".tsv",
    params:
        scripts = scripts,
        mutect = True,
        sample_order = ['n','t'], # sample_order can change for muTect
        dbsnp_freq = True,
        dbsnpFreq = None,
        dbsnpAllele = None,
        code = lambda wildcards: get_code(wildcards.sick),
        vcf = resultdir+"{sick}"+"_mutect.vcf",
    benchmark:
        "benchmarks/benchmark_MakeFinalFileM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_1_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "MakeFinalFile.py"

###########################
#        run_lodn         #
###########################

rule lodn_table:
    """
    This step applies the LODn classifier on table.
    """
    input:
        infile = lambda wildcards: get_lodn_infile(wildcards.sick,'table'),
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
        normal_bai = lambda wildcards: get_bai(wildcards.sick,'N')
    output:
        outfile = resultdir+"{sick}"+"_lodn_table.tsv",
    params:
        fmt = 'table',
        scripts = scripts,
        min_n_cov = min_n_cov,
        min_t_cov = min_t_cov,
    conda:
        "envs/config_conda.yaml"
    threads: LodnT_thrs
    benchmark:
        "benchmarks/benchmark_LODntable_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{LodnT_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, LodnT_thrs=LodnT_thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "LODn_table.py"

rule lodn_vcf:
    """
    This step applies the LODn classifier on vcf.
    """
    input:
        infile = lambda wildcards: get_lodn_infile(wildcards.sick,'vcf'),
        infile_idx = lambda wildcards: get_lodn_infile(wildcards.sick,'vcf','.idx'),
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
        normal_bai = lambda wildcards: get_bai(wildcards.sick,'N')
    output:
        outfile = resultdir+"{sick}"+"_lodn.vcf",
    params:
        fmt = 'vcf',
        scripts = scripts,
        min_n_cov = min_n_cov,
        min_t_cov = min_t_cov,
    conda:
        "envs/config_conda.yaml"
    threads: Lodnv_thrs
    benchmark:
        "benchmarks/benchmark_LODnvcf_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{Lodnv_thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, Lodnv_thrs=Lodnv_thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "LODn_vcf.py"


###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

#########################################################
#                   DOWNLOADING FILES                   #
#########################################################

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

rule download_target:
    """download target from illumina"""
    output:
        target = temp(target + ".bed"),
    shell:
        "wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed && "
        "mv nexterarapidcapture_expandedexome_targetedregions.bed {output.target}"

rule fix_target:
    """fix target: remove prefix chr in first column and last column(name)"""
    input:
        target = target + ".bed",
    output:
        fixed_target = target + "_fixed.bed",
    shell:
        "sed 's/^chr//' {input.target} > no_chr.bed &&"
        "awk \'{{print $1,$2,$3}}\' no_chr.bed > {output.fixed_target} &&"
        "rm no_chr.bed"


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
#                   REFERENCE INDEXING                  #
#########################################################

rule index_bwa:
    """
    Generate the index of the reference genome for the bwa program.
    """
    input:
        hg,
    output:
        bwa_indexes,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    version: 0.1
    shell:
        "bwa index -a bwtsw {hg}"


rule index_picard:
    """
    Generate the index of the reference genome for the picard program.
    """
    input:
        hg=hg,
    output:
        hg.replace('fasta', 'dict'),
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard CreateSequenceDictionary R={input.hg} O={output}"


rule index_samtools:
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        hg=hg,
    output:
        hg+'.fai',
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input.hg} "

#########################################################
#                   CHECKING REQUIREMENTS               #
#########################################################

rule check_GATK:
    """
        This step check if GATK is present in the directory set in config file.
    """
    output:
        gatk,
    priority: 4
    shell:
        "echo 'Error. Genome Analysis ToolKit not found in softwares directory.' && "
        "exit 1"

rule check_muTect:
    """
        This step check if muTect is present in the directory set in config file.
    """
    output:
        muTect,
    priority: 3
    shell:
        "echo 'Error. muTect not found in softwares directory.' && "
        "exit 1"

rule check_Annovar:
    """
        This step check if annovar is present in the directory set in config file.
    """
    output:
        annovar,
    priority: 2
    shell:
        "echo 'Error. Annovar not found in softwares directory.' && "
        "exit 1"

#########################################################
#                        THE END                        #
#########################################################
