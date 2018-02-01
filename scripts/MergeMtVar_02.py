def Merge_MuTect_Varscan(mutecttable, vsntable, logpath, match, outfile, rsID_table):
    import pandas as pd

    ###### MUTECT

    # load lists of somatic filtered variants for mutect and varscan
    mutect = pd.read_table(mutecttable, sep='\t', header=0)

    # specify contig column as 'str'
    mutect = mutect.astype({'contig': 'str'})

    # load annotated file
    mutect_anno_genome = pd.read_table(logpath + match + '_mutect_genome.hg19_multianno.txt', sep='\t', header=0, dtype='str')
    mutect_anno_mt_g = pd.read_table(logpath + match + '_mutect_mt_g.GRCh37_MT_multianno.txt', sep='\t', header=0, dtype='str')

    # select common columns
    mutect_anno_genome = mutect_anno_genome[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                                             'ExonicFunc.refGene', 'AAChange.refGene', 'cytoBand', 'genomicSuperDups',
                                             'esp6500siv2_all', 'snp138', '1000g2014oct_all', 'ExAC_nontcga_ALL',
                                             'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score',
                                             'MutationAssessor_pred']]

    mutect_anno_mt_g = mutect_anno_mt_g[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene',
                                             'ExonicFunc.ensGene', 'AAChange.ensGene']]

    # compute the difference of columns
    cols_diff = len(mutect_anno_genome.columns) - len(mutect_anno_mt_g.columns)

    # add to MT annovar table the missing columns
    mutect_anno_mt_g = pd.concat([mutect_anno_mt_g, pd.DataFrame(columns=['x']*cols_diff)], axis=1)

    # assign same column name
    mutect_anno_mt_g.columns = mutect_anno_genome.columns

    # merge by row
    mutect_anno = pd.concat([mutect_anno_genome, mutect_anno_mt_g], ignore_index=True)

    # merge total annovar output with mutect output by column
    mutect = pd.concat([mutect_anno,
                        pd.DataFrame({'somatic_status': [float('nan') * len(mutect)]}),
                        pd.DataFrame({'n_vaf': mutect['n_alt_count'] / (mutect['n_ref_count'] + mutect['n_alt_count'])}),
                        pd.DataFrame({'n_depth': mutect['n_ref_count'] + mutect['n_alt_count']}),
                        pd.DataFrame({'t_vaf': mutect['t_alt_count'] / (mutect['t_ref_count'] + mutect['t_alt_count'])}),
                        pd.DataFrame({'t_depth': mutect['t_ref_count'] + mutect['t_alt_count']}),
                        pd.DataFrame({'detection_method': ['mutect'] * len(mutect)}),
                       ], axis=1)

    ###### VARSCAN

    ### VARSCAN

    # load lists of somatic filtered variants for mutect and varscan
    varscan = pd.read_table(vsntable, sep='\t', header=0)

    # set conting values to 'str'
    varscan = varscan.astype({'chrom': 'str'})

    # convert 'tumor_var_freq' from percentage to proportion
    varscan['tumor_var_freq'] = [float(freq[:-1])/100 for freq in varscan['tumor_var_freq']]

    # load annotated files
    varscan_anno_genome = pd.read_table(logpath + match + '_varscan_genome.hg19_multianno.txt', sep='\t', header=0, dtype='str')
    varscan_anno_mt_g = pd.read_table(logpath + match + '_varscan_mt_g.GRCh37_MT_multianno.txt', sep='\t', header=0, dtype='str')

    # select common columns
    varscan_anno_genome = varscan_anno_genome[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                                             'ExonicFunc.refGene', 'AAChange.refGene', 'cytoBand', 'genomicSuperDups',
                                             'esp6500siv2_all', 'snp138', '1000g2014oct_all', 'ExAC_nontcga_ALL',
                                             'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score',
                                             'MutationAssessor_pred']]

    varscan_anno_mt_g = varscan_anno_mt_g[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene',
                                             'ExonicFunc.ensGene', 'AAChange.ensGene']]

    # compute the difference of columns
    cols_diff = len(varscan_anno_genome.columns) - len(varscan_anno_mt_g.columns)

    # add to MT annovar table the missing columns
    varscan_anno_mt_g = pd.concat([varscan_anno_mt_g, pd.DataFrame(columns=['x']*cols_diff)], axis=1)

    # assign same column name
    varscan_anno_mt_g.columns = varscan_anno_genome.columns

    # merge by row
    varscan_anno = pd.concat([varscan_anno_genome, varscan_anno_mt_g], ignore_index=True)

    # merge total annovar output with varscan output by column
    varscan = pd.concat([varscan_anno,
                         varscan['somatic_status'],
                         pd.DataFrame({'n_vaf': varscan['normal_reads2']/(varscan['normal_reads1'] + varscan['normal_reads2'])}),
                         pd.DataFrame({'n_depth': varscan['normal_reads1'] + varscan['normal_reads2']}),
                         pd.DataFrame({'t_vaf': varscan['tumor_reads2']/(varscan['tumor_reads1'] + varscan['tumor_reads2'])}),
                         pd.DataFrame({'t_depth': varscan['tumor_reads1'] + varscan['tumor_reads2']}),
                         pd.DataFrame({'detection_method': ['varscan'] * len(varscan)}),
                        ], axis=1)


    # set common column names
    colnames = ['chr', 'start', 'end', 'ref', 'alt', 'exon', 'gene', 'var_type', 'aa_change', 'cytoband', 'segdups', 'esp',
                'dbsnp_id', '1000g', 'exac', 'sift_score', 'sift_pred', 'poly_pred', 'mutasses_pred', 'somatic_status', 'n_vaf',
                'n_depth', 't_vaf', 't_depth', 'detection_method']

    mutect.columns = colnames
    varscan.columns = colnames



    #### MERGE MUTECT WITH VARSCAN

    # set columns to merge
    colnames = ['chr', 'start', 'end', 'ref', 'alt']
    # concatenate the strings in these columns and save into columns 'variant'
    mutect['variant'] = ['_'.join([str(e) for e in mutect.loc[i, colnames]]) for i in mutect.index]
    varscan['variant'] = ['_'.join([str(e) for e in varscan.loc[i, colnames]]) for i in varscan.index]

    # set column 'variant' as index
    mutect.index = mutect.variant
    varscan.index = varscan.variant

    # identify shared and unique variants between mutect and varscan lists
    unique_mutect = [v for v in mutect.variant if v not in varscan.variant]
    common = [v for v in varscan.variant if v in mutect.variant]

    # if a variant in varscan is in 'common', change its 'detection_method' field to 'both'
    varscan.loc[common, 'detection_method'] = 'both'

    # remove common variants from mutect
    mutect = mutect.loc[unique_mutect, :]

    # merge
    all_var = pd.concat([mutect, varscan], ignore_index=True)

    # write to file
    all_var.to_csv(outfile, sep='\t', header=True, index=False)

    # get rsIDs and save to file
    rsID = all_var.dbsnp_id[pd.notnull(all_var.dbsnp_id)]
    rsID.to_csv(rsID_table, sep='\t', header=False, index=False)

Merge_MuTect_Varscan(snakemake.input['m_tsv'],
                     snakemake.input['vsn_tsv'],
                     snakemake.params['workdir'],
                     snakemake.params['name'],
                     snakemake.output['table_out'],
                     snakemake.output['rsID_table'])
