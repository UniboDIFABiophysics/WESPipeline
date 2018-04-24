def varscanSomaticFilterPreAnnovar(genome_snv_infile, genome_indel_infile, mt_snv_infile, mt_indel_infile, varscan, varfiltpath, pair, min_coverage, min_reads2, min_var_freq_sf, outfile):

    import pandas as pd
    import subprocess as sp
    import os
    
    # create log folder if it doesn't already exists
    varfiltpath_logs = varfiltpath + 'logs/'
    os.makedirs(varfiltpath_logs, exist_ok=True)
    
    # load VARSCAN tables
    genome_snv = pd.read_table(genome_snv_infile, sep='\t', header=0, dtype='str')
    genome_indel = pd.read_table(genome_indel_infile, sep='\t', header=0, dtype='str')
    mt_snv = pd.read_table(mt_snv_infile, sep='\t', header=0, dtype='str')
    mt_indel = pd.read_table(mt_indel_infile, sep='\t', header=0, dtype='str')

    # remove 'MT' variants called by 'genomic branch' of the pipeline
    genome_snv = genome_snv.loc[genome_snv.chrom != 'MT',]
    genome_indel = genome_indel.loc[genome_indel.chrom != 'MT',]

    # set variable for somaticFilter (sf)
    varscan = varscan
    sf_infile_snv = varfiltpath + pair + '_varscan_snv.tsv'
    sf_infile_indel = varfiltpath + pair + '_varscan_indel.tsv'
    min_coverage = min_coverage
    min_reads2 = min_reads2
    min_var_freq_sf = min_var_freq_sf
    sf_outfile_snv = varfiltpath + pair + '_varscan_snv_somaticFilter.tsv'
    sf_outfile_indel = varfiltpath + pair + '_varscan_indel_somaticFilter.tsv'
    sf_log_snv = varfiltpath_logs + pair + '_varscan_snv_somaticFilter_err.log'
    sf_log_indel = varfiltpath_logs + pair + '_varscan_indel_somaticFilter_err.log'

    # merge GENOMIC and MT and write to file, to be input of somaticFilter
    pd.concat([genome_snv, mt_snv]).to_csv(sf_infile_snv, sep='\t', index=False)
    pd.concat([genome_indel, mt_indel]).to_csv(sf_infile_indel, sep='\t', index=False)


    # somaticFilter for SNV
    cmd = 'java -jar %s somaticFilter %s --indel-file %s --min-coverage %s --min-reads2 %s --min-var-freq %s --output-file %s 2> %s' %(varscan, sf_infile_snv, sf_infile_indel, min_coverage, min_reads2, min_var_freq_sf, sf_outfile_snv, sf_log_snv)
    sp.run(cmd, shell=True)

    # somaticFilter for INDEL
    cmd = 'java -jar %s somaticFilter %s --min-coverage %s --min-reads2 %s --min-var-freq %s --output-file %s 2> %s' %(varscan, sf_infile_indel, min_coverage, min_reads2, min_var_freq_sf, sf_outfile_indel, sf_log_indel)
    sp.run(cmd, shell=True)

    # load somaticFilter output tables
    snv_sf = pd.read_table(sf_outfile_snv, sep='\t', header=0, dtype='str')
    indel_sf = pd.read_table(sf_outfile_indel, sep='\t', header=0, dtype='str')

    # merge post-somaticFilter tables
    varscan_sf = pd.concat([snv_sf, indel_sf])

    # merge pre-somaticFilter tables
    varscan = pd.concat([genome_snv, mt_snv, genome_indel, mt_indel]).reset_index(drop=True)

    # create list of every row in varscan and varscan_sf, joined by '_'
    joined = ['_'.join(row) for i, row in varscan.iterrows()]
    joined_sf = ['_'.join(row) for i, row in varscan_sf.iterrows()]

    # label variants whether they passed or not the filter
    varscan['somaticFilter'] = ['passed' if v in joined_sf else 'not_passed' for v in joined]

    # label variants whether they have tumor strand bias or not
    varscan['tumor_strand_bias'] = ['yes' if (float(row.tumor_reads2_plus) == 0) | (float(row.tumor_reads2_minus) == 0) else 'no' for i, row in varscan.iterrows()]

    # derive N and T VAF columns
    varscan['n_vaf'] = [n.replace(',', '.') if ',' in n else n for n in varscan['normal_var_freq']]
    varscan['t_vaf'] = [t.replace(',', '.') if ',' in t else t for t in varscan['tumor_var_freq']]
    varscan['n_vaf'] = [float(freq[:-1])/100 for freq in varscan['n_vaf']]
    varscan['t_vaf'] = [float(freq[:-1])/100 for freq in varscan['t_vaf']]

    # drop old columns
    varscan = varscan.drop(['normal_var_freq', 'tumor_var_freq'], axis=1)

    # compute N and T depth
    varscan['n_depth'] = [float(row.normal_reads1) + float(row.normal_reads2) for i, row in varscan.iterrows()]
    varscan['t_depth'] = [float(row.tumor_reads1) + float(row.tumor_reads2) for i, row in varscan.iterrows()]

    # re-organize first columns to fit ANNOVAR
    varscan = pd.concat([varscan.chrom,
                         pd.DataFrame({'start': varscan.position}),
                         pd.DataFrame({'end': varscan.position}),
                         pd.DataFrame({'ref': varscan.ref}),
                         pd.DataFrame({'alt': varscan['var']}),
                         varscan.iloc[:,4:]
                        ], axis=1)
    # convert start and end to 'int'
    varscan = varscan.astype({'start': 'int', 'end': 'int'})


    # get indeces of rows with DELETIONS and INSERTIONS
    deletion = varscan[['-' in v for v in varscan.alt]].index
    insertion = varscan[['+' in v for v in varscan.alt]].index

    ## INSERTIONS
    # change 'ref' column to '-'
    varscan.loc[insertion, 'ref'] = ['-'] * len(insertion)
    # remove '-' from the beginning of strings in 'alt' column
    varscan.loc[insertion, 'alt'] = [v[1:] for v in varscan.loc[insertion, 'alt']]
    # add 1 to 'start' and 'end' chromosomal positions
    varscan.loc[insertion, 'start'] = varscan.loc[insertion, 'start'] + 1
    varscan.loc[insertion, 'end'] = varscan.loc[insertion, 'end'] + 1

    ## DELETIONS
    # substitute strings in 'ref' with strings in 'alt', removing the '-' at the beginning
    varscan.loc[deletion, 'ref'] = [v[1:] for v in varscan.loc[deletion, 'alt']]
    # change 'alt' column to '-'
    varscan.loc[deletion, 'alt'] = ['-'] * len(deletion)
    # add 1 to 'start' chromosomal positions
    varscan.loc[deletion, 'start'] = varscan.loc[deletion, 'start'] + 1
    # add length of string in 'ref' to 'end' chromosomal positions
    varscan.loc[deletion, 'end'] = [varscan.loc[i, 'end'] +
                                         len(varscan.loc[i, 'ref']) for i in deletion]

    # convert start and end to 'str'
    varscan = varscan.astype({'start': 'str', 'end': 'str'})

    # create 'variant' column
    varscan['variant'] = ['_'.join(row[:5].tolist() + [pair]) for i, row in varscan.iterrows()]

    # write to file
    varscan.to_csv(outfile, sep='\t', index=False)



    # remove intermediate files
    sp.run('rm %s%s_varscan_*.tsv' % (varfiltpath, pair), shell=True)


varscanSomaticFilterPreAnnovar(snakemake.input['genome_snv_infile'],
                               snakemake.input['genome_indel_infile'],
                               snakemake.input['mt_snv_infile'],
                               snakemake.input['mt_indel_infile'],
                               snakemake.params['varscan'],
                               snakemake.params['varfiltpath'],
                               snakemake.params['pair'],
                               snakemake.params['min_coverage'],
                               snakemake.params['min_reads2'],
                               snakemake.params['min_var_freq_sf'],
                               snakemake.output['outfile'])
