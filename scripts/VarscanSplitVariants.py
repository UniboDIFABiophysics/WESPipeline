def VarscanSplitVariants(vsn,inoutpath,match):
    import pandas as pd
    ### VARSCAN

    # load lists of somatic filtered variants for mutect and varscan
    varscan = pd.read_table(vsn, sep='\t', header=0)

    # set conting values to 'str'
    varscan = varscan.astype({'chrom': 'str'})

    # convert 'tumor_var_freq' from percentage to proportion
    varscan['tumor_var_freq'] = [float(freq[:-1])/100 for freq in varscan['tumor_var_freq']]

    # extract chr, start, end, ref, alt columns (renaming 'end' one)
    varscan_anno = pd.concat([varscan[['chrom', 'position']],
                     pd.DataFrame({'end': varscan['position']}),
                     varscan[['ref', 'var']]], axis=1)

    # get indeces of rows with DELETIONS and INSERTIONS
    deletion = varscan_anno[['-' in v for v in varscan_anno['var']]].index
    insertion = varscan_anno[['+' in v for v in varscan_anno['var']]].index

    ## INSERTIONS
    # change 'ref' column to '-'
    varscan_anno.loc[insertion, 'ref'] = ['-'] * len(insertion)
    # remove '-' from the beginning of strings in 'var' column
    varscan_anno.loc[insertion, 'var'] = [v[1:] for v in varscan_anno.loc[insertion, 'var']]
    # add 1 to 'start' and 'end' chromosomal positions
    varscan_anno.loc[insertion, 'position'] = varscan_anno.loc[insertion, 'position'] + 1
    varscan_anno.loc[insertion, 'end'] = varscan_anno.loc[insertion, 'end'] + 1

    ## DELETIONS
    # substitute strings in 'ref' with strings in 'var', removing the '-' at the beginning
    varscan_anno.loc[deletion, 'ref'] = [v[1:] for v in varscan_anno.loc[deletion, 'var']]
    # change 'var' column to '-'
    varscan_anno.loc[deletion, 'var'] = ['-'] * len(deletion)
    # add 1 to 'start' chromosomal positions
    varscan_anno.loc[deletion, 'position'] = varscan_anno.loc[deletion, 'position'] + 1
    # add length of string in 'ref' to 'end' chromosomal positions
    varscan_anno.loc[deletion, 'end'] = [varscan_anno.loc[i, 'end'] +
                                         len(varscan_anno.loc[i, 'ref']) for i in deletion]

    # split between genomic variants and MT
    varscan_anno_genome = varscan_anno[varscan_anno.chrom != 'MT'].copy()
    varscan_anno_mt = varscan_anno[varscan_anno.chrom == 'MT'].copy()


    # write to file
    varscan_anno_genome.to_csv(inoutpath + match + '_varscan_annovar_genome.tsv', sep='\t', header=False, index=False)
    varscan_anno_mt.to_csv(inoutpath + match + '_varscan_annovar_mt.tsv', sep='\t', header=False, index=False)

VarscanSplitVariants(snakemake.input['vsn_tsv'],
                     snakemake.params['workdir'],
                     snakemake.params['name'])
