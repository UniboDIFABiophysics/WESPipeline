def SplitVariants(m_som,inoutpath,match):
    import pandas as pd
    ### MUTECT

    # load lists of somatic filtered variants for mutect and varscan
    mutect = pd.read_table(m_som, sep='\t', header=0)

    # specify contig column as 'str'
    mutect = mutect.astype({'contig': 'str'})

    # extract chr, start, end, ref, alt columns
    mutect_anno = pd.concat([mutect[['contig', 'position']],
                             mutect['position'],
                             mutect[['ref_allele', 'alt_allele']]
                            ], axis=1)

    # split genome and mt variants
    mutect_anno_genome = mutect_anno[mutect_anno.contig != 'MT'].copy()
    mutect_anno_mt = mutect_anno[mutect_anno.contig == 'MT'].copy()

    # write to file
    mutect_anno_genome.to_csv(inoutpath + match + '_mutect_annovar_genome.tsv', sep='\t', header=False, index=False)
    mutect_anno_mt.to_csv(inoutpath + match + '_mutect_annovar_mt.tsv', sep='\t', header=False, index=False)

MutectSplitVariants(snakemake.input['m_tsv'],
              snakemake.params['workdir'],
              snakemake.params['name'])
