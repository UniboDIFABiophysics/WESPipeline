def mutectPreAnnovar(genome_infile, mt_infile, pair, outfile):

    import pandas as pd

    # set some needed data types
    dtype = {'contig': 'str',
             'position': 'str',
             'n_ref_count': 'int',
             't_ref_count': 'int'}

    # load MUTECT tables
    genome = pd.read_table(genome_infile, sep='\t', header=1, dtype=dtype)
    mt = pd.read_table(mt_infile, sep='\t', header=1, dtype=dtype)

    # remove 'MT' variants called by 'genomic branch' of the pipeline
    genome = genome.loc[genome.contig != 'MT',]

    # merge
    mutect = pd.concat([genome, mt])

    # re-organize first columns to fit ANNOVAR
    # remove tumor and normal names columns
    mutect = pd.concat([pd.DataFrame({'chrom': mutect.contig}),
                        pd.DataFrame({'start': mutect.position}),
                        pd.DataFrame({'end': mutect.position}),
                        pd.DataFrame({'ref': mutect.ref_allele}),
                        pd.DataFrame({'alt': mutect.alt_allele}),
                        mutect.context,
                        mutect.iloc[:,7:]
                       ], axis=1)

    # compute N and T depth
    mutect['n_depth'] = mutect.n_ref_count + mutect.n_alt_count
    mutect['t_depth'] = mutect.t_ref_count + mutect.t_alt_count

    # compute N and T VAF
    mutect['n_vaf'] = mutect.n_alt_count / mutect.n_depth
    mutect['t_vaf'] = mutect.t_alt_count / mutect.t_depth

    # create 'variant' column
    mutect['variant'] = ['_'.join(row[:5].tolist() + [pair]) for i, row in mutect.iterrows()]

    # write to file
    mutect.to_csv(outfile, sep='\t', index=False)

mutectPreAnnovar(snakemake.input['genome_infile'],
                 snakemake.input['mt_infile'],
                 snakemake.params['pair'],
                 snakemake.output['outfile'])
