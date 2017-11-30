def KeepOnlyVariants(mutect2, varscan_snp2,varscan_indel2,outpath,match):
    import pandas as pd

    mutect = pd.read_table(mutect2, sep='\t', header=1, dtype='str')
    mutect = mutect[(mutect.judgement == 'KEEP') & (mutect.covered == 'COVERED')]
    # laod and merge filtered varscan tables
    varscan_snp = pd.read_table(varscan_snp2, sep='\t', header=0)
    varscan_indel = pd.read_table(varscan_indel2, sep='\t', header=0)

    varscan_all = pd.concat([varscan_snp, varscan_indel], ignore_index=True)


    # keep only variants that varscan called 'Somatic' or 'LOH'
    varscan_all = varscan_all[(varscan_all.somatic_status == 'Somatic') | (varscan_all.somatic_status == 'LOH')]

    # keep only variants with variant reads on both strands (for tumor)
    varscan_all = varscan_all[(varscan_all.tumor_reads2_plus > 0) & (varscan_all.tumor_reads2_minus > 0)]

    # write FILTERED mutect and varscan variants
    mutect.to_csv(outpath + match + '_mutect_somatic.tsv', sep='\t', header=True, index=False)
    varscan_all.to_csv(outpath + match + '_varscan_somatic.tsv', sep='\t', header=True, index=False)


KeepOnlyVariants(snakemake.input['m'],
                 snakemake.input['snv'],
                 snakemake.input['indel'],
                 snakemake.params['outdir'],
                 snakemake.params['name'])
