
def MergeWriteVariants(muTect, muTect_MT, varscansnp,varscansnp_MT,varscanindel,varscanindel_MT, match, outpath):

    # marco

    import pandas as pd
    import os

    if not os.path.exists(outpath):
         os.makedirs(outpath)

    # load mutect and varscan variants, for both genome and MT
    mutect = pd.read_table(muTect, sep='\t', header=1, dtype='str')
    mutect_MT = pd.read_table(muTect_MT, sep='\t', header=1, dtype='str')

    varscan_snp = pd.read_table(varscansnp, sep='\t', header=0, dtype='str')
    varscan_snp_MT = pd.read_table(varscansnp_MT, sep='\t', header=0, dtype='str')

    varscan_indel = pd.read_table(varscanindel, sep='\t', header=0, dtype='str')
    varscan_indel_MT = pd.read_table(varscanindel_MT, sep='\t', header=0,
                                     dtype='str')


    # remove 'MT' variants called by 'genomic branch' of the pipeline
    mutect = mutect.loc[mutect.contig != 'MT',]
    varscan_snp = varscan_snp.loc[varscan_snp.chrom != 'MT',]
    varscan_indel = varscan_indel.loc[varscan_indel.chrom != 'MT',]


    #### IN CASE VAR-FREQ COLUMNS IN VARSCAN HAVE COMMA AS DECIMAL SIGN, INSTEAD OF DOT
    var_file = [varscan_snp, varscan_snp_MT, varscan_indel, varscan_indel_MT]
    for vf in var_file:
        vf['normal_var_freq'] = [n.replace(',', '.') if ',' in n else n for n in vf['normal_var_freq']]
        vf['tumor_var_freq'] = [t.replace(',', '.') if ',' in t else t for t in vf['tumor_var_freq']]

    # merge genome and MT variants
    mutect = pd.concat([mutect, mutect_MT], ignore_index=True)

    varscan_snp = pd.concat([varscan_snp, varscan_snp_MT], ignore_index=True)

    varscan_indel = pd.concat([varscan_indel, varscan_indel_MT], ignore_index=True)


    # merge all varscan variants
    varscan_all = pd.concat([varscan_snp, varscan_indel], ignore_index=True)


    # write to file UNFILTERED variants for mutect and varscan(snv + indel)
    mutect.to_csv(outpath + match + '_mutect.tsv', sep='\t', header=True, index=False)

    varscan_snp.to_csv(outpath + match + '_varscan_snv_temp1.tsv', sep='\t', header=True, index=False)
    varscan_indel.to_csv(outpath + match + '_varscan_indel_temp1.tsv', sep='\t', header=True, index=False)
    varscan_all.to_csv(outpath + match + '_varscan.tsv', sep='\t', header=True, index=False)


MergeWriteVariants(snakemake.input['mutect_genome'],
                   snakemake.input['mutect_MT'],
                   snakemake.input['snv_genome'],
                   snakemake.input['snv_MT'],
                   snakemake.input['indel_genome'],
                   snakemake.input['indel_MT'],
                   snakemake.params['name'],
                   snakemake.params['outdir'])
