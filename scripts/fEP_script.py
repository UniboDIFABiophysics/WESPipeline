def fEP(tsv, output_table, maf):
    # marco
    import pandas as pd

    # load tables

    variants = pd.read_table(tsv, sep='\t', header=0)

    # keep exonic
    variants = variants[variants.exon == 'exonic']

    # remove synonymous
    variants = variants[variants.var_type != 'synonymous SNV']


    # remove polymorphisms according to 1000g, ESP and EXAC
    remove = variants[(variants['1000g'] >= 0.01) | (variants.esp >= 0.01) | (variants.exac >= 0.01)].index
    variants = variants.drop(remove)


    if (maf=="/"):
        # select rsIDs and write to file
        rsID = variants.dbsnp_id[pd.notnull(variants.dbsnp_id)]
        rsID.to_csv(output_table, sep='\t', header=False, index=False)
    else:
        # load output of R script
        rsID = pd.read_table(maf, sep='\t', header=0)

        # if for any rsID a MAF value was retrieved
        if len(rsID) > 0:

            # select rsIDs with MAF >= 0.01
            rsID_remove = rsID[rsID.MAF >= 0.01]['Query']

            # get indeces of those rsIDs in variants table
            remove = [i for i in variants.index if variants.loc[i, 'dbsnp_id'] not in rsID_remove]

            # remove them
            variants = variants.drop(remove)

        # write to file
        variants.to_csv(output_table, sep='\t', header=True, index=False)


fEP(snakemake.input['tsv'],
    snakemake.output['out'],
    snakemake.input['maf'])
