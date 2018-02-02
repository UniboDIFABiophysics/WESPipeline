def labelFilterVariants(annotated, rsID_maf, black_list, processpath, pair, labelled, filtered):

    import pandas as pd
    import subprocess as sp
    import os

    # load annotated variants
    variants = pd.read_table(annotated, sep='\t', header=0)

    # load black list of genes
    with open(black_list) as infile:
        bl = infile.read().splitlines()

    # assign 'yes' label to variants that belong to genes in black list
    variants['black_list'] = ['yes' if row.gene in bl else 'no' for _, row in variants.iterrows()]

    # assign 'yes' label to variants that belong to 3 gene families
    for prefix in ['^MUC', '^OR', '^KRT']:
        sel = [i for i, row in variants.iterrows() if re.search(prefix, row.gene)]
        variants.loc[sel, 'black_list'] = 'yes'


    # get indeces of polymorphisms according to 1000g, ESP and EXAC
    polym_1 = variants[(variants['1000g'] >= 0.01) | (variants.esp >= 0.01) | (variants.exac >= 0.01)].index.tolist()

    # load table of rdIDs and their MAFs (output of R script)
    rsID_maf = pd.read_table(rsID_maf, sep='\t', header=0)

    # if for any rsID a MAF value was retrieved
    if len(rsID_maf) > 0:

        # select rsIDs with MAF >= 0.01
        rsID_remove = rsID_maf[rsID_maf['MAF'] >= 0.01]['Query']

        # select the indeces of their rows
        polym_2 = [i for i, row in variants.iterrows() if row.dbsnp_id in rsID_remove.tolist()]

    # get complete list of indeces of polymorphic variants
    polymorphisms = list(set(polym_1 + polym_2))

    # assign 'yes' label to polymorphisms and 'no' the others
    variants['polymorphism'] = ['yes' if i in polymorphisms else 'no' for i, row in variants.iterrows()]

    # write LABELLED variants to file
    variants.to_csv(labelled, sep='\t', header=True, index=False)


    # keep exonic
    variants = variants[variants.exon == 'exonic']

    # remove synonymous
    variants = variants[variants.var_type != 'synonymous SNV']

    # remove polymorphisms
    variants = variants[variants.polymorphism == 'no']

    # remove genes in black_list
    variants = variants[variants.black_list == 'no']

    # remove variants in SEGDUPS
    to_drop = [i for i, row in variants.iterrows() if pd.notnull(row['segdups'])]
    variants = variants.drop(to_drop)

    # write FILTERED variants to file
    variants.to_csv(filtered, sep='\t', header=True, index=False)



    # delete token file if --keep-going argument was specified in SNAKEMAKE
    if os.path.exists(processpath + 'failed/'):
        failedpath = processpath + 'failed/'
        os.remove(failedpath + pair)




labelFilterVariants(snakemake.input['annotated'],
                    snakemake.input['rsID_maf'],
                    snakemake.params['black_list'],
                    snakemake.params['processpath'],
                    snakemake.params['pair'],
                    snakemake.output['labelled'],
                    snakemake.output['filtered'])
