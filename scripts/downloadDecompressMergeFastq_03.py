def downloadDecompressMergeFastq(name, processpath, pairs, fastq_path, fastq_cs):

    import os
    import re
    import subprocess as sp

    # create token file if --keep-going argument was specified in SNAKEMAKE
    if os.path.exists(processpath + 'failed/'):
        failedpath = processpath + 'failed/'
        for p in pairs:
            if re.search(name, p):
                sp.run('> %s%s' %(failedpath, p), shell=True)

    inpath = fastq_path
    inoutpath = processpath + '01_fastq/01_preprocess/'
    outpath = processpath + '01_fastq/02_postprocess/'

    if not os.path.exists(inoutpath):
         os.makedirs(inoutpath)

    if not os.path.exists(outpath):
         os.makedirs(outpath)

    # create fastq names
    fastq1 = name + '_R1.fastq'
    fastq2 = name + '_R2.fastq'

    # list all fastq.gz files that contain fastq_cs
    gz_list = [gz for gz in os.listdir(inpath) if re.search(fastq_cs, gz)]

    # loop over them
    for gz in gz_list:

        # decompress it while downloading to inoutpath
        sp.run(' '.join(['gzip',
                          '-dc', inpath + gz,
                          '>', inoutpath + gz[:-3]]), shell=True)


    # list all fastq for current sample
    current_sample_fastq = [f for f in os.listdir(inoutpath) if re.search(fastq_cs, f)]

    # set all possible denominations of the paired read
    denominations = ['_R', '_sequence_', '_read']

    # loop over
    for d in denominations:

        # for current sample, list all fastq that contain current denomination
        l = [f for f in current_sample_fastq if re.search(d, f)]

        # if any exist
        if l:

            # list common substrings of each pair of fastq
            common_between_read_pairs = [d.join(re.split(d + '[12]', f)) for f in current_sample_fastq]

            # uniq list
            cbrp = sorted(list(set(common_between_read_pairs)))

            break

    # list all R1 fastq (unordered), concatenating their relative path
    R1 = [inoutpath + f for f in current_sample_fastq if re.search(d + '1', f)]

    # list all R2 fastq (unordered), concatenating their relative path
    R2 = [inoutpath + f for f in current_sample_fastq if re.search(d + '2', f)]

    # create new lists
    R1b = []
    R2b = []

    # loop over common substrings (pairs of fastqs)
    for c in cbrp:

        # re-create original string
        c = (d + '[12]').join(c.split(d))

        # loop over R1 fastq
        for r1 in R1:

            # if common substring matches R1 fastq
            if re.search(c, r1):

                # append it to new list
                R1b.append(r1)

        # loop over R1 fastq
        for r2 in R2:

            # if common substring matches R1 fastq
            if re.search(c, r2):

                # append it to new list
                R2b.append(r2)

    # join elements in new lists
    R1b = ' '.join(R1b)
    R2b = ' '.join(R2b)

    # concatenate multiple fastq (if there are)
    sp.run(' '.join(['cat', R1b, '>', outpath + fastq1]), shell=True)
    sp.run(' '.join(['cat', R2b, '>', outpath + fastq2]), shell=True)

    # remove pre-processed fastq
    R = ' '.join(R1 + R2)
    sp.run('rm %s' %R, shell=True)


# Snakemake call

downloadDecompressMergeFastq(snakemake.params['name'],
                            snakemake.params['processpath'],
                            snakemake.params['pairs'],
                            snakemake.params['fastq_path'],
                            snakemake.params['fastq_cs'])
