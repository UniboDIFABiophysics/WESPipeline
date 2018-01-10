def downloadDecompressMergeFastq(name, all_fastq, fastq_path, fastq_id, homepath, preprocesspath,postprocesspath,logpath,log, processpath ,cadaver=False):

    import os
    import re
    import subprocess as sp

    inpath = fastq_path
    inoutpath = preprocesspath
    outpath = postprocesspath
    logpath = logpath

    if not os.path.exists(inoutpath):
         os.makedirs(inoutpath)

    if not os.path.exists(outpath):
         os.makedirs(outpath)

    if not os.path.exists(logpath):
         os.makedirs(logpath)

    # create fastq names
    fastq1 = name + '_R1.fastq'
    fastq2 = name + '_R2.fastq'


    # in case 'cadaver' option was chosen
    if cadaver == True:

        inpath = inpath.replace('NAS/', '')

        # import list of all fastq
        with open(all_fastq) as infile:
            all_fastq = infile.read().splitlines()

        # list all fastq.gz files that contain fastq_id
        gz_list = [gz for gz in all_fastq if re.search(fastq_id, gz)]

        # open file connection
        with open(processpath + 'cadaver_get.sh', 'w') as outfile:

            # loop over fastq files of current sample
            for gz in gz_list:

                # add a 'get' command for each fastq file, each in a different line
                outfile.write(' '.join(['get',
                                        inpath + gz,
                                        inoutpath + gz + '\n']))

            # add last line to close connection
            outfile.write('quit')

        # run cadaver script
        sp.call(' '.join(['cadaver',
                          'https://ngs-ptl.unibo.it:5006',
                          '-r', processpath + 'cadaver_get.sh',
                          '>', log]), shell=True)

        # delete script
        sp.call(' '.join(['rm', processpath + 'cadaver_get.sh']), shell=True)

        # loop over fastq
        for gz in gz_list:

            # decompress it
            sp.call(' '.join(['gzip',
                              '-d', inoutpath + gz]), shell=True)

    else:

        # complete fastq_path
        fastq_path = homepath + fastq_path


        # list all fastq.gz files that contain fastq_id
        gz_list = [gz for gz in os.listdir(inpath) if re.search(fastq_id, gz)]

        # loop over them
        for gz in gz_list:

            # decompress it while downloading to inoutpath
            index = (os.path.splitext(os.path.basename(inoutpath + gz[:-3]))[0]).split('_')[1]
#% TODO: what is this step for??

            # decompress it while downloading to inoutpath
            sp.call(' '.join(['gzip',
                              '-dc', inpath + gz,
                              '>', inoutpath + gz[:-3]]), shell=True)


    # list all fastq for current sample
    current_sample_fastq = [f for f in os.listdir(inoutpath) if re.search(fastq_id, f)]

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
    sp.call(' '.join(['cat', R1b, '>', outpath + fastq1]), shell=True)
    sp.call(' '.join(['cat', R2b, '>', outpath + fastq2]), shell=True)

    # remove pre-processed fastq
    R = ' '.join(R1 + R2)
    sp.run('rm %s' %R, shell=True)


# Snakemake call

downloadDecompressMergeFastq(snakemake.params['name'],
                            snakemake.input['all_fastq'],
                            snakemake.params['fastq_path'],
                            snakemake.params['fastq_id'],
                            snakemake.params['homepath'],
                            snakemake.params['preprocesspath'],
                            snakemake.params['postprocesspath'],
                            snakemake.params['logpath'],
                            snakemake.log['log'],
                            snakemake.params['processpath'],
                            snakemake.params['cadaver'])
