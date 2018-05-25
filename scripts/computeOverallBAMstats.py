def computeOverallBAMstats(bam_stats, off_target_bed, kit_bed, overall_bam_stats):

    import subprocess as sp
    import pandas as pd

    # load file with general BAM statistics (output of 'samtools stats')
    with open(bam_stats) as infile:
        stats = infile.read().splitlines()

    # many small edits to stats table, to extract only necessary info (fields that start by 'SN')
    stats = [e.replace(':', '').split('\t') for e in stats if 'SN' in e][1:]

    # loop over fields in 'stats' and store label/value for each field
    labels = []
    dic = {}
    for e in stats:
        label = e[1].replace(' ', '_')
        dic[label] = float(e[2])
        labels.append(label)
        if label in ['reads_mapped', 'reads_properly_paired', 'reads_duplicated']:
            dic[label + '_(%)'] = float(e[2])/dic['sequences']*100
            labels.append(label + '_(%)')

    # count OFF-TARGET alignments
    cmd = 'wc -l %s' % off_target_bed
    off_target = sp.run(cmd, shell=True, stdout=sp.PIPE).stdout.decode()
    off_target = float(off_target.split(' ')[0])
    # store OFF-TARGET as new field
    dic['reads_off_target'] = off_target
    labels.append('reads_off_target')
    
    # compute and store percentage of OFF-TARGETs
    dic['reads_off_target_(%)'] = dic['reads_off_target']/dic['sequences']*100
    labels.append('reads_off_target_(%)')

    # compute total EXOME SIZE (bp), by summing up the intervals of the kit BED file
    cmd = "cat %s | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'" %kit_bed
    # store EXOME SIZE as new field
    dic['exome_size'] = float(sp.run(cmd, shell=True, stdout=sp.PIPE).stdout.decode().splitlines()[0])
    labels.append('exome_size')

    # use many info to compute BAM coverage and store it as new field
    ot = dic['reads_off_target']/dic['sequences']   # proportion of OFF-TARGETS
    bm = dic['bases_mapped']                        # total mapped bases
    mapq0 = dic['reads_MQ0']/dic['sequences']       # proportion of alignments with MAPQ=0
    d = dic['reads_duplicated']/dic['sequences']    # proportion of duplicated reads
    es = dic['exome_size']
    dic['coverage'] = bm * (1 - ot) * (1 - d) * (1 - mapq0) / es
    labels.append('coverage')

    # collect every field in a SERIES and write to file
    overall_stats = pd.Series(dic, index=labels)
    overall_stats.to_csv(overall_bam_stats, sep='\t', header=False)


computeOverallBAMstats(snakemake.input['bam_stats'],
                       snakemake.input['off_target_bed'],
                       snakemake.input['kit_bed'],
                       snakemake.output['overall_bam_stats'])
