# The WES reference BED files from Illumina have the chromosome names as 'chr#', 'chrX', 'chrY', 'chrM' and before running the pipeline we must change them to '1', '2', .. , 'Y', 'MT'

import pandas as pd

def editBEDChromField(input_bed, line_to_skip,output_bed):
    # get name of BED file from command line
    bed_name = str(input_bed)

    # some files have a first line that we have to skip
    # if no 2nd argument is provided, skip no lines
    # otherwise skip the line number provided in the 2nd argument
    rows_to_skip = None
    if line_to_skip:
        rows_to_skip = [int(r) for r in line_to_skip]

    # load file
    bed = pd.read_table(bed_name, sep='\t', header=None, dtype='str', skiprows=rows_to_skip)

    # remove MT intervals
    bed = bed.loc[bed[0] != 'chrM',]

    # change 'chr#' to '#'
    bed[0] = [c[3:] for c in bed[0]]

    # write to file
    bed.to_csv(str(output_bed), sep='\t', header=False, index=False)

editBEDChromField(snakemake.input,
                  snakemake.params['line_to_skip'],
                  snakemake.output)
