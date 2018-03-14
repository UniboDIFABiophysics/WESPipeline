import os
import argparse
    
#### Parse command line arguments and load CONFIGFILE ####

parser = argparse.ArgumentParser()
parser.add_argument('-normal', required=True, dest='normal', action='store')
parser.add_argument('-tumor', required=True, dest='tumor', action='store')

# get arguments
args = parser.parse_args()

### get CONFIG file absolute path
normal = args.normal
tumor = args.tumor
    
if os.stat(normal).st_size == 0:
    with open(normal, 'w') as outfile:
        outfile.write('MT\t1\tA\t1\t.\t1')
        
if os.stat(tumor).st_size == 0:
    with open(tumor, 'w') as outfile:
        outfile.write('MT\t1\tA\t1\t.\t1')
