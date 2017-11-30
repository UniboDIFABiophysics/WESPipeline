# 171117

# query a list of rsIDs to ncbi and get their MAF score, then save the scores to file
# usage: Rscript rsIDfilter.R arg1
# 'arg1' is required and must be the whole path of input rsID list file
# input file must be called 'rsID.tsv'

library(rsnps)

# store command line arguments
args = commandArgs(trailingOnly = TRUE)

# load rsID list from file
rsIDs = read.delim(paste0(args[1], 'rsID.tsv'), sep="\t", header=F, quote="", stringsAsFactors=F)

# the number of rsIDs that can be passed in to one request is limited to around 600
# therefore, if there are more than 500 I split them in group of 500 and make multiple queries
if(length(rsIDs)>500){
  
  rsIDs_annotated=NULL
  
  n = length(rsIDs)
  groups = (n%/%500)+1
  
  for(i in 1:groups){
    l1 = (500*(i-1))+1
    if(i<groups){l2 = 500*i}else{l2 = n}
    ids = rsIDs[l1:l2]
    rsIDs_annotated<-rbind(rsIDs_annotated, ncbi_snp_query(ids))
  }
  
}else{rsIDs_annotated<-ncbi_snp_query(rsIDs)}

# write MAF annotation to file
write.table(rsIDs_annotated, file=paste0(args[1], 'rsID_maf.tsv'), sep="\t", row.names=F, col.names=T, dec=".", quote=F)


sessionInfo()
