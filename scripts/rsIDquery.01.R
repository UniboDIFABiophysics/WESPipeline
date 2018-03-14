# 171117


# if rsnps is installed locally, load
if("rsnps" %in% rownames(installed.packages())){
    do.call('library', list("rsnps"))

# if rsnps is not installed locally, download, then load
} else {
    install.packages("rsnps", repos = "https://cloud.r-project.org/")
    do.call("library", list("rsnps"))
 }

# store command line arguments
args = commandArgs(trailingOnly = TRUE)

input_table = args[1]
outtable = args[2]

# load rsID list from file
rsIDs = read.delim(input_table, sep="\t", header=F, quote="", stringsAsFactors=F)[,1]

# the number of rsIDs that can be passed in to one request is limited
# therefore, if there are more than 300 I split them in group of 300 and make multiple queries
q = 300

if(length(rsIDs)>q){

  rsIDs_annotated=NULL

  n = length(rsIDs)
  groups = (n%/%q)+1

  for(i in 1:groups){
    l1 = (q*(i-1))+1
    if(i<groups){l2 = q*i}else{l2 = n}
    ids = rsIDs[l1:l2]
    rsIDs_annotated<-rbind(rsIDs_annotated, ncbi_snp_query(ids))
  }

}else{rsIDs_annotated<-ncbi_snp_query(rsIDs)}

# write MAF annotation to file
write.table(rsIDs_annotated, file=outtable, sep="\t", row.names=F, col.names=T, dec=".", quote=F)


writeLines('\n########################')
writeLines('### SESSION INFO #######')
writeLines('########################')
sessionInfo()
