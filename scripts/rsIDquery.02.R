# 180412

#######################################################################
# This is a modified version of 'rsIDquery.01.R' that can be used     #
# to query very long rsID lists. The problem with long lists is that  #
# ncbi_snp_query sometimes gives random HTTP errors than can be fixed #
# by adjusting the size of the groups in which we split the list of   #
# rsIDs to query.                                                     #
#######################################################################

### USAGE EXAMPLE ###
# 'Rscript arg1 arg2 arg3 arg4 arg5'
# arg1 = input file (the list of rsIDs to query)
# arg2 = desired output file
# arg3 = size of groups to split the list (I suggest to start with 300 and reduce if errors occur)
# arg4 = position in rsID list to start from (0 = start from the beginning)
# if last successful batch was "from x to y", then must be arg4 = y 


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
q = as.numeric(args[3]) # usually around 298~300
start = as.numeric(args[4])

# load rsID list from file
rsIDs = read.delim(input_table, sep="\t", header=F, quote="", stringsAsFactors=F)[,1]

n = length(rsIDs)
cat(n, 'total rsIDs\n')

if(n > q){
  
  rsIDs_annotated=NULL
  n = n - start

  cat(n, 'remaining rsIDs\n')
  cat('starting from', start, '\n\n')
  
  if(start > 0){rsIDs_annotated = read.delim(outtable, sep="\t", header=T, quote="", stringsAsFactors=F)}
  
  groups = (n%/%q)
  if(n%%q > 0){groups = groups + 1}
  
  for(i in 1:groups){
    
    l1 = (q*(i-1))+1
    if(i<groups){l2 = q*i}else{l2 = n}

    l1 = l1 + start
    l2 = l2 + start

    cat('group', i, 'of', groups, '- from', l1, 'to', l2, '\n\n')
    # print(rsIDs[(l1-1):(l1+1)])
    # print(rsIDs[(l2-1):(l2+1)])

    ids = rsIDs[l1:l2]
    rsIDs_annotated<-rbind(rsIDs_annotated, ncbi_snp_query(ids))
    write.table(rsIDs_annotated, file=outtable, sep="\t", row.names=F, col.names=T, dec=".", quote=F)
    
  }

}else{rsIDs_annotated<-ncbi_snp_query(rsIDs)}

# write MAF annotation to file
write.table(rsIDs_annotated, file=outtable, sep="\t", row.names=F, col.names=T, dec=".", quote=F)



writeLines('\n########################')
writeLines('### SESSION INFO #######')
writeLines('########################')
sessionInfo()



# i = 600
# n = length(rsIDs)
# groups = (n%/%q)+1
# l1 = (q*(i-1))+1
# if(i<groups){l2 = q*i}else{l2 = n}
# 
# ids = rsIDs[l1:l2]
# print(ids)
# #ids = rsIDs[(l1+1):(l2+1)]
# #print(ids)
# 
# ids = rsIDs[l1:l2]
# annot = ncbi_snp_query(ids)
# print('annot_done')
# 
# ids = rsIDs[(l1+1):(l2+1)]
# print(ids)
# annot = ncbi_snp_query(ids)
# print('annot_done')

#write.table(ids, file='ids_599.tsv', sep="\t", row.names=F, col.names=F, dec=".", quote=F)
#write.table(annot, file=outtable, sep="\t", row.names=F, col.names=T, dec=".", quote=F)
