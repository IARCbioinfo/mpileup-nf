#! /usr/bin/env Rscript

indiv_run=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
indiv_run[,2]=make.unique(indiv_run[,2],sep="_")

pileups_files=paste(indiv_run[,1],".txt",sep="")
nindiv=nrow(indiv_run)
npos=length(readLines(pileups_files[1]))-1
atcg_matrix=matrix(nrow=npos,ncol=8*nindiv)
coverage_matrix=matrix(nrow=npos,ncol=nindiv + 2)

for (k in 1:nindiv) {
  cur_data=read.table(pileups_files[k],header = T,stringsAsFactors = F,sep="\t",colClasses = c("character","numeric","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character"))
    if (k==1) {
    pos_ref=cur_data[,1:2]
  }
  #ins_coverage = sapply(1:npos, function(pos) sum(as.numeric(unlist(strsplit(cur_data[pos,12], "[ | | :]"))[c(T,F)])) )
  #del_coverage = sapply(1:npos, function(pos) sum(as.numeric(unlist(strsplit(cur_data[pos,13], "[ | | :]"))[c(T,F)])) )
  coverage_matrix[,k + 2] = sapply(1:npos, function(pos) rowSums(as.matrix(cur_data[,4:11]))[pos])
}

colnames(coverage_matrix) = c("chr", "pos", indiv_run[,2])
coverage_matrix[,1] = pos_ref$chr
coverage_matrix[,2] = pos_ref$loc

output_file_name = paste(paste(paste(coverage_matrix[1,1],coverage_matrix[1,2],sep="_"),paste(coverage_matrix[dim(coverage_matrix)[1],1],coverage_matrix[dim(coverage_matrix)[1],2],sep="_"),sep="-"),"_regions_cov.txt",sep="")
write.table(as.data.frame(coverage_matrix), file = output_file_name, quote = F, row.names = F, sep="\t")
