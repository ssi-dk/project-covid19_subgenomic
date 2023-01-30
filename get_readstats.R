library(Rsamtools)
library(IRanges)
library(rjson)
library(tidyverse)
library(foreach)

input <- fromJSON(file ="input.json")

gr <- GRanges('MN908947.3',IRanges(52,67))
p1<- ScanBamParam(which=as(gr, "IntegerRangesList"),what=c("qname","rname","flag","strand","pos","cigar","mpos","isize","seq"))
p2<- ScanBamParam(what=c("qname","rname","flag","strand","pos","cigar","mpos","isize","seq"))
leader_template <- 'GTAGATCTGTTCTCT'
minlen <- 70
N <- 5000

variants<- c('Alpha','Delta','Delta2','BA.1','BA.2','BA.4','BA.5','BA.2.12.1')

for (variant in variants){

workdir <- input$workdir

setwd(paste0(workdir,'/',variant,'/bam/'))

metadatafile <- input$metadatafile

allmeta <- read.table(file=metadatafile,header=T, sep='\t')


allmeta <- allmeta[allmeta$variant == variant & allmeta$readlen > minlen,]



samplist <- list.files(path = '.', pattern = '.bam$')
samplist <- samplist[samplist %in% basename(allmeta$bam_path)]

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)


resall <- c()

#resall <- foreach (i=samplist[1:N], .combine = rbind ) %dopar% {
resall <- foreach (i=samplist[1:length(samplist)], .combine = rbind ) %dopar% {
  library(Rsamtools)
  library(IRanges)
  library(tidyverse)
  mappedreads <- idxstatsBam(i)[1,3]
resall <- rbind(resall,c(as.character(unlist(strsplit(i,'\\.'))[1]),mappedreads))
}
colnames(resall) <- c('sample','mapped')
write.table(file=paste0('idxstatsBam_',variant,'_',length(samplist),'.tsv'),resall,col.names=T, row.names=F, sep='\t',quote=F)

}

