library(Rsamtools)
library(IRanges)
library(rjson)
library(tidyverse)
library(foreach)

# read input parameters from json
input <- fromJSON(file ="input.json")

gr <- GRanges('MN908947.3',IRanges(c(52,28250),c(67,29903)))  
p1<- ScanBamParam(which=as(gr, "IntegerRangesList"),what=c("qname","rname","flag","strand","pos","cigar","mpos","isize","seq"))
p2<- ScanBamParam(what=c("qname","rname","flag","strand","pos","cigar","mpos","isize","seq"))
leader_template <- 'GTAGATCTGTTCTCT'
minlen <- 70
N <- 5000
# which variant sub folder to run 
variant<- input$variant

# work directory
workdir <- input$workdir

# BAM files stored under the structure <workdir>/<variant>/bam/ 
setwd(paste0(workdir,'/',variant,'/bam/'))

#dir.create(paste0(workdir,'/',variant,'/bam/depth'),showWarnings=FALSE)

# read metadata file
# contains the columns sample|variant|readlen|bampath
metadatafile <- input$metadatafile

allmeta <- read.table(file = metadatafile,sep='\t',header=T)
allmeta <- allmeta[allmeta$variant == variant & allmeta$readlen > minlen,]

# get list of bamfiles
samplist <- list.files(path = '.', pattern = '.bam$')
samplist <- samplist[samplist %in% basename(allmeta$bam_path)]

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

resall <- c()

if (N > length(samplist)){
  N <- length(samplist)
}

# loop over BAM files
resall <- foreach (i=samplist[1:N], .combine = rbind ) %dopar% {
  library(Rsamtools)
  library(IRanges)
  library(tidyverse)
  tmp5 <- scanBam(file=i,param=p1)
  tmp5all <- scanBam(file=i,param=p2)

  # get all the reads within [52, 67] range  
  res5 <- tmp5[[1]]
  # get all reads
  res5all <- tmp5all[[1]]
  # get number of mapped reads
  mappedreads <- idxstatsBam(i)[1,3]  
  # get all reads within range that match the template
  qnames <- res5$qname
  mpos <- res5$mpos

  # get all reads with read name in the list of reads matching the range and containing the motif
  res5all <- as.data.frame(res5all)
  res5all <- res5all[res5all[,1] %in% qnames,]

  # order by readname
  res5all <- res5all[order(res5all[,1]),]
  # add sample name
  res5all$sample <- rep(unlist(strsplit(i,'\\.'))[1],nrow(res5all))  
  # return template match positon
  res5all$leaderpos <- res5all$pos + str_locate(res5all$seq,leader_template)[1]
  res5all$mappedreads <- rep(mappedreads,nrow(res5all))

  res5all
  
}
outfile <- paste0(input$outfileprefix,'_',variant,'_',N,'.tsv')
# write read table with columns qname|flag|rname|strand|pos|cigar|mpos|isize|seq|sample|leaderpos|mappedreads|substr_match
write.table(file=outfile,resall,col.names=T, row.names=F, sep='\t',quote=F)


doParallel::stopCluster(cl)



