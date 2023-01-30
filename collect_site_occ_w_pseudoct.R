library(rjson)
library(tidyverse)

# get input parameters from json
input <- fromJSON(file ="input.json")

# list of variant output to process, Delta2='Late Delta'
variants <- c('Alpha','Delta','Delta2','BA.1','BA.2','BA.4','BA.5','BA.2.12.1')
# list of canonic sites to consider
refsites <- c(21552,28255,26237,26468,27041,27385,27884,25382)

outfileprefix<- input$outfileprefix

# expected folder structure is <workdir>/results/
filelist <- list.files(path='results',pattern=outfileprefix,full.names=TRUE)

allres_sites<- c()
for (variant in variants){
    infile <- filelist[grep(paste0(variant,'_'),filelist)] 
M <- c()
for (f in infile){
Mtmp <- read.table(file=f,sep='\t',header=TRUE)
M<- rbind(M,Mtmp[,which(colnames(Mtmp) %in% setdiff(colnames(Mtmp),'mappedreads'))])
}
res <- data.frame(pos=numeric(),count=numeric(),sample=character())
for (s in unique(M$sample)){
  denom <- unique(M[ M$sample ==s & M$pos < 67  & M$isize> 74 & M$isize < 500,c('qname','pos')])
  tmp <- unique(M[M$mpos>21500 & M$sample ==s,c('qname','mpos')])
  T <- table(tmp$mpos)
  T <- rbind(names(T),T)
  T <- t(T)
  T <- cbind(T,rep(as.character(s),nrow(T)),rep(nrow(denom),nrow(T)),rep(as.character(variant),nrow(T)))
  
    colnames(T) <- c('pos','count','sample','leader_depth','lineage')
  res <- rbind(res,T)
    if (nrow(T)>0 & length(setdiff(refsites,unique(T[,1])))>0){
       reszero <- c()
       for (p in setdiff(refsites,unique(T[,1]))){
          reszero <- rbind(reszero,c(as.character(p),'0',unique(T[,3:5])))
       }
       colnames(reszero) <- c('pos','count','sample','leader_depth','lineage')
       res <- rbind(res,reszero)
    }    
  }

# extend sites to the AUG codon position of each ORF  
res_sites <- res[res$pos %in% c(21552:21562,28250:29903,26237:26244,26468:26522,27041:27201,27385:27393,27884:27893,25382:25392),]

res_sites$count <- as.numeric(res_sites$count) +1
res_sites$leader_depth <- as.numeric(res_sites$leader_depth)
res_sites$pos <- as.numeric(res_sites$pos)

allres_sites <- rbind(allres_sites,res_sites)
}

# get mapped reads info
mappedreadstab <- input$mappedreadstab
mappedall <- read.table(file=mappedreadstab,header=TRUE)

allres_sites2 <- merge(allres_sites,mappedall,by='sample',all.x=TRUE)
allres_sites2 <- unique(allres_sites2)

# write output
write.table(file='results/sites_PE_occ_all_allR1_allsites_w_Nmapped_expanded_pseudocount_N_tilAUG.tsv',allres_sites2,col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)


