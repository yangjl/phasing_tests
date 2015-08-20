### Jinliang Yang
### August 7th, 2015

#######################################
FormatData <- function(wgs, cols){
    
    source("lib/load_data.R")
    
    ### load WGS data of 19 teosintes and recoded to `0, 1, 2` format and `3` indicates missing.
    wgs <- recode()
    
    ## PED file
    #'family ID' 'plant ID' 'Paternal ID' 'Maternal ID' 'genotype sequence'
    wgs$snpid2 <- paste0("S", wgs$snpid2)
    #wgs <- wgs[order(wgs$chr, wgs$pos), ]
    
    ### SNP matrix comparison
    ob <- load("largedata/GBS_genomx.RData")
    
    nms <- gsub("_1\\:.*|_mrg\\:.*", "", colnames(genos))
    subgeno <- genos[rownames(genos) %in% wgs$snpid2, ]
    #subgeno <- subgeno[order(rownames(subgeno)), ]
    
    
    #############################################################################################
    #parentage pafile="data/parentage_info.txt"
    pa <- read.table("data/parentage_info.txt", header=T)
    
    pa$parent1 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent1)
    pa$parent2 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent2)
    pa <- subset(pa, parent1==parent2)
    
    
    for(i in cols){
        mid <- names(wgs)[i]
        kids <- subgeno[, subset(pa, parent1==mid)$proid]
        kids[is.na(kids)] <- 3
        mk <- merge(wgs[, c("snpid", "snpid2", "chr", "pos", mid)], as.data.frame(kids), by.x="snpid2", by.y="row.names")
        mk <- mk[order(mk$chr, mk$pos),]
        message(sprintf("###>>> Common SNPs [ %s ] of family [ %s ] with kids [ N=%s]", nrow(subgeno), mid, ncol(kids)))
        save(file=paste0("largedata/sfamdata/", mid, ".RData"), list="mk")
    }
    
}

### write PLINK format file
FormatData(wgs, cols=9:27)

### input files
files <- list.files(pattern="RData", path="largedata/sfamdata", full.names=TRUE)
tem <- data.frame(input=rep(files, each=10), chr=rep(1:10, times=19), output=0)
tem$output <- gsub(".RData", "_chr", tem$input)
tem$output <- paste0(tem$output, tem$chr, ".RData")
write.table(tem, "largedata/sfamdata/all_files.txt", sep="\t", row.names=FALSE, quote=FALSE)

### sbatch codes
tem <- data.frame(input="sbatch -p serial slurm-script/run_phasing.sh", num=1:160)
write.table(tem, "slurm-script/sbatch.txt", sep="\t", row.names=FALSE, quote=FALSE)

