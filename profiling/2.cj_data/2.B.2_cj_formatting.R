### Jinliang Yang
### 8/4/2015

save(pa, file= "largedata/cj_parentage.RData")

#ProgenyArray object: 598043 loci, 70 parents, 4805 progeny
#number of chromosomes: 11
#object size: 11640.408 Mb
#number of progeny: 4805
#number of parents: 70
#proportion missing:
#    progeny: 0.503
#parents:  0.393
#number of complete parental loci: 9202

###########################################################################

map <- as.data.frame(pa@ranges)

geno1 <- pa@parents_geno
geno1[is.na(geno1)] <- 3
geno1 <- t(geno1)

genodf <- data.frame(fid=1:70, iid=row.names(geno1), pid=0, mid=0, as.data.frame(geno1))



##################

getwgs <- function(){
    library("data.table", lib="~/bin/Rlib/")
    genoteo <- fread("/group/jrigrp4/phasing/cj_teosinte/genotypes_teosinte_19_noScaffolds_or_organelles.geno", header=FALSE)
    genoteo <- as.data.frame(genoteo) # 96,908,505
    
    genoteo$snpid <- paste(genoteo$V1, genoteo$V2, sep="_")
    
    #### comparing the alt and ref alleles
    v <- read.table("data/GBS_teos_moms.v3.txt", header=TRUE)
    v <- subset(v, !is.na(v3_coord))
    v$v3_chr <- gsub("chr", "", v$v3_chr)
    v$snpid2 <- paste(v$v2_chr, v$v2_coord, sep="_")
    v$snpid3 <- paste(v$v3_chr, v$v3_coord, sep="_")
    
    steo <- subset(genoteo, snpid %in% v$snpid3)
    message(sprintf("total [ %s ] and GBS set [ %s ]", nrow(genoteo), nrow(steo)))
    ##total [ 96908505 ] and GBS set [ 396818 ]
    steo <- steo[, -24]
    
    ids <- read.csv("data/wgs_teo19_id.csv")
    names(steo) <- c("chr", "pos", "major", "minor", as.character(ids$id), "snpid")
    
    steo <- steo[, c("snpid", "chr", "pos", "major", "minor", as.character(ids$id))]
    save(file="largedata/wgs_teo19.RData", list=c("steo", "v"))
    
}


comp_allels <- function(){
    ob <- load("largedata/wgs_teo19.RData")
    
    info <- read.csv("largedata//teo_info.csv")
    info$snpid <- gsub("S", "", info$snpid)
    info <- merge(info, v[, 5:6], by.x="snpid", by.y="snpid2")
    comp <- merge(steo[, c("snpid", "major", "minor")], info[, c(11, 2:3)], by.x="snpid", by.y="snpid3")
    
    message(sprintf("###>>> WGS [ %s ] | GBS [ %s ] | shared [ %s ]", nrow(steo), nrow(info), nrow(comp)))
    
    ### Teo19 WGS V3 and V4 are major/minor
    idx <- which((comp$major == comp$ref & comp$minor == comp$alt) | (comp$major == comp$alt & comp$minor == comp$ref))
    message(sprintf("###>>> consistent SNP calling [ %s ]", length(idx)))
    
}










