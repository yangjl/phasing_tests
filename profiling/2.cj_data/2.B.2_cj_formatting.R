### Jinliang Yang

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



library("data.table", lib="~/bin/Rlib/")
genoteo <- fread("/group/jrigrp4/phasing/cj_teosinte/genotypes_teosinte_19_noScaffolds_or_organelles.geno", header=FALSE)
genoteo <- as.data.frame(genoteo) # 96,908,505

genoteo$snpid <- paste(genoteo$V1, genoteo$V2, sep="_")


v <- read.table("data/GBS_teos_moms.v3.txt", header=TRUE)
v <- subset(v, !is.na(v3_coord))
v$v3_chr <- gsub("chr", "", v$v3_chr)
v$snpid2 <- paste(v$v2_chr, v$v2_coord, sep="_")
v$snpid3 <- paste(v$v3_chr, v$v3_coord, sep="_")


steo <- subset(genoteo, snpid %in% v$snpid3)


info <- read.csv("largedata//teo_info.csv")
info$snpid <- gsub("S", "", info$snpid)
info <- merge(info, v[, 5:6], by.x="snpid", by.y="snpid2")
comp <- merge(steo[, c("snpid", "V3", "V4")], info[, c(7, 2:3)], by.x="snpid", by.y="snpid3")


### Teo19 WGS V3 and V4 are major/minor
idx <- which(comp$V3 != comp$ref & comp$V4 != comp$ref)







