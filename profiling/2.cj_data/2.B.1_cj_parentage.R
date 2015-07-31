### Jinliang Yang modified from Kate Crosby
# You'll need this amount of memory to open the parentage.Rdata file
# srun -p bigmemh --mem 64000 --pty R


#library(plyr)
ped <- read.table("/group/jrigrp4/phasing/cj_teosinte/parents.txt", header =TRUE)

ped$parent1 <- as.character(ped$parent1)
ped$parent2 <- as.character(ped$parent2)
selfer <- subset(ped, parent1 == parent2)

ox <- subset(ped, parent1 != parent2)


pinfo <- data.frame(table(selfer$parent1))
names(pinfo) <- c("founder", "nselfer")

oxinfo <- data.frame(table(c(ox$parent1, ox$parent2)))
names(oxinfo) <- c("founder", "nox")

pinfo <- merge(pinfo, oxinfo, by="founder", all=TRUE)
pinfo$sid <- gsub("_1.*", "", pinfo$founder)
pinfo$sid <- gsub("_mrg.*", "", pinfo$sid)

wgs <- read.csv("data/wgs_teo19_id.csv")


pinfo <- merge(pinfo, wgs, by.x="sid", by.y="id", all=TRUE)

sub <- subset(pinfo, !is.na(WGS))






