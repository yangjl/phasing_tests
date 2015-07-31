### Jinliang Yang modified from Kate Crosby
# You'll need this amount of memory to open the parentage.Rdata file
# srun -p bigmemh --mem 64000 --pty R

#setwd("/group/jrigrp4/phasing/cj_teosinte")

ob <- load("/group/jrigrp4/phasing/cj_teosinte/parentage.Rdata")

#need the ProgenyArray, ignore
#str(pa)

# parents/kids index
parent.names <- colnames(pa@parents_geno)
progeny.names <- colnames(pa@progeny_geno)


x <- rep(1:70, 1)

parents <- data.frame(parent.names,x)

# index file

idx.all <- data.frame(pa@parents)
progeny.names <- colnames(pa@progeny_geno)

idx.all <- data.frame(progeny.names, idx.all)

idx.all <- idx.all[,1:4]

# Merge

new.df <- merge(idx.all, parents, by.x="parent_1", by.y="x")
colnames(new.df)[colnames(new.df)=="parent.names"] <- "parent1_names"

df <- merge(new.df, parents, by.x="parent_2", by.y="x")
colnames(df)[colnames(df)=="parent.names"] <- "parent2_names"

final.df <- data.frame(df$progeny.names, df$parent1_names, df$parent2_names)
colnames(final.df) <- c("progeny", "parent1", "parent2")
write.table(final.df, "parents.txt")

# Now exit out and open a new session

rm(list=ls())

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



fams[]<- lapply(fams, as.character)
x <- which(fams$parent1==fams$parent2)

new.df <- fams[row.names(fams) %in% x,]

new.df[] <- lapply(new.df, as.factor)


# These below should be identical for the selves
parent1 <- table(new.df$parent1)

sort(parent1)

parent2 <- table(new.df$parent2)

sort(parent2)


# Now get those familes, and make a list of them to retrieve
biggest <- subset(new.df, parent1=="PC_O51_ID2_mrg:250276282")
smallest <- subset(new.df, parent1=="PC_J10_ID1_2:250276214")

biggest.array.kids <- biggest[,1]
smallest.array.kids <- smallest[,1]

#Append the parent name
write.table(biggest.array.kids,"biggest_selfed_family.txt")
write.table(smallest.array.kids, "smallest_selfed_family.txt")

