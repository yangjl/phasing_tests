## convert v2 to v3
## ensembl is the worst and won't even give me an error that means anything
###    they suggest using their batch script that can only be used for humans and mice
### took the coordinates from Jeff Glaubitz's genome annotation db

library(rtracklayer)
library(data.table)
setwd('/group/jrigrp4/phasing/cj_teosinte')
## coordinate conversion
## I don't care about boundary conditions
##   be careful if you worry about such things.
##   They were provided as [1, 96447)
##   I think the reason for a lot of lines in this file is that the 100 N contig boundaries were thrown out of the v2/v3 alignment
conversion=read.table('v2_v3_from_genomeAnnos_db.txt', header=T, sep='\t')
## set up genomic ranges object to do overlaps
vc=GRanges(seqnames=conversion$v2_chr, ranges=IRanges(start=conversion$v2_start, end=conversion$v2_end))
mcols(vc)$v3_chr=conversion$v3_chr
mcols(vc)$v3_start=conversion$v3_start
mcols(vc)$v3_end=conversion$v3_end


#######################################################
###    Switch target file for conversion here       ###
#######################################################
### your file to convert - with v2 coordinates
map=fread('moms_teo20.hmp.txt', header=TRUE)
map <- data.table(map)
# Just selecting out the columns of interest
map <- subset(map, select=c("chrom", "pos"))

map <- data.frame(map)
names(map)=c('v2_chr', 'v2_coord')
### put into a genomic ranges object
## seqnames= will be your v2 chromosome
## start= will be your v2 coordinate
## end= will be your v2 coordinate again, unless you have a range to convert
m=GRanges(seqnames=as.factor(map$v2_chr), ranges=IRanges(start=map$v2_coord, end=map$v2_coord))


### find overlaps of v2 sites with the conversion
ov=findOverlaps(m, vc)
### storing as new column on df of original imported file
map$v3_chr=NA
map$v3_coord=NA
map$v3_chr[queryHits(ov)]=mcols(vc)$v3_chr[subjectHits(ov)] ## should always be the same between v2/v3

##############################
## Do the actual conversion ##
##############################
## v3_position = v3_window_start + (v2_coord - v2_window_start)
map$v3_coord[queryHits(ov)]=mcols(vc)$v3_start[subjectHits(ov)]+(map$v2_coord[queryHits(ov)]-start(vc)[subjectHits(ov)])


### because I need chromosome labels prefaced with 'chr'
map$v3_chr=paste('chr', map$v3_chr, sep='')


### now it's ready, v3 coordinates stored in map$v3_coord and v3 chromosome stored in map$v3_chr
### output in whatever format needed; here, outputting chromosome and v3 position
write.table(map, file='GBS_teos_moms.v3.txt', row.names=F, col.names=F, quote=F, sep='\t')

