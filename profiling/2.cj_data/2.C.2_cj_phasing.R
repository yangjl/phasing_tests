### Jinliang Yang
### August 7th, 2015

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
job <- args[2]

print(c(idx, job))

### loading all the functions in folder "lib"
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)


### start to phasing
files <- read.table("largedata/sfamdata/all_files.txt", header=TRUE)
infile <- as.character(files$input)[idx]
outfile <- as.character(files$output)[idx]
chrnum <- files$chr[idx]

win_length=10
verbose=TRUE
probs <- get_error_mat(0.02, 0.8)[[2]]

#########################
ob <- load(infile)
chr <- subset(mk, chr==chrnum)

progeny <- list()
for(j in 6:ncol(chr)){
    progeny[[j-5]] <- list()
    progeny[[j-5]][[1]] <- chr[, j]
    progeny[[j-5]][[2]] <- chr[, j]
}

message(sprintf("###>>> phasing mom of [ %s ] kids for [ chr%s ]", ncol(chr)-5, i))
newmom <- phasing(estimated_mom=chr[, 5], progeny, win_length, verbose)
#plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)
save(file=outfile, list=c("newmom"))
#########################
