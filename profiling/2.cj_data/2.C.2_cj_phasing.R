### Jinliang Yang
### August 7th, 2015

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
job=args[1]
infile <- read.table(infile, header=TRUE)[args[3], 1]
outfile <- read.table(infile, header=TRUE)[args[3], 2]


### loading all the functions in folder "lib"
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

ob <- load("largedata/sfamdata/PC_I11_ID2.RData")

ob <- load(infile)


win_length=10
verbose=TRUE
probs <- get_error_mat(0.02, 0.8)[[2]]

#########################
newmom_chr <- list()
for(i in 1:10){
    chr <- subset(mk, chr==i)
    
    progeny <- list()
    for(j in 6:ncol(chr)){
        progeny[[j-5]] <- list()
        progeny[[j-5]][[1]] <- chr[, j]
        progeny[[j-5]][[2]] <- chr[, j]
    }
    
    message(sprintf("###>>> phasing mom of [ %s ] kids for [ chr%s ]", ncol(chr)-5), i)
    newmom <- phasing(estimated_mom=chr[, 5], progeny, win_length, verbose)
    #plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)
    newmom_chr[[i]] <- newmom
}
save(file=outfile, list=c("newmom_chr"))
#########################
