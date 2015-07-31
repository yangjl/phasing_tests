### Jinliang Yang modified from VSB
### July 30th, 2015
## phasing.R

library(parallel)
library(devtools)
options(mc.cores=NULL)
load_all("~/bin/tasselr")
load_all("~/bin/ProgenyArray")


#### load h5file
teo <- initTasselHDF5(file="largedata/teo.h5", version="5")
teo <- loadBiallelicGenotypes(teo, verbose = TRUE)

#loading in genotypes from HDF5 file 'teo.h5'... done.
#binding samples together into matrix... done.
#filtering biallelic loci... done.
#encoding genotypes... done.
#Warning message:
#    In loadBiallelicGenotypes(teo, verbose = TRUE) :
#    Removed 357647 loci non-biallelic.

#[1:598043, 1:4875]
save(list="teo", file="largedata/cj_data.Rdata")


save(teo, file=tassel_file)








