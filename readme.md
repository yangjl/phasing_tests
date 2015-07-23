#Phasing Code

Uses selfed progeny from a single mom to 

1) impute mom's genotype
2) phase mom
3) impute self progeny genotypes

Currently assumes data is coded as 0: AA 1: Aa 2: aa and 3 for missing data.


------------------
# split JRI's codes and test independently

```
### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=T)
```


